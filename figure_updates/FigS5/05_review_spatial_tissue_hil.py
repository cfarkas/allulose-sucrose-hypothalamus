#!/usr/bin/env python3
"""DAPI-only HIL for the 3V contours displayed in S5 panels I/J.

This reviewer is separate from the accepted ARC/ME/VMN anatomy review. It
records exactly one closed 3V-cavity polygon per acquisition. The filled polygon
is a negative mask: nuclei and displayed density inside it are excluded from I/J
spatial analyses. It never defines the outer-tissue crop or normalization box.
Revised saves archive prior JSON/TIFF evidence instead of deleting it. The Make
Figure button creates a hash-bound receipt and an absent 600-dpi candidate
directory. The renderer remains backward-compatible with the frozen v1 review
whose region label was historically incorrect.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
import tifffile
from PIL import Image
from skimage.draw import polygon as draw_polygon
from skimage.segmentation import find_boundaries


SCHEMA = "figs5_spatial_3v_hil_v2"
LEGACY_SCHEMAS = {"figs4_spatial_tissue_hil_v1", "figs5_spatial_tissue_hil_v1"}
REGION = "3V"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path, block_size: int = 4 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def maximum_projection(path: Path) -> np.ndarray:
    array = np.squeeze(np.asarray(tifffile.imread(str(path))))
    while array.ndim > 2:
        array = np.max(array, axis=0)
    if array.ndim != 2:
        raise ValueError(f"DAPI image is not 2-D after projection: {path} {array.shape}")
    return array


def load_segmentation(path: Path) -> np.ndarray:
    raw = np.load(str(path), allow_pickle=True)
    if isinstance(raw, np.ndarray) and raw.dtype == object and raw.shape == ():
        raw = raw.item()
    if isinstance(raw, Mapping):
        raw = raw["masks"]
    array = np.squeeze(np.asarray(raw))
    while array.ndim > 2:
        array = array[0]
    if array.ndim != 2:
        raise ValueError(f"DAPI segmentation is not 2-D: {path} {array.shape}")
    return array.astype(np.int32, copy=False)


def normalize_u8(array: np.ndarray) -> np.ndarray:
    values = np.asarray(array, dtype=np.float32)
    finite = values[np.isfinite(values)]
    if not finite.size:
        return np.zeros(values.shape, dtype=np.uint8)
    low, high = np.percentile(finite, (0.2, 99.8))
    if not np.isfinite(high) or high <= low:
        high = low + 1.0
    return np.clip((values - low) / (high - low) * 255.0, 0, 255).astype(np.uint8)


def safe_replace_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    incoming = path.with_name(path.name + ".incoming")
    if incoming.exists():
        raise FileExistsError(f"Stale incoming file blocks safe write: {incoming}")
    incoming.write_bytes(payload)
    os.replace(incoming, path)


def archive_existing(path: Path, review_dir: Path) -> None:
    if not path.exists():
        return
    history = review_dir / "history"
    history.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    target = history / f"{path.stem}_{stamp}{path.suffix}"
    shutil.copy2(path, target)
    if sha256_file(path) != sha256_file(target):
        raise RuntimeError(f"History copy verification failed: {path}")


def load_sheet(workdir: Path) -> pd.DataFrame:
    path = workdir / "analysis_samplesheet.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Figure S5 analysis samplesheet missing: {path}")
    sheet = pd.read_csv(path)
    required = {
        "acquisition_id", "sample_id", "animal_id", "condition",
        "genotype", "dapi_mip", "dapi_seg",
    }
    missing = required - set(sheet.columns)
    if missing:
        raise ValueError(f"Samplesheet missing columns: {sorted(missing)}")
    expected = {
        "2026-07-09_FR4-1_Allulose", "2026-07-09_FR5-1_Allulose",
        "2026-07-09_FR5-2_Allulose", "2026-07-09_FR5-3_Sucrose",
        "2026-07-09_FR5-4_Sucrose", "2026-07-10_FR4-3_Sucrose",
        "2026-07-10_AGUA-M_Water", "2026-07-10_FR5-5_Water",
    }
    npy_acquisition = "2026-07-10_NPY-M_Water"
    if npy_acquisition in set(sheet["acquisition_id"]):
        expected.add(npy_acquisition)
        npy = sheet.loc[sheet["acquisition_id"].eq(npy_acquisition)]
        if len(npy) != 1 or npy.iloc[0]["genotype"] != "NPY-Tg" or npy.iloc[0]["condition"] != "Water":
            raise ValueError("NPY-M must use the author-confirmed NPY-Tg Water assignment")
    wt = sheet.loc[~sheet["acquisition_id"].eq(npy_acquisition)]
    if (set(sheet["acquisition_id"]) != expected or len(sheet) != len(expected)
            or sheet["animal_id"].astype(str).nunique() != len(expected)
            or len(wt) != 8 or not wt["genotype"].eq("WT").all()):
        raise ValueError("Spatial HIL requires the accepted eight WT animals, optionally plus author-confirmed NPY-M")
    order = pd.Categorical(
        sheet["condition"], ["Water", "Sucrose", "Allulose"], ordered=True
    )
    sheet = (
        sheet.assign(_order=order)
        .sort_values(["_order", "acquisition_id"])
        .drop(columns="_order")
        .reset_index(drop=True)
    )
    for index, row in sheet.iterrows():
        for column in ("dapi_mip", "dapi_seg"):
            raw = Path(str(row[column])).expanduser()
            candidates = [raw] if raw.is_absolute() else [workdir / raw, raw]
            for marker in ("prepared", "annotations", "region_masks", "results"):
                if marker in raw.parts:
                    candidates.append(workdir.joinpath(*raw.parts[raw.parts.index(marker):]))
                    break
            target = next((p.resolve() for p in candidates if p.is_file()), None)
            if target is None:
                raise FileNotFoundError(f"{row['acquisition_id']}: missing {column}: {raw}")
            sheet.at[index, column] = str(target)
    return sheet


def annotation_path(review_dir: Path, acquisition: str) -> Path:
    return review_dir / "annotations" / f"{acquisition}.json"


def ventricle_mask_path(review_dir: Path, acquisition: str) -> Path:
    return review_dir / "ventricle_masks" / f"{acquisition}_3v.tif"


def legacy_mask_path(review_dir: Path, acquisition: str) -> Path:
    return review_dir / "tissue_masks" / f"{acquisition}_outer_tissue.tif"


def read_annotation(review_dir: Path, acquisition: str) -> dict[str, Any] | None:
    path = annotation_path(review_dir, acquisition)
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def review_states(review_dir: Path, sheet: pd.DataFrame) -> list[str]:
    states: list[str] = []
    for acquisition in sheet["acquisition_id"].astype(str):
        record = read_annotation(review_dir, acquisition)
        semantics_valid = bool(
            record
            and (
                (
                    record.get("schema") == SCHEMA
                    and record.get("region") == REGION
                )
                or (
                    record.get("schema") in LEGACY_SCHEMAS
                    and record.get("region") == "OUTER_TISSUE"
                )
            )
        )
        states.append(
            "accepted"
            if semantics_valid and record.get("decision") == "include"
            else "pending"
        )
    return states


def rasterize_ventricle(
    points: list[list[float]], shape: tuple[int, int]
) -> np.ndarray:
    if len(points) < 3:
        raise ValueError("The 3V polygon needs at least three vertices")
    xy = np.asarray(points, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all():
        raise ValueError("Invalid 3V polygon coordinates")
    rr, cc = draw_polygon(xy[:, 1], xy[:, 0], shape=shape)
    mask = np.zeros(shape, dtype=np.uint8)
    mask[rr, cc] = 1
    if int(mask.sum()) < 100:
        raise ValueError("The 3V polygon is implausibly small")
    return mask


def save_decision(
    review_dir: Path,
    row: Mapping[str, Any],
    points: list[list[float]],
    shape: tuple[int, int],
    notes: str,
) -> None:
    if (review_dir / "SPATIAL_TISSUE_HIL_RECEIPT.json").exists():
        raise FileExistsError(
            "This HIL set is finalized; start a new review directory for revisions"
        )
    acquisition = str(row["acquisition_id"])
    mask = rasterize_ventricle(points, shape)
    ann_file = annotation_path(review_dir, acquisition)
    ventricle_file = ventricle_mask_path(review_dir, acquisition)
    archive_existing(ann_file, review_dir)
    archive_existing(ventricle_file, review_dir)
    record = {
        "schema": SCHEMA,
        "acquisition_id": acquisition,
        "sample_id": str(row["sample_id"]),
        "animal_id": str(row["animal_id"]),
        "condition": str(row["condition"]),
        "decision": "include",
        "region": REGION,
        "points": [
            [round(float(x), 3), round(float(y), 3)] for x, y in points
        ],
        "image_height": int(shape[0]),
        "image_width": int(shape[1]),
        "channels_used": ["DAPI"],
        "ventricular_contour_available": True,
        "used_for_tissue_geometry": False,
        "used_as_3v_exclusion_mask": True,
        "notes": str(notes),
        "saved_at_utc": utc_now(),
    }
    safe_replace_bytes(
        ann_file, (json.dumps(record, indent=2) + "\n").encode("utf-8")
    )
    ventricle_file.parent.mkdir(parents=True, exist_ok=True)
    incoming = ventricle_file.with_name(
        ventricle_file.name + ".incoming.tif"
    )
    if incoming.exists():
        raise FileExistsError(f"Stale incoming mask blocks safe write: {incoming}")
    tifffile.imwrite(str(incoming), mask, compression="zlib")
    os.replace(incoming, ventricle_file)


def finalize(review_dir: Path, sheet: pd.DataFrame) -> tuple[Path, Path]:
    states = review_states(review_dir, sheet)
    if states.count("pending"):
        raise RuntimeError(
            f"{states.count('pending')} spatial 3V decisions remain"
        )
    rows: list[dict[str, Any]] = []
    for _, row in sheet.iterrows():
        acquisition = str(row["acquisition_id"])
        ann_file = annotation_path(review_dir, acquisition)
        ventricle_file = ventricle_mask_path(review_dir, acquisition)
        record = json.loads(ann_file.read_text(encoding="utf-8"))
        mask = np.squeeze(
            np.asarray(tifffile.imread(str(ventricle_file)), dtype=np.uint8)
        )
        dapi = maximum_projection(Path(str(row["dapi_mip"])))
        if (
            mask.shape != dapi.shape
            or record.get("schema") != SCHEMA
            or record.get("region") != REGION
            or record.get("channels_used") != ["DAPI"]
            or record.get("ventricular_contour_available") is not True
            or record.get("used_for_tissue_geometry") is not False
            or record.get("used_as_3v_exclusion_mask") is not True
        ):
            raise ValueError(f"{acquisition}: HIL mask/DAPI validation failed")
        rows.append(
            {
                "acquisition_id": acquisition,
                "sample_id": str(row["sample_id"]),
                "animal_id": str(row["animal_id"]),
                "condition": str(row["condition"]),
                "decision": "include",
                "annotation_path": str(ann_file.relative_to(review_dir)),
                "annotation_sha256": sha256_file(ann_file),
                "ventricle_mask_path": str(ventricle_file.relative_to(review_dir)),
                "ventricle_mask_sha256": sha256_file(ventricle_file),
                "image_height": int(mask.shape[0]),
                "image_width": int(mask.shape[1]),
                "ventricle_pixels": int(np.count_nonzero(mask)),
                "source_dapi_mip": str(
                    Path(str(row["dapi_mip"])).resolve()
                ),
                "source_dapi_mip_sha256": sha256_file(
                    Path(str(row["dapi_mip"]))
                ),
                "source_dapi_segmentation": str(
                    Path(str(row["dapi_seg"])).resolve()
                ),
                "source_dapi_segmentation_sha256": sha256_file(
                    Path(str(row["dapi_seg"]))
                ),
                "channels_used": "DAPI only",
                "boundary_semantics": (
                    "one closed 3V negative mask; excludes I/J nuclei and density"
                ),
            }
        )
    manifest = review_dir / "SPATIAL_TISSUE_HIL_MANIFEST.csv"
    receipt = review_dir / "SPATIAL_TISSUE_HIL_RECEIPT.json"
    csv_payload = pd.DataFrame(rows).to_csv(index=False).encode("utf-8")
    receipt_payload = {
        "schema": SCHEMA,
        "status": "human_accepted",
        "figure": "S5",
        "panels": ["I", "J"],
        "accepted_acquisitions": len(rows),
        "pending_acquisitions": 0,
        "channels_used": ["DAPI"],
        "boundary_region": REGION,
        "one_ventricle_polygon_per_acquisition": True,
        "ventricular_contour_available": True,
        "used_for_tissue_geometry": False,
        "used_as_3v_exclusion_mask": True,
        "raw_images_modified": False,
        "raw_segmentations_modified": False,
        "manifest": manifest.name,
        "manifest_sha256": hashlib.sha256(csv_payload).hexdigest(),
        "finalized_at_utc": utc_now(),
    }
    receipt_bytes = (
        json.dumps(receipt_payload, indent=2) + "\n"
    ).encode("utf-8")
    if manifest.exists() or receipt.exists():
        raise FileExistsError(
            "This spatial HIL set is already finalized; its receipt is immutable"
        )
    safe_replace_bytes(manifest, csv_payload)
    safe_replace_bytes(receipt, receipt_bytes)
    return manifest, receipt


def display_image(
    row: Mapping[str, Any], mode: str, max_dim: int
) -> tuple[bytes, tuple[int, int], tuple[int, int]]:
    dapi = maximum_projection(Path(str(row["dapi_mip"])))
    shape = tuple(map(int, dapi.shape))
    if mode == "seg":
        segmentation = load_segmentation(Path(str(row["dapi_seg"])))
        if segmentation.shape != shape:
            raise ValueError("DAPI image/segmentation shape mismatch")
        rgb = np.zeros((*shape, 3), dtype=np.uint8)
        rgb[segmentation > 0] = (120, 145, 220)
        rgb[find_boundaries(segmentation, mode="inner")] = (220, 230, 255)
    else:
        u8 = normalize_u8(dapi)
        rgb = np.stack((u8 // 5, u8 // 3, u8), axis=-1)
    scale = min(1.0, float(max_dim) / max(shape))
    display = (
        max(1, int(round(shape[0] * scale))),
        max(1, int(round(shape[1] * scale))),
    )
    image = Image.fromarray(rgb)
    if display != shape:
        image = image.resize(
            (display[1], display[0]), Image.Resampling.LANCZOS
        )
    buffer = io.BytesIO()
    image.save(buffer, format="PNG", optimize=True)
    return buffer.getvalue(), shape, display


HTML = r'''<!doctype html><html><head><meta charset="utf-8"><title>S5 I/J DAPI 3V HIL</title><style>
:root{--ink:#162033;--paper:#f4f7fa;--ok:#167d50;--bad:#a02837;--build:#6650a4}*{box-sizing:border-box}body{margin:0;font:14px Arial;color:#fff;background:#000}button,select,textarea{font:inherit;font-weight:700}button{cursor:pointer;border:0;border-radius:4px;padding:8px 10px}button:disabled{opacity:.55;cursor:not-allowed}#top{height:58px;background:#111;display:flex;align-items:center;gap:8px;padding:6px 10px;border-bottom:2px solid #333}#progress{font-weight:800;white-space:nowrap}#meta{flex:1;font-weight:700;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.build{background:var(--build);color:#fff;border:2px solid #b9abef}#main{display:grid;grid-template-columns:minmax(0,1fr) 340px;height:calc(100vh - 58px)}#wrap{overflow:hidden}svg{width:100%;height:100%;cursor:crosshair;touch-action:none}#side{background:var(--paper);color:var(--ink);padding:13px;overflow:auto}.section{padding-bottom:12px;margin-bottom:12px;border-bottom:1px solid #cad2da}.small{font-size:12px;color:#4a5966;line-height:1.35}.action{width:100%;margin:5px 0}.save{background:var(--ok);color:white}.clear{background:var(--bad);color:white}textarea{width:100%;min-height:55px}#buildStatus{overflow-wrap:anywhere}
</style></head><body><div id="top"><button onclick="loadItem(index-1)">◀</button><button onclick="loadItem(index+1)">▶</button><span id="progress">Loading…</span><select id="mode" onchange="reloadImage()"><option value="raw" selected>Raw DAPI</option><option value="seg">DAPI segmentation</option></select><button onclick="resetView()">Reset view</button><div id="meta"></div><button id="buildTop" class="build" onclick="buildFigure()" disabled>Make S5 candidate · 600 dpi</button></div>
<div id="main"><div id="wrap"><svg id="viewer"><g id="viewport"><image id="bg" x="0" y="0"></image><polygon id="saved" fill="#17d8ff" fill-opacity=".12" stroke="#fff" stroke-width="7" stroke-dasharray="18 12" stroke-linejoin="round" vector-effect="non-scaling-stroke"></polygon><polyline id="working" fill="none" stroke="#17d8ff" stroke-width="7" stroke-linejoin="round" vector-effect="non-scaling-stroke"></polyline><g id="dots"></g></g></svg></div><div id="side">
<div class="section"><b>DAPI-only 3V exclusion review for panels I/J</b><p class="small">Trace exactly one closed polygon along the visible 3V cavity boundary. Nuclei and displayed density inside it are excluded from I/J. It never defines the outer-tissue crop or normalization box.</p></div>
<div class="section"><b>Status:</b> <span id="state"></span><p class="small">Auto-close: after ≥3 vertices, click near the first white dot. Enter or double-click also closes. Shift, middle, or right drag pans; wheel zooms.</p></div>
<button class="action" onclick="finish()">Close polygon</button><button class="action" onclick="undo()">Undo point</button><button class="action clear" onclick="clearBoundary()">Clear boundary</button><textarea id="notes" placeholder="Review notes"></textarea><button class="action save" onclick="saveBoundary(false)">Save accepted</button><button class="action save" onclick="saveBoundary(true)">Save & next</button>
<div class="section"><button id="buildSide" class="action build" onclick="buildFigure()" disabled>Make S5 candidate · 600 dpi</button><div id="buildStatus" class="small">Loading review state…</div></div></div></div>
<script>'use strict';let items=[],index=0,meta=null,points=[],closed=[],zoom=1,panX=0,panY=0,drag=false,lx=0,ly=0,canBuild=false,building=false;const svg=document.getElementById('viewer'),vp=document.getElementById('viewport');
function buttons(){for(const id of ['buildTop','buildSide'])document.getElementById(id).disabled=!(canBuild&&!building)}function update(c){canBuild=c.pending===0;buttons();document.getElementById('buildStatus').textContent=canBuild?'All '+c.accepted+' 3V contours accepted. Ready to render.':c.pending+' DAPI 3V decision'+(c.pending===1?'':'s')+' remain.'}
function resetView(){zoom=1;panX=0;panY=0;transform()}function transform(){vp.setAttribute('transform','translate('+panX+' '+panY+') scale('+zoom+')');render()}function svgPoint(e){let p=svg.createSVGPoint();p.x=e.clientX;p.y=e.clientY;return p.matrixTransform(vp.getScreenCTM().inverse())}function screenPoint(p){let q=svg.createSVGPoint();q.x=p[0];q.y=p[1];return q.matrixTransform(vp.getScreenCTM())}function nearFirst(e){if(points.length<3)return false;let q=screenPoint(points[0]);return Math.hypot(e.clientX-q.x,e.clientY-q.y)<=16}function finish(){if(points.length<3){document.getElementById('state').textContent='At least three vertices are required.';return false}closed=points.slice();points=[];render();return true}function undo(){if(points.length)points.pop();else if(closed.length){points=closed.slice();closed=[];points.pop()}render()}function clearBoundary(){if(confirm('Clear this 3V contour?')){points=[];closed=[];render()}}
function dot(p){let c=document.createElementNS('http://www.w3.org/2000/svg','circle');c.setAttribute('cx',p[0]);c.setAttribute('cy',p[1]);c.setAttribute('r','5');c.setAttribute('fill','#fff');c.setAttribute('stroke','#17d8ff');c.setAttribute('stroke-width','2.5');c.setAttribute('vector-effect','non-scaling-stroke');document.getElementById('dots').appendChild(c)}function render(){document.getElementById('saved').setAttribute('points',closed.map(function(p){return p.join(',')}).join(' '));document.getElementById('working').setAttribute('points',points.map(function(p){return p.join(',')}).join(' '));let d=document.getElementById('dots');d.innerHTML='';for(const p of points)dot(p)}
svg.addEventListener('contextmenu',function(e){e.preventDefault()});svg.addEventListener('wheel',function(e){e.preventDefault();let p=svgPoint(e),f=e.deltaY<0?1.15:1/1.15,o=zoom;zoom=Math.min(15,Math.max(.25,zoom*f));panX+=p.x*(o-zoom);panY+=p.y*(o-zoom);transform()},{passive:false});svg.addEventListener('pointerdown',function(e){if(e.button===1||e.button===2||e.shiftKey){drag=true;lx=e.clientX;ly=e.clientY;svg.setPointerCapture(e.pointerId);return}if(e.button!==0||closed.length)return;if(nearFirst(e)){finish();return}let p=svgPoint(e);if(meta&&p.x>=0&&p.y>=0&&p.x<meta.display_width&&p.y<meta.display_height){points.push([p.x,p.y]);render()}});svg.addEventListener('pointermove',function(e){if(drag){let r=svg.getBoundingClientRect();panX+=(e.clientX-lx)*(meta.display_width/r.width)/zoom;panY+=(e.clientY-ly)*(meta.display_height/r.height)/zoom;lx=e.clientX;ly=e.clientY;transform()}});svg.addEventListener('pointerup',function(){drag=false});svg.addEventListener('dblclick',function(e){e.preventDefault();finish()});document.addEventListener('keydown',function(e){if(e.key==='Enter'){finish();e.preventDefault()}if(e.key==='Backspace'){undo();e.preventDefault()}if(e.key==='Escape'){points=[];render()}});
async function init(){let r=await fetch('/api/list',{cache:'no-store'});items=await r.json();let p=items.findIndex(function(x){return x.status==='pending'});await loadItem(p>=0?p:0)}async function loadItem(i){index=Math.max(0,Math.min(items.length-1,i));let r=await fetch('/api/item/'+index,{cache:'no-store'});meta=await r.json();closed=meta.points||[];points=[];document.getElementById('notes').value=meta.notes||'';document.getElementById('state').textContent=meta.status;document.getElementById('progress').textContent=(index+1)+'/'+items.length+' | accepted '+meta.counts.accepted+', pending '+meta.counts.pending;document.getElementById('meta').textContent=meta.acquisition_id+' | '+meta.condition;update(meta.counts);svg.setAttribute('viewBox','0 0 '+meta.display_width+' '+meta.display_height);let b=document.getElementById('bg');b.setAttribute('width',meta.display_width);b.setAttribute('height',meta.display_height);reloadImage();resetView()}function reloadImage(){document.getElementById('bg').setAttribute('href','/api/image/'+index+'?mode='+document.getElementById('mode').value+'&t='+Date.now())}
async function saveBoundary(next){if(!closed.length&&!finish())return;let body={points:closed,notes:document.getElementById('notes').value,display_width:meta.display_width,display_height:meta.display_height};let r=await fetch('/api/save/'+index,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});let t=await r.text();if(!r.ok){alert(t);return}items=await(await fetch('/api/list',{cache:'no-store'})).json();await loadItem(next?Math.min(index+1,items.length-1):index)}async function buildFigure(){building=true;buttons();let s=document.getElementById('buildStatus');s.textContent='Hashing accepted DAPI 3V HIL and rendering English master plus bilingual subpanels at 600 dpi…';let r=await fetch('/api/build',{method:'POST'}),t=await r.text();if(!r.ok){s.textContent=t;building=false;buttons();return}let d=JSON.parse(t);s.textContent='Built: '+d.output_dir;building=false;buttons()}init().catch(function(e){document.getElementById('buildStatus').textContent=e.message});</script></body></html>'''


def serve(args: argparse.Namespace) -> int:
    from flask import Flask, Response, jsonify, request

    workdir = args.workdir.expanduser().resolve()
    review_dir = args.review_dir.expanduser().resolve()
    review_dir.mkdir(parents=True, exist_ok=True)
    (review_dir / "annotations").mkdir(exist_ok=True)
    if not (review_dir / "SPATIAL_TISSUE_HIL_RECEIPT.json").exists():
        (review_dir / "ventricle_masks").mkdir(exist_ok=True)
    sheet = load_sheet(workdir)
    records = sheet.to_dict("records")
    cache: dict[
        tuple[int, str],
        tuple[bytes, tuple[int, int], tuple[int, int]],
    ] = {}
    app = Flask(__name__)

    def rendered(i: int, mode: str):
        key = (i, mode)
        if key not in cache:
            cache[key] = display_image(records[i], mode, args.web_max_dim)
        return cache[key]

    @app.get("/")
    def root_route():
        response = Response(HTML, mimetype="text/html")
        response.headers["Cache-Control"] = "no-store"
        return response

    @app.get("/api/list")
    def list_route():
        states = review_states(review_dir, sheet)
        return jsonify(
            [
                {
                    "index": i,
                    "acquisition_id": str(row["acquisition_id"]),
                    "status": states[i],
                }
                for i, row in enumerate(records)
            ]
        )

    @app.get("/api/item/<int:i>")
    def item_route(i: int):
        if i < 0 or i >= len(records):
            return jsonify({"error": "index"}), 404
        _, native, display = rendered(i, "raw")
        row = records[i]
        record = read_annotation(
            review_dir, str(row["acquisition_id"])
        ) or {}
        sx, sy = display[1] / native[1], display[0] / native[0]
        points = [
            [float(x) * sx, float(y) * sy]
            for x, y in record.get("points", [])
        ]
        states = review_states(review_dir, sheet)
        counts = {
            state: states.count(state) for state in ("accepted", "pending")
        }
        return jsonify(
            {
                "acquisition_id": str(row["acquisition_id"]),
                "condition": str(row["condition"]),
                "status": states[i],
                "notes": record.get("notes", ""),
                "points": points,
                "display_width": display[1],
                "display_height": display[0],
                "counts": counts,
            }
        )

    @app.get("/api/image/<int:i>")
    def image_route(i: int):
        if i < 0 or i >= len(records):
            return Response(status=404)
        mode = request.args.get("mode", "raw")
        mode = mode if mode in {"raw", "seg"} else "raw"
        return Response(rendered(i, mode)[0], mimetype="image/png")

    @app.post("/api/save/<int:i>")
    def save_route(i: int):
        if i < 0 or i >= len(records):
            return jsonify({"error": "index"}), 404
        incoming = request.get_json(force=True)
        _, native, display = rendered(i, "raw")
        sx = native[1] / max(
            1.0, float(incoming.get("display_width", display[1]))
        )
        sy = native[0] / max(
            1.0, float(incoming.get("display_height", display[0]))
        )
        points = [
            [float(x) * sx, float(y) * sy]
            for x, y in incoming.get("points", [])
        ]
        try:
            save_decision(
                review_dir,
                records[i],
                points,
                native,
                str(incoming.get("notes", "")),
            )
        except Exception as exc:
            return jsonify({"error": str(exc)}), 409
        return jsonify({"ok": True})

    @app.post("/api/build")
    def build_route():
        try:
            receipt = review_dir / "SPATIAL_TISSUE_HIL_RECEIPT.json"
            if not receipt.exists():
                _, receipt = finalize(review_dir, sheet)

            stamp = datetime.now(timezone.utc).strftime(
                "%Y%m%dT%H%M%S%fZ"
            )
            output_dir = (
                args.output_root.expanduser().resolve()
                / f"spatial_hil_candidate_{stamp}"
            )
            command = [
                sys.executable,
                str(args.renderer.expanduser().resolve()),
                "--workdir",
                str(workdir),
                "--spatial-hil-dir",
                str(review_dir),
                "--output-dir",
                str(output_dir),
                "--dpi",
                "600",
            ]
            completed = subprocess.run(
                command,
                text=True,
                capture_output=True,
                timeout=1800,
            )
            if completed.returncode:
                return jsonify(
                    {
                        "error": "S5 render failed",
                        "log": (
                            completed.stdout + completed.stderr
                        )[-12000:],
                    }
                ), 500
            return jsonify(
                {
                    "ok": True,
                    "output_dir": str(output_dir),
                    "receipt": str(receipt),
                    "figure_png": str(output_dir / "Figure_S5.png"),
                }
            )
        except Exception as exc:
            return jsonify({"error": str(exc)}), 409

    @app.post("/api/stop")
    def stop_route():
        def shutdown():
            time.sleep(0.3)
            os.kill(os.getpid(), signal.SIGINT)

        threading.Thread(target=shutdown, daemon=True).start()
        return jsonify({"ok": True})

    print(
        f"S5 panels I/J DAPI-only 3V HIL: "
        f"http://{args.host}:{args.port}",
        flush=True,
    )
    print(f"Review directory: {review_dir}", flush=True)
    app.run(
        host=args.host,
        port=args.port,
        threaded=True,
        use_reloader=False,
    )
    return 0


def main() -> int:
    paper = next(p for p in Path(__file__).resolve().parents if (p / "FigS5").is_dir() and (p / "analyses").is_dir())
    parser = argparse.ArgumentParser(
        description=(
            "DAPI-only 3V exclusion HIL for S5 spatial panels I/J"
        )
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=(
            paper
            / "analyses/FigS5/results/including_NPY_M_20260908"
        ),
    )
    parser.add_argument(
        "--review-dir",
        type=Path,
        default=(
            paper / "FigS5/hil_review/spatial_3v_including_NPY_M_20260908"
        ),
    )
    parser.add_argument(
        "--renderer",
        type=Path,
        default=paper / "FigS5/03_make_figure_s5_acth_clip_cfos.py",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=paper / "analyses/FigS5/figure",
    )
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=34250)
    parser.add_argument("--web-max-dim", type=int, default=1700)
    parser.add_argument("--status", action="store_true")
    args = parser.parse_args()
    sheet = load_sheet(args.workdir.expanduser().resolve())
    if args.status:
        states = review_states(
            args.review_dir.expanduser().resolve(), sheet
        )
        print(
            json.dumps(
                {
                    "accepted": states.count("accepted"),
                    "pending": states.count("pending"),
                    "review_dir": str(
                        args.review_dir.expanduser().resolve()
                    ),
                },
                indent=2,
            )
        )
        return 0
    return serve(args)


if __name__ == "__main__":
    raise SystemExit(main())
