#!/usr/bin/env python3
"""Local human-in-the-loop microglia morphology review server.

This is the second reviewer in the Figure 5 workflow and it is not the same job
as the first. scripts/shared/01_annotate_regions.py delineates ARC, ME and VMN
on whole sections; this one classifies individual Iba1 cells into the four
morphological states, which is what Fig5/README.txt gate 4 asks for.

The analyzer already extracts every cell as a 96x96 patch and writes
microglia_annotation_template.csv with a model proposal in ``pseudo_state`` and
an empty ``state`` column. This serves those patches one at a time with the
proposal preselected, so the human is reviewing a suggestion rather than
labelling from nothing, and writes the accepted labels to a CSV that
Fig5/01_analyze_gfap_iba1_microglia.py consumes through --microglia-labels-csv.

A proposal is never counted as a label. Only a state a human committed is
written, every row records who committed it and when, and the patch SHA-256 is
stored so a label can always be tied back to the exact image that was shown.

The page keeps a live readout of the supervised training design that
--microglia-classifier-mode supervised enforces, so you can see when the gate is
actually satisfied instead of finding out when the analyzer refuses.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import ipaddress
import json
import os
import random
import re
import sys
import threading
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

HERE = Path(__file__).resolve()
PAPER_ROOT = next((path for path in HERE.parents
                   if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                       and (path / "README.txt").is_file()))), None)
if PAPER_ROOT is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")

TEMPLATE_NAME = "microglia_annotation_template.csv"
DEFAULT_OUTPUT_NAME = "microglia_reviewed_labels.csv"
# Kept in step with MICROGLIA_STATES in Fig5/01_analyze_gfap_iba1_microglia.py.
STATES = ("Ramified", "Rod-like", "Activated", "Amoeboid")
STATE_HELP = {
    "Ramified": "small soma, long thin branched processes, surveillant",
    "Rod-like": "elongated bipolar soma, few processes, often aligned",
    "Activated": "thickened retracted processes, enlarged soma, bushy",
    "Amoeboid": "round soma, few or no processes, phagocytic",
}
# The supervised training design Fig5/01 enforces before it will train on
# reviewed labels. Mirrored here so the page can show progress toward it.
GATE = {"min_cells": 120, "min_per_class": 12, "min_classes": 2, "min_animals": 3}
OUTPUT_FIELDS = (
    "cell_uid", "state", "reviewer", "session", "reviewed_at", "patch_sha256",
    "animal_id", "section", "region", "pseudo_state", "pseudo_confidence",
    "agreed_with_model",
)
MAX_BODY = 1 * 1024 * 1024


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class Review:
    """The template, the accepted labels, and the order cells are served in."""

    def __init__(self, analysis_dir: Path, output_csv: Path, seed: int,
                 only_animals=None, only_regions=None):
        self.analysis_dir = analysis_dir.resolve()
        self.output_csv = output_csv.resolve()
        template = self.analysis_dir / TEMPLATE_NAME
        if not template.is_file():
            raise SystemExit(
                f"No {TEMPLATE_NAME} in {self.analysis_dir}. Run "
                "Fig5/01_analyze_gfap_iba1_microglia.py first; --prepare-labels-only "
                "writes the template without training or running statistics."
            )
        with template.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            raise SystemExit(f"{template} is empty")
        for row in rows:
            row["patch_path"] = str(row.get("patch_path", "") or "")
        if only_animals:
            keep = {a.strip() for a in only_animals if a.strip()}
            rows = [r for r in rows if str(r.get("animal_id", "")) in keep]
        if only_regions:
            keep = {r.strip().upper() for r in only_regions if r.strip()}
            rows = [r for r in rows if str(r.get("region", "")).upper() in keep]
        rows = [r for r in rows if r["patch_path"] and Path(r["patch_path"]).is_file()]
        if not rows:
            raise SystemExit(
                "No reviewable cells: every template row is filtered out or its "
                "patch file is missing. Patches live under "
                f"{self.analysis_dir / 'microglia_patches'}."
            )
        self.rows = rows
        self.by_uid = {str(r["cell_uid"]): r for r in rows}
        self.lock = threading.Lock()
        self.labels: dict[str, dict] = {}
        self._load_existing()
        self.order = self._review_order(seed)

    # -- persistence --------------------------------------------------------
    def _load_existing(self) -> None:
        if not self.output_csv.is_file():
            return
        with self.output_csv.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                uid = str(row.get("cell_uid", ""))
                if uid in self.by_uid and str(row.get("state", "")) in STATES:
                    self.labels[uid] = dict(row)

    def _write(self) -> None:
        """Rewrite the whole file atomically; a few thousand rows is nothing."""
        self.output_csv.parent.mkdir(parents=True, exist_ok=True)
        temp = self.output_csv.with_suffix(self.output_csv.suffix + ".tmp")
        with temp.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(OUTPUT_FIELDS))
            writer.writeheader()
            for uid in sorted(self.labels):
                writer.writerow({k: self.labels[uid].get(k, "") for k in OUTPUT_FIELDS})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp, self.output_csv)

    # -- ordering -----------------------------------------------------------
    def _review_order(self, seed: int) -> list:
        """Round-robin over proposed class, then animal, so the gate is reachable.

        Reviewing the template top to bottom means hundreds of Ramified cells
        from one animal before anything else, and the supervised gate needs at
        least twelve cells in two classes across three animals. Interleaving gets
        there in a fraction of the clicks. It is seeded, so the order a reviewer
        sees is reproducible.
        """
        rng = random.Random(seed)
        buckets: dict[tuple, list] = {}
        for row in self.rows:
            key = (str(row.get("pseudo_state", "")), str(row.get("animal_id", "")))
            buckets.setdefault(key, []).append(str(row["cell_uid"]))
        for uids in buckets.values():
            rng.shuffle(uids)
        keys = sorted(buckets)
        rng.shuffle(keys)
        order: list = []
        while any(buckets[k] for k in keys):
            for key in keys:
                if buckets[key]:
                    order.append(buckets[key].pop())
        return order

    # -- state --------------------------------------------------------------
    def counts(self) -> dict:
        per_class = {s: 0 for s in STATES}
        animals = set()
        agreed = 0
        for row in self.labels.values():
            state = str(row.get("state", ""))
            if state in per_class:
                per_class[state] += 1
            animals.add(str(row.get("animal_id", "")))
            if str(row.get("agreed_with_model", "")).lower() == "true":
                agreed += 1
        total = sum(per_class.values())
        usable = [s for s, n in per_class.items() if n >= GATE["min_per_class"]]
        return {
            "labelled": total, "total_cells": len(self.rows),
            "per_class": per_class, "animals": len(animals),
            "usable_classes": len(usable),
            "agreed_with_model": agreed,
            "gate": dict(GATE),
            "gate_met": (
                total >= GATE["min_cells"]
                and len(usable) >= GATE["min_classes"]
                and len(animals) >= GATE["min_animals"]
            ),
            "output_csv": str(self.output_csv),
        }

    def next_uid(self, after: str = "") -> str:
        with self.lock:
            pending = [u for u in self.order if u not in self.labels]
            if not pending:
                return ""
            if after and after in pending:
                index = pending.index(after)
                return pending[(index + 1) % len(pending)]
            return pending[0]

    def cell(self, uid: str) -> dict:
        row = self.by_uid.get(uid)
        if row is None:
            raise KeyError(f"unknown cell_uid: {uid}")
        existing = self.labels.get(uid, {})
        return {
            "cell_uid": uid,
            "animal_id": row.get("animal_id", ""), "section": row.get("section", ""),
            "region": row.get("region", ""), "sample": row.get("sample", ""),
            "pseudo_state": row.get("pseudo_state", ""),
            "pseudo_confidence": row.get("pseudo_confidence", ""),
            "state": existing.get("state", ""),
            "reviewed_at": existing.get("reviewed_at", ""),
            "reviewer": existing.get("reviewer", ""),
        }

    def save(self, payload: dict) -> dict:
        uid = str(payload.get("cell_uid", ""))
        state = str(payload.get("state", ""))
        reviewer = str(payload.get("reviewer", "")).strip()
        session = str(payload.get("session", "")).strip()
        if uid not in self.by_uid:
            raise ValueError(f"unknown cell_uid: {uid}")
        if state not in STATES:
            raise ValueError(f"state must be one of {', '.join(STATES)}")
        if not reviewer or not session:
            raise ValueError("reviewer name and session identifier are required")
        if len(reviewer) > 200 or len(session) > 200:
            raise ValueError("reviewer/session must be 200 characters or fewer")
        row = self.by_uid[uid]
        patch = Path(row["patch_path"])
        if not patch.is_file():
            raise ValueError(f"patch image is missing: {patch}")
        now = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
        with self.lock:
            self.labels[uid] = {
                "cell_uid": uid, "state": state, "reviewer": reviewer,
                "session": session, "reviewed_at": now,
                "patch_sha256": sha256_file(patch),
                "animal_id": row.get("animal_id", ""), "section": row.get("section", ""),
                "region": row.get("region", ""),
                "pseudo_state": row.get("pseudo_state", ""),
                "pseudo_confidence": row.get("pseudo_confidence", ""),
                "agreed_with_model": str(state == str(row.get("pseudo_state", ""))),
            }
            self._write()
        return {"ok": True, "cell_uid": uid, "state": state, **self.counts()}

    def clear(self, uid: str) -> dict:
        with self.lock:
            self.labels.pop(uid, None)
            self._write()
        return {"ok": True, "cell_uid": uid, **self.counts()}

    def patch_bytes(self, uid: str) -> bytes:
        row = self.by_uid.get(uid)
        if row is None:
            raise KeyError(f"unknown cell_uid: {uid}")
        patch = Path(row["patch_path"])
        if not patch.is_file():
            raise KeyError(f"patch image is missing: {patch}")
        return patch.read_bytes()


class Handler(BaseHTTPRequestHandler):
    server_version = "ApotomeMicrogliaReview/1.0"

    @property
    def review(self) -> Review:
        return self.server.review  # type: ignore[attr-defined]

    def log_message(self, *args):  # quiet
        return

    def send_json(self, payload, status=HTTPStatus.OK):
        body = json.dumps(payload).encode()
        self.send_response(int(status))
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def send_bytes(self, body, content_type, status=HTTPStatus.OK):
        self.send_response(int(status))
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        if body:
            self.wfile.write(body)

    def error_json(self, status, message):
        self.send_json({"ok": False, "error": message}, status)

    def do_GET(self):
        parsed = urlparse(self.path)
        query = parse_qs(parsed.query)
        try:
            if parsed.path == "/":
                self.send_bytes(HTML.encode(), "text/html; charset=utf-8")
            elif parsed.path in ("/health", "/api/health"):
                self.send_json({"ok": True, "service": "microglia-review",
                                **self.review.counts()})
            elif parsed.path == "/api/config":
                self.send_json({"ok": True, "states": list(STATES),
                                "state_help": STATE_HELP,
                                "analysis_dir": str(self.review.analysis_dir),
                                **self.review.counts()})
            elif parsed.path == "/api/next":
                uid = self.review.next_uid(query.get("after", [""])[0])
                if not uid:
                    self.send_json({"ok": True, "done": True, **self.review.counts()})
                else:
                    self.send_json({"ok": True, "done": False,
                                    "cell": self.review.cell(uid),
                                    **self.review.counts()})
            elif parsed.path == "/api/cell":
                self.send_json({"ok": True,
                                "cell": self.review.cell(query.get("cell_uid", [""])[0]),
                                **self.review.counts()})
            elif parsed.path == "/api/patch":
                self.send_bytes(self.review.patch_bytes(query.get("cell_uid", [""])[0]),
                                "image/png")
            elif parsed.path == "/favicon.ico":
                self.send_bytes(b"", "image/x-icon", HTTPStatus.NO_CONTENT)
            else:
                self.error_json(HTTPStatus.NOT_FOUND, "endpoint not found")
        except KeyError as exc:
            self.error_json(HTTPStatus.NOT_FOUND, str(exc))
        except Exception as exc:  # pragma: no cover
            self.error_json(HTTPStatus.INTERNAL_SERVER_ERROR, str(exc))

    def do_POST(self):
        endpoint = urlparse(self.path).path
        if endpoint not in ("/api/label", "/api/clear"):
            self.error_json(HTTPStatus.NOT_FOUND, "endpoint not found")
            return
        try:
            if self.headers.get_content_type() != "application/json":
                self.error_json(HTTPStatus.UNSUPPORTED_MEDIA_TYPE,
                                "Content-Type must be application/json")
                return
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > MAX_BODY:
                self.error_json(HTTPStatus.REQUEST_ENTITY_TOO_LARGE,
                                "invalid or oversized request body")
                return
            payload = json.loads(self.rfile.read(length).decode())
            if endpoint == "/api/label":
                self.send_json(self.review.save(payload))
            else:
                self.send_json(self.review.clear(str(payload.get("cell_uid", ""))))
        except (ValueError, TypeError, json.JSONDecodeError) as exc:
            self.error_json(HTTPStatus.UNPROCESSABLE_ENTITY, str(exc))
        except Exception as exc:  # pragma: no cover
            self.error_json(HTTPStatus.INTERNAL_SERVER_ERROR, str(exc))


class Server(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, address, review):
        self.review = review
        super().__init__(address, Handler)


HTML = r'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Microglia morphology review</title><style>
:root{--bg:#12161c;--panel:#1b212b;--line:#2b3441;--ink:#e8edf4;--dim:#93a1b3;--ok:#2ecc71;--warn:#e67e22}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,sans-serif;display:flex;height:100vh}
main{flex:1;display:flex;flex-direction:column;align-items:center;justify-content:center;gap:14px;padding:20px}
aside{width:330px;background:var(--panel);border-left:1px solid var(--line);padding:16px;overflow:auto}
h1{font-size:15px;margin:0 0 12px}
#patch{image-rendering:pixelated;width:384px;height:384px;background:#000;border:1px solid var(--line);border-radius:6px}
.meta{color:var(--dim);font-size:13px;text-align:center}
.states{display:grid;grid-template-columns:1fr 1fr;gap:8px;width:520px}
button.state{padding:11px 8px;border-radius:6px;border:1px solid var(--line);background:#232c38;color:var(--ink);cursor:pointer;text-align:left}
button.state:hover{border-color:#4a5a6f}
button.state .k{display:inline-block;width:18px;color:var(--dim)}
button.state.proposed{border-color:var(--warn);box-shadow:0 0 0 1px var(--warn) inset}
button.state .d{display:block;color:var(--dim);font-size:11px;margin-left:18px}
.row{display:flex;gap:8px;width:520px}
.row button{flex:1;padding:8px;border-radius:6px;border:1px solid var(--line);background:#232c38;color:var(--dim);cursor:pointer}
.group{margin-bottom:16px;padding-bottom:14px;border-bottom:1px solid var(--line)}
label{display:block;color:var(--dim);font-size:12px;margin-bottom:4px}
input{width:100%;padding:7px;border-radius:5px;border:1px solid var(--line);background:#141a22;color:var(--ink)}
.bar{height:6px;background:#141a22;border-radius:3px;overflow:hidden;margin-top:4px}
.bar i{display:block;height:100%;background:var(--ok)}
table{width:100%;border-collapse:collapse;font-size:12px}td{padding:2px 0}td:last-child{text-align:right;color:var(--dim)}
.gate{padding:9px;border-radius:6px;font-size:12px;background:#241d13;border:1px solid var(--warn)}
.gate.met{background:#12291b;border-color:var(--ok)}
.status{font-size:12px;color:var(--dim);min-height:2.2em}
</style></head><body>
<main>
<img id="patch" alt="microglia cell patch">
<div class="meta" id="meta">loading</div>
<div class="states" id="states"></div>
<div class="row"><button id="skip">Skip (s)</button><button id="unlabel">Clear this label</button></div>
<div class="status" id="status"></div>
</main>
<aside>
<h1>Microglia morphology review</h1>
<div class="group">
<label>Reviewer name</label><input id="reviewer" placeholder="required">
<label style="margin-top:8px">Session identifier</label><input id="session" placeholder="required">
</div>
<div class="group">
<div id="progresstext" class="meta" style="text-align:left"></div>
<div class="bar"><i id="progressbar" style="width:0%"></i></div>
<table id="counts"></table>
</div>
<div class="group"><div id="gate" class="gate">gate</div></div>
<div class="group"><div class="meta" style="text-align:left" id="outpath"></div></div>
</aside>
<script>'use strict';
const $=id=>document.getElementById(id);
let cfg={states:[],state_help:{}},cur=null,busy=false;
function ident(){const r=$('reviewer').value.trim(),s=$('session').value.trim();
 if(!r||!s){status('Fill in reviewer name and session identifier first.');return null}return{reviewer:r,session:s}}
function status(m){$('status').textContent=m||''}
function renderStates(){$('states').innerHTML='';cfg.states.forEach((s,i)=>{
 const b=document.createElement('button');b.className='state';b.dataset.state=s;
 b.innerHTML=`<span class="k">${i+1}</span><b>${s}</b><span class="d">${cfg.state_help[s]||''}</span>`;
 b.onclick=()=>label(s);$('states').append(b)})}
function paint(c){cur=c;
 $('patch').src=`/api/patch?cell_uid=${encodeURIComponent(c.cell_uid)}&t=${Date.now()}`;
 const conf=c.pseudo_confidence?` (${Number(c.pseudo_confidence).toFixed(2)})`:'';
 $('meta').textContent=`${c.animal_id} · ${c.section} · ${c.region} · ${c.cell_uid}  —  model proposes ${c.pseudo_state}${conf}`;
 document.querySelectorAll('button.state').forEach(b=>{
  b.classList.toggle('proposed',b.dataset.state===c.pseudo_state)});
 status(c.state?`Currently labelled ${c.state} by ${c.reviewer}.`:'')}
function stats(d){
 $('progresstext').textContent=`${d.labelled} labelled of ${d.total_cells} cells · ${d.animals} animals`;
 $('progressbar').style.width=Math.min(100,100*d.labelled/Math.max(1,d.gate.min_cells))+'%';
 $('counts').innerHTML=Object.entries(d.per_class).map(([k,v])=>
  `<tr><td>${k}</td><td>${v}</td></tr>`).join('')+
  `<tr><td>agreed with model</td><td>${d.agreed_with_model}</td></tr>`;
 const g=d.gate,ok=d.gate_met;const el=$('gate');el.className='gate'+(ok?' met':'');
 el.innerHTML=ok?`<b>Supervised gate met.</b><br>Run the analyzer with --microglia-classifier-mode supervised and --microglia-labels-csv.`
  :`<b>Supervised gate not yet met.</b><br>needs ${g.min_cells} cells (${d.labelled}), ${g.min_per_class}+ in ${g.min_classes} classes (${d.usable_classes}), ${g.min_animals} animals (${d.animals})`;
 $('outpath').textContent='writing '+d.output_csv}
async function next(after){const r=await fetch('/api/next'+(after?`?after=${encodeURIComponent(after)}`:''));
 const d=await r.json();stats(d);
 if(d.done){$('patch').removeAttribute('src');$('meta').textContent='Every cell in this selection has been reviewed.';cur=null}
 else paint(d.cell)}
async function label(s){if(busy||!cur)return;const id=ident();if(!id)return;busy=true;
 try{const r=await fetch('/api/label',{method:'POST',headers:{'Content-Type':'application/json'},
  body:JSON.stringify({...id,cell_uid:cur.cell_uid,state:s})});const d=await r.json();
  if(!r.ok||!d.ok)throw Error(d.error||'could not save');stats(d);await next(cur.cell_uid)}
 catch(e){status(e.message)}finally{busy=false}}
async function unlabel(){if(!cur)return;
 const r=await fetch('/api/clear',{method:'POST',headers:{'Content-Type':'application/json'},
  body:JSON.stringify({cell_uid:cur.cell_uid})});const d=await r.json();stats(d);status('Label cleared.')}
$('skip').onclick=()=>cur&&next(cur.cell_uid);
$('unlabel').onclick=unlabel;
document.addEventListener('keydown',e=>{if(/input/i.test(e.target.tagName))return;
 const i=parseInt(e.key,10);if(i>=1&&i<=cfg.states.length){label(cfg.states[i-1]);e.preventDefault()}
 else if(e.key==='s'){cur&&next(cur.cell_uid);e.preventDefault()}});
(async()=>{const r=await fetch('/api/config');cfg=await r.json();renderStates();stats(cfg);await next('')})();
</script></body></html>'''


def loopback(host: str) -> bool:
    if host.lower() == "localhost":
        return True
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return False


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Browser UI for image-by-image microglia morphology review",
        epilog=("Reviews the cell patches an analysis run already extracted, with the "
                "model's proposal preselected, and writes a labels CSV that "
                "Fig5/01_analyze_gfap_iba1_microglia.py accepts through "
                "--microglia-labels-csv."))
    parser.add_argument("--analysis-dir", type=Path, required=True,
                        help=f"Analysis directory holding {TEMPLATE_NAME} and microglia_patches/")
    parser.add_argument("--output-csv", type=Path, default=None,
                        help=f"Accepted labels. Default: <analysis-dir>/{DEFAULT_OUTPUT_NAME}")
    parser.add_argument("--animals", nargs="*", default=None,
                        help="Restrict review to these animal_id values")
    parser.add_argument("--regions", nargs="*", default=None,
                        help="Restrict review to these regions, e.g. ARC ME VMN")
    parser.add_argument("--seed", type=int, default=20260825,
                        help="Seed for the interleaved review order")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--validate-only", action="store_true",
                        help="Report the label tally and the supervised gate, then exit")
    args = parser.parse_args(argv)
    analysis = args.analysis_dir.expanduser().resolve()
    args.analysis_dir = analysis
    args.output_csv = (args.output_csv.expanduser().resolve() if args.output_csv
                       else analysis / DEFAULT_OUTPUT_NAME)
    return args


def main(argv=None) -> int:
    args = parse_args(argv)
    if not 0 <= args.port <= 65535:
        raise SystemExit("--port must be 0..65535")
    review = Review(args.analysis_dir, args.output_csv, args.seed,
                    only_animals=args.animals, only_regions=args.regions)
    counts = review.counts()
    if args.validate_only:
        print(json.dumps(counts, indent=2, sort_keys=True))
        return 0 if counts["gate_met"] else 1
    if not loopback(args.host):
        print("WARNING: non-loopback binding exposes microscopy and label writes "
              "without authentication or TLS.", file=sys.stderr)
    server = Server((args.host, args.port), review)
    host, port = server.server_address[:2]
    browser = "127.0.0.1" if host in ("0.0.0.0", "::") else host
    print(f"Microglia morphology review: http://{browser}:{port}")
    print(f"Reviewable cells: {counts['total_cells']}   already labelled: {counts['labelled']}")
    print(f"Accepted labels: {review.output_csv}")
    print(f"Supervised gate: {'met' if counts['gate_met'] else 'not yet met'} "
          f"({GATE['min_cells']} cells, {GATE['min_per_class']}+ in "
          f"{GATE['min_classes']} classes, {GATE['min_animals']} animals)")
    print("Keys 1-4 label the shown cell, s skips it. Press Ctrl-C to stop.")
    try:
        server.serve_forever(0.25)
    except KeyboardInterrupt:
        print("\nStopping server.")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
