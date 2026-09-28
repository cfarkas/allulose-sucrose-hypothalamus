"""Post-segmentation POMC size QC against DAPI in the same section/region.

The factor is fixed at 1.0 before inspecting treatment contrasts. This is a
size heuristic, not a biological identification of neuronal somata. Source
images, historical training annotations and HIL anatomical masks stay intact.
"""
from __future__ import annotations

import numpy as np
from skimage.measure import regionprops

SCHEMA = "pomc-local-dapi-size-qc-v1"
AREA_FACTOR = 1.0
REGIONS = {0: "OUTSIDE_DELINEATED_REGIONS", 1: "ARC", 2: "ME", 3: "VMN"}


def _code(regions, centroid):
    y, x = (int(round(float(v))) for v in centroid)
    y = min(max(y, 0), regions.shape[0] - 1)
    x = min(max(x, 0), regions.shape[1] - 1)
    return int(regions[y, x])


def filter_pomc_by_local_dapi(pomc, dapi, regions, *, area_factor=AREA_FACTOR):
    """Remove whole POMC instances smaller than regional mean DAPI area.

    Both marker and DAPI instances are assigned by centroid. Areas use full
    masks, not region-clipped fragments. Equality is retained. Empty regional
    DAPI references cause no exclusion and are recorded for review. No
    treatment, c-FOS status, p value or POMC frequency enters the decision.
    Labels of retained objects are unchanged.
    """
    if not np.isfinite(area_factor) or area_factor <= 0:
        raise ValueError("Area factor must be finite and positive")
    pomc, dapi, regions = map(np.asarray, (pomc, dapi, regions))
    if pomc.ndim != 2 or not (pomc.shape == dapi.shape == regions.shape):
        raise ValueError("POMC, DAPI and accepted regions must share a 2-D native geometry")
    if not np.issubdtype(pomc.dtype, np.integer) or not np.issubdtype(dapi.dtype, np.integer):
        raise ValueError("POMC and DAPI must contain integer instance labels")
    if min(int(pomc.min()), int(dapi.min())) < 0:
        raise ValueError("Instance labels must be nonnegative")
    references = {}
    for obj in regionprops(dapi):
        references.setdefault(_code(regions, obj.centroid), []).append(int(obj.area))
    rows = []
    removed = []
    for obj in regionprops(pomc):
        code = _code(regions, obj.centroid)
        values = references.get(code, [])
        mean = float(np.mean(values)) if values else None
        threshold = mean * area_factor if mean is not None else None
        reject = threshold is not None and obj.area < threshold
        if reject:
            removed.append(int(obj.label))
        rows.append({
            "pomc_label": int(obj.label), "region_code": code,
            "region": REGIONS.get(code, str(code)), "area_px": int(obj.area),
            "centroid_y": float(obj.centroid[0]), "centroid_x": float(obj.centroid[1]),
            "dapi_reference_n": len(values), "mean_dapi_area_px": mean,
            "minimum_pomc_area_px": threshold, "excluded_small": bool(reject),
            "decision": "exclude_below_local_dapi_mean" if reject else
                        "retain" if values else "retain_no_local_dapi_reference",
        })
    cleaned = pomc.copy()
    cleaned[np.isin(cleaned, removed)] = 0
    return cleaned, {
        "schema": SCHEMA, "area_factor": float(area_factor),
        "source_masks_modified": False, "regions_assigned_by": "object centroid",
        "area_definition": "full instance area in native pixels",
        "objects_before": len(rows), "objects_removed": len(removed),
        "objects_after": len(rows) - len(removed),
        "dapi_references": {
            str(code): {"n": len(values), "mean_area_px": float(np.mean(values))}
            for code, values in sorted(references.items())
        }, "objects": rows,
    }
