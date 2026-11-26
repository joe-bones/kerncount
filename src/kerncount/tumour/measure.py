"""Per-tumour measurement with a cache under <project>/_work/tumour/ (re-runs only segment what changed).

Per tumour: `<slug>.json` (readouts + cache key) and `<slug>.labels.tif` (zlib-compressed uint8/uint16 label image:
1 = main mass, 2.. = satellites). The key has two parts: the **labels key** (segmentation settings + the original
file's size / mtime) and the **readouts key** (mid-section / outline scales + readout definitions). When only the
readouts key changed, the readouts are re-derived from the cached labels (seconds) instead of re-segmenting.
Tickets 05 (shape, burden), 06 (truncation, re-picking the main mass) and 07 (Fiji) read the cached labels instead of
re-segmenting.
"""
from __future__ import annotations

import hashlib
import json
import re
import time
from pathlib import Path

import pandas as pd

from .config import tumour_settings
from .segment import measure_file
from . import shape
from .shape import KEYS as SHAPE_KEYS, ShapeOptions

WORK = Path("_work") / "tumour"
SEG_KEYS = ("smooth_sigma_um", "threshold", "min_object_um3")  # change the label image -> re-segment
# Bump when the label image changes (segmentation): every stack is re-segmented.
LABELS_VERSION = 1
# Bump when the readouts derived from the labels change (new columns / definitions): re-derived from cached labels.
# 2 = ticket 05 (mid-section, area, perimeter, shape, burden); 3 = exact-thickness mid-section (review fix 2026-11-20).
READOUTS_VERSION = 3
CACHE_VERSION = 2  # the single version of caches written before the review fix (see `_legacy_key`)


def slug(tumour_id: str) -> str:
    """Unique file-system-safe name for a tumour ID: readable part + short hash of the exact ID
    ("TLR7 KO / HERV / -clodro / n3A" -> "TLR7-KO_HERV_-clodro_n3A_<6 hex>")."""
    parts = [re.sub(r"[^A-Za-z0-9+\-.]+", "-", p.strip()) for p in str(tumour_id).split("/")]
    readable = "_".join(p for p in parts if p) or "tumour"
    return f"{readable}_{hashlib.sha1(str(tumour_id).encode()).hexdigest()[:6]}"


def raw_dir(project_dir: Path, settings: dict) -> Path:
    raw = settings["project"]["raw_dir"]
    if not raw:
        raise SystemExit("settings.toml: project.raw_dir is empty")
    p = Path(raw)
    return p if p.is_absolute() else (Path(project_dir) / p).resolve()


def _hash(blob: dict) -> str:
    return hashlib.sha1(json.dumps(blob, sort_keys=True, default=str).encode()).hexdigest()[:16]


def _file(path: Path) -> list:
    st = path.stat()
    return [path.name, st.st_size, st.st_mtime_ns]


def labels_key(path: Path, ts: dict) -> str:
    """Identity of the cached label image: segmentation settings + the original file."""
    return _hash({"seg": {k: ts[k] for k in SEG_KEYS}, "file": _file(path), "version": LABELS_VERSION})


def readouts_key(ts: dict) -> str:
    """Identity of the readouts derived from a label image: the [tumour] readout scales + definitions."""
    return _hash({"shape": {k: ts[k] for k in SHAPE_KEYS}, "version": READOUTS_VERSION})


def _legacy_key(path: Path, ts: dict) -> str:
    """The single key of caches written before the review fix (segmentation + shape settings, CACHE_VERSION): its
    labels are still valid when it matches the current settings."""
    return _hash({"seg": {k: ts[k] for k in (*SEG_KEYS, *SHAPE_KEYS)}, "file": _file(path), "version": CACHE_VERSION})


def _key(path: Path, ts: dict) -> str:
    """The key of one tumour's measurement (labels + readouts): QC, the sweep and the shared-threshold check build
    on it, so their caches follow any change of either part."""
    return f"{labels_key(path, ts)}:{readouts_key(ts)}"


def labels_path(project_dir: Path, tumour_id: str) -> Path:
    return Path(project_dir) / WORK / f"{slug(tumour_id)}.labels.tif"


def measure(project_dir: Path, layout: pd.DataFrame, settings: dict, resegment: bool = False) -> pd.DataFrame:
    """One row per included tumour: tumour_id + segmentation readouts (segment.segment)."""
    import tifffile

    from .segment import MissingCalibration

    prj = Path(project_dir)
    work = prj / WORK
    work.mkdir(parents=True, exist_ok=True)
    raw = raw_dir(prj, settings)
    ts = tumour_settings(settings)
    included = layout[layout["include"].astype(str).str.strip().str.lower().isin(("yes", "y", "true", "1"))]
    rows = []
    for i, (_, row) in enumerate(included.iterrows(), 1):
        tid = row["tumour_id"]
        src = raw / row["file"]
        if not src.is_file():
            raise SystemExit(f"{tid}: raw file not found: {src}")
        key, lkey = _key(src, ts), labels_key(src, ts)
        meta, lab = work / f"{slug(tid)}.json", labels_path(prj, tid)
        cached = json.loads(meta.read_text()) if meta.exists() and lab.exists() and not resegment else None
        if cached and cached.get("key") == key:
            rows.append({"tumour_id": tid, **cached["readouts"]})
            continue
        if cached and (cached.get("labels_key") == lkey or cached.get("key") == _legacy_key(src, ts)):
            ro = cached["readouts"]  # labels still valid: only the readouts derived from them changed
            vox = (float(ro["voxel_z_um"]), float(ro["voxel_y_um"]), float(ro["voxel_x_um"]))
            ro = {**ro, **shape.label_readouts(tifffile.imread(lab), vox, ShapeOptions.from_settings(ts))}
            meta.write_text(json.dumps({"tumour_id": tid, "file": row["file"], "key": key, "labels_key": lkey,
                                        "readouts": ro}, indent=1, default=float))
            rows.append({"tumour_id": tid, **ro})
            print(f"  [{i}/{len(included)}] {tid}: readouts re-derived from the cached labels", flush=True)
            continue
        t0 = time.time()
        try:
            seg = measure_file(src, ts)
        except MissingCalibration as e:
            raise SystemExit(f"{tid}: {e}") from None
        tifffile.imwrite(lab, seg.labels, compression="zlib", compressionargs={"level": 6})
        meta.write_text(json.dumps({"tumour_id": tid, "file": row["file"], "key": key, "labels_key": lkey,
                                    "readouts": seg.readouts}, indent=1, default=float))
        rows.append({"tumour_id": tid, **seg.readouts})
        print(f"  [{i}/{len(included)}] {tid}: main mass {seg.readouts['main_mass_volume_um3'] / 1e9:.4g} mm³ "
              f"(threshold {seg.readouts['threshold']:.1f}, {time.time() - t0:.0f} s)", flush=True)
    return pd.DataFrame(rows)
