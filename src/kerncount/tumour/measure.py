"""Per-tumour measurement with a cache under <project>/_work/tumour/ (re-runs only segment what changed).

Per tumour: `<slug>.json` (readouts + cache key) and `<slug>.labels.tif` (zlib-compressed uint8/uint16 label image:
1 = main mass, 2.. = satellites). The key covers the segmentation settings and the original file's size / mtime.
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
from .shape import KEYS as SHAPE_KEYS

WORK = Path("_work") / "tumour"
SEG_KEYS = ("smooth_sigma_um", "threshold", "min_object_um3", *SHAPE_KEYS)
# Bump when the readouts computed per tumour change (new columns / definitions): old caches are then re-measured
# instead of silently lacking columns. 2 = ticket 05 (mid-section, area, perimeter, shape, burden).
CACHE_VERSION = 2


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


def _key(path: Path, ts: dict) -> str:
    st = path.stat()
    blob = json.dumps({"seg": {k: ts[k] for k in SEG_KEYS}, "file": [path.name, st.st_size, st.st_mtime_ns],
                       "version": CACHE_VERSION}, sort_keys=True, default=str)
    return hashlib.sha1(blob.encode()).hexdigest()[:16]


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
        key = _key(src, ts)
        meta, lab = work / f"{slug(tid)}.json", labels_path(prj, tid)
        cached = json.loads(meta.read_text()) if meta.exists() and lab.exists() and not resegment else None
        if cached and cached.get("key") == key:
            rows.append({"tumour_id": tid, **cached["readouts"]})
            continue
        t0 = time.time()
        try:
            seg = measure_file(src, ts)
        except MissingCalibration as e:
            raise SystemExit(f"{tid}: {e}") from None
        tifffile.imwrite(lab, seg.labels, compression="zlib", compressionargs={"level": 6})
        meta.write_text(json.dumps({"tumour_id": tid, "file": row["file"], "key": key, "readouts": seg.readouts},
                                   indent=1, default=float))
        rows.append({"tumour_id": tid, **seg.readouts})
        print(f"  [{i}/{len(included)}] {tid}: main mass {seg.readouts['main_mass_volume_um3'] / 1e9:.4g} mm³ "
              f"(threshold {seg.readouts['threshold']:.1f}, {time.time() - t0:.0f} s)", flush=True)
    return pd.DataFrame(rows)
