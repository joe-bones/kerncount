"""Fiji legacy readout per tumour (ADR 0002): exploratory `fiji_*` columns of tumours.csv, always labelled legacy.

Source, in order of preference (recorded per tumour in `fiji_source`):
1. "Fiji statistics" — the lab's own 3D Objects Counter table (`legacy_stats` in layout.csv);
2. "Fiji mask" — objects counted on the lab's Fiji mask (`legacy_mask`) with `fiji.objects_3d`, which reproduces
   3D Objects Counter's numbers exactly;
3. "re-implementation" — the macro re-run on the original stack (`fiji.macro_mask`, voxel-identical to the lab's
   masks on the stacks checked), when neither exists and `[tumour] fiji_fallback = true`.
Unparseable rows of a statistics file are skipped and listed in `fiji_note` (and printed), never guessed.
Per-tumour results are cached in `_work/tumour/<slug>.fiji.json` (key: the files' size / mtime + settings), Fiji's main
object in `<slug>.fiji.labels.tif` (zlib uint8), so a change of the mid-section rule or scales re-derives the
`fiji_*` section readouts from it instead of re-running the macro.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from typing import TYPE_CHECKING

import pandas as pd

from . import fiji
from .config import tumour_settings
from .measure import WORK, raw_dir, readouts_key, slug
from .shape import ShapeOptions

if TYPE_CHECKING:
    import numpy as np

SOURCES = ("Fiji statistics", "Fiji mask", "re-implementation")


def _stat(p: Path) -> list:
    return [p.name, p.stat().st_size, p.stat().st_mtime_ns] if p.is_file() else [p.name, None]


def _key(files: list[Path], ts: dict) -> str:
    """Identity of Fiji's objects (volumes + main-object mask): the files and the object settings."""
    blob = json.dumps({"files": [_stat(p) for p in files], "min": ts["fiji_min_voxels"],
                       "fallback": bool(ts["fiji_fallback"]), "version": 3}, sort_keys=True, default=str)
    return hashlib.sha1(blob.encode()).hexdigest()[:16]


# Fiji counterparts of the confirmatory mid-section readouts (same rule as the primary method, ticket 05)
SECTION_COLUMNS = ("mid_section_centre_layer", "mid_section_layers", "mid_section_um", "mid_section_fits",
                   "cross_sectional_area_um2", "perimeter_um")


def _sections(labels, vox: tuple[float, float, float], opts: ShapeOptions | None, source: str) -> dict:
    """fiji_<mid-section readouts> of Fiji's largest object (label 1 of `fiji.objects_3d`)."""
    import numpy as np

    from . import shape

    main = (labels == 1).astype(np.uint8)
    r = shape.readouts(main, vox, opts)
    return {**{f"fiji_{k}": r[k] for k in SECTION_COLUMNS}, "fiji_section_source": source if main.any() else ""}


def _no_sections() -> dict:
    return {**{f"fiji_{k}": (math.nan if k not in ("mid_section_fits",) else "") for k in SECTION_COLUMNS},
            "fiji_section_source": ""}


def tumour_readouts(original: Path, stats: Path | None, mask: Path | None, min_voxels: int = fiji.MIN_VOXELS,
                    fallback: bool = True, n_layers: int | None = None, shape_opts: ShapeOptions | None = None
                    ) -> dict:
    """Legacy readouts of one tumour from the best available source (see module doc). The mid-section area and
    perimeter of Fiji's main object (`fiji_cross_sectional_area_um2`, `fiji_perimeter_um`) need a mask: the lab's,
    else the re-implementation's (`fiji_section_source`)."""
    return _readouts(original, stats, mask, min_voxels, fallback, n_layers, shape_opts)[0]


def _readouts(original: Path, stats: Path | None, mask: Path | None, min_voxels: int = fiji.MIN_VOXELS,
              fallback: bool = True, n_layers: int | None = None, shape_opts: ShapeOptions | None = None
              ) -> tuple[dict, np.ndarray | None]:
    """(`tumour_readouts`, Fiji's main object as a uint8 mask or None)."""
    import numpy as np
    from .segment import read_stack, voxel_size_um

    notes: list[str] = []
    vol: dict | None = None
    if stats is not None and stats.is_file():
        st = fiji.read_stats(stats)
        notes += st.problems
        if st.objects:
            vol = {**fiji.readouts(st.objects, min_voxels, n_layers), "fiji_source": SOURCES[0]}
    if mask is not None and mask.is_file():
        try:
            m = fiji.read_mask(mask)
            vox = voxel_size_um(original)
            labels, objs = fiji.objects_3d(m, vox, min_voxels)
            if vol is None:
                vol = {**fiji.readouts(objs, min_voxels, n_layers or m.shape[0]), "fiji_source": SOURCES[1]}
            return ({**vol, **_sections(labels, vox, shape_opts, SOURCES[1]), "fiji_note": "; ".join(notes)},
                    (labels == 1).astype(np.uint8))
        except Exception as e:  # unreadable / not a stack: fall through to the re-implementation
            notes.append(f"Fiji mask {mask.name} unreadable ({type(e).__name__})")
    if not fallback:
        notes.append("no Fiji mask; re-implementation switched off ([tumour] fiji_fallback)" if vol else
                     "no Fiji output; re-implementation switched off ([tumour] fiji_fallback)")
        return ({**(vol or {**fiji.readouts([], min_voxels), "fiji_source": ""}), **_no_sections(),
                 "fiji_note": "; ".join(notes)}, None)
    stack = read_stack(original)
    if stack.dtype != "uint8":
        notes.append(f"re-implementation needs an 8-bit stack (got {stack.dtype})")
        return ({**(vol or {**fiji.readouts([], min_voxels), "fiji_source": ""}), **_no_sections(),
                 "fiji_note": "; ".join(notes)}, None)
    m, _ = fiji.macro_mask(stack)
    del stack
    vox = voxel_size_um(original)
    labels, objs = fiji.objects_3d(m, vox, min_voxels)
    if vol is None:
        vol = {**fiji.readouts(objs, min_voxels, n_layers or m.shape[0]), "fiji_source": SOURCES[2]}
    return ({**vol, **_sections(labels, vox, shape_opts, SOURCES[2]), "fiji_note": "; ".join(notes)},
            (labels == 1).astype(np.uint8))


def measure(project_dir: Path, layout: pd.DataFrame, settings: dict, measurements: pd.DataFrame | None = None
            ) -> pd.DataFrame:
    """tumour_id + fiji_* legacy readouts for every included tumour (empty frame when `fiji_legacy = false`)."""
    import tifffile

    from .segment import voxel_size_um

    ts = tumour_settings(settings)
    if not ts["fiji_legacy"]:
        return pd.DataFrame({"tumour_id": []})
    prj = Path(project_dir)
    work = prj / WORK
    work.mkdir(parents=True, exist_ok=True)
    raw = raw_dir(prj, settings)
    layers = {}
    if measurements is not None and "n_layers" in measurements.columns:
        layers = dict(zip(measurements["tumour_id"], measurements["n_layers"]))
    included = layout[layout["include"].astype(str).str.strip().str.lower().isin(("yes", "y", "true", "1"))]
    rows = []
    for _, row in included.iterrows():
        tid = row["tumour_id"]
        src = raw / row["file"]
        stats = raw / row["legacy_stats"] if str(row.get("legacy_stats", "") or "").strip() else None
        mask = raw / row["legacy_mask"] if str(row.get("legacy_mask", "") or "").strip() else None
        nl = layers.get(tid)
        nl = int(nl) if nl is not None and not (isinstance(nl, float) and math.isnan(nl)) else None
        okey = _key([src] + [p for p in (stats, mask) if p is not None], ts)
        key = f"{okey}:{readouts_key(ts)}"
        cache, mpath = work / f"{slug(tid)}.fiji.json", work / f"{slug(tid)}.fiji.labels.tif"
        cached = json.loads(cache.read_text()) if cache.exists() else None
        if cached and cached.get("key") == key and cached.get("n_layers") == nl:
            r = cached["readouts"]
        elif (cached and cached.get("objects_key") == okey and cached.get("n_layers") == nl
              and (mpath.exists() or not cached["readouts"].get("fiji_section_source"))):
            r = dict(cached["readouts"])  # Fiji's objects unchanged: only the section readouts are re-derived
            if r.get("fiji_section_source"):
                r.update(_sections(tifffile.imread(mpath), voxel_size_um(src), ShapeOptions.from_settings(ts),
                                   r["fiji_section_source"]))
            cache.write_text(json.dumps({"tumour_id": tid, "key": key, "objects_key": okey, "n_layers": nl,
                                         "readouts": r}, indent=1, default=float))
        else:
            r, main = _readouts(src, stats, mask, int(ts["fiji_min_voxels"]), bool(ts["fiji_fallback"]), nl,
                                ShapeOptions.from_settings(ts))
            if main is not None:
                tifffile.imwrite(mpath, main, compression="zlib", compressionargs={"level": 6})
            cache.write_text(json.dumps({"tumour_id": tid, "key": key, "objects_key": okey, "n_layers": nl,
                                         "readouts": r}, indent=1, default=float))
            if r["fiji_note"]:
                print(f"  {tid}: Fiji legacy: {r['fiji_note']}", flush=True)
        rows.append({"tumour_id": tid, **r})
    return pd.DataFrame(rows)
