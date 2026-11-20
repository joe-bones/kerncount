"""Threshold sweep (ADR 0002): the readout table recomputed with the primary threshold x each `[tumour] sweep_scales`.

Each original stack is read and smoothed once and thresholded at every scale (`segment.sweep`); only the per-tumour
readouts are kept (no labels, nothing full-resolution on disk), cached in `_work/tumour/<slug>.sweep.json` under the
same key as the primary measurement plus the scales. The per-scale tables go through `table.build`, so every column
derived there (e.g. mm³) is derived identically. A readout is sweepable when `segment.segment` (and `table.build`)
produce it; confirmatory readouts computed elsewhere are reported as "not available in the sweep".

Review corrections are respected: `table` holds only the tumours left after review exclusions, and for a tumour whose
main mass was re-picked (`main_mass_object` != 1) the main mass at every scale is the object overlapping the corrected
x1 main mass (from the cached labels) most; the pick is part of the cache key.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd

from . import table as table_mod
from .config import tumour_settings
from .measure import WORK, _key, labels_path, raw_dir, slug


def scales(settings: dict) -> list[float]:
    """Declared sweep scales (without 1 = the primary result), sorted; validated as positive."""
    out = sorted({float(s) for s in tumour_settings(settings)["sweep_scales"]} - {1.0})
    if any(s <= 0 for s in out):
        raise ValueError(f"[tumour] sweep_scales must be positive, got {out}")
    return out


def _pick(row: pd.Series) -> int:
    """The review's main-mass object of a table row (1 = the automatic largest object)."""
    v = pd.to_numeric(pd.Series([row.get("main_mass_object", 1)]), errors="coerce").iloc[0]
    return int(v) if pd.notna(v) else 1


def measurements(project_dir: Path, table: pd.DataFrame, settings: dict) -> dict[float, pd.DataFrame]:
    """{scale: tumour_id + segmentation readouts at that scale} for every tumour of `table` (cached)."""
    from .segment import sweep_file

    prj = Path(project_dir)
    work = prj / WORK
    work.mkdir(parents=True, exist_ok=True)
    raw = raw_dir(prj, settings)
    ts = tumour_settings(settings)
    sc = scales(settings)
    rows: dict[float, list[dict]] = {s: [] for s in sc}
    for i, (_, row) in enumerate(table.iterrows(), 1):
        tid, src = row["tumour_id"], raw / row["file"]
        pick = _pick(row)
        key = f"{_key(src, ts)}:{','.join(f'{s:g}' for s in sc)}" + (f":main={pick}" if pick != 1 else "")
        path = work / f"{slug(tid)}.sweep.json"
        cached = json.loads(path.read_text()) if path.exists() else None
        if cached and cached.get("key") == key:
            res = {float(k): v for k, v in cached["scales"].items()}
        else:
            t0 = time.time()
            ref = None
            if pick != 1:  # re-picked at review: follow the corrected main mass
                import tifffile

                ref = tifffile.imread(labels_path(prj, tid)) == pick
            res = sweep_file(src, ts, sc, ref)
            del ref
            path.write_text(json.dumps({"tumour_id": tid, "key": key, "scales": {f"{k:g}": v for k, v in res.items()}},
                                       indent=1, default=float))
            print(f"  sweep [{i}/{len(table)}] {tid}: {len(sc)} thresholds ({time.time() - t0:.0f} s)", flush=True)
        for s in sc:
            rows[s].append({"tumour_id": tid, **res[s]})
    return {s: pd.DataFrame(r) for s, r in rows.items()}


def rebuild(table: pd.DataFrame, sweep_meas: pd.DataFrame, settings: dict) -> pd.DataFrame:
    """`table` with the sweep's segmentation readouts in place of the primary ones, derived columns re-derived."""
    from .config import factors

    facs = [f.name for f in factors(settings)]
    ident = [c for c in [*table_mod.LEADING, *facs, *table_mod.IDENTITY, "legacy_mask", "legacy_stats", "note"]
             if c in table.columns]
    layout = table[ident].assign(include="yes")
    meas = table.drop(columns=[c for c in ident if c != "tumour_id"]).set_index("tumour_id")
    new = sweep_meas.set_index("tumour_id")
    for c in new.columns:
        meas[c] = new[c].reindex(meas.index)
    return table_mod.build(layout, meas.reset_index(), settings)


def tables(project_dir: Path, table: pd.DataFrame, settings: dict) -> dict[float, pd.DataFrame]:
    """{scale: readout table at that threshold scale}, scale 1 = `table` itself."""
    out = {s: rebuild(table, m, settings) for s, m in measurements(project_dir, table, settings).items()}
    out[1.0] = table
    return dict(sorted(out.items()))
