"""Fiji legacy readout per tumour (ADR 0002): exploratory `fiji_*` columns of tumours.csv, always labelled legacy.

Source, in order of preference (recorded per tumour in `fiji_source`):
1. "Fiji statistics" — the lab's own 3D Objects Counter table (`legacy_stats` in layout.csv);
2. "Fiji mask" — objects counted on the lab's Fiji mask (`legacy_mask`) with `fiji.objects_3d`, which reproduces
   3D Objects Counter's numbers exactly;
3. "re-implementation" — the macro re-run on the original stack (`fiji.macro_mask`, voxel-identical to the lab's
   masks on the stacks checked), when neither exists and `[tumour] fiji_fallback = true`.
Unparseable rows of a statistics file are skipped and listed in `fiji_note` (and printed), never guessed.
Per-tumour results are cached in `_work/tumour/<slug>.fiji.json` (key: the files' size / mtime + settings).
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import pandas as pd

from . import fiji
from .config import tumour_settings
from .measure import WORK, raw_dir, slug

SOURCES = ("Fiji statistics", "Fiji mask", "re-implementation")


def _stat(p: Path) -> list:
    return [p.name, p.stat().st_size, p.stat().st_mtime_ns] if p.is_file() else [p.name, None]


def _key(files: list[Path], ts: dict) -> str:
    blob = json.dumps({"files": [_stat(p) for p in files], "min": ts["fiji_min_voxels"],
                       "fallback": bool(ts["fiji_fallback"]), "version": 1}, sort_keys=True, default=str)
    return hashlib.sha1(blob.encode()).hexdigest()[:16]


def tumour_readouts(original: Path, stats: Path | None, mask: Path | None, min_voxels: int = fiji.MIN_VOXELS,
                    fallback: bool = True, n_layers: int | None = None) -> dict:
    """Legacy readouts of one tumour from the best available source (see module doc)."""
    from .segment import read_stack, voxel_size_um

    notes: list[str] = []
    if stats is not None and stats.is_file():
        st = fiji.read_stats(stats)
        notes += st.problems
        if st.objects:
            return {**fiji.readouts(st.objects, min_voxels, n_layers), "fiji_source": SOURCES[0],
                    "fiji_note": "; ".join(notes)}
    if mask is not None and mask.is_file():
        try:
            m = fiji.read_mask(mask)
            _, objs = fiji.objects_3d(m, voxel_size_um(original), min_voxels)
            return {**fiji.readouts(objs, min_voxels, n_layers or m.shape[0]), "fiji_source": SOURCES[1],
                    "fiji_note": "; ".join(notes)}
        except Exception as e:  # unreadable / not a stack: fall through to the re-implementation
            notes.append(f"Fiji mask {mask.name} unreadable ({type(e).__name__})")
    if not fallback:
        notes.append("no Fiji output; re-implementation switched off ([tumour] fiji_fallback)")
        return {**fiji.readouts([], min_voxels), "fiji_source": "", "fiji_note": "; ".join(notes)}
    stack = read_stack(original)
    if stack.dtype != "uint8":
        notes.append(f"re-implementation needs an 8-bit stack (got {stack.dtype})")
        return {**fiji.readouts([], min_voxels), "fiji_source": "", "fiji_note": "; ".join(notes)}
    m, _ = fiji.macro_mask(stack)
    del stack
    _, objs = fiji.objects_3d(m, voxel_size_um(original), min_voxels)
    return {**fiji.readouts(objs, min_voxels, n_layers or m.shape[0]), "fiji_source": SOURCES[2],
            "fiji_note": "; ".join(notes)}


def measure(project_dir: Path, layout: pd.DataFrame, settings: dict, measurements: pd.DataFrame | None = None
            ) -> pd.DataFrame:
    """tumour_id + fiji_* legacy readouts for every included tumour (empty frame when `fiji_legacy = false`)."""
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
        key = _key([src] + [p for p in (stats, mask) if p is not None], ts)
        cache = work / f"{slug(tid)}.fiji.json"
        cached = json.loads(cache.read_text()) if cache.exists() else None
        if cached and cached.get("key") == key and cached.get("n_layers") == nl:
            r = cached["readouts"]
        else:
            r = tumour_readouts(src, stats, mask, int(ts["fiji_min_voxels"]), bool(ts["fiji_fallback"]), nl)
            cache.write_text(json.dumps({"tumour_id": tid, "key": key, "n_layers": nl, "readouts": r}, indent=1,
                                        default=float))
            if r["fiji_note"]:
                print(f"  {tid}: Fiji legacy: {r['fiji_note']}", flush=True)
        rows.append({"tumour_id": tid, **r})
    return pd.DataFrame(rows)
