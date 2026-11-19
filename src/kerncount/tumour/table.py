"""Layout + measurements -> the per-tumour readout table (tumours.csv): one row per included tumour = one n."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import factors

FILENAME = "tumours.csv"
LEADING = ["tumour_id", "imaging_session", "date"]
IDENTITY = ["section", "hemisphere", "tumour", "file"]
READOUTS = ["main_mass_volume_um3", "main_mass_volume_mm3"]


def build(layout: pd.DataFrame, measurements: pd.DataFrame, settings: dict) -> pd.DataFrame:
    facs = [f.name for f in factors(settings)]
    inc = layout[layout["include"].astype(str).str.strip().str.lower().isin(("yes", "y", "true", "1"))]
    keep = [c for c in [*LEADING, *facs, *IDENTITY, "legacy_mask", "legacy_stats", "note"] if c in inc.columns]
    table = inc[keep].merge(measurements, on="tumour_id", how="left", validate="one_to_one")
    table["main_mass_volume_mm3"] = table["main_mass_volume_um3"] / 1e9
    first = [c for c in [*LEADING, *facs, *IDENTITY, *READOUTS] if c in table.columns]
    last = [c for c in ("legacy_mask", "legacy_stats", "note") if c in table.columns]
    rest = [c for c in table.columns if c not in first and c not in last]
    return table[first + rest + last].reset_index(drop=True)


def write(project_dir: Path, table: pd.DataFrame, settings: dict) -> Path:
    from ..fileio import atomic_write

    path = Path(project_dir) / FILENAME
    atomic_write(path, table.to_csv(index=False, sep=settings["project"]["csv_delimiter"]))
    return path
