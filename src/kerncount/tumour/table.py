"""Layout + measurements -> the per-tumour readout table (tumours.csv): one row per included tumour = one n."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import batch, factors

FILENAME = "tumours.csv"
LEADING = ["tumour_id", "imaging_session", "session_label", "file_date", "culture_batch", "acquisition_day", "batch"]
NO_DAY = "(day not recorded)"
NO_BATCH = "(culture batch not recorded)"
# [tumour] batch -> (layout column, placeholder for an empty cell)
SOURCES = {"culture_batch": ("culture_batch", NO_BATCH), "session_folder": ("imaging_session", ""),
           "acquisition_day": ("acquisition_day", NO_DAY)}
IDENTITY = ["section", "hemisphere", "tumour", "file"]
READOUTS = ["cross_sectional_area_mm2", "perimeter_mm", "mid_section_centre_layer", "mid_section_layers",
            "mid_section_um", "mid_section_fits", "footprint_area_mm2", "slab_volume_mm3", "main_mass_volume_um3",
            "main_mass_volume_mm3", "main_mass_surface_mm2", "tumour_burden_mm3", "invasion_fraction", "n_satellites",
            "elongation", "solidity", "edge_roughness", "circularity"]
# Derived columns: (new column, measured column, divisor) -- mm / mm² / mm³ for readable plots and statistics.
UNITS = [("cross_sectional_area_mm2", "cross_sectional_area_um2", 1e6), ("perimeter_mm", "perimeter_um", 1e3),
         ("footprint_area_mm2", "footprint_area_um2", 1e6), ("slab_volume_mm3", "slab_volume_um3", 1e9),
         ("main_mass_volume_mm3", "main_mass_volume_um3", 1e9), ("main_mass_surface_mm2", "main_mass_surface_um2", 1e6),
         ("tumour_burden_mm3", "tumour_burden_um3", 1e9),
         ("fiji_cross_sectional_area_mm2", "fiji_cross_sectional_area_um2", 1e6),
         ("fiji_perimeter_mm", "fiji_perimeter_um", 1e3)]


def build(layout: pd.DataFrame, measurements: pd.DataFrame, settings: dict) -> pd.DataFrame:
    facs = [f.name for f in factors(settings)]
    inc = layout[layout["include"].astype(str).str.strip().str.lower().isin(("yes", "y", "true", "1"))]
    keep = [c for c in [*LEADING, *facs, *IDENTITY, "legacy_mask", "legacy_stats", "note"] if c in inc.columns]
    table = inc[keep].merge(measurements, on="tumour_id", how="left", validate="one_to_one")
    # tickets 12, 13: the batch of every tumour ([tumour] batch): block effect + confounding rule (stats.batch_of)
    src, empty = SOURCES[batch(settings)]
    if src in table.columns:
        b = table[src].fillna("").astype(str)
        table["batch"] = b.where(b != "", empty) if empty else b
    for new_col, col, div in UNITS:
        if col in table.columns:
            table[new_col] = pd.to_numeric(table[col], errors="coerce") / div
    first = [c for c in [*LEADING, *facs, *IDENTITY, *READOUTS] if c in table.columns]
    last = [c for c in ("legacy_mask", "legacy_stats", "note") if c in table.columns]
    rest = [c for c in table.columns if c not in first and c not in last]
    return table[first + rest + last].reset_index(drop=True)


def write(project_dir: Path, table: pd.DataFrame, settings: dict) -> Path:
    from ..fileio import atomic_write

    path = Path(project_dir) / FILENAME
    atomic_write(path, table.to_csv(index=False, sep=settings["project"]["csv_delimiter"]))
    return path
