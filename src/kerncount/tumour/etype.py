"""The registered `tumour_volume` experiment type: chains layout -> measure -> readout table -> analyse."""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ..experiment_types import ExperimentType, RunOptions, register

if TYPE_CHECKING:
    import pandas as pd


@register
class TumourVolume(ExperimentType):
    """mCherry tumours in brain sections (ADR 0002-0004): one tumour = one n; main-mass volume from z-stacks."""
    name = "tumour_volume"
    label = "Tumour volume (z-stacks)"

    def settings_template(self) -> str | None:
        from .config import TEMPLATE

        return TEMPLATE

    def scan_layout(self, raw: Path, settings: dict) -> pd.DataFrame:
        from . import layout

        return layout.scan(raw, settings)

    def report_layout(self, layout: pd.DataFrame, settings: dict, raw: Path) -> None:
        import pandas as pd

        from .config import factors

        facs = [f.name for f in factors(settings)]
        with pd.option_context("display.width", 200, "display.max_colwidth", 70, "display.max_rows", 500):
            print(layout[["tumour_id", "imaging_session", "include", "note"]].to_string(index=False))
            inc = layout[layout["include"] == "yes"]
            print(f"\n{len(inc)} tumours included ({len(layout) - len(inc)} excluded); per group:")
            print(inc.groupby(facs).size().rename("tumours").to_string())
        legacy = int((layout["legacy_mask"] != "").sum())
        print(f"Fiji legacy masks found for {legacy} tumours (kept as legacy outputs, not tumours).")

    def measure(self, project_dir: Path, layout: pd.DataFrame, settings: dict, options: RunOptions) -> pd.DataFrame:
        from . import measure

        return measure.measure(project_dir, layout, settings, resegment=options.resegment)

    def readout_table(self, project_dir: Path, layout: pd.DataFrame, measurements: pd.DataFrame,
                      settings: dict) -> pd.DataFrame:
        from . import table

        t = table.build(layout, measurements, settings)
        table.write(project_dir, t, settings)
        return t

    def analyse(self, project_dir: Path, table: pd.DataFrame, settings: dict, options: RunOptions) -> None:
        from . import plots
        from .config import factors

        out = plots.group_plot(table, factors(settings), "main_mass_volume_mm3", "Main-mass volume (mm³)",
                               Path(project_dir) / "plots" / "main_mass_volume_absolute")
        print(f"Wrote {Path(project_dir) / 'tumours.csv'} and {out}")
