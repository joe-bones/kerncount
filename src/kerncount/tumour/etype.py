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

    def merge_layout(self, old: pd.DataFrame, fresh: pd.DataFrame, settings: dict
                     ) -> tuple[pd.DataFrame, list[str]] | None:
        """Ticket 12: files moved (session folder renamed / reorganised), new files or columns added -> the updated
        layout with every edit carried over; None when nothing changed (layout.csv kept, layout.new.csv convention)."""
        from . import layout

        m = layout.migrate(old, fresh, settings, add_new=True)
        return (m.layout, m.describe()) if m.changed else None

    layout_hint = ("Check the tumour IDs, factor columns and notes in layout.csv: fix a value by editing its cell "
                   "(noted as an override), or set include = no and give the reason in note")

    def report_layout(self, layout: pd.DataFrame, settings: dict, raw: Path) -> None:
        import pandas as pd

        from . import layout as layout_mod
        from .config import factors

        facs = [f.name for f in factors(settings)]
        v = layout_mod.validate(layout, raw, settings)
        with pd.option_context("display.width", 200, "display.max_colwidth", 90, "display.max_rows", 500):
            print(layout[["tumour_id", "imaging_session", "include", "note"]].to_string(index=False))
            if "culture_batch" in layout.columns and len(layout):  # ticket 13: culture batches per session folder
                cb = layout.assign(batch=layout["culture_batch"].replace("", "(not recorded)"))
                print("\nTumours per imaging session and culture batch (file-name date = culture preparation):")
                print(pd.crosstab(cb["imaging_session"], cb["batch"]).to_string())
            if "acquisition_day" in layout.columns and len(layout):  # ticket 12: imaging days per session folder
                days = layout.assign(day=layout["acquisition_day"].replace("", "(not recorded)"))
                print("\nTumours per imaging session and acquisition day (OME metadata):")
                print(pd.crosstab(days["imaging_session"], days["day"]).to_string())
            inc = v.layout[(v.layout["include"] == "yes") & (v.layout["tumour_id"] != "")]
            print(f"\n{len(inc)} tumours included ({len(layout) - len(inc)} excluded or unreadable); per group:")
            if len(inc):
                print(inc.groupby(facs).size().rename("tumours").to_string())
        legacy = int((layout["legacy_mask"] != "").sum())
        stats = int((layout["legacy_stats"] != "").sum())
        print(f"Fiji legacy outputs (not tumours): masks for {legacy}, statistics for {stats} tumours.")
        if "czi" in layout.columns:
            print(f"Zeiss .czi twins recorded as provenance: {int((layout['czi'] != '').sum())}.")
        other = layout_mod.unmatched_files(raw, layout)
        if other:
            print(f"Other files, neither tumour nor matched legacy output (ignored): {'; '.join(other)}")
        for w in v.warnings:
            print(f"Note: {w}")
        if v.problems:
            print("\n" + str(layout_mod.LayoutBlocked(v.problems).code))

    def run(self, project_dir: Path, options: RunOptions) -> None:
        """Validate layout.csv against the raw folder first (blocks with every problem listed), then the base chain
        on the validated layout; excluded tumours and their reasons go to excluded_tumours.csv."""
        from .. import discover, settings as settings_mod
        from ..fileio import atomic_write
        from . import layout as layout_mod
        from .measure import raw_dir

        s = settings_mod.load(project_dir)
        from . import lock_log
        from .results import lock

        try:  # review fix: the CLI enforces the analysis-plan lock like the web interface (ADR 0003)
            note = lock_log.enforce(Path(project_dir), lock(s)[0], options.unlock)
        except ValueError as e:  # PlanLocked or a missing reason
            raise SystemExit(str(e)) from None
        if note:
            print(note)
        lp = Path(project_dir) / "layout.csv"
        if not lp.exists():
            raise SystemExit(f"{lp} not found - run 'kerncount init' first")
        v = layout_mod.check(discover.read_table(lp), raw_dir(Path(project_dir), s), s)
        for w in v.warnings:
            print(f"Note: {w}")
        import pandas as pd

        from . import corrections

        # ticket 06: review corrections (blind codes in review/corrections.csv) -> exclusions with reason
        corr = corrections.read(Path(project_dir), s)
        layout, review_excluded = corrections.exclusions(v.layout.drop(columns="exclusion_reason"), corr)
        excluded = pd.concat([v.excluded, review_excluded], ignore_index=True) if len(review_excluded) else v.excluded
        atomic_write(Path(project_dir) / layout_mod.EXCLUDED,
                     excluded.to_csv(index=False, sep=s["project"]["csv_delimiter"]))
        measurements = self.measure(project_dir, layout, s, options)
        corrections.write_used(Path(project_dir), corr, list(layout["tumour_id"]))
        print(f"Review: {corr.describe()}")
        table = self.readout_table(project_dir, layout, measurements, s)
        self.analyse(project_dir, table, s, options)

    def measure(self, project_dir: Path, layout: pd.DataFrame, settings: dict, options: RunOptions) -> pd.DataFrame:
        from . import corrections, legacy, measure, qc

        m = measure.measure(project_dir, layout, settings, resegment=options.resegment)
        # ticket 06: re-picked main masses (relabelled cached labels), QC flags, blinded review sheets
        m = qc.measure(project_dir, layout, settings, m, corrections.read(Path(project_dir), settings),
                       render=not options.no_review)
        fiji = legacy.measure(project_dir, layout, settings, m)  # ticket 07: fiji_* legacy readouts
        return m.merge(fiji, on="tumour_id", how="left") if len(fiji.columns) > 1 else m

    def readout_table(self, project_dir: Path, layout: pd.DataFrame, measurements: pd.DataFrame,
                      settings: dict) -> pd.DataFrame:
        from . import table

        t = table.build(layout, measurements, settings)
        table.write(project_dir, t, settings)
        return t

    def analyse(self, project_dir: Path, table: pd.DataFrame, settings: dict, options: RunOptions) -> None:
        from . import results

        written = results.analyse(Path(project_dir), table, settings)
        from .config import factors
        from .plots import shape_panel  # ticket 05: exploratory shape readouts + area vs perimeter

        written.append(str(shape_panel(table, factors(settings), Path(project_dir) / "plots" / "shape_panel")))
        print(f"Wrote {Path(project_dir) / 'tumours.csv'}, " + ", ".join(written))
