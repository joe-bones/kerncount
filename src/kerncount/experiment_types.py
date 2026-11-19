"""Experiment types (ADR 0004): `[experiment] type` in settings.toml selects how an experiment is

1. discovered    scan_layout:   raw folder -> layout.csv rows (one row per unit the user checks once)
2. measured      measure:       project + layout -> raw measurements per measured unit
3. tabulated     readout_table: measurements -> the readout table (one row per datapoint, with its layout metadata)
4. analysed      analyse:       readout table -> statistics, plots and results for the type's statistical design

Keep this module importable without numpy / pandas (settings_form lists the types).
Everything else (data/ never modified, analysis/<experiment>/, settings.toml, layout.csv, backups, GUI) is shared.
`kerncount init` calls scan_layout (+ report_layout), `kerncount run` calls run(), which chains 2-4.
Omitted type = "cell_count", the original nucleus-count / intensity pipeline, unchanged.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar

if TYPE_CHECKING:  # pandas only at call time: settings_form (stdlib only) imports this registry
    import pandas as pd

DEFAULT = "cell_count"


@dataclass(frozen=True)
class RunOptions:
    """`kerncount run` flags; a type ignores the ones that do not apply to it."""
    rebuild: list[str] = field(default_factory=list)
    resegment: bool = False  # ignore cached segmentations / measurements
    no_review: bool = False  # skip drawing review images


class UnknownExperimentType(ValueError):
    pass


class ExperimentType:
    """Base class: subclass, set `name` / `label`, implement the four seams and decorate with @register."""
    name: ClassVar[str]
    label: ClassVar[str]

    def scan_layout(self, raw: Path, settings: dict) -> pd.DataFrame:
        """Layout parsing: read the raw folder (never modify it) into layout.csv rows."""
        raise NotImplementedError

    def report_layout(self, layout: pd.DataFrame, settings: dict, raw: Path) -> None:
        """Print what `init` found so the user can check it."""
        import pandas as pd

        with pd.option_context("display.width", 200, "display.max_colwidth", 60, "display.max_rows", 500):
            print(layout.to_string(index=False))

    def measure(self, project_dir: Path, layout: pd.DataFrame, settings: dict, options: RunOptions) -> pd.DataFrame:
        """Measurement: one row per measured unit (cached under _work/ where expensive)."""
        raise NotImplementedError

    def readout_table(self, project_dir: Path, layout: pd.DataFrame, measurements: pd.DataFrame,
                      settings: dict) -> pd.DataFrame:
        """Readout table: one row per datapoint with its layout metadata; written into the project."""
        raise NotImplementedError

    def analyse(self, project_dir: Path, table: pd.DataFrame, settings: dict, options: RunOptions) -> None:
        """Statistical design: tests, plots and results from the readout table."""
        raise NotImplementedError

    def run(self, project_dir: Path, options: RunOptions) -> None:
        """`kerncount run`: layout -> measure -> readout table -> analyse."""
        from . import discover, settings as settings_mod

        s = settings_mod.load(project_dir)
        layout_path = Path(project_dir) / "layout.csv"
        if not layout_path.exists():
            raise SystemExit(f"{layout_path} not found - run 'kerncount init' first")
        layout = discover.read_table(layout_path)
        measurements = self.measure(project_dir, layout, s, options)
        table = self.readout_table(project_dir, layout, measurements, s)
        self.analyse(project_dir, table, s, options)


TYPES: dict[str, ExperimentType] = {}


def register(cls: type[ExperimentType]) -> type[ExperimentType]:
    TYPES[cls.name] = cls()
    return cls


def get(name: str) -> ExperimentType:
    key = str(name).strip() or DEFAULT
    if key not in TYPES:
        raise UnknownExperimentType(f"settings.toml [experiment] type = {name!r} is not a known experiment type; "
                                    f"valid types: {', '.join(sorted(TYPES))}")
    return TYPES[key]


def of(settings: dict) -> ExperimentType:
    """The experiment type a project's settings select (omitted = cell_count)."""
    return get((settings.get("experiment") or {}).get("type") or DEFAULT)


@register
class CellCount(ExperimentType):
    """Nucleus count + marker intensity on plates (field ⊂ well ⊂ condition ⊂ n = plate). Its measurement, readout
    table (qc_details.csv / counts.csv with manual-correction generations) and statistics are interleaved inside
    pipeline.run, so run() delegates to it as a whole: wrapped, not split, to keep its behaviour unchanged."""
    name = "cell_count"
    label = "Cell count / intensity (plates)"

    def scan_layout(self, raw: Path, settings: dict) -> pd.DataFrame:
        from . import discover, settings as settings_mod

        return discover.scan(raw, settings["project"]["date_order"], files=settings["channels"]["files"],
                             n_channels=settings_mod.n_channels_expected(settings) or 1)

    def report_layout(self, layout: pd.DataFrame, settings: dict, raw: Path) -> None:
        import pandas as pd

        from . import discover, settings as settings_mod

        with pd.option_context("display.width", 200, "display.max_colwidth", 60, "display.max_rows", 500):
            print(layout[["n", "treatment", "treatment_folder", "images", "note"]].to_string(index=False))
        if settings_mod.interleaved(settings):
            print(f"\nChannels are separate files ({', '.join(settings['channels']['file_order'])} in file order): "
                  "'images' counts fields (file groups).")
            return
        try:
            det = discover.detect_interleaved(raw)
        except Exception:  # detection is only a hint
            det = {"interleaved": False}
        if det["interleaved"]:
            first, second = ("marker", "DAPI") if det["bright_position"] == 1 else ("DAPI", "marker")
            print("\nHint: the files look like ONE channel each, alternating dim / bright (e.g. marker, DAPI, marker, "
                  "DAPI, ...).\nIf so, set in settings.toml and run `kerncount init` again with --force-layout:\n"
                  f"  [channels]\n  files = \"interleaved\"\n  file_order = [\"{first}\", \"{second}\"]   "
                  "# your channel names, in file order\n  segment = \"DAPI\"")

    def run(self, project_dir: Path, options: RunOptions) -> None:
        from . import pipeline

        pipeline.run(project_dir, rebuild=options.rebuild, resegment=options.resegment, no_review=options.no_review)
