"""kerncount: automated nucleus counting and analysis.

Workspace convention (recommended):
    <workspace>/data/<experiment>/       raw images (input, never modified)
    <workspace>/analysis/<experiment>/   project: settings, counts.csv, review images, results
Run commands from the workspace folder and refer to experiments by name:
    kerncount init my_experiment    kerncount run my_experiment    kerncount run --all
Explicit paths also work:  kerncount init RAW_DIR PROJECT_DIR  /  kerncount run PROJECT_DIR
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

DATA, ANALYSIS = "data", "analysis"


def workspace_root(start: Path | None = None) -> Path | None:
    """Nearest folder (cwd or a parent) that contains a data/ folder."""
    env = os.environ.get("KERNCOUNT_WORKSPACE")
    if env:
        return Path(env).resolve()
    p = (start or Path.cwd()).resolve()
    for d in [p, *p.parents]:
        if (d / DATA).is_dir():
            return d
    return None


def _experiments(ws: Path) -> list[str]:
    return sorted(d.name for d in (ws / DATA).iterdir() if d.is_dir() and not d.name.startswith("."))


def _bare(arg: str) -> bool:
    """An experiment name (no folder separators, not absolute) rather than a path."""
    return not Path(arg).is_absolute() and "/" not in arg and "\\" not in arg and arg not in (".", "..")


def resolve_raw(arg: str) -> Path:
    """A bare name is looked up in <workspace>/data/ FIRST, then as a folder in the working directory; a path is
    used as given. Symlinks are kept (absolute path, not resolved): data/<exp> may link to a folder elsewhere and
    the project should keep pointing at the link."""
    ws = workspace_root()
    if _bare(arg) and ws and (ws / DATA / arg).is_dir():
        return Path(os.path.abspath(ws / DATA / arg))
    p = Path(arg)
    if p.is_dir():
        return Path(os.path.abspath(p))
    sys.exit(f"raw folder not found: {arg}" + (f" (also looked in {ws / DATA})" if ws else ""))


def resolve_project(arg: str, must_exist: bool = True) -> Path:
    p = Path(arg)
    if (p / "settings.toml").exists() or (not must_exist and p.parent.is_dir() and len(p.parts) > 1):
        return p.resolve()
    ws = workspace_root()
    if ws and ((ws / ANALYSIS / arg / "settings.toml").exists() or not must_exist):
        return (ws / ANALYSIS / arg).resolve()
    if p.is_dir():
        return p.resolve()
    sys.exit(f"project not found: {arg}" + (f" (also looked in {ws / ANALYSIS})" if ws else "")
             + " - create it with `kerncount init`")


def _experiment_type(s: dict):
    """The project's [experiment] type (ADR 0004); an unknown name stops with the list of valid ones."""
    from . import experiment_types

    try:
        return experiment_types.of(s)
    except experiment_types.UnknownExperimentType as e:
        sys.exit(str(e))


class TypeMismatch(Exception):
    """`init --type X` on a project whose settings.toml has another type (without --force)."""


def _init_one(raw: Path, prj: Path, force: bool, force_layout: bool = False, etype_name: str | None = None):
    from . import discover, experiment_types, settings

    new_type = None
    if etype_name:
        try:
            new_type = experiment_types.get(etype_name)
        except experiment_types.UnknownExperimentType as e:
            sys.exit(str(e).replace("settings.toml [experiment] type", "--type"))
    prj.mkdir(parents=True, exist_ok=True)
    if (prj / "settings.toml").exists() and not force:
        current = _experiment_type(settings.load(prj))
        if new_type is not None and new_type.name != current.name:
            raise TypeMismatch(f"settings.toml in {prj} is for type {current.name!r}; use --force to reset it to "
                               f"--type {new_type.name}")
        print(f"settings.toml already exists in {prj} (kept). Use --force to reset it.")
    else:
        had = (prj / "settings.toml").exists()
        template = new_type.settings_template() if new_type is not None else None
        if new_type is not None and template is None:  # a type without its own template: default + type line
            template = settings.TEMPLATE.replace('type = "cell_count"', f'type = "{new_type.name}"', 1)
        settings.write_template(prj, raw, template)  # backs up an existing settings.toml to _work/backups/ first
        if had:
            print(f"settings.toml reset to the defaults; the previous version is kept in {settings.BACKUPS}/")
    s = settings.load(prj)
    etype = _experiment_type(s)
    layout = etype.scan_layout(raw, s)
    lp = prj / "layout.csv"
    merged = etype.merge_layout(discover.read_table(lp), layout, s) if lp.exists() and not (force or force_layout) \
        else None
    if merged is not None:  # raw folder reorganised: layout.csv updated with the user's edits, old one backed up
        layout, messages = merged
        bk = settings.backup(prj, "layout.csv")
        print("\n".join(messages))
        kept = bk.relative_to(prj).as_posix() if bk is not None else "-"
        print(f"layout.csv updated; the previous version is kept in {kept}")
    elif lp.exists() and not (force or force_layout):
        lp = prj / "layout.new.csv"
        print(f"layout.csv already exists (kept); fresh scan written to {lp.name} for comparison.")
    elif lp.exists():
        bk = settings.backup(prj, "layout.csv")
        print(f"layout.csv overwritten; the previous version is kept in {bk.relative_to(prj).as_posix()}")
    discover.write_layout(layout, lp, s["project"]["csv_delimiter"])
    etype.report_layout(layout, s, raw)
    ws = workspace_root()
    target = prj.name if ws and prj.resolve().parent == (ws / ANALYSIS).resolve() else str(prj)
    hint = getattr(etype, "layout_hint", "Check the n / treatment names (edit the CSV if needed)")
    print(f"\nWrote {lp}\n{hint}, then run:  kerncount run \"{target}\"")


def cmd_init(a):
    if a.all:
        ws = workspace_root() or sys.exit("--all needs a workspace: run from a folder that contains data/")
        skipped = []
        for exp in _experiments(ws):
            print(f"\n=== {exp}")
            try:
                _init_one(ws / DATA / exp, ws / ANALYSIS / exp, a.force, a.force_layout, a.type)
            except TypeMismatch as e:  # report it and go on with the other experiments
                print(f"skipped: {e}")
                skipped.append(exp)
        if skipped:
            print(f"\nNot initialised as --type {a.type} (settings.toml has another type): {', '.join(skipped)}")
        return
    if not a.raw_dir:
        sys.exit("give an experiment name (from data/), a raw folder, or --all")
    raw = resolve_raw(a.raw_dir)
    if a.project_dir:
        prj = Path(a.project_dir).resolve()
    else:
        ws = workspace_root()
        # data/<name> may be a symlink to a folder elsewhere: the project still belongs in this workspace
        in_ws = ws and (raw.parent == (ws / DATA).resolve() or (ws / DATA / raw.name).resolve() == raw)
        prj = (ws / ANALYSIS / raw.name) if in_ws else raw.parent.parent / ANALYSIS / raw.name
    try:
        _init_one(raw, prj, a.force, a.force_layout, a.type)
    except TypeMismatch as e:
        sys.exit(str(e))


def _projects(a) -> list[Path]:
    if a.all:
        ws = workspace_root() or sys.exit("--all needs a workspace: run from a folder that contains data/ and analysis/")
        prjs = [ws / ANALYSIS / e for e in _experiments(ws) if (ws / ANALYSIS / e / "settings.toml").exists()]
        missing = [e for e in _experiments(ws) if not (ws / ANALYSIS / e / "settings.toml").exists()]
        if missing:
            print(f"not initialised yet (run `kerncount init <name>` first): {', '.join(missing)}")
        return prjs
    if not a.project:
        sys.exit("give an experiment/project name or --all")
    return [resolve_project(a.project)]


def cmd_run(a):
    from . import experiment_types, settings

    opts = experiment_types.RunOptions(rebuild=a.rebuild, resegment=a.resegment, no_review=a.no_review,
                                       unlock=a.unlock)
    for prj in _projects(a):
        print(f"\n=== {prj.name}")
        _experiment_type(settings.load(prj)).run(prj, opts)


def cmd_regen(a):
    from . import pipeline

    for prj in _projects(a):
        print(f"\n=== {prj.name}")
        pipeline.regenerate(prj, rebuild=a.rebuild)


def cmd_sampling(a):
    from . import pipeline

    for prj in _projects(a):
        print(f"\n=== {prj.name}")
        pipeline.sampling(prj, a.generation, a.iterations)


def cmd_wells_template(a):
    from . import pipeline

    for prj in _projects(a):
        print(f"\n=== {prj.name}")
        pipeline.wells_template(prj, overwrite=a.force)


def cmd_platemap(a):
    from . import platemap

    for prj in _projects(a):
        print(f"\n=== {prj.name}")
        try:
            out = platemap.plan_project(prj, plates=a.plates, plate_format=a.format, wells=a.wells,
                                        reference_wells=a.reference_wells, avoid_edge=True if a.avoid_edge else None,
                                        seed=a.seed, balance=not a.no_balance, negative_control=a.negative_control,
                                        negative_control_wells=a.control_wells)
        except ValueError as e:
            sys.exit(f"platemap: {e}")
        print((out / "README.txt").read_text(encoding="utf-8"))
        print(f"Written to {out}")


def _ref_wells(v: str):
    return v if v.strip().lower() == "auto" else int(v)


def cmd_gui(a):
    from . import gui

    gui.serve(a.workspace, a.port, not a.no_browser)


def main(argv=None):
    p = argparse.ArgumentParser(prog="kerncount", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("init", help="scan raw images and create a project (analysis/<experiment>)")
    i.add_argument("raw_dir", nargs="?", help="experiment name in data/, or a raw folder path")
    i.add_argument("project_dir", nargs="?", help="optional; default analysis/<experiment>")
    i.add_argument("--all", action="store_true", help="every experiment in data/")
    i.add_argument("--force", action="store_true", help="overwrite settings.toml and layout.csv")
    i.add_argument("--force-layout", action="store_true",
                   help="keep settings.toml, rewrite layout.csv (e.g. after setting [channels] files = \"interleaved\")")
    i.add_argument("--type", default=None, metavar="NAME",
                   help="experiment type of a new project, e.g. tumour_volume (default: cell_count)")
    i.set_defaults(func=cmd_init)
    r = sub.add_parser("run", help="count, check, and analyse (safe to re-run any time)")
    r.add_argument("project", nargs="?", help="experiment name, or a project folder path")
    r.add_argument("--all", action="store_true", help="every initialised experiment")
    r.add_argument("--rebuild", nargs="+", metavar="gN", default=[],
                   help="overwrite existing generation folders, e.g. --rebuild g2  or  --rebuild all")
    r.add_argument("--resegment", action="store_true", help="ignore cached segmentations")
    r.add_argument("--no-review", action="store_true", help="skip drawing review images (faster)")
    r.add_argument("--unlock", default=None, metavar="REASON",
                   help="tumour experiments: unlock a changed analysis plan (logged with this reason in lock_log.csv)")
    r.set_defaults(func=cmd_run)
    rg = sub.add_parser("regen", help="rebuild generations from counts.csv only (fast; after manual corrections)")
    rg.add_argument("project", nargs="?", help="experiment name, or a project folder path")
    rg.add_argument("--all", action="store_true", help="every initialised experiment")
    rg.add_argument("--rebuild", nargs="+", metavar="gN", default=[],
                    help="overwrite existing generation folders, e.g. --rebuild g2 g3  or  --rebuild all")
    rg.set_defaults(func=cmd_regen)
    m = sub.add_parser("sampling", help="how many images per condition are needed? (variance + power + subsampling)")
    m.add_argument("project", nargs="?", help="experiment name, or a project folder path")
    m.add_argument("--all", action="store_true", help="every initialised experiment")
    m.add_argument("--generation", default=None, help="e.g. g2 (default: latest)")
    m.add_argument("--iterations", type=int, default=1000)
    m.set_defaults(func=cmd_sampling)
    wt = sub.add_parser("wells-template", help="write wells.csv (which well / field each image is) to edit by hand")
    wt.add_argument("project", nargs="?", help="experiment name, or a project folder path")
    wt.add_argument("--all", action="store_true", help="every initialised experiment")
    wt.add_argument("--force", action="store_true", help="overwrite an existing wells.csv (default: write wells.new.csv)")
    wt.set_defaults(func=cmd_wells_template)
    pm = sub.add_parser("platemap", help="randomised plate maps for the next experiment (one plate per n)")
    pm.add_argument("project", nargs="?", help="experiment name, or a project folder path")
    pm.add_argument("--all", action="store_true", help="every initialised experiment")
    pm.add_argument("--plates", type=int, default=3, help="number of plates (= n's) to plan (default 3)")
    pm.add_argument("--format", type=int, default=None, choices=[6, 12, 24, 48, 96],
                    help="wells per plate (default: [planner] plate_format)")
    pm.add_argument("--wells", type=int, default=2, help="wells per treatment condition (default 2)")
    pm.add_argument("--reference-wells", type=_ref_wells, default=None,
                    help='wells per reference condition, or "auto" (Dunnett sqrt(k)); default: [planner] reference_wells')
    pm.add_argument("--avoid-edge", action="store_true",
                    help="inner wells first (default: [planner] avoid_edge)")
    pm.add_argument("--seed", type=int, default=None, help="random seed (default: new; printed for reuse)")
    pm.add_argument("--no-balance", action="store_true", help="independent random plates (no rotation across plates)")
    pm.add_argument("--negative-control", default=None, metavar="NAME",
                    help='reserve wells for a negative-control condition (e.g. "secondary only"); '
                         'default: [experiment] negative_control; "" = none')
    pm.add_argument("--control-wells", type=int, default=2, choices=[1, 2],
                    help="wells per plate for the negative control (default 2)")
    pm.set_defaults(func=cmd_platemap)
    g = sub.add_parser("gui",help="open the point-and-click web interface in your browser")
    g.add_argument("workspace", nargs="?", help="workspace folder (default: found automatically)")
    g.add_argument("--port", type=int, default=8765)
    g.add_argument("--no-browser", action="store_true", help="don't open a browser tab")
    g.set_defaults(func=cmd_gui)
    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
