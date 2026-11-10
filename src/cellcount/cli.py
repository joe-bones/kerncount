"""cellcount: automated nucleus counting and analysis.

Workspace convention (recommended):
    <workspace>/data/<experiment>/       raw images (input, never modified)
    <workspace>/analysis/<experiment>/   project: settings, counts.csv, review images, results
Run commands from the workspace folder and refer to experiments by name:
    cellcount init DAPI_output      cellcount run DAPI_output      cellcount run --all
Explicit paths also work:  cellcount init RAW_DIR PROJECT_DIR  /  cellcount run PROJECT_DIR
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

DATA, ANALYSIS = "data", "analysis"


def workspace_root(start: Path | None = None) -> Path | None:
    """Nearest folder (cwd or a parent) that contains a data/ folder."""
    env = os.environ.get("CELLCOUNT_WORKSPACE")
    if env:
        return Path(env).resolve()
    p = (start or Path.cwd()).resolve()
    for d in [p, *p.parents]:
        if (d / DATA).is_dir():
            return d
    return None


def _experiments(ws: Path) -> list[str]:
    return sorted(d.name for d in (ws / DATA).iterdir() if d.is_dir() and not d.name.startswith("."))


def resolve_raw(arg: str) -> Path:
    p = Path(arg)
    if p.is_dir():
        return p.resolve()
    ws = workspace_root()
    if ws and (ws / DATA / arg).is_dir():
        return (ws / DATA / arg).resolve()
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
             + " - create it with `cellcount init`")


def _init_one(raw: Path, prj: Path, force: bool):
    from . import discover, settings

    prj.mkdir(parents=True, exist_ok=True)
    if (prj / "settings.toml").exists() and not force:
        print(f"settings.toml already exists in {prj} (kept). Use --force to reset it.")
    else:
        settings.write_template(prj, raw)
    s = settings.load(prj)
    layout = discover.scan(raw, s["project"]["date_order"])
    lp = prj / "layout.csv"
    if lp.exists() and not force:
        lp = prj / "layout.new.csv"
        print(f"layout.csv already exists (kept); fresh scan written to {lp.name} for comparison.")
    discover.write_layout(layout, lp, s["project"]["csv_delimiter"])
    with __import__("pandas").option_context("display.width", 200, "display.max_colwidth", 60, "display.max_rows", 500):
        print(layout[["n", "treatment", "treatment_folder", "images", "note"]].to_string(index=False))
    ws = workspace_root()
    target = prj.name if ws and prj.resolve().parent == (ws / ANALYSIS).resolve() else str(prj)
    print(f"\nWrote {lp}\nCheck the n / treatment names (edit the CSV if needed), then run:  cellcount run \"{target}\"")


def cmd_init(a):
    if a.all:
        ws = workspace_root() or sys.exit("--all needs a workspace: run from a folder that contains data/")
        for exp in _experiments(ws):
            print(f"\n=== {exp}")
            _init_one(ws / DATA / exp, ws / ANALYSIS / exp, a.force)
        return
    if not a.raw_dir:
        sys.exit("give an experiment name (from data/), a raw folder, or --all")
    raw = resolve_raw(a.raw_dir)
    if a.project_dir:
        prj = Path(a.project_dir).resolve()
    else:
        ws = workspace_root()
        prj = (ws / ANALYSIS / raw.name) if ws and raw.parent == (ws / DATA).resolve() else raw.parent.parent / ANALYSIS / raw.name
    _init_one(raw, prj, a.force)


def _projects(a) -> list[Path]:
    if a.all:
        ws = workspace_root() or sys.exit("--all needs a workspace: run from a folder that contains data/ and analysis/")
        prjs = [ws / ANALYSIS / e for e in _experiments(ws) if (ws / ANALYSIS / e / "settings.toml").exists()]
        missing = [e for e in _experiments(ws) if not (ws / ANALYSIS / e / "settings.toml").exists()]
        if missing:
            print(f"not initialised yet (run `cellcount init <name>` first): {', '.join(missing)}")
        return prjs
    if not a.project:
        sys.exit("give an experiment/project name or --all")
    return [resolve_project(a.project)]


def cmd_run(a):
    from . import pipeline

    for prj in _projects(a):
        print(f"\n=== {prj.name}")
        pipeline.run(prj, rebuild=a.rebuild, resegment=a.resegment, no_review=a.no_review)


def cmd_sampling(a):
    from . import pipeline

    for prj in _projects(a):
        print(f"\n=== {prj.name}")
        pipeline.sampling(prj, a.generation, a.iterations)


def main(argv=None):
    p = argparse.ArgumentParser(prog="cellcount", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("init", help="scan raw images and create a project (analysis/<experiment>)")
    i.add_argument("raw_dir", nargs="?", help="experiment name in data/, or a raw folder path")
    i.add_argument("project_dir", nargs="?", help="optional; default analysis/<experiment>")
    i.add_argument("--all", action="store_true", help="every experiment in data/")
    i.add_argument("--force", action="store_true", help="overwrite settings.toml and layout.csv")
    i.set_defaults(func=cmd_init)
    r = sub.add_parser("run", help="count, check, and analyse (safe to re-run any time)")
    r.add_argument("project", nargs="?", help="experiment name, or a project folder path")
    r.add_argument("--all", action="store_true", help="every initialised experiment")
    r.add_argument("--rebuild", nargs="+", metavar="gN", default=[],
                   help="overwrite existing generation folders, e.g. --rebuild g2  or  --rebuild all")
    r.add_argument("--resegment", action="store_true", help="ignore cached segmentations")
    r.add_argument("--no-review", action="store_true", help="skip drawing review images (faster)")
    r.set_defaults(func=cmd_run)
    m = sub.add_parser("sampling", help="how many images per condition are needed? (variance + power + subsampling)")
    m.add_argument("project", nargs="?", help="experiment name, or a project folder path")
    m.add_argument("--all", action="store_true", help="every initialised experiment")
    m.add_argument("--generation", default=None, help="e.g. g2 (default: latest)")
    m.add_argument("--iterations", type=int, default=1000)
    m.set_defaults(func=cmd_sampling)
    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
