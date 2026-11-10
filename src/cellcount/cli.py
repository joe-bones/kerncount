"""Command line: `cellcount init RAW_DIR PROJECT_DIR` then `cellcount run PROJECT_DIR`."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def cmd_init(a):
    from . import discover, settings

    raw, prj = Path(a.raw_dir).resolve(), Path(a.project_dir).resolve()
    if not raw.is_dir():
        sys.exit(f"raw folder not found: {raw}")
    prj.mkdir(parents=True, exist_ok=True)
    if (prj / "settings.toml").exists() and not a.force:
        print(f"settings.toml already exists in {prj} (kept). Use --force to reset it.")
    else:
        settings.write_template(prj, raw)
    s = settings.load(prj)
    layout = discover.scan(raw, s["project"]["date_order"])
    lp = prj / "layout.csv"
    if lp.exists() and not a.force:
        lp = prj / "layout.new.csv"
        print(f"layout.csv already exists (kept); fresh scan written to {lp.name} for comparison.")
    discover.write_layout(layout, lp, s["project"]["csv_delimiter"])
    with __import__("pandas").option_context("display.width", 200, "display.max_colwidth", 60, "display.max_rows", 500):
        print(layout[["n", "treatment", "treatment_folder", "images", "note"]].to_string(index=False))
    print(f"\nWrote {lp}\nCheck the n / treatment names (edit the CSV if needed), then run:\n  cellcount run \"{prj}\"")


def cmd_run(a):
    from . import pipeline

    pipeline.run(Path(a.project_dir), rebuild=a.rebuild, resegment=a.resegment, no_review=a.no_review)


def cmd_sampling(a):
    from . import pipeline

    pipeline.sampling(Path(a.project_dir), a.generation, a.iterations)


def main(argv=None):
    p = argparse.ArgumentParser(prog="cellcount", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("init", help="scan a raw image folder and create a project folder")
    i.add_argument("raw_dir")
    i.add_argument("project_dir")
    i.add_argument("--force", action="store_true", help="overwrite settings.toml and layout.csv")
    i.set_defaults(func=cmd_init)
    r = sub.add_parser("run", help="count, check, and analyse (safe to re-run any time)")
    r.add_argument("project_dir")
    r.add_argument("--rebuild", nargs="+", metavar="gN", default=[],
                   help="overwrite existing generation folders, e.g. --rebuild g2  or  --rebuild all")
    r.add_argument("--resegment", action="store_true", help="ignore cached segmentations")
    r.add_argument("--no-review", action="store_true", help="skip drawing review images (faster)")
    r.set_defaults(func=cmd_run)
    m = sub.add_parser("sampling", help="how many images per condition are needed? (variance + power + subsampling)")
    m.add_argument("project_dir")
    m.add_argument("--generation", default=None, help="e.g. g2 (default: latest)")
    m.add_argument("--iterations", type=int, default=1000)
    m.set_defaults(func=cmd_sampling)
    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
