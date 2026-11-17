"""Synthetic projects with unusual designs, run through the whole chain (docs/design.md, "Unusual designs").

`project()` writes raw files, layout.csv, qc_details.csv and counts.csv as `kerncount run` would, so that
`chain()` can run everything after detection exactly as the CLI and the browser do: regenerate (design check,
duplicate check, generations: models, stats, readout ladder) -> sampling + planner -> plate maps -> the GUI
endpoints that render the results."""
from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd

from kerncount import discover, gui, master, pipeline, platemap, settings, settings_form

BASE = {"Control": (60, 0.40), "LyoVec": (58, 0.40), "CL075": (45, 0.25), "miR-1": (30, 0.20), "miR-2": (40, 0.30)}
FOLDER = {"Control": "cntrl", "LyoVec": "lyovec", "CL075": "cl075"}
MEASURED = ["readout_frac_positive_Ki67", "intensity_Ki67_field_raw", "intensity_Ki67_nucleus_mean_bgsub"]


def rows(ns=("n1", "n2", "n3"), treatments=("Control", "LyoVec", "miR-1"), images=range(1, 13), seed=0,
         cell_line="") -> list[dict]:
    """A plain design: every n x treatment x image, counts and % Ki-67-positive around BASE."""
    rng = np.random.default_rng(seed)
    out = []
    for k, n in enumerate(ns):
        for t in treatments:
            c, f = BASE.get(t, (50, 0.3))
            for i in images:
                out.append(dict(n=n, treatment=t, image=int(i), count=int(max(1, rng.normal(c * (1 + 0.1 * k), 6))),
                                frac=float(np.clip(rng.normal(f, 0.05), 0, 1)), cell_line=cell_line))
    return out


def _n_folder(n: str) -> str:
    return f"0{int(n[1:]) % 10}.06.26 {n} 72h"


def project(root: Path, data: list[dict], settings_over: dict | None = None, name: str = "E",
            empty: tuple = (), copies: dict | None = None) -> Path:
    """Write a project as `kerncount run` leaves it (without detection).
    data: rows(...) dicts; optional keys qc ("high"/"low"), frac/count = nan for missing readouts.
    empty: (n, treatment) folders with no images; copies: {(n, treatment): (n, treatment)} = byte-identical copy."""
    copies = copies or {}
    raw = root / "data" / name
    prj = root / "analysis" / name
    prj.mkdir(parents=True)
    d = pd.DataFrame(data)
    for r in d.itertuples():
        folder = raw / _n_folder(r.n) / FOLDER.get(r.treatment, r.treatment)
        folder.mkdir(parents=True, exist_ok=True)
        src = copies.get((r.n, r.treatment), (r.n, r.treatment))  # a copied folder: the other condition's bytes
        (folder / f"{r.image}.tif").write_bytes(f"{src[0]}/{src[1]}/{r.image}".encode())
    for n, t in empty:
        (raw / _n_folder(n) / FOLDER.get(t, t)).mkdir(parents=True, exist_ok=True)
    settings.write_template(prj, raw)
    over = {"analysis": {"reference": "LyoVec", "normalise_to": ["LyoVec", "Control"]}}
    for k, v in (settings_over or {}).items():
        over.setdefault(k, {}).update(v)
    sp = prj / "settings.toml"
    sp.write_text(settings_form.update_text(sp.read_text(), settings_form.coerce(over)))
    s = settings.load(prj)
    lay = discover.scan(raw)
    if "cell_line" in d and d["cell_line"].astype(bool).any():
        lines = d.groupby("n")["cell_line"].first()
        lay["cell_line"] = lay["n"].map(lambda x: lines.get(x.split("_")[0], ""))
    lay["n"] = lay["n"].str.split("_").str[0]  # n1_2026-06-01 -> n1 (the n labels of `data`)
    discover.write_layout(lay, prj / "layout.csv")
    qc = pd.DataFrame(dict(
        n=d["n"], treatment=d["treatment"], image=d["image"], auto_count=d["count"], count_low=d["count"],
        count_high=d["count"], qc_level=d["qc"] if "qc" in d else "high", qc_ok="yes", qc_flags="",
        source_file=[str(raw / _n_folder(r.n) / FOLDER.get(r.treatment, r.treatment) / f"{r.image}.tif")
                     for r in d.itertuples()],
        readout_frac_positive_Ki67=d["frac"], n_measured_Ki67=d["count"],
        intensity_Ki67_field_raw=200 + 300 * d["frac"], intensity_Ki67_nucleus_mean_bgsub=50 + 400 * d["frac"]))
    qc = master.order_rows(qc, s["analysis"]["treatment_order"])
    qc.to_csv(prj / "qc_details.csv", index=False)
    master.merge(prj / "counts.csv", qc, s["analysis"]["treatment_order"], ",", prj / "_work" / "backups")
    return prj


def chain(prj: Path) -> dict:
    """Everything after detection, as the CLI / browser run it. Returns every result plus the log text."""
    out = {}
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        pipeline.regenerate(prj)
        try:
            pipeline.sampling(prj, iterations=20)
        except SystemExit as e:  # e.g. no generation
            out["sampling_exit"] = str(e)
        out["platemap"] = platemap.plan_project(prj, plates=2, seed=1)
    out["log"] = buf.getvalue()
    gen = gui.generations(prj)[-1]
    out["gen"] = prj / "generations" / gen
    out["overview"] = gui.overview(prj)
    out["project"] = gui.project_info(prj)
    out["review"] = gui.review_index(prj)
    out["counts"] = gui.counts_data(prj)
    out["manifest"] = gui.manifest(prj, gen)
    out["stats"] = gui.stats_long(prj, gen)
    lines = gui.plan_lines(prj, gen) or [dict(path=None)]
    out["planner"] = [gui.planner_info(prj, None, ln.get("path")) for ln in lines]
    out["plans"] = gui.plans_info(prj)
    out["design"] = pd.read_csv(prj / "design_check.csv", keep_default_na=False)
    out["stats_long"] = pd.read_csv(out["gen"] / "stats_long.csv")
    out["notes"] = out["manifest"].get("notes", [])
    json.dumps({k: v for k, v in out.items() if k in ("overview", "manifest", "stats", "planner", "plans")},
               default=str)  # every endpoint answer is serialisable
    return out
