"""Tumour statistics outputs (ticket 04): stats_long.csv, results.xlsx, Prism CSVs, plots and stats_report.md.

For every declared readout present in the table: fit the N-way log-scale model (stats.fit), run every family
(stats.readout_results), normalise to the focus factor's reference, and write the outputs. The analysis plan
(factors, readout split, families) is **locked**: its hash is written into every stats row, the workbook and the report.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from . import planner, plots, sensitivity, stats
from .config import Factor, Family, Readout, batch, factors, families, normalise_by, readouts, tumour_settings

MODEL = "N-way linear model on log(readout), full interaction, pooled variance"
BLOCK = "; batch as fixed block effect"  # ticket 13: appended when the model estimates block effects


def model_text(model: stats.CellModel) -> str:
    return MODEL + (BLOCK + f" ({model.n_blocks} parameter{'s' if model.n_blocks != 1 else ''})"
                    if model.n_blocks else "")
NAN = float("nan")


# [tumour] settings that define the confirmatory values (segmentation + mid-section rule): part of the locked plan
# (review fix 2026-11-20; ADR 0003 "locked in settings before results are viewed"). Exploratory-only scales
# (roughness, slab, surface grid) are not.
MEASUREMENT_KEYS = ("smooth_sigma_um", "threshold", "min_object_um3", "mid_section_um", "outline_smooth_um")


def _measurement(settings: dict) -> dict:
    from .config import tumour_settings

    ts = tumour_settings(settings)
    return {k: (float(ts[k]) if isinstance(ts[k], (int, float)) and not isinstance(ts[k], bool) else str(ts[k]))
            for k in MEASUREMENT_KEYS}


def lock(settings: dict) -> tuple[str, dict]:
    """(12-hex hash, spec) of the analysis plan: factors (levels, reference), readouts (column, role), families,
    the measurement settings behind the confirmatory readouts (MEASUREMENT_KEYS) and the batch rule."""
    spec = {
        "factors": [{"name": f.name, "levels": list(f.levels), "reference": f.reference} for f in factors(settings)],
        "readouts": [{"name": r.name, "column": r.column, "role": r.role} for r in readouts(settings)],
        "families": [{"id": f.id, "effect": f.effect, "factors": list(f.factors), "role": f.role}
                     for f in families(settings)],
        "normalise_by": (normalise_by(settings, factors(settings)) or Factor("", ("",), "")).name,
        "measurement": _measurement(settings),
    }
    if batch(settings) != "session_folder":  # ticket 12: decides which contrasts are tested; in the hash only when
        spec["batch"] = batch(settings)  # not the session folder, so locks made before ticket 13 stay valid
        spec["batch_effect"] = "fixed block"  # ticket 13 (ADR 0003 amendment 2026-11-20)
    return hashlib.sha1(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:12], spec


def _group_label(cell: tuple) -> str:
    return " / ".join(cell)


def _key(facs: list[Factor], cell: tuple) -> tuple:
    by = {f.name: lv for f, lv in zip(facs, cell)}

    def role(r: str):
        f = next((f for f in facs if f.role == r), None)
        return by[f.name] if f is not None else None
    return role("row"), role("panel"), role("x")


def _stats_rows(res: pd.DataFrame, means: pd.DataFrame, facs: list[Factor], readout: Readout, h: str,
                model_text: str = MODEL) -> list[dict]:
    rows = []
    for _, m in means.iterrows():
        cell = tuple(str(m[f.name]) for f in facs)
        rows.append(dict(level="all_n", n=None, metric=readout.name, family="estimated means",
                         test="model cell mean (geometric mean, back-transformed)", comparison=_group_label(cell),
                         estimate=m["mean"], estimate_type="mean", ci_low=m["ci_low"], ci_high=m["ci_high"], p=NAN,
                         p_adj=NAN, significance="", N=m["n"], model=model_text, note=f"{readout.role}; lock {h}",
                         estimand=f"geometric mean {readout.label}", cell_line=""))
    for _, r in res.iterrows():
        tested = np.isfinite(r["p"])
        test = (f"{r['family_role']} contrast, t (df {r['df']}), Holm within family" if tested else
                "not tested (estimate and CI only)" if np.isfinite(r["estimate"]) else "not estimable")
        rows.append(dict(level="all_n", n=None, metric=readout.name, family=r["family"], test=test,
                         comparison=r["comparison"], estimate=r["estimate"], estimate_type="ratio",
                         ci_low=r["ci_low"], ci_high=r["ci_high"], p=r["p"], p_adj=r["p_adj"],
                         significance=r["significance"], N=r["N"], model=model_text,
                         note="; ".join(x for x in (f"{r['family_role']} family", r["note"], f"lock {h}") if x),
                         estimand=("fold-change of geometric means" if " vs " in r["comparison"] and
                                   not r["comparison"].startswith("(") else "ratio of fold-changes")
                         + f" of {readout.label}", cell_line=""))
    return rows


def _prism(table: pd.DataFrame, facs: list[Factor], values: pd.Series) -> pd.DataFrame:
    from ..analysis import prism_table
    from .families import cells

    keys = stats.cell_of(table, facs)
    groups = {}
    for cell in cells(facs):
        v = values[(keys == cell) & values.notna()]
        if len(v):
            groups[_group_label(cell)] = v.tolist()
    return prism_table(groups)


def analyse(project_dir: Path, table: pd.DataFrame, settings: dict) -> list[str]:
    """Run the pre-registered statistics on every declared readout; write all outputs; return report lines."""
    from ..manifest import write_stats_long

    out = Path(project_dir)
    facs = factors(settings)
    focus = normalise_by(settings, facs)
    h, spec = lock(settings)
    from . import lock_log  # ticket 10: a rebuild with another plan (and any unlock) is logged and reported

    lock_log.record_run(out, h)
    names = {f.name for f in facs}
    fams: list[Family] = []
    notes: list[str] = []
    for fam in families(settings):
        missing = [x for x in fam.factors if x not in names]
        if missing:
            notes.append(f"family {fam.id} skipped: factor(s) {', '.join(missing)} not declared")
        else:
            fams.append(fam)
    delim = settings.get("project", {}).get("csv_delimiter", ",")
    tum = table.copy()
    stats_rows: list[dict] = []
    all_res, all_means, written = [], [], []
    n_perm, seed = sensitivity.options(settings)
    all_sens: list[pd.DataFrame] = []
    plans: list[pd.DataFrame] = []
    pooled_sd: dict[str, float] = {}
    blocks_seen: list[stats.CellModel] = []
    for ro in readouts(settings):
        if ro.column not in table.columns:
            notes.append(f"readout {ro.name}: column {ro.column} not in tumours.csv (not analysed)")
            continue
        model = stats.fit(table, facs, ro.column)
        if model.dropped:
            notes.append(f"readout {ro.name}: {model.dropped} tumour(s) without a positive value or outside the "
                         "declared levels left out")
        res = stats.readout_results(table, facs, fams, ro, model)
        means = stats.cell_means(model)
        pct_col = f"{ro.name}_pct"
        if focus is not None:
            tum[pct_col] = stats.normalised(table, model, focus)
            norm = stats.normalised_groups(model, focus)
            means = means.merge(norm[[f.name for f in facs] + ["pct", "ci_low", "ci_high"]].rename(
                columns={"pct": "pct_of_reference", "ci_low": "pct_ci_low", "ci_high": "pct_ci_high"}),
                on=[f.name for f in facs], how="left")
        stats_rows += _stats_rows(res, means, facs, ro, h, model_text(model))
        if model.n_blocks:
            blocks_seen.append(model)
        adjusted = stats.batch_adjusted(table, model)  # ticket 13: sensitivity on within-batch (adjusted) values
        sens = sensitivity.readout_sensitivity(adjusted, facs, fams, ro, res, n_perm, seed)
        stats_rows += sensitivity.long_rows(sens, ro, h)
        all_sens.append(sens)
        if ro.role == "confirmatory":
            fc, pw, alpha = planner.options(settings)
            plans.append(planner.plan_readout(model, facs, fams, ro, fc, pw, alpha, batch(settings)))
            pooled_sd[ro.name] = float(np.sqrt(model.s2))
        written.append(str(sensitivity.residual_plot(sensitivity.groups(adjusted, facs, ro.column), facs, ro.label,
                                                     out / "plots" / f"{ro.name}_residuals")))
        all_res.append(res)
        all_means.append(means.assign(readout=ro.name))
        # Prism-ready per-tumour tables (one column per group)
        (out / "prism").mkdir(parents=True, exist_ok=True)
        vals = pd.to_numeric(table[ro.column], errors="coerce")
        _prism(table, facs, vals).to_csv(out / "prism" / f"{ro.name}_prism.csv", index=False, sep=delim)
        if focus is not None:
            _prism(table, facs, tum[pct_col]).to_csv(out / "prism" / f"{ro.name}_pct_prism.csv", index=False,
                                                     sep=delim)
        written += _plots(out, table, tum, facs, focus, ro, model, res, means, fams)
    write_stats_long(out, stats_rows)
    lines = _about(settings, h, spec, fams, notes, table, facs, blocks_seen[0] if blocks_seen else None)
    lines += lock_log.report_lines(out)
    _excel(out, lines, tum, stats_rows, all_means, fams)
    report = _report(out, lines, all_res)
    _sensitivity_outputs(out, report, all_sens, n_perm, seed)
    _planning_outputs(out, report, plans, pooled_sd, planner.options(settings))
    from . import concordance  # ticket 07: primary vs Fiji legacy, threshold sweep

    written += concordance.analyse(out, table, settings, facs, fams, report)
    from . import acquisition, shared  # ticket 11: acquisition settings; shared-threshold robustness check

    acq_written, consistency = acquisition.analyse(out, table, settings, facs, report)
    written += acq_written
    written += shared.analyse(out, table, settings, facs, fams, report, consistency)
    from . import qc  # ticket 06: flow table, sensitivity without flagged tumours, A-B hemisphere check

    written += qc.analyse(out, table, settings, facs, fams, report, h)
    print("\n".join(x for x in lines if x.startswith(("Tumour statistics", "Readouts", "Levels without", "NOTE"))))
    return [str(out / "stats_long.csv"), str(out / "results.xlsx"), str(report), *written]


def _sensitivity_outputs(out: Path, report: Path, sens: list[pd.DataFrame], n_perm: int, seed: int) -> None:
    """Ticket 08: 'Sensitivity' sheet in results.xlsx and a '## Sensitivity' section at the end of the report."""
    with pd.ExcelWriter(out / "results.xlsx", engine="openpyxl", mode="a") as xw:
        (pd.concat(sens, ignore_index=True) if sens else pd.DataFrame(columns=sensitivity.SENSITIVITY_COLUMNS)
         ).to_excel(xw, sheet_name="Sensitivity", index=False)
    with report.open("a", encoding="utf-8") as fh:
        fh.write("\n" + "\n".join(sensitivity.report_lines(sens, n_perm, seed)) + "\n")


def _planning_outputs(out: Path, report: Path, plans: list[pd.DataFrame], sd: dict[str, float],
                      opts: tuple[float, float, float]) -> None:
    """Ticket 09: 'Planning aid' sheets (family summary + per contrast) and a separate report section."""
    rows = pd.concat(plans, ignore_index=True) if plans else pd.DataFrame(columns=planner.PLAN_COLUMNS)
    summ = planner.summarise(rows)
    with pd.ExcelWriter(out / "results.xlsx", engine="openpyxl", mode="a") as xw:
        summ.to_excel(xw, sheet_name="Planning aid", index=False)
        rows.to_excel(xw, sheet_name="Planning aid (contrasts)", index=False)
    with report.open("a", encoding="utf-8") as fh:
        fh.write("\n" + "\n".join(planner.report_lines(summ, sd, *opts)) + "\n")


def _plots(out: Path, table: pd.DataFrame, tum: pd.DataFrame, facs: list[Factor], focus: Factor | None,
           ro: Readout, model: stats.CellModel, res: pd.DataFrame, means: pd.DataFrame, fams: list[Family]) -> list[str]:
    placed = all(f.role in ("row", "panel", "x") for f in facs)  # every cell is one plotted group
    ci = {}
    norm_ci = {}
    for _, m in means.iterrows():
        cell = tuple(str(m[f.name]) for f in facs)
        ci[_key(facs, cell)] = (m["mean"], m["ci_low"], m["ci_high"])
        if "pct_of_reference" in m:
            norm_ci[_key(facs, cell)] = (m["pct_of_reference"], m["pct_ci_low"], m["pct_ci_high"])
    # marks on the normalised plot: the simple-effect family of the focus factor (e.g. F1 line vs cntrl)
    marks: dict = {}
    fam = next((f for f in fams if focus is not None and f.effect == "simple" and f.factors == (focus.name,)), None)
    if fam is not None and placed:
        from .families import contrasts

        for c, (_, r) in zip(contrasts(facs, fam), res[res["family"] == fam.id].iterrows()):
            cell = next(k for k, w in c.weights.items() if w > 0)
            marks[_key(facs, cell)] = r["significance"] or ("conf." if r["confounded"] else "")
    pd_ = out / "plots"
    hollow = "flagged" if "flagged" in table.columns else None  # ticket 06: QC-flagged tumours drawn hollow
    w = [plots.group_plot(table, facs, ro.column, ro.label, pd_ / f"{ro.name}_absolute", ci=ci if placed else None,
                          hollow=hollow,
                          suptitle=f"{ro.label} ({ro.role}): every tumour one dot; line + bar = model geometric "
                                   "mean, 95 % CI")]
    if focus is not None and f"{ro.name}_pct" in tum.columns:
        w.append(plots.group_plot(
            tum, facs, f"{ro.name}_pct", f"{ro.label}, % of {focus.reference}", pd_ / f"{ro.name}_normalised",
            ci=norm_ci if placed else None, marks=marks, ref_line=100.0, hollow=hollow,
            suptitle=f"% of model-estimated {focus.reference} in the same cell of the other factors; bars = "
                  f"fold-change 95 % CI" + (f"; stars: {fam.id} (Holm)" if fam and ro.role == "confirmatory"
                                             and fam.role == "confirmatory" else "")))
    w.append(plots.forest_plot(res, ro.label, pd_ / f"{ro.name}_forest"))
    w.append(plots.interaction_plot(means, facs, ro.label, pd_ / f"{ro.name}_interaction"))
    return [str(p) for p in w]


BATCH_TEXT = {
    "culture_batch": "one batch = one culture preparation (the date at the start of the file names); slices were "
                     "fixed after incubation, so the imaging day is recorded but not a batch",
    "session_folder": "one batch = one session folder",
    "acquisition_day": "one batch = one acquisition day from the OME metadata",
}


def _about(settings: dict, h: str, spec: dict, fams: list[Family], notes: list[str], table: pd.DataFrame,
           facs: list[Factor], blocked: stats.CellModel | None = None) -> list[str]:
    lines = [f"Tumour statistics — analysis plan lock {h}", f"Built {datetime.now():%Y-%m-%d %H:%M}",
             f"Model: {MODEL}; one tumour = one n; effects are fold-changes (ratios of geometric means) with 95 % CI.",
             "Holm within each family and readout. Stars only for confirmatory families of confirmatory readouts; "
             "exploratory results carry p-values but no stars; descriptive families are never tested.",
             f"Batch definition: {batch(settings)} ([tumour] batch; {BATCH_TEXT[batch(settings)]}). Batches that "
             "share a group are connected and enter the model as a fixed block effect (contrasts inside a connected "
             "set are within-batch estimates); contrasts whose compared groups come from different connected sets "
             f"(not a ratio of within-set ratios) are labelled '{stats.batch_of(table)[1]}' and never tested.",
             "Factors: " + "; ".join(f"{f['name']} ({', '.join(f['levels'])}; reference {f['reference']})"
                                     for f in spec["factors"]),
             "Readouts: " + "; ".join(f"{r['name']} ({r['role']})" for r in spec["readouts"]),
             f"Normalised values: % of the model-estimated {spec['normalise_by']} reference in the same cell of the "
             "other factors."]
    ts = tumour_settings(settings)
    thick, zrule = str(ts.get("section_thickness_um", "")).strip(), str(ts.get("z_range_rule", "")).strip()
    lines += [f"Mid-section: {float(ts['mid_section_um']):g} µm of depth for every tumour, centred on the layer of "
              "largest main-mass cross-section; the outermost layers count with the fraction of their extent inside "
              "it, so the thickness is the same at every z-step (achieved thickness per tumour: tumours.csv "
              "mid_section_um).",
              f"Section thickness: {thick + ' µm' if thick else 'not recorded'}; z-range at the microscope: "
              f"{zrule or 'not recorded'} ([tumour] section_thickness_um, z_range_rule; experiment metadata)."]
    lines += [f"Family {f.id} [{f.role}]: {f.effect} effect of {' × '.join(f.factors)} — {f.label}" for f in fams]
    present = {f.name: sorted(set(table[f.name].astype(str)) & set(f.levels), key=f.levels.index) for f in facs}
    absent = [f"{f.name}: {', '.join(lv for lv in f.levels if lv not in present[f.name])}" for f in facs
              if len(present[f.name]) < len(f.levels)]
    if absent:
        lines.append("Levels without data (contrasts needing them are listed as 'needs … data'): " + "; ".join(absent))
    if "imaging_session" in table.columns:
        lines.append("Session folders: " + ", ".join(sorted(set(table["imaging_session"].astype(str)))))
    col = stats.batch_of(table)[0]
    if col in table.columns:
        lines.append("Batches: " + ", ".join(sorted(set(table[col].astype(str)))))
    if blocked is not None:  # ticket 13: which batches were connected, block effects of the first such readout
        multi = [bs for bs in blocked.batch_sets if len(bs) > 1]
        lines.append("Block effect: connected batches " + "; ".join(" + ".join(bs) for bs in multi)
                     + " (estimated per readout; sensitivity models and residual plots use the batch-adjusted values)")
    lines += [f"NOTE: {n}" for n in notes]
    return lines


def _excel(out: Path, lines: list[str], tum: pd.DataFrame, stats_rows: list[dict], means: list[pd.DataFrame],
           fams: list[Family]) -> None:
    from ..manifest import STATS_LONG_COLUMNS

    with pd.ExcelWriter(out / "results.xlsx", engine="openpyxl") as xw:
        pd.DataFrame({"": lines}).to_excel(xw, sheet_name="About", index=False, header=False)
        tum.to_excel(xw, sheet_name="Tumours", index=False)
        pd.DataFrame(stats_rows, columns=STATS_LONG_COLUMNS).to_excel(xw, sheet_name="Statistics", index=False)
        (pd.concat(means, ignore_index=True) if means else pd.DataFrame()).to_excel(
            xw, sheet_name="Estimated means", index=False)
        pd.DataFrame([dict(id=f.id, effect=f.effect, factors=" × ".join(f.factors), role=f.role, label=f.label)
                      for f in fams]).to_excel(xw, sheet_name="Families", index=False)


def _fmt(r: pd.Series) -> str:
    if not np.isfinite(r["estimate"]):
        return "—"
    ci = f" [{r['ci_low']:.2f}, {r['ci_high']:.2f}]" if np.isfinite(r["ci_low"]) else ""
    return f"{r['estimate']:.2f}{ci}"


def _report(out: Path, lines: list[str], results: list[pd.DataFrame]) -> Path:
    md = ["# Tumour statistics report", "", *[f"- {x}" for x in lines], ""]
    for res in results:
        if not len(res):
            continue
        md += [f"## {res['readout'].iloc[0]} ({res['readout_role'].iloc[0]})", ""]
        for fid, g in res.groupby("family", sort=False):
            md += [f"### {fid} [{g['family_role'].iloc[0]}] {g['family_label'].iloc[0]}", "",
                   "| comparison | fold-change [95 % CI] | p | p (Holm) | | note |", "|---|---|---|---|---|---|"]
            for _, r in g.iterrows():
                p = f"{r['p']:.3g}" if np.isfinite(r["p"]) else ""
                pa = f"{r['p_adj']:.3g}" if np.isfinite(r["p_adj"]) else ""
                md.append(f"| {sensitivity.md_cell(r['comparison'])} | {_fmt(r)} | {p} | {pa} | {r['significance']} | "
                          f"{r['note']} |")
            md.append("")
    path = out / "stats_report.md"
    path.write_text("\n".join(md), encoding="utf-8")
    return path
