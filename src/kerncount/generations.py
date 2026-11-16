"""Write one generation folder (g1, g2, ...) with outputs per treatment, per n and for all n.

Every metric of the registry (metrics.py) goes through the same statistics. Besides the classical outputs
(Prism tables, % tables, per-n Dunnett, RM ANOVA, ratio t-tests, two-way ANOVA) each metric gets the hierarchical
model across all n (models.py) with its comparison families, equivalence test, leave-one-n-out and per-n effects,
and the generation gets acceptance.csv, stats_long.csv and manifest.json (the contracts for the browser UI,
docs/statistics.md).
"""
from __future__ import annotations

import json
import re
import shutil
import traceback
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd

from . import acceptance, analysis, design, manifest, modelplots, models
from . import metrics as metrics_mod

README = """Generation {name}
==============
Source: {source_text}
Built:  {built}

Layout
  counts_used.csv            the value used for every image in this generation (and where it came from)
  results.xlsx               everything below in one readable workbook
  stats_long.csv             every test result in one long table (filterable)
  manifest.json              list of every plot and table with a one-line description (used by the browser app)
  acceptance.csv             per-n acceptance criteria (pass / fail)
  <n>/<treatment>/           per-field values and plot for one treatment in one n
  <n>/                       one experiment: summary, Prism table, stats, plots per metric
                             (+ <metric>_fields_by_well.png, <metric>_position_effect.png)
  all_n/                     all experiments combined: summary, Prism tables (raw and % of control),
                             stats, plots (dot colour = n, see n_colours.csv), and the hierarchical model:
                             <metric>_model_means.csv / _model_absolute.png   estimates in absolute units
                             <metric>_contrasts.csv / _forest_<family>.png    ratios vs each control etc.
                             <metric>_variance_components.csv / .png          where the variation comes from
                             <metric>_diagnostics.png                         model checks
                             <metric>_loo*.csv, _per_n_ratio*.{{csv,png}}       robustness

Prism import: *_prism*.csv are Column tables (one column per treatment). In Prism use
File > Import, or open the CSV and paste into a Column table.

Statistics (details: docs/statistics.md). Every control group (reference + normalise_to) gets its own tests.
Within one n, fields are technical replicates (descriptive only): Dunnett vs each control.
Across n, each n is one biological replicate.
  Primary test ({primary_test}): {primary_text}
  Always reported as well: RM one-way ANOVA + Dunnett on the per-n means (Prism standard), ratio paired t-tests
  (Holm), two-way ANOVA (fields as replicates; optimistic), the hierarchical model with its comparison families,
  the vehicle equivalence test, leave-one-n-out and the effect in every n.
Absolute values (nuclei per field, intensity) and relative values (ratio / % of each control) are both given.
"""

PRIMARY_TEXT = {
    "mixed": "hierarchical (mixed) model on the log scale - fields in wells in n, a separate n-to-n variance per "
             "treatment class, Satterthwaite t; stars = adjusted p of the model ratio vs each control.",
    "rm_anova": "repeated-measures one-way ANOVA with n as the matched replicate, then Dunnett vs each control with "
                "the pooled error (raw scale).",
}


def _slug(s: str) -> str:
    return re.sub(r"[^\w.\-]+", "_", s).strip("_")


def _autosize(writer: pd.ExcelWriter):
    for ws in writer.book.worksheets:
        ws.freeze_panes = "A2"
        for col in ws.columns:
            width = max((len(str(c.value)) if c.value is not None else 0) for c in col)
            ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 8), 60)


def _fields_plot(sub: pd.DataFrame, out: Path, title: str, positions_ok: bool = True):
    """count per image of one treatment in one n; bar colour = well, position under each bar."""
    modelplots.treatment_fields(sub, out, title, positions_ok)


def _stars_footer(res, ref: str) -> str:
    split = "separate n-to-n variance per treatment class" if "class-specific" in res.model else "one common n-to-n variance"
    return (f"Stars: mixed model, ratio vs {ref} (log scale; fields in wells in n; {split}); p adjusted within the "
            "family (single-step multivariate t), Satterthwaite df. * p<0.05  ** p<0.01  *** p<0.001  **** p<0.0001; "
            "no mark = not significant.")


def _model_marks(res, ref: str) -> dict:
    c = res.contrasts
    fam = c[(c["family_kind"] == "vs_control") & (c["comparator"] == ref)]
    return {r.group: r.significance for r in fam.itertuples() if isinstance(r.significance, str)
            and r.significance not in ("", "ns")}


def _write_model_outputs(res, m: str, met, ad: Path, ncol: dict, cfg, main_ref: str | None, star_source: str):
    """Tables and plots of the hierarchical model for one metric in all_n/."""
    info = res.status_dict()
    pd.DataFrame([dict(key=k, value="; ".join(v) if isinstance(v, list) else v) for k, v in info.items()]).to_csv(
        ad / f"{m}_model_info.csv", index=False)
    if not res.varcomp.empty:
        res.varcomp.to_csv(ad / f"{m}_variance_components.csv", index=False)
    if not res.per_n.empty:
        res.per_n.to_csv(ad / f"{m}_per_n_ratios.csv", index=False)
    if not res.loo.empty:
        res.loo.to_csv(ad / f"{m}_loo.csv", index=False)
    if not res.loo_summary.empty:
        res.loo_summary.to_csv(ad / f"{m}_loo_summary.csv", index=False)
    label = met.label
    order = res.order or (res.prepared.treatments if res.prepared else [])
    if res.ok:
        e = res.emm.copy()
        for ref in cfg.controls:
            if ref in set(e["treatment"]):
                base = float(e.loc[e["treatment"] == ref, "estimate"].iloc[0])
                e[f"pct_of_{ref}"] = 100 * e["estimate"] / base if base else np.nan
        e.to_csv(ad / f"{m}_model_means.csv", index=False)
        c = res.contrasts.copy()
        c.insert(c.columns.get_loc("ci_high") + 1, "pct_of_comparator", 100 * c["estimate"])
        c.to_csv(ad / f"{m}_contrasts.csv", index=False)
        sig = _model_marks(res, main_ref) if (main_ref and star_source == "mixed") else None
        if star_source != "mixed" and main_ref:
            sig = None  # stars of the RM ANOVA are on the bar plots; model plot shows none to avoid mixing tests
        foot = ((_stars_footer(res, main_ref) if sig is not None else
                 "No stars here: the primary test is the RM ANOVA (stars on the bar plots). ") +
                f" Response {res.scale.describe}; zero fields are drawn at the bottom edge.")
        modelplots.model_absolute(res, ad / f"{m}_model_absolute", ncol, sig or {}, main_ref, foot,
                                  f"All n: {label} (mixed model)")
        loo = res.loo if (res.primary_test == "mixed") else None
        for fid, fam in res.contrasts.groupby("family_id", sort=False):
            fam = fam.reset_index(drop=True)
            kind = fam["family_kind"].iloc[0]
            if kind == "vehicle_check":
                foot = (f"Equivalence (TOST): the vehicle counts as equivalent if the whole 90% CI lies within the "
                        f"shaded margin ({1 / (1 + cfg.margin):.2f}-{1 + cfg.margin:.2f}). 'Not significant' alone "
                        "does not show that there is no effect.")
                modelplots.forest(fam, None, ad / f"{m}_forest_{fid}", ncol,
                                  f"All n: {label}, vehicle check (equivalence, 90% CI)", foot,
                                  res.scale.estimate_type, margin=cfg.margin)
                continue
            adj = "simultaneous 95% CI (dashed) = adjusted for the family"
            if star_source == "mixed":
                st_txt = (f"Stars: this model comparison, p adjusted within the family (single-step multivariate t, "
                          f"Satterthwaite df). {adj}.")
            else:
                st_txt = (f"Stars shown are the model's; the primary test (RM ANOVA + Dunnett) is on the bar plots. {adj}.")
            foot = (f"{st_txt} Coloured dots: estimate with that n left out (leave-one-n-out). Model: {res.model}. "
                    "* p<0.05  ** p<0.01  *** p<0.001  **** p<0.0001; ns = not significant.")
            modelplots.forest(fam, loo, ad / f"{m}_forest_{fid}", ncol, f"All n: {label}, {fam['family'].iloc[0]}",
                              foot, res.scale.estimate_type)
        modelplots.variance_components(
            res.varcomp, ad / f"{m}_variance_components",
            f"All n: {label}, where the variation comes from (mixed model)",
            "Variances on the log scale (they add up). SD on the log scale ≈ CV for small values. Lines = approximate "
            "95% CI where available (rough with few n's; the n-level variances rest on very few levels).")
        try:
            modelplots.diagnostics(res, ad / f"{m}_diagnostics", ncol, f"All n: {label}, model diagnostics")
        except Exception:  # noqa: BLE001 - a diagnostic plot must not stop the build
            pass
    if not res.per_n.empty:
        for ref in cfg.controls:
            fam = None
            if res.ok:
                c = res.contrasts
                fam = c[(c["family_kind"] == "vs_control") & (c["comparator"] == ref)]
            modelplots.per_n_ratio(
                res.per_n, fam, ref, order, ad / f"{m}_per_n_ratio_vs_{_slug(ref)}", ncol,
                f"All n: {label}, effect in each n (ratio to {ref})",
                "Dots: ratio of geometric means of the fields within each n (colour = n, joined per n). Black bar: "
                "model estimate across all n. Spread between the dots = how much the effect differs between experiments.")


def build(project: Path, gen: dict, image_table: pd.DataFrame, metrics: dict | None, settings: dict,
          design_table: pd.DataFrame | None = None) -> Path:
    """metrics: the metric registry (metrics.registry) or a legacy {column: label} dict; count is always included.
    design_table: n, treatment, image, well, position, timepoint from design.assign (computed here if None)."""
    out = project / "generations" / gen["name"]
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    a = settings["analysis"]
    reg = metrics_mod.normalise(metrics)
    extra_cols = [m.column for k, m in reg.items() if k != "count"]

    t = gen["table"].copy()
    t = t.merge(image_table, on=["n", "treatment", "image"], how="left", suffixes=("", "_img"))
    if "qc_level" not in t:
        t["qc_level"] = ""
    if a.get("exclude_low_qc"):
        drop = (t["qc_level"] == "low") & (t["source"] == "auto")
        t.loc[drop, "excluded"] = True
        t.loc[drop, "source"] = "auto (excluded: QC low)"
    if design_table is None:
        design_table = design.assign(t[["n", "treatment", "image"]], settings, project)
    dcols = ["well", "position", "timepoint"]
    t = t.drop(columns=[c for c in dcols if c in t]).merge(
        design_table[["n", "treatment", "image"] + dcols].astype({"image": int}).drop_duplicates(["n", "treatment", "image"]),
        on=["n", "treatment", "image"], how="left")
    used_cols = ["n", "treatment", "image", "count", "source", "excluded", "qc_level"] + dcols
    t[used_cols].to_csv(out / "counts_used.csv", index=False)
    try:
        positions_ok = design.positions_meaningful(settings)
    except ValueError:
        positions_ok = False

    # acceptance per n (on the counts); enforce = leave failing n's out of every all-n statistic
    acc = acceptance.evaluate(t[used_cols], settings)
    acc.to_csv(out / "acceptance.csv", index=False)
    acc_sum = acceptance.summary(acc, settings)

    metrics_live = {k: m for k, m in reg.items() if m.column in t and pd.to_numeric(t[m.column], errors="coerce").notna().any()}
    live = t[~t["excluded"].astype(bool)]
    live_all = live[~live["n"].isin(acc_sum["excluded"])]
    per_n_summ, per_n_stats, all_summ, all_stats, results = {}, {}, {}, {}, {}
    long_rows: list[dict] = []
    n_by_slug, tr_by_slug = {}, {}
    ns_all = sorted(t["n"].unique())
    ncol = modelplots.n_colours(ns_all)

    for n in ns_all:
        nd = out / _slug(n)
        nd.mkdir()
        n_by_slug[_slug(n)] = n
        for tr in analysis.treatment_order(t.loc[t["n"] == n, "treatment"], a["treatment_order"]):
            td = nd / _slug(tr)
            td.mkdir()
            tr_by_slug[(_slug(n), _slug(tr))] = tr
            sub = t[(t["n"] == n) & (t["treatment"] == tr)].sort_values("image")
            keep = ["image", "count", "source", "excluded", "qc_level"] + dcols[:2] + [c for c in extra_cols if c in sub]
            sub[[c for c in keep if c in sub]].to_csv(td / "fields.csv", index=False)
            _fields_plot(sub, td / "count_fields.png", f"{n} / {tr}", positions_ok)
        for m, met in metrics_live.items():
            vals = live.assign(value=pd.to_numeric(live[met.column], errors="coerce"))
            if vals.loc[vals["n"] == n, "value"].notna().sum() == 0:
                continue
            st = analysis.per_n(vals, m, met.label, n, nd, settings)
            per_n_stats.setdefault(m, []).append(st)
            summ = pd.read_csv(nd / f"{m}_summary.csv")
            per_n_summ.setdefault(m, []).append(summ)
            long_rows += manifest.rows_classical(st, "n", m, n, summ.set_index("treatment")["mean"])
            if (nd / f"{m}_twoway_anova.csv").exists():
                long_rows += manifest.rows_twoway(pd.read_csv(nd / f"{m}_twoway_pairs.csv"),
                                                  pd.read_csv(nd / f"{m}_twoway_anova.csv"), "n", m, n)
            sub = vals[vals["n"] == n]
            order_n = analysis.treatment_order(sub["treatment"], a["treatment_order"])
            modelplots.fields_by_well(sub, order_n, met.label, nd / f"{m}_fields_by_well",
                                      f"{n}: {met.label}, every field", positions_ok)
            if positions_ok:
                modelplots.position_effect(sub.assign(well_l=sub["well"]), nd / f"{m}_position_effect",
                                           f"{n}: {met.label} by field position (per well)", "well",
                                           {w: modelplots.WELL_COLOURS[i % 6] for i, w in enumerate(sorted(sub["well"].unique()))},
                                           lambda w: f"well {int(w)}")

    ad = out / "all_n"
    ad.mkdir()
    cfg = models.config(settings, live_all["treatment"].unique())
    model_status, fallbacks = {}, []
    for m, met in metrics_live.items():
        vals = live_all.assign(value=pd.to_numeric(live_all[met.column], errors="coerce"))
        vals = vals[vals["value"].notna()]
        if vals.empty:
            continue
        res = models.analyse(vals, met, settings, positions_ok)
        results[m] = res
        means = vals.groupby("treatment")["value"].mean()
        main_ref = cfg.primary_ref if cfg.primary_ref in set(vals["treatment"]) else \
            next((c for c in cfg.controls if c in set(vals["treatment"])), None)
        use_mixed = cfg.primary_test == "mixed" and res.ok
        stars_by_ref, star_note = None, ""
        if use_mixed:
            stars_by_ref = {ref: (_model_marks(res, ref), _stars_footer(res, ref)) for ref in cfg.controls
                            if ref in set(vals["treatment"])}
        elif cfg.primary_test == "mixed":
            star_note = f"Mixed model not used for this metric ({res.model}). "
            fallbacks.append(m)
        star_source = "mixed" if use_mixed else "rm_anova"
        res.primary_test = star_source
        summ, st = analysis.all_n(vals, m, met.label, ad, settings, stars_by_ref=stars_by_ref, main_ref=main_ref,
                                  star_note=star_note)
        all_summ[m], all_stats[m] = summ, st
        long_rows += manifest.rows_classical(st, "all_n", m, None, means)
        if (ad / f"{m}_twoway_anova.csv").exists():
            long_rows += manifest.rows_twoway(pd.read_csv(ad / f"{m}_twoway_pairs.csv"),
                                              pd.read_csv(ad / f"{m}_twoway_anova.csv"), "all_n", m)
        try:
            _write_model_outputs(res, m, met, ad, ncol, cfg, main_ref, star_source)
        except Exception as e:  # noqa: BLE001 - never crash a generation build on a plot
            res.notes.append(f"model outputs incomplete ({type(e).__name__}: {e})")
            traceback.print_exc()
        long_rows += manifest.rows_model(res, m)
        model_status[m] = res.status_dict() | dict(stars_from=star_source)
        if positions_ok and m == "count":
            modelplots.position_effect(vals, ad / f"{m}_position_effect", f"All n: {met.label} by field position (per n)",
                                       "n", ncol, manifest.short)

    stats_long = manifest.write_stats_long(out, long_rows)
    _excel(out, gen, t, used_cols, a, cfg, reg, metrics_live, all_summ, all_stats, per_n_summ, per_n_stats, results,
           acc, acc_sum, model_status, ad)

    src = "fully automatic (auto_count)" if gen["source"] == "auto" else \
        f"previous generation + manual values from column {gen['source']} of counts.csv ({gen['changes']} cells)"
    (out / "README.txt").write_text(README.format(name=gen["name"], source_text=src, built=f"{datetime.now():%Y-%m-%d %H:%M}",
                                                  primary_test=cfg.primary_test, primary_text=PRIMARY_TEXT[cfg.primary_test]),
                                    encoding="utf-8")
    (out / "generation.json").write_text(json.dumps(dict(name=gen["name"], source=gen["source"], fingerprint=gen["fingerprint"],
                                                         changes=gen["changes"], built=datetime.now().isoformat(timespec="seconds")),
                                                    indent=2), encoding="utf-8")
    notes = []
    if acc_sum["excluded"]:
        notes.append(f"acceptance enforced: {', '.join(acc_sum['excluded'])} left out of all-n statistics")
    if fallbacks:
        notes.append(f"mixed model not available for {', '.join(fallbacks)}: RM ANOVA drives the stars there")
    primary = dict(test=cfg.primary_test, family=f"vs {cfg.primary_ref}", reference=cfg.primary_ref,
                   fallback_metrics=fallbacks, controls=cfg.controls,
                   classes={c: [t_ for t_, k in cfg.classes.items() if k == c] for c in cfg.class_order},
                   positive_control=cfg.positive_control or None,
                   vehicle_check=[f"{v} vs {c}" for v, c in cfg.vehicle_pairs],
                   equivalence_margin=cfg.margin)
    manifest.build_manifest(out, gen, metrics_live, n_by_slug, tr_by_slug, primary, model_status, acc_sum, notes)
    return out


def _excel(out, gen, t, used_cols, a, cfg, reg, metrics_live, all_summ, all_stats, per_n_summ, per_n_stats, results,
           acc, acc_sum, model_status, ad):
    other = [m for m in metrics_live if m != "count" and m in all_summ]
    with pd.ExcelWriter(out / "results.xlsx", engine="openpyxl") as xw:
        lines = [f"Generation {gen['name']}", f"Source: {gen['source']}",
                 f"Manual changes in this step: {gen['changes']}",
                 f"Built {datetime.now():%Y-%m-%d %H:%M}",
                 f"Control groups (main first): {', '.join(analysis.control_groups(a))}",
                 f"Primary test (drives the stars): {cfg.primary_test}; main claim: vs {cfg.primary_ref}",
                 "Within-n stats treat fields as technical replicates (descriptive); Dunnett vs each control.",
                 "All-n stats: each n = one biological replicate.",
                 "  Mixed model (sheets 'Model - ...'): log scale (results as ratios), fields in wells in n, n-to-n "
                 "variability estimated separately per treatment class; Satterthwaite t; p adjusted within each family "
                 "(single-step multivariate t). Means are given in absolute units (geometric means) and as % of each control.",
                 "  RM one-way ANOVA + Dunnett vs each control (Prism standard, raw scale); ratio paired t-tests (Holm) as a check.",
                 "  Vehicle check: equivalence test (TOST, 90% CI within the margin) - 'not significant' is not 'no effect'.",
                 "  Robustness: leave-one-n-out (is the conclusion the same whichever n is dropped?) and the ratio in every n.",
                 "Two-way ANOVA (extra): per n treatment x well, all n treatment x n; fields = replicates; Tukey all pairs.",
                 f"Treatment classes: " + "; ".join(f"{c}: {', '.join(x for x, k in cfg.classes.items() if k == c)}"
                                                   for c in cfg.class_order)]
        for m, s in model_status.items():
            lines.append(f"Model for {m}: {s['model']} (stars from {s['stars_from']})" +
                         (f"; notes: {'; '.join(s['notes'])}" if s["notes"] else ""))
        fails = acceptance.failing(acc)
        lines.append("Acceptance: " + ("all n pass" if not fails else f"failing: {', '.join(fails)}") +
                     (" (enforced: left out of all-n statistics)" if acc_sum["enforce"] and fails else
                      " (report only)" if fails else ""))
        pd.DataFrame({"": lines}).to_excel(xw, sheet_name="About", index=False, header=False)
        t[used_cols].to_excel(xw, sheet_name="Counts used", index=False)
        if "count" in all_summ:
            all_summ["count"].to_excel(xw, sheet_name="Count - all n", index=False)
            all_stats["count"].to_excel(xw, sheet_name="Count - stats all n", index=False)
        if "count" in per_n_summ:
            pd.concat(per_n_summ["count"]).to_excel(xw, sheet_name="Count - per n", index=False)
            pd.concat(per_n_stats["count"]).to_excel(xw, sheet_name="Count - stats per n", index=False)
        inten = [m for m in other if m.startswith("intensity_")]
        rest = [m for m in other if not m.startswith("intensity_")]
        for name, ms in (("Intensity", inten), ("Other metrics", rest)):
            if not ms:
                continue
            pd.concat([all_summ[m].assign(metric=m) for m in ms]).to_excel(xw, sheet_name=f"{name} - all n", index=False)
            pd.concat([all_stats[m].assign(metric=m) for m in ms]).to_excel(xw, sheet_name=f"{name} - stats all n", index=False)
            pn = [pd.concat(per_n_summ[m]).assign(metric=m) for m in ms if m in per_n_summ]
            if pn:
                pd.concat(pn).to_excel(xw, sheet_name=f"{name} - per n", index=False)
        for level, files in (("all n", sorted(ad.glob("*_twoway_anova.csv"))),
                             ("per n", sorted(out.glob("*/*_twoway_anova.csv")))):
            files = [f for f in files if f.parent != ad] if level == "per n" else files
            if not files:
                continue
            for kind in ("anova", "pairs"):
                tabs = []
                for f in files:
                    m = f.name[: -len("_twoway_anova.csv")]
                    g = f.with_name(f"{m}_twoway_{kind}.csv")
                    if g.exists():
                        tab = pd.read_csv(g)
                        tab.insert(0, "metric", m)
                        if level == "per n":
                            tab.insert(0, "n", f.parent.name)
                        tabs.append(tab)
                if tabs:
                    pd.concat(tabs).to_excel(xw, sheet_name=f"Two-way {kind} - {level}", index=False)

        def cat(attr, filt=None):
            tabs = []
            for m, r in results.items():
                x = getattr(r, attr)
                if x is None or len(x) == 0:
                    continue
                x = x if filt is None else filt(x)
                tabs.append(x.assign(metric=m, unit=reg[m].unit if m in reg else ""))
            return pd.concat(tabs) if tabs else None

        sheets = [
            ("Model - means (absolute)", cat("emm")),
            ("Model - ratios vs controls", cat("contrasts", lambda c: c[c["family_kind"] == "vs_control"].assign(
                pct_of_comparator=100 * c["estimate"]))),
            ("Comparison families", cat("contrasts", lambda c: c.assign(pct_of_comparator=100 * c["estimate"]))),
            ("Variance components", cat("varcomp")),
            ("Robustness - stability", cat("loo_summary")),
            ("Robustness - leave-one-out", cat("loo")),
            ("Robustness - effect per n", cat("per_n")),
        ]
        for name, tab in sheets:
            if tab is not None and len(tab):
                tab.to_excel(xw, sheet_name=name, index=False)
        ms = pd.DataFrame([dict(metric=m, **{k: ("; ".join(v) if isinstance(v, list) else v) for k, v in s.items()})
                           for m, s in model_status.items()])
        if len(ms):
            ms.to_excel(xw, sheet_name="Model status", index=False)
        acc.to_excel(xw, sheet_name="Acceptance", index=False)
        _autosize(xw)
