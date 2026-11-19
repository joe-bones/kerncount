"""Write one generation folder (g1, g2, ...) with outputs per treatment, per n and for all n.

Every metric of the registry (metrics.py) goes through the same statistics. Besides the classical outputs
(Prism tables, % tables, per-n Dunnett, RM ANOVA, ratio t-tests, two-way ANOVA) each metric gets the hierarchical
model across all n (models.py) with its comparison families, equivalence test, leave-one-n-out and per-n effects;
fraction metrics also get a density-adjusted model. The generation gets acceptance.csv, stats_long.csv and
manifest.json (the contracts for the browser UI, docs/statistics.md).

Rules applied here:
- the negative-control condition ([experiment] negative_control) is left out of every test, plot and summary;
- cell lines are never pooled: with more than one cell line the all-n analysis runs once per line, in
  all_n/<cell line>/ (one line: all_n/ as before);
- with one n the all-n level keeps every file, column and plot, filled with that n's estimates and "needs ≥ 2 n"
  instead of across-n CIs, p-values and stars.
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

from . import acceptance, analysis, design, fileio, manifest, modelplots, models
from . import duplicates as duplicates_mod
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
                             (+ <metric>_fields_by_well.png, <metric>_position_effect.png, <fraction>_density.png)
  all_n/                     all experiments combined{line_text}: summary, Prism tables (raw and % of control),
                             stats, plots (dot colour = n, see n_colours.csv), and the hierarchical model:
                             <metric>_model_means.csv / _model_absolute.png   estimates in absolute units
                             <metric>_contrasts.csv / _forest_<family>.png    ratios (fractions: % points) vs controls
                             <metric>_variance_components.csv / .png          where the variation comes from
                             <metric>_diagnostics.png                         model checks
                             <metric>_loo*.csv, _per_n_ratio*.{{csv,png}}       robustness
                             <fraction>_density_adjusted_*.csv, _density.png  adjusted for local crowding

Prism import: *_prism*.csv are Column tables (one column per treatment). In Prism use
File > Import, or open the CSV and paste into a Column table.

Statistics (details: docs/statistics.md). Every control group (reference + normalise_to) gets its own tests.
Within one n, fields are technical replicates (descriptive only): Dunnett vs each control.
Across n, each n is one biological replicate.
  Model unit: {unit_text}
  Primary test ({primary_test}): {primary_text}
  Always reported as well: RM one-way ANOVA + Dunnett on the per-n means (Prism standard), ratio paired t-tests
  (Holm), two-way ANOVA (fields as replicates; optimistic), the hierarchical model with its comparison families,
  the vehicle equivalence test, leave-one-n-out and the effect in every n.
Absolute values (nuclei per field, % of nuclei, intensity) and relative values (ratio / % of each control,
percentage-point differences for fractions) are both given.
{extra_notes}"""

PRIMARY_TEXT = {
    "mixed": "hierarchical (mixed) model - wells in n, a separate n-to-n variance per treatment class (default: one "
             "per non-control treatment), Satterthwaite t; stars = adjusted p of the model comparison vs each control.",
    "rm_anova": "repeated-measures one-way ANOVA with n as the matched replicate, then Dunnett vs each control with "
                "the pooled error (raw scale).",
}
UNIT_TEXT = {
    "well": "well - each well's mean per field (pooled % for fractions); effects are ratios of mean values "
            "(e.g. 'ratio of mean nuclei per field') or differences in percentage points for fractions.",
    "field": "field - every field is one value; effects are ratios of geometric means of fields.",
}
STARS_LEGEND = "* p<0.05  ** p<0.01  *** p<0.001  **** p<0.0001; no mark = not significant."


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


def _class_text(res) -> str:
    if "class-specific" not in res.model:
        return "one common n-to-n variance"
    return "a separate n-to-n variance per treatment class"


def _stars_footer(res, ref: str) -> str:
    est = res.estimand
    what = "difference in percentage points (test on the logit scale)" if res.scale.kind == "fraction" else est
    return (f"Estimand: {est}. Stars: mixed model, {what} vs {ref} (unit: {res.unit}; {_class_text(res)}); "
            f"p adjusted within the family (single-step multivariate t), Satterthwaite df. {STARS_LEGEND}")


def _one_n_footer(res) -> str:
    return (f"Estimand: {res.estimand}. Only one n: estimates come from that n (N = 1); across-n CIs, p-values and "
            f"stars: {models.NEEDS_2N}. Field and well variation is not used as a biological test.")


def _model_marks(res, ref: str) -> dict:
    c = res.contrasts
    fam = c[(c["family_kind"] == "vs_control") & (c["comparator"] == ref)]
    return {r.group: r.significance for r in fam.itertuples() if isinstance(r.significance, str)
            and r.significance not in ("", "ns")}


def _write_model_outputs(res, m: str, met, ad: Path, ncol: dict, cfg, main_ref: str | None, star_source: str,
                         where: str = "All n"):
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
    frac = res.scale is not None and res.scale.kind == "fraction"
    raw = res.scale is not None and res.scale.kind == "difference"  # well means <= 0: raw-scale difference (I3)
    unit = met.unit or ""
    raw_txt = (models.RAW_SCALE_NOTE + ". ") if raw else ""
    order = res.order or (res.prepared.treatments if res.prepared else [])
    one_n = res.descriptive
    if res.has_estimates:
        e = res.emm.copy()
        for ref in cfg.controls:
            if ref in set(e["treatment"]):
                base = float(e.loc[e["treatment"] == ref, "estimate"].iloc[0])
                e[f"pct_of_{ref}"] = 100 * e["estimate"] / base if base > 0 else np.nan
        e.to_csv(ad / f"{m}_model_means.csv", index=False)
        res.contrasts.to_csv(ad / f"{m}_contrasts.csv", index=False)
        sig = _model_marks(res, main_ref) if (main_ref and star_source == "mixed" and res.ok) else None
        if one_n:
            foot = _one_n_footer(res)
        elif sig is not None:
            foot = _stars_footer(res, main_ref)
        else:
            foot = f"Estimand: {res.estimand}. No stars here: the primary test is the RM ANOVA (stars on the bar plots)."
        foot = raw_txt + foot + f" Response {res.scale.describe}" + (
            "" if (frac or raw) else "; zero fields are drawn at the bottom edge") + "."
        modelplots.model_absolute(res, ad / f"{m}_model_absolute", ncol, sig or {}, main_ref, foot,
                                  f"{where}: {label} (mixed model)" if res.ok else f"{where}: {label} (one n)")
        loo = res.loo if (res.primary_test == "mixed" and res.ok) else None
        for fid, fam in res.contrasts.groupby("family_id", sort=False):
            fam = fam.reset_index(drop=True)
            kind = fam["family_kind"].iloc[0]
            if kind == "vehicle_check":
                what = ("difference in percentage points" if frac else "raw-scale difference" if raw else "ratio")
                mtxt = (f"±{cfg.margin_pp:g} percentage points" if frac else
                        f"±{cfg.margin:.0%} of the comparator mean" if raw else
                        f"{1 / (1 + cfg.margin):.2f}-{1 + cfg.margin:.2f}")
                ci = str(fam["test"].iloc[0]).rsplit("; ", 1)[-1][:-1] if len(fam) else "90% CI"  # drop the closing ")"
                foot = (raw_txt + f"Equivalence (TOST, Holm across pairs): the vehicle counts as equivalent if the "
                        f"Holm-adjusted TOST p < 0.05, i.e. the whole {ci} of the {what} lies within the shaded "
                        f"margin ({mtxt}); not equivalent if it lies wholly outside. 'Not significant' alone does "
                        "not show that there is no effect." + (f" Only one n: {models.NEEDS_2N}." if one_n else ""))
                modelplots.forest(fam, None, ad / f"{m}_forest_{fid}", ncol,
                                  f"{where}: {label}, vehicle check (equivalence)", foot,
                                  res.scale.estimate_type, margin=cfg.margin, pp=frac, unit=unit)
                continue
            adj = "simultaneous 95% CI (dashed) = adjusted for the family"
            if one_n:
                st_txt = (f"Only one n: estimates from that n (N = 1); CIs, p-values and stars: {models.NEEDS_2N}.")
            elif star_source == "mixed":
                st_txt = (f"Stars: this model comparison, p adjusted within the family (single-step multivariate t, "
                          f"Satterthwaite df{'; test on the logit scale' if frac else ''}). {adj}.")
            else:
                st_txt = f"Stars shown are the model's; the primary test (RM ANOVA + Dunnett) is on the bar plots. {adj}."
            extra = (" Odds ratios and % of the comparator are in the contrasts table; a ratio of near-zero fractions "
                     "is not meaningful." if frac else "")
            loo_txt = "" if one_n else " Coloured dots: estimate with that n left out (leave-one-n-out)."
            foot = (raw_txt + f"Estimand: {res.estimand}. {st_txt}{loo_txt} Model: {res.model}.{extra} "
                    + ("" if one_n else STARS_LEGEND))
            modelplots.forest(fam, loo, ad / f"{m}_forest_{fid}", ncol, f"{where}: {label}, {fam['family'].iloc[0]}"
                              + (" (% points)" if frac else " (DIFFERENCE, raw scale)" if raw else ""), foot,
                              res.scale.estimate_type, pp=frac, unit=unit)
        modelplots.variance_components(
            res.varcomp, ad / f"{m}_variance_components",
            f"{where}: {label}, where the variation comes from" + (" (mixed model)" if res.ok else " (one n)"),
            (raw_txt + ("Variances on the raw scale (squared intensity units; they add up). " if raw else
                        "Variances on the analysis scale (log, or logit for fractions; they add up). SD on the log "
                        "scale ≈ CV for small values. ") + "Lines = approximate 95% CI where available (rough with few n's; the n-level variances "
             "rest on very few levels)." if res.ok else
             f"Only one n: the n-level variances: {models.NEEDS_2N}; the well-to-well variance shown is descriptive "
             "(technical replicates)."))
        try:
            modelplots.diagnostics(res, ad / f"{m}_diagnostics", ncol, f"{where}: {label}, model diagnostics")
        except Exception:  # noqa: BLE001 - a diagnostic plot must not stop the build
            traceback.print_exc()
    if not res.per_n.empty:
        for ref in cfg.controls:
            fam = None
            if res.has_estimates:
                c = res.contrasts
                fam = c[(c["family_kind"] == "vs_control") & (c["comparator"] == ref)]
            if raw:
                col, txt = "difference", (raw_txt + "Dots: difference of the mean values within each n (colour = n, "
                                          "joined per n; thin lines = technical 95% interval from the fields, "
                                          "descriptive).")
            elif frac:
                col, txt = "pp_difference", ("Dots: difference of the pooled fractions within each n (colour = n, "
                                             "joined per n; thin lines = technical 95% interval from the fields, "
                                             "descriptive).")
            elif res.unit == "well":
                col, txt = "ratio_of_means", ("Dots: ratio of mean values within each n (colour = n, joined per n; "
                                              "thin lines = technical 95% interval from the fields, descriptive).")
            else:
                col, txt = "ratio", "Dots: ratio of geometric means of the fields within each n (colour = n, joined per n)."
            bar = ("Black bar: model estimate across all n." if res.ok else
                   f"Black bar: estimate from the one n ({models.NEEDS_2N} for an across-n test).")
            modelplots.per_n_ratio(
                res.per_n, fam, ref, order, ad / f"{m}_per_n_ratio_vs_{_slug(ref)}", ncol,
                f"{where}: {label}, effect in each n " + (f"(difference from {ref})" if (frac or raw)
                                                         else f"(ratio to {ref})"),
                f"{txt} {bar} Spread between the dots = how much the effect differs between experiments.",
                col=col, pp=frac or raw)
    if res.density is not None:
        _write_density(res, m, met, ad, cfg, where)


def _write_density(res, m: str, met, ad: Path, cfg, where: str):
    dens = res.density
    em = dens.emm.copy()
    raw = res.emm.set_index("treatment") if len(res.emm) else pd.DataFrame()
    if len(raw):
        em.insert(em.columns.get_loc("percent") + 3, "raw_percent", em["treatment"].map(raw["percent"]))
    em.insert(1, "covariate", dens.covariate)
    em.to_csv(ad / f"{m}_density_adjusted_means.csv", index=False)
    ct = dens.contrasts.copy()
    if len(res.contrasts) and len(ct):
        rc = res.contrasts.assign(family=res.contrasts["family"] + " (density-adjusted)").set_index(["family", "comparison"])
        key = list(zip(ct["family"], ct["comparison"]))
        ct.insert(ct.columns.get_loc("pp_difference") + 1, "raw_pp_difference",
                  [rc["pp_difference"].get(k, np.nan) for k in key])
    ct.insert(0, "covariate", dens.covariate)
    for k_, v_ in dens.slope.items():
        if k_ != "covariate":
            ct[f"slope_{k_}"] = v_
    ct.to_csv(ad / f"{m}_density_adjusted_contrasts.csv", index=False)
    pts, curves = models.density_curves(dens)
    cov_label = metrics_mod.describe_column(dens.covariate)
    xl = f"{cov_label.label} ({cov_label.unit}; well mean)" if cov_label else f"{dens.covariate} (well mean)"
    flagged = sorted({r.comparison for r in ct.itertuples() if isinstance(r.density_flag, str) and r.density_flag})
    flag_txt = (f" Flag: {models.DENSITY_FLAG} for {', '.join(flagged[:6])}{' ...' if len(flagged) > 6 else ''}"
                f" (overlap < {models.OVERLAP_MIN:.0%}: each condition needs at least half of its wells inside the "
                "other's density range)." if flagged else " All compared conditions overlap in density.")
    foot = (f"Dots = wells (marker = n), lines = density-adjusted model per treatment ({dens.model}). The raw "
            f"estimate stays the headline; the adjusted one asks 'at the same crowding?'.{flag_txt}")
    modelplots.density_scatter(pts, curves, res.order or dens.prepared.treatments, ad / f"{m}_density",
                               f"{where}: {met.label} vs local crowding", foot, xl, f"{met.label} (% of nuclei)")


def _density_per_n(sub: pd.DataFrame, met, cov: str, order: list[str], nd: Path, n: str):
    got = models.density_fields(sub, met, cov)
    if got is None:
        return
    pts, curves = got
    cl = metrics_mod.describe_column(cov)
    xl = f"{cl.label} ({cl.unit})" if cl else cov
    modelplots.density_scatter(pts, curves, order, nd / f"{met.id}_density", f"{n}: {met.label} vs local crowding",
                               "Dots = fields; lines = descriptive fit within this experiment (empirical logit on "
                               "treatment + one common density slope). Technical replicates only: describes this n.",
                               xl, f"{met.label} (% of nuclei)", by_n=False)


def _write_layers(out: Path, blocks, results: dict, metrics_live: dict, settings: dict, order: list[str], ncol: dict,
                  n_by_slug: dict) -> dict:
    """Layer comparison (layers.py) for every channel with a readout ladder: all_n[/<line>]/ and <n>/.
    Returns the manifest entry {channel: {layers: [...], extras: [...], references: [...], files: {...}}}."""
    from . import layers as layers_mod

    chans = layers_mod.channels(metrics_live)
    info = {}
    if not chans:
        return info
    slug_of = {v: k for k, v in n_by_slug.items()}
    for ch in chans:
        lays = layers_mod.layer_list(ch, metrics_live)
        refs_all = []
        for line, ad, rows in blocks:
            key = (lambda m: m) if len(blocks) == 1 else (lambda m, line=line: f"{line}/{m}")
            res_b = {L["metric"]: results[key(L["metric"])] for L in lays if key(L["metric"]) in results}
            if not res_b:
                continue
            present = set(rows["treatment"])
            refs = [c for c in models.config(settings, present).controls if c in present]
            refs_all += [r for r in refs if r not in refs_all]
            where = f"All n ({line})" if len(blocks) > 1 else "All n"
            frac_m = f"readout_frac_positive_{ch}"
            dens = (f"{frac_m}_density.png in this folder and intensity/density_[n]__{ch}.png"
                    if frac_m in metrics_live else f"intensity/density_[n]__{ch}.png")
            try:
                layers_mod.write_level(ad, ch, lays, res_b, refs, order, "all_n", where, line, ncol=ncol,
                                       density_hint=dens)
                for n in sorted(rows["n"].unique()):
                    nd = out / slug_of.get(n, _slug(n))
                    refs_n = [r for r in refs if ((rows["n"] == n) & (rows["treatment"] == r)).any()]
                    layers_mod.write_level(nd, ch, lays, res_b, refs_n, order, "n", manifest.short(n), line,
                                           n_filter=n, ncol=ncol, density_hint=f"{frac_m}_density.png in this folder")
            except Exception:  # noqa: BLE001 - never crash a generation build on the layer comparison
                traceback.print_exc()
        info[ch] = dict(layers=[dict(index=L["index"], key=L["key"], label=L["label"], metric=L["metric"],
                                     step=L["step"]) for L in lays if not L["extra"]],
                        extras=[dict(key=L["key"], label=L["label"], metric=L["metric"]) for L in lays if L["extra"]],
                        references=refs_all,
                        files=dict(effects=f"layers_{ch}_effects.csv", absolute=f"layers_{ch}_absolute.csv",
                                   concordance=f"layers_{ch}_concordance.csv",
                                   ladder_plot=f"layers_{ch}_vs_[reference].png",
                                   absolute_plot=f"layers_{ch}_absolute.png"))
    return info


def coverage_notes(t: pd.DataFrame, metrics_live: dict) -> list[str]:
    """Plain-language notes on data that drop out of the statistics without an error (docs/design.md, "Unusual
    designs"): a well or a whole condition with no usable field left (QC low / manual exclusions), and a condition
    that has fields but no value of a readout in an n (e.g. no measurable nuclei for % positive)."""
    notes = []
    exc = t["excluded"].astype(bool)
    for (n, tr), g in t.groupby(["n", "treatment"], sort=True):
        ok = ~exc.loc[g.index]
        if not ok.any():
            notes.append(f"{n} / {tr}: no usable field (all {len(g)} excluded: {_why(g)}); left out of every statistic")
            continue
        if "well" in g:
            for w, gw in g.groupby("well"):
                if not (~exc.loc[gw.index]).any():
                    notes.append(f"{n} / {tr}: well {int(w)} has no usable field (all {len(gw)} excluded: {_why(gw)}); "
                                 f"this condition rests on {g.loc[ok, 'well'].nunique()} well(s) in {n}")
    live = t[~exc]
    for m, met in metrics_live.items():
        if met.column == "count" or met.column not in live:
            continue
        v = pd.to_numeric(live[met.column], errors="coerce")
        has = v.notna().groupby([live["n"], live["treatment"]]).any()
        none = has[~has]
        if none.empty:
            continue
        by_tr: dict[str, list[str]] = {}
        for (n, tr) in none.index:
            by_tr.setdefault(tr, []).append(str(n))
        for tr, ns in by_tr.items():
            everywhere = not has.loc[(slice(None), tr)].any()
            what = "no value in any field" + (" of any n" if everywhere else f" in {', '.join(ns)}")
            why = " (e.g. no measurable nuclei)" if met.kind == "fraction" else ""
            notes.append(f"{met.label}: {tr} has {what}{why}; "
                         + ("not testable for this readout" if everywhere else "those n's are left out of its comparisons"))
    return notes


def _why(g: pd.DataFrame) -> str:
    src = g["source"].astype(str) if "source" in g else pd.Series("", index=g.index)
    kinds = sorted({"QC low" if "QC low" in s else ("negative control" if "negative control" in s else "manual")
                    for s in src})
    return ", ".join(kinds)


def _cell_lines(t: pd.DataFrame) -> list[str]:
    if "cell_line" not in t:
        return []
    return sorted(x for x in t["cell_line"].fillna("").astype(str).str.strip().unique() if x)


def build(project: Path, gen: dict, image_table: pd.DataFrame, metrics: dict | None, settings: dict,
          design_table: pd.DataFrame | None = None, duplicates: pd.DataFrame | None = None) -> Path:
    """metrics: the metric registry (metrics.registry) or a legacy {column: label} dict; count is always included.
    design_table: n, treatment, image, well, position, timepoint (optional cell_line) from design.assign (computed
    here if None). duplicates: duplicates.find() table of byte-identical raw images (None = not checked).

    The generation is built in generations/<gN>.tmp and only then swapped in by renaming (fileio.swap_in): a failed
    build never leaves a half-written gN, and a file in gN held open (Excel on Windows) gives a clear error instead
    of a half-deleted folder. Files in the old gN that kerncount did not write go to gN/_kept/ (never deleted)."""
    final = project / "generations" / gen["name"]
    out = final.with_name(gen["name"] + ".tmp")
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    try:
        _build_into(out, project, gen, image_table, metrics, settings, design_table, duplicates)
    except BaseException:
        shutil.rmtree(out, ignore_errors=True)
        raise
    kept = fileio.swap_in(out, final, _known_files(final))
    if kept:
        print(f"  {gen['name']}: {len(kept)} file(s) not written by kerncount kept in {gen['name']}/{fileio.KEPT}/: "
              + ", ".join(kept[:5]) + (" ..." if len(kept) > 5 else ""))
    return final


def _raw_scale(model_status: dict, m: str) -> bool:
    """Was metric m analysed on the raw scale (well means <= 0, models.needs_raw_scale) in any all-n block?"""
    return any(v.get("raw_scale") for k, v in model_status.items() if k == m or k.endswith("/" + m))


GENERATION_FILES = {"README.txt", "generation.json", "manifest.json"}


def _known_files(folder: Path) -> set[str] | None:
    """Files an earlier build wrote into `folder` (its manifest, PDF twins of every plot, and the fixed files);
    None when there is no manifest (unknown origin: everything not rebuilt is kept)."""
    try:
        man = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    paths = {f["path"] for f in man.get("files", []) if isinstance(f, dict) and f.get("path")}
    return paths | {p.rsplit(".", 1)[0] + ".pdf" for p in paths} | GENERATION_FILES


def _build_into(out: Path, project: Path, gen: dict, image_table: pd.DataFrame, metrics: dict | None,
                settings: dict, design_table: pd.DataFrame | None, duplicates: pd.DataFrame | None) -> Path:
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
    # design + experiment metadata (well, position, timepoint, cell_line, meta_*), in design.output_columns order
    dcols = [c for c in design.output_columns(design_table) if c in design.COLUMNS]
    mcols = [c for c in design.output_columns(design_table) if c not in design.COLUMNS]
    # a metadata value already in the image table (e.g. cell_line) is kept where the design leaves it empty
    prior = {c: t[c].fillna("").astype(str).str.strip().to_numpy() for c in mcols if c in t}
    t = t.drop(columns=[c for c in dcols + mcols if c in t]).merge(
        design_table[["n", "treatment", "image"] + dcols + mcols].astype({"image": int}).drop_duplicates(
            ["n", "treatment", "image"]), on=["n", "treatment", "image"], how="left")
    for c in mcols:
        t[c] = t[c].fillna("").astype(str).str.strip()
        if c in prior:
            t[c] = t[c].where(t[c] != "", prior[c])
    src_cols = [c for c in t.columns if c == "source_file" or c.startswith("source_file_")]
    # negative control: only sets positivity cut-offs, never analysed
    neg = str((settings.get("experiment") or {}).get("negative_control", "") or "").strip()
    neg_present = bool(neg) and (t["treatment"] == neg).any()
    if neg_present:
        t.loc[t["treatment"] == neg, "excluded"] = True
        t.loc[t["treatment"] == neg, "source"] = "negative control (cut-off only; not analysed)"
    used_cols = ["n", "treatment", "image", "count", "source", "excluded", "qc_level"] + dcols + mcols
    t[used_cols + src_cols].to_csv(out / "counts_used.csv", index=False)
    if neg_present:
        t = t[t["treatment"] != neg].copy()
    try:
        positions_ok = design.positions_meaningful(settings)
    except ValueError:
        positions_ok = False

    # acceptance per n (on the counts); enforce = leave failing n's out of every all-n statistic
    acc = acceptance.evaluate(t[used_cols], settings, duplicates)
    acc.to_csv(out / "acceptance.csv", index=False)
    acc_sum = acceptance.summary(acc, settings)
    if duplicates is not None:
        duplicates.to_csv(out / "duplicates.csv", index=False)

    metrics_live = {k: m for k, m in reg.items() if m.column in t and pd.to_numeric(t[m.column], errors="coerce").notna().any()}
    live = t[~t["excluded"].astype(bool)]
    live_all = live[~live["n"].isin(acc_sum["excluded"])]
    per_n_summ, per_n_stats = {}, {}
    long_rows: list[dict] = []
    n_by_slug, tr_by_slug = {}, {}
    ns_all = sorted(t["n"].unique())
    ncol = modelplots.n_colours(ns_all)
    lines = _cell_lines(t)
    line_of_n = {}
    if "cell_line" in t:
        for n, g in t.groupby("n"):
            ls = sorted(set(g["cell_line"]) - {""})
            line_of_n[n] = ls[0] if len(ls) == 1 else ("; ".join(ls) if ls else "")
    cov = models.config(settings, live["treatment"].unique()).density
    data_notes = coverage_notes(t, metrics_live)
    for x in data_notes:
        print(f"  note: {x}")
    order_all = analysis.treatment_order(live["treatment"], a["treatment_order"])  # stable colours across plots

    for n in ns_all:
        nd = out / _slug(n)
        nd.mkdir()
        n_by_slug[_slug(n)] = n
        for tr in analysis.treatment_order(t.loc[t["n"] == n, "treatment"], a["treatment_order"]):
            td = nd / _slug(tr)
            td.mkdir()
            tr_by_slug[(_slug(n), _slug(tr))] = tr
            sub = t[(t["n"] == n) & (t["treatment"] == tr)].sort_values("image")
            keep = (["image", "count", "source", "excluded", "qc_level"] + dcols[:2] + [c for c in extra_cols if c in sub]
                    + src_cols)
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
            long_rows += manifest.rows_classical(st, "n", m, n, summ.set_index("treatment")["mean"], met.unit,
                                                 line_of_n.get(n) or None)
            if (nd / f"{m}_twoway_anova.csv").exists():
                long_rows += manifest.rows_twoway(pd.read_csv(nd / f"{m}_twoway_pairs.csv"),
                                                  pd.read_csv(nd / f"{m}_twoway_anova.csv"), "n", m, n, met.unit,
                                                  line_of_n.get(n) or None)
            sub = vals[vals["n"] == n]
            order_n = analysis.treatment_order(sub["treatment"], a["treatment_order"])
            modelplots.fields_by_well(sub, order_n, met.label, nd / f"{m}_fields_by_well",
                                      f"{n}: {met.label}, every field", positions_ok)
            if positions_ok:
                modelplots.position_effect(sub.assign(well_l=sub["well"]), nd / f"{m}_position_effect",
                                           f"{n}: {met.label} by field position (per well)", "well",
                                           {w: modelplots.WELL_COLOURS[i % 6] for i, w in enumerate(sorted(sub["well"].unique()))},
                                           lambda w: f"well {int(w)}")
            if met.kind == "fraction" and cov and cov in sub and cov != met.column:
                try:
                    _density_per_n(sub, met, cov, order_all, nd, n)
                except Exception:  # noqa: BLE001
                    traceback.print_exc()

    # all n: once, or once per cell line (never pooled)
    blocks = []
    mixed_line_note = ""
    if len(lines) > 1:
        for line in lines + (["unknown"] if (live_all["cell_line"] == "").any() else []):
            rows = live_all[live_all["cell_line"] == ("" if line == "unknown" else line)]
            blocks.append((line, out / "all_n" / _slug(line), rows))
    else:
        label = lines[0] if lines else None
        if lines and "cell_line" in live_all and (live_all["cell_line"].fillna("").astype(str).str.strip() == "").any():
            # some n's are labelled with the one cell line, others are not: do not call the pooled block "U251"
            unl = sorted(set(live_all.loc[live_all["cell_line"].fillna("").astype(str).str.strip() == "", "n"]))
            label = f"{lines[0]} + unlabelled"
            mixed_line_note = (f"cell line: {', '.join(unl)} have no cell line set; they are analysed together with "
                               f"the {lines[0]} n's (labelled '{label}'). Set [experiment] cell_line (or the "
                               "cell_line column of layout.csv) if they are a different line: lines are never pooled")
            print(f"  WARNING: {mixed_line_note}")
        blocks.append((label, out / "all_n", live_all))
    (out / "all_n").mkdir()
    all_summ, all_stats, results, model_status, fallbacks = {}, {}, {}, {}, []
    one_n_blocks = []
    for line, ad, rows in blocks:
        ad.mkdir(parents=True, exist_ok=True)
        where = f"All n ({line})" if len(blocks) > 1 else "All n"
        key = (lambda m: m) if len(blocks) == 1 else (lambda m, line=line: f"{line}/{m}")
        cfg = models.config(settings, rows["treatment"].unique())
        if rows["n"].nunique() < 2:
            one_n_blocks.append(line or "all")
        for m, met in metrics_live.items():
            vals = rows.assign(value=pd.to_numeric(rows[met.column], errors="coerce"))
            vals = vals[vals["value"].notna()]
            if vals.empty:
                continue
            res = models.analyse(vals, met, settings, positions_ok)
            results[key(m)] = res
            for x in models.boundary_notes(vals, met):
                x = f"{met.label}{f' ({line})' if len(blocks) > 1 else ''}: {x}"
                if x not in data_notes:
                    data_notes.append(x)
            means = vals.groupby("treatment")["value"].mean()
            main_ref = cfg.primary_ref if cfg.primary_ref in set(vals["treatment"]) else \
                next((c for c in cfg.controls if c in set(vals["treatment"])), None)
            use_mixed = cfg.primary_test == "mixed" and res.ok
            stars_by_ref, star_note = None, ""
            if use_mixed:
                stars_by_ref = {ref: (_model_marks(res, ref), _stars_footer(res, ref)) for ref in cfg.controls
                                if ref in set(vals["treatment"])}
            elif res.descriptive:
                star_note = ""
            elif cfg.primary_test == "mixed":
                star_note = f"Mixed model not used for this metric ({res.model}). "
                fallbacks.append(key(m))
            star_source = "mixed" if use_mixed else ("none" if res.descriptive else "rm_anova")
            res.primary_test = star_source if not res.descriptive else "none"
            summ, st = analysis.all_n(vals, m, met.label, ad, settings, stars_by_ref=stars_by_ref, main_ref=main_ref,
                                      star_note=star_note)
            all_summ[key(m)], all_stats[key(m)] = summ, st
            long_rows += manifest.rows_classical(st, "all_n", m, None, means, met.unit, line)
            if (ad / f"{m}_twoway_anova.csv").exists():
                long_rows += manifest.rows_twoway(pd.read_csv(ad / f"{m}_twoway_pairs.csv"),
                                                  pd.read_csv(ad / f"{m}_twoway_anova.csv"), "all_n", m, None,
                                                  met.unit, line)
            try:
                _write_model_outputs(res, m, met, ad, ncol, cfg, main_ref, star_source, where)
            except Exception as e:  # noqa: BLE001 - never crash a generation build on a plot
                res.notes.append(f"model outputs incomplete ({type(e).__name__}: {e})")
                traceback.print_exc()
            long_rows += manifest.rows_model(res, m, line, line_of_n)
            model_status[key(m)] = res.status_dict() | dict(stars_from=star_source, cell_line=line, metric=m)
            if positions_ok and m == "count":
                modelplots.position_effect(vals, ad / f"{m}_position_effect",
                                           f"{where}: {met.label} by field position (per n)", "n", ncol, manifest.short)

    # readout ladder: every intensity layer's effect side by side, per channel, at all n and in every n
    layer_info = _write_layers(out, blocks, results, metrics_live, settings, order_all, ncol, n_by_slug)
    stats_long = manifest.write_stats_long(out, long_rows)
    cfg_all = models.config(settings, live_all["treatment"].unique())
    dup_lines = duplicates_mod.summary_lines(duplicates) if duplicates is not None else []
    dup_aff = duplicates_mod.affected(duplicates) if duplicates is not None else {}
    notes = []
    if acc_sum["excluded"]:
        notes.append(f"acceptance enforced: {', '.join(acc_sum['excluded'])} left out of all-n statistics")
    if fallbacks:
        notes.append(f"mixed model not available for {', '.join(fallbacks)}: RM ANOVA drives the stars there")
    if neg_present:
        notes.append(f"negative control '{neg}' left out of every test, plot and summary (it only sets the positivity "
                     "cut-off)")
    elif neg:
        notes.append(f"negative control '{neg}' set but not present in the data")
    if mixed_line_note:
        notes.append(mixed_line_note)
    notes += data_notes
    if len(blocks) > 1:
        notes.append(f"cell lines are never pooled: all-n analysis run separately for {', '.join(b[0] for b in blocks)} "
                     f"(folders {', '.join((b[1].relative_to(out)).as_posix() for b in blocks)})")
    if one_n_blocks:
        notes.append(f"only one n ({', '.join(one_n_blocks)}): all-n estimates from that n, across-n CIs, p-values and "
                     f"stars {models.NEEDS_2N}")
    for line in dup_lines:
        notes.append(f"DUPLICATE RAW IMAGES: {line}")
    if dup_aff:
        notes.append("duplicate images shared across conditions in: " +
                     "; ".join(f"{n}: {', '.join(trs)}" for n, trs in dup_aff.items()) +
                     " - one of the copies is not the image it claims to be")
        # duplicates combined with a missing reference (the real n3 case plus an n without LyoVec): say how many
        # n's are left with a reference that is neither missing nor a copy
        for line, _, rows in blocks:
            for ref in cfg_all.controls:
                have = sorted(set(rows.loc[rows["treatment"] == ref, "n"]))
                if not have or not any(ref in dup_aff.get(n, []) for n in have):
                    continue
                clean = [n for n in have if ref not in dup_aff.get(n, [])]
                lab = f" ({line})" if len(blocks) > 1 and line else ""
                notes.append(f"comparisons vs {ref}{lab}: only {len(clean)} n ({', '.join(clean) or 'none'}) "
                             f"{'has' if len(clean) == 1 else 'have'} a {ref} that is neither missing nor a duplicate"
                             + (f" - {models.NEEDS_2N} for a trustworthy test" if len(clean) < 2 else ""))
    _excel(out, gen, t, used_cols, a, cfg_all, reg, metrics_live, all_summ, all_stats, per_n_summ, per_n_stats,
           results, acc, acc_sum, model_status, out / "all_n", notes, duplicates)

    src = "fully automatic (auto_count)" if gen["source"] == "auto" else \
        f"previous generation + manual values from column {gen['source']} of counts.csv ({gen['changes']} cells)"
    unit = cfg_all.unit
    (out / "README.txt").write_text(README.format(
        name=gen["name"], source_text=src, built=f"{datetime.now():%Y-%m-%d %H:%M}", primary_test=cfg_all.primary_test,
        primary_text=PRIMARY_TEXT[cfg_all.primary_test], unit_text=UNIT_TEXT[unit],
        line_text=(" (one sub-folder per cell line: " + ", ".join(b[1].name for b in blocks) + ")") if len(blocks) > 1 else "",
        extra_notes=("\nNotes\n" + "\n".join(f"  - {x}" for x in notes) + "\n") if notes else ""), encoding="utf-8")
    (out / "generation.json").write_text(json.dumps(dict(name=gen["name"], source=gen["source"], fingerprint=gen["fingerprint"],
                                                         changes=gen["changes"], built=datetime.now().isoformat(timespec="seconds"),
                                                         **({"inputs": gen["inputs"]} if gen.get("inputs") else {}),
                                                         **({"rebuilt_because": gen["rebuilt_because"]}
                                                            if gen.get("rebuilt_because") else {})),
                                                    indent=2, default=str), encoding="utf-8")
    primary = dict(test=cfg_all.primary_test, family=f"vs {cfg_all.primary_ref}", reference=cfg_all.primary_ref,
                   fallback_metrics=fallbacks, controls=cfg_all.controls,
                   classes={c: [t_ for t_, k in cfg_all.classes.items() if k == c] for c in cfg_all.class_order},
                   classes_mode=cfg_all.classes_mode,
                   positive_control=cfg_all.positive_control or None,
                   vehicle_check=[f"{v} vs {c}" for v, c in cfg_all.vehicle_pairs],
                   equivalence_margin=cfg_all.margin, equivalence_margin_pp=cfg_all.margin_pp)
    extra = dict(
        model_unit=unit,
        estimands={m: (lambda met_: dict(effect=metrics_mod.estimand(met_, unit),
                                         mean=metrics_mod.estimand(met_, unit, "mean"),
                                         headline={"fraction": "pp_difference", "difference": "difference"}.get(
                                             met_.kind, "ratio"),
                                         raw_scale=met_.kind == "difference"))(
                       models.raw_scale_metric(met) if _raw_scale(model_status, m) else met)
                   for m, met in metrics_live.items()},
        density_covariate=cfg_all.density or None,
        n_count=int(live_all["n"].nunique()),
        cell_lines=lines,
        all_n_folders=[dict(cell_line=b[0], path=b[1].relative_to(out).as_posix(), n_count=int(b[2]["n"].nunique()))
                       for b in blocks],
        negative_control=dict(treatment=neg or None, present=neg_present),
        layers=layer_info,
        duplicates=dict(checked=duplicates is not None,
                        images=int(len(duplicates)) if duplicates is not None else 0,
                        cross_condition=int(duplicates["cross_condition"].astype(bool).sum()) if duplicates is not None
                        and len(duplicates) else 0,
                        affected=dup_aff, summary=dup_lines))
    line_by_slug = {b[1].name: b[0] for b in blocks} if len(blocks) > 1 else None
    manifest.build_manifest(out, gen, metrics_live, n_by_slug, tr_by_slug, primary, model_status, acc_sum, notes,
                            extra=extra, line_by_slug=line_by_slug, line_of_n=line_of_n or None)
    return out


def _excel(out, gen, t, used_cols, a, cfg, reg, metrics_live, all_summ, all_stats, per_n_summ, per_n_stats, results,
           acc, acc_sum, model_status, ad, notes=(), duplicates=None):
    def base_metric(k):
        return k.split("/", 1)[1] if "/" in k else k

    def line_of(k):
        return k.split("/", 1)[0] if "/" in k else ""

    other = [m for m in all_summ if base_metric(m) != "count"]
    unit = cfg.unit
    with pd.ExcelWriter(out / "results.xlsx", engine="openpyxl") as xw:
        lines = [f"Generation {gen['name']}", f"Source: {gen['source']}",
                 f"Manual changes in this step: {gen['changes']}",
                 f"Built {datetime.now():%Y-%m-%d %H:%M}",
                 f"Control groups (main first): {', '.join(analysis.control_groups(a))}",
                 f"Primary test (drives the stars): {cfg.primary_test}; main claim: vs {cfg.primary_ref}",
                 f"Model unit: {UNIT_TEXT[unit]}",
                 "Estimands: " + "; ".join(f"{m}: {metrics_mod.estimand(met, unit)}" for m, met in metrics_live.items()),
                 "Within-n stats treat fields as technical replicates (descriptive); Dunnett vs each control.",
                 "All-n stats: each n = one biological replicate.",
                 "  Mixed model (sheets 'Model - ...'): wells in n (or fields, see model unit); log scale for counts and "
                 "intensities (results as ratios of means), logit for fractions (headline = difference in percentage "
                 "points; odds ratio and % of reference as secondary); n-to-n variability estimated separately per "
                 "treatment class; Satterthwaite t; p adjusted within each family (single-step multivariate t).",
                 "  RM one-way ANOVA + Dunnett vs each control (Prism standard, raw scale); ratio paired t-tests (Holm) as a check.",
                 "  Vehicle check: equivalence test (TOST, 90% CI within the margin) - 'not significant' is not 'no effect'.",
                 "  Robustness: leave-one-n-out (is the conclusion the same whichever n is dropped?) and the effect in every n.",
                 "  Fractions: a density-adjusted model (covariate " + (cfg.density or "off") + ") is shown beside the raw "
                 "estimate; flagged when the conditions do not overlap in density.",
                 "Two-way ANOVA (extra): per n treatment x well, all n treatment x n; fields = replicates; Tukey all pairs.",
                 "Treatment classes (" + cfg.classes_mode + "): " + "; ".join(
                     f"{c}: {', '.join(x for x, k in cfg.classes.items() if k == c)}" for c in cfg.class_order)]
        for m, s in model_status.items():
            lines.append(f"Model for {m}: {s['model']} (stars from {s['stars_from']})" +
                         (f"; notes: {'; '.join(s['notes'])}" if s["notes"] else ""))
        fails = acceptance.failing(acc)
        lines.append("Acceptance: " + ("all n pass" if not fails else f"failing: {', '.join(fails)}") +
                     (" (enforced: left out of all-n statistics)" if acc_sum["enforce"] and fails else
                      " (report only)" if fails else ""))
        for nt in notes:
            lines.append(f"NOTE: {nt}")
        pd.DataFrame({"": lines}).to_excel(xw, sheet_name="About", index=False, header=False)
        t[used_cols].to_excel(xw, sheet_name="Counts used", index=False)

        def cat_summ(d, ms):
            return pd.concat([d[m].assign(metric=base_metric(m), cell_line=line_of(m)) for m in ms])

        counts = [m for m in all_summ if base_metric(m) == "count"]
        if counts:
            cat_summ(all_summ, counts).to_excel(xw, sheet_name="Count - all n", index=False)
            cat_summ(all_stats, counts).to_excel(xw, sheet_name="Count - stats all n", index=False)
        if "count" in per_n_summ:
            pd.concat(per_n_summ["count"]).to_excel(xw, sheet_name="Count - per n", index=False)
            pd.concat(per_n_stats["count"]).to_excel(xw, sheet_name="Count - stats per n", index=False)
        inten = [m for m in other if base_metric(m).startswith("intensity_")]
        rest = [m for m in other if not base_metric(m).startswith("intensity_")]
        for name, ms in (("Intensity", inten), ("Other metrics", rest)):
            if not ms:
                continue
            cat_summ(all_summ, ms).to_excel(xw, sheet_name=f"{name} - all n", index=False)
            cat_summ(all_stats, ms).to_excel(xw, sheet_name=f"{name} - stats all n", index=False)
            pn = [pd.concat(per_n_summ[base_metric(m)]).assign(metric=base_metric(m)) for m in dict.fromkeys(
                base_metric(x) for x in ms) if m in per_n_summ]
            if pn:
                pd.concat(pn).to_excel(xw, sheet_name=f"{name} - per n", index=False)
        for level, files in (("all n", sorted(ad.rglob("*_twoway_anova.csv"))),
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
                        elif f.parent != ad:
                            tab.insert(0, "cell_line", f.parent.name)
                        tabs.append(tab)
                if tabs:
                    pd.concat(tabs).to_excel(xw, sheet_name=f"Two-way {kind} - {level}", index=False)

        def cat(attr, filt=None, sub=None):
            tabs = []
            for k, r in results.items():
                obj = getattr(r, sub) if sub else r
                if obj is None:
                    continue
                x = getattr(obj, attr)
                if x is None or len(x) == 0:
                    continue
                x = x if filt is None else filt(x)
                m = base_metric(k)
                tabs.append(x.assign(metric=m, cell_line=line_of(k), unit=reg[m].unit if m in reg else ""))
            return pd.concat(tabs) if tabs else None

        sheets = [
            ("Model - means (absolute)", cat("emm")),
            ("Model - ratios vs controls", cat("contrasts", lambda c: c[c["family_kind"] == "vs_control"])),
            ("Comparison families", cat("contrasts")),
            ("Density-adjusted means", cat("emm", sub="density")),
            ("Density-adjusted effects", cat("contrasts", sub="density")),
            ("Variance components", cat("varcomp")),
            ("Robustness - stability", cat("loo_summary")),
            ("Robustness - leave-one-out", cat("loo")),
            ("Robustness - effect per n", cat("per_n")),
        ]
        for name, tab in sheets:
            if tab is not None and len(tab):
                tab.to_excel(xw, sheet_name=name, index=False)
        ms = pd.DataFrame([dict(key=m, **{k: ("; ".join(v) if isinstance(v, list) else v) for k, v in s.items()})
                           for m, s in model_status.items()])
        if len(ms):
            ms.to_excel(xw, sheet_name="Model status", index=False)
        acc.to_excel(xw, sheet_name="Acceptance", index=False)
        if duplicates is not None and len(duplicates):
            duplicates.to_excel(xw, sheet_name="Duplicate images", index=False)
        _autosize(xw)
