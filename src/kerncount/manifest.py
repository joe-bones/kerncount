"""Output contracts of one generation folder for the browser UI: manifest.json and stats_long.csv.

manifest.json registers every plot (.png; the matching .pdf is given as `pdf`) and every table (.csv, .xlsx) in the
generation folder with path, level, n, treatment, metric, kind, title and a one-sentence description. Files are
described from their names (one rule per file type), so anything any module writes is listed. The schema is in
docs/statistics.md.
"""
from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

STATS_LONG_COLUMNS = ["level", "n", "metric", "family", "test", "comparison", "estimate", "estimate_type", "ci_low",
                      "ci_high", "p", "p_adj", "significance", "N", "model", "note", "estimand", "cell_line"]
ESTIMATE_TYPES = ("ratio", "odds_ratio", "difference", "pp_difference", "mean")
LEVELS = ("treatment", "n", "all_n", "generation")
KINDS = ("absolute", "relative", "model", "contrast", "robustness", "diagnostic", "design", "acceptance", "data")


# ---------------------------------------------------------------- describing files by name

def _family_text(fid: str) -> str:
    if fid.startswith("vs_positive_control"):
        return "each treatment compared with the positive control"
    if fid == "vehicle_check":
        return "equivalence test of the vehicle against its control"
    if fid.startswith("within_"):
        return f"all pairs within the class {fid[len('within_'):]}"
    if fid.startswith("vs_"):
        return f"each treatment compared with {fid[3:]}"
    return fid


# (regex on the part after "<metric>_", level-set, kind, title template, description template)
_METRIC_RULES = [
    (r"", ("n",), "absolute", "{label}: per-field values and stars vs the first control",
     "Bars show the mean of the fields in this experiment (dots = fields) with stars from a within-experiment test."),
    (r"", ("all_n",), "absolute", "{label}: mean of the experiments",
     "Bars show the mean of the per-experiment means (dots = experiments) with stars from the primary test."),
    (r"vs_(?P<x>.+)", ("n",), "absolute", "{label}: compared with {x}",
     "Same values as the main plot, with stars for the comparison with {x} in this experiment."),
    (r"pct_of_(?P<x>.+)", ("all_n",), "relative", "{label} as % of {x}",
     "Each experiment's value as a percentage of its own {x}, averaged over experiments."),
    (r"prism_pct_of_(?P<x>.+)", None, "relative", "Prism table: {label} as % of {x}",
     "One column per treatment, one row per experiment, values as % of {x}; paste into a Prism column table."),
    (r"prism", None, "data", "Prism table: {label}",
     "One column per treatment with the replicate values, ready to paste into a Prism column table."),
    (r"summary", None, "absolute", "Summary of {label}",
     "Mean, SD, SEM and count per treatment, with % of each control where available."),
    (r"stats", None, "contrast", "Statistical tests for {label}",
     "Every classical test (ANOVA, Dunnett, ratio t-tests) with raw and adjusted p-values."),
    (r"twoway", None, "contrast", "{label}: two-way ANOVA",
     "Grouped bars and a matrix of all pairwise Tukey comparisons; fields are counted as replicates."),
    (r"twoway_anova", None, "contrast", "Two-way ANOVA table for {label}",
     "Sums of squares, F and p for treatment, the second factor and their interaction."),
    (r"twoway_pairs", None, "contrast", "Tukey all-pairs comparisons for {label}",
     "Every pair of treatments compared within the two-way ANOVA, with adjusted p-values."),
    (r"fields_by_well", None, "design", "{label}: every field by well and position",
     "Each dot is one field, coloured by well and shaped by position, to show where variation within this experiment comes from."),
    (r"position_effect", None, "design", "{label}: effect of field position",
     "Fields divided by their treatment's mean and averaged per position, showing whether some positions systematically give higher values."),
    (r"model_absolute", None, "model", "{label}: model estimates in absolute units",
     "Model-estimated typical value (geometric mean) per treatment with 95% CI, with each experiment's mean as a coloured dot."),
    (r"model_means", None, "model", "Model estimates of {label} per treatment",
     "Model-estimated typical value per treatment in absolute units with 95% confidence interval."),
    (r"model_info", None, "model", "Model used for {label}",
     "Which mixed model was fitted, its formula, degrees-of-freedom rule, convergence and notes."),
    (r"contrasts", None, "contrast", "Model comparisons for {label}",
     "Every model-based comparison (ratios with 95% CI and adjusted p) in every comparison family, including the equivalence test."),
    (r"forest_(?P<x>.+)", None, "contrast", "{label}: {fam}",
     "Ratios with 95% CI on a log axis (1 = no difference); small coloured dots show the estimate when one experiment is left out."),
    (r"variance_components", None, "model", "{label}: sources of variation",
     "How much of the variation comes from experiments, treatment-by-experiment, wells, field positions and single fields."),
    (r"diagnostics", None, "diagnostic", "{label}: model checks",
     "Residual plots and estimated experiment effects to check that the model fits the data reasonably."),
    (r"loo", None, "robustness", "{label}: leave-one-experiment-out",
     "Every comparison re-estimated with each experiment left out in turn."),
    (r"loo_summary", None, "robustness", "{label}: stability of the conclusions",
     "For each comparison, whether direction and significance stay the same whichever experiment is left out."),
    (r"per_n_ratios", None, "robustness", "{label}: effect in each experiment",
     "The ratio of each treatment to each control, separately for every experiment."),
    (r"per_n_ratio_vs_(?P<x>.+)", None, "robustness", "{label}: effect per experiment vs {x}",
     "Each experiment's own effect vs {x} (coloured dots, joined per experiment; thin lines = technical interval from "
     "the fields) next to the model estimate (black bar)."),
    (r"density_adjusted_means", None, "model", "{label}: density-adjusted estimates per treatment",
     "Model estimates at the average local crowding (density covariate), beside the raw estimates; raw stays the headline."),
    (r"density_adjusted_contrasts", None, "contrast", "{label}: density-adjusted comparisons",
     "Comparisons adjusted for local crowding, with the density overlap of each pair and a warning when the "
     "conditions do not overlap in density."),
    (r"density", ("n",), "design", "{label} vs local crowding",
     "Each field's value against its crowding (one colour per treatment, line = descriptive fit within this experiment)."),
    (r"density", ("all_n",), "model", "{label} vs local crowding (all n)",
     "Each well's value against its crowding (one colour per treatment, marker = experiment, lines = density-adjusted "
     "model per treatment)."),
]

_ROOT = {
    "counts_used.csv": ("data", "Counts used in this generation",
                        "The value used for every image, where it came from, and its well, position and timepoint."),
    "results.xlsx": ("data", "Results workbook", "All tables of this generation in one Excel workbook."),
    "stats_long.csv": ("data", "All test results (long table)",
                       "Every statistical result of this generation in one table, one row per test or comparison."),
    "acceptance.csv": ("acceptance", "Acceptance criteria per experiment",
                       "Whether each experiment passes the quality criteria set in the acceptance settings."),
    "duplicates.csv": ("acceptance", "Duplicate raw images",
                       "Raw image files that are byte-identical to another image (e.g. a copied folder); empty if none."),
}


def _metric_split(stem: str, metric_ids: list[str]) -> tuple[str | None, str]:
    for m in sorted(metric_ids, key=len, reverse=True):
        if stem == m:
            return m, ""
        if stem.startswith(m + "_"):
            return m, stem[len(m) + 1:]
    return None, stem


def describe(rel: str, metrics: dict, n_by_slug: dict, tr_by_slug: dict, line_by_slug: dict | None = None,
             line_of_n: dict | None = None) -> dict:
    """line_by_slug: {folder: cell line} of all_n/<cell line>/ sub-folders (several cell lines);
    line_of_n: {n: cell line} for per-n files."""
    parts = rel.split("/")
    name = parts[-1]
    stem, ext = name.rsplit(".", 1) if "." in name else (name, "")
    entry = dict(path=rel, level="generation", n=None, treatment=None, metric=None, kind="data",
                 title=name, description="", cell_line=None)
    if len(parts) == 1:
        k, t, d = _ROOT.get(name, ("data", name, "A file of this generation."))
        entry.update(kind=k, title=t, description=d)
        return entry
    if parts[0] == "all_n":
        entry["level"] = "all_n"
        if len(parts) == 3 and line_by_slug:
            entry["cell_line"] = line_by_slug.get(parts[1], parts[1])
    elif len(parts) == 2:
        entry.update(level="n", n=n_by_slug.get(parts[0], parts[0]))
    else:
        entry.update(level="treatment", n=n_by_slug.get(parts[0], parts[0]),
                     treatment=tr_by_slug.get((parts[0], parts[1]), parts[1]))
    if entry["n"] and line_of_n:
        entry["cell_line"] = line_of_n.get(entry["n"]) or None
    if entry["level"] == "treatment":
        if name == "fields.csv":
            entry.update(kind="data", title=f"Fields of {entry['treatment']} in {entry['n']}",
                         description="One row per image with its count, source, QC level and every measured value.")
        elif name == "count_fields.png":
            entry.update(kind="absolute", metric="count", title=f"Nuclei per image: {entry['treatment']} in {entry['n']}",
                         description="Nuclei counted in each image of this treatment; bar colour = well, label = image and field position.")
        return entry
    if name == "n_colours.csv":
        entry.update(kind="data", title="Colour of each experiment", description="The dot colour used for each experiment in the all-n plots.")
        return entry
    m, rest = _metric_split(stem, list(metrics))
    if m is None:
        entry["description"] = "A file of this generation."
        return entry
    met = metrics[m]
    entry["metric"] = m
    label = met.label if hasattr(met, "label") else str(met)
    for pat, levels, kind, title, desc in _METRIC_RULES:
        if levels is not None and entry["level"] not in levels:
            continue
        mm = re.fullmatch(pat, rest)
        if not mm:
            continue
        x = mm.groupdict().get("x", "") or ""
        fam = _family_text(x) if "forest" in pat else ""
        if pat.startswith("forest") and x == "vehicle_check":
            desc = "Ratio of the vehicle to its control with 90% CI against the equivalence margin (shaded): inside = equivalent."
        entry.update(kind=kind, title=title.format(label=label, x=x, fam=fam), description=desc.format(x=x))
        break
    else:
        entry.update(title=f"{label}: {rest.replace('_', ' ')}", description="A file of this generation.")
    if entry["level"] == "n":
        entry["title"] = f"{short(entry['n'])}: {entry['title']}"
    elif entry["level"] == "all_n":
        entry["title"] = (f"All n ({entry['cell_line']}): {entry['title']}" if entry.get("cell_line")
                          else f"All n: {entry['title']}")
    return entry


def short(n) -> str:
    return str(n).split("_")[0] if n else ""


def build_manifest(out: Path, gen: dict, metrics: dict, n_by_slug: dict, tr_by_slug: dict, primary: dict,
                   model_status: dict, acceptance: dict, notes: list[str] | None = None, extra: dict | None = None,
                   line_by_slug: dict | None = None, line_of_n: dict | None = None) -> dict:
    """extra: further top-level keys (model_unit, estimands, cell_lines, all_n_folders, negative_control,
    duplicates, ...; see docs/statistics.md)."""
    files = []
    for f in sorted(out.rglob("*")):
        if not f.is_file() or f.suffix.lower() not in (".png", ".csv", ".xlsx"):
            continue
        rel = f.relative_to(out).as_posix()
        e = describe(rel, metrics, n_by_slug, tr_by_slug, line_by_slug, line_of_n)
        pdf = f.with_suffix(".pdf")
        if f.suffix.lower() == ".png" and pdf.exists():
            e["pdf"] = pdf.relative_to(out).as_posix()
        files.append(e)
    man = dict(schema="kerncount.generation.manifest/1", generation=gen["name"], source=gen["source"],
               built=datetime.now().isoformat(timespec="seconds"),
               metrics=[dict(id=k, label=m.label, unit=m.unit, kind=m.kind, description=getattr(m, "description", ""))
                        for k, m in metrics.items()],
               primary=primary, model_status=model_status, acceptance=acceptance, notes=list(notes or []),
               **(extra or {}), files=files)
    (out / "manifest.json").write_text(json.dumps(_clean(man), indent=2, ensure_ascii=False), encoding="utf-8")
    return man


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


# ---------------------------------------------------------------- stats_long

def _family_of(test: str) -> str:
    m = re.search(r"vs (\S+)", test)
    if "ANOVA" in test and "Dunnett" not in test:
        return "omnibus"
    if "Tukey" in test:
        return "all pairs (two-way)"
    return f"vs {m.group(1)}" if m else ""


def _classical_estimand(test: str, unit: str) -> str:
    if "ratio paired t-test" in test:
        return f"ratio of per-n mean {unit} (geometric mean over n)"
    if "Dunnett" in test or "Tukey" in test:
        return f"difference in mean {unit}"
    return ""


def rows_classical(st: pd.DataFrame, level: str, metric: str, n=None, means: pd.Series | None = None,
                   unit: str = "", cell_line=None) -> list[dict]:
    """Rows from an analysis.per_n / all_n *_stats table."""
    out = []
    for r in st.itertuples():
        test = r.test
        ttype = getattr(r, "effect_type", "")
        if isinstance(ttype, str) and ttype.startswith("ratio"):
            et = "ratio"
        elif "ANOVA" in test and "Dunnett" not in test:
            et = ""
        else:
            et = "difference"
        est = getattr(r, "effect", np.nan)
        if (est is None or not np.isfinite(est)) and means is not None and " vs " in str(r.comparison):
            a, b = str(r.comparison).split(" vs ", 1)
            if a in means and b in means:
                est = float(means[a] - means[b])
        note = getattr(r, "note", "")
        out.append(dict(level=level, n=n, metric=metric, family=_family_of(test), test=test, comparison=r.comparison,
                        estimate=est, estimate_type=et, ci_low=getattr(r, "ci_low", np.nan),
                        ci_high=getattr(r, "ci_high", np.nan), p=r.p, p_adj=r.p_adj, significance=r.significance,
                        N=getattr(r, "N", np.nan), model="classical", note=note if isinstance(note, str) else "",
                        estimand=_classical_estimand(test, unit), cell_line=cell_line))
    return out


def rows_twoway(pairs: pd.DataFrame, anova: pd.DataFrame, level: str, metric: str, n=None, unit: str = "",
                cell_line=None) -> list[dict]:
    out = []
    for r in anova.itertuples():
        if r.source == "residual":
            continue
        out.append(dict(level=level, n=n, metric=metric, family="omnibus", test=f"two-way ANOVA ({r.model})",
                        comparison=r.source, estimate=np.nan, estimate_type="", ci_low=np.nan, ci_high=np.nan, p=r.p,
                        p_adj=r.p, significance=r.significance if isinstance(r.significance, str) else "", N=np.nan,
                        model="classical (fields as replicates)", note=str(getattr(r, "note", "") or ""),
                        estimand="", cell_line=cell_line))
    for r in pairs.itertuples():
        out.append(dict(level=level, n=n, metric=metric, family="all pairs (two-way)", test=r.test,
                        comparison=r.comparison, estimate=r.mean_diff, estimate_type="difference", ci_low=r.ci_low,
                        ci_high=r.ci_high, p=np.nan, p_adj=r.p_adj,
                        significance=r.significance if isinstance(r.significance, str) else "",
                        N=r.matched_levels, model="classical (fields as replicates)",
                        note=str(getattr(r, "note", "") or ""), estimand=f"difference in mean {unit}",
                        cell_line=cell_line))
    return out


def _f(x):
    return x if (x is not None and np.isfinite(x)) else np.nan


def rows_model(res, metric: str, cell_line=None, line_of_n: dict | None = None) -> list[dict]:
    """Rows of the hierarchical model: estimated means, every comparison (fractions: one row with the difference in
    percentage points, the headline, plus one with the odds ratio), density-adjusted comparisons, leave-one-n-out
    and the per-n effects (descriptive)."""
    from .metrics import estimand

    out = []
    model = res.model
    frac = res.scale is not None and res.scale.kind == "fraction"
    unit = res.unit
    est_mean = estimand(res.metric, unit, "mean")
    est_eff = estimand(res.metric, unit)

    def base(**k):
        return dict(level="all_n", n=None, metric=metric, model=model, cell_line=cell_line) | k

    for r in res.emm.itertuples():
        note = f"{r.mean_type}; df {r.df:.1f}" + (f"; {r.note}" if isinstance(r.note, str) and r.note else "")
        out.append(base(family="estimated means", test="mixed model estimated marginal mean (back-transformed)"
                        if res.ok else "descriptive mean (one n)", comparison=r.treatment, estimate=r.estimate,
                        estimate_type="mean", ci_low=r.ci_low, ci_high=r.ci_high, p=np.nan, p_adj=np.nan,
                        significance="", N=r.N, note=note, estimand=est_mean))
    for r in res.contrasts.itertuples():
        note = (r.note or "") if isinstance(r.note, str) else ""
        dfn = f"df {r.df:.1f}" if np.isfinite(r.df) and r.df >= 1 else ""
        note = "; ".join(x for x in (note, dfn) if x)
        if frac:
            out.append(base(family=r.family, test=r.test, comparison=r.comparison, estimate=r.pp_difference,
                            estimate_type="pp_difference", ci_low=r.pp_ci_low, ci_high=r.pp_ci_high, p=r.p,
                            p_adj=r.p_adj, significance=r.significance, N=r.N,
                            note="; ".join(x for x in (note, f"{r.pct_group:.2f}% vs {r.pct_comparator:.2f}%"
                                                       if np.isfinite(r.pct_group) else "") if x),
                            estimand=est_eff))
            out.append(base(family=r.family, test=r.test, comparison=r.comparison, estimate=r.estimate,
                            estimate_type="odds_ratio", ci_low=r.ci_low, ci_high=r.ci_high, p=r.p, p_adj=r.p_adj,
                            significance=r.significance, N=r.N, note=note,
                            estimand="odds ratio (pooled per well; logit model; same test as the pp row)"))
        else:
            out.append(base(family=r.family, test=r.test, comparison=r.comparison, estimate=r.estimate,
                            estimate_type=r.estimate_type, ci_low=r.ci_low, ci_high=r.ci_high, p=r.p, p_adj=r.p_adj,
                            significance=r.significance, N=r.N, note=note, estimand=est_eff))
    dens = getattr(res, "density", None)
    if dens is not None and len(dens.contrasts):
        for r in dens.contrasts.itertuples():
            flag = r.density_flag if isinstance(r.density_flag, str) else ""
            note = "; ".join(x for x in (f"density overlap {r.density_overlap:.2f}" if np.isfinite(r.density_overlap)
                                         else "", flag, r.note if isinstance(r.note, str) else "") if x)
            out.append(base(family=r.family, test=f"{r.test}; adjusted for {dens.covariate}", comparison=r.comparison,
                            estimate=r.pp_difference, estimate_type="pp_difference", ci_low=r.pp_ci_low,
                            ci_high=r.pp_ci_high, p=r.p, p_adj=r.p_adj, significance=r.significance, N=r.N,
                            model=dens.model, note=note,
                            estimand=est_eff + f" at the average {dens.covariate} (density-adjusted)"))
    for r in res.loo.itertuples():
        if not r.family:
            continue
        pp = frac and np.isfinite(getattr(r, "pp_difference", np.nan))
        out.append(base(family=f"{r.family} (leave-one-n-out)",
                        test=f"{'mixed model' if res.ok and res.primary_test == 'mixed' else 'RM ANOVA + Dunnett'}"
                             f" without {r.dropped_n}", comparison=r.comparison,
                        estimate=r.pp_difference if pp else r.estimate,
                        estimate_type="pp_difference" if pp else r.estimate_type,
                        ci_low=np.nan if pp else r.ci_low, ci_high=np.nan if pp else r.ci_high, p=r.p, p_adj=r.p_adj,
                        significance=_stars(r.p_adj), N=r.N, note=f"without {r.dropped_n}",
                        estimand=est_eff if r.estimate_type != "difference" else "difference in mean (raw scale)"))
    for r in res.per_n.itertuples():
        line = (line_of_n or {}).get(r.n, cell_line)
        if frac:
            out.append(dict(level="n", n=r.n, metric=metric, family=f"vs {r.control} (per-n difference)",
                            test="difference of pooled fractions within this n (descriptive: technical replicates)",
                            comparison=f"{r.treatment} vs {r.control}", estimate=r.pp_difference,
                            estimate_type="pp_difference", ci_low=r.tech_ci_low, ci_high=r.tech_ci_high, p=np.nan,
                            p_adj=np.nan, significance="", N=r.fields, model="descriptive",
                            note=f"{r.interval}; ratio of fractions {r.ratio_of_means:.3g}",
                            estimand="difference in percentage points (pooled over the fields of this n)",
                            cell_line=line))
        else:
            out.append(dict(level="n", n=r.n, metric=metric, family=f"vs {r.control} (per-n ratio)",
                            test="ratio of mean values within this n (descriptive: technical replicates)",
                            comparison=f"{r.treatment} vs {r.control}", estimate=r.ratio_of_means,
                            estimate_type=r.estimate_type, ci_low=r.tech_ci_low, ci_high=r.tech_ci_high, p=np.nan,
                            p_adj=np.nan, significance="", N=r.fields, model="descriptive",
                            note=f"{r.interval}; model-scale ratio {r.ratio:.3g}",
                            estimand=f"ratio of mean {res.metric.unit or res.metric.label} (fields of this n)",
                            cell_line=line))
    return out


def _stars(p):
    from .models import stars
    return stars(p)


def write_stats_long(out: Path, rows: list[dict]) -> pd.DataFrame:
    t = pd.DataFrame(rows, columns=STATS_LONG_COLUMNS)
    t.to_csv(out / "stats_long.csv", index=False)
    return t
