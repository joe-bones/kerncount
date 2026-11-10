"""Summaries, statistics, Prism tables and plots for one metric at three levels:
per treatment (fields), per n (treatments compared within one experiment) and all n combined.

Statistics
    within one n:  fields are technical replicates -> one-way ANOVA + Dunnett vs reference (descriptive only)
    across n:      each n contributes one value per treatment (mean of its fields)
                   - repeated-measures one-way ANOVA on n's that have every treatment
                   - ratio paired t-test of each treatment vs reference (log ratio), Holm-adjusted;
                     uses every n that has both groups, so a missing group does not drop the whole n
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

N_COLOURS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b", "#e377c2", "#17becf", "#7f7f7f", "#bcbd22"]


def stars(p: float) -> str:
    if p is None or np.isnan(p):
        return ""
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def holm(p: list[float]) -> list[float]:
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    idx = np.argsort(p[ok])
    m = ok.sum()
    adj = np.empty(m)
    run = 0.0
    for rank, i in enumerate(idx):
        run = max(run, min(1.0, (m - rank) * p[ok][i]))
        adj[i] = run
    out[ok] = adj
    return out.tolist()


def treatment_order(present, preferred: list[str]) -> list[str]:
    present = list(dict.fromkeys(present))
    return [t for t in preferred if t in present] + sorted(t for t in present if t not in preferred)


def prism_table(groups: dict[str, list[float]], row_titles: list[str] | None = None) -> pd.DataFrame:
    """Prism 'Column' table: one column per group, replicate values down the rows."""
    length = max((len(v) for v in groups.values()), default=0)
    df = pd.DataFrame({k: list(v) + [np.nan] * (length - len(v)) for k, v in groups.items()})
    if row_titles is not None:
        df.insert(0, "", list(row_titles) + [""] * (length - len(row_titles)))
    return df


# ---------------------------------------------------------------- plotting

def _bar_plot(order, values: dict, err: str, ylabel: str, title: str, out: Path, point_colour=None, ref_line=None):
    fig, ax = plt.subplots(figsize=(max(4, 0.8 * len(order) + 1.5), 4))
    rng = np.random.default_rng(0)
    for i, t in enumerate(order):
        v = np.asarray([x for x in values.get(t, []) if not np.isnan(x)])
        if len(v) == 0:
            continue
        m = v.mean()
        e = (v.std(ddof=1) if err == "SD" else v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0
        ax.bar(i, m, color="#d9d9d9", edgecolor="black", width=0.65, zorder=1)
        ax.errorbar(i, m, yerr=e, color="black", capsize=4, zorder=3)
        cols = point_colour(t) if point_colour else ["black"] * len(v)
        ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=16, c=cols, zorder=4, edgecolors="none")
    if ref_line is not None:
        ax.axhline(ref_line, color="grey", lw=0.8, ls="--")
    ax.set_xticks(range(len(order)), order, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=200)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- per n

def per_n(df: pd.DataFrame, metric: str, label: str, n: str, out_dir: Path, settings: dict) -> pd.DataFrame:
    a = settings["analysis"]
    sub = df[df["n"] == n]
    order = treatment_order(sub["treatment"], a["treatment_order"])
    groups = {t: sub.loc[sub["treatment"] == t, "value"].dropna().tolist() for t in order}
    rows = []
    for t in order:
        v = np.asarray(groups[t])
        rows.append(dict(n=n, treatment=t, fields=len(v), mean=v.mean() if len(v) else np.nan,
                         sd=v.std(ddof=1) if len(v) > 1 else np.nan,
                         sem=v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan,
                         median=np.median(v) if len(v) else np.nan))
    summ = pd.DataFrame(rows)
    for ref in a["normalise_to"]:
        if ref in groups and groups[ref]:
            summ[f"pct_of_{ref}"] = 100 * summ["mean"] / np.mean(groups[ref])
    summ.to_csv(out_dir / f"{metric}_summary.csv", index=False)
    prism_table(groups).to_csv(out_dir / f"{metric}_prism.csv", index=False)

    st = []
    valid = {t: v for t, v in groups.items() if len(v) >= 2}
    ref = a["reference"]
    if len(valid) >= 2:
        f, p = stats.f_oneway(*valid.values())
        st.append(dict(test="one-way ANOVA (fields as replicates, descriptive)", comparison="all", statistic=f, p=p, p_adj=p))
        if ref in valid and len(valid) > 1:
            others = [t for t in valid if t != ref]
            res = stats.dunnett(*[valid[t] for t in others], control=valid[ref])
            for t, pv, sv in zip(others, res.pvalue, res.statistic):
                st.append(dict(test="Dunnett vs reference (fields, descriptive)", comparison=f"{t} vs {ref}",
                               statistic=sv, p=pv, p_adj=pv))
    st = pd.DataFrame(st, columns=["test", "comparison", "statistic", "p", "p_adj"])
    st["significance"] = [stars(p) for p in st["p_adj"]]
    st.to_csv(out_dir / f"{metric}_stats.csv", index=False)
    _bar_plot(order, groups, "SD", label, f"{n}: {label} per field (mean ± SD)", out_dir / metric)
    return st.assign(n=n)


# ---------------------------------------------------------------- all n

def _rm_anova(mat: np.ndarray) -> tuple[float, float]:
    n, k = mat.shape
    gm = mat.mean()
    ss_t = n * ((mat.mean(0) - gm) ** 2).sum()
    ss_s = k * ((mat.mean(1) - gm) ** 2).sum()
    ss_e = ((mat - gm) ** 2).sum() - ss_t - ss_s
    df_t, df_e = k - 1, (n - 1) * (k - 1)
    if df_e <= 0 or ss_e <= 0:
        return np.nan, np.nan
    F = (ss_t / df_t) / (ss_e / df_e)
    return F, stats.f.sf(F, df_t, df_e)


def all_n(df: pd.DataFrame, metric: str, label: str, out_dir: Path, settings: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    a = settings["analysis"]
    ns = sorted(df["n"].unique())
    order = treatment_order(df["treatment"], a["treatment_order"])
    means = df.groupby(["n", "treatment"])["value"].mean().unstack("treatment").reindex(index=ns, columns=order)
    ncol = {n: N_COLOURS[i % len(N_COLOURS)] for i, n in enumerate(ns)}

    summ = means.T.copy()
    summ.columns = [f"{c}" for c in summ.columns]
    summ["N"] = means.notna().sum()
    summ["mean"] = means.mean()
    summ["sd"] = means.std(ddof=1)
    summ["sem"] = summ["sd"] / np.sqrt(summ["N"])
    norms = {}
    for ref in a["normalise_to"]:
        if ref in means:
            pct = means.div(means[ref], axis=0) * 100
            norms[ref] = pct
            summ[f"pct_of_{ref}_mean"] = pct.mean()
            summ[f"pct_of_{ref}_sem"] = pct.std(ddof=1) / np.sqrt(pct.notna().sum())
    summ.index.name = "treatment"
    summ.reset_index().to_csv(out_dir / f"{metric}_summary.csv", index=False)

    def to_prism(m: pd.DataFrame) -> pd.DataFrame:
        t = m.copy()
        t.insert(0, "", t.index)
        return t.reset_index(drop=True)

    to_prism(means).to_csv(out_dir / f"{metric}_prism.csv", index=False)
    for ref, pct in norms.items():
        to_prism(pct).to_csv(out_dir / f"{metric}_prism_pct_of_{ref}.csv", index=False)

    # statistics
    st = []
    complete = means.dropna(axis=0, how="any")
    if len(complete) >= 2 and means.shape[1] >= 2:
        F, p = _rm_anova(complete.to_numpy(float))
        st.append(dict(test=f"RM one-way ANOVA (n's with every treatment: {len(complete)})", comparison="all",
                       N=len(complete), effect=np.nan, ci_low=np.nan, ci_high=np.nan, statistic=F, p=p))
    rows = []
    refs = list(dict.fromkeys([a["reference"]] + list(a["normalise_to"])))
    for ref in refs:
        if ref not in means:
            continue
        fam = []
        for t in order:
            if t == ref:
                continue
            pair = means[[t, ref]].dropna()
            pair = pair[(pair > 0).all(axis=1)]
            r = dict(test=f"ratio paired t-test vs {ref} (Holm within this family)", comparison=f"{t} vs {ref}",
                     N=len(pair), effect=np.nan, ci_low=np.nan, ci_high=np.nan, statistic=np.nan, p=np.nan)
            if len(pair) >= 2:
                lr = np.log(pair[t] / pair[ref])
                tt = stats.ttest_1samp(lr, 0)
                ci = tt.confidence_interval()
                r.update(effect=float(np.exp(lr.mean())), ci_low=float(np.exp(ci.low)), ci_high=float(np.exp(ci.high)),
                         statistic=tt.statistic, p=tt.pvalue)
            fam.append(r)
        for r, pa in zip(fam, holm([r["p"] for r in fam])):
            r["p_adj"] = pa
        rows += fam
    for r in st:
        r["p_adj"] = r["p"]
    st = pd.DataFrame(st + rows, columns=["test", "comparison", "N", "effect", "ci_low", "ci_high", "statistic", "p", "p_adj"])
    st = st.rename(columns={"effect": "ratio"})
    st["significance"] = [stars(p) for p in st["p_adj"]]
    st.to_csv(out_dir / f"{metric}_stats.csv", index=False)

    vals = {t: means[t].tolist() for t in order}
    colour = lambda t: [ncol[n] for n in means.index if not np.isnan(means.at[n, t])]
    _bar_plot(order, vals, "SEM", label, f"All n: {label} (mean of n ± SEM; dots = n)", out_dir / metric, colour)
    for ref_, pct in norms.items():
        pv = {t: pct[t].tolist() for t in order}
        colour_p = lambda t, pct=pct: [ncol[n] for n in pct.index if not np.isnan(pct.at[n, t])]
        _bar_plot(order, pv, "SEM", f"{label} (% of {ref_})", f"All n: % of {ref_} (mean ± SEM)",
                  out_dir / f"{metric}_pct_of_{ref_}", colour_p, ref_line=100)
    legend = pd.DataFrame({"n": ns, "colour": [ncol[n] for n in ns]})
    legend.to_csv(out_dir / "n_colours.csv", index=False)
    return summ.reset_index(), st
