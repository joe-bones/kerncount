"""Summaries, statistics, Prism tables and plots for one metric at three levels:
per treatment (fields), per n (treatments compared within one experiment) and all n combined.

Statistics (every control group in `reference` + `normalise_to` gets its own test family; default Control, then LyoVec)
    within one n:  fields are technical replicates -> one-way ANOVA + Dunnett vs each control (descriptive only)
    across n:      each n contributes one value per treatment (mean of its fields)
                   - Prism standard (drives the stars): repeated-measures one-way ANOVA (n = matched replicate) and
                     Dunnett vs each control with the pooled residual error; least-squares fit, so an n missing a
                     group still contributes its other groups
                   - check: ratio paired t-test of each treatment vs each control (log ratio), Holm-adjusted per
                     control; one test per pair, no pooling (more conservative)
    two-way (extra, `two_way_anova`): fields as replicates, Type II sums of squares (handles missing groups)
                   - within one n: treatment x well        - across n: treatment x n
                   - Tukey(-Kramer) all-pairs comparison of treatments, matched within the second factor
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


def control_groups(a: dict) -> list[str]:
    """Control groups in priority order: `reference` first, then `normalise_to`."""
    return list(dict.fromkeys([a["reference"]] + list(a["normalise_to"])))


# ---------------------------------------------------------------- two-way ANOVA

def _design(d: pd.DataFrame, fa: str, fb: str, la: list, lb: list, terms: str) -> np.ndarray:
    A = np.column_stack([(d[fa] == x).to_numpy(float) for x in la[1:]]) if len(la) > 1 else np.empty((len(d), 0))
    B = np.column_stack([(d[fb] == x).to_numpy(float) for x in lb[1:]]) if len(lb) > 1 else np.empty((len(d), 0))
    cols = [np.ones((len(d), 1))]
    if "A" in terms:
        cols.append(A)
    if "B" in terms:
        cols.append(B)
    if "X" in terms and A.shape[1] and B.shape[1]:
        AB = np.column_stack([A[:, i] * B[:, j] for i in range(A.shape[1]) for j in range(B.shape[1])])
        cols.append(AB[:, AB.any(axis=0)])  # drop empty cells (missing groups)
    return np.hstack(cols)


def _fit(X: np.ndarray, y: np.ndarray) -> tuple[float, int]:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(((y - X @ beta) ** 2).sum()), int(np.linalg.matrix_rank(X))


def two_way(df: pd.DataFrame, factor_b: str, order: list[str], factor_a: str = "treatment"):
    """Two-way ANOVA (value ~ A + B + A:B, replicates = rows) with Type II sums of squares, and Tukey(-Kramer)
    all-pairs comparisons of the levels of A. Each pair is compared on the cell means of the B levels where both
    are present (matched), using the model's residual mean square, so a group missing in one B level does not bias
    the others. Without replication the interaction is not estimable and the additive model is used.
    Returns (anova table, pairwise table, cell means) or None if there is too little data."""
    d = df[[factor_a, factor_b, "value"]].dropna()
    la = [t for t in order if t in set(d[factor_a])]
    lb = sorted(d[factor_b].unique(), key=str)
    if len(la) < 2 or len(d) < 3:
        return None
    y = d["value"].to_numpy(float)
    fit = {t: _fit(_design(d, factor_a, factor_b, la, lb, t), y) for t in ("1", "A", "B", "AB", "ABX")}
    N = len(y)
    full = "ABX" if N - fit["ABX"][1] > 0 and fit["ABX"][1] > fit["AB"][1] else "AB"
    rss_res, rank_res = fit[full]
    df_res = N - rank_res
    if df_res <= 0:
        return None
    mse = rss_res / df_res
    ss_tot = float(((y - y.mean()) ** 2).sum())
    rows = []
    terms = [(factor_a, fit["B"], fit["AB"]), (factor_b, fit["A"], fit["AB"])]
    if full == "ABX":
        terms.append((f"{factor_a} x {factor_b} (interaction)", fit["AB"], fit["ABX"]))
    for name, (rss0, r0), (rss1, r1) in terms:
        dfx = r1 - r0
        if dfx <= 0:
            continue
        ss = max(rss0 - rss1, 0.0)
        F = (ss / dfx) / mse if mse > 0 else np.nan
        rows.append(dict(source=name, SS=ss, df=dfx, MS=ss / dfx, F=F, p=stats.f.sf(F, dfx, df_res) if mse > 0 else np.nan,
                         pct_of_total_variation=100 * ss / ss_tot if ss_tot > 0 else np.nan))
    rows.append(dict(source="residual", SS=rss_res, df=df_res, MS=mse, F=np.nan, p=np.nan,
                     pct_of_total_variation=100 * rss_res / ss_tot if ss_tot > 0 else np.nan))
    anova = pd.DataFrame(rows)
    anova["significance"] = [stars(p) for p in anova["p"]]
    anova.insert(0, "model", "interaction" if full == "ABX" else "additive (no replication: interaction not estimable)")

    g = d.groupby([factor_a, factor_b])["value"]
    cm, cn = g.mean().unstack(factor_b), g.count().unstack(factor_b)
    k = len(la)
    q_crit = stats.studentized_range.ppf(0.95, k, df_res) if k >= 2 else np.nan
    pairs = []
    for i in range(k):
        for j in range(i + 1, k):
            ti, tj = la[i], la[j]
            common = [b for b in lb if not (np.isnan(cm.at[ti, b]) or np.isnan(cm.at[tj, b]))]
            r = dict(comparison=f"{ti} vs {tj}", group_1=ti, group_2=tj, matched_levels=len(common),
                     mean_diff=np.nan, ci_low=np.nan, ci_high=np.nan, q=np.nan, p_adj=np.nan)
            if common:
                diff = float(np.mean([cm.at[ti, b] - cm.at[tj, b] for b in common]))
                se = float(np.sqrt(mse * sum(1 / cn.at[ti, b] + 1 / cn.at[tj, b] for b in common)) / len(common))
                if se > 0:
                    q = np.sqrt(2) * abs(diff) / se
                    half = q_crit * se / np.sqrt(2)
                    r.update(mean_diff=diff, ci_low=diff - half, ci_high=diff + half, q=q,
                             p_adj=float(stats.studentized_range.sf(q, k, df_res)))
            pairs.append(r)
    pairs = pd.DataFrame(pairs)
    pairs["significance"] = [stars(p) for p in pairs["p_adj"]]
    pairs.insert(0, "test", f"Tukey all pairs ({factor_a}, matched within {factor_b}; residual df {df_res})")
    return anova, pairs, cm.reindex(index=la)


def _two_way_plot(df: pd.DataFrame, factor_b: str, order: list[str], res, label: str, title: str, out: Path,
                  colours: dict | None = None, b_label: str | None = None, note: str = ""):
    anova, pairs, cm = res
    la = list(cm.index)
    lb = list(cm.columns)
    b_label = b_label or factor_b
    colours = colours or {b: N_COLOURS[i % len(N_COLOURS)] for i, b in enumerate(lb)}
    k = len(la)
    fig, (ax, axm) = plt.subplots(1, 2, figsize=(max(10, 1.1 * k + 6), max(4.8, 0.35 * k + 2.5)),
                                  gridspec_kw=dict(width_ratios=[1.5, 1]))
    # left: grouped bars, one bar per (treatment, second factor), mean ± SEM of fields
    w = 0.8 / max(len(lb), 1)
    g = df.dropna(subset=["value"]).groupby(["treatment", factor_b])["value"]
    mean, sem = g.mean(), g.sem()
    for j, b in enumerate(lb):
        xs, ms, es = [], [], []
        for i, t in enumerate(la):
            if (t, b) in mean.index:
                xs.append(i - 0.4 + w * (j + 0.5))
                ms.append(mean[(t, b)])
                es.append(0 if np.isnan(sem[(t, b)]) else sem[(t, b)])
        ax.bar(xs, ms, width=w, color=colours.get(b, "#999999"), edgecolor="black", lw=0.4, label=str(b), zorder=1)
        ax.errorbar(xs, ms, yerr=es, fmt="none", ecolor="black", capsize=2, lw=0.7, zorder=2)
    ax.set_xticks(range(k), la, rotation=45, ha="right")
    ax.set_ylabel(label)
    ax.set_title(f"{title}\nbars = mean ± SEM of fields per {b_label}", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(title=b_label, fontsize=7, title_fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    # right: Tukey all-pairs matrix (lower triangle), colour = -log10 p, text = stars
    mat = np.full((k, k), np.nan)
    txt = [["" for _ in range(k)] for _ in range(k)]
    for _, r in pairs.iterrows():
        i, j = la.index(r["group_1"]), la.index(r["group_2"])
        if not np.isnan(r["p_adj"]):
            mat[j, i] = -np.log10(max(r["p_adj"], 1e-6))
            txt[j][i] = r["significance"]
    im = axm.imshow(mat, cmap="Blues", vmin=0, vmax=4)
    for i in range(k):
        for j in range(k):
            if txt[i][j]:
                axm.text(j, i, txt[i][j], ha="center", va="center", fontsize=7,
                         color="white" if mat[i, j] > 2.5 else "black")
    axm.set_xticks(range(k), la, rotation=45, ha="right", fontsize=7)
    axm.set_yticks(range(k), la, fontsize=7)
    axm.set_title("Tukey all pairs (p adjusted)", fontsize=9)
    for sp in axm.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=axm, fraction=0.046, pad=0.04)
    cb.set_label("-log10 p", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    parts = [f"{r['source']}: F({int(r['df'])},{int(anova['df'].iloc[-1])})={r['F']:.2f}, p={r['p']:.2g} ({r['significance']})"
             for _, r in anova.iterrows() if r["source"] != "residual"]
    foot = ("Two-way ANOVA (Type II SS; fields = replicates): " + "; ".join(parts) +
            ". * p<0.05  ** p<0.01  *** p<0.001  **** p<0.0001; ns = not significant." + note)
    fig.text(0.01, 0.005, foot, fontsize=6.5, color="#444444", ha="left", va="bottom", wrap=True)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(out.with_suffix(".png"), dpi=200)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def _write_two_way(df: pd.DataFrame, factor_b: str, order: list[str], metric: str, label: str, title: str,
                   out_dir: Path, colours: dict | None = None, b_label: str | None = None, note: str = "") -> None:
    res = two_way(df, factor_b, order)
    if res is None:
        return
    anova, pairs, _ = res
    anova.to_csv(out_dir / f"{metric}_twoway_anova.csv", index=False)
    pairs.to_csv(out_dir / f"{metric}_twoway_pairs.csv", index=False)
    _two_way_plot(df, factor_b, order, res, label, title, out_dir / f"{metric}_twoway", colours, b_label, note)


# ---------------------------------------------------------------- plotting

def significance_marks(st: pd.DataFrame, ref: str, test_contains: str) -> dict[str, str]:
    """{treatment: stars} for significant comparisons '<t> vs <ref>' of one test family."""
    out = {}
    if st is None or st.empty:
        return out
    fam = st[st["test"].str.contains(test_contains, regex=False) & st["comparison"].str.endswith(f" vs {ref}")]
    for _, r in fam.iterrows():
        mark = stars(r["p_adj"])
        if mark and mark != "ns":
            out[r["comparison"][: -len(f" vs {ref}")]] = mark
    return out


def _bar_plot(order, values: dict, err: str, ylabel: str, title: str, out: Path, point_colour=None, ref_line=None,
              sig: dict | None = None, ref: str | None = None, footer: str = "", legend: dict | None = None):
    fig, ax = plt.subplots(figsize=(max(4.5, 0.8 * len(order) + 1.5), 4.4))
    rng = np.random.default_rng(0)
    tops = {}
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
        tops[t] = max(v.max(), m + e)
    if ref_line is not None:
        ax.axhline(ref_line, color="grey", lw=0.8, ls="--")
    lo, hi = ax.get_ylim()
    pad = 0.04 * (hi - lo)
    for i, t in enumerate(order):
        if t not in tops:
            continue
        if sig and t in sig:
            ax.text(i, tops[t] + pad, sig[t], ha="center", va="bottom", fontsize=11, fontweight="bold")
        elif ref and t == ref and sig is not None:
            ax.text(i, tops[t] + pad, "ref", ha="center", va="bottom", fontsize=8, color="grey")
    ax.set_ylim(lo, hi + 4 * pad)
    ax.set_xticks(range(len(order)), order, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    if legend:
        from matplotlib.lines import Line2D
        handles = [Line2D([], [], marker="o", ls="", color=c, markersize=5, label=k) for k, c in legend.items()]
        ax.legend(handles=handles, fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    if footer:
        fig.text(0.01, 0.005, footer, fontsize=6.5, color="#444444", ha="left", va="bottom", wrap=True)
    fig.tight_layout(rect=(0, 0.06 if footer else 0, 1, 1))
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
    refs = [c for c in control_groups(a) if c in valid]  # main reference first; missing controls are skipped
    if len(valid) >= 2:
        f, p = stats.f_oneway(*valid.values())
        st.append(dict(test="one-way ANOVA (fields as replicates, descriptive)", comparison="all", statistic=f, p=p, p_adj=p))
        for ref in refs:
            others = [t for t in valid if t != ref]
            res = stats.dunnett(*[valid[t] for t in others], control=valid[ref])
            for t, pv, sv in zip(others, res.pvalue, res.statistic):
                st.append(dict(test=f"Dunnett vs {ref} (fields, descriptive)", comparison=f"{t} vs {ref}",
                               statistic=sv, p=pv, p_adj=pv))
    st = pd.DataFrame(st, columns=["test", "comparison", "statistic", "p", "p_adj"])
    st["significance"] = [stars(p) for p in st["p_adj"]]
    st.to_csv(out_dir / f"{metric}_stats.csv", index=False)
    title = f"{n}: {label} (mean ± SD; dots = fields)"
    if not refs or len(valid) < 2:
        _bar_plot(order, groups, "SD", label, title, out_dir / metric, ref=a["reference"],
                  footer="No comparison possible: no control group with >= 2 fields in this n.")
    for i, ref in enumerate(refs if len(valid) >= 2 else []):
        foot = (f"Stars: Dunnett vs {ref} (fields = technical replicates; describes this experiment only). "
                "* p<0.05  ** p<0.01  *** p<0.001  **** p<0.0001; no mark = not significant.")
        if i == 0 and ref != a["reference"]:
            foot = f"{a['reference']} missing in this n; compared with {ref}. " + foot
        # the first available control gets the plain file name, every further control its own plot
        stem = metric if i == 0 else f"{metric}_vs_{ref}"
        _bar_plot(order, groups, "SD", label, title, out_dir / stem,
                  sig=significance_marks(st, ref, f"Dunnett vs {ref} ("), ref=ref, footer=foot)
    if a.get("two_way_anova", True) and "image" in sub:
        from .sampling import design, well_of  # local import: sampling imports this module
        W, K, layout, _ = design(settings)
        if W >= 2:
            sw = sub.assign(well=[f"well {well_of(int(i), K, W, layout)}" for i in sub["image"]])
            _write_two_way(sw, "well", order, metric, label, f"{n}: {label}, two-way ANOVA (treatment x well)",
                           out_dir, b_label="well", note=" Describes this experiment only.")
    return st.assign(n=n)


# ---------------------------------------------------------------- all n

def _dunnett_padj(t: np.ndarray, R: np.ndarray, df: int) -> tuple[np.ndarray, float]:
    """Two-sided Dunnett adjusted p-values for statistics t with correlation R (multivariate t), and the 95% critical value."""
    mvt = stats.multivariate_t(shape=R, df=df, allow_singular=True)
    m = len(t)

    def cover(c):
        return float(mvt.cdf(np.full(m, c), lower_limit=np.full(m, -c), random_state=0))

    if m == 1:
        padj = 2 * stats.t.sf(np.abs(t), df)
        return padj, float(stats.t.ppf(0.975, df))
    padj = np.array([min(1.0, max(0.0, 1 - cover(abs(x)))) for x in t])
    from scipy.optimize import brentq
    crit = brentq(lambda c: cover(c) - 0.95, stats.t.ppf(0.975, df), stats.t.ppf(1 - 0.025 / m, df) + 1)
    return padj, crit


def rm_anova_dunnett(means: pd.DataFrame, order: list[str], ref: str) -> tuple[dict | None, list[dict]]:
    """Prism standard for matched data: repeated-measures one-way ANOVA on the per-n means (n = subject), then
    Dunnett's multiple comparisons vs `ref` using the pooled residual error. Fitted as value ~ treatment + n by least
    squares, so an n missing a group still contributes all its other groups (equal to the RM ANOVA when nothing is
    missing; Prism would switch to its mixed-effects model in that case)."""
    long = means.reindex(columns=[t for t in order if t in means]).stack().rename("value").dropna().reset_index()
    long.columns = ["n", "treatment", "value"]
    la = [ref] + [t for t in order if t in set(long["treatment"]) and t != ref]
    lb = sorted(long["n"].unique(), key=str)
    if ref not in set(long["treatment"]) or len(la) < 2 or len(lb) < 2:
        return None, []
    X = _design(long, "treatment", "n", la, lb, "AB")
    y = long["value"].to_numpy(float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    rank = int(np.linalg.matrix_rank(X))
    dfe = len(y) - rank
    if dfe <= 0:
        return None, []
    mse = float(((y - X @ beta) ** 2).sum()) / dfe
    rss_b, r_b = _fit(_design(long, "treatment", "n", la, lb, "B"), y)
    df_t = rank - r_b
    F = ((rss_b - mse * dfe) / df_t) / mse if mse > 0 and df_t > 0 else np.nan
    anova = dict(test=f"RM one-way ANOVA (n = matched replicate; {len(lb)} n's, residual df {dfe})", comparison="all",
                 N=len(lb), effect=np.nan, ci_low=np.nan, ci_high=np.nan, statistic=F,
                 p=stats.f.sf(F, df_t, dfe) if not np.isnan(F) else np.nan)
    cov = mse * np.linalg.pinv(X.T @ X)
    idx = list(range(1, len(la)))  # treatment coefficients = difference vs ref, adjusted for n
    est, se = beta[idx], np.sqrt(np.diag(cov)[idx])
    ok = se > 0
    rows = []
    if ok.any():
        tt = est[ok] / se[ok]
        sub = cov[np.ix_(np.array(idx)[ok], np.array(idx)[ok])]
        R = sub / np.outer(np.sqrt(np.diag(sub)), np.sqrt(np.diag(sub)))
        padj, crit = _dunnett_padj(tt, R, dfe)
    k = 0
    for i, t in enumerate(la[1:]):
        both = int(means[[t, ref]].dropna().shape[0])
        r = dict(test=f"Dunnett vs {ref} (RM ANOVA, pooled error; Prism standard)", comparison=f"{t} vs {ref}",
                 N=both, effect=np.nan, ci_low=np.nan, ci_high=np.nan, statistic=np.nan, p=np.nan, p_adj=np.nan)
        if ok[i]:
            r.update(effect=est[i], ci_low=est[i] - crit * se[i], ci_high=est[i] + crit * se[i], statistic=tt[k],
                     p=2 * stats.t.sf(abs(tt[k]), dfe), p_adj=padj[k])
            k += 1
        rows.append(r)
    return anova, rows


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

    # statistics: Prism standard (RM one-way ANOVA + Dunnett vs each control) drives the stars;
    # ratio paired t-tests (one test per pair, no pooling) are kept as a check
    st = []
    rows = []
    refs = control_groups(a)
    for ref in refs:
        if ref not in means:
            continue
        anova, dun = rm_anova_dunnett(means, order, ref)
        if anova is not None and not st:
            st.append(anova)
        rows += dun
    for ref in refs:
        if ref not in means:
            continue
        fam = []
        for t in order:
            if t == ref:
                continue
            pair = means[[t, ref]].dropna()
            pair = pair[(pair > 0).all(axis=1)]
            r = dict(test=f"ratio paired t-test vs {ref} (Holm within this family; check)", comparison=f"{t} vs {ref}",
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
    st.insert(3, "effect_type", ["mean difference vs control" if "Dunnett" in x else "ratio vs control (geometric mean)"
                                 if "ratio" in x else "" for x in st["test"]])
    st["significance"] = [stars(p) for p in st["p_adj"]]
    st.to_csv(out_dir / f"{metric}_stats.csv", index=False)

    def foot_for(ref_):
        anova = st[st["test"].str.startswith("RM one-way ANOVA")]
        a_txt = (f"RM one-way ANOVA p={anova['p'].iloc[0]:.2g} ({stars(anova['p'].iloc[0])}). "
                 if len(anova) and not np.isnan(anova["p"].iloc[0]) else "")
        return (f"{a_txt}Stars: Dunnett vs {ref_} after RM one-way ANOVA (n = matched replicate, pooled error). "
                "* p<0.05  ** p<0.01  *** p<0.001  **** p<0.0001; no mark = not significant.")

    vals = {t: means[t].tolist() for t in order}
    colour = lambda t: [ncol[n] for n in means.index if not np.isnan(means.at[n, t])]
    main_ref = next((r_ for r_ in refs if r_ in means), None)
    _bar_plot(order, vals, "SEM", label, f"All n: {label} (mean of n ± SEM; dots = n)", out_dir / metric, colour, legend=ncol,
              sig=significance_marks(st, main_ref, f"Dunnett vs {main_ref} (") if main_ref else None, ref=main_ref,
              footer=foot_for(main_ref) if main_ref else "")
    for ref_, pct in norms.items():
        pv = {t: pct[t].tolist() for t in order}
        colour_p = lambda t, pct=pct: [ncol[n] for n in pct.index if not np.isnan(pct.at[n, t])]
        _bar_plot(order, pv, "SEM", f"{label} (% of {ref_})", f"All n: % of {ref_} (mean ± SEM)",
                  out_dir / f"{metric}_pct_of_{ref_}", colour_p, ref_line=100, legend=ncol,
                  sig=significance_marks(st, ref_, f"Dunnett vs {ref_} ("), ref=ref_, footer=foot_for(ref_))
    if a.get("two_way_anova", True) and len(ns) >= 2:
        _write_two_way(df, "n", order, metric, label, f"All n: {label}, two-way ANOVA (treatment x n)", out_dir,
                       colours=ncol, b_label="n",
                       note=" Fields are counted as replicates here, so these p-values are optimistic; the RM ANOVA +"
                            " Dunnett across n (one value per n) remains the main result.")
    legend = pd.DataFrame({"n": ns, "colour": [ncol[n] for n in ns]})
    legend.to_csv(out_dir / "n_colours.csv", index=False)
    return summ.reset_index(), st
