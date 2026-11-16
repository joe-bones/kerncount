"""Plots for the hierarchical model, comparison families, robustness and the plate design.

Same style as analysis.py: white background, no top/right spines, grey footers explaining the statistics,
dot colour = n (analysis.N_COLOURS), PNG at 200 dpi plus a PDF.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

from .analysis import N_COLOURS

FOOT = dict(fontsize=6.5, color="#444444", ha="left", va="bottom", wrap=True)
WELL_COLOURS = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3", "#937860"]
POS_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">", "p"]


def _save(fig, out: Path, footer: str = "", bottom: float = 0.0):
    if footer:
        fig.text(0.01, 0.005, footer, **FOOT)
    fig.tight_layout(rect=(0, bottom if footer else 0, 1, 1))
    fig.savefig(out.with_suffix(".png"), dpi=200)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def _plain(v, _=None) -> str:
    if v <= 0:
        return ""
    if v >= 100:
        return f"{v:.0f}"
    if v >= 1:
        return f"{v:g}"
    return f"{v:.2g}"


def _log_axis(ax, axis="y"):
    a = ax.yaxis if axis == "y" else ax.xaxis
    (ax.set_yscale if axis == "y" else ax.set_xscale)("log")
    a.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 3.0, 5.0)))
    a.set_major_formatter(FuncFormatter(_plain))
    a.set_minor_formatter(NullFormatter())


def _clean(ax):
    ax.spines[["top", "right"]].set_visible(False)


def n_colours(ns) -> dict:
    return {n: N_COLOURS[i % len(N_COLOURS)] for i, n in enumerate(sorted(ns))}


def short_n(n: str) -> str:
    return str(n).split("_")[0]


# ---------------------------------------------------------------- all n: model

def model_absolute(res, out: Path, ncol: dict, sig: dict, ref: str | None, footer: str, title: str):
    """Model means (absolute units: mean per field, or % of nuclei for fractions) ± 95% CI; per-n observed values as
    coloured dots; faint per-field points. With one n the bar is the estimate from that n and there is no CI."""
    p, e = res.prepared, res.emm
    fl = p.fields if p.fields is not None else p.d
    fl = fl[fl["value"].notna()]
    order = [t for t in _order_hint(res) if t in set(p.d["treatment"])]  # user order (the baseline is not shown first)
    frac = p.scale.kind == "fraction"
    k = 100.0 if frac else 1.0
    fig, ax = plt.subplots(figsize=(max(5.5, 0.85 * len(order) + 2.5), 4.8))
    rng = np.random.default_rng(0)
    from .models import _observed
    obs = _observed(p)
    pos = fl.loc[fl["value"] > 0, "value"]
    lows = e["ci_low"].dropna()
    floor = float(min(pos.min() if len(pos) else 1.0, lows.min() if len(lows) else 1.0)) / 1.5
    cap = np.inf
    if frac:  # keep the estimates readable: single extreme fields are drawn at the top edge
        tops = [float(np.nanmax(e[["estimate", "ci_high"]].to_numpy(float))) if len(e) else 0.0,
                float(obs.max()) if len(obs) else 0.0,
                float(np.nanpercentile(fl["value"], 97)) if len(fl) else 0.0]
        cap = max(tops) * k * 1.15 or 1.0
    for i, t in enumerate(order):
        f = fl[fl["treatment"] == t]
        fv = f["value"].to_numpy(float) * k
        ax.scatter(i + rng.uniform(-0.28, 0.28, len(fv)), np.clip(fv, None, cap) if frac else np.clip(fv, floor, None), s=5,
                   color="#bbbbbb", zorder=1, edgecolors="none")
        ns = sorted(set(p.d.loc[p.d["treatment"] == t, "n"]))
        xs = i + np.linspace(-0.12, 0.12, max(len(ns), 1)) if len(ns) > 1 else [i]
        for x, n in zip(xs, ns):
            if (n, t) in obs.index:
                ax.scatter(x, float(obs[(n, t)]) * k, s=26, color=ncol.get(n, "black"), zorder=3,
                           edgecolors="white", linewidths=0.6)
        r = e[e["treatment"] == t]
        if len(r):
            r = r.iloc[0]
            ax.hlines(r["estimate"] * k, i - 0.3, i + 0.3, color="black", lw=2, zorder=4)
            if np.isfinite(r["ci_low"]):
                ax.vlines(i, r["ci_low"] * k, r["ci_high"] * k, color="black", lw=1.2, zorder=4)
                for yv in (r["ci_low"], r["ci_high"]):
                    ax.hlines(yv * k, i - 0.08, i + 0.08, color="black", lw=1.2, zorder=4)
    if not frac:
        _log_axis(ax)
    else:
        ax.set_ylim(bottom=0)
    lo, hi = ax.get_ylim()
    for i, t in enumerate(order):
        r = e[e["treatment"] == t]
        tops = [float(r["ci_high"].iloc[0]) * k if len(r) and np.isfinite(r["ci_high"].iloc[0]) else lo]
        tops += [float(obs[(n, t)]) * k for n in set(p.d.loc[p.d["treatment"] == t, "n"]) if (n, t) in obs.index]
        if len(r):
            tops.append(float(r["estimate"].iloc[0]) * k)
        top = max(tops)
        y = top * 1.12 if not frac else top + 0.04 * (hi - lo)
        txt = sig.get(t) if sig else None
        if txt:
            ax.text(i, y, txt, ha="center", va="bottom", fontsize=11, fontweight="bold")
        elif t == ref and sig is not None:
            ax.text(i, y, "ref", ha="center", va="bottom", fontsize=8, color="grey")
    if not frac:
        ax.set_ylim(max(lo, floor / 1.2), hi * 1.9)
    else:
        ax.set_ylim(0, min(hi, cap) * 1.12)
    ax.set_xticks(range(len(order)), order, rotation=45, ha="right")
    if frac:
        ylab = f"{res.metric.label} (% of nuclei)"
        kind = "% of nuclei (pooled per well)" if p.unit == "well" else "% of nuclei"
    else:
        unit = res.metric.unit or res.metric.label
        ylab = (f"{res.metric.label}" + (f" ({unit})" if unit and unit.lower() not in res.metric.label.lower() else "")
                + " (log scale)")
        kind = "mean per field" if p.unit == "well" else "geometric mean of fields"
    ax.set_ylabel(ylab, fontsize=9)
    bar = (f"bar = model {kind} ± 95% CI" if res.ok else f"bar = estimate from the data ({kind}; CI {NEEDS_2N_TXT})")
    ax.set_title(f"{title}\n{bar}; dots = per-n values; grey = fields", fontsize=9)
    _clean(ax)
    handles = [Line2D([], [], marker="o", ls="", color=c, markersize=5, label=short_n(n)) for n, c in ncol.items()
               if n in set(p.d["n"])]
    ax.legend(handles=handles, title="n", fontsize=7, title_fontsize=7, frameon=False, loc="upper left",
              bbox_to_anchor=(1.0, 1.0))
    _save(fig, out, footer, 0.11)


NEEDS_2N_TXT = "needs ≥ 2 n"


def _order_hint(res) -> list[str]:
    cfg_order = getattr(res, "order", None)
    return list(cfg_order) if cfg_order else list(res.prepared.treatments)


def forest(fam: pd.DataFrame, loo: pd.DataFrame, out: Path, ncol: dict, title: str, footer: str,
           estimate_type: str = "ratio", stars_col: str = "significance", margin: float | None = None,
           pp: bool = False):
    """One row per comparison: estimate with 95% CI (log axis, line at 1; pp=True: difference in percentage points
    on a linear axis, line at 0), leave-one-n-out estimates as small coloured marks (colour = the n that was
    dropped), stars from the primary test on the right. Rows without a CI show the estimate only ('needs ≥ 2 n')."""
    fam = fam.reset_index(drop=True)
    if pp:
        cols = dict(est="pp_difference", lo="pp_ci_low", hi="pp_ci_high", loa="pp_ci_low_adj", hia="pp_ci_high_adj")
    else:
        cols = dict(est="estimate", lo="ci_low", hi="ci_high", loa="ci_low_adj", hia="ci_high_adj")
    get = lambda r, c: float(getattr(r, cols[c], np.nan)) if cols[c] in fam else np.nan  # noqa: E731
    k = len(fam)
    fig, ax = plt.subplots(figsize=(8.2, max(2.4, 0.42 * k + 1.6)))
    ys = np.arange(k)[::-1]
    for y, r in zip(ys, fam.itertuples()):
        if np.isfinite(get(r, "lo")):
            ax.hlines(y, get(r, "lo"), get(r, "hi"), color="black", lw=1.4, zorder=2)
        if np.isfinite(get(r, "loa")):
            ax.hlines(y, get(r, "loa"), get(r, "hia"), color="black", lw=0.6, ls=(0, (2, 2)), zorder=1)
        ax.scatter(get(r, "est"), y, marker="D", s=34, color="black" if np.isfinite(get(r, "lo")) else "white",
                   edgecolors="black", zorder=4)
        if loo is not None and len(loo):
            dr = loo[(loo["family"] == r.family) & (loo["comparison"] == r.comparison)]
            lcol = "pp_difference" if pp else "estimate"
            for q in dr.itertuples():
                v = getattr(q, lcol, np.nan)
                if np.isfinite(v):
                    ax.scatter(v, y + 0.22, s=14, color=ncol.get(q.dropped_n, "grey"), zorder=3,
                               edgecolors="white", linewidths=0.4)
    null = 0.0 if pp else 1.0
    ax.axvline(null, color="grey", lw=0.8, ls="--", zorder=0)
    if margin:
        ax.axvspan(1 / (1 + margin), 1 + margin, color="#c7e9c0", alpha=0.6, zorder=0, lw=0)
        lo, hi = ax.get_xlim()
        ax.set_xlim(min(lo, 1 / (1 + margin) / 1.15), max(hi, (1 + margin) * 1.15))
    if not pp:
        _log_axis(ax, "x")
    ax.set_yticks(ys, [c.replace(" vs ", "  vs  ") for c in fam["comparison"]], fontsize=8)
    ax.set_ylim(-0.7, k - 0.2)
    if pp:
        ax.set_xlabel("difference (percentage points; 0 = no difference)", fontsize=8)
    else:
        lab = "odds ratio" if estimate_type == "odds_ratio" else "ratio"
        ax.set_xlabel(f"{lab} (log scale; 1 = no difference)", fontsize=8)
    ax.set_title(title, fontsize=9)
    _clean(ax)
    trans = ax.get_yaxis_transform()
    fmt = (lambda v: f"{v:+.1f}".replace("-0.0", "0.0").replace("+0.0", "0.0")) if pp else (lambda v: f"{v:.2f}")  # noqa: E731
    for y, r in zip(ys, fam.itertuples()):
        est, lo_, hi_ = get(r, "est"), get(r, "lo"), get(r, "hi")
        txt = f"{fmt(est)} [{fmt(lo_)}, {fmt(hi_)}]" if np.isfinite(lo_) else f"{fmt(est)} ({NEEDS_2N_TXT})"
        ax.text(1.02, y, txt, transform=trans, fontsize=7, va="center", ha="left", color="#333333")
        s = getattr(r, stars_col, "")
        if isinstance(s, str) and s:
            ax.text(1.27, y, s, transform=trans, fontsize=8 if len(s) > 4 else 9, va="center", ha="left",
                    fontweight="bold" if s != "ns" else "normal")
    unit = " pp" if pp else ""
    ax.text(1.02, k - 0.35, ("estimate [90% CI]" if margin else f"estimate{unit} [95% CI]"), transform=trans,
            fontsize=6.5, color="grey", ha="left")
    handles = [Line2D([], [], marker="D", ls="-", color="black", markersize=5,
                      label="all n (90% CI)" if margin else "all n (95% CI)")]
    if margin:
        handles.append(plt.Rectangle((0, 0), 1, 1, color="#c7e9c0", label=f"equivalence margin ±{margin:.0%}"))
    elif cols["loa"] in fam and fam[cols["loa"]].notna().any():
        handles.append(Line2D([], [], ls=(0, (2, 2)), color="black", lw=0.8, label="simultaneous 95% CI"))
    if not fam[cols["lo"]].notna().any():
        handles.append(Line2D([], [], marker="D", ls="", markerfacecolor="white", color="black", markersize=5,
                              label=f"estimate only ({NEEDS_2N_TXT})"))
    handles += [Line2D([], [], marker="o", ls="", color=c, markersize=4, label=f"without {short_n(n)}")
                for n, c in ncol.items() if loo is not None and len(loo) and n in set(loo["dropped_n"])]
    ax.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(1.33, 1.02))
    fig.subplots_adjust(right=0.62)
    if footer:
        fig.text(0.01, 0.005, footer, **FOOT)
    fig.tight_layout(rect=(0, 0.12 if footer else 0, 1, 1))
    fig.savefig(out.with_suffix(".png"), dpi=200)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def per_n_ratio(pn: pd.DataFrame, fam: pd.DataFrame | None, control: str, order: list[str], out: Path, ncol: dict,
                title: str, footer: str, col: str = "ratio", pp: bool = False):
    """Effect of each treatment vs `control` within each n (dots, coloured by n, joined across treatments per n;
    thin vertical lines = technical 95% interval from the fields), with the model estimate (black bar).
    col: the per-n column drawn (ratio, ratio_of_means or pp_difference); pp = linear axis with 0 as no effect."""
    sub = pn[pn["control"] == control]
    trs = [t for t in order if t in set(sub["treatment"])]
    if not trs:
        return False
    fig, ax = plt.subplots(figsize=(max(5.5, 0.8 * len(trs) + 2.5), 4.4))
    for n, g in sub.groupby("n"):
        g = g.set_index("treatment").reindex(trs)
        xs = np.arange(len(trs)) + (sorted(ncol).index(n) - (len(ncol) - 1) / 2) * 0.07 if n in ncol else np.arange(len(trs))
        if "tech_ci_low" in g and col != "ratio":
            ax.vlines(xs, g["tech_ci_low"], g["tech_ci_high"], color=ncol.get(n, "grey"), lw=0.7, alpha=0.5, zorder=1)
        ax.plot(xs, g[col], color=ncol.get(n, "grey"), lw=0.8, alpha=0.6, zorder=2)
        ax.scatter(xs, g[col], s=26, color=ncol.get(n, "grey"), zorder=3, edgecolors="white", linewidths=0.6,
                   label=short_n(n))
    if fam is not None and len(fam):
        mcol = "pp_difference" if pp else "estimate"
        for i, t in enumerate(trs):
            r = fam[fam["group"] == t]
            if len(r) and np.isfinite(r[mcol].iloc[0]):
                ax.hlines(r[mcol].iloc[0], i - 0.3, i + 0.3, color="black", lw=2, zorder=4)
    ax.axhline(0 if pp else 1, color="grey", lw=0.8, ls="--", zorder=0)
    if not pp:
        _log_axis(ax)
    ax.set_xticks(range(len(trs)), trs, rotation=45, ha="right")
    ax.set_ylabel(f"difference from {control} within each n (percentage points)" if pp else
                  f"ratio to {control} within each n (log scale)", fontsize=9)
    ax.set_title(title, fontsize=9)
    _clean(ax)
    ax.legend(title="n", fontsize=7, title_fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    _save(fig, out, footer, 0.09)
    return True


def variance_components(vc: pd.DataFrame, out: Path, title: str, footer: str):
    v = vc[vc["variance"].notna()].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.5, 0.45 * max(len(v), 1) + 1.8))
    if v.empty:
        ax.text(0.5, 0.5, f"No variance component estimable ({NEEDS_2N_TXT} for the n-level components).",
                ha="center", va="center", fontsize=9, color="#555555", transform=ax.transAxes)
        ax.set_axis_off()
        ax.set_title(title, fontsize=9)
        _save(fig, out, footer, 0.1)
        return True
    ys = np.arange(len(v))[::-1]
    ax.barh(ys, v["variance"], color="#9ecae1", edgecolor="black", lw=0.5, height=0.6, zorder=2)
    for y, r in zip(ys, v.itertuples()):
        resid = r.component.startswith("field (residual") or r.component.startswith("well (residual")
        if np.isfinite(r.ci_low) and np.isfinite(r.ci_high) and not resid:
            ax.hlines(y, r.ci_low, min(r.ci_high, v["variance"].max() * 3), color="black", lw=0.9, zorder=3)
        label = f"{r.pct_of_total:.0f}%" if np.isfinite(r.pct_of_total) else ""
        if np.isfinite(r.cv_pct):
            label += f"  (SD {r.sd:.2f} ≈ {r.cv_pct:.0f}% CV)"
        if "descriptive" in str(r.note):
            label += "  (descriptive)"
        right = max(r.variance, r.ci_high if np.isfinite(r.ci_high) and not resid else 0)
        ax.text(min(right, v["variance"].max() * 3), y, "  " + label, va="center", fontsize=7, color="#333333", zorder=4)
    ax.set_yticks(ys, v["component"], fontsize=8)
    ax.set_xlabel("variance on the analysis scale (log or logit)", fontsize=8)
    notres = ~(v["component"].str.startswith("field (residual") | v["component"].str.startswith("well (residual"))
    xmax = max(v["variance"].max(), np.nanmax(np.where(notres, v["ci_high"].fillna(0), 0)))
    ax.set_xlim(0, min(xmax, v["variance"].max() * 3) * 1.6)
    ax.set_title(title, fontsize=9)
    _clean(ax)
    _save(fig, out, footer, 0.1)
    return True


def diagnostics(res, out: Path, ncol: dict, title: str):
    """Residuals vs fitted, normal QQ of residuals, and the estimated random effects per n. Without a model fit
    (one n) the residuals are deviations from the condition means and the random-effects panel says why it is empty."""
    from scipy import stats

    f = res.fit
    d = res.prepared.d
    if f is not None:
        fitted, resid = f.fitted_resid(len(d))
    else:
        fitted = d.groupby(["n", "treatment"])["y"].transform("mean").to_numpy(float)
        resid = d["y"].to_numpy(float) - fitted
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), gridspec_kw=dict(width_ratios=[1, 1, 1.3]))
    ax = axs[0]
    for n in sorted(set(d["n"])):
        m = (d["n"] == n).to_numpy()
        ax.scatter(fitted[m], resid[m], s=9, color=ncol.get(n, "grey"), alpha=0.7, edgecolors="none", label=short_n(n))
    ax.axhline(0, color="grey", lw=0.8, ls="--")
    ax.set_xlabel("fitted (analysis scale)", fontsize=8)
    ax.set_ylabel("residual", fontsize=8)
    ax.set_title("Residuals vs fitted (should be a shapeless band)", fontsize=8.5)
    ax.legend(fontsize=6.5, frameon=False, markerscale=1.5)
    _clean(ax)
    ax = axs[1]
    if len(resid) >= 3:
        (osm, osr), (slope, inter, _) = stats.probplot(resid, dist="norm")
        ax.scatter(osm, osr, s=9, color="#4c72b0", edgecolors="none")
        ax.plot(osm, inter + slope * np.asarray(osm), color="black", lw=0.8)
    ax.set_xlabel("normal quantile", fontsize=8)
    ax.set_ylabel("residual quantile", fontsize=8)
    ax.set_title("Normal QQ plot of residuals (points near the line)", fontsize=8.5)
    _clean(ax)
    ax = axs[2]
    if f is None:
        ax.text(0.5, 0.5, f"No random effects: {NEEDS_2N_TXT}.", ha="center", va="center", fontsize=9,
                color="#555555", transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("Random effects per n", fontsize=8.5)
        _clean(ax)
        fig.suptitle(title, fontsize=9.5)
        _save(fig, out, f"Descriptive only (one n): residuals = deviations from each condition's mean in that n. "
                        f"Response: {res.scale.describe}.", 0.06)
        return
    re = f.random_effects()
    ns = sorted(re)
    tmap = {f"I{i}": t for i, t in enumerate(res.prepared.treatments)}
    for i, n in enumerate(ns):
        terms = re[n]
        ax.scatter(i - 0.15, float(terms["n"][0]), marker="s", s=40, color=ncol.get(n, "grey"), edgecolors="black",
                   lw=0.5, zorder=3)
        vals, labs = [], []
        for k, b in terms.items():
            if k.startswith("nt_"):
                vals += list(b)
                labs += [tmap.get(c, c) for c in f.vc_colnames.get(k, [""] * len(b))]
        vals = np.asarray(vals, float)
        ax.scatter(np.full(len(vals), i + 0.12) + np.linspace(-0.08, 0.08, len(vals)), vals, s=12,
                   color=ncol.get(n, "grey"), alpha=0.8, edgecolors="none", zorder=2)
        if len(vals):
            j = int(np.argmax(np.abs(vals)))
            ax.annotate(labs[j], (i + 0.2, vals[j]), fontsize=6, color="#555555", xytext=(4, 0),
                        textcoords="offset points", va="center")
    ax.axhline(0, color="grey", lw=0.8, ls="--")
    ax.set_xticks(range(len(ns)), [short_n(n) for n in ns])
    ax.set_ylabel("estimated random effect (analysis scale)", fontsize=8)
    ax.set_title("Random effects per n: square = n, dots = n × treatment\n(largest labelled)", fontsize=8.5)
    _clean(ax)
    fig.suptitle(title, fontsize=9.5)
    _save(fig, out, f"Model: {res.model}. Response: {res.scale.describe}. Estimation: {res.fit.method}; "
                    f"converged: {'yes' if res.converged else 'no'}.", 0.06)


TR_COLOURS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
              "#1f78b4", "#b2df8a", "#fb9a99", "#cab2d6"]
N_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]


def density_scatter(pts: pd.DataFrame, curves: pd.DataFrame, order: list[str], out: Path, title: str, footer: str,
                    xlabel: str, ylabel: str, by_n: bool = True):
    """Fraction (y, %) vs density (x) per condition: points (one per well or field; marker = n when by_n) and the
    fitted line of each condition (curves: treatment, x, y). One colour per condition."""
    trs = [t for t in order if t in set(pts["treatment"])]
    col = {t: TR_COLOURS[order.index(t) % len(TR_COLOURS)] for t in trs}  # stable across plots
    ns = sorted(set(pts["n"])) if by_n else []
    mk = {n: N_MARKERS[i % len(N_MARKERS)] for i, n in enumerate(ns)}
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for t in trs:
        g = pts[pts["treatment"] == t]
        if by_n:
            for n, gn in g.groupby("n"):
                ax.scatter(gn["x"], gn["y"], s=22, marker=mk[n], color=col[t], alpha=0.8, edgecolors="white",
                           linewidths=0.4, zorder=3)
        else:
            ax.scatter(g["x"], g["y"], s=14, color=col[t], alpha=0.75, edgecolors="none", zorder=3)
        c = curves[curves["treatment"] == t] if len(curves) else curves
        if len(c):
            ax.plot(c["x"], c["y"], color=col[t], lw=1.6, zorder=4)
    ax.set_xlabel(xlabel, fontsize=8.5)
    ax.set_ylabel(ylabel, fontsize=8.5)
    ax.set_ylim(bottom=0)
    ax.set_title(title, fontsize=9)
    _clean(ax)
    handles = [Line2D([], [], marker="o", ls="-", color=col[t], markersize=5, label=t) for t in trs]
    if by_n and len(ns) > 1:
        handles += [Line2D([], [], marker=mk[n], ls="", color="#666666", markersize=5, label=short_n(n)) for n in ns]
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    _save(fig, out, footer, 0.12)
    return True


# ---------------------------------------------------------------- design: wells and positions

def fields_by_well(sub: pd.DataFrame, order: list[str], label: str, out: Path, title: str, positions_ok: bool):
    """One n: every field as a dot, colour = well, marker = field position."""
    sub = sub.dropna(subset=["value"])
    if sub.empty:
        return False
    trs = [t for t in order if t in set(sub["treatment"])]
    wells = sorted(sub["well"].unique())
    poss = sorted(sub["position"].unique())
    fig, ax = plt.subplots(figsize=(max(5.5, 0.8 * len(trs) + 2.8), 4.4))
    nw = len(wells)
    for i, t in enumerate(trs):
        g = sub[sub["treatment"] == t]
        for j, w in enumerate(wells):
            gw = g[g["well"] == w]
            if gw.empty:
                continue
            x0 = i + (j - (nw - 1) / 2) * 0.28
            xs = x0 + (gw["position"].to_numpy() - np.mean(poss)) * (0.2 / max(len(poss), 1))
            for x, yv, pos in zip(xs, gw["value"], gw["position"]):
                ax.scatter(x, yv, s=20, marker=POS_MARKERS[(int(pos) - 1) % len(POS_MARKERS)] if positions_ok else "o",
                           color=WELL_COLOURS[j % len(WELL_COLOURS)], edgecolors="white", linewidths=0.4, zorder=3)
            ax.hlines(gw["value"].mean(), x0 - 0.12, x0 + 0.12, color=WELL_COLOURS[j % len(WELL_COLOURS)], lw=1.6, zorder=2)
    ax.set_xticks(range(len(trs)), trs, rotation=45, ha="right")
    ax.set_ylabel(label, fontsize=9)
    ax.set_title(title + "\ncolour = well (line = well mean)" + ("; marker = field position" if positions_ok else ""),
                 fontsize=9)
    _clean(ax)
    handles = [Line2D([], [], marker="o", ls="", color=WELL_COLOURS[j % len(WELL_COLOURS)], markersize=5,
                      label=f"well {int(w)}") for j, w in enumerate(wells)]
    if positions_ok:
        handles += [Line2D([], [], marker=POS_MARKERS[(int(pp) - 1) % len(POS_MARKERS)], ls="", color="#666666",
                           markersize=5, label=f"position {int(pp)}") for pp in poss]
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    _save(fig, out, "Fields are technical replicates: this shows where the variation within one experiment comes from.",
          0.06)
    return True


def position_effect(sub: pd.DataFrame, out: Path, title: str, by: str, colours: dict, label_fn=str):
    """Mean of (field / mean of its treatment in that n) by field position; one line per `by` level
    (per n: per well; all n: per n). 1 = average field."""
    s = sub.dropna(subset=["value"])
    s = s[s["value"] > 0]
    if s.empty or s["position"].nunique() < 2:
        return False
    key = ["n", "treatment"]
    s = s.assign(rel=s["value"] / s.groupby(key)["value"].transform("mean"))
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    poss = sorted(s["position"].unique())
    for lv, g in s.groupby(by):
        m = g.groupby("position")["rel"].agg(["mean", "sem", "count"]).reindex(poss)
        ax.errorbar(m.index, m["mean"], yerr=m["sem"].fillna(0), color=colours.get(lv, "grey"), marker="o", ms=4,
                    lw=1.2, capsize=2, label=label_fn(lv))
    allm = s.groupby("position")["rel"].mean().reindex(poss)
    ax.plot(allm.index, allm.values, color="black", lw=2.2, label="all", zorder=5)
    ax.axhline(1, color="grey", lw=0.8, ls="--")
    ax.set_xticks(poss)
    ax.set_xlabel("field position within the well", fontsize=8)
    ax.set_ylabel("field / mean of its treatment (1 = average)", fontsize=8)
    ax.set_title(title, fontsize=9)
    _clean(ax)
    ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    _save(fig, out, "Each field is divided by the mean of its treatment in that n, so treatment effects cancel and "
                    "only the position pattern remains (mean ± SEM).", 0.06)
    return True


def treatment_fields(sub: pd.DataFrame, out: Path, title: str, positions_ok: bool):
    """Per treatment in one n: count per image, bar colour = well, position under each bar; x = excluded."""
    fig, ax = plt.subplots(figsize=(max(5, 0.38 * len(sub) + 1.5), 3.2))
    xs = np.arange(len(sub))
    wells = sorted(sub["well"].dropna().unique()) if "well" in sub else []
    wmap = {w: WELL_COLOURS[j % len(WELL_COLOURS)] for j, w in enumerate(wells)}
    for x, r in zip(xs, sub.itertuples()):
        if r.excluded:
            ax.text(x, 0, "x", ha="center", va="bottom", color="red")
            continue
        ax.bar(x, r.count, color=wmap.get(getattr(r, "well", None), "#9ecae1"), edgecolor="black", lw=0.5)
    labels = [f"{int(r.image)}" + (f"\np{int(r.position)}" if positions_ok and "position" in sub and
                                    pd.notna(r.position) else "") for r in sub.itertuples()]
    ax.set_xticks(xs, labels, fontsize=7)
    ax.set_xlabel("image" + (" / field position" if positions_ok else ""), fontsize=8)
    ax.set_ylabel("nuclei", fontsize=8)
    ax.set_title(title, fontsize=9)
    _clean(ax)
    if wells:
        ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c, label=f"well {int(w)}") for w, c in wmap.items()],
                  fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
