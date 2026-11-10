"""How many images are needed? Variance components, predicted power, and subsampling of real data.

Design assumed: each n x treatment = `plates` plates x `fields_per_plate` fields (images).
Variation in a per-n treatment mean comes from
    field-to-field (within a plate)      -> shrinks with more fields
    plate-to-plate (within an n)         -> shrinks only with more plates
    n-to-n (biological, of the ratio vs reference) -> shrinks only with more n
The across-n test (ratio paired t-test vs reference) has variance per n of
    s2_bio + 2 * (cv2_plate / p + cv2_field / (p * k))     (log scale, delta method)
so the value of extra images depends on how big the measurement part is relative to s2_bio.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from . import analysis


def plate_of(image: int, fields_per_plate: int, plates: int, layout: str) -> int:
    if layout == "alternating":
        return (image - 1) % plates + 1
    return min((image - 1) // fields_per_plate + 1, plates)


def variance_components(df: pd.DataFrame) -> dict:
    """Pooled nested ANOVA on counts scaled by their n x treatment mean (so components are CV^2)."""
    ss_f = df_f = ss_p = df_p = 0.0
    k_list = []
    for _, g in df.groupby(["n", "treatment"]):
        m = g["value"].mean()
        if m <= 0 or len(g) < 2:
            continue
        y = g["value"] / m
        for _, pg in y.groupby(g["plate"]):
            if len(pg) > 1:
                ss_f += ((pg - pg.mean()) ** 2).sum()
                df_f += len(pg) - 1
        pm = y.groupby(g["plate"]).agg(["mean", "size"])
        if len(pm) > 1:
            ss_p += (pm["size"] * (pm["mean"] - y.mean()) ** 2).sum()
            df_p += len(pm) - 1
            k_list.append(pm["size"].mean())
    ms_f = ss_f / df_f if df_f else np.nan
    ms_p = ss_p / df_p if df_p else np.nan
    k = float(np.mean(k_list)) if k_list else np.nan
    cv2_plate = max((ms_p - ms_f) / k, 0.0) if df_p else 0.0
    return dict(cv2_field=ms_f, cv2_plate=cv2_plate, ms_plate=ms_p, df_field=df_f, df_plate=df_p, k=k)


def _power(effect: float, sd: float, n: int, alpha: float = 0.05) -> float:
    if n < 2 or sd <= 0 or np.isnan(sd) or np.isnan(effect):
        return np.nan
    dfree = n - 1
    nc = abs(effect) / (sd / np.sqrt(n))
    tc = stats.t.ppf(1 - alpha / 2, dfree)
    return float(1 - stats.nct.cdf(tc, dfree, nc) + stats.nct.cdf(-tc, dfree, nc))


def _n_needed(effect, sd, target=0.8, alpha=0.05, cap=50):
    for n in range(2, cap + 1):
        p = _power(effect, sd, n, alpha)
        if not np.isnan(p) and p >= target:
            return n
    return np.nan


def _per_n_ratios(df: pd.DataFrame, ref: str) -> pd.DataFrame:
    means = df.groupby(["n", "treatment"])["value"].mean().unstack("treatment")
    if ref not in means:
        return pd.DataFrame()
    return np.log(means.drop(columns=ref).div(means[ref], axis=0)).replace([np.inf, -np.inf], np.nan)


def run(df: pd.DataFrame, settings: dict, out: Path, metric_label: str = "nuclei per field",
        iterations: int = 1000, seed: int = 0) -> pd.DataFrame:
    """df columns: n, treatment, image, value (excluded images already removed)."""
    out.mkdir(parents=True, exist_ok=True)
    d = settings["design"]
    P, K, layout = int(d["plates"]), int(d["fields_per_plate"]), d["plate_layout"]
    ref = settings["analysis"]["reference"]
    df = df.assign(plate=[plate_of(int(i), K, P, layout) for i in df["image"]])

    vc = variance_components(df)
    full_lr = _per_n_ratios(df, ref)
    order = analysis.treatment_order(full_lr.columns, settings["analysis"]["treatment_order"])
    meas = lambda p, k: vc["cv2_plate"] / p + vc["cv2_field"] / (p * k)

    # biological variance of log ratio per treatment = observed variance - measurement part at the full design
    rows = []
    for t in order:
        lr = full_lr[t].dropna()
        if len(lr) < 2:
            continue
        s2_obs = lr.var(ddof=1)
        s2_bio = max(s2_obs - 2 * meas(P, K), 0.0)
        rows.append(dict(treatment=t, n_pairs=len(lr), ratio=float(np.exp(lr.mean())), log_ratio=lr.mean(),
                         s2_observed=s2_obs, s2_bio=s2_bio, s2_meas_full=2 * meas(P, K)))
    eff = pd.DataFrame(rows)

    designs = [(p, k) for p in range(1, P + 1) for k in range(1, K + 1)]
    pw = []
    for p, k in designs:
        for _, e in eff.iterrows():
            sd = np.sqrt(e["s2_bio"] + 2 * meas(p, k))
            sd_full = np.sqrt(e["s2_bio"] + 2 * meas(P, K))
            pw.append(dict(plates=p, fields_per_plate=k, images=p * k, treatment=e["treatment"],
                           ratio=e["ratio"], sd_log_ratio=sd,
                           power_at_current_n=_power(e["log_ratio"], sd, int(e["n_pairs"])),
                           n_for_80pct_power=_n_needed(e["log_ratio"], sd),
                           variance_vs_current=(sd / sd_full) ** 2 if sd_full > 0 else np.nan))
    pw = pd.DataFrame(pw)

    # empirical: subsample real images, redo the across-n test, compare with full data
    rng = np.random.default_rng(seed)
    full_p = {}
    for t in order:
        lr = full_lr[t].dropna()
        full_p[t] = stats.ttest_1samp(lr, 0).pvalue if len(lr) >= 2 else np.nan
    groups = {key: g for key, g in df.groupby(["n", "treatment"])}
    emp = []
    for p, k in designs:
        agree = {t: 0 for t in order}
        err = {t: [] for t in order}
        valid = {t: 0 for t in order}
        for _ in range(iterations):
            parts = []
            for (n, tr), g in groups.items():
                plates = sorted(g["plate"].unique())
                chosen = rng.choice(plates, size=min(p, len(plates)), replace=False)
                for pl in chosen:
                    gg = g[g["plate"] == pl]
                    parts.append(gg.iloc[rng.choice(len(gg), size=min(k, len(gg)), replace=False)])
            sub_lr = _per_n_ratios(pd.concat(parts), ref)
            for t in order:
                if t not in sub_lr:
                    continue
                lr = sub_lr[t].dropna()
                if len(lr) < 2 or np.isnan(full_p[t]):
                    continue
                pv = stats.ttest_1samp(lr, 0).pvalue
                valid[t] += 1
                agree[t] += (pv < 0.05) == (full_p[t] < 0.05)
                err[t].extend((lr - full_lr[t].reindex(lr.index)).dropna().tolist())
        for t in order:
            if valid[t]:
                e = np.asarray(err[t])
                emp.append(dict(plates=p, fields_per_plate=k, images=p * k, treatment=t,
                                same_conclusion_pct=100 * agree[t] / valid[t],
                                per_n_ratio_error_pct=100 * (np.exp(np.sqrt(np.mean(e ** 2))) - 1) if len(e) else np.nan))
    emp = pd.DataFrame(emp)

    # ---- write
    pd.DataFrame([vc]).to_csv(out / "variance_components.csv", index=False)
    eff.to_csv(out / "effects.csv", index=False)
    pw.to_csv(out / "design_power.csv", index=False)
    emp.to_csv(out / "subsampling.csv", index=False)
    summary = (pw.groupby(["plates", "fields_per_plate", "images"])
               .agg(median_variance_vs_current=("variance_vs_current", "median"),
                    median_n_for_80pct_power=("n_for_80pct_power", "median")).reset_index())
    if len(emp):
        summary = summary.merge(emp.groupby(["plates", "fields_per_plate"]).agg(
            same_conclusion_pct=("same_conclusion_pct", "mean"),
            per_n_ratio_error_pct=("per_n_ratio_error_pct", "median")).reset_index(), on=["plates", "fields_per_plate"])
    summary.to_csv(out / "summary_by_design.csv", index=False)
    _plots(summary, pw, out, P)
    _report(vc, eff, summary, out, P, K, layout, ref, metric_label, iterations)
    return summary


def _plots(summary, pw, out, P):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    for p in sorted(summary["plates"].unique()):
        s = summary[summary["plates"] == p]
        lab = f"{p} plate{'s' if p > 1 else ''}"
        axes[0].plot(s["images"], s["median_variance_vs_current"], "o-", label=lab)
        if "same_conclusion_pct" in s:
            axes[1].plot(s["images"], s["same_conclusion_pct"], "o-", label=lab)
            axes[2].plot(s["images"], s["per_n_ratio_error_pct"], "o-", label=lab)
    axes[0].axhline(1, color="grey", ls="--", lw=0.8)
    axes[0].set(xlabel="images per condition per n", ylabel="variance of test / current design",
                title="Predicted noise in the across-n test")
    axes[1].set(xlabel="images per condition per n", ylabel="% same conclusion as 12 images",
                title="Subsampling real data", ylim=(0, 105))
    axes[2].set(xlabel="images per condition per n", ylabel="typical error of per-n ratio (%)",
                title="Per-n estimate vs all images")
    for a in axes:
        a.legend(fontsize=8)
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "images_needed.png", dpi=180)
    fig.savefig(out / "images_needed.pdf")
    plt.close(fig)


def _report(vc, eff, summary, out, P, K, layout, ref, label, iterations):
    meas_full = 2 * (vc["cv2_plate"] / P + vc["cv2_field"] / (P * K))
    bio = eff["s2_bio"].median() if len(eff) else np.nan
    lines = [
        f"How many images are needed? ({label})",
        "=" * 60,
        f"Design analysed: {P} plates x {K} fields (plate layout: {layout}); reference group: {ref}",
        f"Subsampling iterations per design: {iterations}",
        "",
        "Where the variation comes from (as coefficient of variation of one field / one plate):",
        f"  field-to-field within a plate : CV = {np.sqrt(vc['cv2_field']):.0%}",
        f"  plate-to-plate within an n    : CV = {np.sqrt(vc['cv2_plate']):.0%}",
        f"  measurement noise in one per-n ratio at {P*K} images: {np.sqrt(meas_full):.1%} (log scale SD)",
        f"  biological n-to-n noise in the ratio (median over treatments): {np.sqrt(bio):.1%}" if not np.isnan(bio) else "",
        "",
        "Design table (median over treatments):",
        summary.round(3).to_string(index=False),
        "",
        "Reading the table",
        "  median_variance_vs_current : 1.0 = as precise as now; 1.2 = 20% more noise in the across-n test",
        "  median_n_for_80pct_power   : n's needed to detect the observed effects with 80% power",
        "  same_conclusion_pct        : how often a real-data subsample gives the same significant/ns call",
        "  per_n_ratio_error_pct      : typical deviation of one n's treatment/reference ratio from the 12-image value",
        "",
        "Caveats: with few n's the biological variance is itself uncertain, and the observed effects are",
        "used as if they were the true effects. Treat this as guidance for the design of the next experiments.",
    ]
    (out / "REPORT.txt").write_text("\n".join(l for l in lines if l is not None), encoding="utf-8")
