"""How many images are needed? Variance components, batch checks, predicted power, and subsampling of real data.

Experimental layout (see CONTEXT.md):
    n      = one independent experiment = one plate, seeded/treated/imaged on one day
    well   = each condition sits in `wells` wells of that plate
    field  = each well is imaged `fields_per_well` times (one image = one field)

Variation in a per-n condition mean comes from
    field-to-field (within a well)  -> shrinks with more fields
    well-to-well   (within a plate) -> shrinks only with more wells
    n-to-n (biological + plate/day batch, which cannot be separated: one plate per n) -> shrinks only with more n
The across-n test (ratio paired t-test vs a reference, as in analysis.all_n) has variance per n of
    s2_bio + 2 * (cv2_well / w + cv2_field / (w * k))     (log scale, delta method)
so the value of extra images depends on how big the measurement part is relative to s2_bio.
Plate-wide multiplicative batch effects (seeding density, imaging day) cancel in the ratio; batch effects that change
the size of the treatment effect (e.g. transfection efficiency) do not, and are reported in the batch checks.
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

WELL_Z_FLAG = 2.5  # well pair flagged when its difference is this many field-noise SEs (two-sided p ~ 0.012)
SHIFT_FLAG = 1.5  # n flagged when its treatment ratios are on average >1.5x or <1/1.5 of the across-n value...
SHIFT_AGREE = 0.8  # ...and at least 80% of treatments move in the same direction


def design(settings: dict) -> tuple[int, int, str, list[str]]:
    """(wells, fields_per_well, well_layout, notes). Accepts the old names plates / fields_per_plate / plate_layout."""
    d = settings["design"]
    notes = []
    old = {"plates": "wells", "fields_per_plate": "fields_per_well", "plate_layout": "well_layout"}
    vals = {}
    for o, new in old.items():
        if o in d:  # old names are no longer in the defaults, so their presence means the user's file uses them
            vals[new] = d[o]
            notes.append(f"settings.toml [design] uses the old name '{o}'; please rename it to '{new}'")
        else:
            vals[new] = d[new]
    return int(vals["wells"]), int(vals["fields_per_well"]), vals["well_layout"], notes


def well_of(image: int, fields_per_well: int, wells: int, layout: str) -> int:
    if layout == "alternating":
        return (image - 1) % wells + 1
    return min((image - 1) // fields_per_well + 1, wells)


def variance_components(df: pd.DataFrame) -> dict:
    """Pooled nested ANOVA on counts scaled by their n x treatment mean (so components are CV^2)."""
    ss_f = df_f = ss_w = df_w = 0.0
    k_list = []
    for _, g in df.groupby(["n", "treatment"]):
        m = g["value"].mean()
        if m <= 0 or len(g) < 2:
            continue
        y = g["value"] / m
        for _, wg in y.groupby(g["well"]):
            if len(wg) > 1:
                ss_f += ((wg - wg.mean()) ** 2).sum()
                df_f += len(wg) - 1
        wm = y.groupby(g["well"]).agg(["mean", "size"])
        if len(wm) > 1:
            ss_w += (wm["size"] * (wm["mean"] - y.mean()) ** 2).sum()
            df_w += len(wm) - 1
            k_list.append(wm["size"].mean())
    ms_f = ss_f / df_f if df_f else np.nan
    ms_w = ss_w / df_w if df_w else np.nan
    k = float(np.mean(k_list)) if k_list else np.nan
    cv2_well = max((ms_w - ms_f) / k, 0.0) if df_w else 0.0
    return dict(cv2_field=ms_f, cv2_well=cv2_well, ms_well=ms_w, df_field=df_f, df_well=df_w, k=k)


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


def _log_ratios(means: pd.DataFrame, ref: str) -> pd.DataFrame:
    """n x treatment log(mean / reference mean); n's without the reference are dropped."""
    if ref not in means:
        return pd.DataFrame()
    m = means[means[ref].notna() & (means[ref] > 0)]
    with np.errstate(divide="ignore", invalid="ignore"):  # a zero mean gives -inf, mapped to NaN below
        return np.log(m.drop(columns=ref).div(m[ref], axis=0)).replace([np.inf, -np.inf], np.nan)


def _means(df: pd.DataFrame) -> pd.DataFrame:
    return df.groupby(["n", "treatment"])["value"].mean().unstack("treatment")


# ---------------------------------------------------------------- batch checks
def batch_checks(df: pd.DataFrame, refs: list[str], cv2_field: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per condition (n x treatment): well means and a flag for wells that disagree more than field noise explains;
    ratios vs each reference. Per n: reference levels, missing references, systematic well difference, and how far
    this n's treatment effects sit from the across-n value (a batch effect on the effect size)."""
    means = _means(df)
    wm = df.groupby(["n", "treatment", "well"])["value"].agg(["mean", "size"])
    cond = []
    for (n, t), g in wm.groupby(level=["n", "treatment"]):
        g = g.droplevel(["n", "treatment"])
        r = dict(n=n, treatment=t, images=int(g["size"].sum()), mean=float(means.loc[n, t]))
        for w, row in g.iterrows():
            r[f"well{w}_mean"] = float(row["mean"])
        r["well_ratio_last_first"] = r["well_z"] = np.nan
        r["well_flag"] = ""
        if len(g) >= 2:
            lo, hi = g.index.min(), g.index.max()
            a, b = g.loc[lo], g.loc[hi]
            if a["mean"] > 0 and b["mean"] > 0:
                lr = np.log(b["mean"] / a["mean"])
                se = np.sqrt(cv2_field * (1 / a["size"] + 1 / b["size"]))
                r["well_ratio_last_first"] = float(np.exp(lr))
                r["well_z"] = float(lr / se) if se > 0 else np.nan
                r["well_flag"] = "check" if abs(r["well_z"]) > WELL_Z_FLAG else ""
            elif (a["mean"] > 0) != (b["mean"] > 0):
                r["well_flag"] = "check"  # one well empty, the other not
        for ref in refs:
            if ref in means and t != ref and pd.notna(means.loc[n, ref]) and means.loc[n, ref] > 0:
                r[f"ratio_vs_{ref}"] = float(means.loc[n, t] / means.loc[n, ref])
        cond.append(r)
    cond = pd.DataFrame(cond)

    per_n = []
    controls = set(refs)
    for n in means.index:
        r = dict(n=n, images=int(cond.loc[cond["n"] == n, "images"].sum()))
        for ref in refs:
            r[f"{ref}_mean"] = float(means.loc[n, ref]) if ref in means and pd.notna(means.loc[n, ref]) else np.nan
        r["missing_references"] = ", ".join(ref for ref in refs if ref not in means or pd.isna(means.loc[n, ref]))
        wr = cond.loc[cond["n"] == n, "well_ratio_last_first"].dropna()
        r["well_last_vs_first_geo"] = float(np.exp(np.log(wr).mean())) if len(wr) else np.nan
        r["wells_flagged"] = int((cond.loc[cond["n"] == n, "well_flag"] == "check").sum())
        per_n.append(r)
    per_n = pd.DataFrame(per_n).set_index("n")

    # effect-size shift: this n's log ratio minus the across-n mean, averaged over (non-control) treatments
    for ref in refs:
        lr = _log_ratios(means, ref)
        lr = lr[[c for c in lr.columns if c not in controls]]
        if lr.shape[0] < 2 or lr.shape[1] == 0:
            continue
        dev = lr - lr.mean(axis=0)
        shift = dev.mean(axis=1)
        agree = dev.apply(lambda row: (np.sign(row.dropna()) == np.sign(row.dropna().mean())).mean(), axis=1)
        per_n[f"shift_vs_{ref}"] = np.exp(shift)  # index-aligned: n's without this reference stay NaN
        per_n[f"shift_vs_{ref}_agree"] = agree
        per_n[f"shift_vs_{ref}_flag"] = pd.Series(
            ["check" if (abs(s) > np.log(SHIFT_FLAG) and a >= SHIFT_AGREE) else "" for s, a in zip(shift, agree)],
            index=shift.index)
    return cond, per_n.reset_index()


# ---------------------------------------------------------------- main
def run(df: pd.DataFrame, settings: dict, out: Path, metric_label: str = "nuclei per field",
        iterations: int = 1000, seed: int = 0) -> pd.DataFrame:
    """df columns: n, treatment, image, value (excluded images already removed)."""
    out.mkdir(parents=True, exist_ok=True)
    W, K, layout, notes = design(settings)
    a = settings["analysis"]
    df = df.assign(well=[well_of(int(i), K, W, layout) for i in df["image"]])
    means = _means(df)
    refs = [r for r in dict.fromkeys([a["reference"]] + list(a["normalise_to"])) if r in means]

    vc = variance_components(df)
    meas = lambda w, k: vc["cv2_well"] / w + vc["cv2_field"] / (w * k)
    cond, per_n = batch_checks(df, refs, vc["cv2_field"])

    full_lr = {ref: _log_ratios(means, ref) for ref in refs}
    designs = [(w, k) for w in range(1, W + 1) for k in range(1, K + 1)]

    # biological variance of log ratio per treatment = observed variance - measurement part at the full design
    eff, pw, full_p = [], [], {}
    for ref in refs:
        lrs = full_lr[ref]
        for t in analysis.treatment_order(lrs.columns, a["treatment_order"]):
            lr = lrs[t].dropna()
            if len(lr) < 2:
                continue
            s2_obs = lr.var(ddof=1)
            s2_bio = max(s2_obs - 2 * meas(W, K), 0.0)
            p_full = stats.ttest_1samp(lr, 0).pvalue
            full_p[(ref, t)] = p_full
            eff.append(dict(reference=ref, treatment=t, n_pairs=len(lr), ratio=float(np.exp(lr.mean())),
                            log_ratio=lr.mean(), p_full=p_full, s2_observed=s2_obs, s2_bio=s2_bio,
                            s2_meas_full=2 * meas(W, K)))
    eff = pd.DataFrame(eff)
    for w, k in designs:
        for _, e in eff.iterrows():
            sd = np.sqrt(e["s2_bio"] + 2 * meas(w, k))
            sd_full = np.sqrt(e["s2_bio"] + 2 * meas(W, K))
            pw.append(dict(reference=e["reference"], wells=w, fields_per_well=k, images=w * k,
                           treatment=e["treatment"], ratio=e["ratio"], sd_log_ratio=sd,
                           power_at_current_n=_power(e["log_ratio"], sd, int(e["n_pairs"])),
                           n_for_80pct_power=_n_needed(e["log_ratio"], sd),
                           variance_vs_current=(sd / sd_full) ** 2 if sd_full > 0 else np.nan))
    pw = pd.DataFrame(pw)

    # empirical: subsample real images, redo the across-n test, compare with full data (all references per draw)
    rng = np.random.default_rng(seed)
    groups = {key: g for key, g in df.groupby(["n", "treatment"])}
    emp = []
    for w, k in designs:
        agree = {key: 0 for key in full_p}
        valid = {key: 0 for key in full_p}
        err = {key: [] for key in full_p}
        for _ in range(iterations):
            parts = []
            for g in groups.values():
                wells = sorted(g["well"].unique())
                for wl in rng.choice(wells, size=min(w, len(wells)), replace=False):
                    gg = g[g["well"] == wl]
                    parts.append(gg.iloc[rng.choice(len(gg), size=min(k, len(gg)), replace=False)])
            sub_means = _means(pd.concat(parts))
            for ref in refs:
                sub_lr = _log_ratios(sub_means, ref)
                for (r_, t) in full_p:
                    if r_ != ref or t not in sub_lr:
                        continue
                    lr = sub_lr[t].dropna()
                    if len(lr) < 2:
                        continue
                    pv = stats.ttest_1samp(lr, 0).pvalue
                    valid[(ref, t)] += 1
                    agree[(ref, t)] += (pv < 0.05) == (full_p[(ref, t)] < 0.05)
                    err[(ref, t)].extend((lr - full_lr[ref][t].reindex(lr.index)).dropna().tolist())
        for (ref, t) in full_p:
            if valid[(ref, t)]:
                e = np.asarray(err[(ref, t)])
                emp.append(dict(reference=ref, wells=w, fields_per_well=k, images=w * k, treatment=t,
                                same_conclusion_pct=100 * agree[(ref, t)] / valid[(ref, t)],
                                per_n_ratio_error_pct=100 * (np.exp(np.sqrt(np.mean(e ** 2))) - 1) if len(e) else np.nan))
    emp = pd.DataFrame(emp)

    # ---- write
    pd.DataFrame([vc]).to_csv(out / "variance_components.csv", index=False)
    eff.to_csv(out / "effects.csv", index=False)
    pw.to_csv(out / "design_power.csv", index=False)
    emp.to_csv(out / "subsampling.csv", index=False)
    cond.to_csv(out / "batch_by_condition.csv", index=False)
    per_n.to_csv(out / "batch_by_n.csv", index=False)
    keys = ["reference", "wells", "fields_per_well", "images"]
    summary = (pw.groupby(keys)
               .agg(median_variance_vs_current=("variance_vs_current", "median"),
                    median_n_for_80pct_power=("n_for_80pct_power", "median")).reset_index())
    if len(emp):
        summary = summary.merge(emp.groupby(keys).agg(
            same_conclusion_pct=("same_conclusion_pct", "mean"),
            per_n_ratio_error_pct=("per_n_ratio_error_pct", "median")).reset_index(), on=keys)
    summary.to_csv(out / "summary_by_design.csv", index=False)
    for i, ref in enumerate(refs):
        _plots(summary[summary["reference"] == ref], out, ref,
               "images_needed" if i == 0 else f"images_needed_vs_{ref}")
    _report(vc, eff, summary, cond, per_n, out, W, K, layout, refs, means, notes, metric_label, iterations)
    return summary


def _plots(summary, out, ref, stem):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    for w in sorted(summary["wells"].unique()):
        s = summary[summary["wells"] == w]
        lab = f"{w} well{'s' if w > 1 else ''}"
        axes[0].plot(s["images"], s["median_variance_vs_current"], "o-", label=lab)
        if "same_conclusion_pct" in s:
            axes[1].plot(s["images"], s["same_conclusion_pct"], "o-", label=lab)
            axes[2].plot(s["images"], s["per_n_ratio_error_pct"], "o-", label=lab)
    axes[0].axhline(1, color="grey", ls="--", lw=0.8)
    axes[0].set(xlabel="images per condition per n", ylabel="variance of test / current design",
                title=f"Predicted noise in the across-n test (vs {ref})")
    axes[1].set(xlabel="images per condition per n", ylabel="% same conclusion as all images",
                title="Subsampling real data", ylim=(0, 105))
    axes[2].set(xlabel="images per condition per n", ylabel="typical error of per-n ratio (%)",
                title="Per-n estimate vs all images")
    for ax in axes:
        ax.legend(fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / f"{stem}.png", dpi=180)
    fig.savefig(out / f"{stem}.pdf")
    plt.close(fig)


def _fmt_table(df: pd.DataFrame) -> list[str]:
    return ["  " + line for line in df.to_string(index=False).splitlines()]


def _report(vc, eff, summary, cond, per_n, out, W, K, layout, refs, means, notes, label, iterations):
    meas_full = 2 * (vc["cv2_well"] / W + vc["cv2_field"] / (W * K))
    ns = list(means.index)
    lay = (f"images 1-{K} = well 1, {K + 1}-{2 * K} = well 2" if layout == "consecutive" and W == 2
           else f"well layout: {layout}")
    L = [
        f"How many images are needed? ({label})",
        "=" * 70,
        f"Layout: one plate per n; each condition in {W} wells x {K} fields ({lay}).",
        f"n's: {len(ns)} ({', '.join(ns)}); references analysed: {', '.join(refs)}",
        f"Subsampling iterations per design: {iterations}",
    ]
    L += [f"NOTE: {x}" for x in notes]
    L += [
        "",
        "1. Where the variation comes from (coefficient of variation of one field / one well)",
        "-" * 70,
        f"  field-to-field within a well  : CV = {np.sqrt(vc['cv2_field']):.0%}",
        f"  well-to-well within a plate   : CV = {np.sqrt(vc['cv2_well']):.0%}",
        f"  measurement noise in one per-n ratio at {W * K} images: {np.sqrt(meas_full):.1%} (log-scale SD)",
        "  plate, day and passage effects cannot be separated from n-to-n variation (one plate per n).",
        "  Plate-wide effects (seeding density, imaging day) cancel in the treatment/reference ratio;",
        "  effects that change the size of the treatment response (e.g. transfection efficiency) do not.",
        "",
        "2. Batch checks",
        "-" * 70,
        "  Reference groups per n (mean " + label + "):",
    ]
    ref_tab = per_n[["n"] + [f"{r}_mean" for r in refs]].rename(columns={f"{r}_mean": r for r in refs})
    L += _fmt_table(ref_tab.round(1).fillna("-"))
    for _, r in per_n.iterrows():
        if r["missing_references"]:
            for ref in r["missing_references"].split(", "):
                L.append(f"  ! {r['n']} has no {ref}: it is left out of every comparison vs {ref}.")
    L += ["", "  Does one n respond more or less strongly than the others? (ratio of this n's treatment/reference",
          "  ratios to the across-n value, averaged over treatments; controls excluded. 1.00 = typical;",
          "  agree = share of treatments shifted the same way)"]
    for ref in refs:
        col = f"shift_vs_{ref}"
        if col not in per_n:
            continue
        t = per_n[["n", col, f"{col}_agree", f"{col}_flag"]].dropna(subset=[col])
        t = t.rename(columns={col: f"vs {ref}", f"{col}_agree": "agree", f"{col}_flag": "flag"})
        L.append(f"  vs {ref}:")
        L += ["  " + s for s in _fmt_table(t.round(2))]
        for _, r in t.iterrows():
            if r["flag"]:
                direction = "stronger" if r[f"vs {ref}"] < 1 else "weaker"
                L.append(f"    ! {r['n']}: ratios {r[f'vs {ref}']:.2f}x the typical value in {r['agree']:.0%} of "
                         f"treatments -> effects look {direction} (if the treatments lower the count) in this n.")
    L += ["  A flagged n is a batch effect on the effect size: check what differed on that day (transfection, cell",
          "  passage, confluence, timepoint). It is real n-to-n variation unless there is a technical reason to exclude.",
          "", f"  Wells: geometric mean of well {W} / well 1 per n (1.00 = no systematic well difference):"]
    L += _fmt_table(per_n[["n", "well_last_vs_first_geo", "wells_flagged"]].round(2)
                    .rename(columns={"well_last_vs_first_geo": f"well{W}/well1", "wells_flagged": "pairs flagged"}))
    flagged = cond[cond["well_flag"] == "check"]
    if len(flagged):
        L.append(f"  Well pairs that differ more than field noise explains (|z| > {WELL_Z_FLAG}); check the review images:")
        for _, r in flagged.iterrows():
            L.append(f"    {r['n']:<16} {r['treatment']:<12} well{W}/well1 = {r['well_ratio_last_first']:.2f}  "
                     f"(z = {r['well_z']:+.1f})")
        L.append(f"  (with {cond['well_z'].notna().sum()} pairs, about {cond['well_z'].notna().sum() * 0.012:.1f} "
                 "would be flagged by chance)")
    else:
        L.append("  No well pair differs more than field noise explains.")

    for i, ref in enumerate(refs):
        e = eff[eff["reference"] == ref]
        s = summary[summary["reference"] == ref].drop(columns="reference")
        if not len(e):
            continue
        n_used = int(e["n_pairs"].max())
        bio = e["s2_bio"].median()
        L += ["", f"{3 + i}. Design table vs {ref} ({n_used} n's; median over treatments)", "-" * 70,
              f"  biological n-to-n noise in the ratio (median over treatments): {np.sqrt(bio):.1%} (log-scale SD)",
              f"  measurement noise at {W * K} images: {np.sqrt(meas_full):.1%}"]
        L += _fmt_table(s.round(3))
        zero = e.loc[e["s2_bio"] == 0, "treatment"].tolist()
        if zero:
            L.append(f"  ! n-to-n variance estimated as 0 for {', '.join(zero)} (too few n's to estimate it):")
            L.append("    their predictions are over-optimistic about cutting images.")
        border = e.loc[e["p_full"].between(0.01, 0.2), "treatment"].tolist()
        if border:
            L.append(f"  ! p with all images is between 0.01 and 0.2 for {', '.join(border)}: any subsample can flip")
            L.append("    the call, so same_conclusion_pct is low for them whatever the design.")
    L += [
        "",
        "Reading the design table",
        f"  wells x fields_per_well     : the design being simulated ({W} x {K} = now)",
        "  median_variance_vs_current : 1.0 = as precise as now; 1.2 = 20% more noise in the across-n test",
        "  median_n_for_80pct_power   : n's needed to detect the observed effects with 80% power",
        "  same_conclusion_pct        : how often a real-data subsample gives the same significant/ns call",
        "  per_n_ratio_error_pct      : typical deviation of one n's treatment/reference ratio from the all-image value",
        "  At equal image numbers, compare 1 well x k fields with 2 wells x k/2 fields: spreading over wells",
        "  protects against one bad well.",
        "",
        "Files: batch_by_n.csv, batch_by_condition.csv (per-well means, ratios, flags), effects.csv,",
        "design_power.csv, subsampling.csv, summary_by_design.csv, variance_components.csv.",
        "",
        "Caveats: with few n's the biological variance is itself uncertain, and the observed effects are",
        "used as if they were the true effects. Treat this as guidance for the design of the next experiments.",
    ]
    (out / "REPORT.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
