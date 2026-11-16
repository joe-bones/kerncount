"""Readout ladder: the same treatment effect measured at every layer of the intensity analysis, side by side.

ADR docs/adr/0001 ("Readout ladder"), docs/intensity.md section "Readout ladder". For every marker channel the
intensity stage writes per-image columns that each remove one analysis choice:

    1 intensity_<ch>_field_raw           whole field as acquired (no flat field, no background, no segmentation)
    2 intensity_<ch>_field_corrected     whole field, flat-field corrected and background subtracted
    3 intensity_<ch>_nucleus_mean_raw    DAPI nuclei, raw (median of the per-nucleus means; Fiji-like)
    4 intensity_<ch>_nucleus_mean_bgsub  DAPI nuclei, shrunk outline, local background, flat field
    5 readout_frac_positive_<ch>         % positive nuclei (cut-off)
    extras: intensity_<ch>_stained_area_mean (MFI in the stained area, no DAPI), intensity_<ch>_field_total

They are ordinary registry metrics (metrics.py), so the generation analyses each one fully. This module only
collects those results per channel and reference and writes, per level (all n, and every n):

    layers_<ch>_effects.csv            effect of every treatment at every layer (% of reference / % points)
    layers_<ch>_absolute.csv           value of every treatment at every layer
    layers_<ch>_concordance.csv        per treatment x reference: the call at every layer, agreement, and which
                                       step changes the conclusion (plain language)
    layers_<ch>_vs_<ref>.png/.pdf      ladder plot: rows = treatments, columns = layers, effect with its interval
    layers_<ch>_absolute.png/.pdf      the same grid with absolute values

All n: model contrasts (95% CI, adjusted p, stars; with one n the estimates only, "needs >= 2 n").
Per n: the within-n effect with the technical 95% interval from the fields (descriptive; no p-values).
"""
from __future__ import annotations

import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import models

LADDER = [  # (key, short label, what this step changes compared with the previous layer)
    ("field_raw", "1 Field, raw", ""),
    ("field_corrected", "2 Field, corrected", "flat-field and background correction of the whole field"),
    ("nucleus_mean_raw", "3 Nucleus, raw", "measuring per nucleus instead of per field"),
    ("nucleus_mean_bgsub", "4 Nucleus, corrected", "nucleus correction: shrunk outline, local background, flat field"),
    ("frac_positive", "5 % positive", "the positivity cut-off"),
]
EXTRAS = [("stained_area_mean", "Stained-area MFI"), ("field_total", "Total stain / field")]
SIG = 0.05
SAME_PCT, SAME_PP = 2.0, 1.0  # |effect| below 2 % (intensities) / 1 percentage point and not significant: "no change"
TECH = "technical 95% interval from the fields (descriptive)"
REF_COLOUR, FRAC_COLOUR, EXTRA_COLOUR = "#2b5d8a", "#b5532c", "#7a7a7a"


def column(ch: str, key: str) -> str:
    return f"readout_frac_positive_{ch}" if key == "frac_positive" else f"intensity_{ch}_{key}"


def channels(metric_ids) -> list[str]:
    """Channels (slugs) with at least two ladder layers among the metric ids."""
    ids = set(metric_ids)
    out = []
    for m in sorted(ids):
        g = re.fullmatch(r"intensity_(.+)_field_raw", m) or re.fullmatch(r"readout_frac_positive_(?!p90_)(.+)", m)
        if g and g.group(1) not in out:
            ch = g.group(1)
            if sum(column(ch, k) in ids for k, _, _ in LADDER) >= 2:
                out.append(ch)
    return sorted(set(out))


def layer_list(ch: str, metric_ids) -> list[dict]:
    ids = set(metric_ids)
    out = []
    for i, (k, lab, step) in enumerate(LADDER, 1):
        if column(ch, k) in ids:
            out.append(dict(key=k, label=lab, index=i, step=step, metric=column(ch, k), extra=False,
                            frac=k == "frac_positive"))
    for k, lab in EXTRAS:
        if column(ch, k) in ids:
            out.append(dict(key=k, label=lab, index=0, step="", metric=column(ch, k), extra=True, frac=False))
    return out


# ---------------------------------------------------------------- effects and values from the model results

def _call(est: float, lo: float, hi: float, frac: bool, sig: bool | None) -> tuple[str, str]:
    """(direction, call). direction: higher | lower | no change; sig None = significance not available."""
    if not np.isfinite(est):
        return "", "n/a"
    null = 0.0 if frac else 100.0
    small = abs(est - null) < (SAME_PP if frac else SAME_PCT)
    d = "higher" if est > null else "lower"
    if small and not sig:
        d = "no change"
    if sig is None:
        return d, d
    return d, f"{d} *" if sig else (f"{d} (ns)" if d != "no change" else d)


def effects(res, layer: dict, ref: str, level: str, cell_line=None) -> pd.DataFrame:
    """Effect rows of one layer vs `ref` from a ModelResult: all n (contrasts) or every n (per_n, technical)."""
    rows = []
    frac = layer["frac"]
    base = dict(layer=layer["label"], layer_key=layer["key"], layer_index=layer["index"], extra=layer["extra"],
                metric=layer["metric"], reference=ref, cell_line=cell_line,
                effect_type="pp_difference" if frac else "pct_of_reference")
    if level == "all_n":
        c = res.contrasts
        if c is None or not len(c):
            return pd.DataFrame()
        fam = c[(c["family_kind"] == "vs_control") & (c["comparator"] == ref)]
        for r in fam.itertuples():
            if frac:
                est, lo, hi = r.pp_difference, r.pp_ci_low, r.pp_ci_high
            else:
                est, lo, hi = 100 * r.estimate, 100 * r.ci_low, 100 * r.ci_high
            has_ci = bool(np.isfinite(lo) and np.isfinite(hi))
            sig = (bool(np.isfinite(r.p_adj) and r.p_adj < SIG)) if has_ci else None
            d, call = _call(est, lo, hi, frac, sig)
            rows.append(dict(base, level="all_n", n=None, treatment=r.group, estimate=est, ci_low=lo, ci_high=hi,
                             interval="95% CI (mixed model across n)" if has_ci else models.NEEDS_2N,
                             p_adj=r.p_adj, significance=r.significance if isinstance(r.significance, str) else "",
                             N=r.N, direction=d, call=call))
    else:
        pn = res.per_n
        if pn is None or not len(pn):
            return pd.DataFrame()
        col = "pp_difference" if frac else ("ratio_of_means" if res.unit == "well" else "ratio")
        for r in pn[pn["control"] == ref].itertuples():
            v = getattr(r, col)
            k = 1.0 if frac else 100.0
            est, lo, hi = k * v, k * r.tech_ci_low, k * r.tech_ci_high
            has = bool(np.isfinite(lo) and np.isfinite(hi))
            null = 0.0 if frac else 100.0
            sig = (lo > null or hi < null) if has else None
            d, call = _call(est, lo, hi, frac, sig)
            if sig is not None:
                call = call.replace(" *", " (tech.)")
            rows.append(dict(base, level="n", n=r.n, treatment=r.treatment, estimate=est, ci_low=lo, ci_high=hi,
                             interval=TECH if has else "no interval (too few fields)", p_adj=np.nan, significance="",
                             N=r.fields, direction=d, call=call))
    return pd.DataFrame(rows)


def absolute(res, layer: dict, level: str, cell_line=None) -> pd.DataFrame:
    """Value of every treatment at one layer: all n = model means with 95% CI (one n: estimate only); per n = mean of
    the fields (pooled % for fractions) with the technical 95% interval."""
    rows = []
    k = 100.0 if layer["frac"] else 1.0
    base = dict(layer=layer["label"], layer_key=layer["key"], layer_index=layer["index"], extra=layer["extra"],
                metric=layer["metric"], unit="% of nuclei" if layer["frac"] else getattr(res.metric, "unit", ""),
                cell_line=cell_line)
    if level == "all_n":
        for r in res.emm.itertuples():
            rows.append(dict(base, level="all_n", n=None, treatment=r.treatment, value=k * r.estimate,
                             ci_low=k * r.ci_low, ci_high=k * r.ci_high,
                             interval="95% CI (mixed model across n)" if np.isfinite(r.ci_low) else models.NEEDS_2N))
    elif res.prepared is not None:
        fs = models._field_summary(res.prepared)
        for (n, t), r in fs.iterrows():
            v, se = float(r["value"]), float(r["se"])
            if layer["frac"]:
                lo, hi = max(v - 1.96 * se, 0.0), min(v + 1.96 * se, 1.0)
            else:  # se is relative (SE / mean) for intensities
                lo, hi = v * np.exp(-1.96 * se), v * np.exp(1.96 * se)
            rows.append(dict(base, level="n", n=n, treatment=t, value=k * v, ci_low=k * lo, ci_high=k * hi,
                             interval=TECH))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- concordance

def _sig_word(call: str) -> str:
    return "significant" if call.endswith("*") else ("technical interval excludes no effect" if "(tech.)" in call
                                                       else "not significant" if "(ns)" in call else "")


def concordance(eff: pd.DataFrame, level_note: str = "") -> pd.DataFrame:
    """One row per (n, reference, treatment): the call at every ladder layer, whether the layers agree, at which
    steps the conclusion changes, and a plain-language summary. Extras are listed but not part of the agreement."""
    rows = []
    if eff.empty:
        return pd.DataFrame()
    keys = ["n", "reference", "treatment"]
    for key, g in eff.assign(n=eff["n"].fillna("")).groupby(keys, sort=False):
        n, ref, tr = key
        lad = g[~g["extra"]].sort_values("layer_index")
        ext = g[g["extra"]]
        calls = list(zip(lad["layer_index"], lad["layer"], lad["call"], lad["direction"]))
        valid = [c for c in calls if c[2] != "n/a"]
        dirs = {c[3] for c in valid}
        has_sig = any(("*" in c[2]) or ("(ns)" in c[2]) or ("(tech.)" in c[2]) for c in valid)
        changes = []
        for (i0, l0, c0, _), (i1, l1, c1, _) in zip(valid, valid[1:]):
            if c0 != c1:
                step = " and ".join(s for k_, _, s in LADDER[i0:i1] if s)
                changes.append(dict(step=f"{i0}->{i1}", what=step, before=f"{l0}: {c0}", after=f"{l1}: {c1}"))
        agree_dir = not {"higher", "lower"} <= dirs  # "no change" next to a direction is not a flip
        agree = not changes
        if not valid:
            text = "No layer has an estimate."
        elif agree:
            d = valid[0][3]
            if d == "no change":
                what = f"no clear change vs {ref}"
            else:
                what = f"{d} than {ref}"
            if has_sig:
                sw = _sig_word(valid[0][2])
                text = f"All {len(valid)} layers agree: {what}" + (f" ({sw} in every layer)." if sw else ".")
            else:
                text = f"All {len(valid)} layers point the same way: {what}" + (
                    f" ({level_note})." if level_note else ".")
        else:
            parts = [f"step {c['step']}, {c['what']} ({c['before']} -> {c['after']})" for c in changes]
            lead = ("No two layers point in opposite directions, but the call changes" if agree_dir
                    else "The direction flips between layers")
            text = f"{lead}: " + "; ".join(parts) + "."
            if any(c["step"].startswith("2->") and c["step"].endswith("3") for c in changes):
                text += (" Field layers also scale with how many cells are in the field, nuclear layers do not: "
                         "check the count and the density scatter.")
            head = valid[-1]
            text += f" Headline ({head[1]}): {head[2]}."
        if len(ext):
            text += " Extras: " + ", ".join(f"{r.layer} {r.call}" for r in ext.itertuples()) + "."
        row = dict(level=g["level"].iloc[0], n=n or None, cell_line=g["cell_line"].iloc[0], reference=ref,
                   treatment=tr, layers=len(valid), direction_agrees=agree_dir, conclusion_agrees=agree,
                   changes_at="; ".join(c["step"] for c in changes),
                   changing_steps="; ".join(c["what"] for c in changes), summary=text)
        for r in g.itertuples():
            row[f"{r.layer_key}"] = r.call
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- plots

def _fmt(v: float, frac: bool) -> str:
    if not np.isfinite(v):
        return ""
    return f"{v:+.1f} pp" if frac else f"{v:.0f}%"


def _fmt_abs(v: float) -> str:
    if not np.isfinite(v):
        return ""
    a = abs(v)
    if a >= 1e6:
        return f"{v:.2e}".replace("e+0", "e").replace("e+", "e")
    return f"{v:,.0f}" if a >= 100 else f"{v:.3g}"


def ladder_plot(eff: pd.DataFrame, layers: list[dict], order: list[str], out: Path, title: str, footer: str,
                level: str, ncol: dict | None = None):
    """Rows = treatments, columns = layers: effect vs the reference with its interval. Intensity layers share one
    x axis (% of reference) so compression between layers is visible; % positive is in percentage points."""
    trs = [t for t in order if t in set(eff["treatment"])]
    if not trs or not layers:
        return
    per_n = level == "n" and eff["n"].nunique() > 1
    ncols = len(layers)
    h = 1.25 + 0.42 * len(trs) + 0.95
    fig, axes = plt.subplots(1, ncols, figsize=(1.0 + 2.05 * ncols, h), sharey=True, squeeze=False)
    axes = axes[0]
    ints = eff[~eff["layer_key"].eq("frac_positive")]
    vals = pd.concat([ints["estimate"], ints["ci_low"], ints["ci_high"]]).dropna()
    lo_i = min(float(vals.min()) if len(vals) else 90, 95)
    hi_i = max(float(vals.max()) if len(vals) else 110, 105)
    pad = 0.08 * (hi_i - lo_i)
    y = np.arange(len(trs))[::-1]
    ypos = dict(zip(trs, y))
    for ax, L in zip(axes, layers):
        g = eff[eff["layer_key"] == L["key"]]
        frac = L["frac"]
        colour = FRAC_COLOUR if frac else (EXTRA_COLOUR if L["extra"] else REF_COLOUR)
        if L["extra"]:
            ax.set_facecolor("#f6f6f6")
        null = 0.0 if frac else 100.0
        ax.axvline(null, color="#999999", lw=0.8, zorder=0)
        if frac:
            fv = pd.concat([g["estimate"], g["ci_low"], g["ci_high"]]).dropna()
            a_, b_ = (min(float(fv.min()), -2), max(float(fv.max()), 2)) if len(fv) else (-5, 5)
            p_ = 0.12 * (b_ - a_)
            ax.set_xlim(a_ - p_, b_ + p_)
        else:
            ax.set_xlim(lo_i - pad, hi_i + pad)
        for r in g.itertuples():
            if r.treatment not in ypos:
                continue
            yy = ypos[r.treatment]
            if per_n:  # several n's in one panel: small offsets, colour = n
                ns = sorted(eff["n"].dropna().unique())
                yy = yy + 0.25 - 0.5 * ns.index(r.n) / max(len(ns) - 1, 1)
                c = (ncol or {}).get(r.n, colour)
            else:
                c = colour
            if np.isfinite(r.ci_low) and np.isfinite(r.ci_high):
                ax.plot([r.ci_low, r.ci_high], [yy, yy], color=c, lw=1.6 if level == "all_n" else 1.0,
                        solid_capstyle="butt")
            ax.plot([r.estimate], [yy], "o", color=c, ms=5.5 if level == "all_n" else 4.5,
                    mfc=c if (r.call.endswith("*") or "(tech.)" in r.call) else "white", mew=1.3, zorder=3)
            if not per_n:
                star = r.significance if isinstance(r.significance, str) and r.significance not in ("", "ns") else ""
                txt = _fmt(r.estimate, frac) + (f" {star}" if star else "")
                ax.annotate(txt, (1.0, yy), xycoords=("axes fraction", "data"), xytext=(-2, 5),
                            textcoords="offset points", ha="right", va="bottom", fontsize=6.5, color="#333333")
        unit = "% points" if frac else "% of ref."
        ax.set_title(("Extra: " if L["extra"] else "") + L["label"] + f"\n({unit})", fontsize=8,
                     color=colour if not L["extra"] else "#444444")
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color="#eeeeee", lw=0.6)
        ax.set_ylim(-0.7, len(trs) - 0.3)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(trs, fontsize=8)
    if per_n and ncol:
        from matplotlib.lines import Line2D

        hs = [Line2D([], [], marker="o", ls="", color=ncol.get(n, "#333"), label=str(n).split("_")[0])
              for n in sorted(eff["n"].dropna().unique())]
        fig.legend(handles=hs, loc="upper right", fontsize=7, frameon=False, ncol=min(len(hs), 6))
    fig.suptitle(title, fontsize=10, x=0.01, ha="left")
    _save(fig, out, footer, h)


def absolute_plot(ab: pd.DataFrame, layers: list[dict], order: list[str], ref: str | None, out: Path, title: str,
                  footer: str):
    trs = [t for t in order if t in set(ab["treatment"])]
    if not trs or not layers:
        return
    ncols = len(layers)
    h = 1.25 + 0.42 * len(trs) + 0.95
    fig, axes = plt.subplots(1, ncols, figsize=(1.0 + 2.05 * ncols, h), sharey=True, squeeze=False)
    axes = axes[0]
    y = np.arange(len(trs))[::-1]
    ypos = dict(zip(trs, y))
    for ax, L in zip(axes, layers):
        g = ab[ab["layer_key"] == L["key"]]
        colour = FRAC_COLOUR if L["frac"] else (EXTRA_COLOUR if L["extra"] else REF_COLOUR)
        if L["extra"]:
            ax.set_facecolor("#f6f6f6")
        if ref is not None and (g["treatment"] == ref).any():
            rv = float(g.loc[g["treatment"] == ref, "value"].iloc[0])
            if np.isfinite(rv):
                ax.axvline(rv, color="#999999", lw=0.8, ls="--", zorder=0)
        v = pd.concat([g["value"], g["ci_low"], g["ci_high"]]).dropna()
        if len(v):
            a_, b_ = float(v.min()), float(v.max())
            p_ = 0.1 * (b_ - a_) if b_ > a_ else max(abs(a_) * 0.1, 1)
            ax.set_xlim(min(a_ - p_, 0 if a_ >= 0 and a_ < 0.5 * b_ else a_ - p_), b_ + 2.2 * p_)
        for r in g.itertuples():
            if r.treatment not in ypos:
                continue
            yy = ypos[r.treatment]
            if np.isfinite(r.ci_low) and np.isfinite(r.ci_high):
                ax.plot([r.ci_low, r.ci_high], [yy, yy], color=colour, lw=1.4)
            ax.plot([r.value], [yy], "o", color=colour, ms=5, mfc="white" if r.treatment != ref else colour, mew=1.3)
            ax.annotate(_fmt_abs(r.value), (1.0, yy), xycoords=("axes fraction", "data"), xytext=(-2, 5),
                        textcoords="offset points", ha="right", va="bottom", fontsize=6.5, color="#333333")
        unit = g["unit"].iloc[0] if len(g) else ""
        ax.set_title(("Extra: " if L["extra"] else "") + L["label"] + "\n" + "\n".join(textwrap.wrap(f"({unit})", 30)),
                     fontsize=8, color=colour if not L["extra"] else "#444444")
        ax.tick_params(labelsize=7)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: _fmt_abs(x)))
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color="#eeeeee", lw=0.6)
        ax.set_ylim(-0.7, len(trs) - 0.3)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(trs, fontsize=8)
    fig.suptitle(title, fontsize=10, x=0.01, ha="left")
    _save(fig, out, footer, h)


def _save(fig, out: Path, footer: str, height_in: float):
    width = fig.get_size_inches()[0]
    lines = textwrap.wrap(footer, width=int(width * 15.5)) if footer else []
    foot_h = (0.12 * len(lines) + 0.08) / height_in
    if lines:
        fig.text(0.01, 0.01, "\n".join(lines), fontsize=6.5, color="#444444", ha="left", va="bottom")
    fig.tight_layout(rect=(0, foot_h, 1, 0.97))
    fig.savefig(out.with_suffix(".png"), dpi=200)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- writing one level

FOOT_COMMON = ("Raw layers (1, 3) compress effects toward 100%: a fixed baseline (camera offset, background, "
               "autofluorescence) is in both treatment and reference. Field layers (1, 2) and the total also scale "
               "with how many cells are in the field: compare with the density scatter ({density}). Filled dot = "
               "{sigdef}; open = not. Layer 5 is in percentage points (pp).")


def write_level(where: Path, ch: str, layers: list[dict], results: dict, refs: list[str], order: list[str],
                level: str, title_prefix: str, cell_line=None, n_filter: str | None = None, ncol: dict | None = None,
                density_hint: str = "") -> list[Path]:
    """Write the layer outputs of one channel at one level into `where`. results: {metric id: ModelResult} (the
    all-n block's results; per-n effects come from their per_n tables). n_filter: the n for per-n output."""
    eff_l, abs_l = [], []
    for L in layers:
        res = results.get(L["metric"])
        if res is None:
            continue
        for ref in refs:
            e = effects(res, L, ref, level, cell_line)
            if len(e) and n_filter is not None:
                e = e[e["n"] == n_filter]
            if len(e):
                eff_l.append(e)
        a = absolute(res, L, level, cell_line) if (res.has_estimates or level == "n") else pd.DataFrame()
        if len(a) and n_filter is not None:
            a = a[a["n"] == n_filter]
        if len(a):
            abs_l.append(a)
    if not eff_l:
        return []
    eff = pd.concat(eff_l, ignore_index=True)
    ab = pd.concat(abs_l, ignore_index=True) if abs_l else pd.DataFrame()
    one_n = level == "all_n" and eff["interval"].eq(models.NEEDS_2N).all()
    note = f"significance {models.NEEDS_2N}" if one_n else ("technical intervals only, descriptive"
                                                          if level == "n" else "")
    conc = concordance(eff, note)
    written = []
    p = where / f"layers_{ch}_effects.csv"
    eff.to_csv(p, index=False)
    written.append(p)
    if len(ab):
        p = where / f"layers_{ch}_absolute.csv"
        ab.to_csv(p, index=False)
        written.append(p)
    p = where / f"layers_{ch}_concordance.csv"
    conc.to_csv(p, index=False)
    written.append(p)
    present = [L for L in layers if L["key"] in set(eff["layer_key"])]
    dens = density_hint or "the density plots"
    if level == "all_n":
        sigdef = ("significant; with one n every dot is open (significance " + models.NEEDS_2N + ")"
                  if one_n else "adjusted p < 0.05 (mixed model, vs this reference)")
        how = ("Only one n: dots are the estimates from that n; across-n intervals, p-values and stars "
               f"{models.NEEDS_2N}; technical intervals are in the per-n folder." if one_n else
               "Lines: 95% CI of the mixed model across n; stars: adjusted p within the family vs this reference.")
    else:
        sigdef = "the technical interval excludes no effect (descriptive)"
        how = ("Within one n: lines are technical 95% intervals from the fields. Descriptive only: fields are "
               "technical replicates, not independent experiments; no p-values.")
    for ref in refs:
        e = eff[eff["reference"] == ref]
        if e.empty:
            continue
        c = conc[conc["reference"] == ref]
        changed = c[~c["conclusion_agrees"]]
        verdict = ("Layers agree for every treatment." if changed.empty else
                   "Conclusion changes between layers for: " + ", ".join(
                       f"{r.treatment} (step {r.changes_at})" for r in changed.itertuples()) + ".")
        foot = (f"Effect of each treatment vs {ref} at every layer of the {ch} analysis. {how} "
                + FOOT_COMMON.format(density=dens, sigdef=sigdef) + f" {verdict} Details: layers_{ch}_concordance.csv.")
        tr_order = [t for t in order if t != ref]
        p = where / f"layers_{ch}_vs_{re.sub(r'[^\w.\-]+', '_', ref).strip('_')}"
        ladder_plot(e, present, tr_order, p,
                    f"{title_prefix}: {ch} effect vs {ref} at every analysis layer"
                    + (" (descriptive)" if level == "n" else " (one n: estimates only)" if one_n else ""),
                    foot, level, ncol)
        written.append(p.with_suffix(".png"))
    if len(ab):
        how_a = ("Model estimates per treatment" + (" (one n: estimate only, no CI)" if one_n else " with 95% CI")
                 if level == "all_n" else "Mean of the fields (pooled % for layer 5) with technical 95% intervals "
                 "(descriptive)")
        foot = (f"Absolute {ch} values at every layer. {how_a}. Dashed line and filled dot = the main reference "
                f"({refs[0] if refs else '-'}). Units differ between layers (raw / corrected intensity units, % of "
                "nuclei, integrated signal); compare treatments within a column, not across columns.")
        p = where / f"layers_{ch}_absolute"
        absolute_plot(ab, present, order, refs[0] if refs else None, p,
                      f"{title_prefix}: {ch} absolute values at every analysis layer", foot)
        written.append(p.with_suffix(".png"))
    return written
