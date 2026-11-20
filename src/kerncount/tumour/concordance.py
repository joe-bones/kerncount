"""Method concordance (ADR 0002, ticket 07): primary method vs the Fiji legacy readout, and the threshold sweep.

Readout-generic: everything runs over the readouts declared **confirmatory** in settings; the Fiji counterpart of a
readout is `config.legacy_counterpart` (a `legacy = "<column>"` key on the readout, else the default map), so a
change of the confirmatory readout definition needs no code change here.

- Per-tumour scatter new vs Fiji (log axes, identity line) per confirmatory readout with a Fiji counterpart.
- Side-by-side effects: the pre-registered statistics (stats.readout_results, same families, Holm) re-run on the
  Fiji column; per contrast: direction agreement (same side of 1) and significance agreement (Holm p < 0.05 in both
  or in neither); summary per readout x family. The Fiji results are legacy / exploratory: never starred anywhere.
- Threshold sweep: each confirmatory effect at threshold x scale (sweep.tables), plotted against the scale.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import stats
from .config import Factor, Family, Readout, legacy_counterpart, readouts, tumour_settings

ALPHA = 0.05
LEGACY = "Fiji legacy (exploratory)"
COMPARE_COLUMNS = ["readout", "fiji_column", "family", "family_role", "comparison", "estimate", "ci_low", "ci_high",
                   "p_adj", "fiji_estimate", "fiji_ci_low", "fiji_ci_high", "fiji_p_adj", "direction_agree",
                   "significance_agree", "note"]
SUMMARY_COLUMNS = ["readout", "family", "family_role", "contrasts", "direction_agree", "tested_in_both",
                   "significance_agree", "conclusion_differs"]
SWEEP_COLUMNS = ["readout", "scale", "family", "family_role", "comparison", "estimate", "ci_low", "ci_high", "p_adj",
                 "significant", "N"]


def _sign(x: float) -> int:
    return 0 if not np.isfinite(x) or x == 1 else (1 if x > 1 else -1)


def _yes_no(b: bool | None) -> str:
    return "" if b is None else ("yes" if b else "no")


def compare(new: pd.DataFrame, fiji: pd.DataFrame, fiji_column: str = "") -> pd.DataFrame:
    """Per contrast: both methods' estimates + Holm p, direction and significance agreement.
    `new` / `fiji` are stats.readout_results frames of the same families (matched on family + comparison)."""
    f = fiji.set_index(["family", "comparison"])
    rows = []
    for _, r in new.iterrows():
        k = (r["family"], r["comparison"])
        g = f.loc[k] if k in f.index else None
        fe = float(g["estimate"]) if g is not None else np.nan
        fpa = float(g["p_adj"]) if g is not None else np.nan
        direction = None if not (np.isfinite(r["estimate"]) and np.isfinite(fe)) else _sign(r["estimate"]) == _sign(fe)
        tested = np.isfinite(r["p_adj"]) and np.isfinite(fpa)
        sig = None
        if tested:
            sn, sf = r["p_adj"] < ALPHA, fpa < ALPHA
            sig = sn == sf and (not sn or _sign(r["estimate"]) == _sign(fe))
        note = "" if tested else ("not tested (" + (r["note"] or "see statistics") + ")"
                                  if not np.isfinite(r["p_adj"]) else "Fiji side not tested")
        rows.append(dict(readout=r["readout"], fiji_column=fiji_column, family=r["family"],
                         family_role=r["family_role"], comparison=r["comparison"], estimate=r["estimate"],
                         ci_low=r["ci_low"], ci_high=r["ci_high"], p_adj=r["p_adj"], fiji_estimate=fe,
                         fiji_ci_low=float(g["ci_low"]) if g is not None else np.nan,
                         fiji_ci_high=float(g["ci_high"]) if g is not None else np.nan, fiji_p_adj=fpa,
                         direction_agree=_yes_no(direction), significance_agree=_yes_no(sig), note=note))
    return pd.DataFrame(rows, columns=COMPARE_COLUMNS)


def summary(cmp: pd.DataFrame) -> pd.DataFrame:
    """Per readout x family: '<agree>/<n>' counts and whether any tested conclusion differs."""
    rows = []
    for (ro, fam), g in cmp.groupby(["readout", "family"], sort=False):
        d = g[g["direction_agree"] != ""]
        s = g[g["significance_agree"] != ""]
        rows.append(dict(readout=ro, family=fam, family_role=g["family_role"].iloc[0], contrasts=len(g),
                         direction_agree=f"{(d['direction_agree'] == 'yes').sum()}/{len(d)}",
                         tested_in_both=len(s), significance_agree=f"{(s['significance_agree'] == 'yes').sum()}/{len(s)}",
                         conclusion_differs="yes" if (s["significance_agree"] == "no").any() else "no"))
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def sweep_results(tables: dict[float, pd.DataFrame], facs: list[Factor], fams: list[Family],
                  ros: list[Readout]) -> pd.DataFrame:
    """Every contrast of every family per readout at each threshold scale (the same statistics as the primary)."""
    rows = []
    for scale, t in tables.items():
        for ro in ros:
            if ro.column not in t.columns:
                continue
            res = stats.readout_results(t, facs, fams, ro)
            for _, r in res.iterrows():
                rows.append(dict(readout=ro.name, scale=scale, family=r["family"], family_role=r["family_role"],
                                 comparison=r["comparison"], estimate=r["estimate"], ci_low=r["ci_low"],
                                 ci_high=r["ci_high"], p_adj=r["p_adj"],
                                 significant=_yes_no(bool(r["p_adj"] < ALPHA)) if np.isfinite(r["p_adj"]) else "",
                                 N=r["N"]))
    return pd.DataFrame(rows, columns=SWEEP_COLUMNS)


# ---------------------------------------------------------------- plots


def scatter_plot(table: pd.DataFrame, facs: list[Factor], new_col: str, fiji_col: str, label: str, out: Path,
                 x_label: str | None = None, x_name: str = "Fiji") -> Path:
    """Per tumour: primary (y) vs Fiji legacy (x; or another method: `x_label` axis label, `x_name` in the title),
    log axes, identity line, coloured by the x-role factor."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .plots import INK, MUTED, SERIES, _plain_log, _save

    x = pd.to_numeric(table[fiji_col], errors="coerce")
    y = pd.to_numeric(table[new_col], errors="coerce")
    ok = (x > 0) & (y > 0)
    colour_by = next((f for f in facs if f.role == "x"), facs[0] if facs else None)
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    if ok.any():
        lo, hi = float(min(x[ok].min(), y[ok].min())) / 1.5, float(max(x[ok].max(), y[ok].max())) * 1.5
        ax.plot([lo, hi], [lo, hi], color=MUTED, lw=1, ls="--", label="identity (y = x)", zorder=1)
        groups = colour_by.levels if colour_by is not None else (None,)
        for i, lv in enumerate(groups):
            sel = ok & ((table[colour_by.name].astype(str) == lv) if colour_by is not None else True)
            if sel.any():
                ax.scatter(x[sel], y[sel], s=26, color=SERIES[i % len(SERIES)], edgecolor="white", linewidth=0.5,
                           label=str(lv) if lv is not None else None, zorder=3)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        _plain_log(ax, "x")
        _plain_log(ax, "y")
        r = float(np.exp(np.median(np.log(y[ok] / x[ok]))))
        rho = float(pd.Series(np.log(x[ok])).corr(pd.Series(np.log(y[ok])), method="spearman")) if ok.sum() > 2 \
            else float("nan")
        ax.set_title(f"median new / {x_name} = {r:.2f}; Spearman ρ = {rho:.2f}; n = {int(ok.sum())}", fontsize=9,
                     color=INK)
        ax.legend(fontsize=8, frameon=False, loc="upper left")
    else:
        ax.text(0.5, 0.5, "no tumour with both values", ha="center", va="center", transform=ax.transAxes)
    ax.set_xlabel(x_label or f"{LEGACY}: {fiji_col}")
    ax.set_ylabel(f"primary method: {label}")
    fig.tight_layout()
    return _save(fig, out)


def sweep_plot(sw: pd.DataFrame, readout: Readout, out: Path) -> Path:
    """One panel per family: each contrast's fold-change (log y, 95 % CI) vs threshold scale; filled = Holm p < 0.05."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .plots import MUTED, SERIES, _plain_log, _save, _wrap

    d = sw[sw["readout"] == readout.name]
    # only families with an estimate (KO-only data: F4-F7 'need WT data' would be empty panels)
    est = pd.to_numeric(d["estimate"], errors="coerce")
    fams = list(dict.fromkeys(d.loc[np.isfinite(est), "family"])) or list(dict.fromkeys(d["family"]))[:1]
    fig, axes = plt.subplots(1, max(len(fams), 1), figsize=(4.6 * max(len(fams), 1), 4.4), squeeze=False)
    for ax, fam in zip(axes[0], fams):
        g = d[d["family"] == fam]
        for i, (comp, c) in enumerate(g.groupby("comparison", sort=False)):
            c = c.sort_values("scale")
            col = SERIES[i % len(SERIES)]
            ok = np.isfinite(c["estimate"].to_numpy(float))
            if not ok.any():
                continue
            ax.plot(c["scale"], c["estimate"], color=col, lw=1.2, label=comp)
            ax.fill_between(c["scale"], c["ci_low"], c["ci_high"], color=col, alpha=0.08, lw=0)
            sig = (c["significant"] == "yes").to_numpy()
            ax.scatter(c["scale"][sig], c["estimate"][sig], color=col, s=22, zorder=3)
            ax.scatter(c["scale"][~sig], c["estimate"][~sig], facecolor="white", edgecolor=col, s=22, zorder=3)
        ax.axhline(1, color=MUTED, lw=1, ls="--")
        ax.axvline(1, color=MUTED, lw=0.8, ls=":")
        ax.set_yscale("log")
        _plain_log(ax, "y")
        ax.set_xlabel("threshold × scale (1 = primary)")
        ax.set_title(f"{fam} [{g['family_role'].iloc[0]}]", fontsize=10)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=6, frameon=False, loc="best")
    axes[0][0].set_ylabel(f"fold-change, {readout.label}")
    fig.suptitle(_wrap(f"Threshold sweep: {readout.label} — filled = Holm p < 0.05, band = 95 % CI",
                       fig.get_figwidth() * 0.95, 10), fontsize=10)
    fig.tight_layout()
    return _save(fig, out)


# ---------------------------------------------------------------- outputs


def analyse(out: Path, table: pd.DataFrame, settings: dict, facs: list[Factor], fams: list[Family],
            report: Path | None = None) -> list[str]:
    """Concordance + threshold-sweep outputs for the confirmatory readouts; returns the written files.
    Writes concordance.csv, threshold_sweep.csv, plots, a 'Concordance' / 'Threshold sweep' sheet in results.xlsx
    (when it exists) and a report section appended to `report`."""
    ts = tumour_settings(settings)
    conf = [r for r in readouts(settings) if r.role == "confirmatory" and r.column in table.columns]
    written: list[str] = []
    md = ["## Method concordance: primary vs Fiji legacy (ADR 0002)", "",
          "**Fiji volume mostly reflects the acquired stack depth, not tumour size** (ADR 0002 amendment "
          "2026-11-20): only ~80-150 µm of each tumour is imaged, and the lab's per-layer Triangle threshold turns "
          "faint haze in near-empty layers into tumour, so Fiji's main object spans every imaged layer and its "
          "volume is about footprint x acquired depth (an acquisition choice). Compare the mid-section area and "
          "perimeter below, not volumes.", ""]
    cmps = []
    for ro in conf:
        fcol = legacy_counterpart(settings, ro)
        if not fcol or fcol not in table.columns or not pd.to_numeric(table[fcol], errors="coerce").gt(0).any():
            md.append(f"- {ro.name}: no Fiji counterpart in tumours.csv ({fcol or 'none declared'}); not compared.")
            continue
        written.append(str(scatter_plot(table, facs, ro.column, fcol, ro.label,
                                        out / "plots" / f"{ro.name}_vs_fiji")))
        new = stats.readout_results(table, facs, fams, ro)
        legacy = stats.readout_results(table, facs, fams, Readout(f"{ro.name}_fiji", fcol,
                                                                  f"{LEGACY} {ro.label}", ro.role))
        cmps.append(compare(new, legacy, fcol))
    cmp = pd.concat(cmps, ignore_index=True) if cmps else pd.DataFrame(columns=COMPARE_COLUMNS)
    summ = summary(cmp) if len(cmp) else pd.DataFrame(columns=SUMMARY_COLUMNS)
    delim = settings.get("project", {}).get("csv_delimiter", ",")
    if len(cmp):
        cmp.to_csv(out / "concordance.csv", index=False, sep=delim)
        written.append(str(out / "concordance.csv"))
        md += ["Same families and Holm correction on both methods; the Fiji column is a legacy, exploratory readout "
               "(never starred). Direction = same side of 1; significance = Holm p < 0.05 in both or neither "
               "(and the same direction).", "",
               "| readout | family | contrasts | direction agrees | significance agrees (tested in both) | "
               "conclusion differs |", "|---|---|---|---|---|---|"]
        for _, s in summ.iterrows():
            md.append(f"| {s['readout']} | {s['family']} [{s['family_role']}] | {s['contrasts']} | "
                      f"{s['direction_agree']} | {s['significance_agree']} | {s['conclusion_differs']} |")
        differs = cmp[cmp["significance_agree"] == "no"]
        if len(differs):
            md += ["", "**The conclusion differs between the methods for:** " + "; ".join(
                f"{r['readout']} {r['family']} {r['comparison']} (new {r['estimate']:.2f}, p_Holm "
                f"{r['p_adj']:.3g}; Fiji {r['fiji_estimate']:.2f}, p_Holm {r['fiji_p_adj']:.3g})"
                for _, r in differs.iterrows())]
    md.append("")
    sw = pd.DataFrame(columns=SWEEP_COLUMNS)
    md += ["## Threshold sweep (ADR 0002)", ""]
    if ts["threshold_sweep"] and conf and "file" in table.columns:
        from . import sweep

        tabs = sweep.tables(out, table, settings)
        sw = sweep_results(tabs, facs, fams, conf)
        sw.to_csv(out / "threshold_sweep.csv", index=False, sep=delim)
        written.append(str(out / "threshold_sweep.csv"))
        for ro in conf:
            if ro.column not in tabs[min(tabs)].columns:
                md.append(f"- {ro.name}: not available in the sweep (not produced by the segmentation step)")
                continue
            written.append(str(sweep_plot(sw, ro, out / "plots" / f"{ro.name}_threshold_sweep")))
        md += _sweep_lines(sw)
    else:
        md.append("- not run (`[tumour] threshold_sweep = false`)." if not ts["threshold_sweep"]
                  else "- not run (no confirmatory readout / no file column).")
    _excel(out, cmp, summ, sw)
    if report is not None and report.exists():
        with report.open("a", encoding="utf-8") as fh:
            fh.write("\n" + "\n".join(md) + "\n")
    return written


def _sweep_lines(sw: pd.DataFrame) -> list[str]:
    """Per confirmatory-family contrast: range of the fold-change over the sweep and whether significance flips."""
    md = ["Fold-change of each confirmatory-family contrast over the threshold scales (primary = ×1); "
          "'stable' = Holm significance the same at every scale.", "",
          "| readout | family | comparison | fold-change range | significance |", "|---|---|---|---|---|"]
    d = sw[sw["family_role"] == "confirmatory"]
    for (ro, fam, comp), g in d.groupby(["readout", "family", "comparison"], sort=False):
        e = g["estimate"][np.isfinite(g["estimate"].astype(float))]
        if not len(e):
            continue
        flags = set(g["significant"]) - {""}
        sig = "not tested" if not flags else ("stable" if len(flags) == 1 else
                                              "changes: " + ", ".join(f"×{s:g} {v}" for s, v in
                                                                      zip(g["scale"], g["significant"])))
        md.append(f"| {ro} | {fam} | {comp.replace('|', chr(92) + '|')} | {e.min():.2f}–{e.max():.2f} | {sig} |")
    return md + [""]


def _excel(out: Path, cmp: pd.DataFrame, summ: pd.DataFrame, sw: pd.DataFrame) -> None:
    path = out / "results.xlsx"
    if not path.exists():
        return
    with pd.ExcelWriter(path, engine="openpyxl", mode="a", if_sheet_exists="replace") as xw:
        summ.to_excel(xw, sheet_name="Concordance summary", index=False)
        cmp.to_excel(xw, sheet_name="Concordance", index=False)
        sw.to_excel(xw, sheet_name="Threshold sweep", index=False)
