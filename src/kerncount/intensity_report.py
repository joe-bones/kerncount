"""Project-level intensity outputs: <project>/intensity/ (docs/intensity.md, "Outputs").

Every file is registered in intensity/manifest.json as {path, kind, title, description, n, channel} for the
browser app. Plots are descriptive; statistics live in the generations (metrics registry).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

# categorical slots in fixed order (one colour per condition, same across n's and plots)
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
EXTRA = ["#7a7a7a", "#a0522d", "#2f4f4f", "#b8860b"]
INK, MUTED, GRID = "#222222", "#666666", "#e6e6e6"
SCHEMA = "kerncount.intensity.manifest/1"


def _slug(s) -> str:
    import re

    return re.sub(r"[^\w.\-]+", "_", str(s)).strip("_")


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True})
    return plt


def treatment_colours(treatments: list[str]) -> dict:
    cols = PALETTE + EXTRA
    return {t: cols[i % len(cols)] for i, t in enumerate(treatments)}


class Manifest:
    def __init__(self, root: Path):
        self.root = root
        self.files: list[dict] = []

    def add(self, path: Path, kind: str, title: str, description: str, n=None, channel=None):
        self.files.append(dict(path=path.relative_to(self.root).as_posix(), kind=kind, title=title,
                               description=description, n=n, channel=channel))

    def write(self, extra: dict):
        doc = dict(schema=SCHEMA, built=time.strftime("%Y-%m-%dT%H:%M:%S"), **extra, files=self.files)
        (self.root / "manifest.json").write_text(json.dumps(doc, indent=2, default=_jsonable), encoding="utf-8")


def _jsonable(v):
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.bool_):
        return bool(v)
    return str(v)


def _clean(obj):
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    return obj


# ---------------------------------------------------------------- plots

def plot_flatfield(ff, title: str, out: Path, obj_before: dict, obj_after: dict):
    from . import flatfield as ff_mod

    plt = _plt()
    F = ff.flat()
    H, W = F.shape
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw=dict(width_ratios=[1.25, 1]))
    ax = axs[0]
    im = ax.imshow(F, cmap="Blues_r", vmin=min(F.min(), 0.8), vmax=max(F.max(), 1.05))
    cs = ax.contour(F, levels=np.round(np.arange(0.5, 1.21, 0.05), 2), colors=INK, linewidths=0.6)
    ax.clabel(cs, fontsize=7, fmt="%.2f")
    ax.set_title("Flat field F (1 = image average)", fontsize=9, color=INK)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(im, ax=ax, fraction=0.04)
    ax = axs[1]
    ax.plot(np.arange(W) / (W - 1), F[H // 2], color=PALETTE[0], lw=2, label="horizontal (through centre)")
    ax.plot(np.arange(H) / (H - 1), F[:, W // 2], color=PALETTE[1], lw=2, label="vertical (through centre)")
    ax.set_xlabel("position across the image (0 = left / top, 1 = right / bottom)")
    ax.set_ylabel("relative illumination")
    ax.legend(frameon=False, fontsize=8)
    mag = ff.info if ff.info else ff_mod.magnitude(ff)
    sub = ff_mod.describe(mag)
    if obj_before:
        sub += (f"\nCheck from the nuclei: edges {obj_before.get('edge_dim_pct', float('nan')):.0f}% dimmer before, "
                f"{obj_after.get('edge_dim_pct', float('nan')):.0f}% after correction")
    fig.suptitle(f"{title}\n{sub}", fontsize=9.5, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)


def plot_sweep(sw: pd.DataFrame, order: list[str], colours: dict, title: str, out: Path):
    plt = _plt()
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, measure, lab in zip(axs, ("mean", "p90"), ("mean in nucleus", "90th percentile in nucleus")):
        d = sw[sw["measure"] == measure]
        if d.empty:
            ax.set_visible(False)
            continue
        for t in order + ["all"]:
            g = d[d["treatment"] == t].sort_values("cutoff")
            if g.empty:
                continue
            ax.plot(g["cutoff"], 100 * g["frac_positive"], color=colours.get(t, INK) if t != "all" else INK,
                    lw=2.6 if t == "all" else 1.6, ls="--" if t == "all" else "-", label="all conditions" if t == "all" else t)
        chosen = d[np.isclose(d["factor"], 1.0)]["cutoff"]
        rule = d["rule"].iloc[0]
        if len(chosen):
            ax.axvline(chosen.iloc[0], color=MUTED, lw=1)
            ax.text(chosen.iloc[0], 2, f" cut-off used ({rule})", color=MUTED, fontsize=7.5, rotation=90, va="bottom")
        ax.set_xscale("log")
        ax.set_ylim(0, 102)
        ax.set_xlabel(f"cut-off: {lab}, minus local background (corrected units)" + ("\n(x-axis = median per-nucleus cut-off)" if rule == "lod" else ""))
        ax.set_title(f"% positive by {lab}", fontsize=9, color=INK)
    axs[0].set_ylabel("% positive nuclei")
    axs[0].legend(frameon=False, fontsize=7.5, ncol=1, loc="lower left")
    fig.suptitle(title, fontsize=9.5, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)


def plot_distribution(nuc: pd.DataFrame, order: list[str], colours: dict, cut: float, rule: str, title: str, out: Path):
    """Per-nucleus signal per condition on a log axis (violins + quartiles), and overlaid histograms (share of
    nuclei per log bin). Values at or below the floor (no detectable signal, incl. negatives) are piled at the
    floor and counted in the labels."""
    plt = _plt()
    d = nuc[nuc["measurable"]]
    present = [t for t in order if (d["treatment"] == t).any()]
    vals = [d.loc[d["treatment"] == t, "mean_bgsub_ff"].to_numpy(float) for t in present]
    allv = np.concatenate(vals) if vals else np.array([1.0])
    posv = allv[allv > 0]
    floor = max(float(np.percentile(posv, 0.5)) if len(posv) else 1.0, 1e-3)
    if np.isfinite(cut) and cut > 0:
        floor = min(floor, cut / 3)
    top = float(np.percentile(allv, 99.9)) if len(allv) else floor * 10
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw=dict(width_ratios=[1.3, 1]))
    ax = axs[0]
    logs = [np.log10(np.clip(v, floor, top)) for v in vals]
    parts = ax.violinplot(logs, positions=range(len(present)), showextrema=False, widths=0.85)
    for body, t in zip(parts["bodies"], present):
        body.set_facecolor(colours[t])
        body.set_edgecolor(colours[t])
        body.set_alpha(0.45)
    ytop = np.log10(top) + 0.08
    for i, (t, v) in enumerate(zip(present, vals)):
        q1, q2, q3 = np.log10(np.clip(np.percentile(v, [25, 50, 75]), floor, None))
        ax.plot([i, i], [q1, q3], color=INK, lw=3, solid_capstyle="round")
        ax.plot(i, q2, "o", color="white", ms=4, mec=INK, mew=1)
        pos = d.loc[d["treatment"] == t, "positive"]
        if pos.notna().any():
            ax.text(i, ytop, f"{100 * pos.astype(float).mean():.0f}%", ha="center", va="bottom", fontsize=7.5, color=INK)
    if np.isfinite(cut) and cut > 0:
        ax.axhline(np.log10(cut), color=MUTED, lw=1, ls="--")
        ax.text(len(present) - 0.5, np.log10(cut), f"cut-off ({rule}{', median' if rule == 'lod' else ''}) ",
                color=MUTED, fontsize=7.5, va="bottom", ha="right")
    ax.set_ylim(np.log10(floor) - 0.1, ytop + 0.25)
    ticks = [t for t in range(int(np.floor(np.log10(floor))), int(np.ceil(ytop)) + 1)]
    ax.set_yticks(ticks, [f"{10 ** t:g}" for t in ticks])
    ax.set_xticks(range(len(present)), present, rotation=30, ha="right")
    ax.set_ylabel("mean in nucleus - local background\n(flat-field corrected, log scale)")
    ax.set_title(f"Per nucleus, per condition (% positive on top; bottom = {floor:.3g} or less)", fontsize=9, color=INK)
    ax = axs[1]
    bins = np.geomspace(floor, top, 50)
    for t, v in zip(present, vals):
        below = float((v <= floor).mean())
        w = np.full(len(v), 1.0 / max(len(v), 1))
        ax.hist(v[v > floor], bins=bins, weights=w[v > floor], histtype="step", color=colours[t], lw=1.4,
                label=f"{t} ({100 * below:.0f}% at or below {floor:.3g})")
    if np.isfinite(cut) and cut > 0:
        ax.axvline(cut, color=MUTED, lw=1, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("mean - local background (log scale)")
    ax.set_ylabel("share of the condition's nuclei per bin")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    ax.set_title("Distribution", fontsize=9, color=INK)
    fig.suptitle(title, fontsize=9.5, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)


def plot_density(per_image: pd.DataFrame, order: list[str], colours: dict, xcol: str, xlab: str, ycol: str, title: str, out: Path):
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 4.6))
    for t in order:
        g = per_image[per_image["treatment"] == t]
        if g.empty:
            continue
        ax.scatter(g[xcol], 100 * g[ycol], s=36, color=colours[t], edgecolors="white", linewidths=1.2, label=t, zorder=3)
    ax.set_xlabel(xlab)
    ax.set_ylabel("% positive nuclei (per image)")
    ax.set_ylim(0, 102)
    ax.legend(frameon=False, fontsize=7.5, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.set_title(title, fontsize=9.5, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)


def plot_alignment(al: pd.DataFrame, warn: float, title: str, out: Path):
    plt = _plt()
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    ax.scatter(al["dx"], al["dy"], s=30, color=PALETTE[0], edgecolors="white", linewidths=1, zorder=3)
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(warn * np.cos(t), warn * np.sin(t), color=PALETTE[7], lw=1, ls="--", label=f"warning limit ({warn:g} px)")
    lim = max(warn * 1.5, float(np.nanmax(np.abs(al[["dx", "dy"]].to_numpy()))) * 1.2 if len(al) else warn)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("shift x (px)")
    ax.set_ylabel("shift y (px)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(title, fontsize=9, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------- write everything

def write(prj, nuc: pd.DataFrame, finfos: list, ffs: dict, cutoffs: pd.DataFrame, sweep: pd.DataFrame,
          names: list[str], chans: list[int], seg_ch: int, log=print) -> Path:
    from . import analysis, flatfield as ff_mod

    s = prj.settings
    root = prj.path / "intensity"
    root.mkdir(parents=True, exist_ok=True)
    for old in list(root.glob("*.png")) + list(root.glob("*.csv")) + list((root / "flatfield").glob("*")):
        old.unlink()
    man = Manifest(root)
    delim = s["project"].get("csv_delimiter", ",")
    order = analysis.treatment_order(sorted(nuc["treatment"].unique()), s["analysis"]["treatment_order"]) if len(nuc) else []
    colours = treatment_colours(order)

    # flat fields
    (root / "flatfield").mkdir(exist_ok=True)
    ff_rows = []
    for (n, ch), ff in sorted(ffs.items()):
        stem = f"{_slug(n)}__{_slug(ch)}"
        g = nuc[(nuc["n"] == n) & (nuc["channel"] == ch) & nuc["measurable"]]
        shape = ff.shape
        before = ff_mod.object_profile(g["mean_bgsub"].to_numpy(float), g["y"].to_numpy(float), g["x"].to_numpy(float), shape)
        after = ff_mod.object_profile(g["mean_bgsub_ff"].to_numpy(float), g["y"].to_numpy(float), g["x"].to_numpy(float), shape)
        mag = ff.info if ff.source != "off" else {}
        cr = cutoffs[(cutoffs["n"] == n) & (cutoffs["channel"] == ch) & (cutoffs["measure"] == "mean")] if len(cutoffs) else cutoffs
        ff_rows.append(dict(n=n, channel=ch, source=ff.source, images=ff.n_images,
                            edge_dim_pct=mag.get("edge_dim_pct"), darkest_side=mag.get("darkest_side"),
                            darkest_side_dim_pct=mag.get("darkest_side_dim_pct"), range_pct=mag.get("range_pct"),
                            nuclei_edge_dim_pct_before=before.get("edge_dim_pct"),
                            nuclei_edge_dim_pct_after=after.get("edge_dim_pct"),
                            frac_positive_without_flatfield=cr["frac_positive_without_flatfield"].iloc[0] if len(cr) else np.nan,
                            frac_positive_with_flatfield=cr["frac_positive"].iloc[0] if len(cr) else np.nan,
                            summary=ff_mod.describe(mag) if mag else "no correction"))
        if ff.source == "off":
            continue
        np.save(root / "flatfield" / f"{stem}.npy", ff.flat())
        man.add(root / "flatfield" / f"{stem}.npy", "flatfield_map", f"Flat field {n} / {ch} (array)",
                "The illumination map (float32, image size; 1 = image average). Intensities were divided by it.", n, ch)
        png = root / "flatfield" / f"{stem}.png"
        plot_flatfield(ff, f"{n} | {ch}: illumination (flat field, {ff.source} from {ff.n_images} images)", png, before, after)
        man.add(png, "flatfield", f"Flat field {n} / {ch}",
                f"Uneven illumination of {ch} in {n}: {ff_mod.describe(mag)}. Background-subtracted intensities are "
                "divided by this map; counts are not.", n, ch)
    if ff_rows:
        p = root / "flatfield.csv"
        pd.DataFrame(ff_rows).to_csv(p, index=False, sep=delim)
        man.add(p, "flatfield_table", "Flat-field magnitude and effect",
                "Per n and channel: how much dimmer the edges are than the centre (from the map and, as a check, from "
                "the nuclei before and after correction) and % positive with and without the correction.")

    # cut-offs and sweep
    if len(cutoffs):
        p = root / "cutoffs.csv"
        cutoffs.to_csv(p, index=False, sep=delim)
        man.add(p, "cutoffs", "Positivity cut-offs",
                "Per n, channel and measure (mean / 90th percentile): the rule used (negative control, locked manual "
                "value, limit of detection), its value in corrected intensity units, % positive, and where it came from.")
    if len(sweep):
        p = root / "positivity_sweep.csv"
        sweep.to_csv(p, index=False, sep=delim)
        man.add(p, "sweep_table", "Cut-off sweep",
                "% positive per n, condition and channel when the cut-off is scaled from 1/4 to 4 times the one used.")
        for (n, ch), g in sweep.groupby(["n", "channel"]):
            png = root / f"sweep_{_slug(n)}__{_slug(ch)}.png"
            plot_sweep(g, order, colours, f"{n} | {ch}: how % positive depends on the cut-off", png)
            man.add(png, "sweep", f"Cut-off sweep {n} / {ch}",
                    "% positive per condition against the cut-off (log scale); the grey line is the cut-off used. "
                    "Conditions whose curves cross change order depending on the cut-off: interpret with care.", n, ch)

    # distributions, per-condition summary, density
    rows = pd.DataFrame([row for _, row, _, _ in finfos])
    summ = []
    for (n, ch), g in nuc.groupby(["n", "channel"]):
        cr = cutoffs[(cutoffs["n"] == n) & (cutoffs["channel"] == ch) & (cutoffs["measure"] == "mean")] if len(cutoffs) else cutoffs
        cut, rule = (float(cr["value"].iloc[0]), cr["rule"].iloc[0]) if len(cr) else (np.nan, "")
        png = root / f"distribution_{_slug(n)}__{_slug(ch)}.png"
        plot_distribution(g, order, colours, cut, rule, f"{n} | {ch}: per-nucleus intensity "
                          f"({int(g['measurable'].sum())} measurable nuclei; {int((g['counted'] & g['saturated']).sum())} "
                          "saturated, excluded)", png)
        man.add(png, "distribution", f"Per-nucleus {ch} intensity {n}",
                f"Background-subtracted, flat-field-corrected mean {ch} per nucleus, per condition (violin + quartiles) "
                "and as histograms on a log scale; dashed = cut-off.", n, ch)
        m = g[g["measurable"]]
        for t, gt in m.groupby("treatment"):
            per_img = gt.groupby("image")["positive"].apply(lambda x: x.astype(float).mean()) if gt["positive"].notna().any() else pd.Series(dtype=float)
            summ.append(dict(n=n, channel=ch, treatment=t, images=gt["image"].nunique(), nuclei=len(gt),
                             saturated_excluded=int((g[(g["treatment"] == t)]["counted"] & g[(g["treatment"] == t)]["saturated"]).sum()),
                             frac_positive_pooled=round(float(gt["positive"].astype(float).mean()), 4) if gt["positive"].notna().any() else np.nan,
                             frac_positive_mean_of_images=round(float(per_img.mean()), 4) if len(per_img) else np.nan,
                             frac_positive_p90_pooled=round(float(gt["positive_p90"].astype(float).mean()), 4) if gt["positive_p90"].notna().any() else np.nan,
                             median_signal=round(float(gt["mean_bgsub_ff"].median()), 2),
                             median_signal_positive=round(float(gt.loc[gt["positive"].fillna(False).astype(bool), "mean_bgsub_ff"].median()), 2)
                             if gt["positive"].fillna(False).any() else np.nan,
                             cutoff_rule=rule, cutoff=cut))
        sl = _slug(ch)
        ycol = f"readout_frac_positive_{sl}"
        if ycol in rows and rows[ycol].notna().any():
            pi = rows[rows["n"] == n]
            xcol, xlab = ("readout_coverage", "nuclear coverage (fraction of image area covered by nuclei)") \
                if "readout_coverage" in pi and pi["readout_coverage"].notna().any() else ("auto_count", "nuclei per image")
            png = root / f"density_{_slug(n)}__{sl}.png"
            plot_density(pi, order, colours, xcol, xlab, ycol, f"{n} | {ch}: % positive vs local density (one dot per image)", png)
            man.add(png, "density", f"% positive vs density {n} / {ch}",
                    "Each image's % positive against its crowding. Dense fields often proliferate less (contact "
                    "inhibition); if conditions differ in density, part of a % positive difference may be density. "
                    "Descriptive; density-adjusted estimates are in the statistics.", n, ch)
    if summ:
        p = root / "positive_by_condition.csv"
        pd.DataFrame(summ).to_csv(p, index=False, sep=delim)
        man.add(p, "summary", "% positive per condition",
                "Per n, channel and condition: images, measurable nuclei, saturated nuclei left out, % positive (pooled "
                "over nuclei and mean of images), the 90th-percentile variant, median signal, and the cut-off rule.")

    # alignment and channel checks
    al = []
    for r, row, fi, _ in finfos:
        for ch, a in (fi.get("alignment") or {}).items():
            b = fi.get("brightness", {})
            al.append(dict(n=r.n, treatment=r.treatment, image=r.image, channel=ch, files="; ".join(Path(p).name for p in r.sources),
                           dy=a["dy"], dx=a["dx"], shift_px=a["shift_px"], corr=a["corr"], corrected=a.get("corrected", False),
                           nuclear_brightest=b.get("nuclear_brightest"), flags=row.get("channel_check", "")))
    if al:
        ad = pd.DataFrame(al)
        p = root / "alignment.csv"
        ad.to_csv(p, index=False, sep=delim)
        man.add(p, "alignment_table", "Channel alignment per field",
                "Shift between each marker file and its nuclear file (phase correlation), the correlation after "
                "shifting (low = possibly mis-paired files), whether the nuclear channel is the brightest, and flags.")
        warn = float(s["channels"].get("align_warn_px", 2.0))
        for n, g in ad.groupby("n"):
            png = root / f"alignment_{_slug(n)}.png"
            plot_alignment(g, warn, f"{n}: shift of marker vs nuclear image\n({len(g)} fields; "
                           f"{int((g['shift_px'] > warn).sum())} above {warn:g} px)", png)
            man.add(png, "alignment", f"Channel alignment {n}", "One dot per field: shift between the marker and the "
                    "nuclear image. Dots outside the dashed circle are flagged.", n)
    p = prj.work / "nuclei_intensity.csv"
    extra = dict(channels=[names[c] for c in chans], nuclear_channel=names[seg_ch],
                 per_nucleus_table=str(p.relative_to(prj.path).as_posix()),
                 cutoffs=_clean(cutoffs.to_dict("records")) if len(cutoffs) else [],
                 flatfield=_clean(ff_rows))
    man.write(extra)
    log(f"  intensity outputs: {root}")
    for c in (cutoffs[cutoffs["measure"] == "mean"].to_dict("records") if len(cutoffs) else []):
        log(f"  {c['n']} / {c['channel']}: cut-off {c['rule_label']} = {c['value']:.1f}, "
            f"{100 * c['frac_positive']:.1f}% positive ({100 * c['frac_positive_without_flatfield']:.1f}% without flat field)")
    for f in ff_rows:
        log(f"  {f['n']} / {f['channel']}: flat field {f['summary']}")
    return root
