"""Group plots laid out from the factors' display roles (row, panel, x, colour), never from factor names.

Group plots: one row of panels per level of the `row` factor that is present, one panel per level of the `panel`
factor, the `x` factor along the x-axis, every tumour a dot, log y, group n under each group; optionally the model
estimate with its 95 % CI and a mark (stars / label) per group. Used for absolute and normalised values.
Forest plot: every contrast of a readout grouped by family. Interaction plot: model means per cell.
Flagged tumours (QC, ticket 06) are drawn as hollow markers.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import Factor

INK, MUTED = "#1f2328", "#8c959f"
DOT = "#3b5b92"
EST = "#c2410c"
SERIES = ["#3b5b92", "#c2410c", "#2f855a", "#7c3aed", "#b7791f", "#0e7490"]

GroupKey = tuple  # (row level, panel level, x level); None where the role is not declared


def _by_role(facs: list[Factor], role: str) -> Factor | None:
    return next((f for f in facs if f.role == role), None)


def _present(f: Factor | None, table: pd.DataFrame) -> list[str | None]:
    if f is None:
        return [None]
    have = set(table[f.name].astype(str))
    return [lv for lv in f.levels if lv in have] or [None]


def _plain_log(ax, axis: str) -> None:
    """Readable log-axis ticks (1, 2, 5 per decade, plain numbers) instead of 4×10⁻¹ labels."""
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    a = ax.yaxis if axis == "y" else ax.xaxis
    a.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    a.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    a.set_minor_formatter(NullFormatter())


def _save(fig, out: Path) -> Path:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    fig.savefig(out.with_suffix(".pdf"))
    import matplotlib.pyplot as plt

    plt.close(fig)
    return out.with_suffix(".png")


def group_plot(table: pd.DataFrame, facs: list[Factor], value: str, label: str, out: Path, log: bool = True,
               ci: dict[GroupKey, tuple[float, float, float]] | None = None, marks: dict[GroupKey, str] | None = None,
               ref_line: float | None = None, suptitle: str | None = None, hollow: str | None = None) -> Path:
    """`ci` maps (row, panel, x) -> (estimate, low, high) drawn beside the dots instead of the plain geometric-mean
    line; `marks` maps the same key to a text above the group (stars, "confounded", …); `ref_line` draws e.g. 100 %.
    `hollow` names a boolean column (QC `flagged`): those tumours are drawn as hollow markers (ticket 06)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    row_f, panel_f, x_f = _by_role(facs, "row"), _by_role(facs, "panel"), _by_role(facs, "x")
    rows, panels = _present(row_f, table), _present(panel_f, table)
    xs = list(x_f.levels) if x_f else [None]
    fig, axes = plt.subplots(len(rows), len(panels), figsize=(2.2 + 1.0 * len(xs) * len(panels), 0.6 + 3.0 * len(rows)),
                             sharey=True, squeeze=False)
    rng = np.random.default_rng(0)
    vals = pd.to_numeric(table[value], errors="coerce")
    is_hollow = (_as_bool(table[hollow]) if hollow is not None and hollow in table.columns
                 else pd.Series(False, index=table.index))
    for i, r in enumerate(rows):
        for j, p in enumerate(panels):
            ax = axes[i][j]
            sel = pd.Series(True, index=table.index)
            if row_f is not None and r is not None:
                sel &= table[row_f.name].astype(str) == r
            if panel_f is not None and p is not None:
                sel &= table[panel_f.name].astype(str) == p
            for k, xv in enumerate(xs):
                s = sel & (table[x_f.name].astype(str) == xv) if x_f is not None and xv is not None else sel
                ok = s & vals.notna() & ((vals > 0) if log else True)
                v = vals[ok].to_numpy()
                hv = is_hollow[ok].to_numpy()
                key = (r, p, xv)
                if len(v):
                    jit = k + rng.uniform(-0.15, 0.15, len(v))
                    if (~hv).any():
                        ax.scatter(jit[~hv], v[~hv], s=22, color=DOT, edgecolor="white", linewidth=0.6, zorder=3)
                    if hv.any():
                        ax.scatter(jit[hv], v[hv], s=22, facecolors="none", edgecolors=DOT, linewidth=1.0, zorder=3)
                    if ci is not None and key in ci and np.isfinite(ci[key][0]):
                        e, lo, hi = ci[key]
                        ax.hlines(e, k - 0.28, k + 0.28, color=INK, linewidth=1.5, zorder=2)
                        if np.isfinite(lo) and np.isfinite(hi):
                            ax.vlines(k + 0.3, lo, hi, color=EST, linewidth=1.4, zorder=4)
                            ax.plot([k + 0.3], [e], marker="o", ms=3.5, color=EST, zorder=5)
                    else:
                        ax.hlines(np.exp(np.mean(np.log(v))) if log else np.mean(v), k - 0.28, k + 0.28, color=INK,
                                  linewidth=1.5, zorder=2)
                if marks and marks.get(key):
                    ax.annotate(marks[key], (k, 1.0), xycoords=("data", "axes fraction"), xytext=(0, -2),
                                textcoords="offset points", ha="center", va="top", fontsize=7, color=INK)
                ax.annotate(f"n={len(v)}", (k, 0), xycoords=("data", "axes fraction"), xytext=(0, -26),
                            textcoords="offset points", ha="center", va="top", fontsize=7, color=MUTED)
            ax.set_xticks(range(len(xs)), [str(x) if x is not None else "" for x in xs])
            ax.set_xlim(-0.6, len(xs) - 0.4)
            if log:
                ax.set_yscale("log")
                _plain_log(ax, "y")
            if ref_line is not None:
                ax.axhline(ref_line, color=MUTED, linewidth=0.8, linestyle="--", zorder=1)
            ax.grid(axis="y", color="#d0d7de", linewidth=0.5, which="major")
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            title = " · ".join(t for t in (r if len(rows) > 1 or row_f else None, p) if t)
            ax.set_title(title, fontsize=9, color=INK)
            if j == 0:
                ax.set_ylabel(label)
    title = suptitle or f"{label}: every tumour one dot, line = geometric mean"
    if bool(is_hollow.any()):
        title += f"; hollow = flagged by QC ({int(is_hollow.sum())}: truncated, kept in the analysis)"
    fig.suptitle(title, fontsize=9, color=MUTED)
    fig.tight_layout()
    return _save(fig, out)


def _as_bool(col: pd.Series) -> pd.Series:
    """Boolean QC column read back from CSV / Excel (True/False, yes/no, 1/0)."""
    if col.dtype == bool:
        return col
    return col.astype(str).str.strip().str.lower().isin(("true", "yes", "1", "1.0"))


def forest_plot(results: pd.DataFrame, label: str, out: Path) -> Path:
    """One row per contrast, grouped by family (header row with role), fold-change and 95 % CI on a log axis.
    Filled = tested; hollow grey = descriptive / confounded / not tested; rows that need data show only the text."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lines: list[tuple[str, pd.Series | None]] = []
    for fid, g in results.groupby("family", sort=False):
        r0 = g.iloc[0]
        lines.append((f"{fid} — {r0['family_label'] or ''} [{r0['family_role']}]", None))
        lines += [(str(r["comparison"]), r) for _, r in g.iterrows()]
    fig, ax = plt.subplots(figsize=(13, 0.6 + 0.24 * max(len(lines), 1)))
    finite = results[["ci_low", "ci_high", "estimate"]].to_numpy(float)
    finite = finite[np.isfinite(finite) & (finite > 0)]
    lo, hi = (float(finite.min()), float(finite.max())) if finite.size else (0.5, 2.0)
    lo, hi = min(lo, 0.9) / 1.15, max(hi, 1.1) * 1.15
    for i, (text, r) in enumerate(lines):
        y = len(lines) - 1 - i
        if r is None:
            ax.text(lo, y, text, fontsize=7.5, fontweight="bold", color=INK, va="center", ha="left")
            continue
        ax.text(-0.01, y, text, fontsize=7, color=INK, va="center", ha="right", transform=ax.get_yaxis_transform())
        tested = np.isfinite(r["p"])
        e, cl, ch = float(r["estimate"]), float(r["ci_low"]), float(r["ci_high"])
        if np.isfinite(cl) and np.isfinite(ch):
            ax.hlines(y, cl, ch, color=INK if tested else MUTED, linewidth=1.2)
        if np.isfinite(e):
            ax.plot([e], [y], marker="o", ms=4.5, color=INK if tested else "white",
                    markeredgecolor=INK if tested else MUTED)
        if np.isfinite(e):
            txt = f"{e:.2f} [{cl:.2f}, {ch:.2f}]" if np.isfinite(cl) else f"{e:.2f}"
        else:
            txt = ""
        sig = str(r["significance"] or "")
        note = str(r["note"] or "")
        txt = "  ".join(x for x in (txt, sig, note) if x)
        ax.text(1.01, y, txt, fontsize=6.5, color=INK if tested else MUTED, va="center", ha="left",
                transform=ax.get_yaxis_transform())
    ax.axvline(1.0, color=MUTED, linewidth=0.8, linestyle="--")
    ax.set_xscale("log")
    _plain_log(ax, "x")
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.7, len(lines) - 0.3)
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel(f"{label}: fold-change (95 % CI), log scale")
    ax.set_title(f"{label}: contrasts by family (Holm within family; stars only for confirmatory families of "
                 "confirmatory readouts)", fontsize=8.5, color=MUTED)
    fig.subplots_adjust(left=0.28, right=0.55, top=1 - 0.4 / fig.get_figheight(), bottom=0.5 / fig.get_figheight())
    return _save(fig, out)


def interaction_plot(means: pd.DataFrame, facs: list[Factor], label: str, out: Path) -> Path:
    """Model means (95 % CI) per cell: x factor on the x-axis, one line per panel-factor level, one plot per row-factor
    level. Parallel lines on the log axis = no interaction."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    row_f, panel_f, x_f = _by_role(facs, "row"), _by_role(facs, "panel"), _by_role(facs, "x")
    rows, panels = _present(row_f, means), _present(panel_f, means)
    xs = list(x_f.levels) if x_f else [None]
    fig, axes = plt.subplots(1, len(rows), figsize=(0.8 + 3.2 * len(rows), 3.4), sharey=True, squeeze=False)
    for i, r in enumerate(rows):
        ax = axes[0][i]
        for j, p in enumerate(panels):
            sel = pd.Series(True, index=means.index)
            if row_f is not None and r is not None:
                sel &= means[row_f.name].astype(str) == r
            if panel_f is not None and p is not None:
                sel &= means[panel_f.name].astype(str) == p
            pts = []
            for k, xv in enumerate(xs):
                s = sel & (means[x_f.name].astype(str) == xv) if x_f is not None and xv is not None else sel
                if s.sum() == 1:
                    m = means[s].iloc[0]
                    pts.append((k + (j - (len(panels) - 1) / 2) * 0.08, m["mean"], m["ci_low"], m["ci_high"]))
            if not pts:
                continue
            col = SERIES[j % len(SERIES)]
            xk, e, lo, hi = (np.array(v, float) for v in zip(*pts))
            ax.plot(xk, e, marker="o", ms=4, color=col, linewidth=1.3, label=str(p) if p is not None else None)
            ok = np.isfinite(lo) & np.isfinite(hi)
            ax.vlines(xk[ok], lo[ok], hi[ok], color=col, linewidth=1.0)
        ax.set_xticks(range(len(xs)), [str(x) if x is not None else "" for x in xs])
        ax.set_xlim(-0.5, len(xs) - 0.5)
        ax.set_yscale("log")
        _plain_log(ax, "y")
        ax.grid(axis="y", color="#d0d7de", linewidth=0.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_title(str(r) if r is not None else "", fontsize=9, color=INK)
        if i == 0:
            ax.set_ylabel(f"{label}\nmodel mean, 95 % CI")
        if any(p is not None for p in panels):
            ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f"{label}: interaction plot (parallel lines = no interaction on the fold-change scale)", fontsize=9,
                 color=MUTED)
    fig.tight_layout()
    return _save(fig, out)


# Exploratory shape readouts in the shape panel (ticket 05): column, label; linear axes (fractions and counts have 0s).
SHAPE_PANEL = [("elongation", "Elongation (major / minor axis)"), ("solidity", "Solidity (area / convex hull)"),
               ("edge_roughness", "Edge roughness (P / P smoothed)"), ("circularity", "Circularity 4πA/P² (reference)"),
               ("invasion_fraction", "Invasion fraction"), ("n_satellites", "Satellites (count)")]
MARKERS = ["o", "^", "s", "D"]


def shape_panel(table: pd.DataFrame, facs: list[Factor], out: Path) -> Path:
    """Shape panel (exploratory): cross-sectional area vs perimeter of every tumour with the circle reference
    P = 2·sqrt(π·A) (points above it: elongated / lobed / ragged outlines), and the exploratory shape readouts by group
    (every tumour a dot, line = group median). Colour = the `x` factor, marker = the `panel` factor."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    from .families import cells

    x_f, panel_f = _by_role(facs, "x"), _by_role(facs, "panel")
    shown = [(c, lab) for c, lab in SHAPE_PANEL if c in table.columns]
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.2), squeeze=False)
    flat = [ax for row in axes for ax in row]

    def colour(row: pd.Series) -> str:
        return SERIES[x_f.levels.index(str(row[x_f.name])) % len(SERIES)] if x_f is not None and str(
            row[x_f.name]) in x_f.levels else DOT

    def marker(row: pd.Series) -> str:
        return MARKERS[panel_f.levels.index(str(row[panel_f.name])) % len(MARKERS)] if panel_f is not None and str(
            row[panel_f.name]) in panel_f.levels else "o"

    # area vs perimeter (log-log) with the circle reference
    ax = flat[0]
    def column(name: str) -> pd.Series:
        return (pd.to_numeric(table[name], errors="coerce") if name in table.columns
                else pd.Series(np.nan, index=table.index))

    a, p = column("cross_sectional_area_mm2"), column("perimeter_mm")
    for i, row in table.iterrows():
        if np.isfinite(a[i]) and np.isfinite(p[i]) and a[i] > 0 and p[i] > 0:
            ax.scatter(a[i], p[i], s=30, marker=marker(row), color=colour(row), edgecolor="white", linewidth=0.6,
                       zorder=3)
    ok = a[(a > 0) & np.isfinite(a)]
    if len(ok):
        g = np.geomspace(ok.min() / 1.3, ok.max() * 1.3, 50)
        ax.plot(g, 2 * np.sqrt(np.pi * g), color=MUTED, linestyle="--", linewidth=1.0, label="circle: P = 2√(πA)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        _plain_log(ax, "x")
        _plain_log(ax, "y")
    ax.set_xlabel("Cross-sectional area (mm²)")
    ax.set_ylabel("Perimeter (mm)")
    ax.set_title("Mid-section: area vs perimeter", fontsize=9, color=INK)
    handles = [Line2D([], [], color=MUTED, linestyle="--", label="circle: P = 2√(πA)")]
    if x_f is not None:
        handles += [Line2D([], [], marker="o", linestyle="", color=SERIES[k % len(SERIES)], label=lv)
                    for k, lv in enumerate(x_f.levels) if lv in set(table[x_f.name].astype(str))]
    if panel_f is not None:
        handles += [Line2D([], [], marker=MARKERS[k % len(MARKERS)], linestyle="", color=INK, label=lv)
                    for k, lv in enumerate(panel_f.levels) if lv in set(table[panel_f.name].astype(str))]
    ax.legend(handles=handles, fontsize=7, frameon=False)

    # exploratory shape readouts by group
    from .stats import cell_of

    keys = cell_of(table, facs) if facs else pd.Series([()] * len(table), index=table.index)
    groups = [c for c in cells(facs) if (keys == c).any()] if facs else [()]
    rng = np.random.default_rng(0)
    for ax, (col, lab) in zip(flat[1:], shown):
        vals = pd.to_numeric(table[col], errors="coerce")
        for k, cell in enumerate(groups):
            sel = (keys == cell) & vals.notna()
            if not sel.any():
                continue
            rows = table[sel]
            jit = k + rng.uniform(-0.15, 0.15, len(rows))
            for (_, row), xj, v in zip(rows.iterrows(), jit, vals[sel]):
                ax.scatter(xj, v, s=18, marker=marker(row), color=colour(row), edgecolor="white", linewidth=0.5,
                           zorder=3)
            ax.hlines(float(vals[sel].median()), k - 0.28, k + 0.28, color=INK, linewidth=1.3, zorder=2)
        ax.set_xticks(range(len(groups)), [" / ".join(c) for c in groups], rotation=60, ha="right", fontsize=6)
        ax.set_xlim(-0.6, len(groups) - 0.4)
        ax.set_title(lab, fontsize=9, color=INK)
        ax.grid(axis="y", color="#d0d7de", linewidth=0.5)
        if col == "n_satellites":
            from matplotlib.ticker import MaxNLocator

            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    for ax in flat:
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in flat[1 + len(shown):]:
        ax.set_visible(False)
    fig.suptitle("Shape panel (exploratory, never starred): 2D shape on the mid-section; every tumour one dot, "
                 "line = group median", fontsize=9, color=MUTED)
    fig.tight_layout()
    return _save(fig, out)
