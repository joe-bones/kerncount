"""Group plots laid out from the factors' display roles (row, panel, x, colour), never from factor names.

Absolute plot: one row of panels per level of the `row` factor that is present, one panel per level of the `panel`
factor, the `x` factor along the x-axis, every tumour a dot, log y, group n under each group. Statistics (ticket 04)
add normalised plots, forest and interaction plots; QC (ticket 06) adds hollow markers for flagged tumours.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import Factor

INK, MUTED = "#1f2328", "#8c959f"
DOT = "#3b5b92"


def _by_role(facs: list[Factor], role: str) -> Factor | None:
    return next((f for f in facs if f.role == role), None)


def _present(f: Factor | None, table: pd.DataFrame) -> list[str | None]:
    if f is None:
        return [None]
    have = set(table[f.name].astype(str))
    return [lv for lv in f.levels if lv in have] or [None]


def group_plot(table: pd.DataFrame, facs: list[Factor], value: str, label: str, out: Path, log: bool = True) -> Path:
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
                v = vals[s].dropna().to_numpy()
                v = v[v > 0] if log else v
                if len(v):
                    ax.scatter(k + rng.uniform(-0.15, 0.15, len(v)), v, s=22, color=DOT, edgecolor="white",
                               linewidth=0.6, zorder=3)
                    ax.hlines(np.exp(np.mean(np.log(v))) if log else np.mean(v), k - 0.28, k + 0.28, color=INK,
                              linewidth=1.5, zorder=2)
                ax.annotate(f"n={len(v)}", (k, 0), xycoords=("data", "axes fraction"), xytext=(0, -26),
                            textcoords="offset points", ha="center", va="top", fontsize=7, color=MUTED)
            ax.set_xticks(range(len(xs)), [str(x) if x is not None else "" for x in xs])
            ax.set_xlim(-0.6, len(xs) - 0.4)
            if log:
                ax.set_yscale("log")
            ax.grid(axis="y", color="#d0d7de", linewidth=0.5, which="major")
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            title = " · ".join(t for t in (r if len(rows) > 1 or row_f else None, p) if t)
            ax.set_title(title, fontsize=9, color=INK)
            if j == 0:
                ax.set_ylabel(label)
    fig.suptitle(f"{label}: every tumour one dot, line = geometric mean (no statistics yet)", fontsize=9,
                 color=MUTED)
    fig.tight_layout()
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)
    return out.with_suffix(".png")
