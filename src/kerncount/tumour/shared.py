"""Shared-threshold robustness check (ADR 0002, ticket 11): every tumour segmented at ONE experiment-wide threshold.

The primary method thresholds each tumour at its own automatic level (half-way between its background and tumour
signal), which removes brightness differences between tumours. ADR 0002 keeps one experiment-wide threshold as a
robustness check when the acquisition settings are identical across tumours (`acquisition.consistency`): grey values
are then directly comparable, and a confirmatory effect that also holds at one shared threshold does not rest on the
per-tumour adaptation. An effect that appears only at the shared threshold is a brightness difference, not a size
difference (a dimmer tumour shrinks at a fixed threshold).

Rule (fixed, not tuned): shared threshold = median of the per-tumour automatic thresholds of the analysed tumours.
Each stack is re-segmented once at that grey value (`segment.sweep` with scale = shared / own automatic threshold, so
smoothing, fill, components and readouts are those of the primary method); readouts cached in
`_work/tumour/<slug>.shared.json`. Not run when a checked acquisition setting differs (reported instead).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .config import Factor, Family, Readout, readouts, tumour_settings
from .measure import WORK, _key, raw_dir, slug

FILENAME = "shared_threshold.csv"
COLUMNS = ["readout", "family", "family_role", "comparison", "estimate", "ci_low", "ci_high", "p_adj",
           "shared_estimate", "shared_ci_low", "shared_ci_high", "shared_p_adj", "direction_agree",
           "significance_agree", "note"]


def value(table: pd.DataFrame) -> float:
    """The experiment-wide threshold: median of the per-tumour automatic thresholds (`threshold` column)."""
    t = pd.to_numeric(table["threshold"], errors="coerce")
    t = t[np.isfinite(t) & (t > 0)]
    if not len(t):
        raise ValueError("no per-tumour automatic threshold in the table")
    return float(np.median(t))


def measurements(project_dir: Path, table: pd.DataFrame, settings: dict, threshold: float) -> pd.DataFrame:
    """tumour_id + segmentation readouts at the grey value `threshold` for every tumour of `table` (cached)."""
    from .segment import sweep_file

    prj = Path(project_dir)
    work = prj / WORK
    work.mkdir(parents=True, exist_ok=True)
    raw = raw_dir(prj, settings)
    ts = tumour_settings(settings)
    rows = []
    for i, (_, row) in enumerate(table.iterrows(), 1):
        tid, src = row["tumour_id"], raw / row["file"]
        key = f"{_key(src, ts)}:shared={threshold:.4f}"
        path = work / f"{slug(tid)}.shared.json"
        cached = json.loads(path.read_text()) if path.exists() else None
        if cached and cached.get("key") == key:
            res = cached["readouts"]
        else:
            t0 = time.time()
            own = float(row["threshold"])
            res = sweep_file(src, ts, [threshold / own])[threshold / own]
            path.write_text(json.dumps({"tumour_id": tid, "key": key, "readouts": res}, indent=1, default=float))
            print(f"  shared threshold [{i}/{len(table)}] {tid}: own {own:.1f} -> {threshold:.1f} "
                  f"({time.time() - t0:.0f} s)", flush=True)
        rows.append({"tumour_id": tid, **res})
    return pd.DataFrame(rows)


def compare(new: pd.DataFrame, shared: pd.DataFrame) -> pd.DataFrame:
    """Per contrast: primary vs shared-threshold estimate, Holm p, direction / significance agreement."""
    from .concordance import compare as _compare

    c = _compare(new, shared).rename(columns=lambda k: k.replace("fiji_", "shared_"))
    c["note"] = c["note"].str.replace("Fiji side", "shared-threshold side")
    return c[COLUMNS]


def analyse(out: Path, table: pd.DataFrame, settings: dict, facs: list[Factor], fams: list[Family],
            report: Path | None, consistency=None) -> list[str]:
    """Run the check when `[tumour] shared_threshold_check` is on and the acquisition settings allow it; writes
    shared_threshold.csv, plots/<readout>_vs_shared_threshold, an xlsx sheet and a report section."""
    from . import concordance, stats, sweep

    ts = tumour_settings(settings)
    conf = [r for r in readouts(settings) if r.role == "confirmatory" and r.column in table.columns]
    md = ["## Shared-threshold robustness check (ADR 0002)", ""]
    written: list[str] = []
    cmp = pd.DataFrame(columns=COLUMNS)
    # outputs of an earlier run (e.g. before a batch with other settings arrived) must not outlive this one
    for p in [out / FILENAME, *(out / "plots").glob("*_vs_shared_threshold.*")]:
        p.unlink(missing_ok=True)
    if not ts["shared_threshold_check"]:
        md.append("- not run (`[tumour] shared_threshold_check = false`).")
    elif consistency is not None and not consistency.identical:
        md.append(f"- not run: {consistency.describe()}; one grey value is not comparable across tumours.")
    elif not conf or "file" not in table.columns or "threshold" not in table.columns:
        md.append("- not run (no confirmatory readout / no file or threshold column).")
    else:
        t = value(table)
        shared_tab = sweep.rebuild(table, measurements(out, table, settings, t), settings)
        cmps = []
        for ro in conf:
            if ro.column not in shared_tab.columns:
                md.append(f"- {ro.name}: not produced by the segmentation step; not checked.")
                continue
            new = stats.readout_results(table, facs, fams, ro)
            alt = stats.readout_results(shared_tab, facs, fams,
                                        Readout(f"{ro.name}_shared", ro.column, ro.label, ro.role))
            cmps.append(compare(new, alt))
            both = table[["tumour_id", *[f.name for f in facs], ro.column]].merge(
                shared_tab[["tumour_id", ro.column]].rename(columns={ro.column: "shared"}), on="tumour_id")
            written.append(str(concordance.scatter_plot(
                both, facs, ro.column, "shared", ro.label, out / "plots" / f"{ro.name}_vs_shared_threshold",
                x_label=f"shared threshold {t:.1f}: {ro.label}", x_name="shared")))
        cmp = pd.concat(cmps, ignore_index=True) if cmps else cmp
        delim = settings.get("project", {}).get("csv_delimiter", ",")
        cmp.to_csv(out / FILENAME, index=False, sep=delim)
        written.append(str(out / FILENAME))
        own = pd.to_numeric(table["threshold"], errors="coerce")
        md += [f"- Shared threshold = median of the per-tumour automatic thresholds = **{t:.1f}** (grey value of the "
               f"smoothed stack; per-tumour range {own.min():.1f}–{own.max():.1f}). "
               + (consistency.describe().capitalize() + "." if consistency is not None else ""),
               "- Same statistics (families, Holm) on the shared-threshold readouts; a robustness check, never "
               "starred. An effect present only at the shared threshold reflects brightness, not size.", "",
               "| readout | family | comparison | primary [95 % CI] | shared threshold [95 % CI] | direction | "
               "significance agrees |", "|---|---|---|---|---|---|---|"]
        for _, r in cmp[cmp["family_role"] == "confirmatory"].iterrows():
            md.append(f"| {r['readout']} | {r['family']} | {r['comparison'].replace('|', chr(92) + '|')} | "
                      f"{_ci(r['estimate'], r['ci_low'], r['ci_high'])} | "
                      f"{_ci(r['shared_estimate'], r['shared_ci_low'], r['shared_ci_high'])} | "
                      f"{r['direction_agree']} | {r['significance_agree'] or '—'} |")
        differs = cmp[cmp["significance_agree"] == "no"]
        if len(differs):
            md += ["", "**The conclusion differs at the shared threshold for:** " + "; ".join(
                f"{r['readout']} {r['family']} {r['comparison']}" for _, r in differs.iterrows())]
    xlsx = out / "results.xlsx"
    if xlsx.exists():
        with pd.ExcelWriter(xlsx, engine="openpyxl", mode="a", if_sheet_exists="replace") as xw:
            cmp.to_excel(xw, sheet_name="Shared threshold", index=False)
    if report is not None and report.exists():
        with report.open("a", encoding="utf-8") as fh:
            fh.write("\n" + "\n".join(md) + "\n\n")
    return written


def _ci(e: float, lo: float, hi: float) -> str:
    if not np.isfinite(e):
        return "—"
    return f"{e:.2f} [{lo:.2f}, {hi:.2f}]" if np.isfinite(lo) else f"{e:.2f}"
