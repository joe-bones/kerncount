"""Quality control of tumours (ticket 06, ADR 0002 + amendment 2026-11-20): pre-set rules, never judgement calls.

Per tumour, from its (possibly corrected) labels — label 1 = main mass — and the original stack:
- **xy-edge contact**: the main mass touches the image frame in any layer; `edge_outline_frac` = share of the main
  mass's outline (footprint, all layers merged) lying on the frame, `edge_outline_frac_mid` = the same per layer of
  the mid-section, averaged. Outline lengths are measured along sub-pixel contours (marching squares).
- **mid-section fit**: the fixed-thickness mid-section (centred on the layer of largest main-mass cross-section)
  must lie inside the imaged layers. z-extent beyond the stack alone is NOT a flag: every tumour extends beyond it.
- **truncated** = edge contact or a mid-section that does not fit. Truncated tumours are `flagged`: kept in the
  primary analysis, drawn hollow, dropped in the sensitivity analysis.
- **saturation**: share of main-mass voxels in the mid-section at the top grey value; a warning above
  `[tumour] qc_saturation_warn` (default 40 %: the intermedians threshold is protected up to ~50 %). Not a flag.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import tumour_settings

NAN = float("nan")
QC_COLUMNS = ["edge_contact", "edge_outline_frac", "edge_outline_frac_mid", "mid_section_fits",
              "saturated_frac_mid_section", "saturation_warning", "truncated", "flagged", "qc_flags"]
YES, NO = "yes", "no"


def _shape_options(settings: dict):
    """Ticket 05's readout scales (mid-section thickness etc.) from the [tumour] settings."""
    from .shape import ShapeOptions

    return ShapeOptions.from_settings(tumour_settings(settings))


def mid_section(area_per_layer: np.ndarray, voxel_z_um: float, settings: dict) -> tuple[int, int, int, bool]:
    """(centre, first, last, fits), 0-based and clipped to the stack: ticket 05's mid-section rule
    (`shape.mid_section`, [tumour] mid_section_um), so QC and the readouts always use the same layers."""
    from .shape import mid_section as rule

    ms = rule(area_per_layer, voxel_z_um, _shape_options(settings).mid_section_um)
    return ms.centre, ms.first, ms.last, bool(ms.fits)


def edge_outline_frac(mask: np.ndarray) -> float:
    """Share of the outline length of a 2D mask that lies on the image frame (NaN for an empty mask)."""
    from skimage.measure import find_contours

    if not mask.any():
        return NAN
    p = np.pad(mask, 1).astype(np.float32)
    lo, hy, hx = 0.5, p.shape[0] - 1.5, p.shape[1] - 1.5
    total = on_edge = 0.0
    for c in find_contours(p, 0.5):
        y, x = c[:, 0], c[:, 1]
        seg = np.hypot(np.diff(y), np.diff(x))
        total += float(seg.sum())
        for side in (np.isclose(y, lo), np.isclose(y, hy), np.isclose(x, lo), np.isclose(x, hx)):
            on_edge += float(seg[side[:-1] & side[1:]].sum())
    return on_edge / total if total else NAN


def _yn(b: bool) -> str:
    return YES if b else NO


def tumour_qc(labels: np.ndarray, stack: np.ndarray, voxel_um: tuple[float, float, float], settings: dict) -> dict:
    """QC row of one tumour (QC_COLUMNS; flags as "yes"/"no" like the other tumours.csv flags) from its labels
    (1 = main mass) and the original ZYX stack."""
    return assess(labels, stack, voxel_um, settings)[0]


def assess(labels: np.ndarray, stack: np.ndarray, voxel_um: tuple[float, float, float], settings: dict
           ) -> tuple[dict, int]:
    """(QC row, 0-based mid-section centre layer; -1 without a main mass)."""
    ts = tumour_settings(settings)
    warn = float(ts.get("qc_saturation_warn", 0.40))
    from scipy import ndimage as ndi

    from .shape import section_area

    main = labels == 1
    if not main.any():
        return dict(edge_contact=NO, edge_outline_frac=NAN, edge_outline_frac_mid=NAN, mid_section_fits="",
                    saturated_frac_mid_section=NAN, saturation_warning=NO, truncated=NO, flagged=NO,
                    qc_flags="no main mass"), -1
    contact = bool(main[:, 0].any() or main[:, -1].any() or main[:, :, 0].any() or main[:, :, -1].any())
    frac = edge_outline_frac(main.any(axis=0)) if contact else 0.0
    sl = ndi.find_objects(main.any(axis=0).astype(np.uint8))[0]  # areas as ticket 05 measures them (holes filled)
    area = np.array([section_area(m[sl], (voxel_um[1], voxel_um[2])) if m.any() else 0.0 for m in main])
    centre, first, last, fits = mid_section(area, voxel_um[0], settings)
    mid = main[first:last + 1]
    if contact:
        fr = [edge_outline_frac(m) for m in mid if m.any()]
        frac_mid = float(np.mean(fr)) if fr else NAN
    else:
        frac_mid = 0.0
    top = np.iinfo(stack.dtype).max if np.issubdtype(stack.dtype, np.integer) else float(stack.max())
    n_mid = int(mid.sum())
    sat = int(np.count_nonzero(stack[first:last + 1][mid] >= top)) / n_mid if n_mid else NAN
    sat_warn = bool(np.isfinite(sat) and sat > warn)
    flags = []
    if contact:
        share = "<1%" if frac < 0.005 else f"{frac:.0%}"  # a few touching pixels must not read as "0%"
        flags.append(f"main mass touches the image edge ({share} of its outline)")
    if not fits:
        flags.append(f"mid-section does not fit in the {len(area)} imaged layers")
    if sat_warn:
        flags.append(f"{sat:.0%} of the mid-section main mass saturated (warning above {warn:.0%})")
    truncated = contact or not fits
    return dict(edge_contact=_yn(contact), edge_outline_frac=round(frac, 4),
                edge_outline_frac_mid=round(frac_mid, 4) if np.isfinite(frac_mid) else NAN,
                mid_section_fits=_yn(fits),
                saturated_frac_mid_section=round(sat, 4) if np.isfinite(sat) else NAN,
                saturation_warning=_yn(sat_warn), truncated=_yn(truncated), flagged=_yn(truncated),
                qc_flags="; ".join(flags)), centre


FLOW_COLUMNS = ["files", "tumours", "excluded", "exclusion_reasons", "analysed", "flagged", "analysed_without_flagged"]


def flagged_mask(table) -> pd.Series:
    """The QC `flagged` column as booleans (False where absent; tolerant of CSV / Excel round trips)."""

    if "flagged" not in table.columns:
        return pd.Series(False, index=table.index)
    col = table["flagged"]
    if col.dtype == bool:
        return col
    return col.astype(str).str.strip().str.lower().isin(("true", "yes", "1", "1.0"))


def flow_table(table, excluded, facs) -> pd.DataFrame:
    """Per group (every factor level combination present in the layout) + a total row 'all':
    files -> tumours (readable ID) -> excluded (reasons from excluded_tumours.csv, incl. review exclusions)
    -> analysed (primary) -> flagged (QC) -> analysed without flagged (sensitivity)."""

    from .families import cells

    names = [f.name for f in facs]
    flg = flagged_mask(table)

    def key(df: pd.DataFrame) -> pd.Series:
        if not len(df):
            return pd.Series([], dtype=object)
        return df.apply(lambda r: tuple(str(r.get(n, "") if pd.notna(r.get(n, "")) else "") for n in names), axis=1)

    tk, ek = list(key(table)), list(key(excluded))
    present = set(tk) | set(ek)
    groups = [c for c in cells(facs) if c in present]
    groups += sorted(present - set(groups))
    rows = []
    for c in [*groups, None]:
        t_sel = pd.Series([c is None or k == c for k in tk], index=table.index, dtype=bool)
        ex = excluded[[c is None or k == c for k in ek]] if len(excluded) else excluded
        reasons = (ex["exclusion_reason"].fillna("").astype(str).replace("", "no reason given").value_counts()
                   if len(ex) and "exclusion_reason" in ex.columns else pd.Series(dtype=int))
        ids = ex["tumour_id"].fillna("").astype(str) if len(ex) and "tumour_id" in ex.columns else pd.Series([], dtype=str)
        n_t, n_f = int(t_sel.sum()), int((t_sel & flg).sum())
        rows.append({**{n: (c[i] if c is not None else "all") for i, n in enumerate(names)},
                     "files": n_t + len(ex), "tumours": n_t + int((ids != "").sum()), "excluded": len(ex),
                     "exclusion_reasons": "; ".join(f"{r} ({k})" for r, k in reasons.items()),
                     "analysed": n_t, "flagged": n_f, "analysed_without_flagged": n_t - n_f})
    return pd.DataFrame(rows, columns=[*names, *FLOW_COLUMNS])


# ---------------------------------------------------------------- results: flow, sensitivity without flagged, A-B plot

SIDE_COLUMNS = ["readout", "readout_role", "family", "family_role", "comparison", "estimate", "ci_low", "ci_high",
                "p_adj", "significance", "N", "without_flagged_estimate", "without_flagged_ci_low",
                "without_flagged_ci_high", "without_flagged_p_adj", "without_flagged_N", "conclusion_changed"]
FLOW_FILE = "flow.csv"


def without_flagged(table: pd.DataFrame, facs, fams, readout) -> pd.DataFrame:
    """The primary statistics of one readout beside the same statistics on the table without flagged tumours."""
    from . import stats

    prim = stats.readout_results(table, facs, fams, readout)
    sub = table[~flagged_mask(table)]
    alt = stats.readout_results(sub, facs, fams, readout).set_index(["family", "comparison"])
    rows = []
    for _, r in prim.iterrows():
        a = alt.loc[(r["family"], r["comparison"])] if (r["family"], r["comparison"]) in alt.index else None
        ap = a["p_adj"] if a is not None else NAN
        changed = (bool((r["p_adj"] < 0.05) != (ap < 0.05)) if np.isfinite(r["p_adj"]) and np.isfinite(ap)
                   else None)
        rows.append({"readout": readout.name, "readout_role": readout.role, "family": r["family"],
                     "family_role": r["family_role"], "comparison": r["comparison"], "estimate": r["estimate"],
                     "ci_low": r["ci_low"], "ci_high": r["ci_high"], "p_adj": r["p_adj"],
                     "significance": r["significance"], "N": r["N"],
                     "without_flagged_estimate": a["estimate"] if a is not None else NAN,
                     "without_flagged_ci_low": a["ci_low"] if a is not None else NAN,
                     "without_flagged_ci_high": a["ci_high"] if a is not None else NAN,
                     "without_flagged_p_adj": ap, "without_flagged_N": a["N"] if a is not None else NAN,
                     "conclusion_changed": changed})
    return pd.DataFrame(rows, columns=SIDE_COLUMNS)


def _long_rows(side: pd.DataFrame, readout, lock_hash: str) -> list[dict]:
    """stats_long rows, family `<id> (without flagged)` (cf. `(Welch-type)`, `(leave-one-n-out)`); never starred."""
    rows = []
    for _, r in side[side["without_flagged_estimate"].notna()].iterrows():
        note = [f"sensitivity without QC-flagged tumours ({r['family_role']} family)"]
        if r["conclusion_changed"] is True:
            note.append("conclusion differs from the primary analysis")
        rows.append(dict(level="all_n", n=None, metric=readout.name, family=f"{r['family']} (without flagged)",
                         test="sensitivity: primary model without flagged (truncated) tumours, Holm within family",
                         comparison=r["comparison"], estimate=r["without_flagged_estimate"], estimate_type="ratio",
                         ci_low=r["without_flagged_ci_low"], ci_high=r["without_flagged_ci_high"], p=NAN,
                         p_adj=r["without_flagged_p_adj"], significance="", N=r["without_flagged_N"],
                         model="N-way linear model on log(readout) without flagged tumours",
                         note="; ".join([*note, f"lock {lock_hash}"]), estimand=f"fold-change of {readout.label}",
                         cell_line=""))
    return rows


def hemisphere_pairs(table: pd.DataFrame, facs, column: str) -> pd.DataFrame:
    """One row per section with both an A and a B tumour (same group, imaging session and section label)."""
    need = {"section", "hemisphere", column}
    if not need <= set(table.columns):
        return pd.DataFrame(columns=["group", "section", "A", "B", "flagged"])
    keys = [f.name for f in facs] + [c for c in ("imaging_session",) if c in table.columns] + ["section"]
    t = table.assign(_v=pd.to_numeric(table[column], errors="coerce"), _f=flagged_mask(table))
    rows = []
    for k, g in t.groupby(keys, sort=False, dropna=False):
        k = k if isinstance(k, tuple) else (k,)
        a, b = g[g["hemisphere"].astype(str) == "A"], g[g["hemisphere"].astype(str) == "B"]
        if len(a) == 1 and len(b) == 1:
            rows.append({"group": " / ".join(str(x) for x in k[: len(facs)]), "section": str(k[-1]),
                         "A": float(a["_v"].iloc[0]), "B": float(b["_v"].iloc[0]),
                         "flagged": bool(a["_f"].iloc[0] or b["_f"].iloc[0])})
    return pd.DataFrame(rows, columns=["group", "section", "A", "B", "flagged"])


def hemisphere_plot(pairs: pd.DataFrame, label: str, out) -> object:
    """A vs B of each section (log-log, identity line), one colour per group; hollow = a flagged tumour in the pair.
    Exploratory: hemispheres are independent n's; this only shows whether they share a section effect."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .plots import INK, MUTED, SERIES, _plain_log, _save

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ok = pairs[(pairs["A"] > 0) & (pairs["B"] > 0)]
    markers = ["o", "s", "^", "D"]
    for i, (grp, g) in enumerate(ok.groupby("group", sort=False)):
        col, mk = SERIES[i % len(SERIES)], markers[(i // len(SERIES)) % len(markers)]
        f = g["flagged"].to_numpy(bool)
        ax.scatter(g["A"][~f], g["B"][~f], s=26, color=col, marker=mk, label=f"{grp} (n={len(g)})", zorder=3)
        if f.any():
            ax.scatter(g["A"][f], g["B"][f], s=26, facecolors="none", edgecolors=col, marker=mk, zorder=3)
    if len(ok):
        lo, hi = float(ok[["A", "B"]].min().min()) / 1.2, float(ok[["A", "B"]].max().max()) * 1.2
        ax.plot([lo, hi], [lo, hi], color=MUTED, linestyle="--", linewidth=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        _plain_log(ax, "x")
        _plain_log(ax, "y")
        r = float(np.corrcoef(np.log(ok["A"]), np.log(ok["B"]))[0, 1]) if len(ok) > 2 else NAN
        ax.set_title(f"{label}: hemisphere A vs B of the same section ({len(ok)} sections, r(log) = {r:.2f})",
                     fontsize=8.5, color=INK)
    ax.set_xlabel(f"{label}, hemisphere A")
    ax.set_ylabel(f"{label}, hemisphere B")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(fontsize=6.5, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.suptitle("EXPLORATORY check, not a test: each hemisphere is its own n; hollow = flagged tumour in the pair",
                 fontsize=7.5, color=MUTED)
    fig.tight_layout()
    return _save(fig, out)


def _md_table(df: pd.DataFrame) -> list[str]:
    from .sensitivity import md_cell

    lines = ["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
    lines += ["| " + " | ".join(md_cell(str(v)) for v in r) + " |" for r in df.itertuples(index=False)]
    return lines


def _corrections_lines(out) -> list[str]:
    """Report lines on the review-corrections generation (corrections_used.csv written by the run)."""
    from pathlib import Path

    from .corrections import USED_FILE

    p = Path(out) / USED_FILE
    if not p.exists():
        return ["Review corrections: none applied (no corrections_used.csv)."]
    u = pd.read_csv(p, dtype=str, keep_default_na=False)
    gen = u["generation"].iloc[0] if len(u) else "g1"
    changed = u[u["source"] != "auto"]
    lines = [f"Review corrections: generation **{gen}** "
             + ("(automatic; nothing corrected)." if not len(changed) else
                f"— {len(changed)} tumour(s) corrected at blinded review (main mass re-picked or excluded; "
                "readouts recomputed from the cached labels, never re-segmented):")]
    lines += [f"- {r['tumour_id']}: {r['source']}" for _, r in changed.iterrows()]
    return lines


def analyse(out, table: pd.DataFrame, settings: dict, facs, fams, report, lock_hash: str) -> list[str]:
    """QC outputs (ticket 06): flow.csv + 'QC flow' sheet, 'Without flagged' sheet (primary beside the analysis
    without flagged tumours) + stats_long rows, A-B hemisphere plots, and the report section."""
    from pathlib import Path

    from ..manifest import STATS_LONG_COLUMNS
    from .config import readouts
    from .layout import EXCLUDED

    out = Path(out)
    delim = settings.get("project", {}).get("csv_delimiter", ",")
    ex_path = out / EXCLUDED
    excluded = (pd.read_csv(ex_path, sep=delim, dtype=str, keep_default_na=False) if ex_path.exists()
                else pd.DataFrame(columns=["tumour_id", "exclusion_reason"]))
    flow = flow_table(table, excluded, facs)
    flow.to_csv(out / FLOW_FILE, index=False, sep=delim)
    flg = flagged_mask(table)
    ros = [r for r in readouts(settings) if r.column in table.columns]
    sides = [without_flagged(table, facs, fams, r) for r in ros] if flg.any() else []
    side = pd.concat(sides, ignore_index=True) if sides else pd.DataFrame(columns=SIDE_COLUMNS)
    long = [row for r, s in zip(ros, sides) for row in _long_rows(s, r, lock_hash)]
    if long:
        st = pd.read_csv(out / "stats_long.csv")
        pd.concat([st, pd.DataFrame(long, columns=STATS_LONG_COLUMNS)], ignore_index=True).to_csv(
            out / "stats_long.csv", index=False)
    written = []
    for r in ros:
        if r.role != "confirmatory":
            continue
        pairs = hemisphere_pairs(table, facs, r.column)
        if len(pairs):
            written.append(str(hemisphere_plot(pairs, r.label, out / "plots" / f"{r.name}_hemispheres_AB")))
    with pd.ExcelWriter(out / "results.xlsx", engine="openpyxl", mode="a") as xw:
        flow.to_excel(xw, sheet_name="QC flow", index=False)
        side.to_excel(xw, sheet_name="Without flagged", index=False)
        qc_cols = [c for c in ["tumour_id", *QC_COLUMNS] if c in table.columns]
        if len(qc_cols) > 1:
            table[qc_cols].to_excel(xw, sheet_name="QC", index=False)
    md = ["", "## QC and problem tumours", "",
          "Pre-set rules (ADR 0002): a tumour is **truncated** when its main mass touches the xy image edge or its "
          "mid-section does not fit inside the imaged layers (extending beyond the imaged depth is not a flag — every "
          "tumour does). Truncated tumours are **flagged**: kept in the primary analysis, drawn as hollow markers in "
          "every group plot, and dropped in the sensitivity analysis below. Saturation of the mid-section main mass "
          "above the settings threshold is a warning, not a flag. Per-tumour thresholds are never changed.", "",
          *_corrections_lines(out), "",
          "### Flow per group (files → tumours → excluded → analysed → flagged → analysed without flagged)", "",
          *_md_table(flow), ""]
    if flg.any():
        lst = table[flg]
        md += [f"### Flagged tumours ({int(flg.sum())})", ""]
        md += [f"- {r['tumour_id']}: {r.get('qc_flags', '') or 'flagged'}" for _, r in lst.iterrows()]
        md += ["", "### Primary analysis beside the analysis without flagged tumours", "",
               "Sheet *Without flagged* of results.xlsx and stats_long families `<id> (without flagged)`; never "
               "starred. Tested contrasts:", "",
               "| readout | family | comparison | primary [95 % CI] | p (Holm) | without flagged [95 % CI] | "
               "p (Holm) | N → N | conclusion changed |", "|---|---|---|---|---|---|---|---|---|"]
        for _, r in side[side["p_adj"].notna()].iterrows():
            def ci(e, lo, hi):
                return f"{e:.2f} [{lo:.2f}, {hi:.2f}]" if np.isfinite(e) else "—"
            wf_p = f"{r['without_flagged_p_adj']:.3g}" if np.isfinite(r["without_flagged_p_adj"]) else ""
            md.append(f"| {r['readout']} | {r['family']} | {r['comparison'].replace('|', '/')} | "
                      f"{ci(r['estimate'], r['ci_low'], r['ci_high'])} {r['significance']} | {r['p_adj']:.3g} | "
                      f"{ci(r['without_flagged_estimate'], r['without_flagged_ci_low'], r['without_flagged_ci_high'])}"
                      f" | {wf_p} | {r['N']} → {r['without_flagged_N']} | "
                      f"{'**yes**' if r['conclusion_changed'] is True else 'no' if r['conclusion_changed'] is False else ''} |")
    else:
        md += ["no tumour is flagged: the sensitivity analysis without flagged tumours equals the primary analysis."]
    if "saturation_warning" in table.columns:
        sat = table[flagged_mask(table.assign(flagged=table["saturation_warning"]))]
        if len(sat):
            md += ["", f"Saturation warnings ({len(sat)}): " + "; ".join(
                f"{r['tumour_id']} ({r['saturated_frac_mid_section']:.0%})" for _, r in sat.iterrows())]
    if written:
        md += ["", "Hemisphere check (exploratory, not a test): A vs B of each section in "
               + ", ".join(Path(w).name for w in written) + "."]
    with Path(report).open("a", encoding="utf-8") as fh:
        fh.write("\n".join(md) + "\n")
    return [str(out / FLOW_FILE), *written]


# ---------------------------------------------------------------- per-tumour pass: QC, corrections, review sheets

QC_VERSION = 2  # 2: edge-contact wording "<1%" for a few touching pixels (ticket 11)


def measure(project_dir, layout: pd.DataFrame, settings: dict, measurements: pd.DataFrame, corr,
            render: bool = True) -> pd.DataFrame:
    """Per included tumour: apply a re-picked main mass to its cached labels, recompute the main-mass readouts
    (only when re-picked), compute QC and draw the blinded review image. Returns `measurements` with the QC columns
    (+ `main_mass_object`) and corrected readouts. Cached in `_work/tumour/<slug>.qc.json`; the original stack is
    read only when something changed. Never re-segments."""
    import hashlib
    import json
    import time
    from pathlib import Path

    import tifffile

    from . import corrections, review
    from .measure import WORK, labels_path, raw_dir, slug
    from .segment import read_stack

    prj = Path(project_dir)
    ts = tumour_settings(settings)
    included = layout[layout["include"].astype(str).str.strip().str.lower() == "yes"]
    m = measurements.set_index("tumour_id")
    tids = [t for t in included["tumour_id"] if t in m.index]
    n_obj = {t: int(m.at[t, "n_satellites"]) + 1 if "n_satellites" in m.columns and pd.notna(m.at[t, "n_satellites"])
             else 1 for t in tids}
    corrections.check_picks(corr, n_obj)
    code_of = review.codes(prj, list(included["tumour_id"]))
    corr.codes.update(code_of)
    px, quality = review.options(settings)
    rd = review.review_dir(prj)
    opts = {k: ts.get(k) for k in ("qc_saturation_warn", "mid_section_um")}
    rows, overrides = [], {}
    raw = raw_dir(prj, settings)
    for i, (_, row) in enumerate(included.iterrows(), 1):
        tid = row["tumour_id"]
        if tid not in m.index:
            continue
        pick = corr.picks.get(tid, (1, ""))[0]
        meta = prj / WORK / f"{slug(tid)}.json"
        mkey = json.loads(meta.read_text()).get("key", "") if meta.exists() else ""
        qkey = hashlib.sha1(json.dumps({"m": mkey, "pick": pick, "qc": opts, "v": QC_VERSION},
                                       sort_keys=True, default=str).encode()).hexdigest()[:16]
        rsig = f"{qkey}:{code_of[tid]}:{px}:{quality}"
        cache_p = prj / WORK / f"{slug(tid)}.qc.json"
        cached: dict = json.loads(cache_p.read_text()) if cache_p.exists() else {}
        img = rd / f"{code_of[tid]}.jpg"
        fresh = cached.get("key") == qkey
        need_img = render and not (img.exists() and fresh and cached.get("review") == rsig)
        res: dict
        if fresh and not need_img:
            res = cached
        else:
            t0 = time.time()
            vox = (float(m.at[tid, "voxel_z_um"]), float(m.at[tid, "voxel_y_um"]), float(m.at[tid, "voxel_x_um"]))
            labels = corrections.relabel(tifffile.imread(labels_path(prj, tid)), pick)
            stack = read_stack(raw / row["file"])
            q, centre = assess(labels, stack, vox, settings)
            res = {"key": qkey, "qc": q,
                   "readouts": corrections.main_mass_readouts(labels, stack, vox, ts) if pick != 1 else {},
                   "review": cached.get("review", "") if fresh else ""}
            if render:
                review.render(stack, labels, centre, code_of[tid], img, (vox[1], vox[2]), px, quality)
                res["review"] = rsig
            del labels, stack
            cache_p.write_text(json.dumps(res, indent=1, default=float))
            print(f"  QC [{i}/{len(included)}] {tid}: {q['qc_flags'] or 'no flags'} ({time.time() - t0:.0f} s)",
                  flush=True)
        rows.append({"tumour_id": tid, "main_mass_object": pick, **res["qc"]})
        if res["readouts"]:
            overrides[tid] = res["readouts"]
    out = measurements.copy().set_index("tumour_id")
    for tid, ro in overrides.items():
        for k, v in ro.items():
            if k not in out.columns:
                out[k] = np.nan
            out[k] = out[k].astype(object) if isinstance(v, str) else out[k]
            out.at[tid, k] = v
    qc_t = pd.DataFrame(rows).set_index("tumour_id") if rows else pd.DataFrame()
    for c in qc_t.columns:
        out[c] = qc_t[c].reindex(out.index)
    if render and rows:
        review.overviews(prj, [code_of[t] for t in qc_t.index])
    review.write_corrections_template(
        prj, settings, pd.DataFrame({"code": [code_of[t] for t in qc_t.index],
                                     "objects": [n_obj[t] for t in qc_t.index],
                                     "qc_flags": list(qc_t["qc_flags"]) if len(qc_t) else []}))
    return out.reset_index()
