"""Turn segmented objects into a count, a plausible range, and a QC verdict.

Categories (per object):
    nucleus        counted
    merged         counted as 1, but probably >1 touching nuclei (raises upper bound)
    debris         too small (micronuclei, fragments, dust)        - not counted
    out_of_focus   blurry, outside focal plane / floating            - not counted (configurable)
    edge           cut by image border, excluded by the edge rule  - not counted
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

COUNTED = ("nucleus", "merged")


@dataclass
class FolderRef:
    median_area: float
    median_focus: float
    median_texture: float = float("nan")


def folder_reference(objs: list[pd.DataFrame]) -> FolderRef:
    allo = pd.concat([o for o in objs if len(o)], ignore_index=True) if any(len(o) for o in objs) else pd.DataFrame()
    if allo.empty:
        return FolderRef(np.nan, np.nan, np.nan)
    med = float(allo["area"].median())
    # typical focus of plausible nuclei (not debris)
    big = allo[allo["area"] >= 0.5 * med]
    tex = float(big["texture"].median()) if "texture" in big else np.nan
    return FolderRef(float(big["area"].median()), float(big["focus"].median()), tex)


def _edge_excluded(df: pd.DataFrame, rule: str) -> pd.Series:
    if rule == "include":
        return pd.Series(False, index=df.index)
    if rule == "exclude":
        return df[["edge_top", "edge_left", "edge_bottom", "edge_right"]].any(axis=1)
    return df["edge_bottom"] | df["edge_right"]  # unbiased counting frame


def classify(objects: pd.DataFrame, info: dict, ref: FolderRef, settings: dict) -> tuple[pd.DataFrame, dict]:
    c = settings["counting"]
    q = settings["qc"]
    df = objects.copy()
    if df.empty or np.isnan(ref.median_area):
        return df.assign(category=pd.Series(dtype=str)), dict(
            auto_count=0, count_low=0, count_high=0, uncertainty=0.0, disagreement=0.0,
            n_debris=0, n_out_of_focus=0, n_edge=0, n_merged=0, saturated_frac=0.0, qc_level="low",
            qc_ok=False, qc_flags="no nuclei found: check image (blank / defocused field?)")

    area_thr = c["min_area_frac"] * ref.median_area
    focus_thr = c["focus_frac"] * ref.median_focus
    debris = df["area"] < area_thr
    # out of focus = blurred edges AND no chromatin texture (glare/halo lowers edge sharpness of sharp nuclei,
    # but not their texture); saturated nuclei are never called out of focus
    tex_thr = c["texture_frac"] * ref.median_texture if "texture" in df and not np.isnan(ref.median_texture) else np.inf
    low_tex = df["texture"] < tex_thr if "texture" in df else pd.Series(True, index=df.index)
    blurry = ((df["focus"] < focus_thr) & low_tex & (df["saturated_frac"] < 0.2) & ~debris
              if c["exclude_out_of_focus"] else pd.Series(False, index=df.index))
    edge = _edge_excluded(df, c["edge_rule"]) & ~debris & ~blurry
    merged = (df["area"] > c["merge_area_frac"] * ref.median_area) & (df["solidity"] < 0.9) & ~debris & ~blurry & ~edge
    cat = np.select([debris, blurry, edge, merged], ["debris", "out_of_focus", "edge", "merged"], "nucleus")
    df["category"] = cat
    counted = df["category"].isin(COUNTED)
    count = int(counted.sum())

    # ---- plausible range from borderline decisions
    near_area = df["area"].between(0.75 * area_thr, 1.5 * area_thr)
    # borderline focus call: blurred-ish edges (0.75-1.33 x threshold) with texture near its threshold,
    # i.e. an object that a small change in either threshold would move in or out of the count
    if c["exclude_out_of_focus"]:
        tex = df["texture"] if "texture" in df else pd.Series(0.0, index=df.index)
        edge_border = df["focus"].between(0.75 * focus_thr, 1.33 * focus_thr)
        tex_border = tex.between(0.75 * tex_thr, 1.25 * tex_thr) if np.isfinite(tex_thr) else pd.Series(False, index=df.index)
        near_focus = (edge_border & (tex < 1.25 * tex_thr)) | ((df["focus"] < 1.33 * focus_thr) & tex_border)
    else:
        near_focus = pd.Series(False, index=df.index)
    low_prob = df["prob"].fillna(1.0) < 0.6
    doubtful_in = counted & (near_area | near_focus | low_prob)
    doubtful_out = (~counted) & (df["category"] != "edge") & (near_area | near_focus)
    extra_merge = int(np.clip(np.round(df.loc[df["category"] == "merged", "area"] / ref.median_area) - 1, 0, None).sum())
    lo = count - int(doubtful_in.sum())
    hi = count + int(doubtful_out.sum()) + extra_merge
    uncertainty = (hi - lo) / max(count, 5)

    # ---- second opinion: classical count with the same size and edge rules (no focus rule)
    disagreement = np.nan
    if info.get("alt_objects"):
        alt = pd.DataFrame(info["alt_objects"], columns=["area", "edge_top", "edge_left", "edge_bottom", "edge_right"]).astype(
            {"edge_top": bool, "edge_left": bool, "edge_bottom": bool, "edge_right": bool})
        alt_count = int(((alt["area"] >= area_thr) & ~_edge_excluded(alt, c["edge_rule"])).sum())
        ref_count = count + int((df["category"] == "out_of_focus").sum())
        disagreement = abs(alt_count - ref_count) / max(ref_count, 5)
    elif info.get("alt_objects") == []:
        disagreement = count / max(count, 5)

    sat = float(info.get("saturated_frac", 0.0))
    flags = []
    if uncertainty > q["high_max_uncertainty"]:
        flags.append(f"range {lo}-{hi}")
    if extra_merge:
        flags.append(f"{int((df['category'] == 'merged').sum())} possible merged")
    n_blur = int((df["category"] == "out_of_focus").sum())
    if n_blur:
        flags.append(f"{n_blur} out of focus excluded")
    if not np.isnan(disagreement) and disagreement > q["high_max_disagreement"]:
        flags.append(f"methods disagree {disagreement:.0%}")
    if sat > q["max_saturated_frac"]:
        flags.append(f"saturated {sat:.0%}")

    dis = 0.0 if np.isnan(disagreement) else disagreement
    if count == 0:
        flags.append("no in-focus nuclei: check image (blank / defocused field?)")
        level = "low"
    elif uncertainty >= q["low_min_uncertainty"] or dis >= q["low_min_disagreement"]:
        level = "low"
    elif uncertainty <= q["high_max_uncertainty"] and dis <= q["high_max_disagreement"] and sat <= q["max_saturated_frac"]:
        level = "high"
    else:
        level = "medium"

    return df, dict(
        auto_count=count, count_low=max(lo, 0), count_high=hi, uncertainty=round(uncertainty, 3),
        disagreement=None if np.isnan(disagreement) else round(float(disagreement), 3),
        n_debris=int((df["category"] == "debris").sum()), n_out_of_focus=n_blur,
        n_edge=int((df["category"] == "edge").sum()), n_merged=int((df["category"] == "merged").sum()),
        saturated_frac=round(sat, 4), qc_level=level, qc_ok=level != "low", qc_flags="; ".join(flags))
