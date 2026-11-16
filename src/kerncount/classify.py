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


# ---------------------------------------------------------------- per-image readouts
# Computed from the classified objects of one image: no re-segmentation and no pixel intensities, so they are
# fast and never invalidate the segmentation cache. Units are relative wherever possible (fractions, multiples of
# the nucleus diameter) so images taken at different magnifications stay comparable. See docs/readouts.md.

class _Readouts(dict):
    """The static readouts below, plus per-channel intensity columns whose names contain the channel
    (docs/intensity.md): `READOUTS.get(column)` answers for them too, so metrics.registry labels them correctly.
    Iterating lists only the static readouts (the ones `readouts.include` selects from)."""

    def get(self, key, default=None):
        if dict.__contains__(self, key):
            return dict.__getitem__(self, key)
        info = channel_metric_info(str(key))
        return info if info is not None else default


def channel_metric_info(col: str) -> dict | None:
    """Label / kind / unit for the per-channel intensity columns written by intensity.process, else None."""
    import re

    pats = [
        (r"readout_frac_positive_p90_(.+)", lambda c: dict(
            label=f"% {c}-positive nuclei (90th percentile)", kind="fraction", unit="fraction",
            description=f"Share of measurable nuclei whose 90th-percentile {c} signal (minus local background, "
                        "flat-field corrected) exceeds the cut-off: picks up dotted / nucleolar staining that the mean "
                        "dilutes. Same cut-off rule as the primary readout (see intensity/cutoffs.csv).",
            min_nuclei=5, n_column=f"n_measured_{c}")),
        (r"readout_frac_positive_(.+)", lambda c: dict(
            label=f"% {c}-positive nuclei", kind="fraction", unit="fraction",
            description=f"Share of measurable nuclei (counted, not saturated) whose mean {c} signal minus the local "
                        "background, divided by the flat field, exceeds the positivity cut-off. Its rule (negative "
                        "control, locked manual value or limit of detection) is in intensity/cutoffs.csv; without "
                        "negative-control wells 'positive' means detectably stained above background.",
            min_nuclei=5, n_column=f"n_measured_{c}")),
        (r"intensity_(.+)_nucleus_mean_bgsub", lambda c: dict(
            label=f"{c} per nucleus (median, background-subtracted)", kind="positive", unit="intensity (a.u., corrected)",
            description=f"Median over the image's measurable nuclei of the mean {c} inside the shrunk nucleus outline "
                        "minus the local background, divided by the flat field. Ladder layer 4.")),
        # the readout ladder (docs/intensity.md "Readout ladder"): each layer removes one analysis choice
        (r"intensity_(.+)_field_raw", lambda c: dict(
            label=f"{c} field MFI, raw", kind="positive", unit="intensity (a.u., raw)",
            description=f"Mean {c} of the whole field as acquired: no flat-field correction, no background "
                        "subtraction, no segmentation (saturated pixels left out). Ladder layer 1.")),
        (r"intensity_(.+)_field_corrected", lambda c: dict(
            label=f"{c} field MFI, corrected", kind="positive", unit="intensity (a.u., corrected)",
            description=f"Mean {c} of the whole field after subtracting the field's background (histogram mode, "
                        "no DAPI) and dividing by the flat field; no segmentation. Ladder layer 2.")),
        (r"intensity_(.+)_nucleus_mean_raw", lambda c: dict(
            label=f"{c} per nucleus, raw (median; Fiji-like)", kind="positive", unit="intensity (a.u., raw)",
            description=f"Median over the image's measurable nuclei of the mean raw {c} over the whole DAPI nucleus "
                        "outline: no background subtraction, no flat field (comparable with a Fiji ROI mean). "
                        "Ladder layer 3.")),
        (r"intensity_(.+)_stained_area_mean", lambda c: dict(
            label=f"{c} MFI in the {c}-stained area", kind="positive", unit="intensity (a.u., corrected)",
            description=f"Mean corrected {c} over the pixels above an automatic threshold on the {c} image itself "
                        "(background mode + 3 noise SD of the smoothed corrected image; no DAPI).")),
        (r"intensity_(.+)_field_total", lambda c: dict(
            label=f"{c} total signal per field", kind="positive", unit="integrated intensity (a.u. x px, corrected)",
            description=f"Integrated corrected {c} over the {c}-stained area of the field (mean in the stained area "
                        "x its pixels): total stain per field, grows with cell number.")),
        (r"intensity_(.+)_positive_mean_bgsub", lambda c: dict(
            label=f"{c} in positive nuclei (mean)", kind="positive", unit="intensity (a.u., corrected)",
            description=f"Mean background-subtracted, flat-field-corrected {c} of the positive nuclei only: how "
                        "strongly the positive cells stain.")),
        (r"intensity_(.+)_fiji_cell", lambda c: dict(
            label=f"{c} mean, nucleus + 6 px (Fiji legacy)", kind="positive", unit="intensity (a.u.)",
            description=f"Legacy readout of the old Fiji macro: mean {c} of the nucleus enlarged by fiji_ring_px, minus "
                        "the image background, averaged over counted nuclei. No flat-field correction.")),
    ]
    for pat, make in pats:
        m = re.fullmatch(pat, col)
        if m:
            return make(m.group(1))
    return None


READOUTS: dict[str, dict] = _Readouts({
    "readout_frac_out_of_focus": dict(
        label="Out-of-focus fraction", kind="fraction", unit="fraction",
        description="Share of nuclei lying outside the focal plane (rounded-up, floating or detaching cells): "
                    "out of focus / (counted + out of focus).",
        min_nuclei=10),
    "readout_micronuclei_per_nucleus": dict(
        label="Micronucleus candidates per nucleus", kind="positive", unit="per nucleus",
        description="Small DAPI objects (at most 1/9 of a typical nucleus' area) lying right next to a nucleus, per "
                    "counted nucleus: a lower-bound proxy for micronuclei (DNA damage, mitotic errors); fragments "
                    "of a disintegrating (apoptotic) nucleus also score.",
        min_nuclei=10),
    "readout_nuclear_area_median": dict(
        label="Median nuclear area", kind="positive", unit="px²",
        description="Median area of whole, in-focus nuclei; compare it as a ratio to the reference within the same n, "
                    "because pixel sizes can differ between magnifications.",
        min_nuclei=5),
    "readout_nuclear_area_median_um2": dict(
        label="Median nuclear area (µm²)", kind="positive", unit="µm²",
        description="Median area of whole, in-focus nuclei in µm² (only present when the pixel size is known).",
        min_nuclei=5),
    "readout_nuclear_area_cv": dict(
        label="Nuclear size heterogeneity (CV)", kind="positive", unit="SD / mean",
        description="Coefficient of variation of the area of whole, in-focus nuclei: high when the population mixes "
                    "small (condensed, apoptotic) and large (polyploid, arrested) nuclei.",
        min_nuclei=10),
    "readout_frac_irregular": dict(
        label="Irregular-nucleus fraction", kind="fraction", unit="fraction",
        description="Share of whole, in-focus nuclei with a non-convex outline (solidity below 0.9): blebbed, lobed or "
                    "fragmenting nuclei.",
        min_nuclei=10),
    "readout_coverage": dict(
        # a 0-1 share, but small (a few %) and density-like, so it is analysed like a count (ratio to reference)
        label="Nuclear coverage", kind="positive", unit="fraction of image area",
        description="Fraction of the image area covered by in-focus nuclei (including nuclei cut by the border): a "
                    "confluence proxy that does not depend on the counting rules.",
        min_nuclei=0),
    "readout_nn_distance_median": dict(
        label="Nearest-neighbour distance", kind="positive", unit="nucleus diameters",
        description="Median distance from each counted nucleus to its closest neighbour, in multiples of the median "
                    "nucleus diameter: small = crowded or clustered, large = sparse.",
        min_nuclei=5),
})

IRREGULAR_SOLIDITY = 0.9  # same convexity cut-off the classifier uses for suspected merges
MICRONUCLEUS_MAX_AREA_FRAC = 1 / 9  # micronucleus diameter <= 1/3 of the main nucleus (Fenech criteria)
MICRONUCLEUS_MAX_DIST = 1.0  # centroid within 1 nucleus diameter of a nucleus centroid (~half a diameter off its rim)
_EDGES = ("edge_top", "edge_left", "edge_bottom", "edge_right")


def readout_names(include="all") -> list[str]:
    """Readout columns selected by settings `readouts.include` ("all" or a list; the readout_ prefix is optional)."""
    if include is None or include == "all" or (isinstance(include, list) and "all" in include):
        return list(READOUTS)
    want = {str(i) if str(i).startswith("readout_") else f"readout_{i}" for i in include}
    return [k for k in READOUTS if k in want]


def compute_readouts(objects: pd.DataFrame, image_shape: tuple[int, ...], pixel_size_um: float | None = None,
                     include="all") -> tuple[dict, pd.DataFrame]:
    """Per-image readouts from the classified objects of one image (output of `classify`).

    Uses only `category`, `area`, `solidity`, centroid `x`/`y` and the edge flags: no intensities.
    Returns (readouts, objects); objects gains `nn_distance_px` (counted nuclei: distance to the nearest in-focus
    nucleus) and `micronucleus_candidate` (True only for debris objects that pass the micronucleus rule).
    A readout whose `min_nuclei` is not met is NaN. `readout_nuclear_area_median_um2` is only returned with a pixel size.
    """
    names = readout_names(include)
    df = objects.copy()
    out = {k: np.nan for k in names if k != "readout_nuclear_area_median_um2" or pixel_size_um}
    if df.empty or "category" not in df:
        if "readout_coverage" in out:
            out["readout_coverage"] = 0.0
        return out, df.assign(nn_distance_px=pd.Series(dtype=float), micronucleus_candidate=pd.Series(dtype=bool))
    cat = df["category"]
    counted = cat.isin(COUNTED)
    in_focus = cat.isin(("nucleus", "merged", "edge"))  # real nuclei in the focal plane, cut by the border or not
    edge_cols = [c for c in _EDGES if c in df]
    touches = df[edge_cols].astype(bool).any(axis=1) if edge_cols else pd.Series(False, index=df.index)
    whole = (cat == "nucleus") & ~touches  # uncut single nuclei: the only ones whose size and shape mean something
    n_counted, n_whole = int(counted.sum()), int(whole.sum())
    if n_whole:
        area_ref = float(df.loc[whole, "area"].median())
    else:
        area_ref = float(df.loc[counted, "area"].median()) if n_counted else np.nan
    diam_ref = float(np.sqrt(4 * area_ref / np.pi)) if np.isfinite(area_ref) else np.nan

    # nearest neighbour: counted nuclei -> any in-focus nucleus (border-cut nuclei are still somebody's neighbour)
    df["nn_distance_px"] = np.nan
    pts = df.loc[in_focus, ["y", "x"]].to_numpy(float)
    tree = None
    if len(pts) >= 2:
        from scipy.spatial import cKDTree

        tree = cKDTree(pts)
        if n_counted:
            d, _ = tree.query(df.loc[counted, ["y", "x"]].to_numpy(float), k=2)
            df.loc[counted, "nn_distance_px"] = d[:, 1]

    # micronucleus candidates: small debris objects hugging a nucleus
    debris = cat == "debris"
    df["micronucleus_candidate"] = False
    if debris.any() and len(pts) and np.isfinite(diam_ref):
        from scipy.spatial import cKDTree

        dist, _ = cKDTree(pts).query(df.loc[debris, ["y", "x"]].to_numpy(float), k=1)
        df.loc[debris, "micronucleus_candidate"] = (
            (df.loc[debris, "area"].to_numpy(float) <= MICRONUCLEUS_MAX_AREA_FRAC * area_ref)
            & (dist <= MICRONUCLEUS_MAX_DIST * diam_ref))

    def ok(name: str, n: int) -> bool:
        return name in out and n > 0 and n >= READOUTS[name]["min_nuclei"]

    n_oof = int((cat == "out_of_focus").sum())
    if ok("readout_frac_out_of_focus", n_counted + n_oof):
        out["readout_frac_out_of_focus"] = n_oof / (n_counted + n_oof)
    if ok("readout_micronuclei_per_nucleus", n_counted):
        out["readout_micronuclei_per_nucleus"] = int(df["micronucleus_candidate"].sum()) / n_counted
    areas = df.loc[whole, "area"].astype(float)
    if ok("readout_nuclear_area_median", n_whole):
        out["readout_nuclear_area_median"] = float(areas.median())
    if ok("readout_nuclear_area_median_um2", n_whole):
        out["readout_nuclear_area_median_um2"] = float(areas.median()) * float(pixel_size_um) ** 2
    if ok("readout_nuclear_area_cv", n_whole):
        out["readout_nuclear_area_cv"] = float(areas.std(ddof=1) / areas.mean())
    if ok("readout_frac_irregular", n_whole) and "solidity" in df:
        out["readout_frac_irregular"] = float((df.loc[whole, "solidity"] < IRREGULAR_SOLIDITY).mean())
    if "readout_coverage" in out:
        h, w = image_shape[:2]
        out["readout_coverage"] = float(df.loc[in_focus, "area"].sum()) / float(h * w)
    if ok("readout_nn_distance_median", n_counted) and tree is not None and np.isfinite(diam_ref):
        out["readout_nn_distance_median"] = float(df.loc[counted, "nn_distance_px"].median()) / diam_ref
    return {k: round(v, 5) if isinstance(v, float) and np.isfinite(v) else v for k, v in out.items()}, df
