"""Nucleus segmentation (StarDist or classical) and per-object measurement."""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import feature, filters, measure, morphology, segmentation

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")

MODEL_DIR = Path(__file__).parent / "models"


@dataclass
class Result:
    labels: np.ndarray
    objects: pd.DataFrame
    info: dict = field(default_factory=dict)


# ---------------------------------------------------------------- helpers

def _norm(img: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(img, [1, 99.8])
    return (img - lo) / max(hi - lo, 1e-6)


def _threshold_objects(img: np.ndarray, sigma: float) -> np.ndarray:
    g = filters.gaussian(_norm(img), sigma)
    try:
        t = filters.threshold_otsu(g)
    except ValueError:
        return np.zeros(img.shape, bool)
    return ndi.binary_fill_holes(g > t)


def estimate_diameter(img: np.ndarray) -> float | None:
    """Typical nucleus diameter (px) from a quick threshold; robust to clumps via convex objects."""
    lab = measure.label(_threshold_objects(img, 2))
    props = [p for p in measure.regionprops(lab) if p.area > 30]
    areas = [p.area for p in props if p.solidity > 0.9]
    if len(areas) < 3:
        areas = [p.area for p in props]
    if len(areas) < 1:
        return None
    return float(np.sqrt(np.median(areas) * 4 / np.pi))


# ---------------------------------------------------------------- backends

_MODEL = None


def _stardist_model():
    global _MODEL
    if _MODEL is None:
        from stardist.models import StarDist2D

        local = MODEL_DIR / "2D_versatile_fluo"
        if (local / "config.json").exists():
            _MODEL = StarDist2D(None, name=local.name, basedir=str(local.parent))
        else:
            _MODEL = StarDist2D.from_pretrained("2D_versatile_fluo")
    return _MODEL


def segment_stardist(img: np.ndarray, diameter: float, target: float, prob_thresh) -> tuple[np.ndarray, dict]:
    from csbdeep.utils import normalize

    model = _stardist_model()
    scale = float(np.clip(target / diameter, 0.2, 3.0))
    kw = {} if prob_thresh in (None, "", "default") else {"prob_thresh": float(prob_thresh)}
    labels, det = model.predict_instances(normalize(img, 1, 99.8), scale=scale, show_tile_progress=False, verbose=False, **kw)
    labels = labels.astype(np.int32)
    prob = np.asarray(det.get("prob", []), dtype=float)
    # StarDist labels are 1..N in the order of det["prob"].
    probs = {i + 1: float(p) for i, p in enumerate(prob)}
    return labels, {"scale": scale, "probs": probs}


def segment_classical(img: np.ndarray, diameter: float) -> np.ndarray:
    """Fiji-like pipeline, but on full bit depth: blur, background removal, Otsu, distance watershed."""
    sigma = max(diameter / 25, 1.0)
    norm = _norm(img)
    bg = filters.gaussian(norm, diameter * 1.5)
    g = filters.gaussian(norm - bg, sigma)
    try:
        t = filters.threshold_otsu(g)
    except ValueError:
        return np.zeros(img.shape, np.int32)
    mask = ndi.binary_fill_holes(g > t)
    mask = ndi.binary_opening(mask, structure=morphology.disk(max(int(diameter / 15), 1)))
    dist = ndi.distance_transform_edt(mask)
    dist = filters.gaussian(dist, sigma)
    peaks = feature.peak_local_max(dist, min_distance=max(int(diameter / 3), 2), threshold_abs=diameter / 8,
                                   labels=measure.label(mask), exclude_border=False)
    markers = np.zeros(mask.shape, np.int32)
    markers[tuple(peaks.T)] = np.arange(1, len(peaks) + 1)
    labels = segmentation.watershed(-dist, markers, mask=mask)
    # objects without a peak (small blobs) still count
    rest = measure.label(mask & (labels == 0))
    rest[rest > 0] += labels.max()
    return (labels + rest).astype(np.int32)


# ---------------------------------------------------------------- over-split repair

def merge_oversplit(labels: np.ndarray, img: np.ndarray, expected_area: float, probs: dict | None,
                    frag_frac: float = 0.35, seam_ratio: float = 0.8) -> tuple[np.ndarray, dict | None, int]:
    """StarDist assumes star-convex shapes, so lobed, kidney-shaped and fragmented nuclei get split.
    Merge touching objects when
      - one of them is a fragment (< frag_frac x expected nucleus area) and the union is not oversized, or
      - the union is nucleus-sized and the seam between them is about as bright as their interiors
        (two separate touching nuclei show a darker seam)."""
    n = int(labels.max())
    if n < 2:
        return labels, probs, 0
    sm = filters.gaussian(_norm(img), 1)
    idx = np.arange(1, n + 1)
    area = np.bincount(labels.ravel(), minlength=n + 1).astype(float)
    inside = np.zeros(n + 1)
    inside[1:] = np.nan_to_num(ndi.mean(sm, labels=labels, index=idx))

    contact: dict[tuple[int, int], list[float]] = {}
    H, W = labels.shape

    def pair(arr, dy, dx):  # arr[y, x] and arr[y + dy, x + dx] over the valid overlap (dy >= 0)
        xa = slice(0, W - dx) if dx >= 0 else slice(-dx, W)
        xb = slice(dx, W) if dx >= 0 else slice(0, W + dx)
        return arr[0: H - dy, xa], arr[dy:H, xb]

    for dy, dx in ((0, 1), (1, 0), (0, 2), (2, 0), (1, 1), (1, -1)):
        a, b = pair(labels, dy, dx)
        sa, sb = pair(sm, dy, dx)
        m = (a > 0) & (b > 0) & (a != b)
        if not m.any():
            continue
        lo, hi = np.minimum(a[m], b[m]), np.maximum(a[m], b[m])
        val = np.minimum(sa[m], sb[m])
        if (dy, dx) == (0, 2):
            val = np.minimum(val, sm[:, 1: W - 1][m])  # the gap pixel between them
        elif (dy, dx) == (2, 0):
            val = np.minimum(val, sm[1: H - 1, :][m])
        for p, q, v in zip(lo.tolist(), hi.tolist(), val.tolist()):
            contact.setdefault((p, q), []).append(v)
    if not contact:
        return labels, probs, 0

    parent = list(range(n + 1))

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    merged = 0
    for (p, q), vals in sorted(contact.items(), key=lambda kv: -len(kv[1])):
        rp, rq = root(p), root(q)
        if rp == rq:
            continue
        union = area[rp] + area[rq]
        small = min(area[p], area[q]) < frag_frac * expected_area
        seam = float(np.mean(vals)) / max(min(inside[p], inside[q]), 1e-6)
        if (small and union <= 1.6 * expected_area) or (union <= 1.4 * expected_area and seam >= seam_ratio):
            parent[rq] = rp
            area[rp] = union
            inside[rp] = (inside[rp] + inside[rq]) / 2
            merged += 1
    if not merged:
        return labels, probs, 0
    lut = np.array([root(i) for i in range(n + 1)])
    new, fwd, _ = segmentation.relabel_sequential(lut[labels])
    new_probs = None
    if probs:
        new_probs = {}
        for old, pr in probs.items():
            k = int(fwd[lut[old]]) if old <= n else 0
            if k:
                new_probs[k] = max(new_probs.get(k, 0.0), pr)
    return new.astype(np.int32), new_probs, merged


# ---------------------------------------------------------------- measurement

def _focus_scores(norm: np.ndarray, labels: np.ndarray, n: int) -> np.ndarray:
    """Boundary sharpness: mean gradient on the outline divided by nucleus-vs-surround contrast."""
    grad = filters.sobel(filters.gaussian(norm, 1))
    bnd = segmentation.find_boundaries(labels, mode="inner")
    bnd[0, :] = bnd[-1, :] = bnd[:, 0] = bnd[:, -1] = False  # image border is not a nucleus edge
    idx = np.arange(1, n + 1)
    g_edge = ndi.mean(grad, labels=np.where(bnd, labels, 0), index=idx)
    inside = ndi.mean(norm, labels=labels, index=idx)
    ring = segmentation.expand_labels(labels, 4) * (labels == 0)
    outside = ndi.mean(norm, labels=ring, index=idx)
    contrast = np.maximum(np.asarray(inside) - np.nan_to_num(np.asarray(outside)), 1e-3)
    return np.asarray(g_edge) / contrast


def _texture_scores(norm: np.ndarray, labels: np.ndarray, n: int, diameter: float = 30.0) -> np.ndarray:
    """Chromatin texture: SD of a Laplacian-of-Gaussian inside the nucleus (rim excluded), with the
    same statistic measured on empty background (camera noise) subtracted, relative to
    nucleus-vs-surround contrast. Defocus removes chromatin texture; halos, glare and noise do not add it."""
    log = ndi.gaussian_laplace(norm, max(diameter / 20, 1.0))
    rim = ndi.binary_dilation(segmentation.find_boundaries(labels, mode="inner"), iterations=2)
    interior = np.where(rim, 0, labels)
    idx = np.arange(1, n + 1)
    sd = np.nan_to_num(np.asarray(ndi.standard_deviation(log, labels=interior, index=idx)))
    bg = segmentation.expand_labels(labels, int(diameter / 2)) == 0
    noise = float(np.std(log[bg])) if bg.sum() > 1000 else 0.0
    sd = np.sqrt(np.maximum(sd ** 2 - noise ** 2, 0.0))
    inside = np.asarray(ndi.mean(norm, labels=labels, index=idx))
    ring = segmentation.expand_labels(labels, 4) * (labels == 0)
    outside = np.nan_to_num(np.asarray(ndi.mean(norm, labels=ring, index=idx)))
    return sd / np.maximum(inside - outside, 1e-3)


def measure_objects(channels: np.ndarray, labels: np.ndarray, seg_ch: int, ring_px: float,
                    max_value: float, probs: dict | None, diameter: float = 30.0) -> tuple[pd.DataFrame, dict]:
    n = int(labels.max())
    img = channels[seg_ch]
    H, W = labels.shape
    info: dict = {}
    if n == 0:
        return pd.DataFrame(), info
    props = measure.regionprops(labels)
    rows = []
    for p in props:
        r0, c0, r1, c1 = p.bbox
        rows.append(dict(label=p.label, area=p.area, solidity=p.solidity, eccentricity=p.eccentricity,
                         y=p.centroid[0], x=p.centroid[1],
                         edge_top=r0 == 0, edge_left=c0 == 0, edge_bottom=r1 == H, edge_right=c1 == W))
    df = pd.DataFrame(rows).set_index("label")
    lab_idx = df.index.to_numpy()
    df["prob"] = [probs.get(int(i), np.nan) for i in lab_idx] if probs else np.nan
    full_focus = _focus_scores(_norm(img), labels, n)
    df["focus"] = full_focus[lab_idx - 1]
    df["texture"] = _texture_scores(_norm(img), labels, n, diameter)[lab_idx - 1]
    sat = ndi.mean((img >= max_value * 0.999).astype(float), labels=labels, index=lab_idx)
    df["saturated_frac"] = sat
    info["saturated_frac"] = float(((img >= max_value * 0.999) & (labels > 0)).sum() / max((labels > 0).sum(), 1))

    ring_px = max(int(round(ring_px)), 1)
    cell = segmentation.expand_labels(labels, ring_px)
    ring = np.where(labels == 0, cell, 0)
    far = segmentation.expand_labels(labels, ring_px * 3) == 0
    df["cell_area"] = ndi.sum(np.ones(labels.shape), labels=cell, index=lab_idx)
    for c in range(channels.shape[0]):
        ch = channels[c]
        bg = float(np.median(ch[far])) if far.any() else 0.0
        info[f"C{c + 1}_background"] = bg
        df[f"C{c + 1}_nucleus_mean"] = ndi.mean(ch, labels=labels, index=lab_idx)
        df[f"C{c + 1}_ring_mean"] = ndi.mean(ch, labels=ring, index=lab_idx)
        df[f"C{c + 1}_cell_mean"] = ndi.mean(ch, labels=cell, index=lab_idx)
    info["ring_px"] = ring_px
    return df.reset_index(), info


def run(channels: np.ndarray, max_value: float, seg_ch: int, settings: dict,
        diameter: float | None) -> Result:
    s = settings["segmentation"]
    img = channels[seg_ch]
    est = estimate_diameter(img)
    d = float(s["nucleus_diameter_px"]) if str(s["nucleus_diameter_px"]) != "auto" else (diameter or est or 30.0)

    if s["backend"] == "stardist":
        labels, extra = segment_stardist(img, d, float(s["target_diameter_px"]), s.get("prob_thresh"))
        probs, scale = extra["probs"], extra["scale"]
        labels, probs, n_merged = merge_oversplit(labels, img, np.pi * (d / 2) ** 2, probs)
    else:
        labels, probs, scale = segment_classical(img, d), None, 1.0
    # Independent second opinion used for QC (always classical).
    alt = segment_classical(img, d) if s["backend"] == "stardist" else None

    ring = settings["intensity"]["ring_px"]
    ring_px = 0.25 * d if str(ring) == "auto" else float(ring)
    objects, info = measure_objects(channels, labels, seg_ch, ring_px, max_value, probs, d)
    info.update(diameter_used=d, diameter_estimate=est, scale=scale, backend=s["backend"],
                oversplit_merges=n_merged if s["backend"] == "stardist" else 0)
    if alt is not None:
        H, W = alt.shape
        info["alt_objects"] = [
            [int(p.area), int(p.bbox[0] == 0), int(p.bbox[1] == 0), int(p.bbox[2] == H), int(p.bbox[3] == W)]
            for p in measure.regionprops(alt)
        ]
    return Result(labels, objects, info)
