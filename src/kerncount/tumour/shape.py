"""Size, shape and burden readouts of a labelled tumour (ticket 05; ADR 0002/0003 amendments 2026-11-20).

Every tumour is ~0.5-1 mm wide but only ~80-150 µm of it is imaged, so the confirmatory readouts are 2D: the main
mass's **cross-sectional area** and **perimeter**, averaged over the layers of a fixed-thickness **mid-section**
centred on the layer where the main mass has its largest cross-section (docs/tumour.md § Readouts).

- Section of the main mass in one layer = all its pieces in that layer with enclosed holes filled (2D), i.e. the area
  inside its outer outline(s). A dark core that opens to the top/bottom of the stack is still inside the outline.
- Area = pixel count x pixel area (unbiased for a digitised shape).
- Perimeter = length of the outline at one fixed physical smoothing scale: the binary section is Gaussian-smoothed
  (sigma `outline_smooth_um`, default 8 µm) and its 0.5 iso-line traced with marching squares (sub-pixel). Counting
  pixel edges instead overestimates any curved outline by 4/pi (27 %) whatever the resolution; marching squares on
  the raw binary still follows the staircase (~+5 %). Smoothed, digitised discs (r >= 80 µm) and ellipses come within
  0.5 % (tests/test_tumour_shape.py). Why 8 µm: real tumour borders are ragged at the µm scale, so a fine-scale
  perimeter measures border texture and moves with the threshold (+-11-17 % for a +-20 % threshold at 2 µm, <= 10 %
  at 8 µm; docs/tumour.md). The smoothing also makes masks of different origin (this method, Fiji's) comparable.
- Edge roughness = outline length at `roughness_fine_um` (2 µm) / at `roughness_smooth_um` (20 µm): the fine ragged
  border the confirmatory perimeter deliberately smooths away.

Mid-section: exactly `mid_section_um` thick for every z-step (review fix 2026-11-20; the real data mix 1.75, 2.0 and
2.36 µm). Layer k represents the slab from its centre - z/2 to its centre + z/2; the mid-section is the interval of
+-thickness/2 around the centre of the layer of largest cross-section, and each layer counts with the fraction of its
extent inside that interval (`MidSection.weights`: 1 inside, a fraction for the two outermost layers). Readouts are
weighted means, so 20 µm is 20 µm at 1.75 µm (13 layers, the outer two at 21 %) as at 2.36 µm (9 layers, the outer two
at 74 %). The slab volume uses the same rule with `slab_um`. When the interval does not fit inside the stack it is
clipped and `mid_section_fits` = "no"; `mid_section_um` records the thickness actually averaged over.

Pure numpy / scipy / scikit-image; nothing here knows about files or settings beyond `ShapeOptions`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage as ndi

NAN = math.nan


@dataclass(frozen=True)
class ShapeOptions:
    """Physical scales of the readouts ([tumour] settings; identical for every tumour)."""
    mid_section_um: float = 20.0  # thickness of the mid-section (µm)
    outline_smooth_um: float = 8.0  # outline smoothing for perimeters (µm, Gaussian sigma)
    roughness_fine_um: float = 2.0  # edge roughness = outline length at this scale ...
    roughness_smooth_um: float = 20.0  # ... / outline length at this (heavy) scale
    slab_um: float = 40.0  # depth of the slab for the slab volume (µm)
    surface_grid_um: float = 2.0  # isotropic grid for the (exploratory) 3D mesh surface area (µm)

    @classmethod
    def from_settings(cls, ts: dict) -> ShapeOptions:
        return cls(**{k: float(ts[k]) for k in cls.__dataclass_fields__ if k in ts})


# Keys of the [tumour] settings that change readouts (part of the measurement cache key).
KEYS = tuple(ShapeOptions.__dataclass_fields__)


@dataclass(frozen=True)
class MidSection:
    centre: int  # 0-based layer of the largest main-mass cross-section
    first: int  # 0-based, inclusive (clipped to the stack): first layer with a weight > 0
    last: int
    fits: bool  # the full thickness lies inside the imaged layers
    weights: np.ndarray = field(default_factory=lambda: np.zeros(0))  # per layer first..last, in (0, 1]
    z_um: float = 1.0

    @property
    def layers(self) -> range:
        return range(self.first, self.last + 1)

    @property
    def thickness_um(self) -> float:
        """Thickness actually covered (= the set thickness unless clipped by the stack)."""
        return float(np.sum(self.weights)) * self.z_um

    def weight(self, k: int) -> float:
        return float(self.weights[k - self.first]) if self.first <= k <= self.last else 0.0


def layer_weights(n_layers: int, centre: int, z_um: float, thickness_um: float) -> tuple[int, int, np.ndarray, bool]:
    """(first, last, weights, fits): fraction of each layer's extent (centre +- z/2) inside centre +- thickness/2."""
    eps = 1e-9
    lo, hi = centre * z_um - thickness_um / 2, centre * z_um + thickness_um / 2
    k = np.arange(n_layers)
    w = np.clip(np.minimum((k + 0.5) * z_um, hi) - np.maximum((k - 0.5) * z_um, lo), 0, None) / z_um
    w[w < eps] = 0.0
    w[w > 1 - eps] = 1.0
    nz = np.flatnonzero(w)
    first, last = (int(nz[0]), int(nz[-1])) if len(nz) else (centre, centre)
    fits = lo >= -0.5 * z_um - eps and hi <= (n_layers - 0.5) * z_um + eps
    return first, last, w[first:last + 1], bool(fits)


def mid_section(areas: np.ndarray, z_um: float, thickness_um: float) -> MidSection:
    """The mid-section for per-layer main-mass areas: centred on the largest area (on ties: the middle one)."""
    a = np.asarray(areas, float)
    ties = np.flatnonzero(a == a.max())
    k = int(ties[len(ties) // 2])
    first, last, w, fits = layer_weights(len(a), k, z_um, thickness_um)
    return MidSection(k, first, last, fits, w, z_um)


def weighted_mean(values: dict[int, float], ms: MidSection) -> float:
    """Mean of per-layer values over the mid-section, each layer weighted by `ms.weight` (missing layers skipped)."""
    pairs = [(ms.weight(k), v) for k, v in values.items() if ms.weight(k) > 0 and np.isfinite(v)]
    tot = sum(w for w, _ in pairs)
    return sum(w * v for w, v in pairs) / tot if tot > 0 else NAN


def _section(mask2d: np.ndarray) -> np.ndarray:
    return _fill_2d(np.asarray(mask2d, bool))


def section_area(mask2d: np.ndarray, px_um: tuple[float, float]) -> float:
    """Area (µm²) inside the outer outline(s) of a 2D mask (holes filled)."""
    return float(np.count_nonzero(_section(mask2d))) * px_um[0] * px_um[1]


def _contours(mask2d: np.ndarray, px_um: tuple[float, float], smooth_um: float) -> list[np.ndarray]:
    """Closed outlines (µm coordinates, y/x) of the smoothed mask's 0.5 iso-line."""
    from skimage.measure import find_contours

    m = np.asarray(mask2d, bool)
    if not m.any():
        return []
    sy, sx = (smooth_um / p for p in px_um)
    pad = int(math.ceil(4 * max(sy, sx))) + 2
    sl = ndi.find_objects(m.astype(np.uint8))[0]
    crop = np.pad(m[sl], pad).astype(np.float32)
    if smooth_um > 0:
        crop = ndi.gaussian_filter(crop, (sy, sx), mode="constant")
    return [c * np.asarray(px_um) for c in find_contours(crop, 0.5)]


def outline_length(mask2d: np.ndarray, px_um: tuple[float, float], smooth_um: float = 2.0) -> float:
    """Total outline length (µm) of a 2D mask at one smoothing scale (outer outlines only: holes filled)."""
    return float(sum(np.hypot(*np.diff(c, axis=0).T).sum() for c in _contours(_section(mask2d), px_um, smooth_um)))


def _polygon_area(c: np.ndarray) -> float:
    y, x = c[:, 0], c[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))))


def _fill_2d(mask2d: np.ndarray) -> np.ndarray:
    """Holes filled (background pieces, 4-connected, not touching the border) via labelling, fast on big layers."""
    bg, n = ndi.label(~mask2d)
    if n == 0:
        return mask2d.copy()
    border = np.zeros(n + 1, bool)
    for edge in (bg[0], bg[-1], bg[:, 0], bg[:, -1]):
        border[np.unique(edge)] = True
    border[0] = True
    return mask2d | ~border[bg]


def section_shape(section: np.ndarray, px_um: tuple[float, float], opts: ShapeOptions) -> dict:
    """Readouts of one filled 2D section: area, perimeter and the 2D shape descriptors."""
    from scipy.spatial import ConvexHull

    area = float(np.count_nonzero(section)) * px_um[0] * px_um[1]
    fine = _contours(section, px_um, opts.outline_smooth_um)
    perim = float(sum(np.hypot(*np.diff(c, axis=0).T).sum() for c in fine))
    def length(scale: float) -> float:
        return float(sum(np.hypot(*np.diff(c, axis=0).T).sum() for c in _contours(section, px_um, scale)))

    rough_fine, coarse = length(opts.roughness_fine_um), length(opts.roughness_smooth_um)
    pts = np.concatenate(fine) if fine else np.zeros((0, 2))
    hull = ConvexHull(pts).volume if len(pts) >= 3 else NAN
    yy, xx = np.nonzero(section)
    cov = np.cov(np.vstack([yy * px_um[0], xx * px_um[1]])) if len(yy) > 2 else np.full((2, 2), NAN)
    ev = np.sort(np.linalg.eigvalsh(cov)) if np.isfinite(cov).all() else (NAN, NAN)
    return {
        "area": area, "perimeter": perim,
        "elongation": math.sqrt(ev[1] / ev[0]) if ev[0] > 0 else NAN,
        "solidity": sum(_polygon_area(c) for c in fine) / hull if hull and hull > 0 else NAN,
        "edge_roughness": rough_fine / coarse if coarse > 0 else NAN,
        "circularity": 4 * math.pi * area / perim ** 2 if perim > 0 else NAN,
    }


def readouts(labels: np.ndarray, voxel_um: tuple[float, float, float], opts: ShapeOptions | None = None,
             prefix: str = "") -> dict:
    """Mid-section, confirmatory (cross-sectional area, perimeter) and exploratory readouts of a label stack
    (1 = main mass, 2.. = satellites; for a plain binary main-mass mask pass it as uint8)."""
    opts = opts or ShapeOptions()
    vz, vy, vx = voxel_um
    px = (vy, vx)
    main = labels == 1
    nz = main.shape[0]
    present = np.flatnonzero(main.any(axis=(1, 2)))
    if not len(present):
        return {f"{prefix}{k}": v for k, v in _empty().items()}
    # one crop around the main mass (all layers) keeps 2D work small on 2000 x 2000 layers
    sl = ndi.find_objects(main.any(axis=0).astype(np.uint8))[0]
    crop = main[:, sl[0], sl[1]]
    sections = {int(k): _fill_2d(crop[k]) for k in present}
    areas = np.zeros(nz)
    for k, s in sections.items():
        areas[k] = np.count_nonzero(s) * vy * vx
    ms = mid_section(areas, vz, opts.mid_section_um)
    per = {k: section_shape(sections[k], px, opts) for k in ms.layers if k in sections}
    mean = {}
    for key in next(iter(per.values())):  # layers of the mid-section without main mass count as area / perimeter 0
        vals = {k: p[key] for k, p in per.items()}
        if key in ("area", "perimeter"):
            vals.update({k: 0.0 for k in ms.layers if k not in per})
        mean[key] = weighted_mean(vals, ms)
    s0, s1, sw, slab_fits = layer_weights(nz, ms.centre, vz, opts.slab_um)
    slab_vox = float(np.dot(sw, crop[s0:s1 + 1].sum(axis=(1, 2))))
    main_vox = int(np.count_nonzero(crop))
    burden_vox = int(np.count_nonzero(labels))
    out = {
        "mid_section_centre_layer": ms.centre + 1, "mid_section_first_layer": ms.first + 1,
        "mid_section_last_layer": ms.last + 1, "mid_section_layers": len(ms.layers),
        "mid_section_um": round(ms.thickness_um, 4), "mid_section_fits": "yes" if ms.fits else "no",
        "cross_sectional_area_um2": round(mean["area"], 3), "perimeter_um": round(mean["perimeter"], 3),
        "footprint_area_um2": round(float(np.count_nonzero(_fill_2d(crop.any(axis=0)))) * vy * vx, 3),
        "slab_volume_um3": round(slab_vox * vz * vy * vx, 3),
        "slab_fits": "yes" if slab_fits else "no",
        "tumour_burden_um3": round(burden_vox * vz * vy * vx, 3),
        "invasion_fraction": (burden_vox - main_vox) / burden_vox,
        "elongation": mean["elongation"], "solidity": mean["solidity"],
        "edge_roughness": mean["edge_roughness"], "circularity": mean["circularity"],
    }
    return {f"{prefix}{k}": v for k, v in out.items()}


def label_readouts(labels: np.ndarray, voxel_um: tuple[float, float, float], opts: ShapeOptions | None = None
                   ) -> dict:
    """Every readout that depends only on the label image and the readout scales (`readouts` + the exploratory mesh
    surface): what the measurement cache re-derives when only the scales or definitions change."""
    opts = opts or ShapeOptions()
    return {**readouts(labels, voxel_um, opts),
            "main_mass_surface_um2": round(surface_area(labels == 1, voxel_um, opts.surface_grid_um,
                                                        opts.outline_smooth_um), 3)}


def _empty() -> dict:
    return {"mid_section_centre_layer": NAN, "mid_section_first_layer": NAN, "mid_section_last_layer": NAN,
            "mid_section_layers": 0, "mid_section_um": NAN, "mid_section_fits": "", "cross_sectional_area_um2": NAN, "perimeter_um": NAN,
            "footprint_area_um2": NAN, "slab_volume_um3": NAN, "slab_fits": "", "tumour_burden_um3": NAN,
            "invasion_fraction": NAN, "elongation": NAN, "solidity": NAN, "edge_roughness": NAN,
            "circularity": NAN}


def surface_area(main: np.ndarray, voxel_um: tuple[float, float, float], grid_um: float = 2.0,
                 smooth_um: float = 2.0) -> float:
    """Exploratory 3D mesh surface (µm²) of a binary main mass: resampled to an isotropic `grid_um` grid (block
    averages / linear interpolation of the occupancy), Gaussian-smoothed at `smooth_um`, marching cubes at 0.5.
    Includes the faces where the tumour is cut by the first/last imaged layer (it extends beyond the stack)."""
    from skimage.measure import marching_cubes, mesh_surface_area

    if not main.any():
        return NAN
    sl = ndi.find_objects(main.astype(np.uint8))[0]
    m = main[sl].astype(np.float32)
    zoom = [v / grid_um for v in voxel_um]
    iso = ndi.zoom(m, zoom, order=1, prefilter=False, grid_mode=True, mode="grid-constant")
    iso = np.pad(iso, 3)
    if smooth_um > 0:
        iso = ndi.gaussian_filter(iso, smooth_um / grid_um, mode="constant")
    if iso.max() <= 0.5:
        return NAN
    verts, faces, _, _ = marching_cubes(iso, 0.5, spacing=(grid_um,) * 3)
    return float(mesh_surface_area(verts, faces))
