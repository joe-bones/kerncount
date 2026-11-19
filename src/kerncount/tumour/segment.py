"""Primary tumour segmentation (ADR 0002): one z-stack -> labelled objects -> main mass and its volume.

    smoothing (Gaussian, fixed physical scale)
    -> per-tumour automatic threshold WITHOUT rolling-ball subtraction ("intermedians", see below / docs/tumour.md)
    -> 3D hole filling -> connected components (26-connected) -> main mass = largest object

Volumes are in µm³ from the voxel size in each file's OME metadata; a file without it is an error, never pixels.

Memory: one stack at a time (uint8, up to ~300 MB). Peak is about stack + float32 smoothed copy (4x stack), then
mask + int32 labels (4x stack): < 2 GB for the largest stacks. Nothing full-resolution in float is kept.

Threshold "intermedians": iterative intermeans (Ridler-Calvard / IsoData) with class MEDIANS, started at Otsu:
t = (median of voxels <= t + median of voxels > t) / 2, repeated until stable. It is
- scale-equivariant: multiplying the recorded signal (detector gain, mCherry expression) multiplies t, so the
  segmentation does not change (only 8-bit rounding and clipping remain);
- the half-way level between tissue background and tumour signal, i.e. the half-maximum edge of a blurred
  boundary, whatever the share of the stack the tumour fills (Otsu weights classes by size and drifts toward the
  background for small tumours);
- insensitive to saturation while less than half of the tumour voxels clip (a median ignores the 255 spike),
  which keeps Fiji's Triangle failure (threshold jumps to the saturated peak) out.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

THRESHOLDS = ("intermedians", "otsu")
_UNIT_UM = {"µm": 1.0, "um": 1.0, "micron": 1.0, "microns": 1.0, "nm": 1e-3, "mm": 1e3, "cm": 1e4, "m": 1e6}
_HIST_BINS = 1024


class MissingCalibration(ValueError):
    """The stack carries no physical voxel size: volumes would silently be in voxels."""


@dataclass
class Segmentation:
    """labels: 0 = background, 1 = main mass, 2.. = satellites (largest first); readouts: one table row."""
    labels: np.ndarray
    readouts: dict


def voxel_size_um(path: Path) -> tuple[float, float, float]:
    """(z, y, x) voxel size in µm from the OME-XML of `path` (units converted; OME's default unit is µm)."""
    import tifffile

    with tifffile.TiffFile(path) as tf:
        xml = tf.ome_metadata
    if not xml:
        raise MissingCalibration(f"{Path(path).name}: no OME metadata, so no voxel size; "
                                 "volumes need the acquisition voxel size (export the original .ome.tiff)")
    root = ET.fromstring(xml)
    px = next((el for el in root.iter() if el.tag.rsplit("}", 1)[-1] == "Pixels"), None)
    out = []
    for ax in ("Z", "Y", "X"):
        v = px.get(f"PhysicalSize{ax}") if px is not None else None
        unit = (px.get(f"PhysicalSize{ax}Unit") if px is not None else None) or "µm"
        if v in (None, "") or unit not in _UNIT_UM or float(v) <= 0:
            raise MissingCalibration(f"{Path(path).name}: voxel size {ax} missing in the OME metadata "
                                     f"(PhysicalSize{ax}={v!r}, unit {unit!r}); refusing to measure in voxels")
        out.append(float(v) * _UNIT_UM[unit])
    return out[0], out[1], out[2]


def read_stack(path: Path) -> np.ndarray:
    """The z-stack as a ZYX array (single channel)."""
    import tifffile

    with tifffile.TiffFile(path) as tf:
        s = tf.series[0]
        arr = s.asarray()
        axes = s.axes
    keep = [i for i, a in enumerate(axes) if a in "ZYX" or arr.shape[i] > 1]
    if len(keep) != arr.ndim:
        arr = arr.squeeze(axis=tuple(i for i in range(arr.ndim) if i not in keep))
    if arr.ndim == 2:
        arr = arr[None]
    if arr.ndim != 3:
        raise ValueError(f"{Path(path).name}: expected one channel z-stack (ZYX), got axes {axes} {arr.shape}")
    return arr


def _wmedian(counts: np.ndarray, centres: np.ndarray) -> float:
    cs = np.cumsum(counts)
    return float(centres[min(int(np.searchsorted(cs, cs[-1] / 2)), len(centres) - 1)])


def threshold_from_histogram(counts: np.ndarray, centres: np.ndarray, method: str = "intermedians"
                             ) -> tuple[float, float, float]:
    """(threshold, background level, tumour level) from a histogram of the smoothed stack."""
    from skimage.filters import threshold_otsu

    if method not in THRESHOLDS:
        raise ValueError(f"[tumour] threshold = {method!r}; valid: {', '.join(THRESHOLDS)}")
    t = float(threshold_otsu(hist=(counts, centres)))
    for _ in range(100 if method == "intermedians" else 0):
        lo = centres <= t
        if not counts[lo].any() or not counts[~lo].any():
            break
        new = (_wmedian(counts[lo], centres[lo]) + _wmedian(counts[~lo], centres[~lo])) / 2
        if abs(new - t) < 1e-6:
            break
        t = new
    lo = centres <= t
    bg = _wmedian(counts[lo], centres[lo]) if counts[lo].any() else float("nan")
    fg = _wmedian(counts[~lo], centres[~lo]) if counts[~lo].any() else float("nan")
    return t, bg, fg


def _fill_holes_3d(mask: np.ndarray) -> None:
    """Fill enclosed background in place: background components (6-connected) not touching the stack border.
    Labelling instead of ndi.binary_fill_holes (iterative dilation, far too slow on 10^8 voxels)."""
    bg, n = ndi.label(~mask)
    if n == 0:
        return
    border = np.zeros(n + 1, bool)
    for i in (0, -1):
        border[np.unique(bg[i])] = True
        border[np.unique(bg[:, i])] = True
        border[np.unique(bg[:, :, i])] = True
    border[0] = True
    mask |= ~border[bg]


def segment(stack: np.ndarray, voxel_um: tuple[float, float, float], smooth_sigma_um: float = 2.0,
            threshold: str = "intermedians", min_object_um3: float = 1000.0) -> Segmentation:
    vz, vy, vx = voxel_um
    voxel_um3 = vz * vy * vx
    top = np.iinfo(stack.dtype).max if np.issubdtype(stack.dtype, np.integer) else float(stack.max())
    saturated = stack >= top
    n_sat = int(np.count_nonzero(saturated))
    del saturated

    sm = ndi.gaussian_filter(stack, [smooth_sigma_um / v for v in voxel_um], output=np.float32)
    counts, edges = np.histogram(sm, bins=_HIST_BINS, range=(0.0, float(top) + 1.0))
    centres = (edges[:-1] + edges[1:]) / 2
    t, bg_level, fg_level = threshold_from_histogram(counts, centres, threshold)
    mask = sm > np.float32(t)
    del sm
    _fill_holes_3d(mask)

    lab, n = ndi.label(mask, structure=np.ones((3, 3, 3), bool))
    del mask
    sizes = np.bincount(lab.ravel(), minlength=n + 1)
    sizes[0] = 0
    min_vox = min_object_um3 / voxel_um3
    order = [int(i) for i in np.argsort(sizes)[::-1] if sizes[i] > 0 and sizes[i] >= min_vox]
    dtype = np.uint8 if len(order) < 255 else np.uint16
    order = order[: np.iinfo(dtype).max]
    lut = np.zeros(n + 1, dtype)
    lut[order] = np.arange(1, len(order) + 1, dtype=dtype)
    labels = lut[lab]
    del lab

    main_vox = int(sizes[order[0]]) if order else 0
    main = labels == 1
    sat_main = int(np.count_nonzero(main & (stack >= top))) if main_vox else 0
    del main
    readouts = {
        "voxel_z_um": vz, "voxel_y_um": vy, "voxel_x_um": vx,
        "n_layers": int(stack.shape[0]), "rows": int(stack.shape[1]), "columns": int(stack.shape[2]),
        "threshold_method": threshold, "threshold": round(t, 3),
        "background_level": round(bg_level, 3), "tumour_level": round(fg_level, 3),
        "smooth_sigma_um": smooth_sigma_um, "min_object_um3": min_object_um3,
        "main_mass_voxels": main_vox, "main_mass_volume_um3": round(main_vox * voxel_um3, 3),
        "n_satellites": max(len(order) - 1, 0),
        "saturated_frac": n_sat / stack.size,
        "saturated_frac_main_mass": sat_main / main_vox if main_vox else float("nan"),
    }
    return Segmentation(labels, readouts)


def measure_file(path: Path, tumour_settings: dict) -> Segmentation:
    """Read one original stack, check its calibration, segment it with the [tumour] settings."""
    vox = voxel_size_um(path)  # before reading the pixels: an uncalibrated file fails fast
    stack = read_stack(path)
    return segment(stack, vox, smooth_sigma_um=float(tumour_settings["smooth_sigma_um"]),
                   threshold=str(tumour_settings["threshold"]),
                   min_object_um3=float(tumour_settings["min_object_um3"]))
