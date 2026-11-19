"""Read microscope TIFFs of any common layout into a (channels, Y, X) float32 array."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile


@dataclass
class Image:
    data: np.ndarray  # (C, Y, X) float32
    max_value: float  # saturation level of the source dtype (or observed max for float)
    pixel_size_um: float | None


def _pixel_size(tf: tifffile.TiffFile) -> float | None:
    try:
        ij = tf.imagej_metadata or {}
        page = tf.pages[0]
        xres = page.tags.get("XResolution")
        if xres is None:
            return None
        num, den = xres.value
        if not num:
            return None
        per_unit = num / den  # pixels per unit
        unit = (ij.get("unit") or "").lower()
        res_unit = page.tags.get("ResolutionUnit")
        if unit in ("micron", "um", "µm", "\\u00b5m"):
            return 1.0 / per_unit
        if res_unit is not None and res_unit.value == 3:  # centimetre
            return 1e4 / per_unit
        return None
    except Exception:
        return None


def read(path: Path) -> Image:
    with tifffile.TiffFile(path) as tf:
        series = tf.series[0]
        arr = series.asarray()
        axes = series.axes.upper()
        px = _pixel_size(tf)
    dtype = arr.dtype
    max_value = float(np.iinfo(dtype).max) if np.issubdtype(dtype, np.integer) else float(arr.max())

    # Collapse to (C, Y, X): Z/T are max-projected, RGB samples become channels.
    arr = np.asarray(arr)
    keep = []
    for i, ax in enumerate(axes):
        if ax in ("Y", "X"):
            keep.append(ax)
        elif ax in ("C", "S"):
            keep.append("C")
        else:
            keep.append(None)
    for i in reversed(range(len(axes))):
        if keep[i] is None:
            arr = arr.max(axis=i)
    axes2 = [a for a in keep if a is not None]
    if axes2.count("C") > 1:  # e.g. C and S both present: merge
        first = axes2.index("C")
        second = axes2.index("C", first + 1)
        arr = np.moveaxis(arr, second, first + 1)
        shp = list(arr.shape)
        arr = arr.reshape(shp[:first] + [shp[first] * shp[first + 1]] + shp[first + 2 :])
        axes2.pop(second)
    if "C" not in axes2:
        arr = arr[None]
    else:
        arr = np.moveaxis(arr, axes2.index("C"), 0)
    arr = arr.astype(np.float32)

    # Greyscale saved as RGB: three identical channels -> one.
    if arr.shape[0] == 3 and np.array_equal(arr[0], arr[1]) and np.array_equal(arr[0], arr[2]):
        arr = arr[:1]
    return Image(arr, max_value, px)


def read_field(paths) -> Image:
    """One field from one file per channel (interleaved channel files): a (C, Y, X) stack in the order of `paths`.
    Each file must hold one channel; greyscale stored as RGB (three identical channels) is collapsed to one.
    Raises ValueError when a file has several different channels or the image sizes differ."""
    paths = list(paths)
    planes, max_value, px = [], 0.0, None
    for p in paths:
        im = read(p)
        if im.data.shape[0] != 1:
            raise ValueError(f"{p}: holds {im.data.shape[0]} different channels, but [channels] files = "
                             "\"interleaved\" expects one channel per file")
        if planes and im.data.shape[1:] != planes[0].shape:
            raise ValueError(f"{p}: image size {im.data.shape[1:]} differs from {paths[0]} {planes[0].shape}")
        planes.append(im.data[0])
        max_value = max(max_value, im.max_value)
        px = px or im.pixel_size_um
    return Image(np.stack(planes), max_value, px)


def read_record(rec) -> Image:
    """The field of a discover.ImageRecord: one multi-channel file, or (interleaved) one file per channel."""
    return read_field(rec.paths) if getattr(rec, "paths", ()) else read(rec.path)
