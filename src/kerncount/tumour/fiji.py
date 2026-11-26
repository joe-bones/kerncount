"""The lab's Fiji macro (`volume calculator tumor.ijm`) as a legacy readout (ADR 0002): its own outputs, read robustly,
and a faithful Python re-implementation used to validate against Fiji's masks and as a fallback.

Macro, per original z-stack (8 bit):
    Subtract Background (rolling ball 50 px, per layer)  -> `subtract_background`
    Median radius 1 (per layer)                           -> 3 x 3 median (ImageJ radius 1 = 3 x 3 kernel)
    setAutoThreshold("Triangle dark stack")               -> `triangle` on the histogram of the whole stack
    Convert to Mask                                       -> value >= level, level re-computed PER LAYER
                                                             (verified on the lab's masks, see macro_mask)
    Fill Holes, Open (per layer)                          -> 2D fill (4-connected background), 3 x 3 open
    3D Objects Counter (min 200 000 voxels)               -> `objects_3d`: 26-connected objects, Fiji's
                                                             border-voxel surface

Everything here follows the ImageJ / 3D Objects Counter sources (BackgroundSubtracter, AutoThresholder,
Counter3D) rather than a textbook version of each step, because the point is to reproduce Fiji's numbers:
- rolling ball: ImageJ shrinks the layer 4x (minimum of 4 x 4 blocks) for radius 50, rolls a trimmed ball patch
  (radius 12.5, 19 x 19) after a 3 x 3 mean pre-smoothing, and interpolates back linearly. That is a grey-scale
  opening with a non-flat structuring element, done here with scipy (fast: the shrunken layer is ~500 x 500).
- 3D Objects Counter's surface is the sum of exposed voxel faces of every surface voxel (6-neighbourhood; faces at
  the stack border count), accumulated in **float32** in scan order: on a 10^6 µm² object this rounding makes
  Fiji's number ~1 % lower than the exact face area. `objects_3d` reproduces it (bit-exact on the real masks).
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

# 3D Objects Counter statistics columns, in Fiji's order (header text -> short name)
FIJI_HEADER = (
    ("Volume (micron^3)", "volume_um3"), ("Surface (micron^2)", "surface_um2"), ("Nb of obj. voxels", "voxels"),
    ("Nb of surf. voxels", "surface_voxels"), ("IntDen", "intden"), ("Mean", "mean"), ("StdDev", "sd"),
    ("Median", "median"), ("Min", "min"), ("Max", "max"), ("X", "x"), ("Y", "y"), ("Z", "z"),
    ("Mean dist. to surf. (micron)", "mean_dist_um"), ("SD dist. to surf. (micron)", "sd_dist_um"),
    ("Median dist. to surf. (micron)", "median_dist_um"), ("XM", "xm"), ("YM", "ym"), ("ZM", "zm"),
    ("BX", "bx"), ("BY", "by"), ("BZ", "bz"), ("B-width", "b_width"), ("B-height", "b_height"),
    ("B-depth", "b_depth"),
)
COLUMNS = tuple(short for _, short in FIJI_HEADER)
_BY_HEADER = {h.lower(): s for h, s in FIJI_HEADER}
MIN_VOXELS = 200_000  # the macro's 3D Objects Counter "min." size (voxels)
ROLLING_BALL_PX = 50


@dataclass
class FijiStats:
    """Objects parsed from one 3D Objects Counter table (one row per object, `COLUMNS` names) + parse problems."""
    objects: list[dict] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)


_NUM = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$|^[+-]?(Infinity|NaN)$")


def _fields(line: str) -> list[str]:
    return [f.strip() for f in line.split("\t" if "\t" in line else ",")]


def parse_stats(text: str, name: str = "statistics") -> FijiStats:
    """Parse a 3D Objects Counter table: with header (`Statistics for …f.csv`, `…f_3Dresults.csv`, optionally a
    leading row-number column) or header-less tab/comma-separated (`… .txt`, Fiji's column order assumed).

    Robust to what spreadsheet round-trips leave behind: trailing `;;`, CRLF, blank lines, a truncated header (rows
    with all 25 values are read in Fiji's order). Rows that are not a full set of numbers (e.g. sum rows with
    decimal commas `3847005,25;;611.277`) are not guessed at: they are listed in `problems`."""
    out = FijiStats()
    header: list[str] | None = None
    index_col = False
    for ln, raw in enumerate(text.splitlines(), 1):
        line = raw.strip().rstrip(";").strip()
        if not line:
            continue
        f = _fields(line)
        if header is None and not out.objects and not all(_NUM.match(x) for x in f if x):
            if any(x.lower() in _BY_HEADER for x in f):
                index_col = f[0] == "" or f[0] == " " or f[0] in ("#", "Label")
                header = [_BY_HEADER.get(x.lower(), x) for x in (f[1:] if index_col else f)]
                unknown = [x for x in header if x not in COLUMNS]
                if unknown:
                    out.problems.append(f"{name} line {ln}: unknown columns {unknown} (ignored)")
                continue
        vals = f[1:] if index_col and len(f) > 1 else f
        if not vals or not all(_NUM.match(x) for x in vals):
            out.problems.append(f"{name} line {ln}: not a row of numbers, skipped: {raw.strip()[:80]!r}")
            continue
        if header is not None and len(vals) == len(header):
            names = header
        elif len(vals) == len(COLUMNS):
            names = list(COLUMNS)  # header-less (.txt) or truncated header: Fiji's column order
        else:
            out.problems.append(f"{name} line {ln}: {len(vals)} values, expected "
                                f"{len(header) if header else len(COLUMNS)}; skipped")
            continue
        out.objects.append({k: float(v) for k, v in zip(names, vals) if k in COLUMNS})
    if not out.objects:
        out.problems.append(f"{name}: no object rows")
    return out


def read_stats(path: Path) -> FijiStats:
    raw = Path(path).read_bytes()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    return parse_stats(text, Path(path).name)


def readouts(objects: list[dict], min_voxels: int = MIN_VOXELS, n_layers: int | None = None) -> dict:
    """Legacy readouts from Fiji's object table: the largest object ("main mass" by Fiji), its volume, Fiji's
    border-voxel surface, its z-extent (layers) and the number of objects >= `min_voxels`. NaN without objects."""
    objs = [o for o in objects if o.get("voxels", math.inf) >= min_voxels and math.isfinite(o.get("volume_um3", math.nan))]
    nan = math.nan
    if not objs:
        return {"fiji_main_volume_um3": nan, "fiji_main_volume_mm3": nan, "fiji_main_surface_um2": nan,
                "fiji_main_surface_mm2": nan, "fiji_n_objects": 0, "fiji_total_volume_mm3": nan,
                "fiji_main_z_start": nan, "fiji_main_z_layers": nan, "fiji_main_full_depth": ""}
    main = max(objs, key=lambda o: o["volume_um3"])
    z0, depth = main.get("bz", nan), main.get("b_depth", nan)
    full = "" if n_layers is None or not math.isfinite(depth) else ("yes" if z0 == 1 and depth == n_layers else "no")
    return {
        "fiji_main_volume_um3": main["volume_um3"], "fiji_main_volume_mm3": main["volume_um3"] / 1e9,
        "fiji_main_surface_um2": main.get("surface_um2", nan),
        "fiji_main_surface_mm2": main.get("surface_um2", nan) / 1e6,
        "fiji_n_objects": len(objs), "fiji_total_volume_mm3": sum(o["volume_um3"] for o in objs) / 1e9,
        "fiji_main_z_start": z0, "fiji_main_z_layers": depth, "fiji_main_full_depth": full,
    }


# ---------------------------------------------------------------- re-implementation of the macro


def _ball(radius: float) -> tuple[np.ndarray, int]:
    """ImageJ RollingBall: (patch heights float32, shrink factor)."""
    if radius <= 10:
        shrink, trim = 1, 24
    elif radius <= 30:
        shrink, trim = 2, 24
    elif radius <= 100:
        shrink, trim = 4, 32
    else:
        shrink, trim = 8, 40
    r = max(radius / shrink, 1.0)
    xtrim = int(trim * r) // 100
    half = int(math.floor(r - xtrim + 0.5))  # Java Math.round
    y, x = np.mgrid[-half:half + 1, -half:half + 1]
    t = r * r - x * x - y * y
    return np.where(t > 0, np.sqrt(np.maximum(t, 0)), 0).astype(np.float32), shrink


def _mean3(a: np.ndarray) -> np.ndarray:
    """ImageJ filter3x3(MEAN): 3-point means along rows, then columns, edge pixels repeated (float32)."""
    third = np.float32(0.333333333)
    for ax in (1, 0):
        p = np.pad(a, [(1, 1) if i == ax else (0, 0) for i in range(2)], mode="edge")
        lo = p[:-2] if ax == 0 else p[:, :-2]
        mid = p[1:-1] if ax == 0 else p[:, 1:-1]
        hi = p[2:] if ax == 0 else p[:, 2:]
        a = ((lo + mid + hi) * third).astype(np.float32)
    return a


def _interp_arrays(length: int, small: int, shrink: int) -> tuple[np.ndarray, np.ndarray]:
    i = np.arange(length)
    idx = np.trunc((i - shrink // 2) / shrink).astype(int)  # Java int division truncates toward zero
    idx = np.minimum(idx, small - 2)
    w = (1 - ((i + 0.5) / shrink - (idx + 0.5))).astype(np.float32)
    return idx, w


def rolling_ball_background(layer: np.ndarray, radius: float = ROLLING_BALL_PX) -> np.ndarray:
    """ImageJ's rolling-ball background of one 2D layer (dark background, smoothing on), float32."""
    ball, shrink = _ball(radius)
    h, w = layer.shape
    fp = _mean3(layer.astype(np.float32))
    if shrink > 1:
        sh, sw = -(-h // shrink), -(-w // shrink)
        pad = np.full((sh * shrink, sw * shrink), np.inf, np.float32)
        pad[:h, :w] = fp
        small = pad.reshape(sh, shrink, sw, shrink).min(axis=(1, 3))
    else:
        small = fp
    # rollBall: ball centres also outside the image (footprint clipped to the image) -> erosion with +inf outside,
    # then dilation of the centre heights, evaluated on the image only
    r = ball.shape[0] // 2
    big = np.pad(small, r, constant_values=np.inf)
    z = ndi.grey_erosion(big, structure=ball, mode="constant", cval=np.inf)
    bg = ndi.grey_dilation(z, structure=ball, mode="constant", cval=-np.inf)[r:-r, r:-r].astype(np.float32)
    if shrink == 1:
        return bg
    sh, sw = bg.shape
    yi, yw = _interp_arrays(h, sh, shrink)
    xi, xw = _interp_arrays(w, sw, shrink)
    lines = bg[:, xi] * xw + bg[:, xi + 1] * (1 - xw)  # x-interpolation of every small line
    return (lines[yi] * yw[:, None] + lines[yi + 1] * (1 - yw[:, None])).astype(np.float32)


def subtract_background(stack: np.ndarray, radius: float = ROLLING_BALL_PX) -> np.ndarray:
    """Process > Subtract Background (rolling ball, per layer) on an 8-bit stack: round(value - background)."""
    out = np.empty(stack.shape, np.uint8)
    for k, layer in enumerate(stack):
        v = layer.astype(np.float32) - rolling_ball_background(layer, radius) + np.float32(0.5)
        out[k] = np.clip(v, 0, 255).astype(np.uint8)
    return out


def triangle(hist: np.ndarray) -> int:
    """ImageJ AutoThresholder.Triangle (incl. its choice of the longer side of the histogram); returns the level;
    the macro's masks keep values >= level."""
    data = np.asarray(hist, np.int64).copy()
    nz = np.nonzero(data)[0]
    if not len(nz):
        return 0
    n = len(data)
    lo = int(nz[0]) - 1 if nz[0] > 0 else 0
    hi = int(nz[-1]) + 1 if nz[-1] < n - 1 else int(nz[-1])
    mx = int(np.argmax(data))
    inverted = (mx - lo) < (hi - mx)
    if inverted:
        data = data[::-1].copy()
        lo, mx = n - 1 - hi, n - 1 - mx
    if lo == mx:
        return lo
    nx, ny = float(data[mx]), float(lo - mx)
    d = math.sqrt(nx * nx + ny * ny)
    nx, ny = nx / d, ny / d
    d = nx * lo + ny * data[lo]
    split, best = lo, 0.0
    for i in range(lo + 1, mx + 1):
        dist = nx * i + ny * data[i] - d
        if dist > best:
            split, best = i, dist
    split -= 1
    return n - 1 - split if inverted else split


def macro_mask(stack: np.ndarray, radius: float = ROLLING_BALL_PX, per_layer_threshold: bool = True
               ) -> tuple[np.ndarray, int]:
    """(binary mask, stack Triangle level) of the macro up to the saved mask. 8-bit stacks only (as the macro).

    `per_layer_threshold=True` is what produced the lab's masks: Convert to Mask re-computed Triangle on every
    layer ("Calculate threshold for each image"), foreground = value >= the layer's level. With it the saved Fiji
    masks are reproduced voxel for voxel (docs/tumour.md). False = one Triangle level from the stack histogram,
    as the macro text alone suggests (setAutoThreshold "stack")."""
    if stack.dtype != np.uint8:
        raise ValueError(f"the Fiji macro re-implementation expects 8-bit stacks, got {stack.dtype}")
    work = np.empty(stack.shape, np.uint8)  # one layer at a time: peak memory ~ 2 x the 8-bit stack
    for k in range(stack.shape[0]):
        work[k] = ndi.median_filter(subtract_background(stack[k:k + 1], radius)[0], size=3, mode="nearest")
    level = triangle(np.bincount(work.ravel(), minlength=256))
    mask = np.empty(stack.shape, bool)
    square = np.ones((3, 3), bool)
    for k, layer in enumerate(work):
        t = triangle(np.bincount(layer.ravel(), minlength=256)) if per_layer_threshold else level
        m = ndi.binary_fill_holes(layer >= t)
        mask[k] = ndi.binary_opening(m, structure=square, border_value=0)
    return mask, level


def objects_3d(mask: np.ndarray, voxel_um: tuple[float, float, float], min_voxels: int = MIN_VOXELS
               ) -> tuple[np.ndarray, list[dict]]:
    """3D Objects Counter on a binary stack: (labels of the kept objects, largest first; one dict per object with
    volume_um3, surface_um2 (Fiji's float32 face sum), voxels, surface_voxels, bz (1-based), b_depth)."""
    vz, vy, vx = voxel_um
    lab, n = ndi.label(mask, structure=np.ones((3, 3, 3), bool))
    sizes = np.bincount(lab.ravel(), minlength=n + 1)
    sizes[0] = 0
    keep = [int(i) for i in np.argsort(sizes)[::-1] if sizes[i] >= max(min_voxels, 1)]
    lut = np.zeros(n + 1, np.uint16)
    lut[keep] = np.arange(1, len(keep) + 1)
    labels = lut[lab]
    del lab
    objs = []
    slices = ndi.find_objects(labels)
    for k in range(1, len(keep) + 1):
        sl = slices[k - 1]
        pad = tuple(slice(max(s.start - 1, 0), s.stop + 1) for s in sl)
        obj = labels[pad] == k
        # same-object 6-neighbours inside the stack (Counter3D: outside the stack is not a neighbour)
        p = np.pad(obj, 1)
        c = p[1:-1, 1:-1, 1:-1]
        nz = p[:-2, 1:-1, 1:-1].astype(np.int8) + p[2:, 1:-1, 1:-1]
        ny = p[1:-1, :-2, 1:-1].astype(np.int8) + p[1:-1, 2:, 1:-1]
        nx = p[1:-1, 1:-1, :-2].astype(np.int8) + p[1:-1, 1:-1, 2:]
        surf = c & ((nz + ny + nx) < 6)
        face = 2 * (vy * vz + vx * vz + vx * vy) - (nx * (vy * vz) + ny * (vx * vz) + nz * (vx * vy))
        area = float(np.cumsum(face[surf].astype(np.float32), dtype=np.float32)[-1]) if surf.any() else 0.0
        vox = int(c.sum())
        zs = np.nonzero(c.any(axis=(1, 2)))[0] + pad[0].start
        objs.append({"volume_um3": vox * vz * vy * vx, "surface_um2": area, "voxels": vox,
                     "surface_voxels": int(surf.sum()), "bz": int(zs[0]) + 1, "b_depth": int(zs[-1] - zs[0] + 1)})
    return labels, objs


def read_mask(path: Path) -> np.ndarray:
    """A Fiji binary mask (0 / 255) as bool ZYX; 3D Objects Counter's threshold=128."""
    import tifffile

    arr = tifffile.imread(path)
    if arr.ndim == 2:
        arr = arr[None]
    return arr >= 128


def dice(a: np.ndarray, b: np.ndarray) -> float:
    s = int(np.count_nonzero(a)) + int(np.count_nonzero(b))
    return 2 * int(np.count_nonzero(a & b)) / s if s else math.nan
