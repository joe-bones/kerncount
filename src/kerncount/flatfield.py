"""Retrospective flat-field (illumination) correction, BaSiC-style, without extra dependencies.

Model (as in BaSiC, Peng et al. 2017):   observed(x) = signal(x) * F(x) + D(x)
    F  smooth multiplicative flat field (illumination of the specific signal), mean 1 over the image
    D  smooth additive dark field (camera offset + autofluorescence / out-of-cell background)

Estimate from all images of one channel in one n (a plate: same lamp, same day). Every image is reduced to blocks
(default 16 x 16 px):
    block low  = 20th percentile inside the block           -> background level there
    block mean = mean inside the block                       -> background + coverage x signal
D  = median over images of block low, fitted by a smooth surface (or a constant when `dark=False`)
F  = trimmed mean over images of (block mean - D) / (image mean of that), fitted by a smooth surface, mean 1.
Cells are placed at random relative to the optics, so averaging many fields leaves only the illumination pattern.
The smooth surface is a 2-D polynomial (degree 2 by default: higher degrees chase noise into the corners) fitted with robust (Huber) reweighting, which
extrapolates sensibly into the corners where a Gaussian blur would flatten the fall-off.

`F` is applied to background-subtracted intensities only (value / F at the nucleus centre). Counts and detection
never use it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

BLOCK = 16
DEGREE = 2
MIN_IMAGES = 10  # an illumination pattern from fewer fields mostly reflects where the cells happen to be
ITERS_FLAT = 1  # plain least squares for F: the block averages are skewed (sparse nuclei), Huber weights would bias it


@dataclass
class FlatField:
    shape: tuple[int, int]
    flat_coef: np.ndarray  # polynomial coefficients of F (normalised so that mean F over the image = 1)
    dark_coef: np.ndarray | None  # polynomial coefficients of D (None = constant `dark_const`)
    dark_const: float = 0.0
    degree: int = DEGREE
    source: str = "estimated"  # "estimated" | "reference" | "off"
    n_images: int = 0
    info: dict = field(default_factory=dict)

    # ------------------------------------------------------------ evaluation
    def flat(self, y=None, x=None) -> np.ndarray:
        """F at points (y, x) in pixels, or the full image when no points are given."""
        if self.source == "off":
            return np.ones(self.shape, np.float32) if y is None else np.ones(np.shape(y), float)
        return _poly_eval(self.flat_coef, self.degree, self.shape, y, x)

    def dark(self, y=None, x=None) -> np.ndarray:
        if self.dark_coef is None:
            return np.full(self.shape if y is None else np.shape(y), self.dark_const, np.float32)
        return _poly_eval(self.dark_coef, self.degree, self.shape, y, x)

    def correct(self, img: np.ndarray) -> np.ndarray:
        """(img - D) / F: an illumination-corrected image (for display; readouts correct per nucleus)."""
        return (img.astype(np.float32) - self.dark()) / self.flat()

    def to_dict(self) -> dict:
        return dict(shape=list(self.shape), flat_coef=self.flat_coef.tolist(),
                    dark_coef=None if self.dark_coef is None else self.dark_coef.tolist(),
                    dark_const=self.dark_const, degree=self.degree, source=self.source, n_images=self.n_images,
                    info=self.info)

    @classmethod
    def from_dict(cls, d: dict) -> "FlatField":
        return cls(tuple(d["shape"]), np.asarray(d["flat_coef"], float),
                   None if d.get("dark_coef") is None else np.asarray(d["dark_coef"], float),
                   float(d.get("dark_const", 0.0)), int(d.get("degree", DEGREE)), d.get("source", "estimated"),
                   int(d.get("n_images", 0)), dict(d.get("info", {})))

    @classmethod
    def identity(cls, shape) -> "FlatField":
        return cls(tuple(shape), np.zeros(1), None, 0.0, 0, "off", 0, {})


# ---------------------------------------------------------------- polynomial surface

def _terms(degree: int):
    return [(i, j) for i in range(degree + 1) for j in range(degree + 1 - i)]


def _design(yn: np.ndarray, xn: np.ndarray, degree: int) -> np.ndarray:
    return np.stack([yn ** i * xn ** j for i, j in _terms(degree)], axis=-1)


def _norm_coords(y, x, shape):
    H, W = shape
    return 2 * np.asarray(y, float) / max(H - 1, 1) - 1, 2 * np.asarray(x, float) / max(W - 1, 1) - 1


def _poly_eval(coef, degree, shape, y=None, x=None) -> np.ndarray:
    if y is None:  # separable evaluation on the full grid: Py @ C @ Px.T
        H, W = shape
        yn, xn = _norm_coords(np.arange(H), np.arange(W), shape)
        C = np.zeros((degree + 1, degree + 1))
        for c, (i, j) in zip(coef, _terms(degree)):
            C[i, j] = c
        Py = yn[:, None] ** np.arange(degree + 1)[None]
        Px = xn[:, None] ** np.arange(degree + 1)[None]
        return (Py @ C @ Px.T).astype(np.float32)
    yn, xn = _norm_coords(y, x, shape)
    return _design(np.ravel(yn), np.ravel(xn), degree) @ coef if np.ndim(y) else float(_design(yn, xn, degree) @ coef)


def fit_surface(values: np.ndarray, ys: np.ndarray, xs: np.ndarray, shape, degree: int = DEGREE,
                weights: np.ndarray | None = None, iters: int = 6) -> np.ndarray:
    """Robust (Huber IRLS) least-squares 2-D polynomial; returns coefficients for `_poly_eval`."""
    ok = np.isfinite(values)
    v, y, x = values[ok], ys[ok], xs[ok]
    w0 = np.ones_like(v) if weights is None else np.asarray(weights, float)[ok]
    A = _design(*_norm_coords(y, x, shape), degree)
    w = w0.copy()
    coef = np.zeros(A.shape[1])
    for _ in range(iters):
        sw = np.sqrt(w)
        coef, *_ = np.linalg.lstsq(A * sw[:, None], v * sw, rcond=None)
        r = v - A @ coef
        s = 1.4826 * np.median(np.abs(r - np.median(r))) or 1e-12
        u = np.abs(r) / (1.345 * s)
        w = w0 * np.where(u <= 1, 1.0, 1.0 / np.maximum(u, 1e-12))
    return coef


# ---------------------------------------------------------------- estimation

def _blocks(img: np.ndarray, b: int) -> tuple[np.ndarray, np.ndarray]:
    H, W = img.shape
    h, w = H // b, W // b
    a = img[: h * b, : w * b].reshape(h, b, w, b).swapaxes(1, 2).reshape(h, w, b * b)
    return a.mean(axis=2), np.percentile(a, 20, axis=2)


def estimate(images, dark: bool = True, block: int = BLOCK, degree: int = DEGREE, min_images: int = MIN_IMAGES) -> FlatField:
    """Estimate F (and D) from an iterable of 2-D images of one channel (same size). Needs >= MIN_IMAGES images (fewer: no correction, source "off")."""
    means, lows, shape = [], [], None
    for img in images:
        img = np.asarray(img, np.float32)
        if shape is None:
            shape = img.shape
        elif img.shape != shape:
            continue
        m, lo = _blocks(img, block)
        means.append(m)
        lows.append(lo)
    if shape is None or len(means) < min_images:
        return FlatField.identity(shape or (1, 1))
    M, L = np.stack(means), np.stack(lows)
    h, w = M.shape[1:]
    yc = (np.arange(h) + 0.5) * block - 0.5
    xc = (np.arange(w) + 0.5) * block - 0.5
    YY, XX = np.meshgrid(yc, xc, indexing="ij")
    Dmap = np.median(L, axis=0)
    if dark:
        dark_coef = fit_surface(Dmap.ravel(), YY.ravel(), XX.ravel(), shape, degree)
        Dfit = (_design(*_norm_coords(YY.ravel(), XX.ravel(), shape), degree) @ dark_coef).reshape(h, w)
        dark_const = float(np.median(Dmap))
    else:
        dark_coef, dark_const = None, float(np.median(Dmap))
        Dfit = np.full((h, w), dark_const)
    S = M - Dfit[None]
    scale = S.reshape(len(S), -1).mean(axis=1)
    keep = scale > 0
    if keep.sum() < min_images:
        return FlatField.identity(shape)
    R = S[keep] / scale[keep, None, None]


    Fraw = R.mean(axis=0)
    coef = fit_surface(Fraw.ravel(), YY.ravel(), XX.ravel(), shape, degree, iters=ITERS_FLAT)
    ff = FlatField(tuple(shape), coef, dark_coef, dark_const, degree, "estimated", int(keep.sum()))
    mean_f = float(ff.flat()[::4, ::4].mean())
    ff.flat_coef = coef / mean_f
    resid = Fraw / mean_f - ff.flat(YY.ravel(), XX.ravel()).reshape(h, w)
    ff.info = dict(block=block, fit_residual_sd=float(np.std(resid)), **magnitude(ff))
    return ff


def from_reference(files: list[Path], shape, dark_files: list[Path] | None = None, degree: int = DEGREE) -> FlatField:
    """Flat field from lab reference images (e.g. a uniformly fluorescent slide), optionally with dark frames."""
    from . import imageio

    def avg(paths):
        acc = None
        for p in paths:
            a = imageio.read(p).data[0].astype(np.float64)
            acc = a if acc is None else acc + a
        return acc / len(paths)

    ref = avg(files)
    if ref.shape != tuple(shape):
        raise ValueError(f"reference images are {ref.shape}, the data {tuple(shape)}")
    d = avg(dark_files) if dark_files else np.zeros_like(ref)
    m, _ = _blocks(ref - d, BLOCK)
    h, w = m.shape
    YY, XX = np.meshgrid((np.arange(h) + 0.5) * BLOCK - 0.5, (np.arange(w) + 0.5) * BLOCK - 0.5, indexing="ij")
    coef = fit_surface(m.ravel(), YY.ravel(), XX.ravel(), shape, degree)
    ff = FlatField(tuple(shape), coef, None, float(np.median(d)), degree, "reference", len(files))
    ff.flat_coef = coef / float(ff.flat()[::4, ::4].mean())
    ff.info = dict(reference_files=[str(p) for p in files], **magnitude(ff))
    return ff


def reference_files(folder: Path, channel: str) -> tuple[list[Path], list[Path]]:
    """Reference and dark images in `folder`: files whose name contains the channel name (else all TIFFs);
    files with "dark" in the name are dark frames."""
    tifs = sorted(p for p in Path(folder).rglob("*") if p.suffix.lower() in (".tif", ".tiff") and p.is_file())
    darks = [p for p in tifs if "dark" in p.stem.lower()]
    refs = [p for p in tifs if p not in darks]
    named = [p for p in refs if channel.lower() in p.stem.lower()]
    return (named or refs), darks


def magnitude(ff: FlatField) -> dict:
    """How uneven the illumination is: mean F in the centre (middle 20% x 20%) vs the outer 10% band, per side,
    and the range. edge_dim_pct = how much dimmer the edges are than the centre, in %."""
    F = ff.flat()
    H, W = F.shape
    cy, cx = slice(int(0.4 * H), int(0.6 * H) + 1), slice(int(0.4 * W), int(0.6 * W) + 1)
    centre = float(F[cy, cx].mean())
    by, bx = max(int(0.1 * H), 1), max(int(0.1 * W), 1)
    band = np.zeros_like(F, bool)
    band[:by] = band[-by:] = True
    band[:, :bx] = band[:, -bx:] = True
    sides = dict(top=F[:by].mean(), bottom=F[-by:].mean(), left=F[:, :bx].mean(), right=F[:, -bx:].mean())
    darkest = min(sides, key=sides.get)
    return dict(centre=round(centre, 4), edge_dim_pct=round(100 * (1 - float(F[band].mean()) / centre), 1),
                darkest_side=darkest, darkest_side_dim_pct=round(100 * (1 - float(sides[darkest]) / centre), 1),
                min=round(float(F.min()), 4), max=round(float(F.max()), 4),
                range_pct=round(100 * (float(F.max()) - float(F.min())) / centre, 1))


def describe(mag: dict) -> str:
    if not mag:
        return "no flat-field correction"
    return (f"edges {mag['edge_dim_pct']:.0f}% dimmer than the centre (darkest: {mag['darkest_side']} edge, "
            f"{mag['darkest_side_dim_pct']:.0f}%); range {mag['range_pct']:.0f}% of the centre")


def object_profile(values: np.ndarray, ys: np.ndarray, xs: np.ndarray, shape, degree: int = 2) -> dict:
    """Cross-check from the nuclei themselves: robust fit of log(background-subtracted intensity) on position.
    Cells sit at random positions, so a position trend is illumination. Returns magnitude() of that surface."""
    ok = np.isfinite(values) & (values > 0)
    if ok.sum() < 200:
        return {}
    coef = fit_surface(np.log(values[ok]), ys[ok], xs[ok], shape, degree)
    ff = FlatField(tuple(shape), coef, None, 0.0, degree)
    F = np.exp(ff.flat())
    ff2 = FlatField(tuple(shape), np.zeros(1), None, 0.0, 0)
    ff2.flat = lambda y=None, x=None, _F=F / F.mean(): _F  # type: ignore[method-assign]
    return magnitude(ff2)
