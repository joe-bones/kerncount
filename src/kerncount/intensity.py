"""Per-nucleus marker intensity (e.g. Ki-67): measurement, local background, noise, channel alignment.

ADR docs/adr/0001 and docs/intensity.md. Everything here is measured on the raw marker channel with the nuclei
found on the nuclear channel; the flat field (flatfield.py) is applied afterwards per nucleus, and positivity
(positivity.py) on top. Detection and counts are never touched.

Per nucleus (one row per label and marker channel):
    area             nucleus area (px)
    erode_px         how far the outline was shrunk before measuring (rule in `erode_amount`)
    measured_px      pixels left after shrinking (N)
    mean, p90        mean and 90th percentile of the raw marker inside the shrunk outline
    total            integrated intensity over the whole (unshrunk) nucleus minus area x local background
    bg               local background: median of a ring outside the nucleus that avoids every nucleus
    bg_px, bg_source ring pixels (M); "ring", "wide ring" or "image" (fallback when the ring is too small)
    fiji_cell_mean   legacy: raw mean of the nucleus enlarged by fiji_ring_px (old Fiji macro, Enlarge 6 px)
    saturated        any pixel of the nucleus at the camera maximum in this channel (excluded from intensity
                     readouts, still counted)
    noise_sd         SD of (mean - bg) expected from background noise alone for this nucleus (raw units):
                     sigma * sqrt(a_c / N + (pi/2) a_c / M); sigma = robust per-pixel noise SD, a_c = noise
                     correlation area in px (noise is not independent between neighbouring pixels)
    noise_sd_p90     same for (p90 - bg): sigma * sqrt(2.92 a_c / N + (pi/2) a_c / M)
    p90_null         expected (p90 - bg) of a nucleus with no signal: 1.2816 sigma
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import segmentation

from . import settings as settings_mod

MEASURE_VERSION = 2
MISPAIR_CORR = 0.2  # marker vs nuclear image correlation below this: the files do not show the same field
Z90 = 1.2816  # standard-normal 90th percentile
P90_VAR = 0.09 / 0.17550 ** 2  # var of a sample 90th percentile x N / sigma^2 (p(1-p) / phi(z90)^2) = 2.92
MEDIAN_VAR = np.pi / 2  # var of a sample median x M / sigma^2


# ---------------------------------------------------------------- settings helpers

def marker_channels(settings: dict, n_channels: int) -> list[int]:
    """0-based channels measured per nucleus ([intensity] per_nucleus); never the nuclear channel under "auto"."""
    it = settings["intensity"]
    per = it.get("per_nucleus", "auto")
    if not it.get("enabled", True) or per in ("off", "", None, False) or per == []:
        return []
    seg = settings_mod.channel_index(settings, settings["channels"]["segment"], n_channels)
    if str(per).lower() == "auto":
        return [c for c in range(n_channels) if c != seg]
    items = per if isinstance(per, list) else [x.strip() for x in str(per).replace(";", ",").split(",") if x.strip()]
    return sorted({settings_mod.channel_index(settings, x, n_channels) for x in items})


def _num(v, default):
    return default if v in ("", None) or str(v).lower() == "auto" else float(v)


def erode_amount(area: np.ndarray, setting="auto") -> np.ndarray:
    """Pixels to shrink each nucleus outline by. "auto": 2 px for nuclei at least 16 px across (equivalent circle
    diameter), else 1 px; a number = that for every nucleus. In `measure_nuclei` the amount is lowered again when
    shrinking would leave fewer than max(9 px, 25% of the area)."""
    area = np.asarray(area, float)
    if str(setting).lower() != "auto":
        return np.full(area.shape, int(float(setting)), int)
    d = np.sqrt(4 * area / np.pi)
    return np.where(d >= 16, 2, 1)


def ring_width(settings: dict, diameter: float) -> int:
    v = settings["intensity"].get("bg_ring_px", "auto")
    return max(int(round(0.5 * diameter if str(v).lower() == "auto" else float(v))), 2)


# ---------------------------------------------------------------- per-label statistics

def _label_sorted(values: np.ndarray, labels: np.ndarray, n: int):
    """Values grouped by label and sorted within each label: (sorted values, start index, count) per label 0..n."""
    m = labels > 0
    lv, vv = labels[m].ravel(), values[m].ravel()
    order = np.lexsort((vv, lv))
    lv, vv = lv[order], vv[order]
    counts = np.bincount(lv, minlength=n + 1)
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    return vv, starts, counts


def _quantile_per_label(values, labels, n, q: float) -> np.ndarray:
    vv, starts, counts = _label_sorted(values, labels, n)
    out = np.full(n + 1, np.nan)
    ok = counts > 0
    pos = q * (counts[ok] - 1)
    lo = np.floor(pos).astype(int)
    hi = np.minimum(lo + 1, counts[ok] - 1)
    frac = pos - lo
    s = starts[ok]
    out[ok] = vv[s + lo] * (1 - frac) + vv[s + hi] * frac
    return out


def background_noise(img: np.ndarray, bgmask: np.ndarray, diameter: float) -> tuple[float, float]:
    """(sigma, a_c): robust per-pixel noise SD of the background after removing its smooth trend, and the noise
    correlation area in pixels (1 = independent pixels), from 4 x 4 block means of background-only blocks."""
    if bgmask.sum() < 500:
        return float("nan"), 1.0
    s = max(diameter, 4.0)
    w = ndi.gaussian_filter(bgmask.astype(np.float32), s)
    trend = ndi.gaussian_filter(np.where(bgmask, img, 0).astype(np.float32), s) / np.maximum(w, 1e-3)
    r = img - trend
    rb = r[bgmask]
    sigma = float(1.4826 * np.median(np.abs(rb - np.median(rb))))
    H, W = img.shape
    est = {}
    for b in (4, 8):
        h, wd = H // b, W // b
        blk_r = r[: h * b, : wd * b].reshape(h, b, wd, b).mean(axis=(1, 3))
        blk_m = bgmask[: h * b, : wd * b].reshape(h, b, wd, b).all(axis=(1, 3))
        if blk_m.sum() >= 50 and sigma > 0:
            x = blk_r[blk_m]
            sb = float(1.4826 * np.median(np.abs(x - np.median(x))))
            est[b] = float(np.clip(b * b * sb ** 2 / sigma ** 2, 1.0, b * b))
    if 4 in est and 8 in est:
        # a block of side b sees a_c(b) ~ A (1 - c / b) (correlation is cut at its edges): extrapolate to large areas
        a_c = max(est[8], 2 * est[8] - est[4])
    else:
        a_c = est.get(8, est.get(4, 1.0))
    return sigma, float(a_c)


def measure_nuclei(img: np.ndarray, labels: np.ndarray, max_value: float, diameter: float, *, erode="auto",
                   gap: int = 3, ring: int | None = None, bg_min: int = 40, fiji_ring: int = 6) -> tuple[pd.DataFrame, dict]:
    """Per-nucleus measurements of one marker image (see the module docstring). `labels` are all detected objects
    (nuclei, debris, ...): every object is kept out of every background ring."""
    img = np.asarray(img, np.float32)
    labels = np.asarray(labels, np.int32)
    n = int(labels.max())
    info: dict = {}
    obj = labels > 0
    dist_out = ndi.distance_transform_edt(~obj)
    far = dist_out > gap  # background pixels: more than `gap` px away from every object
    img_bg = float(np.median(img[far])) if far.sum() > 100 else float(np.median(img[~obj])) if (~obj).any() else 0.0
    sigma, a_c = background_noise(img, far, diameter)
    info.update(image_bg=img_bg, noise_sigma=sigma, noise_corr_area=a_c, bg_pixels=int(far.sum()))
    if n == 0:
        return pd.DataFrame(), info
    idx = np.arange(1, n + 1)
    area = np.bincount(labels.ravel(), minlength=n + 1).astype(float)

    # shrunk outlines: distance from the nucleus rim (contacts with other nuclei count as rim)
    rim = segmentation.find_boundaries(labels, mode="inner") & obj
    depth = ndi.distance_transform_edt(~rim)  # 0 on the rim, 1 next to it, ...
    want = np.zeros(n + 1, int)
    want[1:] = erode_amount(area[1:], erode)
    e = want.copy()
    for k in range(int(want.max()), 0, -1):  # lower the amount where too little would be left
        left = np.bincount(labels[obj & (depth >= e[labels])].ravel(), minlength=n + 1)
        bad = (e >= k) & (left < np.maximum(9, 0.25 * area)) & (np.arange(n + 1) > 0)
        e[bad] = np.maximum(e[bad] - 1, 0)
    core = np.where(obj & (depth >= e[labels]), labels, 0)
    npx = np.bincount(core.ravel(), minlength=n + 1).astype(float)
    sums = np.bincount(core.ravel(), weights=img.ravel(), minlength=n + 1)
    mean = sums / np.maximum(npx, 1)
    p90 = _quantile_per_label(img, core, n, 0.9)

    # local background ring(s), Voronoi-split between neighbours, never closer than `gap` to any object
    ring = ring or max(int(round(0.5 * diameter)), 2)
    bg = np.full(n + 1, np.nan)
    bgpx = np.zeros(n + 1)
    src = np.array([""] * (n + 1), dtype=object)
    todo = np.ones(n + 1, bool)
    todo[0] = False
    for width, name in ((ring, "ring"), (3 * ring, "wide ring")):
        if not todo.any():
            break
        grown = segmentation.expand_labels(labels, gap + width)
        rl = np.where(far, grown, 0)
        rl = np.where(todo[rl], rl, 0)
        med = _quantile_per_label(img, rl, n, 0.5)
        cnt = np.bincount(rl.ravel(), minlength=n + 1)
        ok = todo & (cnt >= bg_min)
        bg[ok], bgpx[ok], src[ok] = med[ok], cnt[ok], name
        todo &= ~ok
    bg[todo], bgpx[todo], src[todo] = img_bg, float(far.sum()), "image"

    full_sum = np.bincount(labels.ravel(), weights=img.ravel(), minlength=n + 1)
    sat = np.bincount(labels.ravel(), weights=(img >= 0.999 * max_value).ravel().astype(float), minlength=n + 1)
    fl = segmentation.expand_labels(labels, fiji_ring)
    fiji = np.bincount(fl.ravel(), weights=img.ravel(), minlength=n + 1) / np.maximum(np.bincount(fl.ravel(), minlength=n + 1), 1)
    cy = np.bincount(labels.ravel(), weights=np.repeat(np.arange(labels.shape[0]), labels.shape[1]).astype(float),
                     minlength=n + 1) / np.maximum(area, 1)
    cx = np.bincount(labels.ravel(), weights=np.tile(np.arange(labels.shape[1]), labels.shape[0]).astype(float),
                     minlength=n + 1) / np.maximum(area, 1)
    sig = sigma if np.isfinite(sigma) else 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        nsd = sig * np.sqrt(a_c / npx + MEDIAN_VAR * a_c / np.maximum(bgpx, 1))
        nsd90 = sig * np.sqrt(P90_VAR * a_c / npx + MEDIAN_VAR * a_c / np.maximum(bgpx, 1))
    present = area[1:] > 0
    df = pd.DataFrame(dict(
        label=idx, y=cy[1:], x=cx[1:], area=area[1:], erode_px=e[1:], measured_px=npx[1:], mean=mean[1:],
        p90=p90[1:], total=full_sum[1:] - area[1:] * bg[1:], bg=bg[1:], bg_px=bgpx[1:], bg_source=src[1:],
        fiji_cell_mean=fiji[1:], saturated=sat[1:] > 0, noise_sd=nsd[1:], noise_sd_p90=nsd90[1:],
        p90_null=Z90 * sig))[present]
    df.loc[df["measured_px"] == 0, ["mean", "p90"]] = np.nan
    return df.reset_index(drop=True), info


# ---------------------------------------------------------------- channel checks

def alignment(reference: np.ndarray, moving: np.ndarray, max_search: int = 30) -> dict:
    """Shift (dy, dx) that registers `moving` onto `reference` (band-passed cross-correlation within +-max_search
    px, refined to 0.1 px), and the correlation of the two band-passed images after the shift (low = the files may
    not belong together)."""
    from skimage.registration import phase_cross_correlation

    def prep(x):  # band-pass (nucleus-sized structure only, no illumination gradient) and taper the borders
        x = np.asarray(x, np.float32)
        x = ndi.gaussian_filter(x, 2) - ndi.gaussian_filter(x, 15)
        wy, wx = np.hanning(x.shape[0]), np.hanning(x.shape[1])
        return x * np.outer(wy, wx).astype(np.float32), x

    aw, a = prep(reference)
    bw, b = prep(moving)
    # cross-correlation, but the peak is searched only within +-max_search px: regularly spaced nuclei otherwise
    # give equally good matches one nucleus spacing away
    cc = np.fft.fftshift(np.real(np.fft.ifft2(np.fft.fft2(aw) * np.conj(np.fft.fft2(bw)))))
    cy, cx = cc.shape[0] // 2, cc.shape[1] // 2
    R = int(max_search)
    win = cc[cy - R: cy + R + 1, cx - R: cx + R + 1]
    iy, ix = np.unravel_index(int(np.argmax(win)), win.shape)
    sy, sx = iy - R, ix - R
    if 0 < iy < win.shape[0] - 1 and 0 < ix < win.shape[1] - 1:  # sub-pixel: refine around the peak
        shift, _, _ = phase_cross_correlation(aw, np.roll(bw, (sy, sx), axis=(0, 1)), upsample_factor=10,
                                              normalization=None)
        if abs(shift[0]) <= 1 and abs(shift[1]) <= 1:
            sy, sx = sy + float(shift[0]), sx + float(shift[1])
    dy, dx = float(sy), float(sx)
    bs = ndi.shift(b, (dy, dx), order=1, mode="nearest")
    m = 16 + int(np.ceil(max(abs(dy), abs(dx))))
    aa, bb = a[m:-m, m:-m].ravel(), bs[m:-m, m:-m].ravel()
    r = float(np.corrcoef(aa, bb)[0, 1]) if aa.size > 10 and aa.std() > 0 and bb.std() > 0 else float("nan")
    return dict(dy=round(dy, 2), dx=round(dx, 2), shift_px=round(float(np.hypot(dy, dx)), 2), corr=round(r, 3))


def brightness_check(stack: np.ndarray, labels: np.ndarray, seg_ch: int) -> dict:
    """Mean of every channel on the nuclei and off them. The nuclear channel should be the brightest on nuclei and
    have the highest nucleus / background contrast."""
    obj = labels > 0
    out = {}
    for c in range(stack.shape[0]):
        on = float(stack[c][obj].mean()) if obj.any() else float("nan")
        off = float(np.median(stack[c][~obj])) if (~obj).any() else float("nan")
        out[c] = dict(on=on, off=off, contrast=on / off if off > 0 else float("nan"))
    seg = out[seg_ch]
    others = [c for c in out if c != seg_ch]
    ok = all(seg["on"] > out[c]["on"] for c in others)
    return dict(per_channel=out, nuclear_brightest=bool(ok))


# ---------------------------------------------------------------- one field (cached)

def measure_hash(settings: dict) -> str:
    it, ch = settings["intensity"], settings["channels"]
    keys = {k: it.get(k) for k in ("per_nucleus", "erode_px", "bg_gap_px", "bg_ring_px", "bg_min_px", "fiji_ring_px")}
    keys.update(align=ch.get("align"), align_max_correct_px=ch.get("align_max_correct_px"), version=MEASURE_VERSION)
    return hashlib.sha1(json.dumps(keys, sort_keys=True, default=str).encode()).hexdigest()[:12]


def field_cache_paths(labels_path: Path) -> tuple[Path, Path]:
    stem = labels_path.name.replace("_labels.tif", "")
    return labels_path.with_name(f"{stem}_intensity.csv"), labels_path.with_name(f"{stem}_intensity.json")


def _sig(p: Path) -> list:
    st = Path(p).stat()
    return [str(p), st.st_size, st.st_mtime_ns]


def measure_field(rec, labels_path: Path, info: dict, settings: dict, chans: list[int], names: list[str],
                  seg_ch: int, read=None) -> tuple[pd.DataFrame, dict]:
    """Measure the marker channels `chans` of one field (cached next to the detection cache).
    Returns (per-nucleus table with a `channel` column, field info with alignment / brightness / noise)."""
    csv_p, json_p = field_cache_paths(labels_path)
    key = dict(hash=measure_hash(settings), labels=_sig(labels_path), sources=[_sig(p) for p in rec.sources],
               chans=chans)
    if csv_p.exists() and json_p.exists():
        meta = json.loads(json_p.read_text())
        if meta.get("key") == json.loads(json.dumps(key)):
            df = pd.read_csv(csv_p) if csv_p.stat().st_size > 1 else pd.DataFrame()
            return df, meta["info"]
    import tifffile

    from . import imageio

    im = (read or imageio.read_record)(rec)
    labels = tifffile.imread(labels_path).astype(np.int32)
    it, chs = settings["intensity"], settings["channels"]
    diameter = float(info.get("diameter_used") or 30.0)
    finfo: dict = dict(max_value=float(im.max_value), shape=list(labels.shape))
    if len(rec.sources) > 1 or im.data.shape[0] > 1:
        finfo["brightness"] = brightness_check(im.data, labels, seg_ch)
    tables = []
    for c in chans:
        ch = im.data[c]
        a = None
        if rec.paths and str(chs.get("align", "check")).lower() != "off":
            a = alignment(im.data[seg_ch], ch)
            a["corrected"] = False
            lim = float(chs.get("align_max_correct_px", 10))
            if str(chs.get("align")).lower() == "correct" and 1.0 <= a["shift_px"] <= lim:
                ch = ndi.shift(ch, (round(a["dy"]), round(a["dx"])), order=0, mode="nearest")
                a["corrected"] = True
            finfo.setdefault("alignment", {})[names[c]] = a
        df, ninfo = measure_nuclei(ch, labels, im.max_value, diameter, erode=it.get("erode_px", "auto"),
                                   gap=int(_num(it.get("bg_gap_px"), 3)),
                                   ring=ring_width(settings, diameter), bg_min=int(_num(it.get("bg_min_px"), 40)),
                                   fiji_ring=int(_num(it.get("fiji_ring_px"), 6)))
        finfo.setdefault("noise", {})[names[c]] = ninfo
        if len(df):
            tables.append(df.assign(channel=names[c]))
    out = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    out.to_csv(csv_p, index=False)
    json_p.write_text(json.dumps(dict(key=key, info=finfo), default=float))
    return out, finfo


# ---------------------------------------------------------------- readout ladder: whole-field layers

LAYERS_VERSION = 2
FIELD_SAT_MAX = 0.01  # more than this share of a field's pixels saturated: its field layers are left empty
STAIN_K = 3.0  # stained area = corrected image above background mode + STAIN_K x its noise SD, then opened
STAIN_SMOOTH = 0.0  # Gaussian sigma (px) before thresholding (0 = none: smoothing adds a halo around every cell)
STAIN_MIN_PX = 50  # fewer stained pixels: no stained-area mean (total = 0)


def background_mode(values: np.ndarray, bins: int = 512) -> tuple[float, float]:
    """(mode, sd) of the background peak of a fluorescence image without any mask: the mode of a smoothed
    histogram between the 0.5th and 99.5th percentiles, and the SD from the half below the mode (a half-normal:
    median(mode - v) = 0.6745 sd). Robust while the background is the most common intensity (true unless a field is
    almost fully covered by stained cells)."""
    v = np.asarray(values, np.float64).ravel()
    v = v[np.isfinite(v)]
    if v.size < 100:
        return float("nan"), float("nan")
    lo, hi = np.percentile(v, [0.5, 99.5])
    if not hi > lo:
        return float(lo), 0.0
    h, edges = np.histogram(v, bins=bins, range=(lo, hi))
    h = ndi.gaussian_filter1d(h.astype(float), 2)
    i = int(np.argmax(h))
    mode = float(0.5 * (edges[i] + edges[i + 1]))
    below = mode - v[v < mode]
    sd = float(np.median(below) / 0.6745) if below.size >= 20 else float(np.std(v))
    return mode, sd


def field_layers(img: np.ndarray, max_value: float, flat: np.ndarray | None = None, dark: np.ndarray | None = None
                 ) -> dict:
    """Whole-field layers of one marker image (no segmentation, no DAPI; docs/intensity.md "Readout ladder").

        field_raw          mean of the image as acquired (saturated pixels left out)
        field_corrected    mean of C = (image - D) / F - b over the field: flat field F and dark field D of the
                           plate (BaSiC model, F = 1 and D = 0 without a flat field), b = the field's own background
                           (mode of the smoothed corrected image)
        stained_area_mean  mean of C over the stained area: C > STAIN_K x background SD, opened (3 x 3 cross)
        field_total        sum of C over the stained area (0 when nothing is stained)
        stained_frac, saturated_frac, bg_mode, bg_sd  (helpers)
    Fields with more than FIELD_SAT_MAX saturated pixels get NaN layers."""
    img = np.asarray(img, np.float32)
    sat = img >= 0.999 * float(max_value)
    ok = ~sat
    out = dict(saturated_frac=float(sat.mean()))
    nan = dict(field_raw=np.nan, field_corrected=np.nan, stained_area_mean=np.nan, field_total=np.nan,
               stained_frac=np.nan, bg_mode=np.nan, bg_sd=np.nan)
    if out["saturated_frac"] > FIELD_SAT_MAX or ok.sum() < 100:
        return out | nan
    F = np.ones_like(img) if flat is None else np.asarray(flat, np.float32)
    D = 0.0 if dark is None else np.asarray(dark, np.float32)
    c0 = (img - D) / np.maximum(F, 1e-3)
    sm = ndi.gaussian_filter(c0, STAIN_SMOOTH) if STAIN_SMOOTH > 0 else c0
    b, sd = background_mode(sm[ok])
    C = c0 - b
    stained = ok & (sm - b > STAIN_K * sd) if np.isfinite(sd) else np.zeros_like(ok)
    # opening with a 3 x 3 cross: drops isolated noise pixels (about 0.1% of background pixels pass 3 SD) and
    # one-pixel slivers; stained cells (many pixels across) keep their shape
    stained = ndi.binary_opening(stained, structure=ndi.generate_binary_structure(2, 1))
    k = int(stained.sum())
    out.update(field_raw=float(img[ok].mean()), field_corrected=float(C[ok].mean()),
               stained_area_mean=float(C[stained].mean()) if k >= STAIN_MIN_PX else np.nan,
               field_total=float(C[stained].sum()) if k >= STAIN_MIN_PX else 0.0,
               stained_frac=float(k / ok.sum()), bg_mode=float(b), bg_sd=float(sd))
    return out


def _ff_key(ff) -> str:
    return hashlib.sha1(json.dumps(ff.to_dict(), sort_keys=True, default=float).encode()).hexdigest()[:12]


def field_layers_cached(rec, labels_path: Path, c: int, name: str, ff, max_value: float, grids: dict) -> dict:
    """field_layers for channel c of one field, cached next to the detection cache (NN_layers.json): redone when
    the image files, the flat field or LAYERS_VERSION change. `grids` caches the full F / D maps per flat field."""
    stem = labels_path.name.replace("_labels.tif", "")
    p = labels_path.with_name(f"{stem}_layers.json")
    srcs = [rec.paths[c]] if rec.paths else [rec.path]
    key = dict(version=LAYERS_VERSION, sources=[_sig(Path(x)) for x in srcs], ff=_ff_key(ff), channel=int(c),
               sat=FIELD_SAT_MAX, k=STAIN_K, smooth=STAIN_SMOOTH)
    data = json.loads(p.read_text()) if p.exists() else {}
    hit = data.get(name)
    if hit and hit.get("key") == json.loads(json.dumps(key)):
        return hit["values"]
    img = channel_image(rec, c)
    gid = id(ff)
    if gid not in grids:
        on = ff.source != "off" and tuple(ff.shape) == tuple(img.shape)
        grids[gid] = (ff.flat() if on else None, ff.dark() if on else None)
    flat, dark = grids[gid]
    vals = field_layers(img, max_value, flat, dark)
    data[name] = dict(key=key, values=vals)
    p.write_text(json.dumps(data, default=float))
    return vals


# ---------------------------------------------------------------- project stage (called by pipeline.run)

def _slug(s: str) -> str:
    import re

    return re.sub(r"[^\w.\-]+", "_", str(s)).strip("_")


FLAT_VERSION = 1


def channel_image(rec, c: int) -> np.ndarray:
    """Channel c of a field, reading only the file needed."""
    from . import imageio

    if rec.paths:
        return imageio.read(rec.paths[c]).data[0]
    return imageio.read(rec.path).data[c]


def flat_field(prj, n: str, c: int, name: str, recs: list, shape, log=print):
    """The flat field of channel c in n: cached in _work/flatfield/. [intensity] flatfield = "auto" | "off" | folder."""
    from . import flatfield as ff_mod

    it = prj.settings["intensity"]
    mode = str(it.get("flatfield", "auto")).strip()
    if mode.lower() in ("off", "", "false", "none"):
        return ff_mod.FlatField.identity(shape)
    srcs = sorted(str(r.paths[c]) if r.paths else str(r.path) for r in recs)
    key = dict(mode=mode, dark=bool(it.get("flatfield_dark", True)), version=FLAT_VERSION, shape=list(shape),
               sources=[_sig(Path(p)) for p in srcs])
    cache = prj.work / "flatfield" / f"{_slug(n)}__{_slug(name)}.json"
    if cache.exists():
        d = json.loads(cache.read_text())
        if d.get("key") == json.loads(json.dumps(key)):
            return ff_mod.FlatField.from_dict(d["flat"])
    if mode.lower() == "auto":
        log(f"  flat field {n} / {name}: estimating from {len(recs)} images ...")
        ff = ff_mod.estimate((channel_image(r, c) for r in recs), dark=bool(it.get("flatfield_dark", True)))
        if ff.source == "off":
            log(f"  note: flat field {n} / {name}: too few usable images, no correction")
    else:
        folder = Path(mode)
        folder = folder if folder.is_absolute() else (prj.path / folder)
        refs, darks = ff_mod.reference_files(folder, name)
        if not refs:
            raise SystemExit(f"[intensity] flatfield = {mode!r}: no reference images (TIFF) found in {folder}")
        ff = ff_mod.from_reference(refs, shape, darks)
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(dict(key=key, flat=ff.to_dict())))
    return ff


def positivity_channels(settings: dict, names: list[str]) -> list[str]:
    p = settings.get("positivity", {})
    if not p.get("enabled", True):
        return []
    want = p.get("channels", "auto")
    if str(want).lower() == "auto":
        return list(names)
    want = want if isinstance(want, list) else str(want).replace(";", ",").split(",")
    low = {str(w).strip().lower() for w in want}
    return [n for n in names if n.lower() in low]


def derive(nuc: pd.DataFrame, ffs: dict, k: float) -> pd.DataFrame:
    """Add flat field, background-subtracted / corrected signals, counted / measurable flags and LoD cut-offs."""
    from . import classify, positivity

    nuc = nuc.copy()
    ffv = np.ones(len(nuc))
    for (n, ch), ff in ffs.items():
        m = ((nuc["n"] == n) & (nuc["channel"] == ch)).to_numpy()
        if m.any():
            ffv[m] = ff.flat(nuc.loc[m, "y"].to_numpy(float), nuc.loc[m, "x"].to_numpy(float))
    nuc["ff"] = ffv
    nuc["mean_bgsub"] = nuc["mean"] - nuc["bg"]
    nuc["p90_bgsub"] = nuc["p90"] - nuc["bg"]
    nuc["mean_bgsub_ff"] = nuc["mean_bgsub"] / nuc["ff"]
    nuc["p90_bgsub_ff"] = nuc["p90_bgsub"] / nuc["ff"]
    nuc["total_bgsub_ff"] = nuc["total"] / nuc["ff"]
    # raw mean over the whole (unshrunk) nucleus: no background, no flat field (Fiji ROI mean; total = sum - area x bg)
    nuc["mean_raw"] = nuc["total"] / nuc["area"].where(nuc["area"] > 0) + nuc["bg"]
    nuc["counted"] = nuc["category"].isin(classify.COUNTED)
    nuc["measurable"] = nuc["counted"] & ~nuc["saturated"].astype(bool) & (nuc["measured_px"] > 0)
    nuc["thr_lod"], nuc["thr90_lod"] = positivity.lod_thresholds(nuc, k, nuc["ff"].to_numpy(float))
    nuc["thr_lod_raw"], nuc["thr90_lod_raw"] = positivity.lod_thresholds(nuc, k)
    return nuc


SUSPECT_SHARE = 0.5  # at least this share of a folder's fields failing the channel check: the folder is suspect


def field_fails_check(fi: dict) -> bool:
    """A field whose marker file does not look like the partner of its nuclear file."""
    b = fi.get("brightness") or {}
    mis = any(np.isfinite(a.get("corr", np.nan)) and a["corr"] < MISPAIR_CORR
              for a in (fi.get("alignment") or {}).values())
    return bool(mis or (b and not b.get("nuclear_brightest", True)))


def suspect_folders(finfos: list) -> dict:
    """{(n, treatment): (failing fields, fields)} for folders where at least SUSPECT_SHARE of the fields (and at
    least 2) fail the channel check."""
    tally: dict = {}
    for r, _, fi, _ in finfos:
        t = tally.setdefault((r.n, r.treatment), [0, 0])
        t[0] += field_fails_check(fi)
        t[1] += 1
    return {k: (b, n) for k, (b, n) in tally.items() if b >= 2 and b >= SUSPECT_SHARE * n}


def process(prj, fields: list, log=print) -> dict | None:
    """Measure marker channels per nucleus for every field, estimate flat fields, decide cut-offs, add per-image
    columns to the rows (in place) and write the intensity outputs. Does nothing (returns None) when the images
    have no marker channel (single-channel DAPI projects)."""
    from . import positivity

    s = prj.settings
    if not fields:
        return None
    nch = int(fields[0][3].get("n_channels", 1))
    try:
        chans = marker_channels(s, nch)
    except ValueError as e:
        log(f"  WARNING: [intensity] per_nucleus: {e}; no per-nucleus intensities")
        return None
    if not chans:
        return None
    names = settings_mod.channel_names(s, nch)
    seg_ch = prj.seg_channel(nch)
    log(f"Measuring {', '.join(names[c] for c in chans)} per nucleus (local background, flat field, positivity)...")

    tables, finfos = [], []
    for r, lab, cobjs, info, row in fields:
        df, fi = measure_field(r, lab, info, s, chans, names, seg_ch)
        if len(df):
            cat = cobjs.set_index("label")["category"] if len(cobjs) and "category" in cobjs else pd.Series(dtype=str)
            df = df.assign(n=r.n, treatment=r.treatment, image=r.image, category=df["label"].map(cat).fillna(""))
            tables.append(df)
        finfos.append((r, row, fi, info))
    cols = ["n", "treatment", "image", "channel", "label", "category"]
    nuc = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(columns=cols + ["y", "x"])
    nuc = nuc[cols + [c for c in nuc.columns if c not in cols]]

    # flat field per n and channel
    shape = tuple(finfos[0][2].get("shape") or (1, 1))
    by_n: dict[str, list] = {}
    for r, _, _, _ in finfos:
        by_n.setdefault(r.n, []).append(r)
    ffs = {}
    for n, recs in sorted(by_n.items()):
        for c in chans:
            ffs[(n, names[c])] = flat_field(prj, n, c, names[c], recs, shape, log)
    k = float(s["positivity"].get("lod_k", 3.0))
    nuc = derive(nuc, ffs, k)
    # fields whose marker file does not show the same cells as the nuclear file (mis-paired / wrong file): the
    # count stays, the marker readouts of that field are left empty
    nuc["excluded"] = ""
    for r, _, fi, _ in finfos:
        for ch, a in (fi.get("alignment") or {}).items():
            if np.isfinite(a.get("corr", np.nan)) and a["corr"] < MISPAIR_CORR:
                m = (nuc["n"] == r.n) & (nuc["treatment"] == r.treatment) & (nuc["image"] == r.image) & (nuc["channel"] == ch)
                nuc.loc[m, "excluded"] = "marker image does not match the nuclear image"
                nuc.loc[m, "measurable"] = False

    # a whole folder whose channel files do not behave like marker + nuclear pairs (most fields mis-paired, or the
    # nuclear channel not the brightest on its nuclei): the folder may hold only one channel, or a missing / extra
    # file shifted the pairing. Its marker readouts are left empty (docs/design.md, "Unusual designs").
    suspect = suspect_folders(finfos)
    for (n, t), (bad, tot) in suspect.items():
        m = (nuc["n"] == n) & (nuc["treatment"] == t)
        nuc.loc[m, "excluded"] = "channel files of this folder do not pair up"
        nuc.loc[m, "measurable"] = False
        log(f"  WARNING: {n} / {t}: {bad} of {tot} fields fail the channel check (mis-paired files or the nuclear "
            f"channel is not the brightest): the folder may hold only one channel, or a missing / extra file shifted "
            f"the pairing. Its marker readouts are left empty; its counts come from every second file only. Check "
            f"the files (Intensity → Channel checks).")

    # cut-offs per n and channel
    pos_ch = positivity_channels(s, [names[c] for c in chans])
    locks = positivity.load_locks(prj.path)
    nuc["thr"], nuc["thr_p90"] = np.nan, np.nan
    nuc["positive"] = pd.Series(pd.NA, index=nuc.index, dtype="boolean")
    nuc["positive_p90"] = pd.Series(pd.NA, index=nuc.index, dtype="boolean")
    cut_rows, sweeps = [], []
    for (n, ch), g in nuc[nuc["measurable"]].groupby(["n", "channel"], sort=True):
        if ch not in pos_ch:
            continue
        dec = positivity.decide(g, s, ch, n, locks)
        nuc.loc[g.index, "thr"] = dec["mean"]["thr"]
        nuc.loc[g.index, "thr_p90"] = dec["p90"]["thr"]
        nuc.loc[g.index, "positive"] = g["mean_bgsub_ff"].to_numpy() > dec["mean"]["thr"]
        nuc.loc[g.index, "positive_p90"] = g["p90_bgsub_ff"].to_numpy() > dec["p90"]["thr"]
        # the same rule without the flat field (before / after)
        raw = g.assign(mean_bgsub_ff=g["mean_bgsub"], p90_bgsub_ff=g["p90_bgsub"], thr_lod=g["thr_lod_raw"],
                       thr90_lod=g["thr90_lod_raw"])
        dec_raw = positivity.decide(raw, s, ch, n, locks)
        before = {"mean": float((raw["mean_bgsub_ff"].to_numpy() > dec_raw["mean"]["thr"]).mean()),
                  "p90": float((raw["p90_bgsub_ff"].to_numpy() > dec_raw["p90"]["thr"]).mean())}
        for measure, pcol in (("mean", "positive"), ("p90", "positive_p90")):
            d = dec[measure]
            cut_rows.append(dict(n=n, channel=ch, measure=measure, rule=d["rule"],
                                 rule_label=positivity.RULE_LABELS[d["rule"]] + (
                                     " (NOT the manual cut-off: that applies to the mean only)" if d.get("note") else ""),
                                 value=round(d["value"], 3),
                                 nuclei=len(g), frac_positive=round(float(nuc.loc[g.index, pcol].astype(float).mean()), 4),
                                 frac_positive_without_flatfield=round(before[measure], 4),
                                 flatfield=ffs[(n, ch)].source, detail=d["detail"]))
        sweeps.append(positivity.sweep(g, dec["mean"]["thr"]).assign(n=n, channel=ch, measure="mean", rule=dec["mean"]["rule"]))
        sweeps.append(positivity.sweep(g, dec["p90"]["thr"], signal="p90_bgsub_ff").assign(
            n=n, channel=ch, measure="p90", rule=dec["p90"]["rule"]))
    cutoffs = pd.DataFrame(cut_rows)
    sweep_df = pd.concat(sweeps, ignore_index=True) if sweeps else pd.DataFrame()
    if len(sweep_df):
        sweep_df = sweep_df[["n", "treatment", "channel", "measure", "rule", "factor", "cutoff", "frac_positive", "nuclei"]]

    # readout ladder: whole-field layers per field and channel (cheap: no segmentation; cached per field)
    layer_vals: dict = {}
    grids: dict = {}
    for r, row, fi, info in finfos:
        mv = float(fi.get("max_value") or info.get("max_value") or 65535)
        lab = next(f[1] for f in fields if f[0] is r)
        for c in chans:
            try:
                layer_vals[(r.n, r.treatment, r.image, names[c])] = field_layers_cached(
                    r, lab, c, names[c], ffs[(r.n, names[c])], mv, grids)
            except Exception as e:  # noqa: BLE001 - a layer must never stop a run
                log(f"  note: field layers {r.n} / {r.treatment} / {r.image}: {type(e).__name__}: {e}")

    # per-image columns
    minn = int(s["positivity"].get("min_nuclei", 5))
    regions = s["intensity"].get("regions", [])
    groups = {key: g for key, g in nuc.groupby(["n", "treatment", "image", "channel"], sort=False)}
    warn = float(s["channels"].get("align_warn_px", 2.0))
    for r, row, fi, info in finfos:
        flags = []
        b = fi.get("brightness")
        if b and not b.get("nuclear_brightest", True):
            flags.append(f"nuclear channel {names[seg_ch]} is not the brightest on its nuclei (channel order?)")
        for ch, a in (fi.get("alignment") or {}).items():
            if a["shift_px"] > warn:
                flags.append(f"{ch} shifted {a['shift_px']:.1f} px vs {names[seg_ch]}" + (" (corrected)" if a.get("corrected") else ""))
            if np.isfinite(a.get("corr", np.nan)) and a["corr"] < MISPAIR_CORR:
                flags.append(f"{ch} barely resembles {names[seg_ch]} (r = {a['corr']:.2f}): files mis-paired? Marker readouts left empty")
        if (r.n, r.treatment) in suspect:
            bad, tot = suspect[(r.n, r.treatment)]
            flags.append(f"folder fails the channel check ({bad} of {tot} fields): only one channel, or files "
                         "mis-paired? Marker readouts of the whole folder left empty")
        if fi.get("alignment") or b:
            row["channel_check"] = "; ".join(flags)
        for c in chans:
            ch = names[c]
            sl = _slug(ch)
            g = groups.get((r.n, r.treatment, r.image, ch), nuc.iloc[0:0])
            meas = g[g["measurable"]]
            row[f"n_measured_{sl}"] = int(len(meas))
            row[f"n_saturated_{sl}"] = int((g["counted"] & g["saturated"].astype(bool)).sum())
            ok = len(meas) >= minn
            mispaired = bool(len(g)) and (g["excluded"] != "").all()
            lv = {} if mispaired else layer_vals.get((r.n, r.treatment, r.image, ch), {})

            def lay(k, nd=3):
                v = lv.get(k, np.nan)
                return round(float(v), nd) if v is not None and np.isfinite(v) else np.nan

            # the ladder, in order (layers 1-4; layer 5 = % positive below), then the two stain-area readouts
            row[f"intensity_{sl}_field_raw"] = lay("field_raw")
            row[f"intensity_{sl}_field_corrected"] = lay("field_corrected")
            row[f"intensity_{sl}_nucleus_mean_raw"] = round(float(meas["mean_raw"].median()), 3) if ok else np.nan
            row[f"intensity_{sl}_nucleus_mean_bgsub"] = round(float(meas["mean_bgsub_ff"].median()), 3) if ok else np.nan
            row[f"intensity_{sl}_stained_area_mean"] = lay("stained_area_mean")
            row[f"intensity_{sl}_field_total"] = lay("field_total", 1)
            row[f"stained_frac_{sl}"] = lay("stained_frac", 5)
            row[f"saturated_px_frac_{sl}"] = lay("saturated_frac", 5)
            if "fiji_cell" in regions:
                cnt = g[g["counted"]]
                bgv = float(info.get(f"C{c + 1}_background", 0.0)) if s["intensity"].get("subtract_background", True) else 0.0
                row[f"intensity_{sl}_fiji_cell"] = float((cnt["fiji_cell_mean"] - bgv).mean()) if len(cnt) else np.nan
            if ch in pos_ch:
                posm = meas["positive"].astype(bool).to_numpy()
                row[f"readout_frac_positive_{sl}"] = round(float(posm.mean()), 5) if ok else np.nan
                row[f"readout_frac_positive_p90_{sl}"] = round(float(meas["positive_p90"].astype(bool).mean()), 5) if ok else np.nan
                row[f"intensity_{sl}_positive_mean_bgsub"] = round(float(meas.loc[posm, "mean_bgsub_ff"].mean()), 3) \
                    if ok and posm.any() else np.nan
                rule = cutoffs[(cutoffs["n"] == r.n) & (cutoffs["channel"] == ch) & (cutoffs["measure"] == "mean")]["rule"] \
                    if len(cutoffs) else pd.Series(dtype=str)
                row[f"cutoff_rule_{sl}"] = rule.iloc[0] if len(rule) else ""

    flagged = [(r, row["channel_check"]) for r, row, _, _ in finfos if row.get("channel_check")]
    for r, msg in flagged:
        log(f"  channel check {r.n} / {r.treatment} / field {r.image:02d} "
            f"(files {', '.join(str(x) for x in r.numbers) or r.path.name}): {msg}")

    from . import intensity_report

    nuc.to_csv(prj.work / "nuclei_intensity.csv", index=False)
    intensity_report.write(prj, nuc, finfos, ffs, cutoffs, sweep_df, names, chans, seg_ch, log)
    return dict(nuclei=nuc, cutoffs=cutoffs, sweep=sweep_df, flatfields=ffs)
