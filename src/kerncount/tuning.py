"""Tuning API for the browser: set a positivity cut-off by eye, blind to condition (docs/intensity.md, "Tuning").

Pure functions on a project that has been run (`kerncount run`): they read `_work/nuclei_intensity.csv`, the
detection labels and the flat fields from `_work/flatfield/`. Nothing here re-detects or re-measures, and only
`lock` writes (cutoffs.json). Results are cached in memory per project (invalidated when the run's output files
change), so after the first call per field everything is fast enough for a live slider.

Blinding: fields are identified by opaque ids (a salted hash, salt in _work/tuning_salt). No function returns a
condition (treatment) name, an n-specific breakdown by condition, or a file path. `preview` gives one overall %
positive. Lock the cut-off before looking at any statistics.

Units: the cut-off is in corrected intensity units = mean inside the (shrunk) nucleus minus the local background,
divided by the flat field (`mean_bgsub_ff` in the per-nucleus table); `measure="p90"` uses the 90th percentile.
"""
from __future__ import annotations

import hashlib
import io
import json
import secrets
import threading
from pathlib import Path

import numpy as np
import pandas as pd

_LOCK = threading.Lock()
_PROJECTS: dict[str, "_Project"] = {}
_FIELD_CACHE_MAX = 64

COLOUR_POSITIVE = (40, 220, 60)
COLOUR_NEGATIVE = (150, 150, 150)
COLOUR_EXCLUDED = (255, 150, 0)  # counted but saturated: excluded from intensity readouts


class TuningError(ValueError):
    """Bad input (unknown field / channel) or a project that has not been run with a marker channel yet."""


def _sig(p: Path):
    try:
        st = p.stat()
        return (st.st_size, st.st_mtime_ns)
    except OSError:
        return None


class _Project:
    def __init__(self, path: Path):
        from . import discover, pipeline

        self.path = path
        self.prj = pipeline.Project(path)
        table = path / "_work" / "nuclei_intensity.csv"
        if not table.exists():
            raise TuningError("no per-nucleus intensities yet: run the analysis first (kerncount run), with a marker "
                              "channel ([intensity] per_nucleus)")
        self.key = self.signature(path)
        nuc = pd.read_csv(table, low_memory=False)
        for c in ("measurable", "counted", "saturated"):
            nuc[c] = nuc[c].astype(str).str.lower().isin(["true", "1"])
        self.nuc = nuc
        self.channels = list(pd.unique(nuc["channel"].astype(str)))
        salt_p = path / "_work" / "tuning_salt"
        if not salt_p.exists():
            salt_p.write_text(secrets.token_hex(16))
        salt = salt_p.read_text().strip()
        layout = discover.read_table(path / "layout.csv")
        recs, _ = self.prj.records(layout)
        self.fields: dict[str, dict] = {}
        have = set(map(tuple, nuc[["n", "treatment", "image"]].drop_duplicates().astype(str).to_numpy()))
        for r in recs:
            if r.path is None or (str(r.n), str(r.treatment), str(r.image)) not in have:
                continue
            fid = "f" + hashlib.sha1(f"{salt}|{r.n}|{r.treatment}|{r.image}".encode()).hexdigest()[:10]
            lab, _, meta = pipeline._cache_paths(self.prj, r)
            self.fields[fid] = dict(rec=r, n=r.n, treatment=r.treatment, image=int(r.image), labels=lab, meta=meta)
        self.by_key = {(f["n"], f["treatment"], f["image"]): fid for fid, f in self.fields.items()}
        nuc["field_id"] = [self.by_key.get((n, t, int(i)), "") for n, t, i in
                           zip(nuc["n"].astype(str), nuc["treatment"].astype(str), nuc["image"])]
        self.field_cache: dict[tuple, dict] = {}

    @staticmethod
    def signature(path: Path):
        return (_sig(path / "_work" / "nuclei_intensity.csv"), _sig(path / "layout.csv"), _sig(path / "settings.toml"))

    def channel(self, channel: str | None) -> str:
        if channel is None:
            return self.channels[0]
        for c in self.channels:
            if c.lower() == str(channel).lower():
                return c
        raise TuningError(f"unknown channel {channel!r}; measured channels: {', '.join(self.channels)}")

    def field(self, field_id: str) -> dict:
        f = self.fields.get(field_id)
        if f is None:
            raise TuningError(f"unknown field id {field_id!r}")
        return f

    def flat(self, n: str, channel: str):
        from . import flatfield, intensity

        p = self.path / "_work" / "flatfield" / f"{intensity._slug(n)}__{intensity._slug(channel)}.json"
        if self.prj.settings["intensity"].get("flatfield", "auto") != "off" and p.exists():
            return flatfield.FlatField.from_dict(json.loads(p.read_text())["flat"])
        return None


def _project(project_dir) -> _Project:
    path = Path(project_dir).resolve()
    with _LOCK:
        p = _PROJECTS.get(str(path))
        if p is None or p.key != _Project.signature(path):
            p = _Project(path)
            _PROJECTS[str(path)] = p
        return p


def _rows(p: _Project, field_id: str, channel: str) -> pd.DataFrame:
    return p.nuc[(p.nuc["field_id"] == field_id) & (p.nuc["channel"] == channel)]


def _col(measure: str) -> str:
    if measure not in ("mean", "p90"):
        raise TuningError("measure must be 'mean' or 'p90'")
    return "mean_bgsub_ff" if measure == "mean" else "p90_bgsub_ff"


# ---------------------------------------------------------------- public API

def channels(project_dir) -> list[str]:
    """Marker channels with per-nucleus measurements."""
    return list(_project(project_dir).channels)


def sample(project_dir, k: int = 12, seed: int | None = None, channel: str | None = None, n: str | None = None) -> dict:
    """k random fields spread over all conditions (round robin over conditions, random within each), in random
    order, identified only by opaque ids. `n` restricts to one n (cut-offs are usually tuned per n: staining varies
    between plates). Returns {"seed", "channel", "n", "fields": [{"field_id", "nuclei"}]}; `nuclei` = measurable
    nuclei in the field."""
    p = _project(project_dir)
    ch = p.channel(channel)
    if seed is None:
        seed = int(np.random.default_rng().integers(0, 2 ** 31 - 1))
    rng = np.random.default_rng(seed)
    cands = {fid: f for fid, f in p.fields.items() if n is None or f["n"] == n}
    if not cands:
        raise TuningError(f"no fields{' in ' + n if n else ''}")
    by_cond: dict[tuple, list[str]] = {}
    for fid in sorted(cands):
        f = cands[fid]
        by_cond.setdefault((f["n"], f["treatment"]), []).append(fid)
    queues = [list(rng.permutation(v)) for _, v in sorted(by_cond.items())]
    order = rng.permutation(len(queues))
    picked: list[str] = []
    while len(picked) < min(k, len(cands)):
        for i in order:
            if queues[i] and len(picked) < k:
                picked.append(str(queues[i].pop()))
    picked = [picked[i] for i in rng.permutation(len(picked))]
    counts = p.nuc[p.nuc["measurable"] & (p.nuc["channel"] == ch)].groupby("field_id").size()
    return dict(seed=int(seed), channel=ch, n=n, fields=[dict(field_id=f, nuclei=int(counts.get(f, 0))) for f in picked])


def _field_geometry(p: _Project, field_id: str) -> dict:
    key = ("geom", field_id)
    g = p.field_cache.get(key)
    if g is None:
        import tifffile
        from skimage import segmentation

        f = p.field(field_id)
        labels = tifffile.imread(f["labels"]).astype(np.int32)
        bnd = segmentation.find_boundaries(labels, mode="inner") & (labels > 0)
        flat_idx = np.flatnonzero(bnd)
        g = dict(labels=labels, shape=labels.shape, bidx=flat_idx, blab=labels.ravel()[flat_idx])
        _put(p, key, g)
    return g


def _put(p: _Project, key, val):
    if len(p.field_cache) >= _FIELD_CACHE_MAX:
        p.field_cache.pop(next(iter(p.field_cache)))
    p.field_cache[key] = val


def nuclei(project_dir, field_id: str, channel: str | None = None, contours: bool = True) -> dict:
    """Per-nucleus values of one field: {"field_id", "channel", "shape": [H, W], "nuclei": [{label, x, y, area,
    signal (mean - local bg, flat-field corrected), signal_p90, lod_cutoff (this nucleus' automatic cut-off),
    measurable (counted and not saturated), saturated, contour: [[x, y], ...] (pixel coordinates, simplified)}]}.
    Only counted nuclei are listed."""
    p = _project(project_dir)
    ch = p.channel(channel)
    key = ("nuclei", field_id, ch, contours)
    hit = p.field_cache.get(key)
    if hit is not None:
        return hit
    rows = _rows(p, field_id, ch)
    rows = rows[rows["counted"]]
    geo = _field_geometry(p, field_id)
    cont = {}
    if contours:
        from skimage import measure

        labels = geo["labels"]
        objs = __import__("scipy.ndimage", fromlist=["find_objects"]).find_objects(labels)
        for lab in rows["label"].astype(int):
            sl = objs[lab - 1] if lab - 1 < len(objs) else None
            if sl is None:
                continue
            y0, x0 = max(sl[0].start - 1, 0), max(sl[1].start - 1, 0)
            sub = labels[y0:sl[0].stop + 1, x0:sl[1].stop + 1] == lab
            cs = measure.find_contours(np.pad(sub, 1).astype(float), 0.5)
            if cs:
                c = measure.approximate_polygon(max(cs, key=len), tolerance=0.8)
                cont[lab] = [[round(float(x) + x0 - 1, 1), round(float(y) + y0 - 1, 1)] for y, x in c]

    def num(v):
        return None if pd.isna(v) else round(float(v), 3)

    out = dict(field_id=field_id, channel=ch, shape=list(geo["shape"]), nuclei=[
        dict(label=int(r.label), x=round(float(r.x), 1), y=round(float(r.y), 1), area=int(r.area),
             signal=num(r.mean_bgsub_ff), signal_p90=num(r.p90_bgsub_ff), lod_cutoff=num(r.thr_lod),
             measurable=bool(r.measurable), saturated=bool(r.saturated), contour=cont.get(int(r.label), []))
        for r in rows.itertuples()])
    _put(p, key, out)
    return out


def _base_image(p: _Project, field_id: str, ch: str) -> np.ndarray:
    key = ("base", field_id, ch)
    b = p.field_cache.get(key)
    if b is None:
        from . import intensity, settings as settings_mod

        f = p.field(field_id)
        meta = json.loads(f["meta"].read_text())
        nch = int(meta.get("n_channels", 1))
        c = [x.lower() for x in settings_mod.channel_names(p.prj.settings, nch)].index(ch.lower())
        img = intensity.channel_image(f["rec"], c).astype(np.float32)
        ff = p.flat(f["n"], ch)
        if ff is not None and tuple(ff.shape) == img.shape:
            img = ff.correct(img)
        lo, hi = np.percentile(img[::2, ::2], [0.5, 99.8])
        g = np.clip((img - lo) / max(hi - lo, 1e-6), 0, 1)
        g = (np.sqrt(g) * 255).astype(np.uint8)  # gamma 0.5: dim nuclei stay visible next to bright ones
        b = np.repeat(g[:, :, None], 3, axis=2)
        _put(p, key, b)
    return b


def render(project_dir, field_id: str, channel: str | None = None, cutoff: float | None = None,
           measure: str = "mean", outline: bool = True) -> bytes:
    """PNG of the marker image of one field (flat-field corrected, contrast-stretched, gamma 0.5) with every counted
    nucleus outlined: green = positive at `cutoff`, grey = negative, orange = saturated (excluded). cutoff=None
    uses each nucleus' current cut-off (the rule in force after the last run)."""
    from PIL import Image as PILImage

    p = _project(project_dir)
    ch = p.channel(channel)
    img = _base_image(p, field_id, ch).copy()
    if outline:
        rows = _rows(p, field_id, ch)
        rows = rows[rows["counted"]]
        geo = _field_geometry(p, field_id)
        n = int(geo["labels"].max())
        colour = np.zeros((n + 1, 3), np.uint8)
        show = np.zeros(n + 1, bool)
        sig = rows[_col(measure)].to_numpy(float)
        thr = np.full(len(rows), float(cutoff)) if cutoff is not None else \
            rows["thr" if measure == "mean" else "thr_p90"].to_numpy(float)
        pos = sig > thr
        labs = rows["label"].to_numpy(int)
        colour[labs] = np.where(pos[:, None], COLOUR_POSITIVE, COLOUR_NEGATIVE)
        colour[labs[rows["saturated"].to_numpy(bool)]] = COLOUR_EXCLUDED
        show[labs] = True
        m = show[geo["blab"]]
        flat = img.reshape(-1, 3)
        flat[geo["bidx"][m]] = colour[geo["blab"][m]]
    buf = io.BytesIO()
    PILImage.fromarray(img).save(buf, format="PNG", compress_level=1)
    return buf.getvalue()


def preview(project_dir, cutoff: float, channel: str | None = None, n: str | None = None, measure: str = "mean",
            field_ids: list[str] | None = None) -> dict:
    """Overall % positive at `cutoff` over all measurable nuclei (of one n, or of the given fields): no
    per-condition breakdown. Returns {"cutoff", "channel", "measure", "nuclei", "positive", "frac_positive",
    "lod_frac_positive" (the automatic rule, for orientation)}."""
    p = _project(project_dir)
    ch = p.channel(channel)
    d = p.nuc[p.nuc["measurable"] & (p.nuc["channel"] == ch)]
    if n is not None:
        d = d[d["n"] == n]
    if field_ids is not None:
        d = d[d["field_id"].isin(set(field_ids))]
    x = d[_col(measure)].to_numpy(float)
    pos = int((x > float(cutoff)).sum())
    lod = d["thr_lod" if measure == "mean" else "thr90_lod"].to_numpy(float)
    return dict(cutoff=float(cutoff), channel=ch, measure=measure, nuclei=int(len(x)), positive=pos,
                frac_positive=pos / len(x) if len(x) else None,
                lod_frac_positive=float((x > lod).mean()) if len(x) else None)


def suggest_range(project_dir, channel: str | None = None, measure: str = "mean") -> dict:
    """Slider range for the UI: 1st and 99th percentile of the signal, and the median automatic (LoD) cut-off."""
    p = _project(project_dir)
    ch = p.channel(channel)
    d = p.nuc[p.nuc["measurable"] & (p.nuc["channel"] == ch)]
    x = d[_col(measure)].to_numpy(float)
    lod = d["thr_lod" if measure == "mean" else "thr90_lod"].to_numpy(float)
    return dict(channel=ch, measure=measure, min=float(np.nanpercentile(x, 1)), max=float(np.nanpercentile(x, 99)),
                lod_median=float(np.nanmedian(lod)))


def lock(project_dir, channel: str, cutoff: float, field_ids: list[str], measure: str = "mean",
         n: str | None = None) -> dict:
    """Lock a manual cut-off for `channel` (and `n`; None = every n) in cutoffs.json with the time and the fields
    (opaque ids + their image files, for the record) it was tuned on. It takes effect at the next `kerncount run`
    (detections and measurements are reused from the cache, so this is quick). Precedence: negative-control wells
    > a settings.toml manual_cutoff > this lock > automatic rule. Returns the stored lock."""
    from . import positivity

    p = _project(project_dir)
    ch = p.channel(channel)
    _col(measure)
    files = []
    for fid in field_ids:
        f = p.field(fid)
        files.append(dict(field_id=fid, files=[str(x) for x in f["rec"].sources]))
    return positivity.save_lock(p.path, ch, float(cutoff), list(field_ids), measure=measure, n=n,
                                extra=dict(field_files=files, units="mean - local background, / flat field"))
