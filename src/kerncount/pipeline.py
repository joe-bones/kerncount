"""End-to-end run: discover -> segment (cached) -> classify/QC -> review images -> counts.csv -> generations."""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from . import classify, discover, generations, imageio, master, overlay, segment
from . import settings as settings_mod

MEASURE_VERSION = 2  # bump when per-object measurements change, so cached results are recomputed


def _slug(s: str) -> str:
    return re.sub(r"[^\w.\-]+", "_", s).strip("_")


def log(msg: str):
    print(msg, flush=True)


class Project:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.settings = settings_mod.load(self.path)
        self.work = self.path / "_work"
        self.review = self.path / "review"
        self.delim = self.settings["project"]["csv_delimiter"]

    def raw_dir(self) -> Path:
        raw = self.settings["project"]["raw_dir"]
        if not raw:
            raise SystemExit("settings.toml: project.raw_dir is empty")
        p = Path(raw)
        return p if p.is_absolute() else (self.path / p).resolve()

    def channel_name(self, c: int) -> str:
        names = self.settings["channels"]["names"]
        return names[c] if c < len(names) else f"C{c + 1}"


# ---------------------------------------------------------------- segmentation cache

def _cache_paths(prj: Project, rec: discover.ImageRecord) -> tuple[Path, Path, Path]:
    d = prj.work / _slug(rec.n) / _slug(rec.treatment)
    stem = f"{rec.image:02d}"
    return d / f"{stem}_labels.tif", d / f"{stem}_objects.csv", d / f"{stem}_meta.json"


def _source_sig(p: Path) -> dict:
    st = p.stat()
    return {"source": str(p), "size": st.st_size, "mtime": st.st_mtime_ns}


def _cached(prj: Project, rec, seg_hash: str) -> bool:
    lab, obj, meta = _cache_paths(prj, rec)
    if not (lab.exists() and obj.exists() and meta.exists()):
        return False
    m = json.loads(meta.read_text())
    sig = _source_sig(rec.path)
    return m.get("seg_hash") == seg_hash and m.get("size") == sig["size"] and m.get("mtime") == sig["mtime"]


def _segment_folder(prj: Project, recs: list, seg_hash: str, force: bool, counter: list[int], total: int):
    todo = [r for r in recs if force or not _cached(prj, r, seg_hash)]
    if not todo:
        counter[0] += len(recs)
        return
    s = prj.settings
    seg_ch = int(s["channels"]["segment"]) - 1
    folder_meta = prj.work / _slug(recs[0].n) / _slug(recs[0].treatment) / "folder.json"
    folder_meta.parent.mkdir(parents=True, exist_ok=True)
    diameter = None
    if str(s["segmentation"]["nucleus_diameter_px"]) == "auto":
        fm = json.loads(folder_meta.read_text()) if folder_meta.exists() else {}
        if fm.get("seg_hash") == seg_hash and len(todo) < len(recs):
            diameter = fm["diameter"]
        else:
            ests = [segment.estimate_diameter(imageio.read(r.path).data[seg_ch]) for r in recs]
            ests = [e for e in ests if e]
            diameter = float(np.median(ests)) if ests else None
            folder_meta.write_text(json.dumps({"seg_hash": seg_hash, "diameter": diameter, "per_image": ests}))
    done_ids = {id(r) for r in todo}
    for r in recs:
        counter[0] += 1
        if id(r) not in done_ids:
            continue
        t0 = time.time()
        im = imageio.read(r.path)
        if seg_ch >= im.data.shape[0]:
            raise SystemExit(f"{r.path}: has {im.data.shape[0]} channel(s) but settings ask for channel {seg_ch + 1}")
        res = segment.run(im.data, im.max_value, seg_ch, s, diameter)
        lab, obj, meta = _cache_paths(prj, r)
        tifffile.imwrite(lab, res.labels.astype(np.uint16 if res.labels.max() < 65535 else np.int32), compression="zlib")
        res.objects.to_csv(obj, index=False)
        info = {k: v for k, v in res.info.items()}
        info.update(_source_sig(r.path), seg_hash=seg_hash, pixel_size_um=im.pixel_size_um, n_channels=int(im.data.shape[0]))
        meta.write_text(json.dumps(info, default=float))
        log(f"  [{counter[0]:>4}/{total}] {r.n} / {r.treatment} / {r.image:02d}: {len(res.objects)} objects "
            f"(diam {res.info['diameter_used']:.0f}px, {time.time() - t0:.1f}s)")


# ---------------------------------------------------------------- per-image results

def _intensity(prj: Project, objs: pd.DataFrame, info: dict) -> dict:
    s = prj.settings["intensity"]
    if not s["enabled"] or objs.empty:
        return {}
    counted = objs[objs["category"].isin(classify.COUNTED)]
    nch = int(info.get("n_channels", 1))
    chans = range(nch) if s["measure"] == "all" else [int(c) - 1 for c in s["measure"] if int(c) - 1 < nch]
    out = {}
    for c in chans:
        bg = info.get(f"C{c + 1}_background", 0.0) if s["subtract_background"] else 0.0
        for region in s["regions"]:
            col = f"C{c + 1}_{region}_mean"
            if col in counted:
                out[f"intensity_{_slug(prj.channel_name(c))}_{region}"] = float((counted[col] - bg).mean()) if len(counted) else np.nan
    return out


def run(project_dir: Path, rebuild: list[str] | None = None, resegment: bool = False, no_review: bool = False) -> None:
    prj = Project(project_dir)
    s = prj.settings
    layout_path = prj.path / "layout.csv"
    if not layout_path.exists():
        raise SystemExit(f"{layout_path} not found - run 'kerncount init' first")
    layout = discover.read_table(layout_path)
    recs, warns = discover.records(layout, prj.raw_dir())
    for w in warns:
        log(f"  note: {w}")
    imaged = [r for r in recs if r.path is not None]
    # pixel size only converts units; it must not invalidate cached segmentations
    hashable = {**s, "segmentation": {**s["segmentation"], "pixel_size_um": ""}}
    seg_hash = settings_mod.section_hash(hashable, "segmentation", "channels", "intensity") + f"v{MEASURE_VERSION}"

    # 1. segmentation (cached)
    groups: dict[tuple, list] = {}
    for r in imaged:
        groups.setdefault((r.n, r.treatment), []).append(r)
    log(f"Segmenting {len(imaged)} images in {len(groups)} folders (cached results are reused)...")
    counter = [0]
    for key in sorted(groups):
        _segment_folder(prj, groups[key], seg_hash, resegment, counter, len(imaged))

    # 2. classify + QC + intensity + review images
    log("Classifying nuclei, scoring quality, drawing review images...")
    rows, cells = [], []
    alpha = float(s["review"]["overlay_alpha"])
    legacy_fiji = _fiji_counts(layout, prj.raw_dir())
    for key in sorted(groups):
        g = groups[key]
        loaded = []
        for r in g:
            lab, obj, meta = _cache_paths(prj, r)
            objs = pd.read_csv(obj) if obj.stat().st_size > 1 else pd.DataFrame()
            loaded.append((r, lab, objs, json.loads(meta.read_text())))
        ref = classify.folder_reference([o for _, _, o, _ in loaded])
        sheet = []
        for r, lab, objs, info in loaded:
            cobjs, summ = classify.classify(objs, info, ref, s)
            row = dict(n=r.n, treatment=r.treatment, image=r.image, **summ,
                       diameter_px=round(info.get("diameter_used", np.nan), 1), stardist_scale=round(info.get("scale", 1), 3),
                       fiji_count=legacy_fiji.get((r.n, r.treatment, r.image)), source_file=str(r.path))
            row.update(_intensity(prj, cobjs, info))
            px = s["segmentation"]["pixel_size_um"] or info.get("pixel_size_um")
            if px:
                h, w = tifffile.TiffFile(lab).pages[0].shape[:2]
                row["nuclei_per_mm2"] = summ["auto_count"] / (h * w * float(px) ** 2 / 1e6)
            rows.append(row)
            if len(cobjs):
                cells.append(cobjs.assign(n=r.n, treatment=r.treatment, image=r.image))
            if not no_review:
                p = _review_image(prj, r, lab, cobjs, summ, alpha)
                label = f"image {r.image:02d}   count {summ['auto_count']} ({summ['count_low']}-{summ['count_high']})   QC {summ['qc_level']}"
                sheet.append((p, label, summ["qc_level"]))
        if sheet and not no_review:
            overlay.contact_sheet(sheet, prj.review / _slug(key[0]) / f"{_slug(key[1])}_overview.jpg",
                                  title=f"{key[0]}  |  {key[1]}   (border: green = QC high, orange = medium, red = low)")

    for r in recs:
        if r.path is None:
            rows.append(dict(n=r.n, treatment=r.treatment, image=r.image, auto_count=r.legacy_count, qc_level="legacy",
                             qc_ok="unknown", qc_flags="old Fiji count, raw image missing", fiji_count=r.legacy_count))
    table = pd.DataFrame(rows)
    table["qc_ok"] = table["qc_ok"].map({True: "yes", False: "no"}).fillna(table["qc_ok"])
    table = master.order_rows(table, s["analysis"]["treatment_order"])
    table.to_csv(prj.path / "qc_details.csv", index=False, sep=prj.delim)
    if cells:
        pd.concat(cells, ignore_index=True).to_csv(prj.work / "cells.csv", index=False)
    _fiji_comparison(prj, table)

    # 3. master counts.csv
    counts_path = prj.path / "counts.csv"
    df, mwarn = master.merge(counts_path, table, s["analysis"]["treatment_order"], prj.delim, prj.work / "backups")
    for w in mwarn:
        log(f"  WARNING: {w}")

    # 4. generations
    build_generations(prj, df, table, rebuild or [])
    _summary(prj, table)


def build_generations(prj: Project, df: pd.DataFrame, table: pd.DataFrame, rebuild: list[str]):
    gens = master.chain(df)
    inten = {c: c.replace("intensity_", "").replace("_", " ") + " mean intensity"
             for c in table.columns if c.startswith("intensity_")}
    img_cols = ["n", "treatment", "image", "qc_level"] + list(inten)
    img_table = table[[c for c in img_cols if c in table]]
    root = prj.path / "generations"
    for g in gens:
        folder = root / g["name"]
        meta = folder / "generation.json"
        forced = "all" in rebuild or g["name"] in rebuild
        if folder.exists() and not forced:
            old = json.loads(meta.read_text()) if meta.exists() else {}
            if old.get("fingerprint") != g["fingerprint"]:
                log(f"  WARNING: {g['name']} already exists but counts.csv now gives different values for it "
                    f"(automatic counts or {g['source']} changed). It was NOT overwritten. "
                    f"To rebuild it: kerncount run \"{prj.path}\" --rebuild {g['name']}")
            else:
                log(f"  {g['name']}: up to date")
            continue
        generations.build(prj.path, g, img_table, inten, prj.settings)
        log(f"  {g['name']}: built ({'automatic' if g['source'] == 'auto' else g['source']}, {g['changes']} manual cells)")
    extra = sorted(p.name for p in root.glob("g*") if p.is_dir() and p.name not in {g["name"] for g in gens})
    if extra:
        log(f"  note: {', '.join(extra)} exist but no longer have a matching manual column (left untouched)")


def _review_image(prj: Project, r, lab_path: Path, objs: pd.DataFrame, summ: dict, alpha: float) -> Path:
    out = prj.review / _slug(r.n) / _slug(r.treatment) / f"{r.image:02d}.jpg"
    sig = json.dumps([summ, objs["category"].tolist() if len(objs) else [], alpha], sort_keys=True, default=str)
    sig_path = out.with_suffix(".sig")
    if out.exists() and sig_path.exists() and sig_path.read_text() == sig:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    im = imageio.read(r.path)
    labels = tifffile.imread(lab_path).astype(np.int32)
    title = f"{r.n} | {r.treatment} | image {r.image:02d}    count {summ['auto_count']} (plausible {summ['count_low']}-{summ['count_high']})    QC {summ['qc_level']}"
    sub = summ["qc_flags"] or "no issues flagged"
    img = overlay.render(im.data[int(prj.settings["channels"]["segment"]) - 1], labels, objs, title, sub, alpha)
    img.save(out, quality=int(prj.settings["review"]["jpeg_quality"]))
    sig_path.write_text(sig)
    return out


def _fiji_counts(layout: pd.DataFrame, raw: Path) -> dict:
    out = {}
    for _, row in layout.iterrows():
        f = discover._legacy_csv(raw / row["n_folder"] / row["treatment_folder"])
        if f:
            for k, c in discover._read_legacy(f):
                out[(row["n"], row["treatment"], k)] = c
    return out


def _fiji_comparison(prj: Project, table: pd.DataFrame):
    t = table[(table["qc_level"] != "legacy") & table["fiji_count"].notna()]
    if t.empty:
        return
    t = t.assign(ratio_new_over_fiji=t["auto_count"] / t["fiji_count"].replace(0, np.nan))
    t[["n", "treatment", "image", "fiji_count", "auto_count", "ratio_new_over_fiji", "qc_level"]].to_csv(
        prj.path / "fiji_comparison.csv", index=False, sep=prj.delim)


def _summary(prj: Project, table: pd.DataFrame):
    lv = table["qc_level"].value_counts().to_dict()
    log("")
    log(f"Done. {len(table)} images: " + ", ".join(f"{k} {v}" for k, v in lv.items()))
    log(f"  Edit manual corrections in:  {prj.path / 'counts.csv'}")
    log(f"  Check segmentations in:      {prj.review}")
    log(f"  Results:                     {prj.path / 'generations'}")


def sampling(project_dir: Path, generation: str | None = None, iterations: int = 1000) -> Path:
    from . import sampling as sampling_mod

    prj = Project(project_dir)
    root = prj.path / "generations"
    gens = sorted((p for p in root.glob("g*") if (p / "counts_used.csv").exists()), key=lambda p: int(p.name[1:]))
    if not gens:
        raise SystemExit("no generations yet - run `kerncount run` first")
    g = root / generation if generation else gens[-1]
    used = pd.read_csv(g / "counts_used.csv")
    used = used[~used["excluded"].astype(bool)].assign(value=lambda d: d["count"])
    out = prj.path / "sampling" / g.name
    log(f"Sampling analysis on {g.name} ({iterations} subsamples per design)...")
    summary = sampling_mod.run(used[["n", "treatment", "image", "value"]], prj.settings, out, iterations=iterations)
    log(summary.round(3).to_string(index=False))
    log(f"\nReport: {out / 'REPORT.txt'}")
    return out
