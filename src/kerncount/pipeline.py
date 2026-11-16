"""End-to-end run: discover -> segment (cached) -> classify/QC -> review images -> counts.csv -> generations."""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from . import classify, design, discover, generations, imageio, master, overlay, segment
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
        return settings_mod.channel_names(self.settings, c + 1)[c]

    def seg_channel(self, n_channels: int) -> int:
        """0-based index of the nuclear (segmentation) channel; [channels] segment is a number or a name."""
        try:
            return settings_mod.channel_index(self.settings, self.settings["channels"]["segment"], n_channels)
        except ValueError as e:
            raise SystemExit(f"settings.toml [channels] segment: {e}") from None

    def seg_hash(self) -> str:
        # pixel size only converts units, and settings added later (per-nucleus intensity) are measured outside
        # the detection cache: neither may invalidate cached segmentations (settings.seg_hash_input)
        hashable = settings_mod.seg_hash_input(self.settings)
        return settings_mod.section_hash(hashable, "segmentation", "channels", "intensity") + f"v{MEASURE_VERSION}"

    def records(self, layout: pd.DataFrame):
        return discover.records(layout, self.raw_dir(), self.settings)


# ---------------------------------------------------------------- segmentation cache

def _cache_paths(prj: Project, rec: discover.ImageRecord) -> tuple[Path, Path, Path]:
    d = prj.work / _slug(rec.n) / _slug(rec.treatment)
    stem = f"{rec.image:02d}"
    return d / f"{stem}_labels.tif", d / f"{stem}_objects.csv", d / f"{stem}_meta.json"


def _source_sig(p: Path) -> dict:
    st = p.stat()
    return {"source": str(p), "size": st.st_size, "mtime": st.st_mtime_ns}


def _record_sig(rec) -> dict:
    """Source signature of a field. Single file: exactly as before; interleaved: plus every file of the field."""
    sig = _source_sig(rec.path)
    if rec.paths:
        sig["sources"] = [_source_sig(p) for p in rec.paths]
    return sig


def _cached(prj: Project, rec, seg_hash: str) -> bool:
    lab, obj, meta = _cache_paths(prj, rec)
    if not (lab.exists() and obj.exists() and meta.exists()):
        return False
    m = json.loads(meta.read_text())
    sig = _record_sig(rec)
    ok = m.get("seg_hash") == seg_hash and m.get("size") == sig["size"] and m.get("mtime") == sig["mtime"]
    if rec.paths:
        ok = ok and [(x.get("source"), x.get("size"), x.get("mtime")) for x in m.get("sources", [])] == \
            [(x["source"], x["size"], x["mtime"]) for x in sig["sources"]]
    return ok


def _segment_folder(prj: Project, recs: list, seg_hash: str, force: bool, counter: list[int], total: int):
    todo = [r for r in recs if force or not _cached(prj, r, seg_hash)]
    if not todo:
        counter[0] += len(recs)
        return
    s = prj.settings
    seg_ch = prj.seg_channel(settings_mod.n_channels_expected(s) or imageio.read_record(todo[0]).data.shape[0])
    folder_meta = prj.work / _slug(recs[0].n) / _slug(recs[0].treatment) / "folder.json"
    folder_meta.parent.mkdir(parents=True, exist_ok=True)
    diameter = None
    if str(s["segmentation"]["nucleus_diameter_px"]) == "auto":
        fm = json.loads(folder_meta.read_text()) if folder_meta.exists() else {}
        if fm.get("seg_hash") == seg_hash and len(todo) < len(recs):
            diameter = fm["diameter"]
        else:
            ests = [segment.estimate_diameter(imageio.read_record(r).data[seg_ch]) for r in recs]
            ests = [e for e in ests if e]
            diameter = float(np.median(ests)) if ests else None
            folder_meta.write_text(json.dumps({"seg_hash": seg_hash, "diameter": diameter, "per_image": ests}))
    done_ids = {id(r) for r in todo}
    for r in recs:
        counter[0] += 1
        if id(r) not in done_ids:
            continue
        t0 = time.time()
        try:
            im = imageio.read_record(r)
        except ValueError as e:
            raise SystemExit(str(e)) from None
        if seg_ch >= im.data.shape[0]:
            raise SystemExit(f"{r.path}: has {im.data.shape[0]} channel(s) but settings ask for channel {seg_ch + 1}")
        res = segment.run(im.data, im.max_value, seg_ch, s, diameter)
        lab, obj, meta = _cache_paths(prj, r)
        tifffile.imwrite(lab, res.labels.astype(np.uint16 if res.labels.max() < 65535 else np.int32), compression="zlib")
        res.objects.to_csv(obj, index=False)
        info = {k: v for k, v in res.info.items()}
        info.update(_record_sig(r), seg_hash=seg_hash, pixel_size_um=im.pixel_size_um, n_channels=int(im.data.shape[0]),
                    max_value=float(im.max_value))
        if r.paths:
            info["file_numbers"] = list(r.numbers)
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
    recs, warns = prj.records(layout)
    for w in warns:
        log(f"  note: {w}")
    imaged = [r for r in recs if r.path is not None]
    seg_hash = prj.seg_hash()

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
    rows, cells, fields = [], [], []
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
                       fiji_count=legacy_fiji.get((r.n, r.treatment, r.image)),
                       source_file="; ".join(str(p) for p in r.paths) if r.paths else str(r.path))
            row.update(_intensity(prj, cobjs, info))
            px = settings_mod.pixel_size_um(s) or info.get("pixel_size_um")
            h, w = tifffile.TiffFile(lab).pages[0].shape[:2]
            if s["readouts"]["enabled"]:
                ro, cobjs = classify.compute_readouts(cobjs, (h, w), float(px) if px else None, s["readouts"]["include"])
                row.update(ro)
            if px:
                row["nuclei_per_mm2"] = summ["auto_count"] / (h * w * float(px) ** 2 / 1e6)
            rows.append(row)
            fields.append((r, lab, cobjs, info, row))
            if len(cobjs):
                cells.append(cobjs.assign(n=r.n, treatment=r.treatment, image=r.image))
            if not no_review:
                p = _review_image(prj, r, lab, cobjs, summ, alpha)
                label = f"image {r.image:02d}   count {summ['auto_count']} ({summ['count_low']}-{summ['count_high']})   QC {summ['qc_level']}"
                sheet.append((p, label, summ["qc_level"]))
        if sheet and not no_review:
            overlay.contact_sheet(sheet, prj.review / _slug(key[0]) / f"{_slug(key[1])}_overview.jpg",
                                  title=f"{key[0]}  |  {key[1]}   (border: green = QC high, orange = medium, red = low)")

    # 2b. marker intensity per nucleus, flat field, % positive (docs/intensity.md); adds columns to the rows in place.
    # Single-channel projects have no marker channel: nothing is measured and nothing is written.
    from . import intensity as intensity_mod

    intensity_mod.process(prj, fields, log)

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


def regenerate(project_dir: Path, rebuild: list[str] | None = None) -> None:
    """Rebuild generations from counts.csv + qc_details.csv alone (after manual corrections): no detection, no review images."""
    prj = Project(project_dir)
    counts, details = prj.path / "counts.csv", prj.path / "qc_details.csv"
    if not (counts.exists() and details.exists()):
        raise SystemExit("counts.csv / qc_details.csv not found - run the full analysis first")
    df = discover.read_table(counts)
    table = discover.read_table(details)
    table["image"] = table["image"].astype(int)
    from . import metrics as metrics_mod
    for m in metrics_mod.registry(table.columns, prj.settings).values():
        if m.column in table and m.column != "count":
            table[m.column] = pd.to_numeric(table[m.column], errors="coerce")
    log("Rebuilding generations from counts.csv (detection results are reused unchanged)...")
    build_generations(prj, df, table, rebuild or [])
    log("Done.")


def build_generations(prj: Project, df: pd.DataFrame, table: pd.DataFrame, rebuild: list[str]):
    gens = master.chain(df)
    from . import metrics as metrics_mod
    # every registered metric column (intensity_*, nuclei_per_mm2, readout_*) flows through all statistics
    reg = metrics_mod.registry(table.columns, prj.settings)
    measured = {k: m for k, m in reg.items() if m.column != "count" and m.column in table}
    img_cols = ["n", "treatment", "image", "qc_level", "cell_line"] + [m.column for m in measured.values()]
    img_table = table[[c for c in dict.fromkeys(img_cols) if c in table]]
    design_table = check_design(prj, df)
    # byte-identical raw images (e.g. a copied folder): duplicates.csv, a loud warning, acceptance criterion
    from . import duplicates as duplicates_mod
    try:
        dups = duplicates_mod.check(prj.path, table, log=log)
    except Exception as e:  # noqa: BLE001 - the check must never stop a build
        log(f"  note: duplicate image check failed ({type(e).__name__}: {e})")
        dups = None
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
        generations.build(prj.path, g, img_table, measured, prj.settings, design_table, duplicates=dups)
        log(f"  {g['name']}: built ({'automatic' if g['source'] == 'auto' else g['source']}, {g['changes']} manual cells)")
    extra = sorted(p.name for p in root.glob("g*") if p.is_dir() and p.name not in {g["name"] for g in gens})
    if extra:
        log(f"  note: {', '.join(extra)} exist but no longer have a matching manual column (left untouched)")


def check_design(prj: Project, df: pd.DataFrame) -> pd.DataFrame:
    """Assign well / position / timepoint to every image in counts.csv once, write design_check.csv, log a summary."""
    keys = df[["n", "treatment", "image"]].astype({"image": int}).drop_duplicates()
    try:
        assigned = design.assign(keys, prj.settings, prj.path)
    except ValueError as e:
        raise SystemExit(f"Plate design: {e}") from None
    desc = design.describe(assigned)
    desc.to_csv(prj.path / "design_check.csv", index=False, sep=prj.delim)
    for note in design.spec(prj.settings).notes:
        log(f"  note: {note}")
    log(design.summary_line(desc, prj.settings))
    return assigned


def wells_template(project_dir: Path, overwrite: bool = False) -> Path:
    """Write wells.csv (wells.new.csv if one exists) reproducing the current assignment of the images in counts.csv."""
    prj = Project(project_dir)
    counts = prj.path / "counts.csv"
    if not counts.exists():
        raise SystemExit("counts.csv not found - run the analysis first (the template lists the images it found)")
    df = discover.read_table(counts)
    keys = df[["n", "treatment", "image"]].astype({"image": int}).drop_duplicates()
    path = design.write_wells_template(prj.path, keys, prj.settings, overwrite=overwrite)
    log(f"Wrote {path}. Edit well / position, rename it to wells.csv if needed, and set [design] well_layout = \"file\".")
    return path


def _review_image(prj: Project, r, lab_path: Path, objs: pd.DataFrame, summ: dict, alpha: float) -> Path:
    out = prj.review / _slug(r.n) / _slug(r.treatment) / f"{r.image:02d}.jpg"
    sig = json.dumps([summ, objs["category"].tolist() if len(objs) else [], alpha], sort_keys=True, default=str)
    sig_path = out.with_suffix(".sig")
    if out.exists() and sig_path.exists() and sig_path.read_text() == sig:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    im = imageio.read_record(r)
    labels = tifffile.imread(lab_path).astype(np.int32)
    title = f"{r.n} | {r.treatment} | image {r.image:02d}    count {summ['auto_count']} (plausible {summ['count_low']}-{summ['count_high']})    QC {summ['qc_level']}"
    sub = summ["qc_flags"] or "no issues flagged"
    img = overlay.render(im.data[prj.seg_channel(im.data.shape[0])], labels, objs, title, sub, alpha)
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
    summary = sampling_mod.run(used[["n", "treatment", "image", "value"]], prj.settings, out, iterations=iterations,
                               project_dir=prj.path)
    log(summary.round(3).to_string(index=False))
    from . import planner

    log("Planning the next experiment (planner.json, design_forward.csv, ...)")
    planner.run(used[["n", "treatment", "image", "value"]], prj.settings, out, project_dir=prj.path)
    log(f"\nReport: {out / 'REPORT.txt'}")
    return out
