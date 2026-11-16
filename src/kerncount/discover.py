"""Find raw images in a messy folder tree and map them onto clean n / treatment / image names.

Expected raw layout (names may be messy):
    RAW_DIR/<one folder per n>/<one folder per treatment>/**/<number>.tif
Sub-folders called "output" and files containing "mask" are ignored (old Fiji results),
except that an old Nuclei_Count.csv is imported when a treatment has no raw images left.

Channels as separate files ([channels] files = "interleaved"): the numbered files of a treatment folder are taken
in numeric order and grouped into fields of k files (k = len(file_order), e.g. Ki67, DAPI, Ki67, DAPI, ...).
Gaps in the numbering are ignored; the field number is the group index 1..m. A leftover file (count not a
multiple of k) is reported and skipped.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

LAYOUT_COLUMNS = ["n_folder", "treatment_folder", "n", "date", "timepoint", "treatment", "include", "images", "note"]
IMAGE_EXT = {".tif", ".tiff"}


@dataclass(frozen=True)
class ImageRecord:
    n: str
    treatment: str
    image: int  # field number (interleaved: index of the file group, 1..m)
    path: Path | None  # None for legacy (count-only) entries; interleaved: the first file of the field
    legacy_count: int | None = None
    paths: tuple = ()  # interleaved: every file of the field, in file_order; () for single-file images
    numbers: tuple = ()  # interleaved: the file numbers of `paths`

    @property
    def sources(self) -> tuple:
        """Every file this field is read from."""
        return self.paths if self.paths else ((self.path,) if self.path is not None else ())


def group_files(numbered: list[tuple[int, Path]], k: int) -> tuple[list[tuple[int, tuple, tuple]], list[Path]]:
    """Group numbered files (any order) into fields of k consecutive files in numeric order.
    Returns ([(field, paths, numbers)], leftover paths). Gaps in the numbering are ignored."""
    items = sorted(numbered, key=lambda t: (t[0], str(t[1])))
    k = max(int(k), 1)
    m = len(items) // k
    fields = [(i + 1, tuple(p for _, p in items[i * k:(i + 1) * k]), tuple(n for n, _ in items[i * k:(i + 1) * k]))
              for i in range(m)]
    return fields, [p for _, p in items[m * k:]]


# ---------------------------------------------------------------- naming rules

def parse_date(name: str, order: str = "DMY", year: int | None = None) -> date | None:
    """Date in a folder name ("05.06.26", "050626"). With `year` given, a bare 4-digit day+month ("0611" = 6 Nov
    for DMY) is also accepted (`scan` takes the year from the file dates and says so in the note)."""
    d = _parse_full_date(name, order)
    if d is None and year:
        m4 = re.search(r"(?<!\d)(\d{2})(\d{2})(?!\d)", name)
        if m4:
            a, b = (int(x) for x in m4.groups())
            day, mo = (b, a) if order.upper() in ("MDY", "YMD") else (a, b)
            try:
                d = date(year, mo, day)
            except ValueError:
                d = None
    return d


def _parse_full_date(name: str, order: str = "DMY") -> date | None:
    m = re.search(r"(\d{1,2})[.\-_/](\d{1,2})[.\-_/](\d{2,4})", name)
    if not m:
        m6 = re.search(r"(?<!\d)(\d{2})(\d{2})(\d{2})(?!\d)", name)
        if not m6:
            return None
        m = m6
    a, b, c = (int(x) for x in m.groups())
    if order.upper() == "DMY":
        d, mo, y = a, b, c
    elif order.upper() == "MDY":
        mo, d, y = a, b, c
    else:  # YMD
        y, mo, d = a, b, c
    if y < 100:
        y += 2000
    try:
        return date(y, mo, d)
    except ValueError:
        return None


_TP_HOURS = re.compile(r"(?<![A-Za-z0-9.])(\d+(?:[.,]\d+)?)\s*(?:h|hr|hrs|hour|hours)(?![a-z])", re.I)
_TP_DAYS = re.compile(r"(?<![A-Za-z0-9.])(?:days?\s*(\d+)|d(\d+)|(\d+)\s*(?:days?|d))(?![a-z\d])", re.I)


def parse_timepoint(name: str) -> str:
    """Timepoint in an n-folder name, normalised: '72h', '48 h', '24hr', '72 hrs' -> '72h' etc.;
    'day 3', 'd3', '3 days' -> 'd3'. '' if there is none."""
    m = _TP_HOURS.search(name)
    if m:
        return m.group(1).replace(",", ".") + "h"
    m = _TP_DAYS.search(name)
    if m:
        return "d" + next(g for g in m.groups() if g)
    return ""


def canonical_treatment(folder: str) -> tuple[str, str]:
    """Return (canonical name, note) for a messy treatment folder name."""
    s = folder.strip()
    low = s.lower()
    if "lyo" in low:
        return "LyoVec", ""
    if re.search(r"cl\s*-?\s*0?75", low):
        return "CL075", ""
    m = re.search(r"(\d+[a-z]?)\s*[-\s]\s*(\d)\s*p\b", low)
    if m:
        return f"miR-{m.group(1)}-{m.group(2)}p", ""
    m = re.fullmatch(r"(\d+[a-z]?)", low)
    if m:
        return f"miR-{m.group(1)}", "arm (-3p/-5p) inferred from other n's if possible"
    if re.search(r"cntrl|ctrl|contr", low):
        extra = re.sub(r"\d{4,}|cntrl|ctrl|control|contr", "", low).strip(" _-.,")
        note = f"original label '{s}'" if extra else ""
        return "Control", note
    return re.sub(r"\s+", "_", s), "unrecognised name: please check"


def image_number(path: Path) -> int | None:
    m = re.fullmatch(r"[^\w]*0*(\d+)", path.stem)
    return int(m.group(1)) if m else None


def _raw_images(folder: Path) -> list[Path]:
    out = []
    for p in folder.rglob("*"):
        if p.suffix.lower() not in IMAGE_EXT or not p.is_file():
            continue
        rel = p.relative_to(folder)
        if any(part.lower() == "output" for part in rel.parts[:-1]):
            continue
        if "mask" in p.stem.lower() or p.name.startswith("."):
            continue
        out.append(p)
    return sorted(out)


def _legacy_csv(folder: Path) -> Path | None:
    hits = sorted(folder.rglob("Nuclei_Count.csv"))
    return hits[0] if hits else None


# ---------------------------------------------------------------- layout

def _file_year(folder: Path) -> int | None:
    """Most common year among the image files' modification dates (to complete a day+month folder date)."""
    years = []
    for p in folder.rglob("*"):
        if p.suffix.lower() in IMAGE_EXT and p.is_file():
            years.append(date.fromtimestamp(p.stat().st_mtime).year)
            if len(years) >= 50:
                break
    return max(set(years), key=years.count) if years else None


def scan(raw_dir: Path, date_order: str = "DMY", files: str = "single", n_channels: int = 1) -> pd.DataFrame:
    """One row per n x treatment folder. files="interleaved": `images` counts fields of n_channels files each."""
    raw_dir = Path(raw_dir)
    inter = str(files).lower() == "interleaved"
    n_dirs = [d for d in sorted(raw_dir.iterdir()) if d.is_dir() and not d.name.startswith(".")]
    n_info, date_notes = [], {}
    for d in n_dirs:
        dt = parse_date(d.name, date_order)
        if dt is None:
            yr = _file_year(d)
            dt = parse_date(d.name, date_order, year=yr) if yr else None
            date_notes[d.name] = (f"date {dt.isoformat()} = day/month from the folder name + year of the file dates: "
                                  "check it (edit n and date if wrong)") if dt else \
                "no date in the folder name: fill in date (and the n label) here"
        n_info.append((dt, d))
    n_info.sort(key=lambda t: (t[0] is None, t[0] or date.min, t[1].name))

    rows = []
    for i, (dt, d) in enumerate(n_info, start=1):
        n_label = f"n{i}" + (f"_{dt.isoformat()}" if dt else "")
        tp = parse_timepoint(d.name)
        n_note = " | ".join(x for x in (tp, date_notes.get(d.name, "")) if x)
        for t in sorted(p for p in d.iterdir() if p.is_dir() and not p.name.startswith(".")):
            name, note = canonical_treatment(t.name)
            imgs = _raw_images(t)
            numbered = [p for p in imgs if image_number(p) is not None]
            notes = [x for x in (n_note, note) if x]
            skipped = [p.name for p in imgs if image_number(p) is None]
            if skipped:
                notes.append("skipped non-numbered: " + "; ".join(skipped))
            n_images = len(numbered)
            if inter and numbered:
                fields, left = group_files([(image_number(p), p) for p in numbered], n_channels)
                n_images = len(fields)
                if left:
                    notes.append(f"{len(numbered)} files are not a multiple of {n_channels} channels: leftover "
                                 + "; ".join(p.name for p in left) + " skipped")
            if not numbered:
                if _legacy_csv(t):
                    notes.append("NO RAW IMAGES - old Fiji counts imported (legacy)")
                else:
                    notes.append("no images found")
            rows.append(dict(n_folder=d.name, treatment_folder=t.name, n=n_label,
                             date=dt.isoformat() if dt else "", timepoint=tp, treatment=name, include="yes",
                             images=n_images, note=" | ".join(notes)))
    df = pd.DataFrame(rows, columns=LAYOUT_COLUMNS)

    # Bare numbers ("21") borrow their arm from other n's ("miR-21-5p").
    full = df["treatment"][df["treatment"].str.fullmatch(r"miR-\d+[a-z]?-\dp")].unique()
    for i, t in df["treatment"].items():
        if re.fullmatch(r"miR-\d+[a-z]?", t):
            cands = [f for f in full if f.startswith(t + "-")]
            if len(cands) == 1:
                df.at[i, "treatment"] = cands[0]
                df.at[i, "note"] = re.sub(r"arm \(-3p/-5p\) inferred.*?possible", "arm inferred from other n's", df.at[i, "note"])
    return df


def write_layout(df: pd.DataFrame, path: Path, delimiter: str = ",") -> None:
    df.to_csv(path, index=False, sep=delimiter)


def read_table(path: Path) -> pd.DataFrame:
    """Read a CSV that may have been re-saved by Excel with ';' and a BOM."""
    text = Path(path).read_text(encoding="utf-8-sig")
    try:
        delim = csv.Sniffer().sniff(text.splitlines()[0], delimiters=",;\t").delimiter
    except csv.Error:
        delim = ","
    return pd.read_csv(path, sep=delim, dtype=str, keep_default_na=False, encoding="utf-8-sig")


def records(layout: pd.DataFrame, raw_dir: Path, settings: dict | None = None) -> tuple[list[ImageRecord], list[str]]:
    """Image records of the included layout rows. With [channels] files = "interleaved" in `settings`, every record
    is one field of len(file_order) consecutive files (see `group_files`)."""
    raw_dir = Path(raw_dir)
    per_field = 0
    if settings is not None and str(settings.get("channels", {}).get("files", "single")).lower() == "interleaved":
        per_field = max(len(settings["channels"].get("file_order") or []), 1)
    recs: list[ImageRecord] = []
    warnings: list[str] = []
    seen: dict[tuple, Path] = {}
    for _, row in layout.iterrows():
        if str(row["include"]).strip().lower() not in ("yes", "y", "true", "1"):
            continue
        folder = raw_dir / row["n_folder"] / row["treatment_folder"]
        if not folder.is_dir():
            warnings.append(f"missing folder: {folder}")
            continue
        imgs = [(image_number(p), p) for p in _raw_images(folder)]
        imgs = [(k, p) for k, p in imgs if k is not None]
        if imgs and per_field:
            fields, left = group_files(imgs, per_field)
            if left:
                warnings.append(f"{row['n']} / {row['treatment']}: {len(imgs)} files are not a multiple of {per_field} "
                                f"channels; leftover file(s) skipped: {', '.join(p.name for p in left)}")
            for f, paths, nums in fields:
                recs.append(ImageRecord(row["n"], row["treatment"], f, paths[0], paths=paths, numbers=nums))
            continue
        if imgs:
            for k, p in sorted(imgs):
                key = (row["n"], row["treatment"], k)
                if key in seen:
                    warnings.append(f"duplicate image number {k}: {seen[key]} and {p} (second ignored)")
                    continue
                seen[key] = p
                recs.append(ImageRecord(row["n"], row["treatment"], k, p))
            continue
        legacy = _legacy_csv(folder)
        if legacy:
            for k, c in _read_legacy(legacy):
                key = (row["n"], row["treatment"], k)
                if key not in seen:
                    seen[key] = legacy
                    recs.append(ImageRecord(row["n"], row["treatment"], k, None, c))
            warnings.append(f"{row['n']} / {row['treatment']}: no raw images, using old Fiji counts")
    return recs, warnings


def _read_legacy(path: Path) -> list[tuple[int, int]]:
    out = []
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 2 or not parts[1].isdigit():
            continue
        m = re.fullmatch(r"[^\w]*0*(\d+)", parts[0])
        if m:
            out.append((int(m.group(1)), int(parts[1])))
    return sorted(out)


# ---------------------------------------------------------------- interleaved channel files: detection

def detect_interleaved(raw_dir: Path, max_folders: int = 4, max_files: int = 8) -> dict:
    """Guess whether channels were saved as alternating single-channel files (e.g. Ki67, DAPI, Ki67, DAPI, ...).

    Looks at the first files of a few treatment folders: every file must hold one channel (greyscale, or RGB with
    three identical channels) and the brightness must alternate (one position of each pair consistently brighter).
    Returns {"interleaved": bool, "bright_position": 0 | 1 | None, "folders": int, "detail": str}; it only
    suggests: the user decides via [channels] files."""
    from . import imageio

    raw_dir = Path(raw_dir)
    folders = []
    for n in sorted(d for d in raw_dir.iterdir() if d.is_dir() and not d.name.startswith(".")):
        folders += sorted(t for t in n.iterdir() if t.is_dir() and not t.name.startswith("."))
    looked, votes, details = 0, [], []
    for t in folders[:max_folders]:
        files = sorted(((image_number(p), p) for p in _raw_images(t)), key=lambda x: (x[0] is None, x[0] or 0))
        files = [p for k, p in files if k is not None][:max_files]
        if len(files) < 4:
            continue
        looked += 1
        levels, single = [], True
        for p in files:
            try:
                im = imageio.read(p)
            except Exception:
                single = False
                break
            if im.data.shape[0] != 1:
                single = False
                break
            levels.append(float(np.percentile(im.data[0], 99)))
        if not single:
            details.append(f"{t.name}: multi-channel files")
            votes.append(None)
            continue
        a, b = np.array(levels[0::2]), np.array(levels[1::2])
        m = min(len(a), len(b))
        a, b = a[:m], b[:m]
        if (b > 1.5 * a).all():
            votes.append(1)
        elif (a > 1.5 * b).all():
            votes.append(0)
        else:
            votes.append(None)
        details.append(f"{t.name}: 99th percentile by file " + ", ".join(f"{v:.0f}" for v in levels))
    ok = looked > 0 and all(v is not None for v in votes) and len(set(votes)) == 1
    return dict(interleaved=bool(ok), bright_position=votes[0] if ok else None, folders=looked,
                detail="; ".join(details))
