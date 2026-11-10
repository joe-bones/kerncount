"""Find raw images in a messy folder tree and map them onto clean n / treatment / image names.

Expected raw layout (names may be messy):
    RAW_DIR/<one folder per n>/<one folder per treatment>/**/<number>.tif
Sub-folders called "output" and files containing "mask" are ignored (old Fiji results),
except that an old Nuclei_Count.csv is imported when a treatment has no raw images left.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import pandas as pd

LAYOUT_COLUMNS = ["n_folder", "treatment_folder", "n", "date", "treatment", "include", "images", "note"]
IMAGE_EXT = {".tif", ".tiff"}


@dataclass(frozen=True)
class ImageRecord:
    n: str
    treatment: str
    image: int
    path: Path | None  # None for legacy (count-only) entries
    legacy_count: int | None = None


# ---------------------------------------------------------------- naming rules

def parse_date(name: str, order: str = "DMY") -> date | None:
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

def scan(raw_dir: Path, date_order: str = "DMY") -> pd.DataFrame:
    raw_dir = Path(raw_dir)
    n_dirs = [d for d in sorted(raw_dir.iterdir()) if d.is_dir() and not d.name.startswith(".")]
    n_info = []
    for d in n_dirs:
        n_info.append((parse_date(d.name, date_order), d))
    n_info.sort(key=lambda t: (t[0] is None, t[0] or date.min, t[1].name))

    rows = []
    for i, (dt, d) in enumerate(n_info, start=1):
        n_label = f"n{i}" + (f"_{dt.isoformat()}" if dt else "")
        n_note = "72h" if re.search(r"\b72\s*h", d.name, re.I) else ""
        for t in sorted(p for p in d.iterdir() if p.is_dir() and not p.name.startswith(".")):
            name, note = canonical_treatment(t.name)
            imgs = _raw_images(t)
            numbered = [p for p in imgs if image_number(p) is not None]
            notes = [x for x in (n_note, note) if x]
            skipped = [p.name for p in imgs if image_number(p) is None]
            if skipped:
                notes.append("skipped non-numbered: " + "; ".join(skipped))
            if not numbered:
                if _legacy_csv(t):
                    notes.append("NO RAW IMAGES - old Fiji counts imported (legacy)")
                else:
                    notes.append("no images found")
            rows.append(dict(n_folder=d.name, treatment_folder=t.name, n=n_label,
                             date=dt.isoformat() if dt else "", treatment=name, include="yes",
                             images=len(numbered), note=" | ".join(notes)))
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


def records(layout: pd.DataFrame, raw_dir: Path) -> tuple[list[ImageRecord], list[str]]:
    raw_dir = Path(raw_dir)
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
