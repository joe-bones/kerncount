"""Raw experiment folder -> layout.csv rows (one row per tumour).

    <raw>/<imaging session>/<date> <tumour line> <±clodro> n<k> <A|B>.ome.tiff     original z-stack = one tumour
                            ... <A|B>f.tif, Statistics for ... <A|B>f.csv|.txt      Fiji legacy outputs (not tumours)

Factor values come from filename tokens or from the session-folder name, as declared per factor (config.Factor).
This is the happy-path parser; unreadable names are excluded with a note (full validation: ticket 03).
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .. import discover
from .config import Factor, factors as declared_factors

ORIGINAL_SUFFIXES = (".ome.tiff", ".ome.tif")
LEGACY_STATS_PREFIX = "Statistics for "


def is_original(name: str) -> bool:
    return name.lower().endswith(ORIGINAL_SUFFIXES)


def is_legacy(name: str) -> bool:
    """Fiji outputs of `volume calculator tumor.ijm`: binary mask `<stem>f.tif` and `Statistics for <stem>f.csv|txt`."""
    low = name.lower()
    return name.startswith(LEGACY_STATS_PREFIX) or (low.endswith((".tif", ".tiff")) and not is_original(name)
                                                    and Path(name).stem.endswith("f"))


def _stem(name: str) -> str:
    low = name.lower()
    for suf in ORIGINAL_SUFFIXES:
        if low.endswith(suf):
            return name[: -len(suf)]
    return name


def _tokens(text: str) -> list[str]:
    return [t for t in re.split(r"[\s_]+", text.strip()) if t]


def _read_factor(f: Factor, tokens: list[str]) -> tuple[str, str, str]:
    """(level, token, note): the first token naming a level of f; note when it was read through an alias."""
    for t in tokens:
        hit = f.match(t)
        if hit:
            level, via_alias = hit
            return level, t, (f"{f.name}: '{t}' read as '{level}'" if via_alias else "")
    return "", "", ""


def tumour_id(values: list[str], section: str, hemisphere: str) -> str:
    """Unique across the experiment: every factor level + section label + hemisphere (CONTEXT.md)."""
    return " / ".join([*values, f"{section}{hemisphere}"])


def _sessions(raw: Path) -> list[tuple[str, Path]]:
    subs = sorted(d for d in raw.iterdir() if d.is_dir() and not d.name.startswith((".", "_")))
    out = [(d.name, d) for d in subs]
    if any(p.is_file() and is_original(p.name) for p in raw.iterdir()):
        out.insert(0, (raw.name, raw))  # flat files directly in the experiment folder: one unnamed session
    return out


def scan(raw: Path, settings: dict) -> pd.DataFrame:
    raw = Path(raw)
    facs = declared_factors(settings)
    order = settings.get("project", {}).get("date_order", "DMY")
    rows = []
    for session, folder in _sessions(raw):
        names = sorted(p.name for p in folder.iterdir() if p.is_file() and not p.name.startswith("."))
        session_tokens = _tokens(session)
        rel = "" if folder == raw else f"{session}/"
        for name in names:
            if not is_original(name):
                continue  # Fiji legacy files and anything else are not tumours
            stem = _stem(name)
            tokens = _tokens(stem)
            values, notes, missing, used = [], [], [], set()
            for f in facs:
                level, token, note = _read_factor(f, session_tokens if f.source == "session" else tokens)
                values.append(level)
                if f.source == "filename":
                    used.add(token)
                if note and f.source == "filename":
                    notes.append(note)
                if not level:
                    missing.append(f.name)
            rest = [t for t in tokens if t not in used]  # a factor level (tumour line "N4") is not a section label
            section = next((t.lower() for t in rest if re.fullmatch(r"[nN]\d+", t)), "")
            hemisphere = next((t.upper() for t in rest if re.fullmatch(r"[ABab]", t)), "")
            if not section:
                missing.append("section label (n<k>)")
            if not hemisphere:
                missing.append("hemisphere (A/B)")
            d = discover.parse_date(stem, order)
            legacy_mask = f"{stem}f.tif" if f"{stem}f.tif" in names else ""
            legacy_stats = next((f"{LEGACY_STATS_PREFIX}{stem}f{ext}" for ext in (".csv", ".txt")
                                 if f"{LEGACY_STATS_PREFIX}{stem}f{ext}" in names), "")
            if missing:
                notes.insert(0, "could not read " + ", ".join(missing))
            rows.append({
                "tumour_id": tumour_id([v or "?" for v in values], section or "n?", hemisphere or "?"),
                "imaging_session": session,
                "date": d.isoformat() if d else "",
                **{f.name: v for f, v in zip(facs, values)},
                "section": section, "hemisphere": hemisphere, "tumour": f"{section}{hemisphere}",
                "file": rel + name,
                "legacy_mask": rel + legacy_mask if legacy_mask else "",
                "legacy_stats": rel + legacy_stats if legacy_stats else "",
                "include": "no" if missing else "yes",
                "note": "; ".join(notes),
            })
    cols = ["tumour_id", "imaging_session", "date", *[f.name for f in facs], "section", "hemisphere", "tumour",
            "file", "legacy_mask", "legacy_stats", "include", "note"]
    return pd.DataFrame(rows, columns=cols)
