"""Raw experiment folder -> layout.csv rows (one row per tumour), and the check that blocks a run on a bad layout.

    <raw>/<imaging session>/<date> <tumour line> <±clodro> n<k> <A|B>.ome.tiff     original z-stack = one tumour
                            ... <A|B>f.tif                                         Fiji legacy binary mask
                            Statistics for ... <A|B>f.csv|.txt, ... <A|B>f_3Dresults.csv, ... <A|B>.txt
                                                                                    Fiji legacy statistics

Factor values come from filename tokens or from the session-folder name, as declared per factor (config.Factor).

`scan` reads the folder: known naming variants (`-clodron`, `.ome.tiff.tif`, upper/lower case, `n3A`, `- clodro`,
underscores, extra spaces, typographic minus) are normalised and noted; names it cannot read keep empty cells and a
note. `validate` checks a (possibly hand-edited) layout.csv against the folder: values typed into layout.csv override
what the names say (noted per row), `include = no` + `note` excludes a tumour with that reason, and every problem that
would turn a file into a wrong datapoint (unreadable values, duplicate tumour IDs, a session whose genotype is unknown,
missing files) is collected. `check` raises them all at once as `LayoutBlocked`.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .. import discover
from .config import Factor, factors as declared_factors

ORIGINAL_SUFFIXES = (".ome.tiff.tif", ".ome.tiff", ".ome.tif")  # longest first; .ome.tiff.tif = re-saved by Fiji
STANDARD_SUFFIXES = (".ome.tiff", ".ome.tif")
LEGACY_STATS_PREFIX = "Statistics for "
YES = ("yes", "y", "true", "1", "ja", "j")
NO = ("no", "n", "false", "0", "nein")
EXCLUDED = "excluded_tumours.csv"
DASHES = str.maketrans({"−": "-", "–": "-", "—": "-", "‐": "-", "‑": "-"})


def _split_original(name: str) -> tuple[str, str] | None:
    """(stem, suffix as written) when `name` is an original z-stack, else None (case-insensitive)."""
    low = name.lower()
    for suf in ORIGINAL_SUFFIXES:
        if low.endswith(suf):
            return name[: -len(suf)], name[-len(suf):]
    return None


def is_original(name: str) -> bool:
    return _split_original(name) is not None


def legacy_candidates(stem: str) -> tuple[list[str], list[str]]:
    """File names the Fiji macro (`volume calculator tumor.ijm`) wrote for the original `<stem>.ome.tiff`:
    (binary masks, statistics tables), most standard first.
    The `f` is the macro stripping `.ome.tif` off `.ome.tiff`."""
    masks = [f"{stem}f.tif", f"{stem}f.tiff"]
    stats = [f"{LEGACY_STATS_PREFIX}{stem}f.csv", f"{LEGACY_STATS_PREFIX}{stem}f.txt", f"{stem}f_3Dresults.csv",
             f"{stem}.txt", f"{stem}f.txt"]
    return masks, stats


def is_legacy(name: str) -> bool:
    """Fiji outputs by their name pattern alone (`scan` additionally matches them to a tumour)."""
    low = name.lower()
    return (name.startswith(LEGACY_STATS_PREFIX) or low.endswith("_3dresults.csv")
            or (low.endswith((".tif", ".tiff")) and not is_original(name) and Path(name).stem.endswith("f")))


# ------------------------------------------------------------------------------------------------ name reading


def _tokens(text: str) -> list[str]:
    return [t for t in re.split(r"[\s_]+", text.strip()) if t]


def _normalise_stem(stem: str) -> tuple[list[str], list[str]]:
    """Tokens of a file stem with spacing / sign / section-label variants undone, plus a note per variant."""
    notes = []
    text = stem.translate(DASHES)
    if text != stem:
        notes.append("typographic dash read as '-'")
    if "_" in text or re.search(r"\s{2,}", text) or text != text.strip():
        notes.append("spacing normalised")
    text = re.sub(r"[\s_]+", " ", text).strip()
    joined = re.sub(r"(?<![^\s])([+-])\s+(?=\S)", r"\1", text)  # "- clodro" -> "-clodro"
    if joined != text:
        notes.append("sign joined to its word ('- clodro' -> '-clodro')")
    tokens = []
    for t in _tokens(joined):
        m = re.fullmatch(r"([nN])(\d+)([ABab])", t)  # "n3A" -> section n3 + hemisphere A
        if m:
            tokens += [m.group(1) + m.group(2), m.group(3)]
            notes.append(f"'{t}' read as section {m.group(1).lower()}{m.group(2)} + hemisphere {m.group(3).upper()}")
        else:
            tokens.append(t)
    return tokens, notes


def _read_factor(f: Factor, tokens: list[str]) -> tuple[str, str, str, str]:
    """(level, token, note, problem) for factor f from tokens: a note when read through an alias or in other case,
    a problem when tokens name two different levels."""
    hits: list[tuple[str, str, bool]] = []
    for t in tokens:
        hit = f.match(t)
        if hit and hit[0] not in [h[0] for h in hits]:
            hits.append((hit[0], t, hit[1]))
    if not hits:
        return "", "", "", ""
    if len(hits) > 1:  # "HERV ... n4": the exact spelling wins over a case variant (section n4 vs tumour line N4)
        exact = [h for h in hits if not h[2] and h[1] == h[0]]
        hits = exact if len(exact) == 1 else hits
    if len(hits) > 1:
        return "", "", "", f"{f.name} ambiguous ({' or '.join(repr(h[1]) for h in hits)})"
    level, token, via_alias = hits[0]
    note = f"{f.name}: '{token}' read as '{level}'" if via_alias or token != level else ""
    return level, token, note, ""


def tumour_id(values: list[str], section: str, hemisphere: str) -> str:
    """Unique across the experiment: every factor level + section label + hemisphere (CONTEXT.md)."""
    return " / ".join([*values, f"{section}{hemisphere}"])


def _sessions(raw: Path) -> list[tuple[str, Path]]:
    subs = sorted(d for d in raw.iterdir() if d.is_dir() and not d.name.startswith((".", "_")))
    out = [(d.name, d) for d in subs]
    if any(p.is_file() and is_original(p.name) for p in raw.iterdir()):
        out.insert(0, (raw.name, raw))  # flat files directly in the experiment folder: one unnamed session
    return out


def _columns(facs: list[Factor]) -> list[str]:
    return ["tumour_id", "imaging_session", "date", *[f.name for f in facs], "section", "hemisphere", "tumour",
            "file", "legacy_mask", "legacy_stats", "include", "note"]


def _row(name: str, session: str, session_tokens: list[str], names: set[str], rel: str, facs: list[Factor],
         order: str) -> dict:
    split = _split_original(name)
    assert split is not None
    stem, suffix = split
    tokens, notes = _normalise_stem(stem)
    if suffix not in STANDARD_SUFFIXES:
        notes.insert(0, f"file ending '{suffix}' read as an original stack")
    values, missing, used = [], [], set()
    for f in facs:
        level, token, note, problem = _read_factor(f, session_tokens if f.source == "session" else tokens)
        values.append(level)
        if f.source == "filename":
            used.add(token)
            if note:
                notes.append(note)
        if problem:
            missing.append(problem)
        elif not level:
            missing.append(f"{f.name} (from the {'session folder' if f.source == 'session' else 'file'} name)")
    rest = [t for t in tokens if t not in used]  # a factor level (tumour line "N4") is not a section label
    sections = [t for t in rest if re.fullmatch(r"[nN]\d+", t)]
    hemis = [t for t in rest if re.fullmatch(r"[ABab]", t)]
    section = sections[0].lower() if len(sections) == 1 else ""
    hemisphere = hemis[0].upper() if len(hemis) == 1 else ""
    if len(sections) == 1 and sections[0] != section:
        notes.append(f"section label '{sections[0]}' read as '{section}'")
    if len(hemis) == 1 and hemis[0] != hemisphere:
        notes.append(f"hemisphere '{hemis[0]}' read as '{hemisphere}'")
    if not section:
        missing.append("section label (n<k>)" + (f" ambiguous ({', '.join(sections)})" if sections else ""))
    if not hemisphere:
        missing.append("hemisphere (A/B)" + (f" ambiguous ({', '.join(hemis)})" if hemis else ""))
    d = discover.parse_date(stem, order)
    masks, stats = legacy_candidates(stem)
    lower = {n.lower(): n for n in names}
    legacy_mask = next((lower[m.lower()] for m in masks if m.lower() in lower), "")
    legacy_stats = next((lower[s.lower()] for s in stats if s.lower() in lower), "")
    if missing:
        notes.insert(0, "could not read " + ", ".join(missing) + ": fill the empty cells or set include = no")
    return {
        "tumour_id": tumour_id(values, section, hemisphere) if not missing else "",
        "imaging_session": session,
        "date": d.isoformat() if d else "",
        **{f.name: v for f, v in zip(facs, values)},
        "section": section, "hemisphere": hemisphere, "tumour": f"{section}{hemisphere}" if section and hemisphere
        else "",
        "file": rel + name,
        "legacy_mask": rel + legacy_mask if legacy_mask else "",
        "legacy_stats": rel + legacy_stats if legacy_stats else "",
        "include": "yes",  # unreadable rows stay included: the run blocks until they are filled in or excluded
        "note": "; ".join(notes),
    }


def scan(raw: Path, settings: dict) -> pd.DataFrame:
    raw = Path(raw)
    facs = declared_factors(settings)
    order = settings.get("project", {}).get("date_order", "DMY")
    rows = []
    for session, folder in _sessions(raw):
        names = sorted(p.name for p in folder.iterdir() if p.is_file() and not p.name.startswith("."))
        rel = "" if folder == raw else f"{session}/"
        session_tokens, _ = _normalise_stem(session)
        rows += [_row(n, session, session_tokens, set(names), rel, facs, order) for n in names if is_original(n)]
    return pd.DataFrame(rows, columns=_columns(facs))


def unmatched_files(raw: Path, layout: pd.DataFrame) -> list[str]:
    """Files in the raw folder that are neither a tumour nor a legacy output matched to one (shown, never used)."""
    raw = Path(raw)
    known = set(layout["file"]) | set(layout["legacy_mask"]) | set(layout["legacy_stats"])
    out = []
    for session, folder in _sessions(raw):
        rel = "" if folder == raw else f"{session}/"
        out += [rel + p.name for p in sorted(folder.iterdir())
                if p.is_file() and not p.name.startswith(".") and rel + p.name not in known]
    return out


# ------------------------------------------------------------------------------------------------ validation


class LayoutBlocked(SystemExit):
    """layout.csv cannot be run: every problem, listed at once."""

    def __init__(self, problems: list[str], where: str = "layout.csv"):
        self.problems = list(problems)
        lines = [f"{where} has {len(problems)} problem{'s' if len(problems) != 1 else ''} that must be fixed before "
                 "the run (no file is measured until every included file maps to exactly one unique tumour ID):"]
        lines += [f"  - {p}" for p in problems]
        lines.append("Fix the values in layout.csv, or set include = no and give the reason in note, then run again.")
        super().__init__("\n".join(lines))


@dataclass
class Validation:
    layout: pd.DataFrame  # every row; include normalised to yes/no, tumour IDs rebuilt from the (overridden) values,
    # exclusion_reason = the user's part of note on excluded rows (scan notes removed)
    problems: list[str] = field(default_factory=list)  # each blocks the run
    warnings: list[str] = field(default_factory=list)  # shown, do not block

    @property
    def included(self) -> pd.DataFrame:
        return self.layout[self.layout["include"] == "yes"].reset_index(drop=True)

    @property
    def excluded(self) -> pd.DataFrame:
        """Excluded tumours with `exclusion_reason` (the row's note), for flow tables and reports."""
        ex = self.layout[self.layout["include"] == "no"]
        keep = [c for c in self.layout.columns if c not in ("include", "note", "legacy_mask", "legacy_stats")]
        return ex[keep].reset_index(drop=True)


def _cell(row: pd.Series, col: str) -> str:
    v = row.get(col, "")
    return "" if v is None or (isinstance(v, float) and pd.isna(v)) else str(v).strip()


def _add_note(note: str, extra: str) -> str:
    parts = [p.strip() for p in note.split(";") if p.strip()]
    return "; ".join(parts + ([extra] if extra not in parts else []))


NO_REASON = "excluded in layout.csv (no reason given)"


def _reason(note: str, parsed: pd.Series | None) -> str:
    """The user's exclusion reason: the row's note without what `scan` wrote there itself."""
    scan_notes = {p.strip() for p in str(parsed["note"]).split(";")} if parsed is not None else set()
    own = [p.strip() for p in note.split(";") if p.strip() and p.strip() not in scan_notes]
    return "; ".join(own) or NO_REASON


def validate(layout: pd.DataFrame, raw: Path | None, settings: dict) -> Validation:
    """Check a layout (as edited by the user) against the raw folder; see the module docstring.
    `raw` = None skips the comparison with the folder (no override notes, no missing-file check)."""
    facs = declared_factors(settings)
    lay = layout.copy().fillna("").astype(str).reset_index(drop=True)
    for c in _columns(facs):
        if c not in lay.columns:
            lay[c] = ""
    lay["exclusion_reason"] = ""
    problems: list[str] = []
    warnings: list[str] = []
    scanned = scan(raw, settings).set_index("file") if raw is not None else None
    session_missing: dict[tuple[str, str], list[str]] = {}

    for i, row in lay.iterrows():
        file = _cell(row, "file") or f"row {int(str(i)) + 2}"
        parsed = scanned.loc[file] if scanned is not None and file in scanned.index else None
        inc = _cell(row, "include").lower()
        if inc in NO:
            lay.at[i, "include"] = "no"
            lay.at[i, "exclusion_reason"] = _reason(_cell(row, "note"), parsed)
            continue
        if inc not in YES:
            problems.append(f"{file}: include = {_cell(row, 'include')!r} must be yes or no")
            continue
        lay.at[i, "include"] = "yes"
        note = _cell(row, "note")
        if scanned is not None and parsed is None:
            problems.append(f"{file}: not found in the raw folder")
        bad = []
        values = []
        for f in facs:
            v = _cell(row, f.name)
            hit = f.match(v) if v else None
            if not v:
                if f.source == "session" and _cell(row, "imaging_session"):
                    session_missing.setdefault((_cell(row, "imaging_session"), f.name), []).append(file)
                else:
                    bad.append(f"{f.name} is empty")
            elif hit is None:
                bad.append(f"{f.name} {v!r} is not one of {', '.join(f.levels)}")
            else:
                lay.at[i, f.name] = v = hit[0]
            values.append(v)
            if parsed is not None and v and v != parsed[f.name]:
                origin = "session folder" if f.source == "session" else "file name"
                was = f"from the {origin}: {parsed[f.name]!r}" if parsed[f.name] else f"not readable from the {origin}"
                note = _add_note(note, f"{f.name} set in layout.csv to {v!r} ({was})")
        section, hemi = _cell(row, "section").lower(), _cell(row, "hemisphere").upper()
        if not re.fullmatch(r"n\d+", section):
            bad.append(f"section {section!r} is not a section label like n3" if section else "section is empty")
        if hemi not in ("A", "B"):
            bad.append(f"hemisphere {hemi!r} is not A or B" if hemi else "hemisphere is empty")
        for col, v in (("section", section), ("hemisphere", hemi)):
            if parsed is not None and v and v != parsed[col]:
                was = f"from the file name: {parsed[col]!r}" if parsed[col] else "not readable from the file name"
                note = _add_note(note, f"{col} set in layout.csv to {v!r} ({was})")
        lay.at[i, "section"], lay.at[i, "hemisphere"] = section, hemi
        if bad:
            problems.append(f"{file}: " + "; ".join(bad))
        lay.at[i, "tumour"] = f"{section}{hemi}"
        lay.at[i, "tumour_id"] = tumour_id(values, section, hemi) if not bad and all(values) else ""
        lay.at[i, "note"] = note

    for (session, fname), files in session_missing.items():
        problems.append(f"imaging session {session!r}: {fname} cannot be determined from the folder name "
                        f"({len(files)} tumour{'s' if len(files) != 1 else ''}); fill the {fname} column for its rows "
                        f"in layout.csv (or declare the folder's spelling as an alias in settings.toml)")

    inc = lay[(lay["include"] == "yes") & (lay["tumour_id"] != "")]
    for tid, grp in inc.groupby("tumour_id", sort=True):
        if len(grp) > 1:
            problems.append(f"tumour ID {tid!r} is given by {len(grp)} files: " + "; ".join(grp["file"])
                            + " - exclude all but one (include = no + note) or correct their values")

    if scanned is not None:
        new = sorted(set(scanned.index) - set(lay["file"]))
        if new:
            warnings.append(f"{len(new)} tumour file(s) in the raw folder are not in layout.csv and are not measured: "
                            + "; ".join(new) + " (run `kerncount init` again: your layout.csv is kept, the fresh scan "
                            "goes to layout.new.csv)")
    no_reason = [_cell(r, "file") for _, r in lay.iterrows()
                 if r["include"] == "no" and r["exclusion_reason"] == NO_REASON]
    if no_reason:
        warnings.append("excluded without a reason in note: " + "; ".join(no_reason))
    return Validation(lay, problems, warnings)


def check(layout: pd.DataFrame, raw: Path | None, settings: dict) -> Validation:
    """validate(), raising LayoutBlocked with every problem when there is any."""
    v = validate(layout, raw, settings)
    if v.problems:
        raise LayoutBlocked(v.problems)
    return v
