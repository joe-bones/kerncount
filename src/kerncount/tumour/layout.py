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

Ticket 12: any number of session folders (hidden / system folders ignored); session-folder names are read word by
word through levels and aliases (`read_session`; the other words are the `session_label`); `file_date` (optional) and
`acquisition_day` (OME metadata) per tumour; IDs clashing between session folders get the folder appended
(`disambiguate`); `migrate` carries a layout.csv over a renamed / reorganised raw folder (rows matched by file name).

Ticket 13: the structure below a session folder is irrelevant. Factors come from the session-folder name and the
file name only; stacks are found at any depth (intermediate folders are provenance: they appear in `file`). A Zeiss
`.czi` of the same tumour (matched by parsed identity, not by base name) is recorded in `czi`; a tumour present only as
`.czi` is reported (`czi_only`), never silently dropped. `_s<k>` scene suffixes and `1n` section labels are read. The
file-name date is the `culture_batch` (culture preparation, the biological batch); section labels repeat per culture
batch, so IDs clashing between batches get the batch appended. A session folder may be a symbolic link (data on an
external drive); a dangling link (drive not mounted) blocks the run with one clear message per folder.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .. import discover
from .config import Factor, factors as declared_factors

ORIGINAL_SUFFIXES = (".ome.tiff.tif", ".ome.tiff", ".ome.tif")  # longest first; .ome.tiff.tif = re-saved by Fiji
STANDARD_SUFFIXES = (".ome.tiff", ".ome.tif")
CZI_SUFFIX = ".czi"  # Zeiss raw file: provenance only (no CZI reader installed)
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


def _split_czi(name: str) -> tuple[str, str] | None:
    return (name[: -len(CZI_SUFFIX)], name[-len(CZI_SUFFIX):]) if name.lower().endswith(CZI_SUFFIX) else None


def is_czi(name: str) -> bool:
    return _split_czi(name) is not None


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
    scene = re.search(r"[\s_]+[sS]\d+$", text)  # "n1 A_s1": Zeiss scene suffix of an exported tile / position
    if scene:
        notes.append(f"scene suffix '{scene.group(0).strip()}' ignored")
        text = text[: scene.start()]
    if "_" in text or re.search(r"\s{2,}", text) or text != text.strip():
        notes.append("spacing normalised")
    text = re.sub(r"[\s_]+", " ", text).strip()
    joined = re.sub(r"(?<![^\s])([+-])\s+(?=\S)", r"\1", text)  # "- clodro" -> "-clodro"
    if joined != text:
        notes.append("sign joined to its word ('- clodro' -> '-clodro')")
    tokens = []
    for t in _tokens(joined):
        m = re.fullmatch(r"([nN])(\d+)([ABab])", t)  # "n3A" -> section n3 + hemisphere A
        r = re.fullmatch(r"(\d+)([nN])([ABab]?)", t)  # "1n" (ticket 13) -> section n1
        if m:
            tokens += [m.group(1) + m.group(2), m.group(3)]
            notes.append(f"'{t}' read as section {m.group(1).lower()}{m.group(2)} + hemisphere {m.group(3).upper()}")
        elif r:
            tokens += ["n" + r.group(1)] + ([r.group(3)] if r.group(3) else [])
            notes.append(f"section label '{r.group(1)}{r.group(2)}' read as 'n{r.group(1)}'"
                         + (f" + hemisphere {r.group(3).upper()}" if r.group(3) else ""))
        else:
            tokens.append(t)
    return tokens, notes


def _spans(tokens: list[str], longest: int) -> list[str]:
    """Every token and every run of up to `longest` neighbouring tokens joined by a space (multi-word levels and
    aliases such as "TLR7 KO" are matched as a phrase)."""
    out: list[str] = []
    for k in range(longest, 1, -1):  # longest first: "TLR7 KO" is read as one phrase, not as "KO" alone
        out += [" ".join(tokens[i:i + k]) for i in range(len(tokens) - k + 1)]
    return out + list(tokens)


def _read_factor(f: Factor, tokens: list[str]) -> tuple[str, str, str, str]:
    """(level, token, note, problem) for factor f from tokens: a note when read through an alias or in other case,
    a problem when tokens name two different levels. `token` may be a phrase of several tokens ("TLR7 KO")."""
    longest = max([len(x.split()) for x in (*f.levels, *f.aliases)] or [1])
    hits: list[tuple[str, str, bool]] = []
    for t in _spans(tokens, longest):
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


def _used(tokens: list[str], phrase: str) -> set[str]:
    """The tokens a matched phrase consists of (a phrase found through `_spans`)."""
    parts = phrase.split()
    return set(parts) if all(p in tokens for p in parts) else {phrase}


def read_session(name: str, facs: list[Factor]) -> tuple[dict[str, str], str, list[str]]:
    """Session-folder name -> ({factor: level} for factors with source = session ('' = not readable), session label =
    the folder's other words ("5d K07" -> "5d"), the words read as no level)."""
    tokens, _ = _normalise_stem(name)
    values, used = {}, set()
    for f in facs:
        if f.source != "session":
            continue
        level, token, _note, _problem = _read_factor(f, tokens)
        values[f.name] = level
        if level:
            used |= _used(tokens, token)
    rest = [t for t in tokens if t not in used]
    return values, " ".join(rest), rest


def tumour_id(values: list[str], section: str, hemisphere: str) -> str:
    """Unique across the experiment: every factor level + section label + hemisphere (CONTEXT.md)."""
    return " / ".join([*values, f"{section}{hemisphere}"])


SYSTEM_PREFIXES = (".", "_", "@", "$", "#", "~")  # hidden / system / NAS (".snapshots", "@eaDir", "$RECYCLE.BIN")
SYSTEM_FOLDERS = {"system volume information", "lost+found"}


def _is_session_folder(d: Path) -> bool:
    return d.is_dir() and not d.name.startswith(SYSTEM_PREFIXES) and d.name.lower() not in SYSTEM_FOLDERS


def unavailable_sessions(raw: Path) -> dict[str, str]:
    """Session folders that are symbolic links to a place that is not there (external drive not mounted):
    {folder name: link target}."""
    raw = Path(raw)
    if not raw.is_dir():
        return {}
    return {d.name: os.readlink(d) for d in sorted(raw.iterdir())
            if d.is_symlink() and not d.exists() and not d.name.startswith(SYSTEM_PREFIXES)}


def _walk(folder: Path, recursive: bool = True) -> list[tuple[str, str]]:
    """(relative directory with trailing '/' or '', file name) of every visible file below `folder`, at any depth;
    symbolic links are followed (each real folder visited once), hidden / system folders skipped."""
    out: list[tuple[str, str]] = []
    if not recursive:
        return [("", p.name) for p in sorted(folder.iterdir()) if p.is_file() and not p.name.startswith(".")]
    seen: set[str] = set()
    for here, dirs, files in os.walk(folder, followlinks=True):
        real = os.path.realpath(here)
        if real in seen:
            dirs[:] = []
            continue
        seen.add(real)
        dirs[:] = sorted(d for d in dirs if not d.startswith(SYSTEM_PREFIXES) and d.lower() not in SYSTEM_FOLDERS)
        rel = Path(here).relative_to(folder).as_posix()
        rel = "" if rel == "." else rel + "/"
        out += [(rel, f) for f in sorted(files) if not f.startswith(".")]
    return sorted(out)


def _sessions(raw: Path) -> list[tuple[str, Path]]:
    """Every sub-folder of the experiment folder is one imaging session (hidden / system folders ignored)."""
    subs = sorted(d for d in raw.iterdir() if _is_session_folder(d))
    out = [(d.name, d) for d in subs]
    if any(p.is_file() and is_original(p.name) for p in raw.iterdir()):
        out.insert(0, (raw.name, raw))  # flat files directly in the experiment folder: one unnamed session
    return out


def _columns(facs: list[Factor]) -> list[str]:
    return ["tumour_id", "imaging_session", "session_label", "file_date", "culture_batch", "acquisition_day",
            *[f.name for f in facs], "section", "hemisphere", "tumour", "file", "legacy_mask", "legacy_stats", "czi",
            "include", "note"]


# read from the folder / file, filled into older layouts
AUTO = ("session_label", "file_date", "culture_batch", "acquisition_day", "czi")


def _row(name: str, session: str, session_tokens: list[str], names: set[str], rel: str, facs: list[Factor],
         order: str, session_label: str = "", inner: str = "") -> dict:
    """One layout row for the stack `rel + name`; `inner` = the folders between the session folder and the file
    (provenance; only read for a culture-batch date when the file name has none)."""
    split = _split_original(name) or _split_czi(name)
    assert split is not None
    stem, suffix = split
    tokens, notes = _normalise_stem(stem)
    if suffix not in STANDARD_SUFFIXES and suffix != CZI_SUFFIX:
        notes.insert(0, f"file ending '{suffix}' read as an original stack")
    values, missing, used = [], [], set()
    for f in facs:
        level, token, note, problem = _read_factor(f, session_tokens if f.source == "session" else tokens)
        values.append(level)
        if f.source == "filename":
            used |= _used(tokens, token)
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
    batch = d.isoformat() if d else ""
    if not d:  # culture batch: the nearest folder carrying a date (the file name has none)
        for folder in reversed([x for x in inner.split("/") if x]):
            fd = discover.parse_date(folder, order)
            if fd:
                batch = fd.isoformat()
                notes.append(f"no date in the file name: culture batch {batch} read from folder '{folder}'")
                break
    masks, stats = legacy_candidates(stem)
    lower = {n.lower(): n for n in names}
    legacy_mask = next((lower[m.lower()] for m in masks if m.lower() in lower), "")
    legacy_stats = next((lower[s.lower()] for s in stats if s.lower() in lower), "")
    if missing:
        notes.insert(0, "could not read " + ", ".join(missing) + ": fill the empty cells or set include = no")
    return {
        "tumour_id": tumour_id(values, section, hemisphere) if not missing else "",
        "imaging_session": session,
        "session_label": session_label,
        "file_date": d.isoformat() if d else "",  # optional; not the imaging date (acquisition_day is)
        "culture_batch": batch,  # ticket 13: culture preparation (file-name date) = the biological batch
        "acquisition_day": "",
        **{f.name: v for f, v in zip(facs, values)},
        "section": section, "hemisphere": hemisphere, "tumour": f"{section}{hemisphere}" if section and hemisphere
        else "",
        "file": rel + name,
        "legacy_mask": rel + legacy_mask if legacy_mask else "",
        "legacy_stats": rel + legacy_stats if legacy_stats else "",
        "czi": "",
        "include": "yes",  # unreadable rows stay included: the run blocks until they are filled in or excluded
        "note": "; ".join(notes),
    }


def acquisition_day(path: Path) -> str:
    """Day the stack was acquired, from its OME metadata (header only); '' when not recorded or unreadable."""
    from .acquisition import read

    try:
        return read(path)["acquisition_day"]
    except Exception:  # not a TIFF / no OME-XML: the day is simply not known
        return ""


NO_BATCH = "no culture batch"


def disambiguate(lay: pd.DataFrame, rows: pd.Series | None = None) -> None:
    """In place: a tumour ID given by tumours of several imaging sessions (repeat sessions of one genotype) gets the
    session appended ("WT / HERV / -clodro / n1A / 5d WT repeat"); then an ID still given by tumours of several
    culture batches (section labels restart per batch, ticket 13) gets the batch appended
    ("WT / cntrl / -clodro / n1B / 2026-07-10"). IDs that are unique stay as they are, so adding a session or batch
    only changes the IDs of tumours that now clash. Duplicates within one session and batch stay duplicates."""
    for col in ("imaging_session", "culture_batch"):
        if col not in lay.columns:
            continue
        sel = lay["tumour_id"] != ""
        if rows is not None:
            sel &= rows
        for tid, grp in lay[sel].groupby("tumour_id", sort=False):
            if grp[col].nunique() > 1:
                for i in grp.index:
                    lay.at[i, "tumour_id"] = f"{tid} / {lay.at[i, col] or NO_BATCH}"


def base_id(tid: str, n_factors: int) -> str:
    """A tumour ID without the session / culture-batch suffixes `disambiguate` added."""
    return " / ".join(str(tid).split(" / ")[: n_factors + 1])


def _identity(row: dict, facs: list[Factor]) -> tuple[str, ...] | None:
    """What makes a tumour within its session folder: file-name factor levels + section + hemisphere."""
    key = (*[row[f.name] for f in facs if f.source == "filename"], row["section"], row["hemisphere"])
    return key if all(key) else None


def _pair_czi(originals: list[dict], czis: list[dict], facs: list[Factor]) -> None:
    """In place: each .czi row of one session folder is recorded on the original (.ome.tiff) row of the same tumour:
    same identity, compatible culture batch (equal, or one side without a date); among several, the one in the same
    folder. A .czi without exactly one free match stays unpaired (`czi_only` reports it)."""
    for c in czis:
        key = _identity(c, facs)
        if key is None:
            continue
        cands = [o for o in originals if not o["czi"] and _identity(o, facs) == key
                 and (not c["culture_batch"] or not o["culture_batch"] or c["culture_batch"] == o["culture_batch"])]
        if len(cands) > 1:
            here = c["file"].rsplit("/", 1)[0]
            cands = [o for o in cands if o["file"].rsplit("/", 1)[0] == here]
        if len(cands) == 1:
            cands[0]["czi"] = c["file"]


def scan(raw: Path, settings: dict) -> pd.DataFrame:
    raw = Path(raw)
    facs = declared_factors(settings)
    order = settings.get("project", {}).get("date_order", "DMY")
    rows = []
    for session, folder in _sessions(raw):
        files = _walk(folder, recursive=folder != raw)  # flat files in the experiment folder: that folder only
        names: dict[str, set[str]] = {}
        for d, n in files:
            names.setdefault(d, set()).add(n)
        top = "" if folder == raw else f"{session}/"
        session_tokens, _ = _normalise_stem(session)
        _, label, _ = read_session(session, facs)
        originals, czis = [], []
        for d, n in files:
            if is_original(n) or is_czi(n):
                row = _row(n, session, session_tokens, names[d], top + d, facs, order, label, inner=d)
                if is_original(n):
                    row["acquisition_day"] = acquisition_day(folder / d / n)
                    originals.append(row)
                else:
                    czis.append(row)
        _pair_czi(originals, czis, facs)
        rows += originals
    out = pd.DataFrame(rows, columns=_columns(facs))
    disambiguate(out)
    return out


def _all_files(raw: Path) -> list[str]:
    raw = Path(raw)
    out = []
    for session, folder in _sessions(raw):
        top = "" if folder == raw else f"{session}/"
        out += [top + d + n for d, n in _walk(folder, recursive=folder != raw)]
    return out


def unmatched_files(raw: Path, layout: pd.DataFrame) -> list[str]:
    """Files in the raw folder that are neither a tumour nor a legacy output matched to one (shown, never used).
    .czi files are either a tumour's twin (`czi` column) or reported by `czi_only`."""
    known = {str(x) for c in ("file", "legacy_mask", "legacy_stats", "czi") if c in layout.columns
             for x in layout[c]}
    return [f for f in _all_files(raw) if f not in known and not is_czi(f)]


def czi_only(raw: Path, layout: pd.DataFrame) -> list[str]:
    """.czi files that are no tumour's twin: tumours present only as .czi (no CZI reader installed: not measured)."""
    known = set(layout["czi"].astype(str)) if "czi" in layout.columns else set()
    return [f for f in _all_files(raw) if is_czi(f) and f not in known]


# ------------------------------------------------------------------------------------------------ migration


@dataclass
class Migration:
    """An older layout.csv brought in line with the current raw folder (`migrate`)."""
    layout: pd.DataFrame
    moved: list[tuple[str, str]] = field(default_factory=list)  # (old file, new file): matched by file name
    added: list[str] = field(default_factory=list)  # files new in the raw folder (rows added; add_new only)
    gone: list[str] = field(default_factory=list)  # rows whose file is no longer found (kept: the check names them)
    ambiguous: dict[str, list[str]] = field(default_factory=dict)  # old file -> several same-named candidates
    filled: list[str] = field(default_factory=list)  # columns added to an older layout.csv

    @property
    def changed(self) -> bool:
        return bool(self.moved or self.added or self.filled)

    def describe(self) -> list[str]:
        out = []
        if self.moved:
            out.append(f"{len(self.moved)} tumour file(s) moved since layout.csv was written (session folder renamed "
                       "or reorganised); rows matched by file name, your edits (include, note, values you typed) "
                       "carried over:")
            out += [f"  {a} -> {b}" for a, b in self.moved]
        if self.added:
            out.append(f"{len(self.added)} new tumour file(s) added: " + "; ".join(self.added))
        if self.filled:
            out.append("columns added: " + ", ".join(self.filled))
        for f, cands in self.ambiguous.items():
            out.append(f"{f}: not found; several files of that name ({'; '.join(cands)}) - set its row by hand")
        if self.gone:
            out.append(f"{len(self.gone)} row(s) whose file is no longer in the raw folder (kept; exclude them or "
                       "remove the rows): " + "; ".join(self.gone))
        return out


def _parts(note: str) -> list[str]:
    return [p.strip() for p in str(note).split(";") if p.strip()]


def _parse_path(file: str, session: str, facs: list[Factor], order: str) -> dict:
    """What a scan says about `file` (relative path) in imaging session `session` (the file need not exist)."""
    name = file.rsplit("/", 1)[-1]
    tokens, _ = _normalise_stem(session)
    _, label, _ = read_session(session, facs)
    return _row(name, session, tokens, set(), "", facs, order, label)


def migrate(old: pd.DataFrame, fresh: pd.DataFrame, settings: dict, add_new: bool = False) -> Migration:
    """Bring an older layout.csv in line with a fresh scan of the raw folder.

    Rows whose file is still there keep everything; missing columns (session_label, file_date, acquisition_day; the
    old `date` column is renamed to file_date) are filled from the scan. A row whose file is gone is matched to a
    scanned file of the same name (one candidate only): folder-derived cells (paths, session, session-sourced factor
    levels) come from the new place, while values the user had typed (cells differing from what the old path reads
    as), include and the user's part of note are carried over. `add_new` appends rows for files not yet listed."""
    facs = declared_factors(settings)
    order = settings.get("project", {}).get("date_order", "DMY")
    cols = _columns(facs)
    lay = old.copy().fillna("").astype(str).reset_index(drop=True)
    m = Migration(lay)
    if "date" in lay.columns and "file_date" not in lay.columns:
        lay = lay.rename(columns={"date": "file_date"})
        m.filled.append("file_date (was: date)")
    m.filled += [c for c in AUTO if c not in lay.columns]
    for c in cols:
        if c not in lay.columns:
            lay[c] = ""
    fr = fresh.copy().fillna("").astype(str).reset_index(drop=True)
    at = {f: i for i, f in enumerate(fr["file"])}
    exact = set(lay["file"]) & set(at)
    by_name: dict[str, list[str]] = {}
    for f in fr["file"]:
        if f not in exact:
            by_name.setdefault(f.rsplit("/", 1)[-1], []).append(f)
    claimed: set[str] = set()
    values = [f.name for f in facs] + ["section", "hemisphere"]
    rows = []
    for _, r in lay.iterrows():
        row = r.to_dict()
        f = row["file"]
        if f in at:
            new = fr.iloc[at[f]]
            for c in AUTO:
                if not row.get(c, ""):
                    row[c] = new[c]
        else:
            cands = [x for x in by_name.get(f.rsplit("/", 1)[-1], []) if x not in claimed]
            if len(cands) == 1:
                new = fr.iloc[at[cands[0]]]
                claimed.add(cands[0])
                was = _parse_path(f, row.get("imaging_session", ""), facs, order)
                moved = {c: new[c] for c in fr.columns}
                for c in values:  # a typed value (differs from what the old place reads as) is the user's
                    if row.get(c, "") != was.get(c, ""):
                        moved[c] = row.get(c, "")
                moved["include"] = row.get("include", "yes")
                own = [p for p in _parts(row.get("note", "")) if p not in _parts(was["note"])]
                moved["note"] = "; ".join(dict.fromkeys(_parts(new["note"]) + own))
                vals = [moved[x.name] for x in facs]
                moved["tumour"] = moved["section"] + moved["hemisphere"]
                moved["tumour_id"] = (tumour_id(vals, moved["section"], moved["hemisphere"])
                                      if all(vals) and moved["section"] and moved["hemisphere"] else "")
                row.update(moved)
                m.moved.append((f, cands[0]))
            else:
                if len(cands) > 1:
                    m.ambiguous[f] = cands
                m.gone.append(f)
        rows.append(row)
    if add_new:
        for f in fr["file"]:
            if f not in exact and f not in claimed:
                rows.append(fr.iloc[at[f]].to_dict())
                m.added.append(f)
    extra = [c for c in lay.columns if c not in cols]
    out = pd.DataFrame(rows, columns=cols + extra).fillna("")
    if m.moved or m.added:  # moved / added rows may create or resolve clashes between sessions: re-derive suffixes
        out["tumour_id"] = [base_id(t, len(facs)) if t else "" for t in out["tumour_id"]]
        disambiguate(out)
    m.layout = out
    return m


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
    problems: list[str] = []
    warnings: list[str] = []
    if raw is not None and not Path(raw).is_dir():  # ticket 13: data on an external drive that is not mounted
        problems.append(f"the raw folder {str(raw)!r} is not available (external drive not mounted?): connect it and "
                        "run again; existing results stay viewable")
        raw = None
    gone_sessions = unavailable_sessions(raw) if raw is not None else {}
    fresh = scan(raw, settings) if raw is not None else None
    if fresh is not None:  # ticket 12: a renamed / reorganised session folder; older layout.csv columns
        mig = migrate(lay, fresh, settings)
        lay = mig.layout
        if mig.moved:
            warnings.append(f"{len(mig.moved)} tumour file(s) moved since layout.csv was written (session folder "
                            "renamed or reorganised); rows matched by file name, edits carried over: "
                            + "; ".join(f"{a} -> {b}" for a, b in mig.moved)
                            + " - run `kerncount init` to update layout.csv")
        for fn, cands in mig.ambiguous.items():
            warnings.append(f"{fn}: several files of that name in the raw folder ({'; '.join(cands)}): set its row by "
                            "hand")
    for c in _columns(facs):
        if c not in lay.columns:
            lay[c] = ""
    lay["exclusion_reason"] = ""
    scanned = fresh.set_index("file") if fresh is not None else None
    session_missing: dict[tuple[str, str], list[str]] = {}
    unavailable: dict[str, list[str]] = {}

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
            if _cell(row, "imaging_session") in gone_sessions:
                unavailable.setdefault(_cell(row, "imaging_session"), []).append(file)
            else:
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
        day = _cell(row, "acquisition_day")
        if parsed is not None and day and parsed["acquisition_day"] and day != parsed["acquisition_day"]:
            note = _add_note(note, f"acquisition_day set in layout.csv to {day!r} (OME metadata: "
                                   f"{parsed['acquisition_day']!r})")
        if bad:
            problems.append(f"{file}: " + "; ".join(bad))
        lay.at[i, "tumour"] = f"{section}{hemi}"
        lay.at[i, "tumour_id"] = tumour_id(values, section, hemi) if not bad and all(values) else ""
        lay.at[i, "note"] = note

    for (session, fname), files in session_missing.items():
        _, _, unread = read_session(session, facs)
        words = f" (words not read as a level: {', '.join(repr(w) for w in unread)})" if unread else ""
        problems.append(f"imaging session {session!r}: {fname} cannot be determined from the folder name{words} "
                        f"({len(files)} tumour{'s' if len(files) != 1 else ''}); declare the folder's spelling as an "
                        f"alias of a level in settings.toml ([[tumour.factors]] name = \"{fname}\": aliases = "
                        f"{{ \"<word>\" = \"<level>\" }}; a new {fname} = a new level), or fill the {fname} column "
                        f"for its rows in layout.csv")

    for session, files in unavailable.items():
        problems.append(f"session folder {session!r} is not available: it links to {gone_sessions[session]!r}, which "
                        f"is not there (external drive not mounted?); {len(files)} tumour{'s' if len(files) != 1 else ''}"
                        " cannot be measured. Connect the drive and run again, or set include = no for its rows; "
                        "existing results stay viewable")

    disambiguate(lay, lay["include"] == "yes")
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
    if raw is not None:
        only = czi_only(raw, lay)
        if only:
            warnings.append(f"{len(only)} .czi file(s) without an OME-TIFF of the same tumour (tumours present only as "
                            "CZI; no CZI reader is installed, so they are not measured - export them as OME-TIFF next "
                            "to the .czi and run `kerncount init` again): " + "; ".join(only))
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
