"""Duplicate raw images: byte-identical files used for different images (e.g. a copied treatment folder).

`check(project_dir, table)` hashes every raw file listed in the `source_file` column of qc_details.csv (SHA-256 of
the file contents). Hashes are cached in `_work/file_hashes.json` by path, size and modification time, so a re-run
reads only new or changed files. Writes `duplicates.csv` in the project folder (one row per image that shares its
content with another image; empty table with the header when there are none) and returns that table.

Columns of duplicates.csv:
    group           1, 2, ... one number per set of identical files
    n, treatment, image, source_file
    hash            SHA-256 of the file contents
    copies          how many images share this content
    cross_condition True if the copies belong to more than one condition (n / treatment): a real data problem
    shared_with     the other images of the group as "n / treatment / image"
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

COLUMNS = ["group", "n", "treatment", "image", "source_file", "hash", "copies", "cross_condition", "shared_with"]
CHUNK = 1 << 20


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(CHUNK)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def file_hashes(paths, cache_path: Path | None) -> dict[str, str]:
    """{path: sha256} for the existing files among `paths`, using and refreshing the cache."""
    cache = {}
    if cache_path is not None and cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text())
        except (OSError, ValueError):
            cache = {}
    out, changed = {}, False
    for p in dict.fromkeys(str(x) for x in paths if isinstance(x, str) and x.strip()):
        fp = Path(p)
        try:
            st = fp.stat()
        except OSError:
            continue
        c = cache.get(p)
        if c and c.get("size") == st.st_size and c.get("mtime_ns") == st.st_mtime_ns:
            out[p] = c["sha256"]
            continue
        try:
            h = _sha256(fp)
        except OSError:
            continue
        cache[p] = dict(size=st.st_size, mtime_ns=st.st_mtime_ns, sha256=h)
        out[p] = h
        changed = True
    if cache_path is not None and changed:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache, indent=0))
    return out


def find(table: pd.DataFrame, hashes: dict[str, str]) -> pd.DataFrame:
    """Images (rows of `table` with n, treatment, image, source_file) whose file content occurs more than once."""
    if "source_file" not in table or table.empty:
        return pd.DataFrame(columns=COLUMNS)
    t = table[["n", "treatment", "image", "source_file"]].copy()
    t["hash"] = t["source_file"].map(lambda s: hashes.get(str(s)))
    t = t[t["hash"].notna()]
    sizes = t.groupby("hash")["source_file"].transform("size")
    dup = t[sizes > 1].copy()
    if dup.empty:
        return pd.DataFrame(columns=COLUMNS)
    rows = []
    for gi, (h, g) in enumerate(sorted(dup.groupby("hash"), key=lambda x: (str(x[1]["n"].iloc[0]),
                                                                          str(x[1]["treatment"].iloc[0]),
                                                                          int(x[1]["image"].iloc[0]))), 1):
        conds = set(zip(g["n"].astype(str), g["treatment"].astype(str)))
        for r in g.itertuples():
            others = [f"{o.n} / {o.treatment} / {int(o.image)}" for o in g.itertuples() if o.Index != r.Index]
            rows.append(dict(group=gi, n=r.n, treatment=r.treatment, image=int(r.image), source_file=r.source_file,
                             hash=h, copies=len(g), cross_condition=len(conds) > 1, shared_with="; ".join(others)))
    return pd.DataFrame(rows, columns=COLUMNS)


def check(project_dir: Path, table: pd.DataFrame, log=print) -> pd.DataFrame:
    """Hash the raw files of qc_details (`table`), write duplicates.csv in the project folder, log a warning."""
    project_dir = Path(project_dir)
    if "source_file" not in table:
        dup = pd.DataFrame(columns=COLUMNS)
    else:
        files = table["source_file"].dropna().astype(str).tolist()
        dup = find(table, file_hashes(files, project_dir / "_work" / "file_hashes.json"))
    dup.to_csv(project_dir / "duplicates.csv", index=False)
    if len(dup):
        cross = dup[dup["cross_condition"]]
        for line in summary_lines(dup):
            log(f"  WARNING: duplicate raw images: {line}")
        if len(cross):
            log("  WARNING: identical image files appear under different conditions. One of the copies is not the "
                "image it claims to be; check against the original microscope files (details: duplicates.csv).")
    return dup


def summary_lines(dup: pd.DataFrame) -> list[str]:
    """One line per pair of conditions sharing images, e.g. 'n3 / Control 1-12 = n3 / LyoVec 1-12 (12 images)'."""
    if dup is None or len(dup) == 0:
        return []
    pairs: dict[tuple, list] = {}
    for _, g in dup.groupby("group"):
        conds = sorted(set(zip(g["n"].astype(str), g["treatment"].astype(str))))
        imgs = {c: sorted(g.loc[(g["n"].astype(str) == c[0]) & (g["treatment"].astype(str) == c[1]), "image"].astype(int))
                for c in conds}
        key = tuple(conds)
        pairs.setdefault(key, []).append(imgs)
    lines = []
    for conds, groups in pairs.items():
        parts = []
        for c in conds:
            ims = sorted({i for g in groups for i in g.get(c, [])})
            parts.append(f"{c[0]} / {c[1]} images {_ranges(ims)}")
        kind = "same condition" if len(conds) == 1 else "different conditions"
        lines.append(" = ".join(parts) + f" ({len(groups)} identical file(s); {kind})")
    return lines


def _ranges(xs: list[int]) -> str:
    out, start, prev = [], None, None
    for x in xs:
        if start is None:
            start = prev = x
        elif x == prev + 1:
            prev = x
        else:
            out.append(f"{start}-{prev}" if prev != start else f"{start}")
            start = prev = x
    if start is not None:
        out.append(f"{start}-{prev}" if prev != start else f"{start}")
    return ", ".join(out)


def affected(dup: pd.DataFrame) -> dict[str, list[str]]:
    """{n: [treatments]} of images shared across conditions."""
    out: dict[str, list[str]] = {}
    if dup is None or len(dup) == 0:
        return out
    for r in dup[dup["cross_condition"]].itertuples():
        out.setdefault(str(r.n), [])
        if r.treatment not in out[str(r.n)]:
            out[str(r.n)].append(r.treatment)
    return out
