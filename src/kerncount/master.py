"""The master counts.csv (one file for all n) and the generation chain built from it.

counts.csv columns:
    n, treatment, image, auto_count, qc_ok, qc_level, manual_1, manual_2, ...

Generation rules:
    g1        = auto_count for every image (always fully automatic)
    g(k+1)    = g(k), with every non-empty cell of manual_k replacing that image's value
    the chain stops at the first manual_k column that is completely empty.
A manual cell may be a whole number (corrected count) or "x" (exclude this image from then on).
"""
from __future__ import annotations

import hashlib
import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd

from .discover import read_table

KEY = ["n", "treatment", "image"]
FIXED = KEY + ["auto_count", "qc_ok", "qc_level"]
EXCLUDE_TOKENS = {"x", "exclude", "excl", "excluded", "-", "na", "n/a"}
MIN_MANUAL_COLS = 3


def manual_columns(df: pd.DataFrame) -> list[str]:
    cols = [c for c in df.columns if c.lower().startswith("manual_") and c[7:].isdigit()]
    return sorted(cols, key=lambda c: int(c[7:]))


def order_rows(df: pd.DataFrame, treatment_order: list[str]) -> pd.DataFrame:
    present = list(dict.fromkeys(df["treatment"]))
    order = [t for t in treatment_order if t in present] + sorted(t for t in present if t not in treatment_order)
    rank = {t: i for i, t in enumerate(order)}
    out = df.assign(_t=df["treatment"].map(rank), _i=df["image"].astype(int))
    return out.sort_values(["n", "_t", "_i"]).drop(columns=["_t", "_i"]).reset_index(drop=True)


def merge(path: Path, auto: pd.DataFrame, treatment_order: list[str], delimiter: str, backup_dir: Path) -> tuple[pd.DataFrame, list[str]]:
    """Refresh automatic columns in counts.csv, keeping every manual entry."""
    warnings: list[str] = []
    auto = auto[FIXED].copy()
    auto["image"] = auto["image"].astype(int)
    if path.exists():
        old = read_table(path)
        old["image"] = old["image"].astype(int)
        mcols = manual_columns(old)
        keep = old[KEY + mcols]
        if mcols:
            has_manual = keep[mcols].astype(str).apply(lambda s: s.str.strip() != "").any(axis=1)
            lost = keep[has_manual].merge(auto[KEY], on=KEY, how="left", indicator=True)
            for _, r in lost[lost["_merge"] == "left_only"].iterrows():
                warnings.append(f"counts.csv row {r['n']}/{r['treatment']}/{r['image']} has manual values but the "
                                "image is no longer found; row dropped (old file kept in _work/backups)")
        df = auto.merge(keep, on=KEY, how="left")
    else:
        mcols = []
        df = auto
    # pad manual columns: at least 3, and always one empty column at the end
    mcols = manual_columns(df)
    while len(mcols) < MIN_MANUAL_COLS or (mcols and df[mcols[-1]].fillna("").astype(str).str.strip().ne("").any()):
        new = f"manual_{len(mcols) + 1}"
        df[new] = ""
        mcols.append(new)
    df[mcols] = df[mcols].fillna("")
    df = order_rows(df[FIXED + mcols], treatment_order)

    if path.exists():
        text_new = df.to_csv(index=False, sep=delimiter)
        if path.read_text(encoding="utf-8-sig") != text_new:
            backup_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, backup_dir / f"counts_{datetime.now():%Y%m%d-%H%M%S}.csv")
    df.to_csv(path, index=False, sep=delimiter)
    return df, warnings


def _parse(v: str, where: str):
    s = str(v).strip()
    if s == "":
        return None
    if s.lower() in EXCLUDE_TOKENS:
        return "exclude"
    try:
        f = float(s.replace(",", "."))
    except ValueError:
        raise ValueError(f"counts.csv {where}: cannot read manual value '{s}' (use a whole number or x)")
    if f < 0 or f != int(f):
        raise ValueError(f"counts.csv {where}: manual value '{s}' must be a whole number >= 0")
    return int(f)


def chain(df: pd.DataFrame) -> list[dict]:
    """Return generations as dicts: {name, source, table, fingerprint}."""
    base = df[KEY + ["auto_count", "qc_level"]].copy()
    base["image"] = base["image"].astype(int)
    base["count"] = pd.to_numeric(base["auto_count"], errors="coerce")
    base["source"] = "auto"
    base["excluded"] = False
    fp = hashlib.sha1(base[KEY + ["count"]].to_csv(index=False).encode()).hexdigest()[:12]
    gens = [dict(name="g1", source="auto", table=base.copy(), fingerprint=fp, changes=0)]
    cur = base
    for col in manual_columns(df):
        vals = [_parse(v, f"row {i + 2}, {col}") for i, v in enumerate(df[col])]
        if all(v is None for v in vals):
            break
        cur = cur.copy()
        changes = 0
        for i, v in enumerate(vals):
            if v is None:
                continue
            changes += 1
            if v == "exclude":
                cur.at[i, "excluded"] = True
                cur.at[i, "source"] = f"{col}: excluded"
            else:
                cur.at[i, "count"] = v
                cur.at[i, "excluded"] = False
                cur.at[i, "source"] = col
        fp = hashlib.sha1((fp + col + ",".join(map(str, vals))).encode()).hexdigest()[:12]
        gens.append(dict(name=f"g{len(gens) + 1}", source=col, table=cur, fingerprint=fp, changes=changes))
    return gens
