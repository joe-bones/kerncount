"""Review corrections (ticket 06, ADR 0002): re-pick the main mass or exclude a tumour — nothing else.

The reviewer edits `review/corrections.csv` (blind: rows are review codes, see review.py). Like the cell-count
`counts.csv` (manual_1 -> g2, manual_2 -> g3, ...), each `manual_<k>` column is one round of corrections:
- a whole number = that object (as numbered on the review image) becomes the main mass ("1" = back to automatic);
- `x` (optionally `x: reason`) = exclude the tumour; the reason goes to excluded_tumours.csv and the flow table;
- empty = no change in this round. Later rounds override earlier ones for the same tumour.
The chain stops at the first completely empty round; generation g<k+1> = automatic + k rounds.

Applying corrections never re-segments: the cached labels (`_work/tumour/<slug>.labels.tif`) are relabelled (the
picked object and object 1 swap numbers) and the readouts recomputed from them with ticket 05's readout functions.
Per-tumour threshold changes are not possible (one rule for every tumour, ADR 0002).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .review import CORRECTIONS, read_key, review_dir

USED_FILE = "corrections_used.csv"
_PICK = re.compile(r"^\s*(\d+)(?:\.0+)?\s*$")
_EXCLUDE = re.compile(r"^\s*x\b[\s:;,\-–]*(.*)$", re.IGNORECASE | re.DOTALL)


class CorrectionsError(SystemExit):
    """Unreadable corrections: the run stops and lists every problem."""


@dataclass
class Corrections:
    picks: dict[str, tuple[int, str]] = field(default_factory=dict)  # tumour ID -> (object, round)
    excluded: dict[str, str] = field(default_factory=dict)  # tumour ID -> reason (with its round)
    rounds: list[str] = field(default_factory=list)  # manual_k columns applied
    codes: dict[str, str] = field(default_factory=dict)  # tumour ID -> code

    @property
    def generation(self) -> str:
        return f"g{len(self.rounds) + 1}"

    def describe(self) -> str:
        if not self.rounds:
            return "generation g1 (automatic; no review corrections)"
        n_pick = sum(1 for k, _ in self.picks.values() if k != 1)
        return (f"generation {self.generation} (automatic + review rounds {', '.join(self.rounds)}): "
                f"{n_pick} main mass(es) re-picked, {len(self.excluded)} tumour(s) excluded")


def read(project_dir: Path, settings: dict) -> Corrections:
    """Parse review/corrections.csv against review_key.csv (blocking on anything unreadable)."""
    path = review_dir(project_dir) / CORRECTIONS
    key = read_key(project_dir)
    tid_of = dict(zip(key["code"], key["tumour_id"]))
    out = Corrections(codes=dict(zip(key["tumour_id"], key["code"])))
    if not path.exists():
        return out
    delim = settings.get("project", {}).get("csv_delimiter", ",")
    t = pd.read_csv(path, sep=delim, dtype=str, keep_default_na=False)
    if "code" not in t.columns:
        raise CorrectionsError(f"{path}: no 'code' column")
    rounds = sorted((c for c in t.columns if re.fullmatch(r"manual_\d+", c)), key=lambda c: int(c.split("_")[1]))
    for r in rounds:
        if not t[r].str.strip().any():
            break
        out.rounds.append(r)
    problems = []
    for _, row in t.iterrows():
        code = str(row["code"]).strip()
        for r in out.rounds:
            cell = str(row[r]).strip()
            if not cell:
                continue
            tid = tid_of.get(code)
            if tid is None:
                problems.append(f"{code}: unknown review code (not in {Path(project_dir) / 'review_key.csv'})")
                break
            m, x = _PICK.match(cell), _EXCLUDE.match(cell)
            if m and int(m.group(1)) >= 1:
                out.picks[tid] = (int(m.group(1)), r)
                out.excluded.pop(tid, None)
            elif x:
                reason = x.group(1).strip() or "no reason given"
                out.excluded[tid] = f"excluded at review ({r}): {reason}"
                out.picks.pop(tid, None)
            else:
                problems.append(f"{code} {r} = {cell!r}: write an object number shown on the review image "
                                "(re-pick the main mass) or 'x: reason' (exclude)")
    if problems:
        raise CorrectionsError("review/corrections.csv cannot be applied:\n" + "\n".join(f"  - {p}" for p in problems)
                               + "\nFix the cells and run again.")
    return out


def check_picks(corr: Corrections, n_objects: dict[str, int]) -> None:
    """Every re-picked object must exist in that tumour's labels (n_objects = objects incl. the main mass)."""
    bad = [f"{corr.codes.get(t, t)} {r}: object {k} does not exist (the review image shows objects 1-{n_objects[t]})"
           for t, (k, r) in corr.picks.items() if t in n_objects and k > n_objects[t]]
    if bad:
        raise CorrectionsError("review/corrections.csv cannot be applied:\n" + "\n".join(f"  - {p}" for p in bad))


def relabel(labels: np.ndarray, pick: int) -> np.ndarray:
    """Object `pick` becomes the main mass (label 1); the former main mass takes number `pick`."""
    if pick == 1:
        return labels
    out = labels.copy()
    out[labels == 1] = pick
    out[labels == pick] = 1
    return out


def main_mass_readouts(labels: np.ndarray, stack: np.ndarray, voxel_um: tuple[float, float, float],
                       tumour_settings: dict) -> dict:
    """The readouts that depend on which object is the main mass, recomputed from (relabelled) labels with the same
    functions the segmentation uses (`shape.readouts`, `shape.surface_area`; ticket 05)."""
    from . import shape

    vz, vy, vx = voxel_um
    main = labels == 1
    n = int(np.count_nonzero(main))
    top = np.iinfo(stack.dtype).max if np.issubdtype(stack.dtype, np.integer) else float(stack.max())
    opts = shape.ShapeOptions.from_settings(tumour_settings)
    return {
        "main_mass_voxels": n, "main_mass_volume_um3": round(n * vz * vy * vx, 3),
        "saturated_frac_main_mass": int(np.count_nonzero(stack[main] >= top)) / n if n else float("nan"),
        **shape.label_readouts(labels, voxel_um, opts),
    }


def exclusions(layout: pd.DataFrame, corr: Corrections) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(layout with review-excluded tumours set to include = no, their rows with `exclusion_reason`)."""
    lay = layout.copy()
    hit = lay["tumour_id"].isin(list(corr.excluded)) & (lay["include"].astype(str) == "yes")
    lay.loc[hit, "include"] = "no"
    ex = layout[hit].assign(exclusion_reason=lambda d: d["tumour_id"].map(corr.excluded))
    return lay, ex.drop(columns=[c for c in ("include", "note", "legacy_mask", "legacy_stats") if c in ex.columns])


def write_used(project_dir: Path, corr: Corrections, tumour_ids: list[str]) -> Path:
    """corrections_used.csv (project folder, NOT blind): code, tumour ID, main-mass object, source per tumour."""
    from ..fileio import atomic_write

    rows = []
    for tid in tumour_ids:
        k, r = corr.picks.get(tid, (1, ""))
        if tid in corr.excluded:
            src = corr.excluded[tid]
        elif r and k != 1:
            src = f"{r}: main mass = object {k}"
        elif r:
            src = f"{r}: automatic main mass confirmed"
        else:
            src = "auto"
        rows.append({"code": corr.codes.get(tid, ""), "tumour_id": tid, "main_mass_object": k, "source": src,
                     "generation": corr.generation})
    path = Path(project_dir) / USED_FILE
    atomic_write(path, pd.DataFrame(rows, columns=["code", "tumour_id", "main_mass_object", "source",
                                                   "generation"]).to_csv(index=False))
    return path
