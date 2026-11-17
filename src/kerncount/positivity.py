"""Which nuclei are positive for a marker (e.g. Ki-67): the cut-off rules (ADR 0001, docs/intensity.md).

Signal of a nucleus = (mean inside the shrunk outline - local background) / flat field at the nucleus, in
"corrected intensity units". The secondary readout uses the 90th percentile instead of the mean (dotted patterns).

Cut-off precedence, decided per n and channel:
    1. negative_control  a negative-control condition (settings [experiment] negative_control) has nuclei in that
                         n: cut-off = control_percentile (default 99th) of their signal
    2. manual            a locked value: [positivity] manual_cutoff in settings.toml, else a lock in cutoffs.json
                         (written by tuning.lock: value, time, the fields it was tuned on); an n-specific lock wins
                         over an all-n lock
    3. lod               limit of detection: signal > k x (background noise SD for that nucleus), noise from the
                         background pixels of the same image (intensity.measure_nuclei: noise_sd). The cut-off
                         differs per nucleus (smaller nuclei average fewer pixels, so need more signal); the median
                         over nuclei is reported as its value.
Without negative-control wells, "positive" means "detectably stained above background", not "above non-specific
staining".
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .fileio import atomic_write, file_lock  # noqa: F401 (re-exported for tuning / gui)

LOCK_FILE = "cutoffs.json"
MIN_CONTROL_NUCLEI = 20
RULE_LABELS = {"negative_control": "negative-control percentile", "manual": "locked manual value",
               "lod": "limit of detection (k x background noise)"}


def _key(s) -> str:
    return re.sub(r"[\s_\-.]+", "", str(s)).lower()


def is_negative_control(treatment: str, settings: dict) -> bool:
    nc = str(settings.get("experiment", {}).get("negative_control", "") or "").strip()
    return bool(nc) and _key(treatment) == _key(nc)


# ---------------------------------------------------------------- locks (cutoffs.json)

def load_locks(project_dir: Path) -> list[dict]:
    p = Path(project_dir) / LOCK_FILE
    if not p.exists():
        return []
    try:
        return list(json.loads(p.read_text(encoding="utf-8")).get("locks", []))
    except (ValueError, OSError):
        return []


def read_lock_file(p: Path) -> dict:
    try:
        return json.loads(Path(p).read_text(encoding="utf-8")) if Path(p).exists() else {}
    except (ValueError, OSError):
        return {}


def lock_fields(lk: dict) -> int:
    """Number of fields a lock was tuned on (stored as a count; older files kept a list of ids)."""
    f = lk.get("fields", 0)
    return len(f) if isinstance(f, list) else int(f or 0)


def save_lock(project_dir: Path, channel: str, cutoff: float, fields: list[str], measure: str = "mean",
              n: str | None = None, extra: dict | None = None) -> dict:
    """Add (replace) a locked cut-off for `channel` (and `n`, None = every n) in cutoffs.json; returns the lock.
    Only the NUMBER of fields is stored (field ids and their files would unblind; tuning.lock keeps the files in
    _work/cutoff_fields.json). Locked read-modify-write with an atomic replace."""
    p = Path(project_dir) / LOCK_FILE
    with file_lock(p):
        data = read_lock_file(p)
        locks = [lk for lk in data.get("locks", [])
                 if not (_key(lk.get("channel")) == _key(channel) and lk.get("n") == n and lk.get("measure", "mean") == measure)]
        lock = dict(channel=channel, n=n, measure=measure, cutoff=float(cutoff),
                    locked_at=time.strftime("%Y-%m-%dT%H:%M:%S"), fields=len(fields), **(extra or {}))
        locks.append(lock)
        data.update(schema="kerncount.cutoffs/1", locks=locks,
                    note="Locked positivity cut-offs (corrected intensity units: mean inside the nucleus minus local "
                         "background, divided by the flat field). Delete an entry to unlock it.")
        atomic_write(p, json.dumps(data, indent=2))
    return lock


def find_lock(locks: list[dict], channel: str, n: str, measure: str = "mean") -> dict | None:
    hits = [lk for lk in locks if _key(lk.get("channel")) == _key(channel) and lk.get("measure", "mean") == measure]
    for lk in hits:
        if lk.get("n") == n:
            return lk
    for lk in hits:
        if lk.get("n") in (None, "", "all"):
            return lk
    return None


def settings_manual(settings: dict, channel: str) -> float | None:
    v = settings.get("positivity", {}).get("manual_cutoff", "")
    if isinstance(v, dict):
        for k, x in v.items():
            if _key(k) == _key(channel) and str(x).strip() != "":
                return float(x)
        return None
    if v in ("", None):
        return None
    return float(v)


# ---------------------------------------------------------------- the decision

def decide(nuc: pd.DataFrame, settings: dict, channel: str, n: str, locks: list[dict] | None = None,
           signal: str = "mean_bgsub_ff", signal90: str = "p90_bgsub_ff",
           lod: str = "thr_lod", lod90: str = "thr90_lod") -> dict:
    """Cut-offs for one n and channel. `nuc`: measurable nuclei (counted, not saturated) with columns treatment,
    the two signals and the two per-nucleus LoD thresholds. Returns
    {"mean": {...}, "p90": {...}} each with rule, value (scalar for reporting), thr (per-nucleus array, aligned
    with nuc), detail."""
    p = settings.get("positivity", {})
    out = {}
    ctrl = nuc["treatment"].map(lambda t: is_negative_control(t, settings)).to_numpy(bool) if len(nuc) else np.zeros(0, bool)
    for measure, col, lcol in (("mean", signal, lod), ("p90", signal90, lod90)):
        x = nuc[col].to_numpy(float)
        if ctrl.sum() >= MIN_CONTROL_NUCLEI:
            pct = float(p.get("control_percentile", 99.0))
            v = float(np.nanpercentile(x[ctrl], pct))
            nc = settings["experiment"]
            out[measure] = dict(rule="negative_control", value=v, thr=np.full(len(x), v),
                                detail=f"{pct:g}th percentile of {int(ctrl.sum())} nuclei in '{nc['negative_control']}' "
                                       f"({nc.get('negative_control_type', '')})")
            continue
        manual = settings_manual(settings, channel) if measure == "mean" else None
        lk = None if manual is not None else find_lock(locks or [], channel, n, measure)
        if manual is not None:
            out[measure] = dict(rule="manual", value=manual, thr=np.full(len(x), manual),
                                detail="[positivity] manual_cutoff in settings.toml")
            continue
        if lk is not None:
            v = float(lk["cutoff"])
            out[measure] = dict(rule="manual", value=v, thr=np.full(len(x), v),
                                detail=f"locked {lk.get('locked_at', '?')} in {LOCK_FILE} "
                                       f"({'n ' + lk['n'] if lk.get('n') else 'all n'}; tuned on "
                                       f"{lock_fields(lk)} fields)")
            continue
        thr = nuc[lcol].to_numpy(float)
        k = float(p.get("lod_k", 3.0))
        out[measure] = dict(rule="lod", value=float(np.nanmedian(thr)) if len(thr) else float("nan"), thr=thr,
                            detail=f"k = {k:g} x background noise SD for each nucleus' area (median shown; "
                                   f"range {np.nanmin(thr):.1f} to {np.nanmax(thr):.1f})" if len(thr) else "no nuclei")
    if out["mean"]["rule"] == "manual" and out["p90"]["rule"] == "lod":
        # a manual value (settings.toml or a tuning lock) is in the units of the MEAN signal and fixes only the
        # mean-based % positive; the p90 readout keeps the automatic rule unless a p90 value is locked separately
        out["p90"]["note"] = P90_NOT_LOCKED
        out["p90"]["detail"] = P90_NOT_LOCKED + "; " + out["p90"]["detail"]
    return out


P90_NOT_LOCKED = ("AUTOMATIC RULE, NOT the manual cut-off: the manual value applies only to the mean signal; lock a "
                  "90th-percentile value in the tuning view (measure p90) to fix this readout too")


def lod_thresholds(nuc: pd.DataFrame, k: float, ff: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Per-nucleus LoD cut-offs (corrected units when the flat field `ff` at the nucleus is given)."""
    f = np.ones(len(nuc)) if ff is None else np.asarray(ff, float)
    thr = k * nuc["noise_sd"].to_numpy(float) / f
    thr90 = (nuc["p90_null"].to_numpy(float) + k * nuc["noise_sd_p90"].to_numpy(float)) / f
    return thr, thr90


SWEEP_FACTORS = np.round(np.geomspace(0.25, 4.0, 41), 4)


def sweep(nuc: pd.DataFrame, thr: np.ndarray, signal: str = "mean_bgsub_ff", factors=SWEEP_FACTORS) -> pd.DataFrame:
    """% positive per treatment (and 'all') when every nucleus' cut-off is multiplied by each factor.
    cutoff = factor x median per-nucleus cut-off (= the cut-off itself for a single-value rule)."""
    x = nuc[signal].to_numpy(float)
    thr = np.asarray(thr, float)
    med = float(np.nanmedian(thr)) if len(thr) else float("nan")
    rows = []
    groups = [("all", np.ones(len(x), bool))] + [(t, (nuc["treatment"] == t).to_numpy()) for t in pd.unique(nuc["treatment"])]
    for f in factors:
        pos = x > f * thr
        for t, m in groups:
            if m.sum():
                rows.append(dict(treatment=t, factor=float(f), cutoff=f * med, frac_positive=float(pos[m].mean()),
                                 nuclei=int(m.sum())))
    return pd.DataFrame(rows)
