"""The analysis-plan lock of a tumour project (ADR 0003): once results exist, the plan behind them stays fixed.

The plan = factors (levels, reference), readouts (column, confirmatory / exploratory), families and the focus factor
of normalised values; `results.lock(settings)` hashes it, and every stats_long.csv row carries `lock <hash>`.
Results exist = `stats_long.csv` in the project folder; their hash is read back from it.

While results exist, a change of the plan (a settings.toml whose hash differs from the results') needs an explicit
**unlock with a reason**. Every unlock, re-lock, plan change and the run that rebuilds the results with a different
hash is appended to `lock_log.csv` in the project folder; the run copies the log into the report (`report_lines`), so a
changed lock is visible next to the results. A run ends an unlock: the new results are locked again.

Enforced by the web interface (`gui_tumour`); a settings.toml edited in a text editor is not blocked, but the next run
records "without a logged unlock" in the log and the report.

Standard library only (the web interface imports it without the scientific stack).
"""
from __future__ import annotations

import csv
import io
import re
import time
from pathlib import Path

LOG = "lock_log.csv"
STATS = "stats_long.csv"
COLUMNS = ["at", "action", "results_lock", "settings_lock", "reason"]
_HASH = re.compile(r"\block ([0-9a-f]{12})\b")
_ENDS_UNLOCK = ("lock", "run")


class PlanLocked(ValueError):
    """A change of the locked analysis plan without an unlock."""


def results_lock(project_dir: Path) -> str:
    """The plan hash the current results were built with ("" = no results)."""
    p = Path(project_dir) / STATS
    if not p.is_file():
        return ""
    with p.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            m = _HASH.search(row.get("note") or "")
            if m:
                return m.group(1)
    return ""


def entries(project_dir: Path) -> list[dict]:
    p = Path(project_dir) / LOG
    if not p.is_file():
        return []
    with p.open(encoding="utf-8-sig", newline="") as fh:
        return [{c: (r.get(c) or "") for c in COLUMNS} for r in csv.DictReader(fh)]


def _append(project_dir: Path, action: str, results: str, settings: str, reason: str) -> dict:
    from ..fileio import atomic_write

    row = dict(at=time.strftime("%Y-%m-%d %H:%M:%S"), action=action, results_lock=results, settings_lock=settings,
               reason=reason)
    rows = entries(project_dir) + [row]
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=COLUMNS, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    atomic_write(Path(project_dir) / LOG, buf.getvalue())
    return row


def _active_unlock(project_dir: Path, current: str) -> dict | None:
    """The unlock in force: the latest unlock / re-lock / run entry is an unlock of the current results."""
    for e in reversed(entries(project_dir)):
        if e["action"] == "unlock":
            return e if e["results_lock"] == current else None
        if e["action"] in _ENDS_UNLOCK:
            return None
    return None


def status(project_dir: Path, settings_lock: str) -> dict:
    """{results, results_lock, settings_lock, changed (settings differ from the results' plan), locked, unlocked
    (the unlock entry in force or None), entries}."""
    cur = results_lock(project_dir)
    unl = _active_unlock(project_dir, cur) if cur else None
    return dict(results=bool(cur), results_lock=cur, settings_lock=settings_lock,
                changed=bool(cur) and cur != settings_lock, locked=bool(cur) and unl is None, unlocked=unl,
                entries=entries(project_dir))


def check_change(project_dir: Path, old_lock: str, new_lock: str) -> None:
    """Raise PlanLocked when settings with plan hash `new_lock` may not replace those with `old_lock`."""
    if old_lock == new_lock:
        return
    st = status(project_dir, new_lock)
    if st["locked"]:
        raise PlanLocked(
            f"the analysis plan is locked: results exist (lock {st['results_lock']}), so factors, levels, reference "
            "levels, readout roles and contrast families cannot change. Unlock it first (Settings -> Analysis plan -> "
            "Unlock, with a reason); the unlock and the change are logged in lock_log.csv and printed in the report.")


def unlock(project_dir: Path, reason: str, settings_lock: str) -> dict:
    reason = " ".join(str(reason or "").split())
    if len(reason) < 3:
        raise ValueError("give a reason for unlocking the analysis plan (it is logged and printed in the report)")
    cur = results_lock(project_dir)
    if not cur:
        raise ValueError("nothing to unlock: there are no results yet, the plan can be edited freely")
    if _active_unlock(project_dir, cur):
        raise ValueError("the analysis plan is already unlocked")
    return _append(project_dir, "unlock", cur, settings_lock, reason)


def relock(project_dir: Path, settings_lock: str) -> dict:
    cur = results_lock(project_dir)
    if not cur or not _active_unlock(project_dir, cur):
        raise ValueError("the analysis plan is not unlocked")
    return _append(project_dir, "lock", cur, settings_lock, "locked again without a run")


def record_change(project_dir: Path, old_lock: str, new_lock: str) -> dict | None:
    """Log a saved change of the plan while results exist (call after check_change)."""
    cur = results_lock(project_dir)
    if old_lock == new_lock or not cur:
        return None
    return _append(project_dir, "plan changed", cur, new_lock, f"plan {old_lock} -> {new_lock}")


def record_run(project_dir: Path, new_lock: str) -> dict | None:
    """Called by the run before the results are rewritten: logs a rebuild with a different plan (noting whether it
    was unlocked) and ends a pending unlock."""
    cur = results_lock(project_dir)
    unl = _active_unlock(project_dir, cur) if cur else None
    if cur and cur != new_lock:
        how = (f"after the unlock of {unl['at']} ({unl['reason']})" if unl else
               "without a logged unlock (settings.toml changed outside the web interface?)")
        return _append(project_dir, "run", cur, new_lock, f"results rebuilt with plan {new_lock}, {how}")
    if unl:
        return _append(project_dir, "run", cur, new_lock, "results rebuilt; plan unchanged")
    return None


def report_lines(project_dir: Path) -> list[str]:
    """The lock history for the report / About sheet ([] when nothing was ever unlocked or changed)."""
    es = entries(project_dir)
    if not es:
        return []
    out = [f"Analysis plan lock history ({LOG}, {len(es)} entr{'y' if len(es) == 1 else 'ies'}):"]
    for e in es:
        locks = e["results_lock"] + (f" -> {e['settings_lock']}" if e["settings_lock"] and
                                     e["settings_lock"] != e["results_lock"] else "")
        out.append(f"Lock log {e['at']}: {e['action']} (lock {locks}): {e['reason']}")
    return out
