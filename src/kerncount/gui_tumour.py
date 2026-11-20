"""Web-interface views of tumour experiments (`[experiment] type = "tumour_volume"`; docs/tumour.md, docs/ui.md).

Every function takes the project folder and returns plain JSON-able data for one page:
- `overview`: experiment home (factors, groups, results, review generation, lock, attention items);
- `layout_info` / `save_layout`: layout.csv with the validation of `tumour.layout.validate` (problems block the run);
- `plan_info` / `save_plan` / `unlock` / `relock`: the analysis plan (factors, readouts, families) and its lock;
  `guard` / `record`: the lock check every settings write of a tumour project goes through;
- `review_info` / `save_corrections`: the **blinded** review sheets and `review/corrections.csv` rounds. These two
  read only the `review/` folder: never `review_key.csv` (code -> tumour ID), `corrections_used.csv`, the layout or
  the results, so the Review and Corrections pages cannot show a group label;
- `results_info` / `stats_info`: plots, tables and the statistics; headline views leave out the sensitivity rows
  (`is_sensitivity`).

Standard library at import time; the tumour modules that need pandas are imported on first use.
"""
from __future__ import annotations

import copy
import re
import time
import tomllib
from pathlib import Path

from . import settings as settings_mod
from .fileio import atomic_write, file_lock
from .tumour import config, lock_log

TYPE = "tumour_volume"
# stats_long families of the sensitivity analyses (tickets 06 and 08): shown on request, never as the headline
SENSITIVITY_SUFFIXES = (" (Welch-type)", " (permutation)", " (without flagged)")
REVIEW = "review"
CORRECTIONS = "review/corrections.csv"
_CODE = re.compile(r"^[A-Za-z0-9_-]{1,40}$")
_MANUAL = re.compile(r"^manual_(\d+)$")
YES = ("yes", "y", "true", "1", "ja", "j")
NO = ("no", "n", "false", "0", "nein")

# per readout, the plots results.analyse / qc / concordance write as plots/<readout>_<kind>.png|pdf
PLOT_KINDS = [
    ("absolute", "Absolute", "Every tumour one dot (hollow = QC-flagged); line and bar = model geometric mean, 95 % CI."),
    ("normalised", "Normalised", "% of the model-estimated reference in the same cell of the other factors; stars = "
                                 "confirmatory family of the focus factor (Holm)."),
    ("forest", "Effects", "Fold-changes with 95 % CIs, per family."),
    ("interaction", "Interaction", "Model geometric means per group, lines joining levels."),
    ("hemispheres_AB", "A vs B", "Hemispheres A and B of the same section side by side (a check, not a test)."),
    ("vs_fiji", "vs Fiji", "This method vs the lab's Fiji legacy value, per tumour (method concordance)."),
    ("threshold_sweep", "Threshold sweep", "Each confirmatory contrast at the threshold × each sweep scale."),
    ("residuals", "Diagnostics", "Residuals by group and QQ plot of the log-scale model."),
]
RESULT_TABLES = ["tumours.csv", "excluded_tumours.csv", "flow.csv", "concordance.csv", "threshold_sweep.csv",
                 "corrections_used.csv", "lock_log.csv"]


def is_tumour(s: dict) -> bool:
    return str((s.get("experiment") or {}).get("type") or "") == TYPE


def is_sensitivity(family: str) -> bool:
    return str(family).endswith(SENSITIVITY_SUFFIXES)


def headline(columns: list[str], rows: list[list[str]]) -> list[list[str]]:
    """stats_long rows without the sensitivity families (Welch-type, permutation, without flagged)."""
    if "family" not in columns:
        return rows
    i = columns.index("family")
    return [r for r in rows if not is_sensitivity(r[i] if i < len(r) else "")]


def _gui():
    from . import gui  # csv helpers (gui imports this module: import on use)

    return gui


def _delim(s: dict) -> str:
    return str((s.get("project") or {}).get("csv_delimiter") or ",")


# ---------------------------------------------------------------- analysis plan and its lock (ADR 0003)

def plan_hash(s: dict) -> str:
    """The lock hash `results.lock` records with the results (raises ValueError on an invalid plan)."""
    from .tumour.results import lock

    return lock(s)[0]


def guard(prj: Path, old: dict, new: dict) -> tuple[str, str] | None:
    """Before a settings write: (old hash, new hash), or PlanLocked when the plan is locked and would change.
    None for projects that are not tumour experiments (before and after)."""
    if not (is_tumour(old) or is_tumour(new)):
        return None
    try:
        old_h = plan_hash(old)
    except (ValueError, KeyError, TypeError):
        old_h = ""
    if is_tumour(new):
        try:
            new_h = plan_hash(new)
        except (KeyError, TypeError) as e:
            raise ValueError(f"the tumour analysis plan in these settings cannot be read: {e}") from None
    else:
        new_h = ""
    lock_log.check_change(prj, old_h, new_h)
    return old_h, new_h


def record(prj: Path, hashes: tuple[str, str] | None) -> None:
    """After a settings write that passed `guard`: log a plan change while results exist."""
    if hashes:
        lock_log.record_change(prj, *hashes)


def write_settings(prj: Path, text: str) -> tuple[str, str] | None:
    """Write settings.toml text through the lock check (tumour projects); returns the (old, new) hashes."""
    st = prj / "settings.toml"
    new = settings_mod._merge(settings_mod.DEFAULTS, tomllib.loads(text))
    with file_lock(st):
        hashes = guard(prj, settings_mod.load(prj), new)
        atomic_write(st, text)
    record(prj, hashes)
    return hashes


def lock_status(prj: Path, s: dict | None = None) -> dict:
    s = settings_mod.load(prj) if s is None else s
    try:
        h, err = plan_hash(s), ""
    except (ValueError, KeyError, TypeError) as e:
        h, err = "", str(e)
    return dict(lock_log.status(prj, h), error=err)


def plan_info(prj: Path) -> dict:
    s = settings_mod.load(prj)
    ts = config.tumour_settings(s)
    facs = [dict(name=str(f.get("name", "")), levels=[str(x) for x in f.get("levels") or []],
                 reference=str(f.get("reference", "")), role=str(f.get("role", "")),
                 source=str(f.get("source", "filename")),
                 aliases={str(k): str(v) for k, v in (f.get("aliases") or {}).items()}) for f in ts["factors"]]
    reads = [dict(name=str(r.get("name") or r.get("column")), column=str(r.get("column", "")),
                  label=str(r.get("label", "")), role=str(r.get("role", "exploratory")),
                  legacy=str(r.get("legacy", ""))) for r in ts["readouts"]]
    fams = [dict(id=str(f.get("id", "")), effect=str(f.get("effect", "")), factors=[str(x) for x in f.get("factors") or []],
                 role=str(f.get("role", "")), label=str(f.get("label", ""))) for f in ts["families"]]
    return dict(factors=facs, readouts=reads, families=fams, normalise_by=str(ts.get("normalise_by") or ""),
                roles=dict(factor=list(config.ROLES), readout=list(config.READOUT_ROLES),
                           family=list(config.FAMILY_ROLES)),
                lock=lock_status(prj, s))


def _levels(v) -> list[str]:
    items = v if isinstance(v, list) else str(v).split(",")
    out = [str(x).strip() for x in items if str(x).strip()]
    if len(set(out)) != len(out):
        raise ValueError(f"levels {out}: a level is listed twice")
    return out


def save_plan(prj: Path, body: dict) -> dict:
    """Apply the plan editor: factor levels / reference / display role, readout roles, family roles. Unchanged
    arrays are not rewritten. The lock check runs before anything is written."""
    from . import settings_form

    st = prj / "settings.toml"
    with file_lock(st):
        text = st.read_text(encoding="utf-8")
        ts = config.tumour_settings(tomllib.loads(text))
        edits = {
            "factors": ({str(f.get("name")): f for f in body.get("factors") or []}, "name"),
            "readouts": ({str(r.get("name")): r for r in body.get("readouts") or []}, "name"),
            "families": ({str(f.get("id")): f for f in body.get("families") or []}, "id"),
        }
        new_text = text
        for key, (by, idk) in edits.items():
            rows = copy.deepcopy(ts[key])
            for r in rows:
                e = by.get(str(r.get(idk) or (r.get("column") if key == "readouts" else "")))
                if e is None:
                    continue
                if key == "factors":
                    if "levels" in e:
                        r["levels"] = _levels(e["levels"])
                    for k in ("reference", "role"):
                        if k in e:
                            r[k] = str(e[k]).strip()
                elif "role" in e:
                    r["role"] = str(e["role"]).strip()
            if rows != ts[key]:
                new_text = settings_form.replace_array(new_text, f"tumour.{key}", rows)
        if new_text == text:
            return dict(ok=True, message="nothing changed", changed=False)
        new = settings_mod._merge(settings_mod.DEFAULTS, tomllib.loads(new_text))
        config.factors(new)  # validated like the run does (reference is a level, roles, aliases)
        config.readouts(new)
        config.families(new)
        hashes = guard(prj, settings_mod.load(prj), new)
        atomic_write(st, new_text)
    record(prj, hashes)
    return dict(ok=True, changed=True, message="analysis plan saved - run the analysis to apply it",
                lock=lock_status(prj))


def unlock(prj: Path, reason: str) -> dict:
    with file_lock(prj / lock_log.LOG):
        lock_log.unlock(prj, reason, lock_status(prj)["settings_lock"])
    return dict(ok=True, message="analysis plan unlocked (logged in lock_log.csv; the next run prints it in the report)",
                lock=lock_status(prj))


def relock(prj: Path) -> dict:
    with file_lock(prj / lock_log.LOG):
        lock_log.relock(prj, lock_status(prj)["settings_lock"])
    return dict(ok=True, message="analysis plan locked again", lock=lock_status(prj))


# ---------------------------------------------------------------- settings: what needs re-running

def rerun_advice(before: dict, after: dict) -> dict:
    """Tumour projects have one step, a run: segmentation, readouts and QC are cached per stack and re-done only when
    their settings changed (the cache key notices). Settings marked reseg make it measure every stack again."""
    from . import settings_form

    changed = []
    for sec in ("project", "experiment", "tumour"):
        b, a = before.get(sec) or {}, after.get(sec) or {}
        for k in sorted(set(b) | set(a)):
            if b.get(k) != a.get(k):
                meta = settings_form.FORM_ITEMS.get(sec, {}).get(k, {})
                changed.append(dict(section=sec, key=k, label=meta.get("label", k), reseg=bool(meta.get("reseg"))))
    if not changed:
        return dict(changed=[], steps=[])
    reseg = [c for c in changed if c["reseg"]]
    why = ", ".join(c["label"] for c in (reseg or changed))
    label = ("Run the analysis: measures every stack again (slow, about 20-60 s per stack)" if reseg else
             "Run the analysis (cached segmentations are reused: quick)")
    return dict(changed=changed, steps=[dict(action="run", options={}, why=why, slow=bool(reseg), label=label)])


# ---------------------------------------------------------------- layout.csv

def _raw(prj: Path, s: dict) -> Path | None:
    from .tumour.measure import raw_dir

    try:
        p = raw_dir(prj, s)
    except SystemExit:
        return None
    return p if p.is_dir() else None


def layout_info(prj: Path) -> dict:
    """layout.csv as {columns, rows} plus the validation (problems block the run, warnings do not), the rebuilt
    tumour ID per row and the number of included tumours per group."""
    from . import discover
    from .tumour import layout as lay_mod

    gui = _gui()
    lp = prj / "layout.csv"
    if not lp.exists():
        return dict(columns=[], rows=[], factors=[], problems=["layout.csv not found: set the project up first"],
                    warnings=[], ids=[], groups=[], unmatched=[])
    s = settings_mod.load(prj)
    cols, rows, _ = gui.read_csv(lp)
    raw = _raw(prj, s)
    facs = config.factors(s)
    v = lay_mod.validate(discover.read_table(lp), raw, s)
    lay = v.layout
    ok = lay[(lay["include"] == "yes") & (lay["tumour_id"] != "")]
    groups = ([dict(group=" / ".join(str(x) for x in (k if isinstance(k, tuple) else (k,))), tumours=int(n))
               for k, n in ok.groupby([f.name for f in facs], sort=True).size().items()] if len(ok) else [])
    problems = list(v.problems)
    if raw is None:
        problems.insert(0, f"the raw folder {s['project']['raw_dir']!r} (settings.toml [project] raw_dir) was not "
                           "found: files cannot be checked")
    return dict(columns=cols, rows=rows, factors=[dict(name=f.name, levels=list(f.levels), source=f.source)
                                                  for f in facs],
                problems=problems, warnings=list(v.warnings), ids=[str(x) for x in lay["tumour_id"]],
                included=int((lay["include"] == "yes").sum()), excluded=int((lay["include"] == "no").sum()),
                groups=groups, unmatched=lay_mod.unmatched_files(raw, discover.read_table(lp)) if raw else [])


def save_layout(prj: Path, rows: list[dict]) -> dict:
    """Write edited cells (factor values, section, hemisphere, include, note) back to layout.csv, row by `file`;
    then validate (the answer lists what still blocks the run)."""
    gui = _gui()
    s = settings_mod.load(prj)
    editable = [f.name for f in config.factors(s)] + ["section", "hemisphere", "include", "note"]
    path = prj / "layout.csv"
    with file_lock(path):
        cols, old, delim = gui.read_csv(path)
        ix = {c: i for i, c in enumerate(cols)}
        if len(rows) != len(old) or "file" not in ix:
            raise ValueError("layout.csv changed on disk - reload the page")
        for r, new in zip(old, rows):
            r += [""] * (len(cols) - len(r))
            if r[ix["file"]] != str(new.get("file", "")):
                raise ValueError("layout.csv changed on disk - reload the page")
            for c in editable:
                if c in ix and c in new:
                    v = " ".join(str(new[c]).split())
                    if c == "include":
                        v = "yes" if v.lower() in YES else "no" if v.lower() in NO else v
                    r[ix[c]] = v
        gui.write_csv(path, cols, old, delim)
    info = layout_info(prj)
    n = len(info["problems"])
    return dict(ok=True, message="layout saved" + (f" - {n} problem(s) still block the run" if n else
                                                   " - run the analysis again to use it"), layout=info)


# ---------------------------------------------------------------- blinded review and corrections

def _manual_cols(cols: list[str]) -> list[str]:
    return sorted((c for c in cols if _MANUAL.match(c)), key=lambda c: int(c.split("_")[1]))


def review_info(prj: Path) -> dict:
    """BLIND: the review sheets and corrections.csv by code only (reads nothing outside review/)."""
    gui = _gui()
    rd = prj / REVIEW
    path = prj / CORRECTIONS
    overviews = sorted((p.name for p in rd.glob("overview_*.jpg")),
                       key=lambda n: int(re.sub(r"\D", "", n) or 0)) if rd.is_dir() else []
    if not path.exists():
        return dict(rows=[], manual=[], overviews=[f"{REVIEW}/{n}" for n in overviews], mtime=None)
    cols, rows, _ = gui.read_csv(path)
    ix = {c: i for i, c in enumerate(cols)}
    manual = _manual_cols(cols)
    get = lambda r, c: r[ix[c]].strip() if c in ix and ix[c] < len(r) else ""  # noqa: E731
    items = []
    for r in rows:
        code = get(r, "code")
        if not _CODE.match(code):
            continue
        img = f"{REVIEW}/{code}.jpg"
        items.append(dict(code=code, objects=get(r, "objects"), qc_flags=get(r, "qc_flags"), note=get(r, "note"),
                          manual={c: get(r, c) for c in manual}, image=img if (prj / img).is_file() else None))
    return dict(rows=items, manual=manual, overviews=[f"{REVIEW}/{n}" for n in overviews],
                mtime=path.stat().st_mtime)


def check_value(v: str, objects: str, where: str) -> str:
    """One correction cell, normalised: '' (no change), an object number >= 1 (re-pick the main mass) or
    'x' / 'x: reason' (exclude). Same rules as `tumour.corrections.read`."""
    from .tumour.corrections import _EXCLUDE, _PICK

    s = " ".join(str(v or "").split())
    if not s:
        return ""
    m = _PICK.match(s)
    if m:
        k = int(m.group(1))
        if k < 1:
            raise ValueError(f"{where}: object numbers start at 1 (1 = the automatic main mass)")
        if objects.strip().isdigit() and k > int(objects):
            raise ValueError(f"{where}: object {k} does not exist (the review image shows objects 1-{objects})")
        return str(k)
    x = _EXCLUDE.match(s)
    if x:
        reason = x.group(1).strip()
        return f"x: {reason}" if reason else "x"
    raise ValueError(f"{where}: '{s}' - type an object number from the review image (re-pick the main mass), "
                     "'x: reason' (exclude), or leave it empty")


def save_corrections(prj: Path, changes: list[dict], mtime: float | None) -> dict:
    """Write review corrections by code. A completely new manual_k (next round) is added if needed, and one empty
    manual column is always kept at the end. Notes are free text."""
    gui = _gui()
    path = prj / CORRECTIONS
    if not path.exists():
        raise ValueError("no review/corrections.csv yet: run the analysis (with review images) first")
    with file_lock(path):
        if mtime is not None and abs(path.stat().st_mtime - float(mtime)) > 1e-3:
            raise ValueError("corrections.csv was changed by another program since this page loaded it (Excel? a "
                             "run?). Reload the page and enter your changes again.")
        cols, rows, delim = gui.read_csv(path)
        if "code" not in cols:
            raise ValueError("review/corrections.csv has no 'code' column")
        for row in rows:
            row += [""] * (len(cols) - len(row))
        ci = cols.index("code")
        by = {r[ci].strip(): r for r in rows}
        for ch in changes:
            col, code = str(ch.get("col", "")), str(ch.get("code", "")).strip()
            r = by.get(code)
            if r is None:
                raise ValueError(f"no review code {code!r} in corrections.csv")
            if col != "note" and not _MANUAL.match(col):
                raise ValueError(f"only manual_N columns and note can be edited, not {col!r}")
            if col not in cols:
                k = int(_MANUAL.match(col).group(1)) if _MANUAL.match(col) else 0  # type: ignore[union-attr]
                have = [int(c.split("_")[1]) for c in _manual_cols(cols)]
                if col == "note" or k != max(have, default=0) + 1:
                    raise ValueError(f"cannot create {col}")
                at = cols.index("note") if "note" in cols else len(cols)
                cols.insert(at, col)
                for row in rows:
                    row.insert(at, "")
            objects = r[cols.index("objects")] if "objects" in cols else ""
            value = str(ch.get("value", ""))
            r[cols.index(col)] = (" ".join(value.split()) if col == "note"
                                  else check_value(value, objects, f"{code} {col}"))
        manual = _manual_cols(cols)
        if manual and any(len(r) > cols.index(manual[-1]) and r[cols.index(manual[-1])].strip() for r in rows):
            nxt = f"manual_{int(manual[-1].split('_')[1]) + 1}"
            at = cols.index(manual[-1]) + 1
            cols.insert(at, nxt)
            for row in rows:
                row.insert(at, "")
        gui.write_csv(path, cols, rows, delim)
    return dict(ok=True, message="corrections saved - run the analysis to apply them", mtime=path.stat().st_mtime)


# ---------------------------------------------------------------- results and statistics

def _exists(prj: Path, rel: str) -> str | None:
    return rel if (prj / rel).is_file() else None


def results_info(prj: Path) -> dict:
    s = settings_mod.load(prj)
    reads = []
    try:
        declared = config.readouts(s)
    except ValueError:
        declared = []
    for ro in sorted(declared, key=lambda r: r.role != "confirmatory"):
        plots = []
        for kind, title, desc in PLOT_KINDS:
            png = _exists(prj, f"plots/{ro.name}_{kind}.png")
            if png:
                plots.append(dict(kind=kind, title=f"{ro.label}: {title}", description=desc, path=png,
                                  pdf=_exists(prj, f"plots/{ro.name}_{kind}.pdf")))
        reads.append(dict(name=ro.name, label=ro.label, role=ro.role, column=ro.column, plots=plots))
    extra = [dict(kind="shape", title="Shape panel (exploratory)", path=p, pdf=_exists(prj, "plots/shape_panel.pdf"),
                  description="Area vs perimeter (circle reference) and the mid-section shape readouts by group.")
             for p in [_exists(prj, "plots/shape_panel.png")] if p]
    extra += [dict(kind="concordance", title=p.stem.replace("_", " "), path=f"plots/{p.name}", pdf=None,
                   description="") for p in sorted((prj / "plots").glob("*.png"))
              if p.stem.startswith(("concordance", "flow"))] if (prj / "plots").is_dir() else []
    prism = sorted(f"prism/{p.name}" for p in (prj / "prism").glob("*.csv")) if (prj / "prism").is_dir() else []
    st = prj / lock_log.STATS
    return dict(readouts=reads, extra=extra, prism=prism,
                tables=[t for t in RESULT_TABLES if (prj / t).is_file()],
                xlsx=_exists(prj, "results.xlsx"), report=_exists(prj, "stats_report.md"),
                built=time.strftime("%Y-%m-%d %H:%M", time.localtime(st.stat().st_mtime)) if st.is_file() else "",
                lock=lock_status(prj, s))


def stats_info(prj: Path, sensitivity: bool = False) -> dict:
    """stats_long.csv: the headline rows (sensitivity families left out unless asked for) and how many were hidden."""
    gui = _gui()
    p = prj / lock_log.STATS
    if not p.is_file():
        return dict(columns=[], rows=[], hidden=0)
    cols, rows, _ = gui.read_csv(p)
    rows = [[("" if v == "nan" else v) for v in r] for r in rows]
    head = headline(cols, rows)
    return dict(columns=cols, rows=rows if sensitivity else head, hidden=0 if sensitivity else len(rows) - len(head))


# ---------------------------------------------------------------- experiment home

def _csv_rows(prj: Path, name: str) -> list[dict]:
    return _gui().read_dicts(prj / name)


def overview(prj: Path) -> dict:
    s = settings_mod.load(prj)
    ts = config.tumour_settings(s)
    lay = layout_info(prj)
    tum = _csv_rows(prj, "tumours.csv")
    excl = _csv_rows(prj, "excluded_tumours.csv")
    used = _csv_rows(prj, "corrections_used.csv")
    lock = lock_status(prj, s)
    yes = lambda v: str(v).strip().lower() in YES  # noqa: E731
    flagged = [r for r in tum if yes(r.get("flagged", ""))]
    satur = [r for r in tum if yes(r.get("saturation_warning", ""))]
    items: list[dict] = []
    ev = lambda label, **k: dict(label=label, **k)  # noqa: E731
    if lay["problems"]:
        items.append(dict(level="error", key="layout", title=f"{len(lay['problems'])} layout problem(s) block the run",
                          text="; ".join(lay["problems"][:6]) + (" …" if len(lay["problems"]) > 6 else ""),
                          evidence=[ev("Layout", type="go", section="layout")]))
    if lay["warnings"]:
        items.append(dict(level="warn", key="layoutwarn", title="Layout notes", text="; ".join(lay["warnings"]),
                          evidence=[ev("Layout", type="go", section="layout")]))
    if lock["error"]:
        items.append(dict(level="error", key="plan", title="The analysis plan in settings.toml cannot be read",
                          text=lock["error"], evidence=[ev("Settings", type="go", section="settings")]))
    if lock["unlocked"]:
        u = lock["unlocked"]
        items.append(dict(level="warn", key="unlocked", title="The analysis plan is unlocked",
                          text=f"Unlocked {u['at']}: {u['reason']}. Changes to factors, readout roles or families are "
                               "possible and logged; the next run locks the plan again and prints the history in "
                               "the report.", evidence=[ev("Analysis plan", type="go", section="settings")]))
    if lock["changed"]:
        items.append(dict(level="warn", key="planchanged", title="The analysis plan changed after the results",
                          text=f"The results were built with plan lock {lock['results_lock']}; settings.toml now has "
                               f"{lock['settings_lock']}. Run the analysis to rebuild the results with the new plan.",
                          evidence=[ev("Run", type="go", section="home", params={"anchor": "run"})]))
    st = prj / lock_log.STATS
    if st.is_file() and (prj / "settings.toml").stat().st_mtime > st.stat().st_mtime + 1 and not lock["changed"]:
        items.append(dict(level="info", key="stale", title="Settings changed after the last run",
                          text="The results may not reflect the current settings: run the analysis again.",
                          evidence=[ev("Settings", type="go", section="settings")]))
    if flagged:
        items.append(dict(level="info", key="flagged", title=f"{len(flagged)} QC-flagged (truncated) tumour(s)",
                          text="Cut by the image edge or with a mid-section outside the imaged layers. They stay in "
                               "the analysis (drawn hollow); the sensitivity analysis without them is in the "
                               "statistics (families '… (without flagged)').",
                          evidence=[ev("tumours.csv", type="table", path="tumours.csv")]))
    if satur:
        items.append(dict(level="info", key="saturation", title=f"{len(satur)} tumour(s) with a saturation warning",
                          text=f"More than {ts.get('qc_saturation_warn')} of the mid-section main mass is saturated "
                               "(8-bit 255). A warning only; see the QC sheet of results.xlsx.",
                          evidence=[ev("tumours.csv", type="table", path="tumours.csv")]))
    if (prj / "layout.new.csv").exists():
        items.append(dict(level="warn", key="layoutnew", title="The raw folder changed since set-up",
                          text="A fresh scan was saved as layout.new.csv. Compare it with layout.csv and use it if "
                               "tumours were added or renamed.", evidence=[ev("Layout", type="go", section="layout")]))
    order = {"error": 0, "warn": 1, "info": 2}
    return dict(
        meta=dict(notes=(s.get("experiment") or {}).get("notes", ""), raw_dir=s["project"]["raw_dir"]),
        factors=lay["factors"], groups=lay["groups"], included=lay.get("included", 0), excluded_layout=lay.get("excluded", 0),
        tumours=len(tum), flagged=len(flagged), saturation=len(satur), excluded=len(excl),
        generation=(used[0].get("generation", "") if used else ""), review=(prj / CORRECTIONS).is_file(),
        threshold_sweep=bool(ts.get("threshold_sweep")), lock=lock,
        results=dict(xlsx=_exists(prj, "results.xlsx"), report=_exists(prj, "stats_report.md"),
                     built=time.strftime("%Y-%m-%d %H:%M", time.localtime(st.stat().st_mtime)) if st.is_file() else ""),
        attention=sorted(items, key=lambda i: order[i["level"]]))
