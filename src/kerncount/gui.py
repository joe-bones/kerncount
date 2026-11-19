"""Point-and-click front end: a local web app that runs, configures and reviews kerncount projects.

Standard library only, so it starts instantly. Listens on 127.0.0.1 (this computer only). Every API call needs a
per-session token embedded in the page, and requests must be addressed to 127.0.0.1/localhost (blocks other
websites and DNS rebinding). Long commands (init, run, regen, sampling) run as a subprocess
(`python -m kerncount.cli ...`) one at a time, with the log streamed to the page. Small edits (settings,
layout, manual corrections) are written directly, with the same validation the command line uses.
"""
from __future__ import annotations

import csv
import io
import json
import mimetypes
import os
import re
import secrets
import subprocess
import sys
import threading
import time
import tomllib
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from . import settings as settings_mod
from . import settings_form
from .cli import ANALYSIS, DATA, workspace_root
from .fileio import atomic_write, file_lock

CONFIG = Path.home() / ".kerncount_gui.json"
TOKEN = secrets.token_urlsafe(16)
WEB = Path(__file__).with_name("web")
_GEN = re.compile(r"^g\d+$")
_MANUAL = re.compile(r"^manual_(\d+)$", re.I)
MAX_ROWS = 20000
EXCLUDE_TOKENS = {"x", "exclude", "excl", "excluded", "-", "na", "n/a"}  # = master.EXCLUDE_TOKENS
mimetypes.add_type("text/plain", ".toml")
mimetypes.add_type("text/plain", ".csv")
mimetypes.add_type("text/plain", ".md")


def slug(s: str) -> str:  # = pipeline._slug (folder names of review images and generation outputs)
    return re.sub(r"[^\w.\-]+", "_", s).strip("_")


def gen_key(name: str) -> int:
    return int(name[1:]) if _GEN.match(name) else -1


# ---------------------------------------------------------------- CSV helpers (stdlib; Excel-saved files too)

def read_csv(path: Path) -> tuple[list[str], list[list[str]], str]:
    text = path.read_text(encoding="utf-8-sig")
    first = text.splitlines()[0] if text else ""
    try:
        delim = csv.Sniffer().sniff(first, delimiters=",;\t").delimiter
    except csv.Error:
        delim = ","
    rows = list(csv.reader(io.StringIO(text), delimiter=delim))
    return (rows[0] if rows else []), rows[1:], delim


def write_csv(path: Path, columns: list[str], rows: list[list[str]], delim: str):
    buf = io.StringIO()
    w = csv.writer(buf, delimiter=delim, lineterminator=os.linesep)
    w.writerow(columns)
    w.writerows(rows)
    try:
        _atomic(path, buf.getvalue(), newline="")
    except PermissionError:
        raise ValueError(f"cannot write {path.name}: is it open in Excel? Close it and try again.") from None


def _atomic(path: Path, text: str, newline: str | None = None):
    atomic_write(path, text, newline=newline)  # temp file + os.replace


def parse_manual(v: str, where: str):
    """Same rules as master._parse: '' = no change, x = exclude, else a whole number >= 0."""
    s = str(v).strip()
    if s == "":
        return None
    if s.lower() in EXCLUDE_TOKENS:
        return "exclude"
    try:
        f = float(s.replace(",", "."))
    except ValueError:
        raise ValueError(f"{where}: cannot read '{s}' (use a whole number, or x to exclude)") from None
    if f < 0 or f != int(f):
        raise ValueError(f"{where}: '{s}' must be a whole number >= 0")
    return int(f)


# ---------------------------------------------------------------- jobs

class Job:
    def __init__(self, label: str, args: list[str], cwd: Path):
        self.label, self.log, self.returncode = label, [], None
        env = dict(os.environ, KERNCOUNT_WORKSPACE=str(cwd), PYTHONIOENCODING="utf-8", PYTHONUNBUFFERED="1")
        flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        self.proc = subprocess.Popen([sys.executable, "-u", "-m", "kerncount.cli", *args], cwd=cwd, env=env,
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                     encoding="utf-8", errors="replace", creationflags=flags)
        threading.Thread(target=self._pump, daemon=True).start()

    def _pump(self):
        for line in self.proc.stdout:  # universal newlines: progress-bar '\r' becomes a line break
            self.log.append(line.rstrip("\n"))
        self.returncode = self.proc.wait()
        self.log.append("\n--- finished OK ---" if self.returncode == 0 else f"\n--- stopped / failed (code {self.returncode}) ---")

    @property
    def running(self) -> bool:
        return self.returncode is None


def build_args(action: str, exp: str, o: dict) -> list[str]:
    """Command line for a job: options first, then "--" and the experiment, so a name starting with "-" can never
    be read as an option (option values that are free text use the --opt=value form)."""
    return _options(action, o) + ["--", exp]


def _options(action: str, o: dict) -> list[str]:
    gens = [g for g in re.split(r"[\s,;]+", str(o.get("rebuild", ""))) if _GEN.match(g) or g == "all"]
    if action == "init":
        return ["init"] + (["--force"] if o.get("force") else []) + (["--force-layout"] if o.get("force_layout") else [])
    if action == "run":
        args = ["run"]
        if o.get("no_review"):
            args.append("--no-review")
        if o.get("resegment"):
            args.append("--resegment")
        return args + (["--rebuild", *gens] if gens else [])
    if action == "regen":
        return ["regen"] + (["--rebuild", *gens] if gens else [])
    if action == "sampling":
        gen = str(o.get("generation", "")).strip()
        args = ["sampling"] + (["--generation", gen] if _GEN.match(gen) else [])
        it = str(o.get("iterations", "")).strip()
        return args + (["--iterations", it] if it.isdigit() and int(it) > 0 else [])
    if action == "wells-template":
        return ["wells-template"] + (["--force"] if o.get("force") else [])
    if action == "platemap":
        p = platemap_options(o)
        args = ["platemap", "--plates", str(p["plates"]), "--wells", str(p["wells"])]
        if p["plate_format"]:
            args += ["--format", str(p["plate_format"])]
        if p["reference_wells"] is not None:
            args += ["--reference-wells", str(p["reference_wells"])]
        if p["avoid_edge"]:
            args.append("--avoid-edge")
        if p["seed"] is not None:
            args += ["--seed", str(p["seed"])]
        if not p["balance"]:
            args.append("--no-balance")
        if p["negative_control"] is not None:
            args += [f"--negative-control={p['negative_control']}", "--control-wells", str(p["negative_control_wells"])]
        return args
    raise ValueError(f"unknown action {action}")


def _int(v, lo: int, hi: int, what: str) -> int:
    try:
        i = int(str(v).strip())
    except ValueError:
        raise ValueError(f"{what}: give a whole number") from None
    if not lo <= i <= hi:
        raise ValueError(f"{what}: must be between {lo} and {hi}")
    return i


def platemap_options(o: dict) -> dict:
    """Validated plate-map options from the page (same meaning as `kerncount platemap`)."""
    fmt = str(o.get("plate_format", "") or "").strip()
    if fmt and fmt not in ("6", "12", "24", "48", "96"):
        raise ValueError("plate format must be 6, 12, 24, 48 or 96")
    rw = str(o.get("reference_wells", "") or "").strip().lower()
    seed = str(o.get("seed", "") or "").strip()
    nc = o.get("negative_control")
    if nc is not None:
        nc = str(nc).strip()
        if len(nc) > 80 or re.search(r"[\x00-\x1f]", nc) or nc.startswith("-"):
            raise ValueError("negative control: give a short condition name (not starting with -)")
    return dict(plates=_int(o.get("plates", 3), 1, 24, "plates"), plate_format=int(fmt) if fmt else None,
                wells=_int(o.get("wells", 2), 1, 24, "wells per condition"),
                reference_wells=None if rw in ("", "default") else "auto" if rw == "auto" else _int(rw, 1, 48, "reference wells"),
                avoid_edge=bool(o.get("avoid_edge")), seed=_int(seed, 0, 2 ** 31 - 1, "seed") if seed else None,
                balance=o.get("balance", True) is not False, negative_control=nc,
                negative_control_wells=_int(o.get("negative_control_wells", 2), 1, 2, "negative-control wells"))


# ---------------------------------------------------------------- workspace / project state

class State:
    def __init__(self, workspace: Path | None):
        self.lock = threading.Lock()
        self.job: Job | None = None
        self.workspace = workspace
        self.httpd: ThreadingHTTPServer | None = None

    def busy(self) -> bool:
        return bool(self.job and self.job.running)

    def experiments(self) -> list[dict]:
        ws = self.workspace
        if not ws or not (ws / DATA).is_dir():
            return []
        out = []
        for d in sorted(p for p in (ws / DATA).iterdir() if p.is_dir() and not p.name.startswith(".")):
            prj = ws / ANALYSIS / d.name
            out.append(dict(name=d.name, initialised=(prj / "settings.toml").exists(),
                            counted=(prj / "counts.csv").exists(), generations=generations(prj),
                            sampling=gen_dirs(prj / "sampling"),
                            intensity=(prj / "intensity" / "manifest.json").exists()))
        return out

    def project(self, exp: str) -> Path:
        if not exp or exp not in {e["name"] for e in self.experiments()}:
            raise ValueError(f"unknown experiment: {exp}")
        return self.workspace / ANALYSIS / exp


def gen_dirs(root: Path) -> list[str]:
    """Sub-folders named g1, g2, ... (anything else is ignored)."""
    return sorted((g.name for g in root.glob("g*") if g.is_dir() and _GEN.match(g.name)), key=gen_key) if root.is_dir() else []


def generations(prj: Path) -> list[str]:
    return gen_dirs(prj / "generations")


def inside(root: Path, rel: str) -> Path:
    """Resolve a relative path and make sure it stays inside root."""
    p = (root / rel).resolve()
    if p != root.resolve() and root.resolve() not in p.parents:
        raise ValueError("path outside the project")
    return p


def default_workspace(arg: str | None) -> Path | None:
    if arg:
        return Path(arg).expanduser().resolve()
    ws = workspace_root()
    if ws:
        return ws
    for d in Path(__file__).resolve().parents:  # e.g. workspace/kerncount/src/kerncount/gui.py
        if (d / DATA).is_dir():
            return d
    try:
        p = Path(json.loads(CONFIG.read_text())["workspace"])
        return p if (p / DATA).is_dir() else None
    except (OSError, ValueError, KeyError):
        return None


def open_path(p: Path):
    if os.name == "nt":
        os.startfile(p)  # noqa: S606 - opens with the default app (Explorer, Excel, Notepad)
    elif sys.platform == "darwin":
        subprocess.Popen(["open", str(p)])
    else:
        subprocess.Popen(["xdg-open", str(p)])


def pick_folder(start: str) -> str:
    """Native folder picker, in a subprocess (Tk must own the main thread on macOS)."""
    code = ("import sys,tkinter as t,tkinter.filedialog as f;r=t.Tk();r.withdraw();r.attributes('-topmost',True);"
            "print(f.askdirectory(initialdir=sys.argv[1] or None,title='Choose the kerncount workspace folder'))")
    try:
        r = subprocess.run([sys.executable, "-c", code, start], capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return ""
    if r.returncode != 0:
        raise ValueError("no folder picker available here - type or paste the path instead")
    return r.stdout.strip()


# ---------------------------------------------------------------- project data for the page

def project_info(prj: Path) -> dict:
    s = settings_mod.load(prj)
    st = prj / "settings.toml"
    gens = []
    for g in generations(prj):
        meta = prj / "generations" / g / "generation.json"
        try:
            gens.append({**json.loads(meta.read_text(encoding="utf-8")), "name": g})
        except (OSError, ValueError):
            gens.append(dict(name=g))
    return dict(
        path=str(prj), settings=s, settings_text=st.read_text(encoding="utf-8") if st.exists() else "",
        schema=settings_form.schema(), generations=gens,
        sampling=gen_dirs(prj / "sampling"), planner=planner_gens(prj),
        colours=n_colours(prj, gens[-1]["name"] if gens else None),
        intensity=(prj / "intensity" / "manifest.json").exists(),
        files={n: (prj / n).exists() for n in ("layout.csv", "layout.new.csv", "counts.csv", "qc_details.csv",
                                                "fiji_comparison.csv", "review", "_work/cells.csv")})


def review_index(prj: Path) -> list[dict]:
    """n -> treatment -> images, from qc_details.csv (count, QC, flags) joined with the review JPEGs on disk."""
    det = prj / "qc_details.csv"
    if not det.exists():
        return []
    cols, rows, _ = read_csv(det)
    ix = {c: i for i, c in enumerate(cols)}
    get = lambda r, c: r[ix[c]] if c in ix and ix[c] < len(r) else ""  # noqa: E731
    src_cols = [c for c in cols if c == "source_file" or c.startswith("source_file_")]
    markers = [c.removeprefix("readout_frac_positive_") for c in cols
               if c.startswith("readout_frac_positive_") and not c.startswith("readout_frac_positive_p90_")]
    dups = {(r.get("n"), r.get("treatment"), r.get("image")): r.get("shared_with", "") for r in read_dicts(prj / "duplicates.csv")}
    out: dict[str, dict[str, list]] = {}
    for r in rows:
        n, t, img = get(r, "n"), get(r, "treatment"), get(r, "image")
        if not img.isdigit():
            continue
        rel = f"review/{slug(n)}/{slug(t)}/{int(img):02d}.jpg"
        out.setdefault(n, {}).setdefault(t, []).append(dict(
            image=int(img), count=get(r, "auto_count"), low=get(r, "count_low"), high=get(r, "count_high"),
            qc=get(r, "qc_level"), flags=get(r, "qc_flags"), path=rel if (prj / rel).exists() else None,
            sources={c.removeprefix("source_file").lstrip("_") or "all": get(r, c) for c in src_cols if get(r, c)},
            marker={ch: dict(positive=get(r, f"readout_frac_positive_{ch}"), rule=get(r, f"cutoff_rule_{ch}"),
                             measured=get(r, f"n_measured_{ch}")) for ch in markers},
            check=get(r, "channel_check"), duplicate=dups.get((n, t, str(int(img))), "")))
    res = []
    for n, ts in out.items():
        res.append(dict(n=n, treatments=[dict(
            treatment=t, images=imgs,
            overview=(f"review/{slug(n)}/{slug(t)}_overview.jpg" if (prj / "review" / slug(n) / f"{slug(t)}_overview.jpg").exists() else None))
            for t, imgs in ts.items()]))
    return res


def save_layout(prj: Path, rows: list[dict]):
    with file_lock(prj / "layout.csv"):
        _save_layout(prj, rows)


def _save_layout(prj: Path, rows: list[dict]):
    path = prj / "layout.csv"
    cols, old, delim = read_csv(path)
    ix = {c: i for i, c in enumerate(cols)}
    if len(rows) != len(old):
        raise ValueError("layout.csv changed on disk - reload the page")
    for r, new in zip(old, rows):
        if (r[ix["n_folder"]], r[ix["treatment_folder"]]) != (new.get("n_folder"), new.get("treatment_folder")):
            raise ValueError("layout.csv changed on disk - reload the page")
        for c in ("n", "treatment", "include", "date"):
            if c in ix and c in new:
                v = str(new[c]).strip()
                if c in ("n", "treatment") and not v:
                    raise ValueError(f"{r[ix['n_folder']]} / {r[ix['treatment_folder']]}: {c} cannot be empty")
                r[ix[c]] = v
    write_csv(path, cols, old, delim)


def save_counts(prj: Path, changes: list[dict], mtime: float | None) -> dict:
    """Apply in-browser corrections to counts.csv. A completely new manual_k column (next round) is added if needed,
    and one empty manual column is always kept at the end (as `kerncount run` does). Locked read-modify-write with
    an atomic replace (two tabs saving at once cannot interleave)."""
    with file_lock(prj / "counts.csv"):
        return _save_counts(prj, changes, mtime)


def _save_counts(prj: Path, changes: list[dict], mtime: float | None) -> dict:
    path = prj / "counts.csv"
    if mtime is not None and abs(path.stat().st_mtime - mtime) > 1e-3:
        raise ValueError("counts.csv was changed by another program since this page loaded it "
                         "(Excel? a run?). Reload the corrections tab and enter your changes again.")
    cols, rows, delim = read_csv(path)
    n_i, t_i, i_i = cols.index("n"), cols.index("treatment"), cols.index("image")
    key = {(r[n_i], r[t_i], str(int(r[i_i])) if r[i_i].strip().isdigit() else r[i_i]): r for r in rows}
    for ch in changes:
        col = str(ch.get("col", ""))
        if not _MANUAL.match(col):
            raise ValueError(f"only manual_N columns can be edited, not {col}")
        if col not in cols:
            k = int(_MANUAL.match(col).group(1))
            have = [int(_MANUAL.match(c).group(1)) for c in cols if _MANUAL.match(c)]
            if k != max(have, default=0) + 1:
                raise ValueError(f"cannot create {col}")
            cols.append(col)
            for r in rows:
                r.append("")
        r = key.get((str(ch.get("n")), str(ch.get("treatment")), str(ch.get("image"))))
        if r is None:
            raise ValueError(f"no row {ch.get('n')} / {ch.get('treatment')} / image {ch.get('image')} in counts.csv")
        where = f"{ch.get('n')} / {ch.get('treatment')} / image {ch.get('image')}"
        v = parse_manual(ch.get("value", ""), where)
        while len(r) < len(cols):
            r.append("")
        r[cols.index(col)] = "" if v is None else ("x" if v == "exclude" else str(v))
    mcols = sorted((c for c in cols if _MANUAL.match(c)), key=lambda c: int(c[7:]))
    if mcols and any(len(r) > cols.index(mcols[-1]) and r[cols.index(mcols[-1])].strip() for r in rows):
        cols.append(f"manual_{int(mcols[-1][7:]) + 1}")
        for r in rows:
            r.append("")
    write_csv(path, cols, rows, delim)
    return dict(mtime=path.stat().st_mtime)


def counts_data(prj: Path) -> dict:
    path = prj / "counts.csv"
    if not path.exists():
        return dict(columns=[], rows=[], mtime=None)
    cols, rows, _ = read_csv(path)
    return dict(columns=cols, rows=rows, mtime=path.stat().st_mtime, generations=generations(prj))


def table_data(path: Path) -> dict:
    cols, rows, _ = read_csv(path)
    return dict(columns=cols, rows=rows[:MAX_ROWS], total=len(rows), truncated=len(rows) > MAX_ROWS)


def list_files(root: Path, rel: str, recursive: bool) -> list[dict]:
    d = inside(root, rel)
    if not d.is_dir():
        raise ValueError(f"not a folder: {rel}")
    if d != root.resolve() and private(root, d):
        raise ValueError("this folder holds internal files and is not shown")
    it = d.rglob("*") if recursive else d.iterdir()
    out = []
    for p in it:
        if p.name.startswith(".") or p.suffix in (".sig",):
            continue
        if private(root, p) and p.name != "cutoffs.json":  # cutoffs.json is listed; opening it shows the summary
            continue
        out.append(dict(path=p.relative_to(root).as_posix(), name=p.name, dir=p.is_dir(),
                        size=p.stat().st_size if p.is_file() else None))
    return sorted(out, key=lambda f: (not f["dir"], f["path"].lower()))


# ---------------------------------------------------------------- cached readers (results, manifests)

_CACHE: dict = {}
_CACHE_LOCK = threading.Lock()
_PLOT_LOCK = threading.Lock()  # matplotlib (plate maps) is not thread-safe


def _sig(p: Path):
    try:
        st = p.stat()
        return (st.st_mtime_ns, st.st_size)
    except OSError:
        return None


def cached(key, paths: list[Path], build):
    """Value of build(), recomputed when any of `paths` changes (mtime / size)."""
    sig = tuple(_sig(p) for p in paths)
    with _CACHE_LOCK:
        hit = _CACHE.get(key)
        if hit is not None and hit[0] == sig:
            return hit[1]
    val = build()
    with _CACHE_LOCK:
        _CACHE[key] = (sig, val)
    return val


def read_dicts(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    cols, rows, _ = read_csv(path)
    return [dict(zip(cols, r + [""] * (len(cols) - len(r)))) for r in rows]


def read_json(path: Path, default=None):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default


def truthy(v) -> bool:
    return str(v).strip().lower() in ("true", "1", "yes", "y")


def fnum(v):
    try:
        f = float(str(v).strip())
    except ValueError:
        return None
    return None if f != f else f


def gen_dir(prj: Path, gen: str) -> Path:
    if not _GEN.match(str(gen or "")):
        raise ValueError(f"unknown generation {gen!r}")
    d = prj / "generations" / gen
    if not d.is_dir():
        raise ValueError(f"no generation {gen} yet")
    return d


def manifest(prj: Path, gen: str) -> dict:
    p = gen_dir(prj, gen) / "manifest.json"
    return cached(("manifest", str(p)), [p], lambda: read_json(p, {}) or {})


def n_colours(prj: Path, gen: str | None) -> dict:
    """n -> plot colour, from the generation's n_colours.csv files (the plots use the same colours)."""
    if not gen:
        return {}
    out = {}
    for f in sorted((prj / "generations" / gen / "all_n").rglob("n_colours.csv")):
        for r in read_dicts(f):
            if r.get("n") and r.get("colour"):
                out.setdefault(r["n"], r["colour"])
    return out


def stats_long(prj: Path, gen: str) -> dict:
    p = gen_dir(prj, gen) / "stats_long.csv"
    if not p.exists():
        return dict(columns=[], rows=[])

    def build():
        cols, rows, _ = read_csv(p)
        return dict(columns=cols, rows=[[("" if v == "nan" else v) for v in r] for r in rows])
    return cached(("stats", str(p)), [p], build)


def generation_summaries(prj: Path) -> list[dict]:
    out = []
    for g in generations(prj):
        d = prj / "generations" / g
        meta = read_json(d / "generation.json", {}) or {}
        man = (read_json(d / "manifest.json", {}) if (d / "manifest.json").exists() else {}) or {}
        out.append(dict(name=g, source=meta.get("source", ""), changes=meta.get("changes"), built=meta.get("built", ""),
                        rebuilt_because=meta.get("rebuilt_because"), has_manifest=bool(man),
                        n_count=man.get("n_count"), primary=man.get("primary", {}).get("test"),
                        reference=man.get("primary", {}).get("reference"), files=len(man.get("files", [])),
                        has_stats=(d / "stats_long.csv").exists(), has_xlsx=(d / "results.xlsx").exists()))
    return out


# ---------------------------------------------------------------- experiment home: metadata, n's, attention

def attention_items(prj: Path, s: dict, man: dict, gen: str | None, qc: list[dict], design: list[dict]) -> list[dict]:
    """Every issue a biologist should know before reading the results, in plain language, with a link to the
    evidence. level: error (results affected) | warn (check this) | info (good to know)."""
    items = []
    ev = lambda label, **k: dict(label=label, **k)  # noqa: E731
    gpath = f"generations/{gen}" if gen else ""
    # duplicate raw images
    dup = [r for r in read_dicts(prj / "duplicates.csv") if truthy(r.get("cross_condition"))]
    if dup:
        groups: dict = {}
        for r in dup:
            groups.setdefault(r.get("n", ""), set()).add(r.get("treatment", ""))
        where = "; ".join(f"{n}: {', '.join(sorted(t))}" for n, t in groups.items())
        summ = (man.get("duplicates") or {}).get("summary") or []
        items.append(dict(level="error", key="duplicates", title=f"Duplicate images in different conditions ({len(dup)} files)",
                          text=f"Some raw images are byte-for-byte identical although they belong to different conditions "
                               f"({where}). One folder is probably a copy of the other, so the results for these "
                               "conditions in that n cannot be trusted until the original microscope files are checked."
                               + (" Details: " + " ".join(summ[:3]) if summ else ""),
                          evidence=[ev("duplicates.csv", type="table", path="duplicates.csv")]))
    # acceptance failures
    acc = read_dicts(prj / gpath / "acceptance.csv") if gen else []
    fails = [r for r in acc if r.get("result") == "fail" and r.get("criterion") != "overall"]
    if fails:
        enforce = (man.get("acceptance") or {}).get("enforce")
        items.append(dict(level="error" if enforce else "warn", key="acceptance",
                          title=f"{len({r['n'] for r in fails})} n failed the acceptance criteria",
                          text="; ".join(f"{r['n']}: {r['criterion'].replace('_', ' ')} ({r.get('detail', '')})" for r in fails)
                               + (". These n's are left out of the all-n statistics (acceptance enforce = true)." if enforce else
                                  ". They are still included in the statistics (acceptance enforce = false): decide whether to keep them."),
                          evidence=[ev("acceptance.csv", type="table", path=f"{gpath}/acceptance.csv")]))
    # channel alignment / mis-paired files
    al = read_dicts(prj / "intensity" / "alignment.csv")
    mis = [r for r in al if "mis-paired" in r.get("flags", "") or (fnum(r.get("corr")) is not None and fnum(r.get("corr")) < 0.2)]
    shift = [r for r in al if r not in mis and "shifted" in r.get("flags", "")]
    if mis:
        folders = sorted({f"{r['n']} / {r['treatment']}" for r in mis if "folder fails the channel check" in r.get("flags", "")})
        items.append(dict(level="error", key="mispaired", title=f"{len(mis)} field(s) with mis-paired channel files",
                          text=(f"Whole folders fail the channel check ({'; '.join(folders)}): they may hold only one "
                                "channel, or a missing or extra file shifted the pairing; their marker readouts are left "
                                "empty and their counts come from every second file. " if folders else "")
                               + "The marker image does not show the same cells as its DAPI image, so these fields keep their "
                               "count but get no marker readouts: " + "; ".join(
                                   f"{r['n']} / {r['treatment']} / field {r['image']} (files {r.get('files', '')}, r = {r.get('corr', '')})" for r in mis)
                               + ". Check the files: a missing or duplicated file shifts the pairing.",
                          evidence=[ev("alignment table", type="table", path="intensity/alignment.csv"),
                                    ev("Intensity checks", type="go", section="intensity", params={"tab": "checks"})]))
    if shift:
        items.append(dict(level="info", key="shift", title=f"{len(shift)} field(s) with shifted channels",
                          text="The marker and the DAPI image of these fields are offset by more than the warning limit "
                               f"({s.get('channels', {}).get('align_warn_px', 2)} px). Small shifts barely change "
                               "per-nucleus values; large ones do (channels.align = \"correct\" shifts them back).",
                          evidence=[ev("alignment table", type="table", path="intensity/alignment.csv")]))
    # design flags
    # files the folder scan could not use (interleaved: a leftover unpaired file; non-numbered files; empty folders)
    lay = [r for r in read_dicts(prj / "layout.csv")
           if any(k in r.get("note", "") for k in ("leftover", "skipped non-numbered", "no images found"))
           and r.get("include", "yes").strip().lower() in ("yes", "y", "true", "1", "")]
    if lay:
        items.append(dict(level="warn", key="files", title=f"{len(lay)} folder(s) with files that were not used",
                          text="; ".join(f"{r.get('n', '')} / {r.get('treatment', '')} ({r.get('treatment_folder', '')}): "
                                         + " | ".join(x for x in r.get("note", "").split(" | ")
                                                      if any(k in x for k in ("leftover", "skipped", "no images")))
                                         for r in lay[:8])
                               + (f"; and {len(lay) - 8} more" if len(lay) > 8 else "")
                               + ". With channels in separate files, an unpaired file means one channel file is missing "
                                 "or extra: check that the pairs (marker, then DAPI) are still in step.",
                          evidence=[ev("Layout", type="go", section="layout")]))
    gone = [r for r in design if r.get("flags", "").strip() and str(r.get("images", "")).strip() == "0"]
    if gone:
        items.append(dict(level="warn", key="absent", title=f"{len(gone)} condition(s) missing from an n",
                          text="; ".join(f"{r['n']} / {r['treatment']}: {r['flags']}" for r in gone[:8])
                               + (f"; and {len(gone) - 8} more" if len(gone) > 8 else "")
                               + ". The statistics use the n's that have each condition (a comparison needs ≥ 2 n "
                                 "with both groups); check that no folder was left out by mistake.",
                          evidence=[ev("design_check.csv", type="table", path="design_check.csv"),
                                    ev("Layout", type="go", section="layout")]))
    flagged = [r for r in design if r.get("flags", "").strip() and r not in gone]
    if flagged:
        items.append(dict(level="warn", key="design", title=f"{len(flagged)} condition(s) with an unusual design",
                          text="; ".join(f"{r['n']} / {r['treatment']}: {r['flags']}" for r in flagged[:8])
                               + (f"; and {len(flagged) - 8} more" if len(flagged) > 8 else "")
                               + ". Missing images make wells unequal; the statistics handle this, but check that no image "
                                 "was lost by mistake.",
                          evidence=[ev("design_check.csv", type="table", path="design_check.csv")]))
    # QC-low share
    low: dict = {}
    for r in qc:
        d = low.setdefault(r.get("n", ""), [0, 0])
        d[0] += 1
        d[1] += r.get("qc_level") == "low"
    lim = fnum(s.get("acceptance", {}).get("max_low_qc_fraction", 0.5)) or 0.5
    high = {n: v[1] / v[0] for n, v in low.items() if v[0] and v[1] / v[0] > 0.25}
    if high:
        items.append(dict(level="warn" if max(high.values()) > lim else "info", key="qc",
                          title="Many images with low quality (QC low)",
                          text="Share of images rated QC low: " + ", ".join(f"{n} {100 * f:.0f}%" for n, f in high.items())
                               + f" (acceptance limit {100 * lim:.0f}%). Low-quality images have an uncertain count "
                                 "(blur, glare, clumps). Look at them and correct or exclude where needed.",
                          evidence=[ev("Review images (QC low)", type="go", section="review", params={"q": "low"})]))
    # model fallbacks
    fb = [(k, v) for k, v in (man.get("model_status") or {}).items() if v.get("status") == "fallback_rm_anova"]
    if fb:
        items.append(dict(level="warn", key="fallback", title=f"Mixed model failed for {len(fb)} metric(s)",
                          text="For " + ", ".join(k for k, _ in fb) + " the hierarchical model could not be fitted, so "
                               "the stars come from the repeated-measures ANOVA (Prism-like; overconfident when effects "
                               "vary between n's).",
                          evidence=[ev("Results: diagnostics", type="go", section="results", params={"view": "diagnostic"})]))
    # one n
    one = [f for f in (man.get("all_n_folders") or []) if (f.get("n_count") or 0) < 2]
    if man and (man.get("n_count", 0) < 2 or one):
        which = ", ".join(f.get("cell_line") or "all" for f in one) if one else ""
        items.append(dict(level="info", key="n1", title="Only one n" + (f" ({which})" if which else ""),
                          text="With a single independent experiment there is no test across experiments: the all-n "
                               "results show estimates only, and every confidence interval, p-value and star says "
                               "\"needs ≥ 2 n\". Fields and wells are technical replicates; they describe how well this "
                               "experiment was measured, not whether the effect is reproducible.",
                          evidence=[ev("Statistics table", type="go", section="stats")]))
    # positivity
    im = read_json(prj / "intensity" / "manifest.json", {}) or {}
    for c in im.get("cutoffs", []):
        fp = c.get("frac_positive")
        if c.get("rule") == "lod" and c.get("measure") == "mean" and fp is not None and fp > 0.95:
            items.append(dict(level="warn", key=f"cutoff_{c.get('n')}_{c.get('channel')}",
                              title=f"Almost every nucleus is {c.get('channel')}-positive in {c.get('n')} ({100 * fp:.1f}%)",
                              text="The cut-off comes from the automatic detection limit (3 x background noise), which only "
                                   "says the stain is detectable, not that a nucleus is truly positive. Without "
                                   "negative-control wells, set a cut-off by eye in the blinded tuning view and lock it "
                                   "before looking at the statistics.",
                              evidence=[ev("Tune the cut-off (blinded)", type="go", section="intensity", params={"tab": "tune"}),
                                        ev("cut-off table", type="table", path="intensity/cutoffs.csv")]))
    if im and not (s.get("experiment", {}).get("negative_control") or "").strip():
        items.append(dict(level="info", key="negctrl", title="No negative-control wells",
                          text="\"Positive\" therefore means \"detectably stained above background\", not \"above "
                               "non-specific staining\". Reserve 1-2 secondary-only wells per plate next time "
                               "(Plate maps page).",
                          evidence=[ev("Plate maps", type="go", section="plates")]))
    locks = (read_json(prj / "cutoffs.json", {}) or {}).get("locks", [])
    if locks and im and (prj / "cutoffs.json").stat().st_mtime > (prj / "intensity" / "manifest.json").stat().st_mtime:
        items.append(dict(level="warn", key="lockpending", title="A locked cut-off is not applied yet",
                          text="A cut-off was locked (or unlocked) after the last run. Run the analysis to apply it "
                               "(detections are reused, so this is quick); the generations are rebuilt automatically.",
                          evidence=[ev("Run", type="go", section="home", params={"anchor": "run"})]))
    # settings / layout changed since the results were built
    if gen and (prj / "settings.toml").exists():
        g = prj / "generations" / gen / "generation.json"
        if g.exists() and (prj / "settings.toml").stat().st_mtime > g.stat().st_mtime + 1:
            items.append(dict(level="info", key="stale", title=f"Settings changed after {gen} was built",
                              text="The results may not reflect the current settings. Run the analysis again "
                                   "(Settings shows which steps are needed after saving).",
                              evidence=[ev("Settings", type="go", section="settings")]))
    if (prj / "layout.new.csv").exists():
        items.append(dict(level="warn", key="layoutnew", title="The raw folders changed since set-up",
                          text="A fresh scan was saved as layout.new.csv. Compare it with layout.csv and use it if "
                               "folders were added or renamed.",
                          evidence=[ev("Layout", type="go", section="layout")]))
    order = {"error": 0, "warn": 1, "info": 2}
    return sorted(items, key=lambda i: order[i["level"]])


def overview(prj: Path) -> dict:
    s = settings_mod.load(prj)
    gens = generation_summaries(prj)
    gen = gens[-1]["name"] if gens else None
    man = manifest(prj, gen) if gen and gens[-1]["has_manifest"] else {}
    layout = read_dicts(prj / "layout.csv")
    qc = read_dicts(prj / "qc_details.csv")
    design = read_dicts(prj / "design_check.csv")
    acc = read_dicts(prj / "generations" / gen / "acceptance.csv") if gen else []
    colours = n_colours(prj, gen)
    ns: dict = {}
    for r in layout:
        if r.get("include", "yes").strip().lower() not in ("yes", "y", "true", "1", ""):
            continue
        d = ns.setdefault(r.get("n", ""), dict(n=r.get("n", ""), date=r.get("date", ""), timepoints=[], conditions=[],
                                                images=0, cell_lines=[]))
        for k, v in (("timepoints", r.get("timepoint", "")), ("conditions", r.get("treatment", "")),
                     ("cell_lines", r.get("cell_line", "") or s.get("experiment", {}).get("cell_line", ""))):
            if v and v not in d[k]:
                d[k].append(v)
        d["images"] += int(fnum(r.get("images")) or 0)
    for r in qc:
        d = ns.get(r.get("n", ""))
        if d is not None:
            d["qc_images"] = d.get("qc_images", 0) + 1
            d["qc_low"] = d.get("qc_low", 0) + (r.get("qc_level") == "low")
    for r in acc:
        if r.get("criterion") == "overall" and r.get("n") in ns:
            ns[r["n"]]["acceptance"] = r.get("result", "")
            ns[r["n"]]["acceptance_detail"] = r.get("detail", "")
    for n, d in ns.items():
        d["colour"] = colours.get(n)

    def mode(key):
        vals = [r.get(key, "") for r in design if r.get(key, "")]
        return max(set(vals), key=vals.count) if vals else ""
    dsum = dict(conditions=len({r.get("treatment") for r in design}), rows=len(design), wells=mode("wells_expected"),
                fields_per_well=mode("fields_per_well"), fields_source=mode("fields_source"),
                flagged=sum(bool(r.get("flags", "").strip()) for r in design),
                unbalanced=sum(truthy(r.get("unbalanced")) for r in design),
                images=sum(int(fnum(r.get("images")) or 0) for r in design),
                positions=s.get("design", {}).get("field_positions", ""), layout=s.get("design", {}).get("well_layout", ""))
    im = read_json(prj / "intensity" / "manifest.json", {}) or {}
    ch = s.get("channels", {})
    meta = dict(s.get("experiment", {}))
    meta.update(channel_files=ch.get("files", "single"), file_order=ch.get("file_order", []), segment=ch.get("segment", ""),
                channel_names=ch.get("names", []), marker_channels=im.get("channels", []),
                nuclear_channel=im.get("nuclear_channel"),
                reference=s.get("analysis", {}).get("reference"), controls=s.get("analysis", {}).get("normalise_to"),
                cell_lines=man.get("cell_lines", []))
    return dict(meta=meta, ns=list(ns.values()), design=dsum, generations=gens, latest=gen,
                attention=attention_items(prj, s, man, gen, qc, design),
                notes=man.get("notes", []), n_count=man.get("n_count"), primary=man.get("primary", {}),
                intensity=bool(im), cutoffs=im.get("cutoffs", []), sampling=gen_dirs(prj / "sampling"),
                plans=len(list((prj / "plans").glob("platemap_*"))) if (prj / "plans").is_dir() else 0,
                wells_csv=(prj / "wells.csv").exists())


# ---------------------------------------------------------------- intensity page and blinded tuning

LOCK_KEYS = ("channel", "n", "measure", "cutoff", "locked_at", "units")
HISTORY_KEYS = ("action", "channel", "n", "measure", "cutoff", "at", "fields", "previous")


def _nfields(lk: dict) -> int:
    f = lk.get("fields", 0)
    return len(f) if isinstance(f, list) else int(fnum(f) or 0)


def cutoff_locks(prj: Path) -> dict:
    """Sanitized locks and history: only whitelisted keys, `fields` as a count. Field ids or the files a lock was
    tuned on (kept by older versions in cutoffs.json) would unblind, so they never leave the server."""
    data = read_json(prj / "cutoffs.json", {}) or {}
    locks = [{k: lk.get(k) for k in LOCK_KEYS if k in lk} | {"fields": _nfields(lk)} for lk in data.get("locks", [])
             if isinstance(lk, dict)]
    hist = [{k: h.get(k) for k in HISTORY_KEYS if k in h} | ({"fields": _nfields(h)} if "fields" in h else {})
            for h in data.get("history", []) if isinstance(h, dict)]
    return dict(locks=locks, history=hist)


def cutoffs_summary(prj: Path) -> str:
    """What /file and /api/text show for cutoffs.json: the sanitized locks, never the raw file."""
    return json.dumps(dict(schema="kerncount.cutoffs/1 (sanitized view: the files and field ids a lock was tuned "
                                  "on are not shown, to keep tuning blinded)", **cutoff_locks(prj)), indent=1)


def _write_cutoffs(prj: Path, data: dict):
    atomic_write(prj / "cutoffs.json", json.dumps(data, indent=1))


PRIVATE_DIRS = ("_work",)  # detection caches, salts / audit trails, backups: never served to the page


def private(prj: Path, p: Path) -> bool:
    """True for paths the page must not read raw: anything in _work/ (audit trail of the files a cut-off was
    tuned on, caches) and cutoffs.json (served only as a sanitized summary)."""
    rel = _rel(prj, p)
    if rel is None:
        return True
    # case-insensitive: on macOS / Windows "_WORK/x" and "CUTOFFS.JSON" open the same files
    return bool(rel) and (rel.split("/")[0] in PRIVATE_DIRS or rel == "cutoffs.json")


def _rel(prj: Path, p: Path) -> str | None:
    try:
        return p.resolve().relative_to(prj.resolve()).as_posix().casefold()
    except ValueError:
        return None


def is_cutoffs(prj: Path, p: Path) -> bool:
    """The lock file itself (any spelling): served only as the sanitized summary."""
    return _rel(prj, p) == "cutoffs.json"


def intensity_info(prj: Path) -> dict:
    im = read_json(prj / "intensity" / "manifest.json", {}) or {}
    s = settings_mod.load(prj)
    return dict(manifest=im, **cutoff_locks(prj), negative_control=s.get("experiment", {}).get("negative_control", ""),
                manual_cutoff=s.get("positivity", {}).get("manual_cutoff", ""))


def _tuning():
    from . import tuning  # heavy (numpy, pandas, scikit-image): imported on first use only
    return tuning


def tuning_info(prj: Path) -> dict:
    t = _tuning()
    chans = t.channels(prj)
    t.audit(prj)  # moves file names an older version left in cutoffs.json to _work/cutoff_fields.json
    ns = sorted({r.get("n", "") for r in read_dicts(prj / "layout.csv") if r.get("n")})
    return dict(channels=chans, scopes=ns, **cutoff_locks(prj))


def tuning_range(prj: Path, arg) -> dict:
    r = _tuning().suggest_range(prj, arg("channel") or None, arg("measure", "mean"))
    return {k: r[k] for k in ("channel", "measure", "min", "max", "lod_median")}


def tuning_sample(prj: Path, body: dict) -> dict:
    t = _tuning()
    seed = body.get("seed")
    seed = None if seed in (None, "") else _int(seed, 0, 2 ** 31 - 1, "seed")
    r = t.sample(prj, k=_int(body.get("k", 12), 1, 60, "fields"), seed=seed, channel=body.get("channel") or None,
                 n=body.get("n") or None)
    # only the opaque ids: no nucleus count (it tracks the treatment), no n, no condition
    return dict(seed=r["seed"], channel=r["channel"], fields=[dict(field_id=f["field_id"]) for f in r["fields"]])


MIN_SAMPLE_PREVIEW = 3


def tuning_preview(prj: Path, arg) -> dict:
    t = _tuning()
    cut = float(arg("cutoff"))
    ch, measure, n = arg("channel") or None, arg("measure", "mean"), arg("n") or None
    ids = [f for f in arg("fields").split(",") if f]
    keep = ("cutoff", "channel", "measure", "nuclei", "positive", "frac_positive", "lod_frac_positive")
    scope = t.preview(prj, cut, ch, n=n, measure=measure)
    out = {k: scope[k] for k in keep}
    if len(ids) >= MIN_SAMPLE_PREVIEW:  # never for one or two fields: their nucleus count tracks the treatment
        smp = t.preview(prj, cut, ch, measure=measure, field_ids=ids)
        out["sample"] = {k: smp[k] for k in keep}
    return out


def tuning_lock(prj: Path, body: dict) -> dict:
    t = _tuning()
    ids = [str(f) for f in body.get("field_ids", [])]
    if not ids:
        raise ValueError("tune on at least one field before locking")
    cut = float(body.get("cutoff"))
    measure = body.get("measure", "mean")
    n = body.get("n") or None
    with file_lock(prj / "cutoffs.json"):
        before = {(lk["channel"], lk.get("n"), lk.get("measure")): lk["cutoff"] for lk in cutoff_locks(prj)["locks"]}
        lk = t.lock(prj, body.get("channel"), cut, ids, measure=measure, n=n)
        data = read_json(prj / "cutoffs.json", {}) or {}
        data.setdefault("history", []).append(dict(action="lock", channel=lk["channel"], n=n, measure=measure,
                                                   cutoff=cut, at=lk["locked_at"], fields=len(ids),
                                                   previous=before.get((lk["channel"], n, measure))))
        _write_cutoffs(prj, data)
    shown = {k: lk.get(k) for k in LOCK_KEYS if k in lk} | {"fields": len(ids)}
    return dict(ok=True, message=f"cut-off {cut:g} locked", lock=shown)


def tuning_unlock(prj: Path, body: dict) -> dict:
    with file_lock(prj / "cutoffs.json"):
        data = read_json(prj / "cutoffs.json", {}) or {}
        ch, n, measure = str(body.get("channel", "")), body.get("n") or None, body.get("measure", "mean")
        keep, gone = [], []
        for lk in data.get("locks", []):
            hit = str(lk.get("channel")).lower() == ch.lower() and lk.get("n") == n and lk.get("measure", "mean") == measure
            (gone if hit else keep).append(lk)
        if not gone:
            raise ValueError("no such lock")
        data["locks"] = keep
        data.setdefault("history", []).append(dict(action="unlock", channel=ch, n=n, measure=measure,
                                                   cutoff=gone[0].get("cutoff"), at=time.strftime("%Y-%m-%dT%H:%M:%S")))
        _write_cutoffs(prj, data)
    return dict(ok=True, message="lock removed: the next run uses the automatic rule again")


def field_overlay(prj: Path, arg) -> bytes:
    """Marker image of one KNOWN field with positive / negative outlines at the cut-off in force. Not blinded: only
    the review page uses it, and that page already shows the condition."""
    t = _tuning()
    fid = t.field_id(prj, arg("n"), arg("treatment"), int(arg("image")))
    return t.render(prj, fid, arg("channel") or None, cutoff=None, measure=arg("measure", "mean"))


# ---------------------------------------------------------------- planner

def _plan_file(prj: Path, gen: str, line: str | None = None) -> Path:
    from . import planner
    if not _GEN.match(str(gen)):
        raise ValueError(f"no planner results for {gen}: run the sampling analysis first")
    try:
        p = planner.plan_path(prj / "sampling" / str(gen), line or None)
    except (FileNotFoundError, ValueError, KeyError):
        p = None
    if p is None or not p.exists():
        raise ValueError(f"no planner results for {gen}" + (f" / {line}" if line else "") +
                         ": run the sampling analysis first")
    return p


def _plan(prj: Path, gen: str, line: str | None = None):
    from . import planner
    p = _plan_file(prj, gen, line)
    return cached(("plan", str(p)), [p], lambda: planner.load(prj, gen, line or None))


def plan_lines(prj: Path, gen: str) -> list[dict]:
    """Cell lines planned separately in sampling/<gen>/ (cell_lines.json); [] = one plan for the generation."""
    data = read_json(prj / "sampling" / gen / "cell_lines.json", {}) or {}
    return [dict(cell_line=x.get("cell_line"), path=x.get("path"), n_count=x.get("n_count"))
            for x in data.get("cell_lines", []) if isinstance(x, dict)]


def planner_gens(prj: Path) -> list[str]:
    return [g for g in gen_dirs(prj / "sampling") if (prj / "sampling" / g / "planner.json").exists()
            or (prj / "sampling" / g / "cell_lines.json").exists()]


def recommendations_text(prj: Path, gen: str, sub: str = "") -> str:
    p = prj / "sampling" / gen / sub / "REPORT.txt"
    if not p.exists():
        return ""
    t = p.read_text(encoding="utf-8", errors="replace")
    if t.startswith("Recommendations for the next experiment"):
        cut = t.find("\n\n\n")
        return t[:cut] if cut > 0 else t
    return ""


def planner_info(prj: Path, gen: str | None, line: str | None = None) -> dict:
    gens = planner_gens(prj)
    if not gens:
        return dict(generations=[], gen=None)
    gen = gen if gen in gens else gens[-1]
    lines = plan_lines(prj, gen)
    cur = next((x for x in lines if line in (x["cell_line"], x["path"])), lines[0] if lines else None)
    sub = cur["path"] if cur else ""
    plan = _plan(prj, gen, cur["cell_line"] if cur else None)
    man = read_json(prj / "sampling" / gen / sub / "manifest.json", {}) or {}
    out = dict(generations=gens, gen=gen, cell_lines=lines, cell_line=cur["cell_line"] if cur else plan.get("cell_line"),
               base=f"sampling/{gen}/" + (f"{sub}/" if sub else ""), excluded_ns=plan.get("excluded_ns", []),
               design=plan["design"], options=plan["options"], costs=plan.get("costs"),
               components=plan["components"], classes=plan.get("classes"), default_class=plan.get("default_class"),
               notes=plan.get("notes", []), effects=[dict(reference=e["reference"], treatment=e["treatment"],
                                                          ratio=e.get("ratio"), n_pairs=e.get("n_pairs"), cls=e.get("cls"))
                                                     for e in plan.get("effects", [])],
               bio=sorted(plan.get("bio", {})), metric=plan.get("metric"),
               recommendations=recommendations_text(prj, gen, sub), files=man.get("files", []),
               family_size=max(int(plan["design"]["conditions"]) - 1, 1))
    try:  # can it plan at all? (needs >= 2 n's with the reference)
        from . import planner
        planner.evaluate(plan, max(plan["design"]["n_count"], 2), plan["design"]["wells"], plan["design"]["fields_per_well"], 0.7)
        out["ok"] = True
    except (ValueError, KeyError) as e:
        out["ok"], out["error"] = False, str(e)
    return out


def _flag(v, default: bool) -> bool:
    if v in (None, ""):
        return default
    if isinstance(v, str):
        return v.strip().lower() not in ("0", "false", "no", "off")
    return bool(v)


def planner_evaluate(prj: Path, body: dict) -> dict:
    from . import planner
    plan = _plan(prj, str(body.get("gen", "")), body.get("cell_line") or None)
    d = plan["design"]
    n = _int(body.get("n", max(d["n_count"], 2)), 2, 200, "n")
    wells = _int(body.get("wells", d["wells"]), 1, 48, "wells")
    fields = _int(body.get("fields", d["fields_per_well"]), 1, 200, "fields")
    eff = body.get("effect", 0.7)
    if isinstance(eff, str) and fnum(eff) is not None:
        eff = fnum(eff)
    if isinstance(eff, (int, float)) and not 0.01 <= eff <= 100:
        raise ValueError("effect: a ratio such as 0.7 (30% fewer)")
    rw = body.get("reference_wells")
    rw = None if rw in (None, "", "same") else "auto" if rw == "auto" else _int(rw, 1, 48, "reference wells")
    kw = dict(alpha=float(body.get("alpha", plan["options"]["alpha"])),
              target_power=float(body.get("target_power", plan["options"]["target_power"])),
              reference=body.get("reference") or None, positions=body.get("positions") or None,
              adjust_for_family=_flag(body.get("adjust_for_family"), plan["options"].get("adjust_for_family", True)))
    if not 0 < kw["alpha"] < 0.5 or not 0.05 < kw["target_power"] < 1:
        raise ValueError("alpha must be between 0 and 0.5, target power between 0.05 and 1")
    res = planner.evaluate(plan, n, wells, fields, eff, reference_wells=rw, **kw)
    vary = body.get("vary", "n")
    if vary not in ("n", "wells", "fields"):
        raise ValueError("vary must be n, wells or fields")
    values = dict(n=range(2, max(plan["options"]["max_n"], n) + 1), wells=range(1, max(6, wells) + 1),
                  fields=range(1, max(3 * d["fields_per_well"], 12, fields) + 1))[vary]
    fixed = dict(n=n, wells=wells, fields=fields, reference_wells=rw)
    fixed.pop(vary)
    pts = planner.curve(plan, vary=vary, fixed=fixed, effect=eff, values=values, **kw)
    keep = ("x", "power", "ci_halfwidth_pct", "mde_ratio", "cost")
    return dict(result=res, curve=[{k: p[k] for k in keep} for p in pts], vary=vary)


# ---------------------------------------------------------------- plate maps

def plans_info(prj: Path) -> dict:
    root = prj / "plans"
    plans = []
    for d in (sorted(root.glob("platemap_*"), reverse=True) if root.is_dir() else []):
        if not d.is_dir():
            continue
        files = sorted(p.name for p in d.iterdir() if p.is_file())
        rel = f"plans/{d.name}"
        eu = d / "edge_usage.csv"
        plans.append(dict(name=d.name, path=rel,
                          plates=[f"{rel}/{f}" for f in files if f.startswith("plate_") and f.endswith(".png")],
                          files=[f"{rel}/{f}" for f in files],
                          readme=(d / "README.txt").read_text(encoding="utf-8", errors="replace") if (d / "README.txt").exists() else "",
                          edge=table_data(eu) if eu.exists() else None))
    s = settings_mod.load(prj)
    try:
        from . import platemap
        conds = platemap.project_conditions(prj, s)
    except Exception:  # noqa: BLE001 - only used to pre-fill the form
        conds = []
    a = s.get("analysis", {})
    refs = [r for r in dict.fromkeys([a.get("reference")] + list(a.get("normalise_to") or [])) if r and r in conds]
    pl = s.get("planner", {})
    return dict(plans=plans, conditions=conds, references=refs,
                defaults=dict(plate_format=pl.get("plate_format", 24), avoid_edge=pl.get("avoid_edge", True),
                              reference_wells=pl.get("reference_wells", "auto"),
                              negative_control=s.get("experiment", {}).get("negative_control", "")))


def make_platemap(prj: Path, body: dict) -> dict:
    from . import platemap
    o = platemap_options(body)
    with _PLOT_LOCK:
        out = platemap.plan_project(prj, plates=o["plates"], plate_format=o["plate_format"], wells=o["wells"],
                                    reference_wells=o["reference_wells"], avoid_edge=o["avoid_edge"] or None,
                                    seed=o["seed"], balance=o["balance"], negative_control=o["negative_control"],
                                    negative_control_wells=o["negative_control_wells"])
    return dict(ok=True, message=f"plate maps written to plans/{out.name}", name=out.name)


# ---------------------------------------------------------------- settings: what needs re-running

RUN_SECTIONS = {"project", "experiment", "channels", "segmentation", "counting", "qc", "intensity", "positivity",
                "readouts", "design", "review"}
REGEN_SECTIONS = {"analysis", "acceptance"}


def rerun_advice(before: dict, after: dict) -> dict:
    """Which steps a settings change needs: re-detection (slow), a run, a rebuild of the generations, or the
    sampling analysis (planner)."""
    changed = []
    for sec in sorted(set(before) | set(after)):
        b, a = before.get(sec, {}), after.get(sec, {})
        if not isinstance(b, dict) or not isinstance(a, dict):
            continue
        for k in sorted(set(b) | set(a)):
            if b.get(k) != a.get(k):
                meta = settings_form.FORM_ITEMS.get(sec, {}).get(k, {})
                changed.append(dict(section=sec, key=k, label=meta.get("label", k), reseg=bool(meta.get("reseg"))))
    lab = lambda cs: ", ".join(c["label"] for c in cs)  # noqa: E731
    reseg = [c for c in changed if c["reseg"]]
    run = [c for c in changed if c["section"] in RUN_SECTIONS]
    regen = [c for c in changed if c["section"] in REGEN_SECTIONS]
    plan = [c for c in changed if c["section"] == "planner"]
    steps = []
    if reseg:
        steps.append(dict(action="run", options=dict(rebuild="all"), why=lab(reseg), slow=True,
                          label="Run the analysis: re-detects every nucleus (slow, about 2-3 s per image), then rebuilds the generations"))
    elif run:
        steps.append(dict(action="run", options=dict(rebuild="all"), why=lab(run),
                          label="Run the analysis and rebuild the generations (detections are reused: quick)"))
    elif regen:
        steps.append(dict(action="regen", options=dict(rebuild="all"), why=lab(regen),
                          label="Rebuild the generations (statistics only, no re-detection)"))
    if plan:
        steps.append(dict(action="sampling", options={}, why=lab(plan),
                          label="Re-run the sampling analysis (planner); plate maps use the new values directly"))
    return dict(changed=changed, steps=steps)


# ---------------------------------------------------------------- HTTP

# The page: scripts only from /static (no inline script, no inline event handlers: app.js dispatches data-on*
# attributes), styles may be inline, images also from data:/blob:. Never framed (clickjacking).
CSP_APP = ("default-src 'self'; img-src 'self' data: blob:; style-src 'self' 'unsafe-inline'; script-src 'self'; "
           "object-src 'none'; base-uri 'none'; form-action 'self'; frame-ancestors 'none'")
# Project files (/file): whatever a file contains, it can never run script in this origin.
CSP_FILE = "default-src 'none'; img-src 'self' data:; style-src 'unsafe-inline'; sandbox; frame-ancestors 'none'"
# project files that browsers would treat as active content: always downloaded, never rendered
ACTIVE_SUFFIXES = {".html", ".htm", ".xhtml", ".svg", ".svgz", ".xml", ".xsl", ".xslt", ".js", ".mjs", ".pdf"}
SAFE_TYPES = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg", ".gif": "image/gif",
              ".webp": "image/webp", ".tif": "image/tiff", ".tiff": "image/tiff", ".csv": "text/plain",
              ".txt": "text/plain", ".toml": "text/plain", ".md": "text/plain", ".log": "text/plain",
              ".json": "application/json",
              ".xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"}


def file_type(name: str) -> tuple[str, bool]:
    """(Content-Type, force download) for a project file. Only a whitelist of passive types is shown inline;
    HTML, SVG, XML, scripts (and anything unknown) are sent as a download (application/octet-stream)."""
    suf = Path(name).suffix.lower()
    if suf in ACTIVE_SUFFIXES:
        return ("application/pdf" if suf == ".pdf" else "application/octet-stream"), True
    t = SAFE_TYPES.get(suf)
    return (t, False) if t else ("application/octet-stream", True)


def content_disposition(name: str) -> str:
    """attachment with an ASCII fallback name (quotes, backslashes, control and non-ASCII characters replaced) and
    the exact name as RFC 5987 filename*=UTF-8''..., so any file name gives a valid latin-1 header."""
    from urllib.parse import quote

    fallback = "".join(c if 32 <= ord(c) < 127 and c not in '"\\;' else "_" for c in name) or "download"
    return f"attachment; filename=\"{fallback}\"; filename*=UTF-8''{quote(name, safe='')}"


def response_headers(ctype: str, length: int, extra: dict | None) -> dict:
    h = {"Content-Type": ctype + ("; charset=utf-8" if ctype.startswith(("text/", "application/json")) else ""),
         "Content-Length": str(length), "Cache-Control": "no-store", "X-Content-Type-Options": "nosniff",
         "X-Frame-Options": "DENY", "Referrer-Policy": "no-referrer", "Content-Security-Policy": CSP_APP}
    h.update(extra or {})
    return h


def make_handler(state: State):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *a):  # keep the console quiet
            pass

        def _send(self, code: int, body, ctype="application/json", extra: dict | None = None):
            data = body if isinstance(body, bytes) else body.encode() if isinstance(body, str) else json.dumps(body, default=str).encode()
            # every header is built (and checked to be encodable) BEFORE the status line goes out, so a bad header
            # becomes a clean error response instead of a half-written one
            headers = response_headers(ctype, len(data), extra)
            try:
                for v in headers.values():
                    v.encode("latin-1")
            except UnicodeEncodeError:
                data = json.dumps({"error": "cannot send this file name"}).encode()
                code, headers = 500, response_headers("application/json", len(data), None)
            self.send_response(code)
            for k, v in headers.items():
                self.send_header(k, v)
            self.end_headers()
            self.wfile.write(data)

        def _local(self) -> bool:
            """Only answer requests addressed to this computer (blocks DNS-rebinding attacks)."""
            port = self.server.server_address[1]
            ok = {f"127.0.0.1:{port}", f"localhost:{port}"}
            origin = self.headers.get("Origin")
            if self.headers.get("Host") not in ok or (origin and origin.removeprefix("http://") not in ok):
                self._send(403, {"error": "forbidden"})
                return False
            return True

        def _token_ok(self, q: dict) -> bool:
            return secrets.compare_digest(self.headers.get("X-Token") or q.get("t", [""])[0], TOKEN)

        def do_GET(self):
            if not self._local():
                return
            u = urlparse(self.path)
            q = parse_qs(u.query)
            if u.path == "/":
                # the token goes into a <meta> tag (read by app.js): no inline script, and a cross-origin page can
                # neither read this HTML nor include it as a script
                html = (WEB / "index.html").read_text(encoding="utf-8").replace("__TOKEN__", TOKEN)
                return self._send(200, html, "text/html")
            if u.path.startswith("/static/"):
                f = WEB / u.path.removeprefix("/static/")
                if f.parent == WEB and f.is_file():
                    return self._send(200, f.read_bytes(), mimetypes.guess_type(f.name)[0] or "application/octet-stream")
                return self._send(404, {"error": "not found"})
            if not self._token_ok(q):
                return self._send(403, {"error": "bad token - reload the page"})
            arg = lambda k, d="": q.get(k, [d])[0]  # noqa: E731
            try:
                if u.path == "/file":
                    prj = state.project(arg("exp"))
                    p = inside(prj, arg("path"))
                    if is_cutoffs(prj, p):
                        return self._send(200, cutoffs_summary(prj), "application/json")
                    if private(prj, p):
                        return self._send(403, {"error": "internal file: not shown"})
                    if not p.is_file():
                        return self._send(404, {"error": "not found"})
                    ctype, attach = file_type(p.name)
                    extra = {"Content-Security-Policy": CSP_FILE}
                    if arg("dl") or attach:
                        extra["Content-Disposition"] = content_disposition(p.name)
                    return self._send(200, p.read_bytes(), ctype, extra)
                if u.path == "/api/tuning/render":  # blinded: the field is known only by its opaque id
                    cut = arg("cutoff")
                    png = _tuning().render(state.project(arg("exp")), arg("field"), arg("channel") or None,
                                           cutoff=float(cut) if cut != "" else None, measure=arg("measure", "mean"))
                    return self._send(200, png, "image/png")
                if u.path == "/api/field_overlay":
                    return self._send(200, field_overlay(state.project(arg("exp")), arg), "image/png")
                return self._send(200, self._get(u.path, arg))
            except Exception as e:
                return self._send(400, {"error": str(e)})

        def _get(self, path, arg):
            if path == "/api/state":
                j = state.job
                return dict(workspace=str(state.workspace or ""), experiments=state.experiments(),
                            job=j and dict(label=j.label, running=j.running, returncode=j.returncode))
            if path == "/api/log":
                since = int(arg("since", "0") or 0)
                j = state.job
                lines = j.log[since:] if j else []
                return dict(lines=lines, next=since + len(lines))
            prj = state.project(arg("exp"))
            if path == "/api/project":
                return project_info(prj)
            if path == "/api/review":
                return review_index(prj)
            if path == "/api/counts":
                return counts_data(prj)
            if path == "/api/table":
                p = inside(prj, arg("path"))
                if private(prj, p):
                    raise ValueError("internal file: not shown")
                return table_data(p) if p.exists() else dict(columns=[], rows=[], total=0)
            if path == "/api/text":
                p = inside(prj, arg("path"))
                if is_cutoffs(prj, p):
                    return dict(text=cutoffs_summary(prj))
                if private(prj, p):
                    raise ValueError("internal file: not shown")
                return dict(text=p.read_text(encoding="utf-8-sig", errors="replace")[:500_000] if p.is_file() else "")
            if path == "/api/files":
                return list_files(prj, arg("dir"), arg("recursive") == "1")
            if path == "/api/overview":
                return overview(prj)
            if path == "/api/manifest":
                gen = arg("gen")
                return dict(manifest=manifest(prj, gen), colours=n_colours(prj, gen))
            if path == "/api/stats":
                return stats_long(prj, arg("gen"))
            if path == "/api/intensity":
                return intensity_info(prj)
            if path == "/api/tuning/info":
                return tuning_info(prj)
            if path == "/api/tuning/range":
                return tuning_range(prj, arg)
            if path == "/api/tuning/preview":
                return tuning_preview(prj, arg)
            if path == "/api/planner":
                return planner_info(prj, arg("gen") or None, arg("cl") or None)
            if path == "/api/plans":
                return plans_info(prj)
            raise ValueError(f"unknown request {path}")

        def do_POST(self):
            if not self._local():
                return
            if not self._token_ok({}):
                return self._send(403, {"error": "bad token - reload the page"})
            try:
                body = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
                res = self._post(urlparse(self.path).path, body)
                self._send(200, res if isinstance(res, dict) else {"ok": True, "message": res})
            except Exception as e:  # report every problem on the page
                self._send(400, {"error": str(e)})

        def _idle(self):
            if state.busy():
                raise ValueError("a job is running - wait for it to finish (or press Stop) before changing files")

        def _post(self, path, body):
            if path == "/api/workspace":
                p = Path(str(body.get("path", "")).strip().strip('"')).expanduser()
                if body.get("create"):
                    if not p.is_absolute():
                        raise ValueError("give a full folder path")
                    (p / DATA).mkdir(parents=True, exist_ok=True)
                if not (p / DATA).is_dir():
                    raise ValueError(f"{p} has no 'data' folder inside it (tick 'create' to make one)")
                state.workspace = p.resolve()
                CONFIG.write_text(json.dumps({"workspace": str(state.workspace)}))
                return "workspace set"
            if path == "/api/browse":
                return dict(path=pick_folder(str(state.workspace or "")))
            if path == "/api/new_experiment":
                name = str(body.get("name", "")).strip()
                if (not state.workspace or not name or name in (".", "..") or re.search(r'[\\/:*?"<>|]', name)
                        or name.startswith("-")):
                    raise ValueError("give a simple folder name (not starting with -; no / \\ : * ? \" < > |)")
                d = state.workspace / DATA / name
                d.mkdir(parents=True, exist_ok=True)
                open_path(d)
                return f"created data/{name} - copy one folder per n into it (opened in your file browser), then press Set up"
            if path == "/api/stop":
                if state.busy():
                    state.job.proc.terminate()
                return "stopping"
            if path == "/api/quit":
                if state.busy():
                    state.job.proc.terminate()
                threading.Thread(target=state.httpd.shutdown, daemon=True).start()
                return "kerncount has stopped - you can close this tab"
            if path == "/api/start":
                exp = body.get("experiment", "")
                state.project(exp)
                args = build_args(body.get("action", ""), exp, body.get("options", {}))
                with state.lock:
                    if state.busy():
                        raise ValueError("another job is still running - wait or press Stop")
                    state.job = Job(f"{args[0]} {exp}", args, state.workspace)
                return "started"
            prj = state.project(body.get("experiment", ""))
            if path == "/api/open":
                what = body.get("what")
                gens = generations(prj)
                target = {"folder": prj, "data": state.workspace / DATA / body["experiment"]}.get(what)
                if target is None and what == "results" and gens:
                    target = prj / "generations" / gens[-1] / "results.xlsx"
                if target is None and body.get("path"):
                    target = inside(prj, body["path"])
                    if private(prj, target) and not is_cutoffs(prj, target):
                        raise ValueError("internal file: not shown")
                if target is None or not target.exists():
                    raise ValueError("not there yet (run the earlier steps first)")
                open_path(target)
                return f"opened {target.name}"
            if path == "/api/settings":
                self._idle()
                st = prj / "settings.toml"
                before = settings_mod.load(prj)
                if "text" in body:
                    text = str(body["text"])
                    try:
                        tomllib.loads(text)
                    except tomllib.TOMLDecodeError as e:
                        raise ValueError(f"settings.toml is not valid: {e}") from None
                else:
                    text = settings_form.update_text(st.read_text(encoding="utf-8"), settings_form.coerce(body.get("values", {})))
                with file_lock(st):
                    atomic_write(st, text)
                return dict(ok=True, message="settings saved", **rerun_advice(before, settings_mod.load(prj)))
            if path == "/api/settings_reset":
                self._idle()
                raw = settings_mod.load(prj)["project"]["raw_dir"]
                raw_abs = Path(raw) if Path(raw).is_absolute() else (prj / raw)
                with file_lock(prj / "settings.toml"):
                    bk = settings_mod.write_template(prj, raw_abs)  # backs the old file up first
                return ("settings reset to defaults" +
                        (f"; the previous settings.toml is kept in {bk.relative_to(prj).as_posix()}" if bk else ""))
            if path == "/api/layout":
                self._idle()
                save_layout(prj, body.get("rows", []))
                return "layout saved"
            if path == "/api/layout_use_new":
                self._idle()
                new = prj / "layout.new.csv"
                if not new.exists():
                    raise ValueError("there is no fresh scan (layout.new.csv)")
                with file_lock(prj / "layout.csv"):
                    bk = settings_mod.backup(prj, "layout.csv")
                    new.replace(prj / "layout.csv")
                return ("layout.csv replaced by the fresh scan" +
                        (f"; the previous layout.csv is kept in {bk.relative_to(prj).as_posix()}" if bk else ""))
            if path == "/api/counts":
                self._idle()
                return dict(ok=True, message="corrections saved", **save_counts(prj, body.get("changes", []), body.get("mtime")))
            if path == "/api/tuning/sample":
                return tuning_sample(prj, body)
            if path == "/api/tuning/lock":
                self._idle()
                return tuning_lock(prj, body)
            if path == "/api/tuning/unlock":
                self._idle()
                return tuning_unlock(prj, body)
            if path == "/api/planner/evaluate":
                return planner_evaluate(prj, body)
            if path == "/api/platemap":
                return make_platemap(prj, body)
            raise ValueError(f"unknown request {path}")

    return Handler


def _warm_up():
    """Import the scientific modules in the background, so the first planner / tuning request is fast."""
    try:
        from . import planner, tuning  # noqa: F401
    except Exception:  # noqa: BLE001 - only a speed-up
        pass


def serve(workspace: str | None = None, port: int = 8765, open_browser: bool = True):
    state = State(default_workspace(workspace))
    for p in (port, 0):  # fall back to any free port
        try:
            httpd = ThreadingHTTPServer(("127.0.0.1", p), make_handler(state))
            break
        except OSError:
            continue
    state.httpd = httpd
    threading.Thread(target=_warm_up, daemon=True).start()
    url =f"http://127.0.0.1:{httpd.server_address[1]}/"
    print(f"kerncount is open in your browser at {url}\n"
          "Keep this window open while you work; close it (or press Quit on the page) to stop.", flush=True)
    if open_browser:
        threading.Timer(0.5, webbrowser.open, (url,)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        if state.busy():
            state.job.proc.terminate()
