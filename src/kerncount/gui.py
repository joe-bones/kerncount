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
import tomllib
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from . import settings as settings_mod
from . import settings_form
from .cli import ANALYSIS, DATA, workspace_root

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
        path.write_text(buf.getvalue(), encoding="utf-8", newline="")
    except PermissionError:
        raise ValueError(f"cannot write {path.name}: is it open in Excel? Close it and try again.") from None


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
    gens = [g for g in re.split(r"[\s,;]+", str(o.get("rebuild", ""))) if _GEN.match(g) or g == "all"]
    if action == "init":
        return ["init", exp] + (["--force"] if o.get("force") else [])
    if action == "run":
        args = ["run", exp]
        if o.get("no_review"):
            args.append("--no-review")
        if o.get("resegment"):
            args.append("--resegment")
        return args + (["--rebuild", *gens] if gens else [])
    if action == "regen":
        return ["regen", exp] + (["--rebuild", *gens] if gens else [])
    if action == "sampling":
        gen = str(o.get("generation", "")).strip()
        args = ["sampling", exp] + (["--generation", gen] if _GEN.match(gen) else [])
        it = str(o.get("iterations", "")).strip()
        return args + (["--iterations", it] if it.isdigit() and int(it) > 0 else [])
    raise ValueError(f"unknown action {action}")


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
                            sampling=sorted((p.name for p in (prj / "sampling").glob("g*") if p.is_dir()), key=gen_key)))
        return out

    def project(self, exp: str) -> Path:
        if not exp or exp not in {e["name"] for e in self.experiments()}:
            raise ValueError(f"unknown experiment: {exp}")
        return self.workspace / ANALYSIS / exp


def generations(prj: Path) -> list[str]:
    root = prj / "generations"
    return sorted((g.name for g in root.glob("g*") if g.is_dir() and _GEN.match(g.name)), key=gen_key) if root.is_dir() else []


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
        sampling=sorted((p.name for p in (prj / "sampling").glob("g*") if p.is_dir()), key=gen_key),
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
    out: dict[str, dict[str, list]] = {}
    for r in rows:
        n, t, img = get(r, "n"), get(r, "treatment"), get(r, "image")
        if not img.isdigit():
            continue
        rel = f"review/{slug(n)}/{slug(t)}/{int(img):02d}.jpg"
        out.setdefault(n, {}).setdefault(t, []).append(dict(
            image=int(img), count=get(r, "auto_count"), low=get(r, "count_low"), high=get(r, "count_high"),
            qc=get(r, "qc_level"), flags=get(r, "qc_flags"), path=rel if (prj / rel).exists() else None))
    res = []
    for n, ts in out.items():
        res.append(dict(n=n, treatments=[dict(
            treatment=t, images=imgs,
            overview=(f"review/{slug(n)}/{slug(t)}_overview.jpg" if (prj / "review" / slug(n) / f"{slug(t)}_overview.jpg").exists() else None))
            for t, imgs in ts.items()]))
    return res


def save_layout(prj: Path, rows: list[dict]):
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
    and one empty manual column is always kept at the end (as `kerncount run` does)."""
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
    it = d.rglob("*") if recursive else d.iterdir()
    out = []
    for p in it:
        if p.name.startswith(".") or p.suffix in (".sig",):
            continue
        out.append(dict(path=p.relative_to(root).as_posix(), name=p.name, dir=p.is_dir(),
                        size=p.stat().st_size if p.is_file() else None))
    return sorted(out, key=lambda f: (not f["dir"], f["path"].lower()))


# ---------------------------------------------------------------- HTTP

def make_handler(state: State):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *a):  # keep the console quiet
            pass

        def _send(self, code: int, body, ctype="application/json", extra: dict | None = None):
            data = body if isinstance(body, bytes) else body.encode() if isinstance(body, str) else json.dumps(body, default=str).encode()
            self.send_response(code)
            self.send_header("Content-Type", ctype + ("; charset=utf-8" if ctype.startswith(("text/", "application/json")) else ""))
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            for k, v in (extra or {}).items():
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
                    p = inside(state.project(arg("exp")), arg("path"))
                    if not p.is_file():
                        return self._send(404, {"error": "not found"})
                    extra = {"Content-Disposition": f'attachment; filename="{p.name}"'} if arg("dl") else {}
                    return self._send(200, p.read_bytes(), mimetypes.guess_type(p.name)[0] or "application/octet-stream", extra)
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
                return table_data(p) if p.exists() else dict(columns=[], rows=[], total=0)
            if path == "/api/text":
                p = inside(prj, arg("path"))
                return dict(text=p.read_text(encoding="utf-8-sig", errors="replace")[:500_000] if p.is_file() else "")
            if path == "/api/files":
                return list_files(prj, arg("dir"), arg("recursive") == "1")
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
                if not state.workspace or not name or name in (".", "..") or re.search(r'[\\/:*?"<>|]', name):
                    raise ValueError("give a simple folder name (no / \\ : * ? \" < > |)")
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
                if target is None or not target.exists():
                    raise ValueError("not there yet (run the earlier steps first)")
                open_path(target)
                return f"opened {target.name}"
            if path == "/api/settings":
                self._idle()
                st = prj / "settings.toml"
                if "text" in body:
                    text = str(body["text"])
                    try:
                        tomllib.loads(text)
                    except tomllib.TOMLDecodeError as e:
                        raise ValueError(f"settings.toml is not valid: {e}") from None
                else:
                    text = settings_form.update_text(st.read_text(encoding="utf-8"), settings_form.coerce(body.get("values", {})))
                st.write_text(text, encoding="utf-8")
                return "settings saved"
            if path == "/api/settings_reset":
                self._idle()
                raw = settings_mod.load(prj)["project"]["raw_dir"]
                raw_abs = Path(raw) if Path(raw).is_absolute() else (prj / raw)
                settings_mod.write_template(prj, raw_abs)
                return "settings reset to defaults"
            if path == "/api/layout":
                self._idle()
                save_layout(prj, body.get("rows", []))
                return "layout saved"
            if path == "/api/layout_use_new":
                self._idle()
                new = prj / "layout.new.csv"
                if not new.exists():
                    raise ValueError("there is no fresh scan (layout.new.csv)")
                new.replace(prj / "layout.csv")
                return "layout.csv replaced by the fresh scan"
            if path == "/api/counts":
                self._idle()
                return dict(ok=True, message="corrections saved", **save_counts(prj, body.get("changes", []), body.get("mtime")))
            raise ValueError(f"unknown request {path}")

    return Handler


def serve(workspace: str | None = None, port: int = 8765, open_browser: bool = True):
    state = State(default_workspace(workspace))
    for p in (port, 0):  # fall back to any free port
        try:
            httpd = ThreadingHTTPServer(("127.0.0.1", p), make_handler(state))
            break
        except OSError:
            continue
    state.httpd = httpd
    url = f"http://127.0.0.1:{httpd.server_address[1]}/"
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
