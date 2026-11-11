"""Point-and-click front end: a small local web page that runs the kerncount commands.

Standard library only. Listens on 127.0.0.1 (this computer only); every action needs a per-session
token embedded in the page, so other websites cannot trigger runs. Commands run as a subprocess
(`python -m kerncount.cli ...`) in the workspace folder, one at a time, with the log streamed to the page.
"""
from __future__ import annotations

import json
import os
import re
import secrets
import subprocess
import sys
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from .cli import ANALYSIS, DATA, workspace_root

CONFIG = Path.home() / ".kerncount_gui.json"
TOKEN = secrets.token_urlsafe(16)
_GEN = re.compile(r"^g\d+$")


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


class State:
    def __init__(self, workspace: Path | None):
        self.lock = threading.Lock()
        self.job: Job | None = None
        self.workspace = workspace

    def experiments(self) -> list[dict]:
        ws = self.workspace
        if not ws or not (ws / DATA).is_dir():
            return []
        out = []
        for d in sorted(p for p in (ws / DATA).iterdir() if p.is_dir() and not p.name.startswith(".")):
            prj = ws / ANALYSIS / d.name
            gens = sorted((g.name for g in (prj / "generations").glob("g*") if _GEN.match(g.name)),
                          key=lambda g: int(g[1:])) if (prj / "generations").is_dir() else []
            out.append(dict(name=d.name, initialised=(prj / "settings.toml").exists(),
                            counted=(prj / "counts.csv").exists(), generations=gens,
                            sampling=(prj / "sampling").is_dir()))
        return out


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


def build_args(action: str, exp: str, o: dict) -> list[str]:
    if action == "init":
        return ["init", exp] + (["--force"] if o.get("force") else [])
    if action == "run":
        args = ["run", exp]
        if o.get("no_review"):
            args.append("--no-review")
        if o.get("resegment"):
            args.append("--resegment")
        rebuild = [g for g in str(o.get("rebuild", "")).replace(",", " ").split() if _GEN.match(g) or g == "all"]
        if rebuild:
            args += ["--rebuild", *rebuild]
        return args
    if action == "sampling":
        gen = str(o.get("generation", "")).strip()
        return ["sampling", exp] + (["--generation", gen] if _GEN.match(gen) else [])
    raise ValueError(f"unknown action {action}")


def make_handler(state: State):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *a):  # keep the console quiet
            pass

        def _send(self, code: int, body, ctype="application/json"):
            data = body.encode() if isinstance(body, str) else json.dumps(body).encode()
            self.send_response(code)
            self.send_header("Content-Type", ctype + "; charset=utf-8")
            self.send_header("Cache-Control", "no-store")
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

        def do_GET(self):
            if not self._local():
                return
            if self.path == "/":
                return self._send(200, PAGE.replace("__TOKEN__", TOKEN), "text/html")
            if self.path == "/api/state":
                j = state.job
                return self._send(200, dict(
                    workspace=str(state.workspace or ""), experiments=state.experiments(),
                    job=j and dict(label=j.label, running=j.running, returncode=j.returncode)))
            if self.path.startswith("/api/log"):
                since = int(self.path.partition("since=")[2] or 0)
                j = state.job
                lines = j.log[since:] if j else []
                return self._send(200, dict(lines=lines, next=since + len(lines)))
            self._send(404, {"error": "not found"})

        def do_POST(self):
            if not self._local():
                return
            if self.headers.get("X-Token") != TOKEN:
                return self._send(403, {"error": "bad token - reload the page"})
            try:
                body = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
                msg = self._post(self.path, body)
                self._send(200, {"ok": True, "message": msg})
            except Exception as e:  # report every problem on the page
                self._send(400, {"error": str(e)})

        def _exp(self, body) -> str:
            name = body.get("experiment", "")
            if name not in {e["name"] for e in state.experiments()}:
                raise ValueError(f"unknown experiment: {name}")
            return name

        def _post(self, path, body):
            if path == "/api/workspace":
                p = Path(body.get("path", "").strip().strip('"')).expanduser()
                if not (p / DATA).is_dir():
                    raise ValueError(f"{p} has no data folder inside it")
                state.workspace = p.resolve()
                CONFIG.write_text(json.dumps({"workspace": str(state.workspace)}))
                return "workspace set"
            if path == "/api/stop":
                if state.job and state.job.running:
                    state.job.proc.terminate()
                return "stopping"
            if path == "/api/start":
                exp = self._exp(body)
                args = build_args(body.get("action", ""), exp, body.get("options", {}))
                with state.lock:
                    if state.job and state.job.running:
                        raise ValueError("another job is still running - wait or press Stop")
                    state.job = Job(f"{args[0]} {exp}", args, state.workspace)
                return "started"
            if path == "/api/open":
                prj = state.workspace / ANALYSIS / self._exp(body)
                what = body.get("what")
                gens = sorted((prj / "generations").glob("g*"), key=lambda g: int(g.name[1:]) if _GEN.match(g.name) else -1)
                target = {"folder": prj, "counts": prj / "counts.csv", "settings": prj / "settings.toml",
                          "layout": prj / "layout.csv", "review": prj / "review",
                          "results": gens[-1] / "results.xlsx" if gens else None,
                          "data": state.workspace / DATA / body["experiment"]}.get(what)
                if target is None or not target.exists():
                    raise ValueError(f"not there yet: {what} (run the earlier steps first)")
                open_path(target)
                return f"opened {target.name}"
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
    url = f"http://127.0.0.1:{httpd.server_address[1]}/"
    print(f"kerncount is open in your browser at {url}\nKeep this window open while you work; close it to quit.", flush=True)
    if open_browser:
        threading.Timer(0.5, webbrowser.open, (url,)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        if state.job and state.job.running:
            state.job.proc.terminate()


PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>kerncount</title>
<style>
 body{font:15px system-ui,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;color:#222}
 h1{font-size:22px;margin:0 0 4px} .muted{color:#666;font-size:13px}
 table{border-collapse:collapse;width:100%;margin:12px 0}
 td,th{border-bottom:1px solid #ddd;padding:8px 6px;text-align:left;vertical-align:top}
 button{font:inherit;font-size:13px;padding:4px 10px;margin:2px;border:1px solid #888;border-radius:4px;background:#f6f6f6;cursor:pointer}
 button.main{background:#2563eb;color:#fff;border-color:#2563eb} button:disabled{opacity:.4;cursor:default}
 .tag{display:inline-block;font-size:12px;padding:1px 6px;border-radius:8px;background:#eee;margin-right:4px}
 .ok{background:#dcfce7} .no{background:#fee2e2}
 pre{background:#111;color:#ddd;padding:10px;height:340px;overflow:auto;font-size:12px;white-space:pre-wrap}
 details{margin:8px 0} input[type=text]{font:inherit;padding:3px 6px}
 #msg{min-height:20px;color:#b45309}
</style></head><body>
<h1>kerncount</h1>
<div class="muted">Steps: <b>1 Set up</b> (once) → check <i>layout</i> → <b>2 Run</b> → look at <i>review images</i> →
optionally correct <i>counts.csv</i> in Excel (save &amp; close it) → <b>Run</b> again → <i>results</i>.</div>
<p>Workspace: <input type="text" id="ws" size="70"> <button onclick="setWs()">Use this folder</button>
<span class="muted">(the folder containing <code>data</code>)</span></p>
<details><summary>Options</summary>
 <p><b>Run:</b> <label><input type="checkbox" id="no_review"> skip review images (faster)</label>
 <label><input type="checkbox" id="resegment"> re-detect every nucleus (ignore cache)</label>
 <label>overwrite generation(s): <input type="text" id="rebuild" size="8" placeholder="e.g. g2"></label></p>
 <p><b>Sampling:</b> <label>generation <input type="text" id="generation" size="5" placeholder="latest"></label></p>
 <p><b>Set up:</b> <label><input type="checkbox" id="force"> reset settings.toml and layout.csv</label></p>
</details>
<table><thead><tr><th>Experiment</th><th>Status</th><th>Actions</th><th>Open</th></tr></thead><tbody id="rows"></tbody></table>
<div id="msg"></div>
<p><b id="joblabel">Log</b> <button id="stop" onclick="post('/api/stop',{})" disabled>Stop</button></p>
<pre id="log"></pre>
<script>
const TOKEN="__TOKEN__"; let next=0, running=false, lastLabel=null;
const $=id=>document.getElementById(id);
async function post(url,body){
  const r=await fetch(url,{method:"POST",headers:{"X-Token":TOKEN,"Content-Type":"application/json"},body:JSON.stringify(body)});
  const j=await r.json(); $("msg").textContent=j.error||""; refresh(); return j;
}
function opts(){return{no_review:$("no_review").checked,resegment:$("resegment").checked,rebuild:$("rebuild").value,
  generation:$("generation").value,force:$("force").checked}}
function start(action,exp){next=0;$("log").textContent="";post("/api/start",{action,experiment:exp,options:opts()})}
function openIt(exp,what){post("/api/open",{experiment:exp,what})}
function setWs(){post("/api/workspace",{path:$("ws").value})}
const esc=s=>s.replace(/[&<>"']/g,c=>"&#"+c.charCodeAt(0)+";");
function row(e){
  const n=esc(JSON.stringify(e.name)), d=running?"disabled":"";
  const st=(e.initialised?'<span class="tag ok">set up</span>':'<span class="tag no">not set up</span>')+
    (e.generations.length?`<span class="tag ok">${e.generations.length} generation(s)</span>`:'');
  const b=(lbl,act,cls,need)=>`<button class="${cls||""}" ${d} ${need===false?"disabled":""} onclick='start("${act}",${n})'>${lbl}</button>`;
  const o=(lbl,w,need)=>`<button ${need?"":"disabled"} onclick='openIt(${n},"${w}")'>${lbl}</button>`;
  return `<tr><td><b>${esc(e.name)}</b></td><td>${st}</td><td>`+
    b("1 Set up","init",e.initialised?"":"main")+b("2 Run","run",e.initialised?"main":"",e.initialised)+
    b("Sampling","sampling","",e.generations.length>0)+`</td><td>`+
    o("Images","data",true)+o("Project folder","folder",e.initialised)+o("layout","layout",e.initialised)+
    o("settings","settings",e.initialised)+o("Review images","review",e.counted)+o("counts.csv","counts",e.counted)+
    o("Results (Excel)","results",e.generations.length>0)+`</td></tr>`;
}
async function refresh(){
  const s=await (await fetch("/api/state")).json();
  if(document.activeElement!==$("ws")) $("ws").value=s.workspace;
  running=!!(s.job&&s.job.running); $("stop").disabled=!running;
  $("joblabel").textContent=s.job?`Log: ${s.job.label} ${running?"(running…)":s.job.returncode===0?"(done)":"(stopped/failed)"}`:"Log";
  $("rows").innerHTML=s.experiments.length?s.experiments.map(row).join(""):
    '<tr><td colspan="4">No experiments found. Enter the workspace folder above (it must contain a <code>data</code> folder with one folder per experiment).</td></tr>';
}
async function poll(){
  try{
    const j=await (await fetch("/api/log?since="+next)).json();
    if(j.lines.length){const L=$("log"),end=L.scrollTop+L.clientHeight>=L.scrollHeight-5;
      L.textContent+=j.lines.join("\n")+"\n"; next=j.next; if(end) L.scrollTop=L.scrollHeight;}
    const was=running; await refresh(); if(was&&!running) refresh();
  }catch(e){$("msg").textContent="Lost contact with kerncount - is its window still open?"}
  setTimeout(poll,1000);
}
poll();
</script></body></html>
"""
