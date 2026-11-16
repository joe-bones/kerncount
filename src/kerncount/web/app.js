"use strict";
// kerncount web interface. Plain JavaScript, no build step, no external libraries (works offline).
// Talks to gui.py over /api/* (token in X-Token / ?t=). Files: app.js (shell, experiment home, set-up, images),
// results.js (results explorer, statistics table, layer comparison), intensity.js (intensity page, blinded
// cut-off tuning), plan.js (planner, plate maps). Navigation state lives in the URL hash:
//   #<experiment>/<section>?key=value&...   (deep links; also used by headless screenshots)

const S = { state: null, exp: null, sec: "home", params: {}, project: null, ov: null, logNext: 0, running: false,
            counts: null, round: null, pending: {}, sel: null, review: null, reviewN: null, reviewT: null, reviewQ: "all",
            filters: { n: "", t: "", q: "", only: false, text: "" }, dir: "", man: {}, stats: {}, renderTok: 0,
            lbOverlay: false };
const $ = id => document.getElementById(id);
const esc = s => String(s ?? "").replace(/[&<>"']/g, c => "&#" + c.charCodeAt(0) + ";");
const js = v => esc(JSON.stringify(v));  // a value inside an onclick='...' attribute
const qs = o => Object.entries(o).map(([k, v]) => encodeURIComponent(k) + "=" + encodeURIComponent(v)).join("&");
const fileUrl = (path, dl) => "/file?" + qs({ exp: S.exp, path, t: TOKEN, ...(dl ? { dl: 1 } : {}) });
const genNum = g => parseInt(String(g).slice(1), 10);
const shortN = n => String(n || "").split("_")[0];
const baseName = p => String(p || "").split(/[\\/]/).pop();
const tail2 = p => String(p || "").split(/[\\/]/).slice(-2).join("/");
const pct = (x, d = 1) => x == null || x === "" || isNaN(+x) ? "–" : (100 * +x).toFixed(d) + "%";
const N_FALLBACK = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f"];

async function get(path, params = {}) {
  const r = await fetch(path + "?" + qs({ ...params, t: TOKEN }));
  const j = await r.json();
  if (!r.ok) throw new Error(j.error || r.statusText);
  return j;
}
async function post(path, body = {}, quiet = false) {
  const r = await fetch(path, { method: "POST", headers: { "X-Token": TOKEN, "Content-Type": "application/json" }, body: JSON.stringify(body) });
  const j = await r.json();
  if (!r.ok) { toast(j.error || r.statusText, true); throw new Error(j.error); }
  if (!quiet && j.message) toast(j.message);
  return j;
}
function toast(msg, err = false) {
  const d = document.createElement("div");
  d.className = "t" + (err ? " err" : "");
  d.textContent = msg;
  $("toast").appendChild(d);
  setTimeout(() => d.remove(), err ? 9000 : 3500);
}
// in-page confirmation (the browser dialogs are never used)
function confirmBox(html, onYes, yesLabel = "Yes, go ahead", noLabel = "Cancel") {
  openModal(`<div style="margin-bottom:14px">${html}</div>
    <div class="row"><button class="primary" id="cfYes">${esc(yesLabel)}</button><button onclick="closeModal()">${esc(noLabel)}</button></div>`, "narrow");
  $("cfYes").onclick = () => { closeModal(); onYes(); };
}

// ---------------------------------------------------------------- navigation: experiment -> sections, URL hash

const SECTIONS = [
  ["Overview", [["home", "Experiment home"]]],
  ["Results", [["results", "Results explorer"], ["stats", "Statistics table"], ["layers", "Layer comparison", e => e.intensity],
               ["intensity", "Intensity & cut-off", e => e.intensity]]],
  ["Images", [["review", "Review images"], ["corrections", "Corrections"]]],
  ["Plan the next experiment", [["planner", "Planner"], ["plates", "Plate maps"]]],
  ["Set-up", [["layout", "Layout"], ["settings", "Settings"], ["files", "All files"]]],
];
const SEC_TITLE = Object.fromEntries(SECTIONS.flatMap(([, l]) => l.map(([k, t]) => [k, t])));
const OLD_TABS = { run: "home", sampling: "planner" };

function hashFor(exp, sec, params) {
  if (!exp) return "#";
  const p = Object.entries(params || {}).filter(([, v]) => v !== "" && v != null);
  return "#" + encodeURIComponent(exp) + "/" + sec + (p.length ? "?" + qs(Object.fromEntries(p)) : "");
}
function parseHash() {
  const h = location.hash.slice(1), [path, query] = h.split("?");
  const [exp, sec] = path.split("/").map(decodeURIComponent);
  const params = {};
  if (query) for (const kv of query.split("&")) { const [k, v = ""] = kv.split("="); if (k) params[decodeURIComponent(k)] = decodeURIComponent(v.replace(/\+/g, " ")); }
  return { exp: exp || null, sec: OLD_TABS[sec] || sec || "home", params };
}
function nav(sec, params = {}, exp = S.exp) {
  if (S.sec === "corrections" && sec !== "corrections" && Object.keys(S.pending).length) {
    return confirmBox("You have unsaved corrections. Leave them unsaved?", () => { S.pending = {}; nav(sec, params, exp); }, "Leave without saving");
  }
  const h = hashFor(exp, sec, params);
  if (location.hash === h) applyHash(); else location.hash = h;
}
function setP(obj, rerender = true) {  // change parameters of the current section (no new history entry)
  for (const [k, v] of Object.entries(obj)) { if (v === "" || v == null) delete S.params[k]; else S.params[k] = String(v); }
  history.replaceState(null, "", hashFor(S.exp, S.sec, S.params));
  if (rerender) render();
}
async function applyHash() {
  const h = parseHash();
  if (h.exp !== S.exp) resetExp();
  S.exp = h.exp; S.sec = SEC_TITLE[h.sec] ? h.sec : "home"; S.params = h.params;
  if (S.sec === "review" && h.params.q) S.reviewQ = h.params.q;
  await refreshState();
  render();
}
function resetExp() {
  S.project = null; S.ov = null; S.counts = null; S.review = null; S.pending = {}; S.round = null; S.dir = ""; S.man = {}; S.stats = {};
  S.reviewN = null; S.reviewT = null; S.tune = null; S.plan = null; S.plans = null; S.intens = null;
}
function go(exp, sec) { nav(OLD_TABS[sec] || sec || "home", {}, exp); }

async function refreshState() {
  S.state = await get("/api/state");
  $("wsLabel").textContent = S.state.workspace || "no workspace chosen";
  const j = S.state.job, was = S.running;
  S.running = !!(j && j.running);
  $("stopBtn").disabled = !S.running;
  $("jobLabel").innerHTML = j ? `<span class="dot ${S.running ? "run" : j.returncode === 0 ? "" : "fail"}"></span>` +
    esc(j.label) + (S.running ? " — running…" : j.returncode === 0 ? " — done" : " — stopped / failed") : "Log";
  renderSide();
  if (was && !S.running) onJobFinished();
  return S.state;
}
function renderSide() {
  if (!S.state) return;
  const html = S.state.experiments.map(e => `<button class="exp ${e.name === S.exp ? "active" : ""}" onclick='go(${js(e.name)})'>
      <span class="nm">${esc(e.name)}</span>
      <span class="st">${!e.initialised ? "not set up" : !e.generations.length ? "set up, not run" : e.generations.length + " generation(s)"}${e.intensity ? " · intensity" : ""}</span></button>
      ${e.name === S.exp ? `<div class="secnav">${SECTIONS.map(([g, items]) => {
        const its = items.filter(([, , show]) => !show || show(e));
        return its.length ? `<div class="grp">${esc(g)}</div>` + its.map(([k, t]) => `<a class="${k === S.sec ? "active" : ""} ${!e.initialised && k !== "home" ? "dis" : ""}" onclick="nav('${k}')">${esc(t)}</a>`).join("") : "";
      }).join("")}</div>` : ""}`).join("")
    || `<div class="muted small" style="padding:0 8px">No experiments yet.</div>`;
  if ($("expList").dataset.h !== html) { $("expList").innerHTML = html; $("expList").dataset.h = html; }
}
function expInfo() { return S.state && S.state.experiments.find(e => e.name === S.exp); }

const RENDER = {};  // section -> async (body, exp info) ; filled below and in the other files
async function render() {
  const tok = ++S.renderTok, m = $("main");
  renderSide();
  if (!S.exp) return renderWelcome(m);
  const e = expInfo();
  if (!e) { S.exp = null; return renderWelcome(m); }
  if (!e.initialised && S.sec !== "home") S.sec = "home";
  if (!S.project && e.initialised) {
    try { S.project = await get("/api/project", { exp: S.exp }); } catch (err) { toast(err.message, true); }
    if (tok !== S.renderTok) return;
  }
  m.innerHTML = `<div class="crumb">${esc(S.exp)} › <b>${esc(SEC_TITLE[S.sec])}</b></div><div id="tabBody"></div>`;
  const body = $("tabBody");
  try { await RENDER[S.sec](body, e, tok); }
  catch (err) { if (tok === S.renderTok) body.innerHTML = `<div class="notice err">${esc(err.message)}</div>`; console.error(err); }
}

function onJobFinished() {
  S.project = null; S.counts = null; S.review = null; S.ov = null; S.man = {}; S.stats = {}; S.intens = null; S.plan = null;
  if (S.sec === "corrections" && Object.keys(S.pending).length) return; // don't wipe unsaved edits
  if (!$("modal").classList.contains("hidden")) return;
  if (S.sec === "intensity" && S.params.tab === "tune") return; // keep a tuning session
  render();
}

async function pollLog() {
  try {
    const j = await get("/api/log", { since: S.logNext });
    if (j.lines.length) {
      const L = $("log"), end = L.scrollTop + L.clientHeight >= L.scrollHeight - 5;
      L.textContent += j.lines.join("\n") + "\n";
      S.logNext = j.next;
      if (end) L.scrollTop = L.scrollHeight;
    }
    await refreshState();
  } catch (e) {
    $("jobLabel").innerHTML = `<span class="dot fail"></span>Lost contact with kerncount — is its window still open?`;
  }
  setTimeout(pollLog, 1000);
}
function toggleLog(ev) { if (ev && ev.target.tagName === "BUTTON") return; $("logDock").classList.toggle("collapsed"); $("logChev").textContent = $("logDock").classList.contains("collapsed") ? "▲" : "▼"; }
function showLog() { $("logDock").classList.remove("collapsed"); $("logChev").textContent = "▼"; }
async function stopJob(ev) { ev.stopPropagation(); await post("/api/stop"); }
async function startJob(action, options = {}) {
  S.logNext = 0; $("log").textContent = "";
  await post("/api/start", { action, experiment: S.exp, options }, true);
  showLog();
  toast(`${action} started — progress in the log below`);
  await refreshState();
}
async function quitApp() {
  confirmBox("Stop kerncount? Any running job is cancelled. You can start it again from the kerncount program.", async () => {
    await post("/api/quit");
    document.body.innerHTML = `<div style="padding:60px;font:16px system-ui">kerncount has stopped. You can close this tab.</div>`;
  }, "Quit kerncount");
}
async function newExperiment() {
  const name = $("newExp").value.trim();
  if (!name) return toast("Type a name for the experiment first", true);
  await post("/api/new_experiment", { name });
  $("newExp").value = "";
  await refreshState();
  go(name, "home");
}
async function openSys(what, path) { await post("/api/open", { experiment: S.exp, what, path }); }

// colours of the n's: the same as in the plots (n_colours.csv of the generation)
function nColour(n) {
  if (S.project && S.project.colours && S.project.colours[n]) return S.project.colours[n];
  const c = (S.ov && S.ov.ns.find(x => x.n === n)) || null;
  if (c && c.colour) return c.colour;
  for (const g of Object.values(S.man)) if (g.colours && g.colours[n]) return g.colours[n];
  const all = (S.ov ? S.ov.ns.map(x => x.n) : []).sort(), i = all.indexOf(n);
  return N_FALLBACK[(i < 0 ? 0 : i) % N_FALLBACK.length];
}
const nChipShort = n => n ? `<span class="nchip" title="${esc(n)}"><span class="ndot" style="background:${esc(nColour(n))}"></span>${esc(shortN(n))}</span>` : "";
const nChip = n => n ? `<span class="nchip" title="${esc(n)}"><span class="ndot" style="background:${esc(nColour(n))}"></span>${esc(n)}</span>` : "";
async function loadOverview(force) {
  if (!S.ov || force) S.ov = await get("/api/overview", { exp: S.exp });
  return S.ov;
}

// ---------------------------------------------------------------- welcome: workspace

function renderWelcome(m) {
  const ws = S.state ? S.state.workspace : "";
  m.innerHTML = `<div class="hero">
    <h1>Welcome to kerncount</h1>
    <p class="muted">Automated nucleus counting, fluorescence readouts, review and statistics. Everything happens in this browser tab; nothing leaves this computer.</p>
    <div class="card"><h3>Workspace folder</h3>
      <p class="muted small">The folder that holds <code>data\\</code> (your raw images, one folder per experiment) and
      <code>analysis\\</code> (all results). Results are written next to your data, never into it.</p>
      <div class="row"><input id="wsIn" style="flex:1;min-width:300px" value="${esc(ws)}" placeholder="e.g. C:\\Users\\me\\Documents\\kerncount-work">
      <button onclick="browseWs()">Browse…</button><button class="primary" onclick="setWs(false)">Use this folder</button></div>
      <div class="small" style="margin-top:6px"><a href="#" onclick="setWs(true);return false">Create a new workspace here</a> (makes an empty <code>data</code> folder)</div>
    </div>
    ${ws ? `<div class="card"><h3>How it works</h3><ol>
      <li>Put each experiment's images in <code>data\\&lt;experiment&gt;\\&lt;one folder per n&gt;\\&lt;one folder per treatment&gt;\\1.tif …</code>
          (or use <b>+ Add</b> on the left to create the folder).</li>
      <li>Pick the experiment on the left → <b>Set up</b>, check <b>Layout</b> and <b>Settings</b>.</li>
      <li><b>Run</b> the analysis (first run ~2–3 s per image; later runs reuse the results).</li>
      <li>Look through <b>Review images</b>; type corrections in <b>Corrections</b> → a new <i>generation</i> (g2, g3, …) is built.</li>
      <li>Read plots and statistics in the <b>Results explorer</b> and the <b>Statistics table</b>; plan the next experiment in the <b>Planner</b>.</li></ol></div>` : ""}
  </div>`;
}
async function browseWs() {
  toast("A folder picker opened (it may be behind this window)");
  const r = await post("/api/browse", {}, true);
  if (r.path) $("wsIn").value = r.path;
}
async function setWs(create) {
  await post("/api/workspace", { path: $("wsIn").value, create });
  await refreshState(); render();
}

// ---------------------------------------------------------------- experiment home

function evidenceBtn(ev) {
  if (ev.type === "table") return `<button class="small" onclick='viewTable(${js(ev.path)})'>${esc(ev.label)}</button>`;
  if (ev.type === "image") return `<button class="small" onclick='openImage(${js(ev.path)})'>${esc(ev.label)}</button>`;
  if (ev.type === "go") return `<button class="small" onclick='nav(${js(ev.section)}, ${js(ev.params || {})})'>${esc(ev.label)} →</button>`;
  return "";
}
function metaVal(v, unit) {
  if (v === "" || v == null || (Array.isArray(v) && !v.length)) return `<div class="v unset">not set</div>`;
  return `<div class="v">${esc(Array.isArray(v) ? v.join(", ") : v)}${unit ? " " + esc(unit) : ""}</div>`;
}
RENDER.home = async function (b, e, tok) {
  if (!e.initialised) { b.innerHTML = runPanel(e); return; }
  const ov = await loadOverview();
  if (tok !== S.renderTok) return;
  const m = ov.meta, d = ov.design;
  const chans = m.channel_files === "interleaved"
    ? `separate files, alternating ${(m.file_order || []).join(" → ")}; nuclei found on ${m.segment}`
    : (m.channel_names && m.channel_names.length ? m.channel_names.join(", ") : "one file per field") + (m.segment ? `; nuclei on channel ${m.segment}` : "");
  const att = ov.attention;
  b.innerHTML = `
    <div class="card exp-head">
      <div class="title"><h1>${esc(S.exp)}</h1>
        <div class="muted small">${esc(S.project ? S.project.path : "")}</div>
        <div style="margin-top:6px">${(m.cell_lines.length ? m.cell_lines : [m.cell_line]).filter(Boolean).map(c => `<span class="tag blue">${esc(c)}</span>`).join("")}
          ${ov.intensity ? `<span class="tag purple">fluorescence: ${esc((m.marker_channels || []).join(", "))}</span>` : `<span class="tag">nucleus counts</span>`}
          ${ov.n_count != null ? `<span class="tag">${ov.n_count} n</span>` : ""}</div></div>
      <div style="flex:2;min-width:420px"><div class="kv">
        <div><div class="k">Cell line</div>${metaVal(m.cell_line)}</div>
        <div><div class="k">Magnification</div>${metaVal(m.magnification ? m.magnification + "x" : "")}</div>
        <div><div class="k">Seeding density</div>${metaVal(m.seeding_density, m.seeding_density ? m.seeding_unit : "")}</div>
        <div><div class="k">Stain / antibodies</div>${metaVal(m.stain)}</div>
        <div><div class="k">Channels</div>${metaVal(chans)}</div>
        <div><div class="k">Negative control</div>${metaVal(m.negative_control ? `${m.negative_control} (${m.negative_control_type.replace("_", " ")})` : "")}</div>
        <div><div class="k">Main reference</div>${metaVal(m.reference)}</div>
        <div><div class="k">Controls</div>${metaVal(m.controls)}</div>
      </div><div class="small" style="margin-top:8px"><a class="link" onclick="nav('settings',{sec:'experiment'})">Edit experiment details →</a>${m.notes ? ` · <span class="muted">${esc(m.notes)}</span>` : ""}</div></div>
    </div>
    <h2>Attention ${att.length ? `<span class="tag ${att.some(a => a.level === "error") ? "fail" : att.some(a => a.level === "warn") ? "medium" : "blue"}">${att.length}</span>` : ""}</h2>
    ${att.length ? att.map(a => `<div class="att ${a.level}"><div class="ico">${a.level === "error" ? "!" : a.level === "warn" ? "!" : "i"}</div>
        <div><div class="t1">${esc(a.title)}</div><div class="t2">${esc(a.text)}</div><div class="ev">${a.evidence.map(evidenceBtn).join("")}</div></div></div>`).join("")
      : `<div class="notice ok">Nothing needs your attention: no duplicates, failed checks or design problems were found.</div>`}
    <div class="cols2" style="margin-top:14px">
      <div class="card"><h3>Experiments (n) and timepoints</h3>
        <table class="grid"><thead><tr><th>n</th><th>date</th><th>timepoint</th><th class="num">conditions</th><th class="num">images</th><th class="num">QC low</th><th>acceptance</th></tr></thead>
        <tbody>${ov.ns.map(n => `<tr><td>${nChip(n.n)}</td><td class="nw">${esc(n.date)}</td><td>${esc(n.timepoints.join(", ")) || `<span class="muted">–</span>`}</td>
          <td class="num">${n.conditions.length}</td><td class="num">${n.qc_images ?? n.images}</td>
          <td class="num">${n.qc_images ? pct(n.qc_low / n.qc_images, 0) : "–"}</td>
          <td>${n.acceptance ? `<span class="tag ${esc(n.acceptance)}" title="${esc(n.acceptance_detail)}">${esc(n.acceptance)}</span>` : `<span class="muted small">not checked</span>`}</td></tr>`).join("")}</tbody></table>
        <p class="small muted" style="margin-bottom:0">Dot colour = the colour of that n in every plot. Timepoints come from the n-folder names (edit them in <a class="link" onclick="nav('layout')">Layout</a>).</p></div>
      <div class="card"><h3>Design</h3>
        <p style="margin:0 0 6px"><b>${d.conditions}</b> conditions · <b>${esc(d.wells || "?")}</b> wells × <b>${esc(d.fields_per_well || "?")}</b> fields per condition per n
          <span class="muted small">(${esc(d.fields_source || "")}; ${esc(d.layout)} numbering; ${esc(d.positions)} field positions)</span></p>
        <p style="margin:0 0 8px">${d.images} images in total. ${d.flagged ? `<span class="tag medium">${d.flagged} condition(s) flagged</span>` : `<span class="tag ok">no design flags</span>`}
          ${d.unbalanced ? `<span class="tag medium">${d.unbalanced} unbalanced</span>` : ""}</p>
        <div class="row"><button class="small" onclick="viewTable('design_check.csv')">design_check.csv</button>
          ${ov.wells_csv ? `<button class="small" onclick="viewTable('wells.csv')">wells.csv</button>` : ""}
          <button class="small" onclick="nav('settings',{sec:'design'})">Design settings</button></div>
        ${ov.intensity && ov.cutoffs.length ? `<h3>Positivity cut-offs</h3>${ov.cutoffs.filter(c => c.measure === "mean").map(c => `<div class="small">${nChip(c.n)} ${esc(c.channel)}: <b>${esc(c.rule_label)}</b> = ${(+c.value).toFixed(1)} → ${pct(c.frac_positive)} positive</div>`).join("")}
          <div style="margin-top:6px"><button class="small" onclick="nav('intensity')">Intensity & cut-off →</button></div>` : ""}
      </div>
    </div>
    <div class="card"><h3>Generations</h3>
      ${ov.generations.length ? `<table class="grid"><thead><tr><th>generation</th><th>built from</th><th>built</th><th class="num">n</th><th>stars from</th><th></th></tr></thead><tbody>
        ${ov.generations.map(g => `<tr><td><b>${esc(g.name)}</b> ${g.name === ov.latest ? `<span class="tag blue">latest</span>` : ""}</td>
          <td>${g.source === "auto" ? "automatic counts" : `previous generation + ${esc(g.changes)} manual value(s) in <code>${esc(g.source)}</code>`}
            ${g.rebuilt_because ? `<div class="small muted">rebuilt because: ${esc([].concat(g.rebuilt_because).join("; "))}</div>` : ""}</td>
          <td class="small">${esc((g.built || "").replace("T", " "))}</td><td class="num">${g.n_count ?? "–"}</td>
          <td class="small">${esc(g.primary === "mixed" ? "mixed model" : g.primary === "rm_anova" ? "RM ANOVA" : g.primary || "–")}${g.reference ? ` vs ${esc(g.reference)}` : ""}</td>
          <td style="white-space:nowrap">${g.has_manifest ? `<button class="small" onclick="nav('results',{gen:'${g.name}'})">Results</button> <button class="small" onclick="nav('stats',{gen:'${g.name}'})">Statistics</button>` : `<span class="muted small">built by an older version</span>`}
            ${g.has_xlsx ? `<a class="btn small" href="${fileUrl(`generations/${g.name}/results.xlsx`, 1)}">Excel</a>` : ""}</td></tr>`).join("")}</tbody></table>`
        : `<div class="muted">No generations yet: run the analysis.</div>`}
      ${ov.notes.length ? `<div class="small muted" style="margin-top:8px">Notes of ${esc(ov.latest)}: ${ov.notes.map(esc).join(" · ")}</div>` : ""}
    </div>
    <h2 id="run">Run</h2>${runPanel(e)}`;
  if (S.params.anchor) { const a = $(S.params.anchor); if (a) a.scrollIntoView(); }
};

function runPanel(e) {
  const gens = e.generations, last = gens[gens.length - 1];
  const step = (n, done, title, text, controls) => `<div class="step ${done ? "done" : ""}"><div><span class="num">${done ? "✓" : n}</span><b>${title}</b></div><p>${text}</p>${controls}</div>`;
  const d = S.running ? "disabled" : "";
  return `
    ${S.running ? `<div class="notice info">A job is running — see the log at the bottom. Buttons that change files are paused until it finishes.</div>` : ""}
    <div class="steps">
      ${step(1, e.initialised, "Set up", "Scan the raw image folders and create the project (layout + settings).",
        `<label class="small"><input type="checkbox" id="o_force"> start over: reset settings and layout</label>
         <label class="small"><input type="checkbox" id="o_forcelayout"> rewrite the layout only (after changing how channels are saved)</label>
         <button class="${e.initialised ? "" : "primary"}" ${d} onclick="startJob('init',{force:$('o_force').checked,force_layout:$('o_forcelayout').checked})">${e.initialised ? "Re-scan folders" : "Set up"}</button>`)}
      ${step(2, gens.length > 0, "Run the analysis", "Count nuclei, measure intensities, score quality, draw review images and build statistics. Always safe to re-run.",
        `<label><input type="checkbox" id="o_noreview"> skip review images (faster)</label>
         <label><input type="checkbox" id="o_reseg"> re-detect every nucleus (ignore cache)</label>
         <label class="small">overwrite generation(s): <input id="o_rebuild" size="8" placeholder="e.g. g2 or all"></label>
         <button class="primary" ${e.initialised ? "" : "disabled"} ${d} onclick="startJob('run',{no_review:$('o_noreview').checked,resegment:$('o_reseg').checked,rebuild:$('o_rebuild').value})">Run analysis</button>`)}
      ${step(3, false, "Rebuild generations", "Statistics only, from counts.csv (after corrections or statistics settings). Fast; no re-detection.",
        `<label class="small">overwrite generation(s): <input id="o_rebuild2" size="8" placeholder="e.g. g2 or all"></label>
         <button ${e.counted ? "" : "disabled"} ${d} onclick="startJob('regen',{rebuild:$('o_rebuild2').value})">Rebuild</button>
         <button class="small" ${e.counted ? "" : "disabled"} onclick="nav('corrections')">Corrections</button>`)}
      ${step(4, e.sampling.length > 0, "Sampling and planner", "Variance components, power and subsampling on a generation: feeds the Planner.",
        `<label class="small">generation <select id="o_sgen">${gens.map(g => `<option ${g === last ? "selected" : ""}>${esc(g)}</option>`).join("")}</select>
           subsamples <input id="o_sit" value="1000" size="5"></label>
         <button ${last ? "" : "disabled"} ${d} onclick="startJob('sampling',{generation:$('o_sgen').value,iterations:$('o_sit').value})">Run sampling</button>
         <button class="small" ${e.sampling.length ? "" : "disabled"} onclick="nav('planner')">Planner</button>`)}
      ${step(5, false, "Well template", "Write wells.csv (which well and position each image is) to edit by hand for unusual plate designs.",
        `<label class="small"><input type="checkbox" id="o_wforce"> overwrite an existing wells.csv</label>
         <button ${e.initialised ? "" : "disabled"} ${d} onclick="startJob('wells-template',{force:$('o_wforce').checked})">Write wells template</button>`)}
      ${step(6, false, "Plate maps", "Randomised, balanced plate maps for the next experiment (one plate per n).",
        `<button ${e.initialised ? "" : "disabled"} onclick="nav('plates')">Plate maps →</button>`)}
    </div>
    <h3>On this computer</h3>
    <div class="row"><button onclick="openSys('data')">Open raw images folder</button>
      <button ${e.initialised ? "" : "disabled"} onclick="openSys('folder')">Open project folder</button>
      <button ${last ? "" : "disabled"} onclick="openSys('results')">Open latest results in Excel</button></div>`;
}

// ---------------------------------------------------------------- Layout

RENDER.layout = async function (b) {
  const t = await get("/api/table", { exp: S.exp, path: "layout.csv" });
  const ix = Object.fromEntries(t.columns.map((c, i) => [c, i]));
  const hasNew = S.project && S.project.files["layout.new.csv"];
  b.innerHTML = `
    ${hasNew ? `<div class="notice">A fresh scan of the raw folders was saved as <code>layout.new.csv</code> (the folders changed since set-up).
       <button class="small" onclick="viewTable('layout.new.csv')">Compare</button>
       <button class="small" onclick="useNewLayout()">Use the fresh scan</button></div>` : ""}
    <p class="muted">One row per raw folder. Fix the <b>n</b>, <b>timepoint</b> or <b>treatment</b> name if it was read wrongly; untick <b>include</b> to leave a folder out.
      Then run the analysis again.</p>
    <div class="bar"><button class="primary" id="saveLayout" onclick="saveLayout()">Save layout</button><span class="muted small" id="layoutState"></span></div>
    <div class="tablewrap"><table class="grid" id="layoutTbl"><thead><tr><th>n folder</th><th>treatment folder</th><th>n</th><th>date</th>${"timepoint" in ix ? "<th>timepoint</th>" : ""}<th>treatment</th><th>include</th><th class="num">images</th><th>note</th></tr></thead>
    <tbody>${t.rows.map((r, i) => `<tr data-i="${i}">
      <td>${esc(r[ix.n_folder])}</td><td>${esc(r[ix.treatment_folder])}</td>
      <td><input type="text" data-c="n" value="${esc(r[ix.n])}"></td>
      <td><input type="text" data-c="date" value="${esc(r[ix.date] ?? "")}" size="11"></td>
      ${"timepoint" in ix ? `<td><input type="text" data-c="timepoint" value="${esc(r[ix.timepoint] ?? "")}" size="6"></td>` : ""}
      <td><input type="text" data-c="treatment" value="${esc(r[ix.treatment])}"></td>
      <td><input type="checkbox" data-c="include" ${/^(yes|y|true|1)$/i.test(r[ix.include]) ? "checked" : ""}></td>
      <td class="num">${esc(r[ix.images])}</td><td class="small muted">${esc(r[ix.note])}</td></tr>`).join("")}</tbody></table></div>`;
  S.layout = t;
  b.querySelectorAll("#layoutTbl input").forEach(el => el.addEventListener("input", () => { $("layoutState").textContent = "unsaved changes"; el.closest("tr").classList.add("changed"); }));
};
async function saveLayout() {
  const t = S.layout, ix = Object.fromEntries(t.columns.map((c, i) => [c, i]));
  const rows = [...document.querySelectorAll("#layoutTbl tbody tr")].map(tr => {
    const r = t.rows[+tr.dataset.i], v = c => tr.querySelector(`[data-c=${c}]`);
    return { n_folder: r[ix.n_folder], treatment_folder: r[ix.treatment_folder], n: v("n").value, date: v("date").value,
             ...(v("timepoint") ? { timepoint: v("timepoint").value } : {}),
             treatment: v("treatment").value, include: v("include").checked ? "yes" : "no" };
  });
  await post("/api/layout", { experiment: S.exp, rows });
  $("layoutState").textContent = "saved — run the analysis again to use it";
  S.ov = null;
  document.querySelectorAll("#layoutTbl tr.changed").forEach(tr => tr.classList.remove("changed"));
}
async function useNewLayout() {
  confirmBox("Replace <code>layout.csv</code> with the fresh scan? Any names you edited in the current layout are lost.", async () => {
    await post("/api/layout_use_new", { experiment: S.exp }); S.project = null; render();
  });
}

// ---------------------------------------------------------------- Settings (every form-item kind)

function tomlTable(v) {  // kind "table": { name = ["a", "b"] } shown as TOML inline-table text; "auto" stays a word
  if (typeof v === "string") return v;
  const keys = Object.keys(v || {});
  return keys.length ? `{ ${keys.map(k => `${/^[A-Za-z0-9_-]+$/.test(k) ? k : JSON.stringify(k)} = [${[].concat(v[k]).map(x => JSON.stringify(x)).join(", ")}]`).join(", ")} }` : "{}";
}
function itemInput(sec, f, v) {
  const id = `f_${sec}_${f.key}`;
  const setBtn = s => `<button class="small" type="button" onclick="$('${id}').value=${esc(JSON.stringify(s))};$('${id}').dispatchEvent(new Event('input',{bubbles:true}))">${s === "" ? "clear" : esc(s)}</button>`;
  switch (f.kind) {
    case "bool": return `<label><input type="checkbox" id="${id}" ${v ? "checked" : ""}> on</label>`;
    case "choice": return `<select id="${id}">${f.options.map(o => `<option ${o === v ? "selected" : ""} value="${esc(o)}">${esc(o === "" ? "(empty)" : o)}</option>`).join("")}</select>`;
    case "multi": return f.options.map(o => `<label class="chk"><input type="checkbox" data-multi="${id}" value="${esc(o)}" ${(v || []).includes(o) ? "checked" : ""}> ${esc(o)}</label>`).join("");
    case "names": {  // "all", or some of the options
      const all = v === "all" || v == null;
      return `<label class="chk"><input type="checkbox" id="${id}_all" ${all ? "checked" : ""} onchange="document.querySelectorAll('[data-names=${id}]').forEach(x=>{x.disabled=this.checked})"> <b>all</b></label><br>` +
        f.options.map(o => `<label class="chk"><input type="checkbox" data-names="${id}" value="${esc(o)}" ${all || (v || []).includes(o) ? "checked" : ""} ${all ? "disabled" : ""}> ${esc(o)}</label>`).join("");
    }
    case "list": return `<input type="text" id="${id}" value="${esc((v || []).join(", "))}" style="width:100%" placeholder="comma-separated">`;
    case "measure": return `<input type="text" id="${id}" value="${esc(Array.isArray(v) ? v.join(", ") : v)}" size="14"> ${setBtn("all")}`;
    case "num_or": return `<input type="text" id="${id}" value="${esc(v)}" size="12" placeholder="${esc(f.special.filter(Boolean).join(" / ") || "not set")}">
      ${f.special.map(setBtn).join(" ")}`;
    case "int": case "float": return `<input type="text" id="${id}" value="${esc(v)}" size="10" inputmode="decimal">${f.min != null || f.max != null ? ` <span class="small muted">${f.min != null ? "min " + f.min : ""}${f.min != null && f.max != null ? ", " : ""}${f.max != null ? "max " + f.max : ""}</span>` : ""}`;
    case "table": return `<textarea id="${id}" rows="2" spellcheck="false" placeholder='{ mirna = ["miR-*"] }'>${esc(tomlTable(v))}</textarea>
      <div class="row" style="margin-top:4px">${(f.special || []).map(setBtn).join("")}${setBtn("{}")}<span class="small muted">TOML table: class = ["name", "pattern*"]</span></div>`;
    default: return `<input type="text" id="${id}" value="${esc(v)}" style="width:100%">`;
  }
}
function itemValue(sec, f) {
  const id = `f_${sec}_${f.key}`;
  if (f.kind === "bool") return $(id).checked;
  if (f.kind === "multi") return [...document.querySelectorAll(`[data-multi="${id}"]:checked`)].map(x => x.value);
  if (f.kind === "names") return $(id + "_all").checked ? "all" : [...document.querySelectorAll(`[data-names="${id}"]:checked`)].map(x => x.value);
  if (f.kind === "list") return $(id).value.split(",").map(s => s.trim()).filter(Boolean);
  return $(id).value;
}
RENDER.settings = async function (b) {
  const p = S.project;
  b.innerHTML = `<p class="muted" style="margin-top:0">Every setting of this experiment. Changes are saved to <code>settings.toml</code>; after saving, kerncount tells you which steps to re-run.
      Settings marked <span class="tag fail">re-detects</span> make the next run re-detect every nucleus (slow); the others reuse the detections.</p>
    <div class="bar"><button class="primary" onclick="saveSettings()">Save settings</button>
      <span class="muted small" id="setState"></span><span class="spacer"></span>
      <select onchange="document.getElementById('sec_'+this.value).scrollIntoView()"><option>Jump to section…</option>${p.schema.map(s => `<option value="${s.section}">${esc(s.title)}</option>`).join("")}</select>
      <button onclick="resetSettings()">Reset to defaults…</button></div>
    <div id="rerun"></div>
    ${p.schema.map(sec => `<div class="form-sec" id="sec_${sec.section}"><h3>${esc(sec.title)} <span class="key small muted">[${sec.section}]</span></h3>${sec.items.map(f => {
      const v = (p.settings[sec.section] || {})[f.key];
      return `<div class="form-item" data-sec="${sec.section}" data-key="${f.key}"><div><div class="lbl">${esc(f.label || f.key)} ${f.reseg ? `<span class="tag fail">re-detects</span>` : ""}</div>
        <div class="key">${sec.section}.${f.key} · default ${esc(JSON.stringify(f.default))}</div></div>
        <div>${itemInput(sec.section, f, v)}</div>${f.help ? `<div class="help">${esc(f.help)}</div>` : ""}</div>`; }).join("")}</div>`).join("")}
    <details class="card"><summary>Advanced: edit settings.toml as text</summary>
      <p class="muted small">For settings not shown above. The file must stay valid TOML.</p>
      <textarea id="rawToml" rows="24" spellcheck="false">${esc(p.settings_text)}</textarea>
      <div class="row" style="margin-top:8px"><button onclick="saveRawSettings()">Save text</button></div></details>`;
  const mark = el => { el.classList.add("changed"); $("setState").textContent = "unsaved changes"; };
  b.querySelectorAll(".form-item").forEach(el => { el.addEventListener("input", () => mark(el)); el.addEventListener("change", () => mark(el)); el.addEventListener("click", ev => { if (ev.target.tagName === "BUTTON") mark(el); }); });
  if (S.params.sec && $("sec_" + S.params.sec)) $("sec_" + S.params.sec).scrollIntoView();
  if (S.rerun) showRerun(S.rerun);
};
function showRerun(r) {
  const box = $("rerun"); if (!box) return;
  if (!r.changed.length) { box.innerHTML = `<div class="notice info">Saved. Nothing changed that needs a re-run.</div>`; return; }
  box.innerHTML = `<div class="notice ${r.steps.some(s => s.slow) ? "" : "info"}"><b>Saved ${r.changed.length} change(s).</b> To use them:
    ${r.steps.map(s => `<div class="row" style="margin-top:6px"><button class="small primary" ${S.running ? "disabled" : ""} onclick='startJob(${js(s.action)}, ${js(s.options)})'>Start</button>
      <span>${esc(s.label)} <span class="muted small">(because of: ${esc(s.why)})</span></span></div>`).join("")}</div>`;
}
async function saveSettings() {
  const values = {}, p = S.project;
  for (const sec of p.schema) for (const f of sec.items) {
    const el = document.querySelector(`.form-item[data-sec="${sec.section}"][data-key="${f.key}"]`);
    if (el && el.classList.contains("changed")) (values[sec.section] ||= {})[f.key] = itemValue(sec.section, f);
  }
  if (!Object.keys(values).length) return toast("Nothing changed");
  const r = await post("/api/settings", { experiment: S.exp, values });
  S.project = null; S.ov = null; S.rerun = r;
  await render();
  S.rerun = null;
  window.scrollTo(0, 0);
}
async function saveRawSettings() {
  const r = await post("/api/settings", { experiment: S.exp, text: $("rawToml").value });
  S.project = null; S.ov = null; S.rerun = r; await render(); S.rerun = null; window.scrollTo(0, 0);
}
function resetSettings() {
  confirmBox("Reset <b>all</b> settings of this experiment to the defaults? (The raw images folder is kept.)", async () => {
    await post("/api/settings_reset", { experiment: S.exp }); S.project = null; S.ov = null; render();
  }, "Reset settings");
}

// ---------------------------------------------------------------- counts / generations model (shared by Review + Corrections)

async function loadCounts(force) {
  if (!S.counts || force) {
    S.counts = await get("/api/counts", { exp: S.exp });
    const c = S.counts.columns;
    S.counts.ix = Object.fromEntries(c.map((x, i) => [x, i]));
    S.counts.manual = c.filter(x => /^manual_\d+$/i.test(x)).sort((a, b) => +a.slice(7) - +b.slice(7));
    S.counts.rowByKey = {};
    S.counts.rows.forEach(r => { S.counts.rowByKey[rowKey(r)] = r; });
  }
  const M = S.counts.manual;
  // the chain uses manual_1, manual_2, ... up to the first completely empty column
  let firstEmpty = M.findIndex(col => !S.counts.rows.some(r => (r[S.counts.ix[col]] || "").trim()));
  if (firstEmpty < 0) firstEmpty = M.length;
  S.counts.active = M.slice(0, firstEmpty);
  S.counts.newCol = `manual_${firstEmpty + 1}`;
  S.counts.ignored = M.slice(firstEmpty + 1).filter(col => S.counts.rows.some(r => (r[S.counts.ix[col]] || "").trim()));
  if (!S.round || (S.round !== S.counts.newCol && !S.counts.active.includes(S.round))) S.round = S.counts.newCol;
  return S.counts;
}
const rowKey = r => `${r[S.counts.ix.n]}\u0001${r[S.counts.ix.treatment]}\u0001${parseInt(r[S.counts.ix.image], 10)}`;
const isExcl = v => /^(x|exclude|excl|excluded|-|na|n\/a)$/i.test(String(v).trim());
function cell(r, col) {
  const k = rowKey(r);
  if (col === S.round && k in S.pending) return S.pending[k];
  return S.counts.ix[col] !== undefined ? (r[S.counts.ix[col]] || "") : "";
}
// value of a row in generation g(k): auto, then manual_1..manual_{k-1}
function valueAt(r, k) {
  let v = r[S.counts.ix.auto_count], ex = false;
  for (let i = 1; i < k; i++) {
    const m = String(cell(r, `manual_${i}`)).trim();
    if (!m) continue;
    if (isExcl(m)) ex = true; else { v = m; ex = false; }
  }
  return ex ? "x" : v;
}
function validManual(v) {
  const s = String(v).trim();
  return s === "" || isExcl(s) || /^\d+([.,]0*)?$/.test(s);
}
function roundOptions() {
  const c = S.counts;
  return c.active.map(col => `<option value="${col}" ${col === S.round ? "selected" : ""}>Revise round ${col.slice(7)} (${col} → g${+col.slice(7) + 1})</option>`).join("") +
    `<option value="${c.newCol}" ${c.newCol === S.round ? "selected" : ""}>New round ${c.newCol.slice(7)} (${c.newCol} → builds g${+c.newCol.slice(7) + 1})</option>`;
}
async function saveCorrections(thenBuild) {
  const changes = Object.entries(S.pending).map(([k, value]) => {
    const [n, treatment, image] = k.split("\u0001");
    return { n, treatment, image, col: S.round, value };
  });
  const bad = changes.filter(c => !validManual(c.value));
  if (bad.length) return toast(`${bad.length} value(s) are not a whole number or x — fix the red cells`, true);
  if (changes.length) {
    const r = await post("/api/counts", { experiment: S.exp, changes, mtime: S.counts.mtime });
    S.pending = {};
    await loadCounts(true);
    S.counts.mtime = r.mtime;
  }
  if (thenBuild) buildGeneration();
  else rerender();
}
function rerender() { if (!$("modal").classList.contains("hidden") && S.lb != null) drawLightbox(); render(); }
function buildGeneration() {
  if (Object.keys(S.pending).length) return toast("Save your corrections first", true);
  const k = +S.round.slice(7), target = `g${k + 1}`;
  const hasValues = S.counts.rows.some(r => String(cell(r, S.round)).trim());
  if (!hasValues) return toast(`${S.round} is empty — type at least one correction first`, true);
  const existing = expInfo().generations.filter(g => genNum(g) >= k + 1);
  if (existing.length) {
    confirmBox(`This rebuilds <b>${existing.join(", ")}</b> from the corrections in <code>${S.round}</code>${existing.length > 1 ? " and the rounds after it" : ""}.
      The old folder${existing.length > 1 ? "s are" : " is"} overwritten. Detection results are reused, so this takes about a minute.`,
      () => startJob("regen", { rebuild: existing.join(" ") }), `Rebuild ${existing.join(", ")}`);
  } else { startJob("regen", {}); toast(`building ${target}`); }
}
function changeRound(sel) {
  const apply = () => { S.pending = {}; S.round = sel.value; render(); };
  if (!Object.keys(S.pending).length) return apply();
  const want = sel.value; sel.value = S.round;
  confirmBox("Discard your unsaved corrections in this round?", () => { sel.value = want; apply(); }, "Discard");
}

// ---------------------------------------------------------------- Corrections

RENDER.corrections = async function (b) {
  if (!expInfo().counted) { b.innerHTML = `<div class="notice">Run the analysis first.</div>`; return; }
  await loadCounts();
  if (!S.review) S.review = await get("/api/review", { exp: S.exp });
  const c = S.counts, ix = c.ix, k = +S.round.slice(7), F = S.filters;
  const ns = [...new Set(c.rows.map(r => r[ix.n]))], ts = [...new Set(c.rows.map(r => r[ix.treatment]))];
  const rows = c.rows.filter(r => (!F.n || r[ix.n] === F.n) && (!F.t || r[ix.treatment] === F.t) && (!F.q || r[ix.qc_level] === F.q)
    && (!F.only || String(cell(r, S.round)).trim()));
  const nPending = Object.keys(S.pending).length;
  const nThis = c.rows.filter(r => String(cell(r, S.round)).trim()).length;
  const gensExist = expInfo().generations;
  b.innerHTML = `
    <p class="muted" style="margin-top:0">Generations: <b>g1</b> = automatic counts. Each round of corrections makes the next generation:
      g1 + round 1 → <b>g2</b>, g2 + round 2 → <b>g3</b>, … Type a corrected count, or <b>x</b> to leave an image out. Empty = no change.
      Existing: ${gensExist.map(g => `<span class="tag">${g}</span>`).join("") || "none"}</p>
    ${c.ignored.length ? `<div class="notice">${c.ignored.join(", ")} ha${c.ignored.length > 1 ? "ve" : "s"} values but come${c.ignored.length > 1 ? "" : "s"} after an empty round, so ${c.ignored.length > 1 ? "they are" : "it is"} ignored. Move them into ${c.newCol} (in Excel) to use them.</div>` : ""}
    <div class="bar">
      <select id="roundSel" onchange="changeRound(this)">${roundOptions()}</select>
      <button class="primary" ${nPending ? "" : "disabled"} onclick="saveCorrections(false)">Save ${nPending || ""} change${nPending === 1 ? "" : "s"}</button>
      <button ${S.running ? "disabled" : ""} onclick="${nPending ? "saveCorrections(true)" : "buildGeneration()"}">${nPending ? "Save &amp; " : ""}build g${k + 1}</button>
      <span class="muted small">${nThis} image(s) corrected in ${S.round}</span>
      <span class="spacer"></span>
      <select onchange="S.filters.n=this.value;render()"><option value="">all n</option>${ns.map(n => `<option ${n === F.n ? "selected" : ""}>${esc(n)}</option>`).join("")}</select>
      <select onchange="S.filters.t=this.value;render()"><option value="">all treatments</option>${ts.map(t => `<option ${t === F.t ? "selected" : ""}>${esc(t)}</option>`).join("")}</select>
      <select onchange="S.filters.q=this.value;render()"><option value="">all QC</option>${["low", "medium", "high"].map(q => `<option ${q === F.q ? "selected" : ""}>${q}</option>`).join("")}</select>
      <label class="small"><input type="checkbox" ${F.only ? "checked" : ""} onchange="S.filters.only=this.checked;render()"> only corrected</label>
    </div>
    <div class="split">
      <div class="tablewrap"><table class="grid" id="corrTbl"><thead><tr><th>n</th><th>treatment</th><th>field</th><th class="num">auto</th><th>QC</th>
        <th class="num" title="value going into this round">g${k}</th><th style="width:110px">${S.round}</th><th class="num" title="result of this round">→ g${k + 1}</th></tr></thead>
      <tbody>${rows.map(r => { const key = rowKey(r), v = cell(r, S.round), after = valueAt(r, k + 1), src = reviewItem(key);
        return `<tr data-k="${esc(key)}" class="${key in S.pending ? "changed" : ""} ${S.sel === key ? "sel" : ""}">
          <td class="small">${nChipShort(r[ix.n])}</td><td>${esc(r[ix.treatment])}</td><td class="num" title="${esc(srcTitle(src))}">${esc(r[ix.image])}</td>
          <td class="num">${esc(r[ix.auto_count])}</td><td><span class="tag ${esc(r[ix.qc_level])}">${esc(r[ix.qc_level])}</span></td>
          <td class="num ${valueAt(r, k) === "x" ? "muted" : ""}">${esc(valueAt(r, k) === "x" ? "excluded" : valueAt(r, k))}</td>
          <td><input type="text" class="${validManual(v) ? "" : "bad"}" value="${esc(v)}" data-k="${esc(key)}"></td>
          <td class="num ${after === "x" ? "muted" : ""}" data-after>${esc(after === "x" ? "excluded" : after)}</td></tr>`; }).join("")}</tbody></table></div>
      <div class="preview" id="corrPreview"><div class="muted">Click a row to see its review image.</div></div>
    </div>`;
  const tbl = $("corrTbl");
  tbl.addEventListener("click", ev => { const tr = ev.target.closest("tr[data-k]"); if (tr) selectRow(tr.dataset.k); });
  tbl.addEventListener("input", ev => {
    const inp = ev.target; if (!inp.dataset.k) return;
    const r = c.rowByKey[inp.dataset.k], orig = r[c.ix[S.round]] ?? "";
    if (inp.value.trim() === String(orig).trim()) delete S.pending[inp.dataset.k]; else S.pending[inp.dataset.k] = inp.value.trim();
    inp.classList.toggle("bad", !validManual(inp.value));
    const tr = inp.closest("tr"); tr.classList.toggle("changed", inp.dataset.k in S.pending);
    const a = valueAt(r, k + 1); tr.querySelector("[data-after]").textContent = a === "x" ? "excluded" : a;
    const n = Object.keys(S.pending).length, btn = b.querySelector(".bar .primary");
    btn.disabled = !n; btn.textContent = `Save ${n || ""} change${n === 1 ? "" : "s"}`;
  });
  tbl.addEventListener("keydown", ev => {
    if (!ev.target.dataset.k || !["Enter", "ArrowDown", "ArrowUp"].includes(ev.key)) return;
    ev.preventDefault();
    const trs = [...tbl.querySelectorAll("tbody tr")], i = trs.indexOf(ev.target.closest("tr"));
    const nx = trs[i + (ev.key === "ArrowUp" ? -1 : 1)];
    if (nx) { const inp = nx.querySelector("input"); inp.focus(); inp.select(); selectRow(nx.dataset.k); }
  });
  tbl.addEventListener("focusin", ev => { if (ev.target.dataset.k) selectRow(ev.target.dataset.k); });
  if (S.sel) selectRow(S.sel, true);
};
function reviewItem(key) {  // the review entry (sources, marker values) of an n / treatment / field
  if (!S.review) return null;
  const [n, t, img] = key.split("\u0001");
  const N = S.review.find(x => x.n === n), T = N && N.treatments.find(x => x.treatment === t);
  return T ? T.images.find(i => String(i.image) === String(img)) : null;
}
function srcTitle(i) { return i && i.sources ? Object.entries(i.sources).map(([k, v]) => (k === "all" ? "" : k + ": ") + v).join("\n") : ""; }
function srcShort(i) {  // the file name(s) of a field, e.g. "cntrl/16.tif" or "Ki67 cntrl/16.tif · DAPI cntrl/17.tif"
  if (!i || !i.sources) return "";
  const per = Object.entries(i.sources).filter(([k]) => k !== "all");
  if (per.length) return per.map(([k, v]) => `${k} ${tail2(v)}`).join(" · ");
  return String(i.sources.all || "").split(";").map(s => tail2(s.trim())).join(" · ");
}
function selectRow(key, scroll) {
  if (!S.counts.rowByKey[key]) return;
  S.sel = key;
  document.querySelectorAll("#corrTbl tr.sel").forEach(tr => tr.classList.remove("sel"));
  const tr = document.querySelector(`#corrTbl tr[data-k="${CSS.escape(key)}"]`);
  if (tr) { tr.classList.add("sel"); if (scroll) tr.scrollIntoView({ block: "center" }); }
  const r = S.counts.rowByKey[key], ix = S.counts.ix, [n, t, img] = key.split("\u0001"), it = reviewItem(key);
  const path = `review/${slugify(n)}/${slugify(t)}/${String(img).padStart(2, "0")}.jpg`;
  $("corrPreview").innerHTML = `<div class="row" style="margin-bottom:6px"><b>${esc(shortN(n))} · ${esc(t)} · field ${esc(img)}</b>
      <span class="tag ${esc(r[ix.qc_level])}">QC ${esc(r[ix.qc_level])}</span><span class="spacer"></span>
      <span class="small muted">auto ${esc(r[ix.auto_count])}</span></div>
    ${it && srcShort(it) ? `<div class="small muted" style="margin-bottom:6px" title="${esc(srcTitle(it))}">file: ${esc(srcShort(it))}</div>` : ""}
    <img src="${fileUrl(path)}" onclick='openLightboxFor(${js(key)})' onerror="this.replaceWith(Object.assign(document.createElement('div'),{className:'muted',textContent:'No review image (run without “skip review images”).'}))">
    <div class="small muted" style="margin-top:6px">Left: original. Right: what was counted. Click to enlarge.</div>`;
}
const slugify = s => String(s).replace(/[^\p{L}\p{N}_.\-]+/gu, "_").replace(/^_+|_+$/g, "");

// ---------------------------------------------------------------- Review images

RENDER.review = async function (b, e) {
  if (!S.review) S.review = await get("/api/review", { exp: S.exp });
  const R = S.review;
  if (!R.length) { b.innerHTML = `<div class="notice">No review images yet — run the analysis (without “skip review images”).</div>`; return; }
  // deep link: #exp/review?n=..&t=..&field=..&overlay=1 (also used by headless screenshots)
  const P = S.params, deep = P.field != null && P.t, want = P.field;
  if (P.n && R.find(x => x.n === P.n)) S.reviewN = P.n;
  if (P.t) S.reviewT = P.t;
  if (P.n || P.t || deep) {  // applied once; afterwards the pills / gallery decide
    for (const k of ["n", "t", "field"]) delete P[k];
    history.replaceState(null, "", hashFor(S.exp, S.sec, P));
  }
  if (!S.reviewN || !R.find(x => x.n === S.reviewN)) S.reviewN = R[0].n;
  const N = R.find(x => x.n === S.reviewN);
  const T = S.reviewT && N.treatments.find(t => t.treatment === S.reviewT);
  const markers = e.intensity ? Object.keys((N.treatments[0].images[0] || {}).marker || {}) : [];
  const tally = imgs => ["high", "medium", "low"].map(q => [q, imgs.filter(i => i.qc === q).length]).filter(x => x[1]).map(([q, c]) => `<span class="tag ${q}">${c} ${q}</span>`).join("");
  let html = `<div class="bar"><div class="pills">${R.map(x => `<button class="pill ${x.n === S.reviewN ? "active" : ""}" onclick='S.reviewN=${js(x.n)};S.reviewT=null;render()'><span class="ndot" style="background:${esc(nColour(x.n))}"></span>${esc(x.n)}</button>`).join("")}</div>
    <span class="spacer"></span><div class="pills">${["all", "low", "medium", "high"].map(q => `<button class="pill ${S.reviewQ === q ? "active" : ""}" onclick="S.reviewQ='${q}';render()">${q === "all" ? "all QC" : "QC " + q}</button>`).join("")}</div></div>`;
  if (!T) {
    const ofQ = imgs => imgs.filter(i => S.reviewQ === "all" || i.qc === S.reviewQ).length;
    html += `<p class="muted">One sheet per treatment; border colour = quality (<span class="tag high">high</span><span class="tag medium">medium</span><span class="tag low">low</span>). Click a treatment to go through its fields.
        ${markers.length ? `Fluorescence experiments: open a field and press <b>Show ${esc(markers[0])} positivity</b> to see which nuclei are counted as positive.` : ""}</p>
      <div class="gallery">${N.treatments.map(t => `<div class="thumb" onclick='S.reviewT=${js(t.treatment)};render()'>
        ${t.overview ? `<img loading="lazy" src="${fileUrl(t.overview)}">` : `<div class="muted" style="padding:30px">no overview sheet</div>`}
        <div class="cap"><b>${esc(t.treatment)}</b><span>${tally(t.images)}</span>${S.reviewQ !== "all" ? `<span class="small muted">${ofQ(t.images)} QC ${S.reviewQ}</span>` : ""}</div></div>`).join("")}</div>`;
  } else {
    const imgs = T.images.filter(i => S.reviewQ === "all" || i.qc === S.reviewQ);
    S.lbList = imgs.map(i => ({ ...i, n: N.n, treatment: T.treatment }));
    html += `<div class="row" style="margin-bottom:10px"><button onclick="S.reviewT=null;render()">← all treatments</button>
        <b style="font-size:15px">${esc(T.treatment)}</b> ${tally(T.images)}<span class="spacer"></span>
        ${T.overview ? `<button onclick='openImage(${js(T.overview)},${js(T.treatment + " overview")})'>Overview sheet</button>` : ""}</div>
      ${imgs.length ? "" : `<div class="muted">No QC-${esc(S.reviewQ)} fields in this treatment.</div>`}
      <div class="gallery small">${imgs.map((i, j) => `<div class="thumb q-${esc(i.qc)}" onclick="openLightbox(${j})">
        ${i.path ? `<img loading="lazy" src="${fileUrl(i.path)}">` : `<div class="muted" style="padding:30px">no image</div>`}
        <div class="cap"><b>field ${i.image}</b><span>${esc(i.count)} <span class="muted">(${esc(i.low)}–${esc(i.high)})</span></span><span class="tag ${esc(i.qc)}">${esc(i.qc)}</span></div>
        ${i.duplicate ? `<div class="src" style="color:var(--red)" title="byte-identical raw file">identical to ${esc(i.duplicate)}</div>` : ""}
        <div class="src" title="${esc(srcTitle(i))}">${esc(srcShort(i))}${markers.map(m => i.marker[m] && i.marker[m].positive !== "" ? ` · ${esc(m)}+ ${pct(i.marker[m].positive, 0)}` : "").join("")}</div></div>`).join("")}</div>`;
  }
  b.innerHTML = html;
  if (deep && S.lbList) {
    const j = S.lbList.findIndex(i => String(i.image) === String(want));
    S.lb = j < 0 ? 0 : j; await loadCounts(); S.lbOverlay = P.overlay === "1"; drawLightbox();
  }
};

// lightbox: one field at a time, arrow keys, a correction box, and (fluorescence) the marker overlay
async function openLightbox(j) { S.lb = j; S.lbOverlay = false; await loadCounts(); drawLightbox(); }
async function openLightboxFor(key) {
  await loadCounts();
  const [n, t, img] = key.split("\u0001");
  if (!S.review) S.review = await get("/api/review", { exp: S.exp });
  const N = S.review.find(x => x.n === n), T = N && N.treatments.find(x => x.treatment === t);
  if (!T) return;
  S.lbList = T.images.map(i => ({ ...i, n, treatment: t }));
  S.lb = S.lbList.findIndex(i => String(i.image) === String(img));
  S.lbOverlay = false;
  drawLightbox();
}
function drawLightbox() {
  const i = S.lbList[S.lb]; if (!i) return closeModal();
  const key = `${i.n}\u0001${i.treatment}\u0001${i.image}`, r = S.counts.rowByKey[key], k = +S.round.slice(7);
  const v = r ? cell(r, S.round) : "";
  const markers = Object.keys(i.marker || {}), mk = markers[0];
  const ovUrl = mk ? "/api/field_overlay?" + qs({ exp: S.exp, n: i.n, treatment: i.treatment, image: i.image, channel: mk, t: TOKEN }) : "";
  const mv = mk ? i.marker[mk] : null;
  openModal(`<div class="modal-head"><span class="title">${nChip(i.n)} ${esc(i.treatment)} · field ${i.image}</span>
      <span class="tag ${esc(i.qc)}">QC ${esc(i.qc)}</span><span>count <b>${esc(i.count)}</b> <span class="muted">(plausible ${esc(i.low)}–${esc(i.high)})</span></span>
      <span class="spacer"></span><span class="small muted"><span class="kbd">←</span> <span class="kbd">→</span> browse · <span class="kbd">Esc</span> close</span>
      <button onclick="lbStep(-1)" ${S.lb ? "" : "disabled"}>‹ Prev</button><button onclick="lbStep(1)" ${S.lb < S.lbList.length - 1 ? "" : "disabled"}>Next ›</button><button onclick="closeModal()">Close</button></div>
    <div class="small" style="margin-bottom:6px"><span class="muted">file${Object.keys(i.sources || {}).length > 1 ? "s" : ""}:</span> <span title="${esc(srcTitle(i))}">${esc(srcShort(i)) || "–"}</span>
      ${i.flags ? ` · <span class="muted">flags:</span> ${esc(i.flags)}` : ""}${i.duplicate ? ` · <span class="tag fail">identical file: ${esc(i.duplicate)}</span>` : ""}${i.check ? ` · <span class="tag medium">${esc(i.check)}</span>` : ""}</div>
    ${mk ? `<div class="row small" style="margin-bottom:6px">
        <button class="small ${S.lbOverlay ? "primary" : ""}" onclick="S.lbOverlay=!S.lbOverlay;drawLightbox()">${S.lbOverlay ? "Show review image" : `Show ${esc(mk)} positivity`}</button>
        <span>${esc(mk)}-positive: <b>${pct(mv.positive)}</b> of ${esc(mv.measured || "–")} measured nuclei (cut-off rule: ${esc(mv.rule || "–")})</span>
        ${S.lbOverlay ? `<span class="legend"><span><i style="border-color:#28dc3c"></i>positive</span><span><i style="border-color:#969696"></i>negative</span><span><i style="border-color:#ff9600"></i>saturated (excluded)</span></span>` : ""}</div>` : ""}
    ${S.lbOverlay && mk ? `<img class="full" src="${ovUrl}" onerror="this.replaceWith(Object.assign(document.createElement('div'),{className:'notice',textContent:'No marker measurements for this field (mis-paired files?).'}))">`
      : i.path ? `<img class="full" src="${fileUrl(i.path)}">` : `<div class="muted">no review image</div>`}
    ${r ? `<div class="row" style="margin-top:10px"><span>In g${k}: <b>${esc(valueAt(r, k))}</b></span>
      <label>Correction for <select id="lbRound">${roundOptions()}</select>
      <input id="lbVal" size="8" value="${esc(v)}" placeholder="count or x"></label>
      <button class="primary" id="lbSave">Save</button><span class="muted small">Enter = save &amp; next. Build the generation in Corrections.</span></div>` : ""}`);
  if (!r) return;
  $("lbRound").onchange = e => { if (Object.keys(S.pending).length) { toast("Save or discard your corrections first", true); e.target.value = S.round; return; } S.round = e.target.value; drawLightbox(); };
  const save = async next => {
    const val = $("lbVal").value.trim();
    if (!validManual(val)) return toast("Type a whole number, x, or leave it empty", true);
    if (val !== String(r[S.counts.ix[S.round]] ?? "").trim()) {
      const keep = { ...S.pending }; S.pending = { [key]: val };
      try { await saveCorrections(false); } finally { S.pending = { ...keep }; delete S.pending[key]; }
      toast(`saved ${S.round} = ${val || "(empty)"}`);
    }
    if (next && S.lb < S.lbList.length - 1) { S.lb++; } drawLightbox();
  };
  $("lbSave").onclick = () => save(false);
  $("lbVal").onkeydown = e => { if (e.key === "Enter") { e.preventDefault(); save(true); } e.stopPropagation(); };
}
function lbStep(d) { S.lb = Math.max(0, Math.min(S.lbList.length - 1, S.lb + d)); drawLightbox(); }

// ---------------------------------------------------------------- All files

const fmtSize = n => n > 1e6 ? (n / 1e6).toFixed(1) + " MB" : n > 1e3 ? Math.round(n / 1e3) + " kB" : n + " B";
function fileTable(files, base) {
  return `<table class="grid filelist"><tbody>${files.map(f => `<tr><td><a onclick='previewFile(${js(f.path)})'>${esc(f.path.slice(base.length).replace(/^\//, ""))}</a></td>
    <td class="num muted small">${f.size != null ? fmtSize(f.size) : ""}</td><td style="width:1%"><a class="small" href="${fileUrl(f.path, 1)}">download</a></td></tr>`).join("")}</tbody></table>`;
}
RENDER.files = async function (b) {
  const files = await get("/api/files", { exp: S.exp, dir: S.dir });
  const parts = S.dir ? S.dir.split("/") : [];
  b.innerHTML = `<div class="row" style="margin-bottom:10px"><b>Path:</b>
      <a href="#" onclick="S.dir='';render();return false">${esc(S.exp)}</a>${parts.map((p, i) => ` / <a href="#" onclick='S.dir=${js(parts.slice(0, i + 1).join("/"))};render();return false'>${esc(p)}</a>`).join("")}
      <span class="spacer"></span><button onclick="openSys(null, S.dir || '.')">Open in file browser</button></div>
    <table class="grid filelist"><tbody>${S.dir ? `<tr><td><a onclick='S.dir=${js(parts.slice(0, -1).join("/"))};render()'>⬑ up</a></td><td></td><td></td></tr>` : ""}
    ${files.map(f => f.dir ? `<tr><td>📁 <a onclick='S.dir=${js(f.path)};render()'>${esc(f.name)}</a></td><td></td><td></td></tr>`
      : `<tr><td>📄 <a onclick='previewFile(${js(f.path)})'>${esc(f.name)}</a></td><td class="num small muted">${fmtSize(f.size)}</td>
         <td style="width:1%"><a class="small" href="${fileUrl(f.path, 1)}">download</a></td></tr>`).join("")}</tbody></table>`;
};

// ---------------------------------------------------------------- viewers

function openModal(html, cls = "") { $("modalBox").className = "modal-box " + cls; $("modalBox").innerHTML = html; $("modal").classList.remove("hidden"); }
function closeModal() { $("modal").classList.add("hidden"); $("modalBox").innerHTML = ""; S.lb = null; }
function openImage(path, title, desc, pdf) {
  openModal(`<div class="modal-head"><span class="title">${esc(title || path)}</span><span class="spacer"></span>
    ${pdf ? `<a class="btn" href="${fileUrl(pdf, 1)}">PDF</a>` : ""}<a class="btn" href="${fileUrl(path, 1)}">Download</a><button onclick="closeModal()">Close</button></div>
    ${desc ? `<div class="small muted" style="margin:-4px 0 8px">${esc(desc)}</div>` : ""}<img class="full" src="${fileUrl(path)}">`);
}
function fmtCell(v) {
  if (v === "" || v == null || v === "nan") return "";
  if (!isNaN(+v) && String(v).includes(".") && String(v).length > 7) { const x = +v; return Math.abs(x) >= 1000 ? x.toFixed(1) : x.toPrecision(4); }
  return v;
}
async function viewTable(path) {
  const t = await get("/api/table", { exp: S.exp, path });
  const numeric = t.columns.map((_, j) => t.rows.slice(0, 50).every(r => r[j] === "" || r[j] === undefined || r[j] === "nan" || !isNaN(+r[j])));
  openModal(`<div class="modal-head"><span class="title">${esc(path)}</span><span class="muted small">${t.total} rows${t.truncated ? " (first " + t.rows.length + " shown)" : ""}</span>
      <span class="spacer"></span><a class="btn" href="${fileUrl(path, 1)}">Download</a><button onclick="closeModal()">Close</button></div>
    ${t.total ? `<div class="tablewrap"><table class="grid"><thead><tr>${t.columns.map((c, j) => `<th class="${numeric[j] ? "num" : ""}">${esc(c)}</th>`).join("")}</tr></thead>
    <tbody>${t.rows.map(r => `<tr>${r.map((v, j) => `<td class="${numeric[j] ? "num" : ""}" title="${esc(v)}">${esc(fmtCell(v))}</td>`).join("")}</tr>`).join("")}</tbody></table></div>`
      : `<div class="muted">This table is empty.</div>`}`);
}
async function viewText(path) {
  const t = await get("/api/text", { exp: S.exp, path });
  openModal(`<div class="modal-head"><span class="title">${esc(path)}</span><span class="spacer"></span>
    <a class="btn" href="${fileUrl(path, 1)}">Download</a><button onclick="closeModal()">Close</button></div><pre class="text">${esc(t.text)}</pre>`);
}
function previewFile(path) {
  const ext = path.split(".").pop().toLowerCase();
  if (["png", "jpg", "jpeg", "gif", "svg"].includes(ext)) return openImage(path, path);
  if (["csv", "tsv"].includes(ext)) return viewTable(path);
  if (["txt", "toml", "json", "md", "log", "ijm"].includes(ext)) return viewText(path);
  if (ext === "pdf") return openModal(`<div class="modal-head"><span class="title">${esc(path)}</span><span class="spacer"></span>
      <a class="btn" href="${fileUrl(path, 1)}">Download</a><button onclick="closeModal()">Close</button></div><iframe src="${fileUrl(path)}" style="width:100%;height:78vh;border:0"></iframe>`);
  location.href = fileUrl(path, 1);
}

document.addEventListener("keydown", e => {
  if ($("modal").classList.contains("hidden")) return;
  if (e.key === "Escape") closeModal();
  if (S.lb != null && ["INPUT", "SELECT", "TEXTAREA"].indexOf(document.activeElement.tagName) < 0) {
    if (e.key === "ArrowLeft") lbStep(-1);
    if (e.key === "ArrowRight") lbStep(1);
  }
});
window.addEventListener("beforeunload", e => { if (Object.keys(S.pending).length) { e.preventDefault(); e.returnValue = ""; } });
window.addEventListener("hashchange", applyHash);

// start after every script file has loaded (results.js, intensity.js, plan.js add their sections to RENDER)
window.addEventListener("DOMContentLoaded", async () => {
  await refreshState();
  const h = parseHash(), ex = S.state.experiments;
  if (!(h.exp && ex.find(e => e.name === h.exp)) && ex.length === 1) history.replaceState(null, "", hashFor(ex[0].name, "home", {}));
  await applyHash();
  pollLog();
});
