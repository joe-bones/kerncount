"use strict";
// kerncount web interface. Plain JavaScript, no build step. Talks to gui.py over /api/* (token in X-Token / ?t=).

const S = { state: null, exp: null, tab: "run", project: null, logNext: 0, running: false,
            counts: null, round: null, pending: {}, sel: null, review: null, reviewN: null, reviewT: null, reviewQ: "all",
            filters: { n: "", t: "", q: "", only: false, text: "" }, gen: null, samp: null, dir: "" };
const $ = id => document.getElementById(id);
const esc = s => String(s ?? "").replace(/[&<>"']/g, c => "&#" + c.charCodeAt(0) + ";");
const qs = o => Object.entries(o).map(([k, v]) => encodeURIComponent(k) + "=" + encodeURIComponent(v)).join("&");
const fileUrl = (path, dl) => "/file?" + qs({ exp: S.exp, path, t: TOKEN, ...(dl ? { dl: 1 } : {}) });
const genNum = g => parseInt(String(g).slice(1), 10);

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
function confirmBox(html, onYes, yesLabel = "Yes, go ahead") {
  openModal(`<div style="max-width:560px"><div style="margin-bottom:14px">${html}</div>
    <div class="row"><button class="primary" id="cfYes">${esc(yesLabel)}</button><button onclick="closeModal()">Cancel</button></div></div>`);
  $("cfYes").onclick = () => { closeModal(); onYes(); };
}

// ---------------------------------------------------------------- shell: experiments, navigation, log

async function refreshState() {
  S.state = await get("/api/state");
  $("wsLabel").textContent = S.state.workspace || "no workspace chosen";
  const j = S.state.job, was = S.running;
  S.running = !!(j && j.running);
  $("stopBtn").disabled = !S.running;
  $("jobLabel").innerHTML = j ? `<span class="dot ${S.running ? "run" : j.returncode === 0 ? "" : "fail"}"></span>` +
    esc(j.label) + (S.running ? " — running…" : j.returncode === 0 ? " — done" : " — stopped / failed") : "Log";
  $("expList").innerHTML = S.state.experiments.map(e => `<button class="exp ${e.name === S.exp ? "active" : ""}" onclick='go(${esc(JSON.stringify(e.name))})'>
      <span class="nm">${esc(e.name)}</span>
      <span class="st">${!e.initialised ? "not set up" : !e.generations.length ? "set up, not run" : e.generations.length + " generation(s)"}</span></button>`).join("")
    || `<div class="muted small" style="padding:0 8px">No experiments yet.</div>`;
  if (was && !S.running) onJobFinished();
  return S.state;
}
function expInfo() { return S.state && S.state.experiments.find(e => e.name === S.exp); }

async function go(exp, tab) {
  S.exp = exp;
  if (tab) S.tab = tab;
  S.project = null; S.counts = null; S.review = null; S.pending = {}; S.round = null; S.gen = null; S.samp = null; S.dir = "";
  await refreshState();
  render();
}
async function setTab(t) {
  if (Object.keys(S.pending).length && S.tab === "corrections" && t !== "corrections") {
    if (!window.confirm("You have unsaved corrections. Leave them unsaved?")) return;
    S.pending = {};
  }
  S.tab = t; render();
}

async function render() {
  const m = $("main");
  history.replaceState(null, "", S.exp ? "#" + encodeURIComponent(S.exp) + "/" + S.tab : "#");
  if (!S.exp) return renderHome(m);
  const e = expInfo();
  if (!e) { S.exp = null; return renderHome(m); }
  if (!S.project && e.initialised) {
    try { S.project = await get("/api/project", { exp: S.exp }); } catch (err) { toast(err.message, true); }
  }
  const tabs = [["run", "Run"], ["layout", "Layout"], ["settings", "Settings"], ["review", "Review images"],
                ["corrections", "Corrections"], ["results", "Results"], ["sampling", "Sampling"], ["files", "All files"]];
  m.innerHTML = `<h1>${esc(S.exp)}</h1>
    <div class="muted small">${e.initialised ? esc(S.project ? S.project.path : "") : "Not set up yet"}</div>
    <div class="tabs">${tabs.map(([k, l]) => `<button class="tab ${S.tab === k ? "active" : ""}" ${!e.initialised && k !== "run" ? "disabled" : ""} onclick="setTab('${k}')">${l}</button>`).join("")}</div>
    <div id="tabBody"></div>`;
  const body = $("tabBody");
  if (!e.initialised && S.tab !== "run") S.tab = "run";
  try {
    await ({ run: renderRun, layout: renderLayout, settings: renderSettings, review: renderReview, corrections: renderCorrections,
             results: renderResults, sampling: renderSampling, files: renderFiles })[S.tab](body, e);
  } catch (err) { body.innerHTML = `<div class="notice">${esc(err.message)}</div>`; }
}

function onJobFinished() {
  S.project = null; S.counts = null; S.review = null;
  if (["corrections"].includes(S.tab) && Object.keys(S.pending).length) return; // don't wipe unsaved edits
  if (!$("modal").classList.contains("hidden")) return;
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
  go(name, "run");
}
async function openSys(what, path) { await post("/api/open", { experiment: S.exp, what, path }); }

// ---------------------------------------------------------------- home: workspace

function renderHome(m) {
  const ws = S.state ? S.state.workspace : "";
  m.innerHTML = `<div class="hero">
    <h1>Welcome to kerncount</h1>
    <p class="muted">Automated nucleus counting, review and statistics. Everything happens in this browser tab; nothing leaves this computer.</p>
    <div class="card"><h3 style="margin-top:0">Workspace folder</h3>
      <p class="muted small">The folder that holds <code>data\\</code> (your raw images, one folder per experiment) and
      <code>analysis\\</code> (all results). Results are written next to your data, never into it.</p>
      <div class="row"><input id="wsIn" style="flex:1;min-width:300px" value="${esc(ws)}" placeholder="e.g. C:\\Users\\me\\Documents\\kerncount-work">
      <button onclick="browseWs()">Browse…</button><button class="primary" onclick="setWs(false)">Use this folder</button></div>
      <div class="small" style="margin-top:6px"><a href="#" onclick="setWs(true);return false">Create a new workspace here</a> (makes an empty <code>data</code> folder)</div>
    </div>
    ${ws ? `<div class="card"><h3 style="margin-top:0">How it works</h3><ol>
      <li>Put each experiment's images in <code>data\\&lt;experiment&gt;\\&lt;one folder per n&gt;\\&lt;one folder per treatment&gt;\\1.tif …</code>
          (or use <b>+ Add</b> on the left to create the folder).</li>
      <li>Pick the experiment on the left → <b>Set up</b>, check the <b>Layout</b> and <b>Settings</b> tabs.</li>
      <li><b>Run</b> the analysis (first run ~2–3 s per image; later runs reuse the results).</li>
      <li>Look through <b>Review images</b>; type corrections in <b>Corrections</b> → a new <i>generation</i> (g2, g3, …) is built.</li>
      <li>Read plots, tables and statistics in <b>Results</b>; download the Excel workbook.</li></ol></div>` : ""}
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

// ---------------------------------------------------------------- Run tab

async function renderRun(b, e) {
  const gens = e.generations, last = gens[gens.length - 1];
  const step = (n, done, title, text, controls) => `<div class="step ${done ? "done" : ""}"><div><span class="num">${done ? "✓" : n}</span><b>${title}</b></div><p>${text}</p>${controls}</div>`;
  const d = S.running ? "disabled" : "";
  b.innerHTML = `
    ${S.running ? `<div class="notice info">A job is running — see the log at the bottom. Buttons that change files are paused until it finishes.</div>` : ""}
    <div class="steps">
      ${step(1, e.initialised, "Set up", "Scan the raw image folders and create the project (layout + settings).",
        `<label class="small"><input type="checkbox" id="o_force"> start over: reset settings and layout</label>
         <button class="${e.initialised ? "" : "primary"}" ${d} onclick="startJob('init',{force:$('o_force').checked})">${e.initialised ? "Re-scan folders" : "Set up"}</button>`)}
      ${step(2, e.initialised, "Check layout & settings", "Make sure every folder maps to the right n and treatment; adjust toggles if needed.",
        `<button ${e.initialised ? "" : "disabled"} onclick="setTab('layout')">Layout</button> <button ${e.initialised ? "" : "disabled"} onclick="setTab('settings')">Settings</button>`)}
      ${step(3, gens.length > 0, "Run the analysis", "Count nuclei, score quality, draw review images and build statistics. Always safe to re-run.",
        `<label><input type="checkbox" id="o_noreview"> skip review images (faster)</label>
         <label><input type="checkbox" id="o_reseg"> re-detect every nucleus (ignore cache)</label>
         <label class="small">overwrite generation(s): <input id="o_rebuild" size="8" placeholder="e.g. g2 or all"></label>
         <button class="primary" ${e.initialised ? "" : "disabled"} ${d} onclick="startJob('run',{no_review:$('o_noreview').checked,resegment:$('o_reseg').checked,rebuild:$('o_rebuild').value})">Run analysis</button>`)}
      ${step(4, false, "Review & correct", "Check the review images; type corrected counts (or x to exclude) to build the next generation.",
        `<button ${e.counted ? "" : "disabled"} onclick="setTab('review')">Review images</button> <button ${e.counted ? "" : "disabled"} onclick="setTab('corrections')">Corrections</button>`)}
      ${step(5, false, "Results", last ? `Latest: <b>${last}</b>. Plots, tables, statistics and the Excel workbook.` : "Appears after the first run.",
        `<button ${last ? "" : "disabled"} onclick="setTab('results')">View results</button>
         <button ${last ? "" : "disabled"} onclick="location.href=fileUrl('generations/${last}/results.xlsx',1)">Download Excel</button>`)}
      ${step(6, e.sampling.length > 0, "How many images do we need?", "Variance components, power and subsampling on a generation.",
        `<button ${last ? "" : "disabled"} onclick="setTab('sampling')">Sampling</button>`)}
    </div>
    <h2>On this computer</h2>
    <div class="row"><button onclick="openSys('data')">Open raw images folder</button>
      <button ${e.initialised ? "" : "disabled"} onclick="openSys('folder')">Open project folder</button>
      <button ${last ? "" : "disabled"} onclick="openSys('results')">Open latest results in Excel</button></div>`;
}

// ---------------------------------------------------------------- Layout tab

async function renderLayout(b) {
  const t = await get("/api/table", { exp: S.exp, path: "layout.csv" });
  const ix = Object.fromEntries(t.columns.map((c, i) => [c, i]));
  const hasNew = S.project && S.project.files["layout.new.csv"];
  b.innerHTML = `
    ${hasNew ? `<div class="notice">A fresh scan of the raw folders was saved as <code>layout.new.csv</code> (the folders changed since set-up).
       <button class="small" onclick="viewTable('layout.new.csv')">Compare</button>
       <button class="small" onclick="useNewLayout()">Use the fresh scan</button></div>` : ""}
    <p class="muted">One row per raw folder. Fix the <b>n</b> or <b>treatment</b> name if it was read wrongly; untick <b>include</b> to leave a folder out.
      Then run the analysis again.</p>
    <div class="bar"><button class="primary" id="saveLayout" onclick="saveLayout()">Save layout</button><span class="muted small" id="layoutState"></span></div>
    <div class="tablewrap"><table class="grid" id="layoutTbl"><thead><tr><th>n folder</th><th>treatment folder</th><th>n</th><th>date</th><th>treatment</th><th>include</th><th>images</th><th>note</th></tr></thead>
    <tbody>${t.rows.map((r, i) => `<tr data-i="${i}">
      <td>${esc(r[ix.n_folder])}</td><td>${esc(r[ix.treatment_folder])}</td>
      <td><input type="text" data-c="n" value="${esc(r[ix.n])}"></td>
      <td><input type="text" data-c="date" value="${esc(r[ix.date] ?? "")}" size="10"></td>
      <td><input type="text" data-c="treatment" value="${esc(r[ix.treatment])}"></td>
      <td><input type="checkbox" data-c="include" ${/^(yes|y|true|1)$/i.test(r[ix.include]) ? "checked" : ""}></td>
      <td class="num">${esc(r[ix.images])}</td><td class="small muted">${esc(r[ix.note])}</td></tr>`).join("")}</tbody></table></div>`;
  S.layout = t;
  b.querySelectorAll("#layoutTbl input").forEach(el => el.addEventListener("input", () => { $("layoutState").textContent = "unsaved changes"; el.closest("tr").classList.add("changed"); }));
}
async function saveLayout() {
  const t = S.layout, ix = Object.fromEntries(t.columns.map((c, i) => [c, i]));
  const rows = [...document.querySelectorAll("#layoutTbl tbody tr")].map(tr => {
    const r = t.rows[+tr.dataset.i], v = c => tr.querySelector(`[data-c=${c}]`);
    return { n_folder: r[ix.n_folder], treatment_folder: r[ix.treatment_folder], n: v("n").value, date: v("date").value,
             treatment: v("treatment").value, include: v("include").checked ? "yes" : "no" };
  });
  await post("/api/layout", { experiment: S.exp, rows });
  $("layoutState").textContent = "saved — run the analysis again to use it";
  document.querySelectorAll("#layoutTbl tr.changed").forEach(tr => tr.classList.remove("changed"));
}
async function useNewLayout() {
  confirmBox("Replace <code>layout.csv</code> with the fresh scan? Any names you edited in the current layout are lost.", async () => {
    await post("/api/layout_use_new", { experiment: S.exp }); S.project = null; render();
  });
}

// ---------------------------------------------------------------- Settings tab

function tomlTable(v) {  // kind "table": { name = ["a", "b"] } shown as TOML inline-table text
  if (typeof v === "string") return v;
  const keys = Object.keys(v || {});
  return keys.length ? `{ ${keys.map(k => `${k} = [${[].concat(v[k]).map(x => JSON.stringify(x)).join(", ")}]`).join(", ")} }` : "{}";
}
function fieldInput(sec, f, v) {
  const id = `f_${sec}_${f.key}`;
  switch (f.kind) {
    case "bool": return `<label><input type="checkbox" id="${id}" ${v ? "checked" : ""}> on</label>`;
    case "choice": return `<select id="${id}">${f.options.map(o => `<option ${o === v ? "selected" : ""} value="${esc(o)}">${esc(o === "" ? "(empty)" : o)}</option>`).join("")}</select>`;
    case "multi": return f.options.map(o => `<label class="chk"><input type="checkbox" data-multi="${id}" value="${esc(o)}" ${(v || []).includes(o) ? "checked" : ""}> ${esc(o)}</label>`).join("");
    case "list": return `<input type="text" id="${id}" value="${esc((v || []).join(", "))}" style="width:100%" placeholder="comma-separated">`;
    case "measure": return `<input type="text" id="${id}" value="${esc(Array.isArray(v) ? v.join(", ") : v)}" size="14">`;
    case "num_or": return `<input type="text" id="${id}" value="${esc(v)}" size="12" placeholder="${esc(f.special.filter(Boolean).join(" / ") || "empty")}">
      ${f.special.filter(Boolean).map(s => `<button class="small" type="button" onclick="$('${id}').value='${s}';$('${id}').dispatchEvent(new Event('input',{bubbles:true}))">${s}</button>`).join("")}`;
    case "int": case "float": return `<input type="text" id="${id}" value="${esc(v)}" size="10" inputmode="decimal">`;
    case "table": return `<input type="text" id="${id}" value="${esc(tomlTable(v))}" style="width:100%" placeholder='{ mirna = ["miR-*"] }'>`;
    default: return `<input type="text" id="${id}" value="${esc(v)}" style="width:100%">`;
  }
}
function fieldValue(sec, f) {
  const id = `f_${sec}_${f.key}`;
  if (f.kind === "bool") return $(id).checked;
  if (f.kind === "multi") return [...document.querySelectorAll(`[data-multi="${id}"]:checked`)].map(x => x.value);
  if (f.kind === "list") return $(id).value.split(",").map(s => s.trim()).filter(Boolean);
  return $(id).value;
}
async function renderSettings(b) {
  const p = S.project;
  b.innerHTML = `<p class="muted">Every toggle for this experiment. Changes are saved to <code>settings.toml</code> and used on the next run.
      Settings marked <span class="tag blue">re-detects</span> make the next run re-detect every nucleus (slow); all others only redo the fast steps.</p>
    <div class="bar"><button class="primary" onclick="saveSettings()">Save settings</button>
      <button onclick="saveSettings(true)">Save &amp; run analysis</button>
      <span class="muted small" id="setState"></span><span class="spacer"></span>
      <button onclick="resetSettings()">Reset to defaults…</button></div>
    ${p.schema.map(sec => `<div class="form-sec"><h3>${esc(sec.title)}</h3>${sec.fields.map(f => {
      const v = (p.settings[sec.section] || {})[f.key];
      return `<div class="field" data-sec="${sec.section}" data-key="${f.key}"><div><div class="lbl">${esc(f.label || f.key)} ${f.reseg ? `<span class="tag blue">re-detects</span>` : ""}</div>
        <div class="key">${sec.section}.${f.key} · default ${esc(JSON.stringify(f.default))}</div></div>
        <div>${fieldInput(sec.section, f, v)}</div>${f.help ? `<div class="help">${esc(f.help)}</div>` : ""}</div>`; }).join("")}</div>`).join("")}
    <details class="card"><summary>Advanced: edit settings.toml as text</summary>
      <p class="muted small">For settings not shown above. The file must stay valid TOML.</p>
      <textarea id="rawToml" rows="24" spellcheck="false">${esc(p.settings_text)}</textarea>
      <div class="row" style="margin-top:8px"><button onclick="saveRawSettings()">Save text</button></div></details>`;
  b.querySelectorAll(".field").forEach(el => el.addEventListener("input", () => { el.classList.add("changed"); $("setState").textContent = "unsaved changes"; }));
  b.querySelectorAll(".field").forEach(el => el.addEventListener("change", () => { el.classList.add("changed"); $("setState").textContent = "unsaved changes"; }));
}
async function saveSettings(thenRun) {
  const values = {}, p = S.project;
  for (const sec of p.schema) for (const f of sec.fields) {
    const el = document.querySelector(`.field[data-sec="${sec.section}"][data-key="${f.key}"]`);
    if (el && el.classList.contains("changed")) (values[sec.section] ||= {})[f.key] = fieldValue(sec.section, f);
  }
  if (!Object.keys(values).length && !thenRun) return toast("Nothing changed");
  if (Object.keys(values).length) await post("/api/settings", { experiment: S.exp, values });
  S.project = null;
  if (thenRun) await startJob("run");
  render();
}
async function saveRawSettings() { await post("/api/settings", { experiment: S.exp, text: $("rawToml").value }); S.project = null; render(); }
function resetSettings() {
  confirmBox("Reset <b>all</b> settings of this experiment to the defaults? (The raw images folder is kept.)", async () => {
    await post("/api/settings_reset", { experiment: S.exp }); S.project = null; render();
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
function rerender() { if (!$("modal").classList.contains("hidden") && S.lb) drawLightbox(); render(); }
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
  } else startJob("regen", {});
}

// ---------------------------------------------------------------- Corrections tab

async function renderCorrections(b) {
  if (!expInfo().counted) { b.innerHTML = `<div class="notice">Run the analysis first.</div>`; return; }
  await loadCounts();
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
      <select id="roundSel" onchange="if(Object.keys(S.pending).length&&!confirm('Discard unsaved corrections?')){this.value=S.round;return}S.pending={};S.round=this.value;render()">${roundOptions()}</select>
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
      <div class="tablewrap"><table class="grid" id="corrTbl"><thead><tr><th>n</th><th>treatment</th><th>img</th><th>auto</th><th>QC</th>
        <th title="value going into this round">g${k}</th><th style="width:110px">${S.round}</th><th title="result of this round">→ g${k + 1}</th></tr></thead>
      <tbody>${rows.map(r => { const key = rowKey(r), v = cell(r, S.round), after = valueAt(r, k + 1);
        return `<tr data-k="${esc(key)}" class="${key in S.pending ? "changed" : ""} ${S.sel === key ? "sel" : ""}">
          <td class="small">${esc(r[ix.n])}</td><td>${esc(r[ix.treatment])}</td><td class="num">${esc(r[ix.image])}</td>
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
}
function selectRow(key, scroll) {
  if (!S.counts.rowByKey[key]) return;
  S.sel = key;
  document.querySelectorAll("#corrTbl tr.sel").forEach(tr => tr.classList.remove("sel"));
  const tr = document.querySelector(`#corrTbl tr[data-k="${CSS.escape(key)}"]`);
  if (tr) { tr.classList.add("sel"); if (scroll) tr.scrollIntoView({ block: "center" }); }
  const r = S.counts.rowByKey[key], ix = S.counts.ix, [n, t, img] = key.split("\u0001");
  const path = `review/${slugify(n)}/${slugify(t)}/${String(img).padStart(2, "0")}.jpg`;
  $("corrPreview").innerHTML = `<div class="row" style="margin-bottom:6px"><b>${esc(n)} · ${esc(t)} · image ${esc(img)}</b>
      <span class="tag ${esc(r[ix.qc_level])}">QC ${esc(r[ix.qc_level])}</span><span class="spacer"></span>
      <span class="small muted">auto ${esc(r[ix.auto_count])}</span></div>
    <img src="${fileUrl(path)}" onclick='openLightboxFor(${esc(JSON.stringify(key))})' onerror="this.replaceWith(Object.assign(document.createElement('div'),{className:'muted',textContent:'No review image (run without “skip review images”).'}))">
    <div class="small muted" style="margin-top:6px">Left: original. Right: what was counted. Click to enlarge.</div>`;
}
const slugify = s => String(s).replace(/[^\p{L}\p{N}_.\-]+/gu, "_").replace(/^_+|_+$/g, "");

// ---------------------------------------------------------------- Review tab

async function renderReview(b) {
  if (!S.review) S.review = await get("/api/review", { exp: S.exp });
  const R = S.review;
  if (!R.length) { b.innerHTML = `<div class="notice">No review images yet — run the analysis (without “skip review images”).</div>`; return; }
  if (!S.reviewN || !R.find(x => x.n === S.reviewN)) S.reviewN = R[0].n;
  const N = R.find(x => x.n === S.reviewN);
  const T = S.reviewT && N.treatments.find(t => t.treatment === S.reviewT);
  const tally = imgs => ["high", "medium", "low"].map(q => [q, imgs.filter(i => i.qc === q).length]).filter(x => x[1]).map(([q, c]) => `<span class="tag ${q}">${c} ${q}</span>`).join("");
  let html = `<div class="bar"><div class="pills">${R.map(x => `<button class="pill ${x.n === S.reviewN ? "active" : ""}" onclick='S.reviewN=${esc(JSON.stringify(x.n))};S.reviewT=null;render()'>${esc(x.n)}</button>`).join("")}</div></div>`;
  if (!T) {
    html += `<p class="muted">One sheet per treatment; border colour = quality (<span class="tag high">high</span><span class="tag medium">medium</span><span class="tag low">low</span>). Click a treatment to go through its images.</p>
      <div class="gallery">${N.treatments.map(t => `<div class="thumb" onclick='S.reviewT=${esc(JSON.stringify(t.treatment))};render()'>
        ${t.overview ? `<img loading="lazy" src="${fileUrl(t.overview)}">` : `<div class="muted" style="padding:30px">no overview sheet</div>`}
        <div class="cap"><b>${esc(t.treatment)}</b><span>${tally(t.images)}</span></div></div>`).join("")}</div>`;
  } else {
    const imgs = T.images.filter(i => S.reviewQ === "all" || i.qc === S.reviewQ);
    S.lbList = imgs.map(i => ({ ...i, n: N.n, treatment: T.treatment }));
    html += `<div class="row" style="margin-bottom:10px"><button onclick="S.reviewT=null;render()">← all treatments</button>
        <b style="font-size:15px">${esc(T.treatment)}</b> ${tally(T.images)}<span class="spacer"></span>
        <div class="pills">${["all", "low", "medium", "high"].map(q => `<button class="pill ${S.reviewQ === q ? "active" : ""}" onclick="S.reviewQ='${q}';render()">${q}</button>`).join("")}</div>
        ${T.overview ? `<button onclick='openImage(${esc(JSON.stringify(T.overview))},${esc(JSON.stringify(T.treatment + " overview"))})'>Overview sheet</button>` : ""}</div>
      <div class="gallery small">${imgs.map((i, j) => `<div class="thumb q-${esc(i.qc)}" onclick="openLightbox(${j})">
        ${i.path ? `<img loading="lazy" src="${fileUrl(i.path)}">` : `<div class="muted" style="padding:30px">no image</div>`}
        <div class="cap"><b>#${i.image}</b><span>${esc(i.count)} <span class="muted">(${esc(i.low)}–${esc(i.high)})</span></span><span class="tag ${esc(i.qc)}">${esc(i.qc)}</span></div></div>`).join("")}</div>`;
  }
  b.innerHTML = html;
}

// lightbox: one image at a time, arrow keys, and a correction box writing into the selected round
async function openLightbox(j) { S.lb = j; await loadCounts(); drawLightbox(); }
async function openLightboxFor(key) {
  await loadCounts();
  const [n, t, img] = key.split("\u0001");
  if (!S.review) S.review = await get("/api/review", { exp: S.exp });
  const N = S.review.find(x => x.n === n), T = N && N.treatments.find(x => x.treatment === t);
  if (!T) return;
  S.lbList = T.images.map(i => ({ ...i, n, treatment: t }));
  S.lb = S.lbList.findIndex(i => String(i.image) === String(img));
  drawLightbox();
}
function drawLightbox() {
  const i = S.lbList[S.lb]; if (!i) return closeModal();
  const key = `${i.n}\u0001${i.treatment}\u0001${i.image}`, r = S.counts.rowByKey[key], k = +S.round.slice(7);
  const v = r ? cell(r, S.round) : "";
  openModal(`<div class="modal-head"><span class="title">${esc(i.n)} · ${esc(i.treatment)} · image ${i.image}</span>
      <span class="tag ${esc(i.qc)}">QC ${esc(i.qc)}</span><span>count <b>${esc(i.count)}</b> <span class="muted">(plausible ${esc(i.low)}–${esc(i.high)})</span></span>
      <span class="spacer"></span><span class="small muted"><span class="kbd">←</span> <span class="kbd">→</span> browse · <span class="kbd">Esc</span> close</span>
      <button onclick="lbStep(-1)" ${S.lb ? "" : "disabled"}>‹ Prev</button><button onclick="lbStep(1)" ${S.lb < S.lbList.length - 1 ? "" : "disabled"}>Next ›</button><button onclick="closeModal()">Close</button></div>
    ${i.flags ? `<div class="small" style="margin-bottom:6px">Flags: ${esc(i.flags)}</div>` : ""}
    ${i.path ? `<img class="full" src="${fileUrl(i.path)}">` : `<div class="muted">no review image</div>`}
    ${r ? `<div class="row" style="margin-top:10px"><span>In g${k}: <b>${esc(valueAt(r, k))}</b></span>
      <label>Correction for <select id="lbRound">${roundOptions()}</select>
      <input id="lbVal" size="8" value="${esc(v)}" placeholder="count or x"></label>
      <button class="primary" id="lbSave">Save</button><span class="muted small">Enter = save &amp; next. Build the generation in the Corrections tab.</span></div>` : ""}`);
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

// ---------------------------------------------------------------- Results / Sampling / Files

function fileTable(files, base) {
  return `<table class="grid filelist"><tbody>${files.map(f => `<tr><td><a onclick='previewFile(${esc(JSON.stringify(f.path))})'>${esc(f.path.slice(base.length).replace(/^\//, ""))}</a></td>
    <td class="num muted small">${f.size != null ? fmtSize(f.size) : ""}</td><td style="width:1%"><a class="small" href="${fileUrl(f.path, 1)}">download</a></td></tr>`).join("")}</tbody></table>`;
}
const fmtSize = n => n > 1e6 ? (n / 1e6).toFixed(1) + " MB" : n > 1e3 ? Math.round(n / 1e3) + " kB" : n + " B";
function plotGrid(files) {
  return `<div class="gallery">${files.map(f => `<div class="thumb" onclick='openImage(${esc(JSON.stringify(f.path))},${esc(JSON.stringify(f.name))})'>
    <img loading="lazy" src="${fileUrl(f.path)}"><div class="cap">${esc(f.name.replace(/\.png$/, ""))}</div></div>`).join("")}</div>`;
}

async function renderResults(b) {
  const gens = S.project.generations;
  if (!gens.length) { b.innerHTML = `<div class="notice">No results yet — run the analysis first.</div>`; return; }
  if (!S.gen || !gens.find(g => g.name === S.gen)) S.gen = gens[gens.length - 1].name;
  const g = gens.find(x => x.name === S.gen), base = `generations/${S.gen}`;
  const files = await get("/api/files", { exp: S.exp, dir: base, recursive: "1" });
  const F = files.filter(f => !f.dir), rel = f => f.path.slice(base.length + 1);
  const top = F.filter(f => !rel(f).includes("/")), alln = F.filter(f => rel(f).startsWith("all_n/"));
  const nDirs = [...new Set(F.map(rel).filter(p => p.includes("/") && !p.startsWith("all_n/")).map(p => p.split("/")[0]))];
  const png = l => l.filter(f => f.name.endsWith(".png")), csvs = l => l.filter(f => !f.name.endsWith(".png"));
  b.innerHTML = `<div class="bar"><div class="pills">${gens.map(x => `<button class="pill ${x.name === S.gen ? "active" : ""}" onclick='S.gen=${esc(JSON.stringify(x.name))};render()'>${esc(x.name)}</button>`).join("")}</div>
      <span class="spacer"></span><a class="btn" href="${fileUrl(base + "/results.xlsx", 1)}">Download results.xlsx</a>
      <button onclick="openSys(null,'${base}/results.xlsx')">Open in Excel</button>
      <button onclick="viewTable('${base}/counts_used.csv')">Counts used</button><button onclick="viewText('${base}/README.txt')">How to read this</button></div>
    <div class="card small">${g.source === "auto" ? "Fully automatic counts." : `Built from the previous generation + <b>${esc(g.changes)}</b> manual value(s) in <code>${esc(g.source)}</code>.`}
      ${g.built ? ` Built ${esc(g.built.replace("T", " "))}.` : ""}</div>
    <h2>All n combined</h2>${plotGrid(png(alln))}
    <details><summary>Tables (all n)</summary>${fileTable(csvs(alln), base + "/all_n")}</details>
    <h2>Per n</h2>
    ${nDirs.map(n => { const L = F.filter(f => rel(f).startsWith(n + "/")), direct = L.filter(f => rel(f).split("/").length === 2), deep = L.filter(f => rel(f).split("/").length > 2);
      return `<details class="card"><summary>${esc(n)}</summary>${plotGrid(png(direct))}
        <h3>Tables</h3>${fileTable(csvs(direct), `${base}/${n}`)}
        <details><summary class="small">Per treatment: field values and plots (${deep.length} files)</summary>${plotGrid(png(deep))}${fileTable(csvs(deep), `${base}/${n}`)}</details></details>`; }).join("")}
    <details><summary>Other files</summary>${fileTable(top, base)}</details>`;
}

async function renderSampling(b) {
  const e = expInfo(), gens = e.generations, done = S.project.sampling;
  if (!S.samp || !done.includes(S.samp)) S.samp = done[done.length - 1] || null;
  let html = `<p class="muted" style="margin-top:0">How many images per condition are needed? Splits the variation into field / well / n,
      estimates power for other designs, and subsamples your data. Uses the plate design from <b>Settings</b>.</p>
    <div class="card row"><label>Generation <select id="sGen">${gens.map(g => `<option ${g === gens[gens.length - 1] ? "selected" : ""}>${esc(g)}</option>`).join("")}</select></label>
      <label>Subsamples per design <input id="sIt" value="1000" size="6"></label>
      <button class="primary" ${S.running || !gens.length ? "disabled" : ""} onclick="startJob('sampling',{generation:$('sGen').value,iterations:$('sIt').value})">Run sampling analysis</button></div>`;
  if (!done.length) { b.innerHTML = html + `<div class="muted">No sampling results yet.</div>`; return; }
  const base = `sampling/${S.samp}`;
  const files = (await get("/api/files", { exp: S.exp, dir: base, recursive: "1" })).filter(f => !f.dir);
  const rep = files.find(f => f.name === "REPORT.txt");
  const report = rep ? (await get("/api/text", { exp: S.exp, path: rep.path })).text : "";
  html += `<div class="bar"><span class="muted">Results for</span><div class="pills">${done.map(g => `<button class="pill ${g === S.samp ? "active" : ""}" onclick='S.samp=${esc(JSON.stringify(g))};render()'>${esc(g)}</button>`).join("")}</div></div>
    ${plotGrid(files.filter(f => f.name.endsWith(".png")))}
    ${report ? `<h2>Report</h2><pre class="text">${esc(report)}</pre>` : ""}
    <h2>Tables and files</h2>${fileTable(files.filter(f => !f.name.endsWith(".png")), base)}`;
  b.innerHTML = html;
}

async function renderFiles(b) {
  const files = await get("/api/files", { exp: S.exp, dir: S.dir });
  const parts = S.dir ? S.dir.split("/") : [];
  b.innerHTML = `<div class="row" style="margin-bottom:10px"><b>Path:</b>
      <a href="#" onclick="S.dir='';render();return false">${esc(S.exp)}</a>${parts.map((p, i) => ` / <a href="#" onclick='S.dir=${esc(JSON.stringify(parts.slice(0, i + 1).join("/")))};render();return false'>${esc(p)}</a>`).join("")}
      <span class="spacer"></span><button onclick="openSys(null, S.dir || '.')">Open in file browser</button></div>
    <table class="grid filelist"><tbody>${S.dir ? `<tr><td><a onclick='S.dir=${esc(JSON.stringify(parts.slice(0, -1).join("/")))};render()'>⬑ up</a></td><td></td><td></td></tr>` : ""}
    ${files.map(f => f.dir ? `<tr><td>📁 <a onclick='S.dir=${esc(JSON.stringify(f.path))};render()'>${esc(f.name)}</a></td><td></td><td></td></tr>`
      : `<tr><td>📄 <a onclick='previewFile(${esc(JSON.stringify(f.path))})'>${esc(f.name)}</a></td><td class="num small muted">${fmtSize(f.size)}</td>
         <td style="width:1%"><a class="small" href="${fileUrl(f.path, 1)}">download</a></td></tr>`).join("")}</tbody></table>`;
}

// ---------------------------------------------------------------- viewers

function openModal(html) { $("modalBox").innerHTML = html; $("modal").classList.remove("hidden"); }
function closeModal() { $("modal").classList.add("hidden"); $("modalBox").innerHTML = ""; S.lb = null; }
function openImage(path, title) {
  openModal(`<div class="modal-head"><span class="title">${esc(title || path)}</span><span class="spacer"></span>
    <a class="btn" href="${fileUrl(path, 1)}">Download</a><button onclick="closeModal()">Close</button></div><img class="full" src="${fileUrl(path)}">`);
}
async function viewTable(path) {
  const t = await get("/api/table", { exp: S.exp, path });
  const numeric = t.columns.map((_, j) => t.rows.slice(0, 50).every(r => r[j] === "" || r[j] === undefined || !isNaN(+r[j])));
  const fmt = v => (v !== "" && !isNaN(+v) && String(v).includes(".") && String(v).length > 8) ? (+v).toPrecision(5) : v;
  openModal(`<div class="modal-head"><span class="title">${esc(path)}</span><span class="muted small">${t.total} rows${t.truncated ? " (first " + t.rows.length + " shown)" : ""}</span>
      <span class="spacer"></span><a class="btn" href="${fileUrl(path, 1)}">Download</a><button onclick="closeModal()">Close</button></div>
    <div class="tablewrap"><table class="grid"><thead><tr>${t.columns.map(c => `<th>${esc(c)}</th>`).join("")}</tr></thead>
    <tbody>${t.rows.map(r => `<tr>${r.map((v, j) => `<td class="${numeric[j] ? "num" : ""}" title="${esc(v)}">${esc(fmt(v))}</td>`).join("")}</tr>`).join("")}</tbody></table></div>`);
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

(async () => {
  await refreshState();
  const ex = S.state.experiments, [hExp, hTab] = decodeURIComponent(location.hash.slice(1)).split("/");
  if (hExp && ex.find(e => e.name === hExp)) { S.exp = hExp; if (hTab) S.tab = hTab; }
  else if (ex.length === 1) S.exp = ex[0].name;
  render();
  pollLog();
})();
