"use strict";
// Results explorer, statistics table and layer comparison: all driven by generations/gN/manifest.json and
// stats_long.csv (docs/statistics.md, sections 5-6).

const KINDS = [  // view -> manifest kinds; colours in app.css (.k-*)
  ["absolute", "Absolute", "Values in their own units (nuclei per field, intensity, %)."],
  ["relative", "Relative", "Values as % of each control."],
  ["model", "Model", "Estimates of the hierarchical model, variance components."],
  ["contrast", "Comparisons", "Each treatment compared with a control (forest plots, tests)."],
  ["robustness", "Robustness", "Effect in every n, leave-one-n-out stability."],
  ["diagnostic", "Diagnostics", "Model checks (residuals, normality)."],
  ["design", "Design", "Field positions, wells, density."],
  ["data", "Tables", "Prism tables and raw data."],
];
const KIND_LABEL = Object.fromEntries(KINDS.map(([k, l]) => [k, l]));

async function loadManifest(gen) {
  if (!S.man[gen]) S.man[gen] = await get("/api/manifest", { exp: S.exp, gen });
  return S.man[gen];
}
function genList() { return (S.project ? S.project.generations : []).map(g => g.name); }
function pickGen() {
  const gens = genList();
  return gens.includes(S.params.gen) ? S.params.gen : gens[gens.length - 1];
}
function genPills(gen, onclick = "setP({gen:'$'})") {
  const gens = genList();
  return gens.length > 1 ? `<div class="pills">${gens.map(g => `<button class="pill ${g === gen ? "active" : ""}" onclick="${onclick.replace("$", g)}">${g}</button>`).join("")}</div>` : `<span class="tag">${esc(gen)}</span>`;
}

// metric groups: Counts; Readouts; Intensity per channel (with % positive)
function metricChannel(id) {
  let m = /^readout_frac_positive_(?:p90_)?(.+)$/.exec(id);
  if (m) return m[1];
  m = /^intensity_([^_]+)_/.exec(id);
  return m ? m[1] : null;
}
function metricGroups(man) {
  const g = new Map([["Counts", []], ["Readouts", []]]);
  for (const m of man.metrics || []) {
    const ch = metricChannel(m.id);
    const name = m.id === "count" || m.id.startsWith("nuclei_per") ? "Counts" : ch ? `Intensity: ${ch}` : "Readouts";
    if (!g.has(name)) g.set(name, []);
    g.get(name).push(m);
  }
  return [...g.entries()].filter(([, l]) => l.length);
}
function metricInfo(man, id) { return (man.metrics || []).find(m => m.id === id) || { id, label: id, unit: "", kind: "" }; }
const dirOf = p => p.includes("/") ? p.slice(0, p.lastIndexOf("/")) : "";

// a plot from a manifest; `base` = the folder its manifest paths are relative to ("generations/g1/", "intensity/", ...)
function plotCard(f, base = "", extra = "") {
  const path = base + f.path, pdf = f.pdf ? base + f.pdf : "";
  return `<div class="plot k-${esc(f.kind)}">
    <div class="ph"><div class="pt">${esc(f.title)}</div><div class="pd">${esc(f.description || "")}</div></div>
    <img loading="lazy" src="${fileUrl(path)}" onclick='openImage(${js(path)},${js(f.title)},${js(f.description || "")},${js(pdf)})' alt="${esc(f.title)}">
    <div class="pf"><span class="tag" style="background:var(--kcs);color:var(--kc)">${esc(KIND_LABEL[f.kind] || f.kind)}</span>${extra}<span class="spacer"></span>
      ${pdf ? `<a class="btn small" href="${fileUrl(pdf, 1)}">PDF</a>` : ""}<a class="btn small" href="${fileUrl(path, 1)}">PNG</a></div></div>`;
}
function tableButtons(files, base = "") {
  return files.length ? `<div class="tbllist">${files.map(f => `<button class="btn" title="${esc(f.description || "")}" onclick='viewTable(${js(base + f.path)})'>${esc(f.title.replace(/^(All n|n\d+)[^:]*: /, ""))}</button>`).join("")}</div>` : "";
}
const isPng = f => f.path.endsWith(".png"), isCsv = f => f.path.endsWith(".csv");

// ---------------------------------------------------------------- Results explorer

RENDER.results = async function (b, e, tok) {
  const gens = genList();
  if (!gens.length) { b.innerHTML = `<div class="notice">No results yet — run the analysis first.</div>`; return; }
  const gen = pickGen();
  const M = await loadManifest(gen);
  if (tok !== S.renderTok) return;
  const man = M.manifest;
  if (!man.files) { b.innerHTML = `<div class="notice">${esc(gen)} was built by an older kerncount (no manifest.json). Rebuild it (Experiment home → Rebuild generations, overwrite ${esc(gen)}).</div>`; return; }
  if (!S.ov) S.ov = await get("/api/overview", { exp: S.exp }).catch(() => null);
  const P = S.params, F = man.files;
  const folders = man.all_n_folders && man.all_n_folders.length ? man.all_n_folders : [{ cell_line: null, path: "all_n", n_count: man.n_count }];
  const ns = [...new Set(F.filter(f => f.level === "n").map(f => f.n))].sort();
  const level = ["all_n", "n", "treatment"].includes(P.level) ? P.level : "all_n";
  const folder = folders.find(f => f.path === P.folder) || folders[0];
  const n = ns.includes(P.n) ? P.n : ns[0];
  const ts = [...new Set(F.filter(f => f.level === "treatment" && f.n === n).map(f => f.treatment))];
  const t = ts.includes(P.t) ? P.t : ts[0];
  const atLevel = f => level === "all_n" ? f.level === "all_n" && dirOf(f.path) === folder.path
    : level === "n" ? f.level === "n" && f.n === n : f.level === "treatment" && f.n === n && f.treatment === t;
  const here = F.filter(atLevel);
  const metricsHere = new Set(here.map(f => f.metric).filter(Boolean));
  const groups = metricGroups(man).map(([g, l]) => [g, l.filter(m => metricsHere.has(m.id))]).filter(([, l]) => l.length);
  const metric = metricsHere.has(P.metric) ? P.metric : (metricsHere.has("count") ? "count" : [...metricsHere][0]);
  const mf = here.filter(f => f.metric === metric);
  const counts = Object.fromEntries(KINDS.map(([k]) => [k, mf.filter(f => f.kind === k).length]));
  let view = KIND_LABEL[P.view] && counts[P.view] ? P.view : (counts.absolute ? "absolute" : KINDS.map(k => k[0]).find(k => counts[k]));
  const mi = metricInfo(man, metric), est = (man.estimands || {})[metric] || {};
  const status = (man.model_status || {})[(folder.cell_line && folders.length > 1 ? folder.cell_line + "/" : "") + metric] || (man.model_status || {})[metric];
  const one = level === "all_n" && (folder.n_count || 0) < 2;
  const shown = mf.filter(f => f.kind === view);
  // absolute values next to relative ones: for relative / comparison / robustness views, the absolute plots of the same metric
  const companion = ["relative", "contrast", "robustness"].includes(view)
    ? mf.filter(f => isPng(f) && (f.kind === "absolute" || (f.kind === "model" && /_model_absolute\.png$/.test(f.path)))) : [];
  const colours = M.colours, G = `generations/${gen}/`;
  b.innerHTML = `
    <div class="selrow">
      <label class="f">Generation ${genPills(gen)}</label>
      <label class="f">Level <div class="pills">
        <button class="pill ${level === "all_n" ? "active" : ""}" onclick="setP({level:'all_n'})">All n${folders.length > 1 ? " (per cell line)" : ""}</button>
        <button class="pill ${level === "n" ? "active" : ""}" ${ns.length ? "" : "disabled"} onclick="setP({level:'n'})">Per n</button>
        <button class="pill ${level === "treatment" ? "active" : ""}" ${F.some(f => f.level === "treatment") ? "" : "disabled"} onclick="setP({level:'treatment'})">Per treatment</button></div></label>
      ${level === "all_n" && folders.length > 1 ? `<label class="f">Cell line <select onchange="setP({folder:this.value})">${folders.map(f => `<option value="${esc(f.path)}" ${f === folder ? "selected" : ""}>${esc(f.cell_line || "unknown")} (${f.n_count} n)</option>`).join("")}</select></label>` : ""}
      ${level !== "all_n" ? `<label class="f">n <div class="pills">${ns.map(x => `<button class="pill ${x === n ? "active" : ""}" onclick='setP({n:${js(x)}})'><span class="ndot" style="background:${esc(colours[x] || nColour(x))}"></span>${esc(shortN(x))}</button>`).join("")}</div></label>` : ""}
      ${level === "treatment" ? `<label class="f">Treatment <select onchange="setP({t:this.value})">${ts.map(x => `<option ${x === t ? "selected" : ""}>${esc(x)}</option>`).join("")}</select></label>` : ""}
      <label class="f">Metric <select onchange="setP({metric:this.value})" style="max-width:340px">${groups.map(([g, l]) => `<optgroup label="${esc(g)}">${l.map(m => `<option value="${esc(m.id)}" ${m.id === metric ? "selected" : ""}>${esc(m.label)}</option>`).join("")}</optgroup>`).join("")}</select></label>
    </div>
    <div class="kindbar" style="margin-bottom:12px">${KINDS.map(([k, l, d]) => `<button class="kpill k-${k} ${k === view ? "active" : ""}" title="${esc(d)}" ${counts[k] ? "" : "disabled"} onclick="setP({view:'${k}'})"><i></i>${l}<span class="cnt">${counts[k] || ""}</span></button>`).join("")}</div>
    <div class="metricbox"><b>${esc(mi.label)}</b> <span class="muted">· ${esc(mi.unit)}${mi.kind && mi.kind !== mi.unit ? " · " + esc(mi.kind) : ""}</span>
      <div class="small">${esc(mi.description || "")}</div>
      ${est.effect ? `<div class="small"><b>Effect:</b> ${esc(est.effect)} · <b>Means:</b> ${esc(est.mean)}${est.headline === "pp_difference" ? ` · <span class="tag purple">headline: percentage points</span>` : ""}</div>` : ""}
      ${level === "all_n" && status ? `<div class="small"><b>Model:</b> ${esc(status.model)} · stars from: ${esc(status.stars_from === "mixed" ? "mixed model" : status.stars_from === "rm_anova" ? "RM ANOVA" : "none")}
        ${status.status === "fallback_rm_anova" ? `<span class="tag medium">mixed model failed: RM ANOVA used</span>` : ""}</div>` : ""}
      ${level === "n" ? `<div class="small muted">Within one n the fields are technical replicates: tests here describe how well this experiment was measured, not whether the effect is reproducible.</div>` : ""}</div>
    ${one ? `<div class="notice n1"><b>Only one n${folder.cell_line ? ` (${esc(folder.cell_line)})` : ""}: needs ≥ 2 n.</b> The all-n results show estimates from that n; confidence intervals, p-values and stars across experiments need at least two independent experiments.</div>` : ""}
    ${!shown.length ? `<div class="muted">Nothing of this kind for this selection.</div>` : ""}
    ${companion.length && shown.some(isPng) ? `<div class="cols2">
        <div class="k-absolute"><div class="colhead">Absolute values</div><div class="plots one">${companion.map(f => plotCard(f, G)).join("")}</div></div>
        <div class="k-${view}"><div class="colhead">${esc(KIND_LABEL[view])}</div><div class="plots one">${shown.filter(isPng).map(f => plotCard(f, G)).join("")}</div></div></div>`
      : `<div class="plots">${shown.filter(isPng).map(f => plotCard(f, G)).join("")}</div>`}
    ${shown.filter(isCsv).length ? `<h3>Tables</h3>${tableButtons(shown.filter(isCsv), G)}` : ""}
    ${level !== "all_n" ? `<h3>${level === "treatment" ? "This treatment" : "This n"}</h3><div class="tbllist">${tableButtons(here.filter(f => !f.metric && isCsv(f)), G)}
      <button class="btn" onclick='nav("review",${js(level === "treatment" ? { n, t } : { n })})'>Review images${level === "treatment" ? " of " + esc(t) : ""} →</button></div>` : ""}
    <details style="margin-top:14px"><summary class="small">Generation files: workbook, long statistics table, counts used</summary>
      <div class="row" style="margin-top:6px"><a class="btn" href="${fileUrl(`generations/${gen}/results.xlsx`, 1)}">results.xlsx</a>
        <button onclick="nav('stats',{gen:'${gen}'})">Statistics table</button>
        <button onclick="viewTable('generations/${gen}/counts_used.csv')">Counts used</button>
        <button onclick="viewText('generations/${gen}/README.txt')">How to read this</button></div></details>`;
};

// ---------------------------------------------------------------- Statistics table (stats_long.csv)

const SIG_ORDER = { "****": 4, "***": 3, "**": 2, "*": 1 };
function fmtNum(x, sig = 3) {
  if (x === "" || x == null || isNaN(+x)) return "";
  const v = +x, a = Math.abs(v);
  if (a === 0) return "0";
  if (a >= 1000) return v.toFixed(0);
  if (a >= 100) return v.toFixed(1);
  return (+v.toPrecision(sig)).toString();
}
function fmtEst(v, type) {
  if (v === "" || v == null) return "";
  if (type === "ratio" || type === "odds_ratio") return fmtNum(v) + "×";
  if (type === "pp_difference") return (+v > 0 ? "+" : "") + (+v).toFixed(2) + " pp";
  if (type === "difference") return (+v > 0 ? "+" : "") + fmtNum(v);
  return fmtNum(v, 4);
}
function fmtP(p) {
  if (p === "" || p == null || isNaN(+p)) return "";
  const v = +p;
  return v < 0.0001 ? "<0.0001" : v < 0.001 ? v.toFixed(4) : v < 0.01 ? v.toFixed(3) : v.toFixed(2);
}
function sigBadge(s) {
  if (!s) return "";
  if (SIG_ORDER[s]) return `<span class="sig">${esc(s)}</span>`;
  if (s === "ns") return `<span class="sig ns">ns</span>`;
  if (s === "equivalent") return `<span class="sig eq">equivalent</span>`;
  return `<span class="sig inc">${esc(s)}</span>`;
}
const needs2 = r => /needs ≥ 2 n/.test(r.note || "") || /needs ≥ 2 n/.test(r.model || "");
const isModel = r => /n×treatment|mixed/.test(r.model || "") && !/descriptive/.test(r.model || "");
const isSignif = r => !!SIG_ORDER[r.significance] || r.significance === "equivalent";

RENDER.stats = async function (b, e, tok) {
  const gens = genList();
  if (!gens.length) { b.innerHTML = `<div class="notice">No results yet — run the analysis first.</div>`; return; }
  const gen = pickGen();
  const M = await loadManifest(gen);
  if (!S.stats[gen]) {
    const t = await get("/api/stats", { exp: S.exp, gen });
    t.objs = t.rows.map((r, i) => Object.assign(Object.fromEntries(t.columns.map((c, j) => [c, r[j] ?? ""])), { _i: i }));
    S.stats[gen] = t;
  }
  if (tok !== S.renderTok) return;
  const T = S.stats[gen], man = M.manifest, P = S.params;
  if (!T.objs.length) { b.innerHTML = `<div class="notice">${esc(gen)} has no stats_long.csv (built by an older kerncount). Rebuild it.</div>`; return; }
  const uniq = c => [...new Set(T.objs.map(r => r[c]).filter(v => v !== ""))].sort();
  const filt = { level: P.level || "", n: P.n || "", cell_line: P.cl || "", metric: P.metric || "", family: P.fam || "", estimate_type: P.et || "" };
  const q = (P.q || "").toLowerCase(), sig = P.sig || "", src = P.src || "";
  const srcOf = r => isModel(r) ? "model" : /^descriptive/.test(r.model || "") ? "descriptive" : r.model ? "classical" : "";
  let rows = T.objs.filter(r => Object.entries(filt).every(([c, v]) => !v || r[c] === v)
    && (!q || (r.comparison + " " + r.test).toLowerCase().includes(q)) && (!src || srcOf(r) === src)
    && (!sig || (sig === "sig" ? isSignif(r) : sig === "ns" ? r.significance === "ns" : sig === "need2" ? needs2(r) : sig === "none" ? !r.significance : true)));
  const sortCol = P.sort || "", dir = P.dir === "desc" ? -1 : 1;
  if (sortCol) {
    const num = ["estimate", "ci_low", "ci_high", "p", "p_adj", "N"].includes(sortCol);
    rows = [...rows].sort((a, z) => {
      if (sortCol === "significance") return dir * ((SIG_ORDER[z.significance] || 0) - (SIG_ORDER[a.significance] || 0));
      const x = a[sortCol], y = z[sortCol];
      if (num) { const xa = x === "" ? Infinity : +x, ya = y === "" ? Infinity : +y; return dir * (xa - ya); }
      return dir * String(x).localeCompare(String(y));
    });
  }
  S.statsView = rows;
  const limit = +(P.limit || 300);
  const mlabel = id => metricInfo(man, id).label;
  const sel = (key, col, label, opts, lab = x => x) => `<label class="f">${label}<select onchange="setP({${key}:this.value,limit:''})"><option value="">all</option>${opts.map(o => `<option value="${esc(o)}" ${String(P[key] || "") === o ? "selected" : ""}>${esc(lab(o))}</option>`).join("")}</select></label>`;
  const th = (col, label, cls = "") => `<th class="sort ${cls}" onclick="setP({sort:'${col}',dir:'${sortCol === col && dir === 1 ? "desc" : "asc"}'})">${label}${sortCol === col ? (dir === 1 ? " ▲" : " ▼") : ""}</th>`;
  const one = (man.n_count || 0) < 2;
  b.innerHTML = `
    <p class="muted" style="margin-top:0">Every statistical result of ${esc(gen)} in one table (<code>stats_long.csv</code>). Filter, sort, click a row to open the plot it belongs to.
      Estimates are shown with what they estimate (the <i>estimand</i>); significant results are highlighted.</p>
    ${one ? `<div class="notice n1"><b>Only one n:</b> all-n rows carry estimates only; their CIs, p-values and stars say “needs ≥ 2 n”.</div>` : ""}
    <div class="selrow">
      <label class="f">Generation ${genPills(gen)}</label>
      ${sel("level", "level", "Level", uniq("level"), x => x === "all_n" ? "all n" : "within one n")}
      ${sel("n", "n", "n", uniq("n"))}
      ${uniq("cell_line").length > 1 ? sel("cl", "cell_line", "Cell line", uniq("cell_line")) : ""}
      ${sel("metric", "metric", "Metric", uniq("metric"), mlabel)}
      ${sel("fam", "family", "Family", uniq("family"))}
      ${sel("et", "estimate_type", "Estimate type", uniq("estimate_type"))}
      <label class="f">Significance<select onchange="setP({sig:this.value,limit:''})">${[["", "all"], ["sig", "significant (* or equivalent)"], ["ns", "not significant"], ["need2", "needs ≥ 2 n"], ["none", "no test"]].map(([v, l]) => `<option value="${v}" ${sig === v ? "selected" : ""}>${l}</option>`).join("")}</select></label>
      <label class="f">Test<select onchange="setP({src:this.value,limit:''})">${[["", "all tests"], ["model", "mixed model (headline stars)"], ["classical", "classical (ANOVA, t-tests)"], ["descriptive", "descriptive (no test)"]].map(([v, l]) => `<option value="${v}" ${src === v ? "selected" : ""}>${l}</option>`).join("")}</select></label>
      <label class="f">Comparison contains<input id="stQ" value="${esc(P.q || "")}" size="14" onchange="setP({q:this.value,limit:''})" placeholder="e.g. miR-132"></label>
      <button onclick="setP({level:'',n:'',cl:'',metric:'',fam:'',et:'',sig:'',q:'',src:'',sort:'',dir:'',limit:''})">Clear filters</button>
    </div>
    <div class="bar"><b>${rows.length}</b> <span class="muted">of ${T.objs.length} results</span>
      <span class="small muted">· ${rows.filter(isSignif).length} significant · ${rows.filter(needs2).length} need ≥ 2 n</span>
      <span class="spacer"></span><button onclick="exportStats('${gen}')">Export this view (CSV)</button></div>
    <div class="tablewrap" style="max-height:none"><table class="grid" id="statsTbl"><thead><tr>
      ${th("level", "Level · n")}${filt.metric ? "" : th("metric", "Metric")}${th("family", "Family · test")}${th("comparison", "Comparison")}
      ${th("estimate", "Estimate", "num")}<th class="num">95% CI</th>${th("p", "p", "num")}${th("p_adj", "p adj.", "num")}${th("significance", "Sig.")}${th("N", "N", "num")}<th>Estimand</th><th>Note</th></tr></thead>
      <tbody>${rows.slice(0, limit).map(r => `<tr class="click ${isSignif(r) ? "signif" : ""}" data-i="${r._i}" title="${esc(r.test)}${r.model ? "\nmodel: " + esc(r.model) : ""}">
        <td class="small">${r.level === "all_n" ? `all n${r.cell_line && uniq("cell_line").length > 1 ? " · " + esc(r.cell_line) : ""}` : nChipShort(r.n)}</td>
        ${filt.metric ? "" : `<td class="small">${esc(mlabel(r.metric))}</td>`}
        <td class="small" style="min-width:200px"><b style="white-space:nowrap">${esc(r.family)}</b>${isModel(r) ? ` <span class="tag purple">mixed model</span>` : ""}<div class="muted">${esc(r.test)}</div></td><td>${esc(r.comparison)}</td>
        <td class="num est">${esc(fmtEst(r.estimate, r.estimate_type))}</td>
        <td class="num ci">${r.ci_low !== "" ? `[${esc(fmtEst(r.ci_low, r.estimate_type))}, ${esc(fmtEst(r.ci_high, r.estimate_type))}]` : ""}</td>
        <td class="num">${fmtP(r.p)}</td><td class="num">${fmtP(r.p_adj)}</td><td>${sigBadge(r.significance)}</td><td class="num">${r.N !== "" && !isNaN(+r.N) ? +(+r.N).toFixed(1) : esc(r.N)}</td>
        <td class="small muted">${esc(r.estimand)}</td>
        <td class="small">${needs2(r) ? `<span class="need2">needs ≥ 2 n</span> ` : ""}<span class="muted">${esc((r.note || "").replace(/;?\s*needs ≥ 2 n( containing both groups)?/, ""))}</span></td></tr>`).join("")}</tbody></table></div>
    ${rows.length > limit ? `<div style="margin:10px 0"><button onclick="setP({limit:${limit + 300}})">Show ${Math.min(300, rows.length - limit)} more</button> <span class="muted small">${rows.length - limit} not shown</span></div>` : ""}`;
  $("statsTbl").addEventListener("click", ev => { const tr = ev.target.closest("tr[data-i]"); if (tr) openRelated(gen, T.objs[+tr.dataset.i]); });
};
function exportStats(gen) {
  const T = S.stats[gen], rows = S.statsView || [];
  const q = v => /[",\n;]/.test(v) ? `"${String(v).replace(/"/g, '""')}"` : v;
  const text = [T.columns.join(","), ...rows.map(r => T.columns.map(c => q(r[c] ?? "")).join(","))].join("\n") + "\n";
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], { type: "text/csv" }));
  a.download = `${S.exp}_${gen}_statistics_filtered.csv`;
  document.body.appendChild(a); a.click(); a.remove();
  toast(`exported ${rows.length} rows`);
}
// the plot a statistics row belongs to (manifest naming rules, docs/statistics.md section 6)
function relatedPlot(man, r) {
  const fam = r.family || "", ctrl = (/^vs (.+?)(?: \(|$)/.exec(fam) || [])[1];
  const pngs = man.files.filter(f => f.metric === r.metric && isPng(f));
  const folders = man.all_n_folders || [];
  const allN = pngs.filter(f => f.level === "all_n" && (folders.length < 2 || !r.cell_line || f.cell_line === r.cell_line));
  const find = (list, re) => list.find(f => re.test(f.path));
  const m = r.metric.replace(/[.*+?^${}()|[\]\\]/g, "\\$&"), c = (ctrl || "").replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
  if (/per-n (ratio|difference)/.test(fam) || /leave-one-n-out/.test(fam)) return find(allN, new RegExp(`${m}_per_n_ratio_vs_${c}\\.png$`)) || find(allN, /_per_n_ratio_vs_/);
  if (r.level === "all_n") {
    if (/density-adjusted/.test(fam)) return find(allN, new RegExp(`${m}_density\\.png$`));
    if (fam === "vehicle check (equivalence)") return find(allN, /_forest_vehicle_check\.png$/);
    if (ctrl) return find(allN, new RegExp(`${m}_forest_vs_${c}\\.png$`));
    if (fam === "estimated means") return find(allN, /_model_absolute\.png$/);
    if (/two-way/.test(fam)) return find(allN, /_twoway\.png$/);
    return find(allN, new RegExp(`/${m}\\.png$`)) || allN[0];
  }
  const inN = pngs.filter(f => f.level === "n" && f.n === r.n);
  if (/two-way/.test(fam)) return find(inN, /_twoway\.png$/);
  if (ctrl) return find(inN, new RegExp(`${m}_vs_${c}\\.png$`)) || find(inN, new RegExp(`/${m}\\.png$`));
  return find(inN, new RegExp(`/${m}\\.png$`)) || inN[0];
}
function openRelated(gen, r) {
  const f = relatedPlot(S.man[gen].manifest, r);
  if (!f) return toast("No plot for this row (it is in the tables only)", true);
  const G = `generations/${gen}/`;
  openImage(G + f.path, f.title, `${r.comparison} · ${r.family} — ${f.description || ""}`, f.pdf ? G + f.pdf : "");
}

// ---------------------------------------------------------------- Layer comparison (readout ladder)

RENDER.layers = async function (b, e, tok) {
  const gens = genList();
  if (!gens.length) { b.innerHTML = `<div class="notice">No results yet — run the analysis first.</div>`; return; }
  const gen = pickGen();
  const M = await loadManifest(gen);
  if (tok !== S.renderTok) return;
  const man = M.manifest, L = man.layers || {}, chans = Object.keys(L);
  if (!chans.length) { b.innerHTML = `<div class="notice info">No layer comparison in ${esc(gen)}: it needs a marker channel measured per nucleus (fluorescence experiments). Rebuild the generation if this is one.</div>`; return; }
  const P = S.params, ch = chans.includes(P.ch) ? P.ch : chans[0], lay = L[ch];
  const lf = man.files.filter(f => f.kind === "layers" && f.channel === ch);
  const levels = [...new Set(lf.map(f => f.level === "all_n" ? "all_n|" + dirOf(f.path) : "n|" + f.n))];
  const lvl = levels.includes(P.lvl) ? P.lvl : levels[0];
  const [lk, lv] = lvl.split("|");
  const here = lf.filter(f => lk === "all_n" ? f.level === "all_n" && dirOf(f.path) === lv : f.level === "n" && f.n === lv);
  const refs = lay.references || [];
  const ref = refs.includes(P.ref) ? P.ref : refs[0];
  const ladder = here.find(f => f.reference === ref && isPng(f)), absolute = here.find(f => /_absolute\.png$/.test(f.path));
  const conc = here.find(f => /_concordance\.csv$/.test(f.path));
  const folder = (man.all_n_folders || []).find(f => f.path === lv);
  const one = lk === "all_n" && folder && folder.n_count < 2;
  let table = null;
  if (conc) table = await get("/api/table", { exp: S.exp, path: `generations/${gen}/${conc.path}` });
  if (tok !== S.renderTok) return;
  // density scatters: intensity/density_<n>__<ch>.png and the generation's % positive density plots
  const im = S.intens || (S.intens = await get("/api/intensity", { exp: S.exp }).catch(() => ({ manifest: {} })));
  const dens = (im.manifest.files || []).filter(f => f.kind === "density" && f.channel === ch && (lk === "all_n" || f.n === lv));
  const fracM = (lay.layers.find(x => x.key === "frac_positive") || {}).metric;
  const gdens = man.files.filter(f => isPng(f) && f.metric === fracM && /_density\.png$/.test(f.path) && (lk === "all_n" ? f.level === "all_n" : f.n === lv));
  const G = `generations/${gen}/`;
  const keys = [...lay.layers.map(x => [x.key, x.label]), ...(lay.extras || []).map(x => [x.key, x.label])];
  const callCls = v => /higher/.test(v) ? "higher" : /lower/.test(v) ? "lower" : "nochange";
  const cIx = table ? Object.fromEntries(table.columns.map((c, i) => [c, i])) : {};
  const crow = table ? table.rows.filter(r => r[cIx.reference] === ref) : [];
  b.innerHTML = `
    <p class="muted" style="margin-top:0">The same ${esc(ch)} result, measured five ways: each layer removes one analysis choice. If a treatment's effect changes direction or significance between layers, the step that changes it is named below.
      Raw layers compress effects towards 100% (a fixed background sits in both conditions); field layers also scale with how many cells are in the field.</p>
    <div class="selrow">
      <label class="f">Generation ${genPills(gen)}</label>
      ${chans.length > 1 ? `<label class="f">Channel <div class="pills">${chans.map(c => `<button class="pill ${c === ch ? "active" : ""}" onclick='setP({ch:${js(c)}})'>${esc(c)}</button>`).join("")}</div></label>` : ""}
      <label class="f">Level <div class="pills">${levels.map(l => { const [a, v] = l.split("|"); return `<button class="pill ${l === lvl ? "active" : ""}" onclick='setP({lvl:${js(l)}})'>${a === "all_n" ? "All n" + (v !== "all_n" ? " · " + esc(v.split("/").pop()) : "") : `<span class="ndot" style="background:${esc(M.colours[v] || nColour(v))}"></span>${esc(shortN(v))}`}</button>`; }).join("")}</div></label>
      <label class="f">Compared with <div class="pills">${refs.map(r => `<button class="pill ${r === ref ? "active" : ""}" onclick='setP({ref:${js(r)}})'>${esc(r)}</button>`).join("")}</div></label>
    </div>
    <div class="ladder">${lay.layers.map((x, i) => `${i ? `<div class="arrow" title="${esc(x.step)}">→</div>` : ""}<div class="lstep"><b>${esc(x.label)}</b>${x.step ? `<span class="muted">${esc(x.step)}</span>` : `<span class="muted">as acquired</span>`}</div>`).join("")}
      ${(lay.extras || []).length ? `<div class="arrow">·</div>${lay.extras.map(x => `<div class="lstep" style="background:#fafafa"><b>${esc(x.label)}</b><span class="muted">extra</span></div>`).join("")}` : ""}</div>
    ${one ? `<div class="notice n1"><b>Only one n: needs ≥ 2 n.</b> All-n layer effects are estimates without intervals or stars; the per-n view shows technical (field) intervals, which are descriptive.</div>` : ""}
    ${lk === "n" ? `<div class="notice info">Within one n: intervals come from the fields (technical replicates) and are descriptive.</div>` : ""}
    <div class="k-layers"><div class="colhead">Effect vs ${esc(ref)} at every layer</div>${ladder ? plotCard(ladder, G) : `<div class="muted">no plot</div>`}</div>
    <div class="k-absolute" style="margin-top:14px"><div class="colhead">Absolute values at every layer (always read the effects next to these)</div>${absolute ? plotCard(absolute, G) : `<div class="muted">no plot</div>`}</div>
    <h2>Does the conclusion depend on the layer?</h2>
    ${crow.length ? `<div class="tablewrap" style="max-height:none"><table class="grid"><thead><tr><th>treatment</th>${keys.map(([k, l]) => `<th title="${esc(l)}">${esc(l)}</th>`).join("")}<th>changes at</th><th>in plain language</th></tr></thead>
      <tbody>${crow.map(r => `<tr class="${r[cIx.direction_agrees] === "False" ? "bad" : r[cIx.conclusion_agrees] === "False" ? "warn" : ""}"><td><b>${esc(r[cIx.treatment])}</b></td>
        ${keys.map(([k]) => `<td>${cIx[k] != null && r[cIx[k]] ? `<span class="call ${callCls(r[cIx[k]])}">${esc(r[cIx[k]])}</span>` : ""}</td>`).join("")}
        <td class="small nw">${esc((r[cIx.changes_at] || "").replace(/->/g, "→").replace(/; /g, ", ")) || `<span class="tag ok">none</span>`}</td><td class="small" style="min-width:320px">${esc(r[cIx.summary])}</td></tr>`).join("")}</tbody></table></div>
      <p class="small muted">Red rows: the direction flips between layers. Yellow rows: same direction, but the call (significance or "no change") changes.</p>`
      : `<div class="muted">No concordance table for this selection.</div>`}
    <div class="row" style="margin:8px 0 0">${here.filter(isCsv).map(f => `<button class="small" onclick='viewTable(${js("generations/" + gen + "/" + f.path)})'>${esc(f.title.replace(/^[^:]*: /, ""))}</button>`).join("")}</div>
    <h2>Density: is a difference just crowding?</h2>
    <p class="small muted" style="margin-top:0">Field layers and % positive can follow how crowded a field is. Compare the conditions at similar density before trusting a difference.</p>
    ${dens.length || gdens.length ? `<div class="plots">${gdens.map(f => plotCard(f, G)).join("")}${dens.map(f => plotCard({ ...f, kind: "design" }, "intensity/")).join("")}</div>` : `<div class="muted">No density plots.</div>`}
    <div class="row" style="margin-top:10px"><button onclick="nav('results',{gen:'${gen}',metric:'count',level:'${lk === "all_n" ? "all_n" : "n"}'${lk === "n" ? `,n:${js(lv)}` : ""}})">Cell counts for comparison →</button></div>`;
};
