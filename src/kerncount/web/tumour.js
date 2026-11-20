"use strict";
// Tumour experiments ([experiment] type = "tumour_volume"; docs/tumour.md, docs/ui.md "Tumour experiments").
// The sections of app.js / results.js stay as they are for cell-count projects; for a tumour project the home,
// layout, review, corrections, results and statistics sections are drawn here, and the settings page gets the
// analysis-plan card (factors, readouts, families) with its lock.
//
// BLINDING: the Review and Corrections sections use only /api/tumour/review (codes, object counts, QC flags, the
// review/ images) and /api/tumour/corrections. They never request the layout, the results, the overview or
// review_key.csv, and show no factor level, session, date or tumour ID.

const isTumour = e => !!e && e.type === "tumour_volume";
const T = {};

function lockLine(L) {
  if (!L) return "";
  if (L.error) return `<span class="tag fail">plan cannot be read</span> <span class="small">${esc(L.error)}</span>`;
  if (!L.results) return `<span class="tag">not locked yet</span> <span class="small muted">plan ${esc(L.settings_lock)} — locked as soon as results exist</span>`;
  const changed = L.changed ? ` <span class="tag medium" title="settings.toml has another plan than the results">settings: plan ${esc(L.settings_lock)} — run to apply</span>` : "";
  if (L.unlocked) return `<span class="tag medium">unlocked</span> <span class="small">${esc(L.unlocked.at)}: ${esc(L.unlocked.reason)}</span>${changed}`;
  return `<span class="tag ok">locked</span> <span class="small muted">results built with plan ${esc(L.results_lock)}</span>${changed}`;
}
function attentionList(att) {
  return att.length ? att.map(a => `<div class="att ${a.level}"><div class="ico">${a.level === "info" ? "i" : "!"}</div>
      <div><div class="t1">${esc(a.title)}</div><div class="t2">${esc(a.text)}</div><div class="ev">${a.evidence.map(evidenceBtn).join("")}</div></div></div>`).join("")
    : `<div class="notice ok">Nothing needs your attention: the layout is complete and the plan is locked as recorded.</div>`;
}

// ---------------------------------------------------------------- experiment home

T.home = async function (b, e, tok) {
  if (!e.initialised) { b.innerHTML = runPanel(e); return; }
  const ov = await get("/api/tumour/overview", { exp: S.exp });
  if (tok !== S.renderTok) return;
  S.tov = ov;
  const R = ov.results;
  b.innerHTML = `
    <div class="card exp-head">
      <div class="title"><h1>${esc(S.exp)}</h1><div class="muted small">${esc(S.project ? S.project.path : "")}</div>
        <div style="margin-top:6px"><span class="tag purple">tumour experiment</span><span class="tag">${ov.included} tumour(s) in the layout</span>
          ${ov.tumours ? `<span class="tag blue">${ov.tumours} analysed</span>` : ""}</div></div>
      <div style="flex:2;min-width:420px"><div class="kv">
        ${ov.factors.map(f => `<div><div class="k">${esc(f.name)} (${esc(f.source)})</div>${metaVal(f.levels)}</div>`).join("")}
        <div><div class="k">Raw folder</div>${metaVal(ov.meta.raw_dir)}</div>
        <div><div class="k">Review generation</div>${metaVal(ov.generation)}</div>
      </div><div class="small" style="margin-top:8px"><a class="link" ${on("click", function () { nav("settings") })}>Settings and analysis plan →</a>${ov.meta.notes ? ` · <span class="muted">${esc(ov.meta.notes)}</span>` : ""}</div></div>
    </div>
    <h2>Attention ${ov.attention.length ? `<span class="tag ${ov.attention.some(a => a.level === "error") ? "fail" : ov.attention.some(a => a.level === "warn") ? "medium" : "blue"}">${ov.attention.length}</span>` : ""}</h2>
    ${attentionList(ov.attention)}
    <div class="cols2" style="margin-top:14px">
      <div class="card"><h3>Groups (included tumours)</h3>
        ${ov.groups.length ? `<table class="grid"><thead><tr><th>group</th><th class="num">tumours</th></tr></thead><tbody>${ov.groups.map(g => `<tr><td>${esc(g.group)}</td><td class="num">${g.tumours}</td></tr>`).join("")}</tbody></table>`
          : `<div class="muted">No readable tumours yet.</div>`}
        <p class="small muted" style="margin-bottom:0">One tumour (hemisphere A or B) = one n. ${ov.excluded_layout} excluded in the layout. <a class="link" ${on("click", function () { nav("layout") })}>Layout →</a></p></div>
      <div class="card"><h3>Results</h3>
        ${R.built ? `<p style="margin:0 0 6px">Built ${esc(R.built)} · <b>${ov.tumours}</b> tumours analysed · ${ov.flagged} QC-flagged · ${ov.excluded} excluded (layout or review)</p>
          <div class="row"><button class="small" ${on("click", function () { nav("results") })}>Results →</button><button class="small" ${on("click", function () { nav("stats") })}>Statistics →</button>
            ${R.xlsx ? `<a class="btn small" href="${fileUrl(R.xlsx, 1)}">results.xlsx</a>` : ""}${R.report ? `<button class="small" ${on("click", function () { viewText(R.report) })}>Report</button>` : ""}</div>`
          : `<div class="muted">No results yet: run the analysis.</div>`}
        <h3>Analysis plan</h3><div>${lockLine(ov.lock)}</div>
        <p class="small muted" style="margin-bottom:0">Factors, readout roles (confirmatory / exploratory) and contrast families are fixed before the results are viewed (ADR 0003). <a class="link" ${on("click", function () { nav("settings") })}>Analysis plan →</a></p></div>
    </div>
    <h2 id="run">Run</h2>${tumourRunPanel(e, ov)}`;
  if (S.params.anchor) { const a = $(S.params.anchor); if (a) a.scrollIntoView(); }
};

function tumourRunPanel(e, ov) {
  const step = (n, done, title, text, controls) => `<div class="step ${done ? "done" : ""}"><div><span class="num">${done ? "✓" : n}</span><b>${title}</b></div><p>${text}</p>${controls}</div>`;
  const d = S.running ? "disabled" : "";
  return `${S.running ? `<div class="notice info">A job is running — see the log at the bottom. Buttons that change files are paused until it finishes.</div>` : ""}
    <div class="steps">
      ${step(1, e.initialised, "Set up", "Scan the raw folder (its session folders, as delivered) into layout.csv: one row per tumour.",
        `<label class="small"><input type="checkbox" id="o_force"> start over: reset settings and layout</label>
         <label class="small"><input type="checkbox" id="o_forcelayout"> rewrite the layout only</label>
         <button ${d} ${on("click", function () { setUp($("o_force").checked, $("o_forcelayout").checked); })}>Re-scan folder</button>`)}
      ${step(2, !!ov.results.built, "Run the analysis", "Segment every stack, measure the readouts, QC, draw the blinded review sheets and build the statistics. Segmentations are cached: re-runs are quick.",
        `<label><input type="checkbox" id="o_noreview"> skip review images (faster)</label>
         <label><input type="checkbox" id="o_reseg"> segment every stack again (ignore the cache)</label>
         <label><input type="checkbox" id="o_sweep" ${ov.threshold_sweep ? "checked" : ""}> threshold sweep (robustness of every confirmatory effect; slow first run)</label>
         <button class="primary" ${d} ${on("click", function () { runTumour() })}>Run analysis</button>`)}
      ${step(3, ov.generation && ov.generation !== "g1", "Review and correct (blinded)", "Look at the review sheets by code; re-pick a main mass or exclude a tumour in Corrections. The next run applies the corrections (g2, g3, …).",
        `<button ${ov.review ? "" : "disabled"} ${on("click", function () { nav("review") })}>Review sheets →</button>
         <button ${ov.review ? "" : "disabled"} ${on("click", function () { nav("corrections") })}>Corrections →</button>`)}
    </div>
    <h3>On this computer</h3>
    <div class="row"><button ${on("click", function () { openSys("data") })}>Open raw folder</button>
      <button ${on("click", function () { openSys("folder") })}>Open project folder</button>
      <button ${ov.results.xlsx ? "" : "disabled"} ${on("click", function () { openSys("results") })}>Open results in Excel</button></div>`;
}
async function runTumour() {
  const sweep = $("o_sweep").checked;
  if (S.tov && sweep !== !!S.tov.threshold_sweep) {  // a setting (not part of the locked plan), saved before the run
    await post("/api/settings", { experiment: S.exp, values: { tumour: { threshold_sweep: sweep } } }, true);
    S.project = null;
  }
  startJob("run", { no_review: $("o_noreview").checked, resegment: $("o_reseg").checked });
}

// ---------------------------------------------------------------- Layout

T.layout = async function (b) {
  const L = await get("/api/tumour/layout", { exp: S.exp });
  S.tlay = L;
  const ix = Object.fromEntries(L.columns.map((c, i) => [c, i]));
  const hasNew = S.project && S.project.files["layout.new.csv"];
  const sel = (c, i, opts, v) => `<select data-c="${esc(c)}">${(opts.includes(v) ? [""] : ["", v]).concat(opts).filter((o, k, a) => a.indexOf(o) === k)
    .map(o => `<option value="${esc(o)}" ${o === v ? "selected" : ""}>${o === "" ? "–" : esc(o) + (opts.includes(o) ? "" : " ?")}</option>`).join("")}</select>`;
  b.innerHTML = `
    ${hasNew ? `<div class="notice">A fresh scan of the raw folder was saved as <code>layout.new.csv</code>.
       <button class="small" ${on("click", function () { viewTable("layout.new.csv") })}>Compare</button>
       <button class="small" ${on("click", function () { useNewLayout() })}>Use the fresh scan</button></div>` : ""}
    ${L.problems.length ? `<div class="notice err"><b>${L.problems.length} problem(s) block the run</b> (nothing is measured until every included file maps to exactly one tumour ID):<ul>${L.problems.map(p => `<li>${esc(p)}</li>`).join("")}</ul></div>`
      : `<div class="notice ok">The layout is complete: ${L.included} tumour(s) included, ${L.excluded} excluded.</div>`}
    ${L.warnings.length ? `<div class="notice"><ul style="margin:0">${L.warnings.map(w => `<li>${esc(w)}</li>`).join("")}</ul></div>` : ""}
    <p class="muted">One row per original z-stack. Factor values come from the file name or the session-folder name; fix a value that was read wrongly
      (it is noted as an override), or untick <b>include</b> and give the reason in <b>note</b>. The tumour ID is rebuilt from the values when you save.</p>
    <div class="bar"><button class="primary" ${on("click", function () { saveTumourLayout() })}>Save layout</button><span class="muted small" id="layoutState"></span></div>
    <div class="tablewrap"><table class="grid" id="tlayTbl"><thead><tr><th>file</th><th>session folder</th>${L.factors.map(f => `<th>${esc(f.name)}</th>`).join("")}
      <th>section</th><th>hemi&shy;sphere</th><th>tumour ID</th><th>include</th><th>note</th></tr></thead>
    <tbody>${L.rows.map((r, i) => `<tr data-i="${i}">
      <td class="small" title="${esc(r[ix.file])}">${esc(r[ix.file])}</td><td class="small">${esc(r[ix.imaging_session])}</td>
      ${L.factors.map(f => `<td>${sel(f.name, i, f.levels, r[ix[f.name]] ?? "")}</td>`).join("")}
      <td><input type="text" data-c="section" size="4" value="${esc(r[ix.section] ?? "")}"></td>
      <td>${sel("hemisphere", i, ["A", "B"], r[ix.hemisphere] ?? "")}</td>
      <td class="small">${L.ids[i] ? esc(L.ids[i]) : `<span class="tag fail">none</span>`}</td>
      <td><input type="checkbox" data-c="include" ${/^(yes|y|true|1|ja|j)$/i.test(r[ix.include] || "") ? "checked" : ""}></td>
      <td><input type="text" data-c="note" value="${esc(r[ix.note] ?? "")}" style="width:100%;min-width:220px"></td></tr>`).join("")}</tbody></table></div>
    ${L.unmatched.length ? `<details class="card"><summary>${L.unmatched.length} other file(s) in the raw folder (neither a tumour nor a matched Fiji output; ignored)</summary><ul>${L.unmatched.map(f => `<li class="small">${esc(f)}</li>`).join("")}</ul></details>` : ""}`;
  b.querySelectorAll("#tlayTbl input, #tlayTbl select").forEach(el => el.addEventListener("input", () => { $("layoutState").textContent = "unsaved changes"; el.closest("tr").classList.add("changed"); }));
};
async function saveTumourLayout() {
  const L = S.tlay, ix = Object.fromEntries(L.columns.map((c, i) => [c, i]));
  const rows = [...document.querySelectorAll("#tlayTbl tbody tr")].map(tr => {
    const r = L.rows[+tr.dataset.i], out = { file: r[ix.file] };
    tr.querySelectorAll("[data-c]").forEach(el => { out[el.dataset.c] = el.type === "checkbox" ? (el.checked ? "yes" : "no") : el.value; });
    return out;
  });
  const r = await post("/api/tumour/layout", { experiment: S.exp, rows });
  toast(r.message); S.tov = null; render();
}

// ---------------------------------------------------------------- blinded review and corrections (codes only)

function tRounds(R) {  // like counts.csv: the chain uses manual_1, manual_2, ... up to the first empty round
  const M = R.manual;
  let first = M.findIndex(c => !R.rows.some(r => (r.manual[c] || "").trim()));
  if (first < 0) first = M.length;
  return { active: M.slice(0, first), newCol: `manual_${first + 1}` };
}
function tValue(r, col) { const k = `${r.code}\u0001${col}`; return k in S.pending ? S.pending[k] : col === "note" ? r.note : (r.manual[col] || ""); }
function tCurrent(r) {  // what the next run applies: the last non-empty round of the chain
  let v = "";
  for (const c of tRounds(S.trev).active) { const x = String(tValue(r, c)).trim(); if (x) v = x; }
  return v;
}
const tValid = v => /^\s*$|^\s*\d+(\.0+)?\s*$|^\s*x\b/i.test(String(v));
function tShow(v) { return !v ? "auto" : /^x/i.test(v) ? "excluded" + (v.replace(/^x[\s:;,\-–]*/i, "") ? ": " + v.replace(/^x[\s:;,\-–]*/i, "") : "") : "object " + v; }
async function loadTReview(force) {
  if (!S.trev || force) S.trev = await get("/api/tumour/review", { exp: S.exp });
  const rr = tRounds(S.trev);
  if (!S.tround || (S.tround !== rr.newCol && !rr.active.includes(S.tround))) S.tround = rr.newCol;
  return S.trev;
}
const BLIND_TXT = `<div class="notice info"><b>Blinded.</b> Tumours are shown by a random code only, in code order, each image contrast-stretched on its own.
  Nothing here names a group, session folder, date or file (the key <code>review_key.csv</code> stays out of these pages). Orange = main mass (object 1), cyan = satellites; the numbers are the objects you can pick.</div>`;

T.review = async function (b) {
  const R = await loadTReview(true);
  if (!R.rows.length) { b.innerHTML = `<div class="notice">No review sheets yet — run the analysis (without “skip review images”).</div>`; return; }
  const only = S.params.q === "flagged";
  const rows = R.rows.filter(r => !only || r.qc_flags);
  S.tlist = rows;
  b.innerHTML = `${BLIND_TXT}
    <div class="bar"><div class="pills">
      <button class="pill ${only ? "" : "active"}" ${on("click", function () { setP({ q: "" }) })}>all (${R.rows.length})</button>
      <button class="pill ${only ? "active" : ""}" ${on("click", function () { setP({ q: "flagged" }) })}>QC-flagged (${R.rows.filter(r => r.qc_flags).length})</button></div>
      <span class="spacer"></span>${R.overviews.map((p, i) => `<button class="small" ${on("click", function () { openImage(p, "Overview sheet " + (i + 1)) })}>Overview ${i + 1}</button>`).join("")}
      <button class="small" ${on("click", function () { nav("corrections") })}>Corrections →</button></div>
    <div class="gallery small">${rows.map((r, j) => `<div class="thumb ${r.qc_flags ? "q-medium" : ""}" ${on("click", function () { tLightbox(j) })}>
      ${r.image ? `<img loading="lazy" src="${fileUrl(r.image)}">` : `<div class="muted" style="padding:30px">no image</div>`}
      <div class="cap"><b>${esc(r.code)}</b><span class="small">${esc(r.objects)} object(s)</span>${tCurrent(r) ? `<span class="tag medium">${esc(tShow(tCurrent(r)))}</span>` : ""}</div>
      ${r.qc_flags ? `<div class="src">${esc(r.qc_flags)}</div>` : ""}</div>`).join("")}</div>`;
};
function tLightbox(j) {
  S.tlb = j;
  const r = S.tlist[j]; if (!r) return closeModal();
  const rr = tRounds(S.trev), col = S.tround;
  openModal(`<div class="modal-head"><span class="title">${esc(r.code)}</span><span>${esc(r.objects)} object(s)</span>${r.qc_flags ? `<span class="tag medium">${esc(r.qc_flags)}</span>` : ""}
      <span class="spacer"></span><button ${j ? "" : "disabled"} id="tlbPrev">‹ Prev</button><button ${j < S.tlist.length - 1 ? "" : "disabled"} id="tlbNext">Next ›</button><button id="tlbClose">Close</button></div>
    ${r.image ? `<img class="full" src="${fileUrl(r.image)}">` : `<div class="muted">no review image</div>`}
    <div class="row" style="margin-top:10px"><span>Next run: <b>${esc(tShow(tCurrent(r)))}</b></span>
      <label>Round <select id="tlbRound">${rr.active.concat([rr.newCol]).map(c => `<option ${c === col ? "selected" : ""}>${esc(c)}</option>`).join("")}</select></label>
      <input id="tlbVal" size="14" value="${esc(tValue(r, col))}" placeholder="object no. or x: reason">
      <input id="tlbNote" size="24" value="${esc(tValue(r, "note"))}" placeholder="note">
      <button class="primary" id="tlbSave">Save</button><span class="muted small">Enter = save &amp; next. The next run applies it.</span></div>`);
  $("tlbPrev").onclick = () => tLightbox(j - 1);
  $("tlbNext").onclick = () => tLightbox(j + 1);
  $("tlbClose").onclick = () => closeModal();
  $("tlbRound").onchange = ev => { S.tround = ev.target.value; tLightbox(j); };
  const save = async next => {
    const v = $("tlbVal").value.trim(), n = $("tlbNote").value.trim();
    if (!tValid(v)) return toast("Type an object number, x: reason, or leave it empty", true);
    const changes = [];
    if (v !== String(r.manual[col] || "").trim()) changes.push({ code: r.code, col, value: v });
    if (n !== String(r.note || "").trim()) changes.push({ code: r.code, col: "note", value: n });
    if (changes.length) {
      await post("/api/tumour/corrections", { experiment: S.exp, changes, mtime: S.trev.mtime });
      await loadTReview(true);
      S.tlist = S.tlist.map(x => S.trev.rows.find(y => y.code === x.code) || x);
      toast(`saved ${r.code}`);
    }
    tLightbox(next && j < S.tlist.length - 1 ? j + 1 : j);
  };
  $("tlbSave").onclick = () => save(false);
  $("tlbVal").onkeydown = ev => { if (ev.key === "Enter") { ev.preventDefault(); save(true); } ev.stopPropagation(); };
  $("tlbNote").onkeydown = ev => { if (ev.key === "Enter") { ev.preventDefault(); save(true); } ev.stopPropagation(); };
  $("tlbVal").focus();
}

T.corrections = async function (b) {
  const R = await loadTReview(false);
  if (!R.rows.length) { b.innerHTML = `<div class="notice">No review/corrections.csv yet — run the analysis (with review images) first.</div>`; return; }
  const rr = tRounds(R), col = S.tround, k = +col.slice(7);
  const prev = rr.active.filter(c => +c.slice(7) < k);
  const nPending = Object.keys(S.pending).length;
  b.innerHTML = `${BLIND_TXT}
    <p class="muted" style="margin-top:0">Each round (<code>manual_1</code>, <code>manual_2</code>, …) builds the next generation: g1 = automatic, g1 + round 1 → <b>g2</b>, …
      Type the <b>number of an object</b> on the review image to make it the main mass (1 = back to automatic), or <b>x: reason</b> to exclude the tumour (the reason goes to
      <code>excluded_tumours.csv</code>); empty = no change. Per-tumour thresholds are not possible (one rule for every tumour). The next run applies the corrections without re-segmenting.</p>
    <div class="bar">
      <select ${on("change", function () { tChangeRound(this) })}>${rr.active.map(c => `<option value="${c}" ${c === col ? "selected" : ""}>Revise round ${c.slice(7)} (${c} → g${+c.slice(7) + 1})</option>`).join("")}
        <option value="${rr.newCol}" ${rr.newCol === col ? "selected" : ""}>New round ${rr.newCol.slice(7)} (${rr.newCol} → g${+rr.newCol.slice(7) + 1})</option></select>
      <button class="primary" ${nPending ? "" : "disabled"} ${on("click", function () { saveTCorrections(false) })}>Save ${nPending || ""} change${nPending === 1 ? "" : "s"}</button>
      <button ${S.running ? "disabled" : ""} ${on("click", function () { saveTCorrections(true) })}>${nPending ? "Save &amp; " : ""}run the analysis</button>
      <span class="spacer"></span>
      <label class="small"><input type="checkbox" ${S.params.q === "flagged" ? "checked" : ""} ${on("change", function () { setP({ q: this.checked ? "flagged" : "" }) })}> only QC-flagged</label></div>
    <div class="split">
      <div class="tablewrap"><table class="grid" id="tcorrTbl"><thead><tr><th>code</th><th class="num">objects</th><th>QC flags</th>
        ${prev.map(c => `<th class="small">${esc(c)}</th>`).join("")}<th style="width:140px">${esc(col)}</th><th>note</th><th>→ g${k + 1}</th></tr></thead>
      <tbody>${R.rows.filter(r => S.params.q !== "flagged" || r.qc_flags).map(r => { const v = tValue(r, col), key = `${r.code}\u0001${col}`;
        return `<tr data-code="${esc(r.code)}" class="${key in S.pending || `${r.code}\u0001note` in S.pending ? "changed" : ""}">
          <td><b>${esc(r.code)}</b></td><td class="num">${esc(r.objects)}</td><td class="small">${esc(r.qc_flags)}</td>
          ${prev.map(c => `<td class="small">${esc(tValue(r, c))}</td>`).join("")}
          <td><input type="text" data-col="${esc(col)}" data-code="${esc(r.code)}" class="${tValid(v) ? "" : "bad"}" value="${esc(v)}" size="12"></td>
          <td><input type="text" data-col="note" data-code="${esc(r.code)}" value="${esc(tValue(r, "note"))}" style="width:100%;min-width:140px"></td>
          <td class="small" data-after>${esc(tShow(tAfter(r, k)))}</td></tr>`; }).join("")}</tbody></table></div>
      <div class="preview" id="tcorrPreview"><div class="muted">Click a row to see its review image.</div></div></div>`;
  const tbl = $("tcorrTbl");
  tbl.addEventListener("click", ev => { const tr = ev.target.closest("tr[data-code]"); if (tr) tPreview(tr.dataset.code); });
  tbl.addEventListener("focusin", ev => { if (ev.target.dataset.code) tPreview(ev.target.dataset.code); });
  tbl.addEventListener("input", ev => {
    const inp = ev.target; if (!inp.dataset.code) return;
    const r = R.rows.find(x => x.code === inp.dataset.code), c = inp.dataset.col, key = `${r.code}\u0001${c}`;
    const orig = c === "note" ? r.note : (r.manual[c] || "");
    if (inp.value.trim() === String(orig).trim()) delete S.pending[key]; else S.pending[key] = inp.value.trim();
    if (c !== "note") inp.classList.toggle("bad", !tValid(inp.value));
    const tr = inp.closest("tr");
    tr.classList.toggle("changed", `${r.code}\u0001${col}` in S.pending || `${r.code}\u0001note` in S.pending);
    tr.querySelector("[data-after]").textContent = tShow(tAfter(r, k));
    const n = Object.keys(S.pending).length, btn = b.querySelector(".bar .primary");
    btn.disabled = !n; btn.textContent = `Save ${n || ""} change${n === 1 ? "" : "s"}`;
  });
  if (S.tsel) tPreview(S.tsel);
};
function tAfter(r, k) {  // value after rounds 1..k (pending edits included)
  let v = "";
  for (let i = 1; i <= k; i++) { const x = String(tValue(r, `manual_${i}`)).trim(); if (x) v = x; }
  return v;
}
function tPreview(code) {
  S.tsel = code;
  document.querySelectorAll("#tcorrTbl tr.sel").forEach(tr => tr.classList.remove("sel"));
  const tr = document.querySelector(`#tcorrTbl tr[data-code="${CSS.escape(code)}"]`); if (tr) tr.classList.add("sel");
  const r = S.trev.rows.find(x => x.code === code); if (!r || !$("tcorrPreview")) return;
  $("tcorrPreview").innerHTML = `<div class="row" style="margin-bottom:6px"><b>${esc(r.code)}</b><span class="small muted">${esc(r.objects)} object(s)</span>${r.qc_flags ? `<span class="tag medium">${esc(r.qc_flags)}</span>` : ""}</div>
    ${r.image ? `<img src="${fileUrl(r.image)}" ${on("click", function () { openImage(r.image, r.code) })}>` : `<div class="muted">No review image (run without “skip review images”).</div>`}
    <div class="small muted" style="margin-top:6px">Left: maximum projection over all layers. Right: the mid-section centre layer. Click to enlarge.</div>`;
}
function tChangeRound(sel) {
  const apply = () => { S.pending = {}; S.tround = sel.value; render(); };
  if (!Object.keys(S.pending).length) return apply();
  const want = sel.value; sel.value = S.tround;
  confirmBox("Discard your unsaved corrections in this round?", () => { sel.value = want; apply(); }, "Discard");
}
async function saveTCorrections(thenRun) {
  const changes = Object.entries(S.pending).map(([k, value]) => { const [code, col] = k.split("\u0001"); return { code, col, value }; });
  if (changes.some(c => c.col !== "note" && !tValid(c.value))) return toast("Fix the red cells: an object number, x: reason, or empty", true);
  if (changes.length) {
    await post("/api/tumour/corrections", { experiment: S.exp, changes, mtime: S.trev.mtime });
    S.pending = {};
    await loadTReview(true);
    toast(`${changes.length} change(s) saved`);
  }
  if (thenRun) return startJob("run", {});
  render();
}

// ---------------------------------------------------------------- results and statistics

async function loadTStats(sens) {
  const key = sens ? "all" : "head";
  S.tstats = S.tstats || {};
  if (!S.tstats[key]) S.tstats[key] = await get("/api/tumour/stats", { exp: S.exp, sensitivity: sens ? "1" : "0" });
  return S.tstats[key];
}
function tStatsTable(st, rows, withMetric) {
  const c = Object.fromEntries(st.columns.map((x, i) => [x, i]));
  const ci = r => r[c.ci_low] !== "" && r[c.ci_low] != null ? `${fmtNum(r[c.ci_low])}–${fmtNum(r[c.ci_high])}` : "";
  return `<div class="tablewrap" style="max-height:none"><table class="grid"><thead><tr>${withMetric ? "<th>readout</th>" : ""}<th>family</th><th>comparison</th><th class="num">estimate</th><th class="num">95% CI</th>
      <th class="num">p</th><th class="num">p (Holm)</th><th></th><th class="num">N</th><th>note</th></tr></thead>
    <tbody>${rows.map(r => `<tr class="${r[c.significance] && r[c.significance] !== "ns" ? "sig" : ""}">${withMetric ? `<td class="small">${esc(r[c.metric])}</td>` : ""}
      <td class="small">${esc(r[c.family])}</td><td>${esc(r[c.comparison])}</td><td class="num">${esc(fmtEst(r[c.estimate], r[c.estimate_type]))}</td>
      <td class="num small">${esc(ci(r))}</td><td class="num">${esc(fmtP(r[c.p]))}</td><td class="num">${esc(fmtP(r[c.p_adj]))}</td>
      <td>${sigBadge(r[c.significance])}</td><td class="num">${esc(r[c.N])}</td><td class="small muted">${esc(r[c.note])}</td></tr>`).join("")}</tbody></table></div>`;
}

T.results = async function (b, e, tok) {
  const R = await get("/api/tumour/results", { exp: S.exp });
  if (tok !== S.renderTok) return;
  if (!R.built) { b.innerHTML = `<div class="notice">No results yet — run the analysis first.</div>`; return; }
  const ro = R.readouts.find(r => r.name === S.params.ro) || R.readouts[0];
  const st = await loadTStats(false);
  if (tok !== S.renderTok) return;
  const mi = st.columns.indexOf("metric"), rows = ro ? st.rows.filter(r => r[mi] === ro.name) : [];
  b.innerHTML = `
    <div class="row" style="margin-bottom:8px">${lockLine(R.lock)}<span class="small muted">built ${esc(R.built)}</span><span class="spacer"></span>
      ${R.xlsx ? `<a class="btn primary" href="${fileUrl(R.xlsx, 1)}">Download results.xlsx</a>` : ""}
      ${R.report ? `<button ${on("click", function () { viewText(R.report) })}>Report</button>` : ""}</div>
    <div class="bar"><div class="pills">${R.readouts.map(r => `<button class="pill ${ro && r.name === ro.name ? "active" : ""}" title="${esc(r.column)}" ${on("click", function () { setP({ ro: r.name }) })}>${esc(r.label)}${r.role === "confirmatory" ? ` <span class="cnt">confirmatory</span>` : ""}</button>`).join("")}</div></div>
    ${ro ? `<p class="muted small">${ro.role === "confirmatory" ? "<b>Confirmatory readout</b>: tested with Holm correction within each family; the only readout that carries stars (confirmatory families)."
      : "<b>Exploratory readout</b>: effect sizes with CIs (p-values labelled exploratory), never stars."} One tumour = one n.</p>
      <div class="plots">${ro.plots.map(f => plotCard(f)).join("") || `<div class="muted">No plots for this readout.</div>`}</div>
      <h3>Statistics of this readout <span class="small muted">(sensitivity analyses: Statistics page)</span></h3>${rows.length ? tStatsTable(st, rows, false) : `<div class="muted">No statistics rows.</div>`}` : ""}
    ${R.extra.length ? `<h3>More plots</h3><div class="plots">${R.extra.map(f => plotCard(f)).join("")}</div>` : ""}
    <h3>Tables</h3><div class="tbllist">${R.tables.map(t => `<button class="btn" ${on("click", function () { viewTable(t) })}>${esc(t)}</button>`).join("")}</div>
    ${R.prism.length ? `<details class="card"><summary>Prism-ready tables (${R.prism.length})</summary>${R.prism.map(p => `<div class="small"><a class="link" ${on("click", function () { viewTable(p) })}>${esc(p)}</a> · <a href="${fileUrl(p, 1)}">download</a></div>`).join("")}</details>` : ""}`;
};

T.stats = async function (b, e, tok) {
  const P = S.params, sens = P.sens === "1";
  const st = await loadTStats(sens);
  if (tok !== S.renderTok) return;
  if (!st.columns.length) { b.innerHTML = `<div class="notice">No statistics yet — run the analysis first.</div>`; return; }
  const c = Object.fromEntries(st.columns.map((x, i) => [x, i]));
  const metrics = [...new Set(st.rows.map(r => r[c.metric]))], fams = [...new Set(st.rows.map(r => r[c.family]))];
  const q = (P.q || "").toLowerCase();
  const rows = st.rows.filter(r => (!P.metric || r[c.metric] === P.metric) && (!P.fam || r[c.family] === P.fam)
    && (P.sig !== "1" || (r[c.significance] && r[c.significance] !== "ns")) && (!q || r.join(" ").toLowerCase().includes(q)));
  b.innerHTML = `<p class="muted" style="margin-top:0">Every estimate of <code>stats_long.csv</code>. Headline view: the pre-registered model only;
      the sensitivity analyses (families ending in <i>(Welch-type)</i>, <i>(permutation)</i>, <i>(without flagged)</i>) are hidden${st.hidden ? ` (${st.hidden} rows)` : ""} unless you show them.</p>
    <div class="bar">
      <select ${on("change", function () { setP({ metric: this.value }) })}><option value="">all readouts</option>${metrics.map(m => `<option ${m === P.metric ? "selected" : ""}>${esc(m)}</option>`).join("")}</select>
      <select ${on("change", function () { setP({ fam: this.value }) })}><option value="">all families</option>${fams.map(f => `<option ${f === P.fam ? "selected" : ""}>${esc(f)}</option>`).join("")}</select>
      <label class="small"><input type="checkbox" ${P.sig === "1" ? "checked" : ""} ${on("change", function () { setP({ sig: this.checked ? "1" : "" }) })}> significant only</label>
      <label class="small"><input type="checkbox" ${sens ? "checked" : ""} ${on("change", function () { setP({ sens: this.checked ? "1" : "", fam: "" }) })}> show sensitivity rows</label>
      <input placeholder="search" value="${esc(P.q || "")}" ${on("change", function () { setP({ q: this.value }) })}>
      <span class="spacer"></span><span class="small muted">${rows.length} row(s)</span><a class="btn small" href="${fileUrl("stats_long.csv", 1)}">stats_long.csv</a></div>
    ${tStatsTable(st, rows, true)}`;
};

// ---------------------------------------------------------------- settings: the analysis plan and its lock

async function tumourPlanCard(b) {
  const P = await get("/api/tumour/plan", { exp: S.exp });
  S.tplan = P;
  const L = P.lock, dis = L.locked ? "disabled" : "";
  const opts = (list, v, empty) => (empty ? [""] : []).concat(list).map(o => `<option value="${esc(o)}" ${o === v ? "selected" : ""}>${o === "" ? "–" : esc(o)}</option>`).join("");
  const box = document.createElement("div");
  box.className = "form-sec"; box.id = "sec_plan";
  box.innerHTML = `<h3>Analysis plan <span class="key small muted">[[tumour.factors]] · [[tumour.readouts]] · [[tumour.families]]</span></h3>
    <p class="small muted">Fixed <b>before the results are viewed</b> (ADR 0003): the factors with their levels and reference levels, which readouts are confirmatory
      (tested with Holm correction, stars) or exploratory (effect sizes, labelled), and the contrast families. The run records a hash of this plan in every
      statistics row and the report. Once results exist the plan is <b>locked</b>; changing it needs an explicit unlock with a reason, which is logged in
      <code>lock_log.csv</code> and printed in the report next to the results. Display roles only place factors in the plots and can change at any time.</p>
    <div class="row" style="margin-bottom:8px">${lockLine(L)}<span class="spacer"></span>
      ${L.locked ? `<button ${on("click", function () { unlockPlan() })}>Unlock…</button>` : L.unlocked ? `<button ${on("click", function () { relockPlan() })}>Lock again</button>` : ""}</div>
    <h4>Factors</h4>
    <table class="grid" id="planFactors"><thead><tr><th>factor</th><th>read from</th><th>levels (comma-separated)</th><th>reference level</th><th>display role</th><th>other spellings</th></tr></thead>
    <tbody>${P.factors.map(f => `<tr data-name="${esc(f.name)}"><td><b>${esc(f.name)}</b></td><td class="small">${esc(f.source === "session" ? "session folder" : "file name")}</td>
      <td><input type="text" data-k="levels" value="${esc(f.levels.join(", "))}" ${dis} style="width:100%"></td>
      <td><select data-k="reference" ${dis}>${opts(f.levels, f.reference)}</select></td>
      <td><select data-k="role">${opts(P.roles.factor, f.role, true)}</select></td>
      <td class="small muted">${esc(Object.entries(f.aliases).map(([k, v]) => `${k} → ${v}`).join(", "))}</td></tr>`).join("")}</tbody></table>
    <p class="small muted">Display role: row / panel / x = where the factor goes in the group plots. Levels and reference levels are part of the plan.</p>
    <h4>Readouts</h4>
    <table class="grid" id="planReadouts"><thead><tr><th>readout</th><th>column in tumours.csv</th><th>role</th><th>Fiji counterpart</th></tr></thead>
    <tbody>${P.readouts.map(r => `<tr data-name="${esc(r.name)}"><td>${esc(r.label)} <span class="key small muted">${esc(r.name)}</span></td><td class="small"><code>${esc(r.column)}</code></td>
      <td><select data-k="role" ${dis}>${opts(P.roles.readout, r.role)}</select></td><td class="small muted">${esc(r.legacy)}</td></tr>`).join("")}</tbody></table>
    <p class="small muted">Confirmatory = main-mass cross-sectional area and perimeter over the mid-section (ADR 0002/0003 amendments): only ~80–150 µm of each
      tumour is imaged, so a 3D volume depends on the acquired depth and is exploratory.</p>
    <h4>Contrast families <span class="small muted">(Holm within each family and readout)</span></h4>
    <table class="grid" id="planFamilies"><thead><tr><th>id</th><th>effect</th><th>factors</th><th>what it compares</th><th>role</th></tr></thead>
    <tbody>${P.families.map(f => `<tr data-name="${esc(f.id)}"><td><b>${esc(f.id)}</b></td><td class="small">${esc(f.effect)}</td><td class="small">${esc(f.factors.join(" × "))}</td>
      <td class="small">${esc(f.label)}</td><td><select data-k="role" ${dis}>${opts(P.roles.family, f.role)}</select></td></tr>`).join("")}</tbody></table>
    <p class="small muted">confirmatory = stars (confirmatory readouts only); exploratory = p-values, no stars; descriptive = estimate and CI only, never tested.
      Batches (setting <code>batch</code>: culture batch = file-name date by default, or session folder, or acquisition day) that share a group enter the model as a fixed block effect; contrasts across unconnected batches (e.g. WT vs KO, prepared in different culture batches) are marked confounded and never tested. Add or remove factors, readouts or
      families in the settings.toml text below.</p>
    <div class="row" style="margin-top:8px"><button class="primary" ${on("click", function () { savePlan() })}>Save analysis plan</button><span class="muted small" id="planState"></span></div>
    ${L.entries.length ? `<h4>Lock history (lock_log.csv)</h4><table class="grid"><thead><tr><th>when</th><th>what</th><th>results plan</th><th>settings plan</th><th>reason</th></tr></thead>
      <tbody>${L.entries.map(x => `<tr><td class="small nw">${esc(x.at)}</td><td>${esc(x.action)}</td><td class="small"><code>${esc(x.results_lock)}</code></td><td class="small"><code>${esc(x.settings_lock)}</code></td><td class="small">${esc(x.reason)}</td></tr>`).join("")}</tbody></table>` : ""}`;
  b.prepend(box);
  box.querySelectorAll("input, select").forEach(el => el.addEventListener("input", () => { $("planState").textContent = "unsaved changes"; }));
}
async function savePlan() {
  const rows = id => [...document.querySelectorAll(`#${id} tbody tr`)].map(tr => {
    const o = { name: tr.dataset.name };
    tr.querySelectorAll("[data-k]").forEach(el => { o[el.dataset.k] = el.value; });
    return o;
  });
  const fams = rows("planFamilies").map(f => ({ id: f.name, role: f.role }));
  const r = await post("/api/tumour/plan", { experiment: S.exp, factors: rows("planFactors"), readouts: rows("planReadouts"), families: fams });
  toast(r.message); S.project = null; S.tov = null; render();
}
function unlockPlan() {
  openModal(`<div class="modal-head"><span class="title">Unlock the analysis plan?</span><span class="spacer"></span><button id="ulCancel">Cancel</button></div>
    <p>Results exist, so the plan they were built with is locked. Unlock it only for a reason that does not come from the results
      (e.g. a decision recorded before the data were seen, or a typo in a level). The unlock, your reason and every change of the plan are written to
      <code>lock_log.csv</code>, and the next run prints the history (old → new lock) in the report and the About sheet of results.xlsx.</p>
    <label>Reason<br><textarea id="ulReason" rows="3" style="width:100%" placeholder="why the plan changes"></textarea></label>
    <div class="row" style="margin-top:8px"><button class="primary" id="ulYes">Unlock and log</button></div>`);
  $("ulCancel").onclick = () => closeModal();
  $("ulYes").onclick = async () => {
    const r = await post("/api/tumour/unlock", { experiment: S.exp, reason: $("ulReason").value });
    closeModal(); toast(r.message); S.project = null; render();
  };
  $("ulReason").focus();
}
async function relockPlan() {
  const r = await post("/api/tumour/relock", { experiment: S.exp });
  toast(r.message); S.project = null; render();
}

// ---------------------------------------------------------------- dispatch: tumour projects use the views above

for (const sec of ["home", "layout", "review", "corrections", "results", "stats"]) {
  const base = RENDER[sec];
  RENDER[sec] = (b, e, tok) => (isTumour(e) ? T[sec] : base)(b, e, tok);
}
{
  const base = RENDER.settings;
  RENDER.settings = async (b, e, tok) => { await base(b, e, tok); if (isTumour(e)) await tumourPlanCard(b); };
}
