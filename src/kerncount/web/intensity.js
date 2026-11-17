"use strict";
// Intensity page (intensity/manifest.json: flat field, cut-offs, sweep, distributions, density, alignment) and the
// BLINDED cut-off tuning view (kerncount.tuning through /api/tuning/*). The tuning view never shows a condition,
// an n of a field or a file name: fields are opaque ids in random order (docs/ui.md, "Blinding").

const ITABS = [["overview", "Cut-offs & flat field"], ["dist", "Distributions & density"], ["checks", "Channel checks"], ["tune", "Tune the cut-off (blinded)"]];

async function loadIntensity(force) {
  if (!S.intens || force) S.intens = await get("/api/intensity", { exp: S.exp });
  return S.intens;
}

RENDER.intensity = async function (b, e, tok) {
  const I = await loadIntensity();
  if (tok !== S.renderTok) return;
  const im = I.manifest;
  if (!im.files) { b.innerHTML = `<div class="notice info">No intensity results: this experiment has no marker channel measured per nucleus (or it has not been run yet).</div>`; return; }
  const tab = ITABS.find(t => t[0] === S.params.tab) ? S.params.tab : "overview";
  b.innerHTML = `<div class="tabs">${ITABS.map(([k, l]) => `<button class="tab ${k === tab ? "active" : ""}" onclick="setP({tab:'${k}'})">${l}</button>`).join("")}</div><div id="itab"></div>`;
  const body = $("itab");
  if (tab === "tune") return renderTune(body, tok);
  if (!S.ov) S.ov = await get("/api/overview", { exp: S.exp }).catch(() => null);
  if (tok !== S.renderTok) return;
  const F = im.files, I0 = "intensity/";
  const pngs = k => F.filter(f => f.kind === k && f.path.endsWith(".png"));
  if (tab === "overview") {
    const lockFor = (c) => I.locks.find(l => l.channel.toLowerCase() === c.channel.toLowerCase() && (l.measure || "mean") === c.measure && (!l.n || l.n === c.n));
    body.innerHTML = `
      <p class="muted" style="margin-top:0">Which nuclei count as <b>${esc(im.channels.join(", "))}-positive</b>, and how uneven illumination was corrected. Cut-off precedence:
        negative-control wells → locked manual value → automatic detection limit (3 × background noise).</p>
      ${!I.negative_control ? `<div class="notice">No negative-control wells: “positive” means <b>detectably stained above background</b>, not “above non-specific staining”.
        ${I.locks.length ? "" : `Set a cut-off by eye in the <a class="link" onclick="setP({tab:'tune'})">blinded tuning view</a> and lock it before looking at the statistics.`}</div>` : ""}
      <div class="card"><h3>Cut-off per n</h3>
        <table class="grid"><thead><tr><th>n</th><th>channel</th><th>measure</th><th>rule used</th><th class="num">cut-off</th><th class="num">% positive</th><th class="num">without flat field</th><th class="num">nuclei</th><th>detail</th></tr></thead>
        <tbody>${im.cutoffs.map(c => { const lk = lockFor(c); return `<tr><td>${nChip(c.n)}</td><td>${esc(c.channel)}</td><td>${c.measure === "p90" ? "90th percentile" : "mean"}</td>
          <td><span class="tag ${c.rule === "lod" ? "medium" : c.rule === "manual" ? "purple" : "ok"}">${esc(c.rule_label)}</span>${lk && c.rule !== "manual" ? `<div class="small" style="color:var(--k-model)">locked ${(+lk.cutoff).toFixed(1)}: applied at the next run</div>` : ""}</td>
          <td class="num">${(+c.value).toFixed(1)}</td><td class="num"><b>${pct(c.frac_positive)}</b></td><td class="num">${pct(c.frac_positive_without_flatfield)}</td>
          <td class="num">${esc(c.nuclei)}</td><td class="small muted">${esc(c.detail)}</td></tr>`; }).join("")}</tbody></table>
        <div class="row" style="margin-top:8px"><button class="small" onclick="viewTable('intensity/cutoffs.csv')">cutoffs.csv</button>
          <button class="small" onclick="viewTable('intensity/positive_by_condition.csv')">% positive per condition</button>
          <button class="small primary" onclick="setP({tab:'tune'})">Tune the cut-off (blinded) →</button></div>
        ${lockHistory(I)}</div>
      <h2>Flat-field correction</h2>
      <p class="small muted" style="margin-top:0">The microscope lights the centre of the field more than the edges. kerncount estimates that pattern per n and channel from all fields and divides the background-subtracted intensities by it; counts are never corrected.</p>
      <div class="cols2">${im.flatfield.map(ff => { const f = F.find(x => x.kind === "flatfield" && x.n === ff.n && x.channel === ff.channel);
        return `<div class="card"><div class="row"><b>${nChip(ff.n)} ${esc(ff.channel)}</b><span class="tag">${esc(ff.source)}${ff.images ? ` from ${ff.images} fields` : ""}</span></div>
          <div class="kv" style="margin:8px 0"><div><div class="k">edges vs centre</div><div class="v">${ff.edge_dim_pct != null ? (+ff.edge_dim_pct).toFixed(0) + "% dimmer" : "–"}</div></div>
            <div><div class="k">darkest side</div><div class="v">${esc(ff.darkest_side || "–")}${ff.darkest_side_dim_pct != null ? ` (${(+ff.darkest_side_dim_pct).toFixed(0)}%)` : ""}</div></div>
            <div><div class="k">range</div><div class="v">${ff.range_pct != null ? (+ff.range_pct).toFixed(0) + "%" : "–"}</div></div>
            <div><div class="k">check from nuclei: edge dimming</div><div class="v">${ff.nuclei_edge_dim_pct_before != null ? `${(+ff.nuclei_edge_dim_pct_before).toFixed(0)}% → ${(+ff.nuclei_edge_dim_pct_after).toFixed(0)}%` : "–"}</div></div>
            <div><div class="k">% positive with / without</div><div class="v">${pct(ff.frac_positive_with_flatfield)} / ${pct(ff.frac_positive_without_flatfield)}</div></div></div>
          ${f ? plotCard(f, I0) : ""}</div>`; }).join("") || `<div class="muted">No flat-field maps.</div>`}</div>
      <div class="row"><button class="small" onclick="viewTable('intensity/flatfield.csv')">flatfield.csv</button></div>
      <h2>Cut-off sweep</h2>
      <p class="small muted" style="margin-top:0">How % positive per condition changes if the cut-off were ¼ to 4 × the one used. Curves that cross mean the ranking of conditions depends on the cut-off.</p>
      <div class="plots">${pngs("sweep").map(f => plotCard({ ...f, kind: "robustness" }, I0)).join("") || `<div class="muted">No sweep plots.</div>`}</div>`;
  } else if (tab === "dist") {
    body.innerHTML = `
      <h3 style="margin-top:0">Per-nucleus intensity</h3>
      <div class="plots">${pngs("distribution").map(f => plotCard({ ...f, kind: "absolute" }, I0)).join("") || `<div class="muted">No distribution plots.</div>`}</div>
      <h3>% positive vs density</h3>
      <div class="plots">${pngs("density").map(f => plotCard({ ...f, kind: "design" }, I0)).join("") || `<div class="muted">No density plots.</div>`}</div>
      <h3>Tables</h3><div class="row"><button onclick="viewTable('intensity/positive_by_condition.csv')">% positive per condition</button><button onclick="viewTable('intensity/positivity_sweep.csv')">Cut-off sweep table</button></div>`;
  } else {
    const al = await get("/api/table", { exp: S.exp, path: "intensity/alignment.csv" });
    if (tok !== S.renderTok) return;
    const ix = Object.fromEntries(al.columns.map((c, i) => [c, i]));
    const bad = r => /mis-paired/.test(r[ix.flags] || "") || +r[ix.corr] < 0.2;
    const flagged = al.rows.filter(r => (r[ix.flags] || "").trim());
    body.innerHTML = `
      <p class="muted" style="margin-top:0">Each field is a pair of files (marker, then DAPI). kerncount measures how far the two images are shifted and whether they show the same cells.
        A field whose files do not match keeps its count but gets no marker readouts.</p>
      <div class="plots">${pngs("alignment").map(f => plotCard({ ...f, kind: "diagnostic" }, I0)).join("")}</div>
      <h3>Flagged fields (${flagged.length} of ${al.rows.length})</h3>
      ${flagged.length ? `<table class="grid"><thead><tr><th>n</th><th>condition</th><th class="num">field</th><th>files</th><th class="num">shift (px)</th><th class="num">correlation</th><th>flags</th></tr></thead>
        <tbody>${flagged.map(r => `<tr class="${bad(r) ? "bad" : "warn"}"><td>${nChip(r[ix.n])}</td><td>${esc(r[ix.treatment])}</td><td class="num">${esc(r[ix.image])}</td><td>${esc(r[ix.files])}</td>
          <td class="num">${esc(r[ix.shift_px])}</td><td class="num">${esc(r[ix.corr])}</td><td class="small">${esc(r[ix.flags])}</td></tr>`).join("")}</tbody></table>`
        : `<div class="notice ok">No field is flagged.</div>`}
      <div class="row" style="margin-top:8px"><button class="small" onclick="viewTable('intensity/alignment.csv')">Full alignment table</button></div>`;
  }
};

function lockHistory(I) {
  const H = I.history || [];
  return `<h3>Locked cut-offs</h3>
    ${I.locks.length ? I.locks.map(l => `<div class="lockbox" style="margin-bottom:6px"><b>${esc(l.channel)} ${l.measure === "p90" ? "(90th percentile)" : "(mean)"} = ${(+l.cutoff).toFixed(1)}</b>
        · ${l.n ? "for " + esc(l.n) : "for every n"} · locked ${esc((l.locked_at || "").replace("T", " "))} on ${l.fields} field(s)</div>`).join("")
      : `<div class="muted small">No locked cut-off: the automatic rule (or negative-control wells) decides.</div>`}
    ${I.manual_cutoff !== "" && I.manual_cutoff != null ? `<div class="small">settings.toml also sets <code>manual_cutoff = ${esc(JSON.stringify(I.manual_cutoff))}</code>, which takes precedence over a lock.</div>` : ""}
    ${H.length ? `<details><summary class="small">History (${H.length})</summary><table class="grid small"><tbody>${H.slice().reverse().map(h => `<tr><td>${esc((h.at || "").replace("T", " "))}</td>
      <td>${h.action === "unlock" ? "unlocked" : "locked"} ${esc(h.channel)} ${h.measure === "p90" ? "(p90)" : ""}</td><td class="num">${h.cutoff != null ? (+h.cutoff).toFixed(1) : ""}</td>
      <td>${h.n ? esc(h.n) : "every n"}</td><td class="muted">${h.previous != null ? "was " + (+h.previous).toFixed(1) : ""}${h.fields ? ` · ${h.fields} fields` : ""}</td></tr>`).join("")}</tbody></table></details>` : ""}`;
}

// ---------------------------------------------------------------- blinded tuning

// slider position (0..1000) <-> cut-off on a log scale between the 1st and 99th percentile of the signal
function sliderToCut(pos, r) { const lo = Math.log(Math.max(r.min, 1e-3) / 2), hi = Math.log(r.max * 1.5); return Math.exp(lo + (hi - lo) * pos / 1000); }
function cutToSlider(c, r) { const lo = Math.log(Math.max(r.min, 1e-3) / 2), hi = Math.log(r.max * 1.5); return Math.round(1000 * (Math.log(Math.max(c, 1e-3)) - lo) / (hi - lo)); }

async function renderTune(body, tok) {
  const P = S.params;
  if (!S.tune || S.tune.exp !== S.exp) {
    body.innerHTML = `<div class="muted">Loading the per-nucleus measurements… (a few seconds the first time)</div>`;
    const info = await get("/api/tuning/info", { exp: S.exp });
    S.tune = { exp: S.exp, info, sample: null, idx: 0, ranges: {} };
  }
  if (tok !== S.renderTok) return;
  const T = S.tune, info = T.info;
  const ch = info.channels.includes(P.ch) ? P.ch : info.channels[0];
  const measure = P.measure === "p90" ? "p90" : "mean";
  const scope = info.scopes.includes(P.scope) ? P.scope : "";
  const key = ch + "|" + measure;
  if (!T.ranges[key]) T.ranges[key] = await get("/api/tuning/range", { exp: S.exp, channel: ch, measure });
  const k = +(P.k || 12), seed = P.seed || "";
  if (!T.sample || T.sample.channel !== ch || T.sampleKey !== `${scope}|${k}|${seed}`) {
    T.sample = await post("/api/tuning/sample", { experiment: S.exp, channel: ch, n: scope, k, seed }, true);
    T.sampleKey = `${scope}|${k}|${seed}`;
    T.idx = 0;
  }
  if (tok !== S.renderTok) return;
  const R = T.ranges[key], fields = T.sample.fields;
  T.idx = Math.max(0, Math.min(fields.length - 1, P.field != null ? +P.field : T.idx));
  T.cut = P.cutoff != null && !isNaN(+P.cutoff) ? +P.cutoff : (T.cut && T.cutKey === key ? T.cut : R.lod_median);
  T.cutKey = key; T.ch = ch; T.measure = measure; T.scope = scope;
  const lockNow = info.locks.find(l => l.channel === ch && (l.measure || "mean") === measure && (l.n || "") === scope);
  body.innerHTML = `
    <div class="notice info"><b>Blinded.</b> Fields appear in random order under opaque ids; their condition, n and file names are hidden here.
      Set the cut-off by eye so that green outlines mark the nuclei you would call ${esc(ch)}-positive, step through several fields, then <b>lock</b> it —
      before looking at any statistics. A run applies the lock.</div>
    <div class="row" style="margin-bottom:10px">
      ${info.channels.length > 1 ? `<label>Channel <select onchange="setP({ch:this.value,cutoff:''})">${info.channels.map(c => `<option ${c === ch ? "selected" : ""}>${esc(c)}</option>`).join("")}</select></label>` : `<span class="tag purple">${esc(ch)}</span>`}
      <div class="pills"><button class="pill ${measure === "mean" ? "active" : ""}" onclick="setP({measure:'mean',cutoff:''})">mean</button><button class="pill ${measure === "p90" ? "active" : ""}" onclick="setP({measure:'p90',cutoff:''})">90th percentile</button></div>
      <label class="small">Lock for <select onchange="setP({scope:this.value,field:''})"><option value="">every n</option>${info.scopes.map(n => `<option value="${esc(n)}" ${n === scope ? "selected" : ""}>${esc(n)} only</option>`).join("")}</select></label>
      <span class="spacer"></span>
      <label class="small">fields <input size="3" value="${k}" onchange="setP({k:this.value,field:''})"></label>
      <button onclick="setP({seed:String(Math.floor(Math.random()*1e9)),field:''})">New random sample</button>
    </div>
    <div class="tune">
      <div class="stage ${T.zoom ? "zoom" : ""}">
        <div class="imgwrap"><img id="tuneImg" alt="marker image with nucleus outlines" ondblclick="openTuneFull()"></div>
        <span class="busy" id="tuneBusy"></span>
        <div class="fieldbar"><button class="small" onclick="tuneStep(-1)" ${T.idx ? "" : "disabled"}>‹ Prev</button>
          <span>Field <b>${T.idx + 1}</b> of ${fields.length} · id <code>${esc(fields[T.idx].field_id)}</code></span>
          <button class="small" onclick="tuneStep(1)" ${T.idx < fields.length - 1 ? "" : "disabled"}>Next ›</button>
          <button class="small" onclick="S.tune.zoom=!S.tune.zoom;document.querySelector('.tune .stage').classList.toggle('zoom')">zoom 1:1</button><span class="spacer"></span>
          <span class="legend"><span><i style="border-color:#28dc3c"></i>positive</span><span><i style="border-color:#969696"></i>negative</span><span><i style="border-color:#ff9600"></i>saturated</span></span></div>
      </div>
      <div>
        <div class="card"><h3>Cut-off (${measure === "p90" ? "90th percentile" : "mean"} signal)</h3>
          <input type="range" class="slider" id="cutSlider" min="0" max="1000" value="${cutToSlider(T.cut, R)}">
          <div class="row" style="justify-content:space-between"><span class="small muted">${R.min.toFixed(0)}</span>
            <span><input id="cutIn" size="8" value="${T.cut.toFixed(1)}"> <span class="small muted">corrected units</span></span><span class="small muted">${R.max.toFixed(0)}</span></div>
          <div class="small muted" style="margin-top:4px">Automatic rule (median): ${R.lod_median.toFixed(1)} <button class="small" onclick="setTuneCut(${R.lod_median})">use</button></div>
          <div class="kv" style="margin-top:12px">
            <div><div class="k">% positive, these ${fields.length} fields</div><div class="big" id="pvSample">…</div></div>
            <div><div class="k">% positive, ${scope ? "this n" : "all fields"}</div><div class="big" id="pvAll">…</div></div>
          </div>
          <div class="small muted" id="pvLod" style="margin-top:4px"></div></div>
        <div class="lockbox"><b>Lock this cut-off</b>
          <p class="small" style="margin:4px 0 8px">Fixes ${esc(ch)} ${measure === "p90" ? "(90th percentile)" : "(mean signal)"} for ${scope ? esc(scope) : "every n"}.
            ${measure === "p90" ? "" : `<b>Only the mean-based % positive:</b> the 90th-percentile readout keeps the automatic rule unless you lock it too (switch to “90th percentile”).`} ${lockNow ? `Currently locked: <b>${(+lockNow.cutoff).toFixed(1)}</b> (${esc((lockNow.locked_at || "").replace("T", " "))}).` : "Nothing is locked yet."}</p>
          <div class="row"><button class="primary" onclick="tuneLock()" ${S.running ? "disabled" : ""}>Lock ${esc(T.cut.toFixed(1))}…</button>
            ${lockNow ? `<button class="small" onclick="tuneUnlock()" ${S.running ? "disabled" : ""}>Unlock</button>` : ""}</div></div>
        <div class="card" style="margin-top:12px">${lockHistory({ ...info, manual_cutoff: "" })}</div>
      </div>
    </div>`;
  const sl = $("cutSlider");
  sl.addEventListener("input", () => setTuneCut(sliderToCut(+sl.value, R), true));
  sl.addEventListener("change", () => setP({ cutoff: T.cut.toFixed(2), field: T.idx }, false));
  $("cutIn").addEventListener("change", ev => { const v = parseFloat(ev.target.value.replace(",", ".")); if (!isNaN(v)) setTuneCut(v); });
  tuneDraw();
}
let tuneTimer = null, tuneSeq = 0;
function setTuneCut(v, fromSlider) {
  const T = S.tune, R = T.ranges[T.ch + "|" + T.measure];
  T.cut = v;
  if (!fromSlider && $("cutSlider")) $("cutSlider").value = cutToSlider(v, R);
  if ($("cutIn")) $("cutIn").value = v.toFixed(1);
  const lb = document.querySelector(".lockbox .primary"); if (lb) lb.textContent = `Lock ${v.toFixed(1)}…`;
  clearTimeout(tuneTimer);
  tuneTimer = setTimeout(() => { tuneDraw(); if (!fromSlider) setP({ cutoff: v.toFixed(2) }, false); }, 70);
}
function tuneImgUrl() {
  const T = S.tune, f = T.sample.fields[T.idx];
  return "/api/tuning/render?" + qs({ exp: S.exp, field: f.field_id, channel: T.ch, cutoff: T.cut.toFixed(3), measure: T.measure, t: TOKEN });
}
async function tuneDraw() {
  const T = S.tune, seq = ++tuneSeq, img = $("tuneImg");
  if (!img) return;
  $("tuneBusy").textContent = "drawing…";
  const pre = new Image();
  pre.onload = () => { if (seq === tuneSeq && $("tuneImg")) { $("tuneImg").src = pre.src; $("tuneBusy").textContent = ""; } };
  pre.onerror = () => { if ($("tuneBusy")) $("tuneBusy").textContent = "could not draw this field"; };
  pre.src = tuneImgUrl();
  try {
    const p = await get("/api/tuning/preview", { exp: S.exp, cutoff: T.cut, channel: T.ch, measure: T.measure, n: T.scope, fields: T.sample.fields.map(f => f.field_id).join(",") });
    if (seq !== tuneSeq || !$("pvAll")) return;
    $("pvAll").textContent = pct(p.frac_positive);
    $("pvSample").textContent = p.sample ? pct(p.sample.frac_positive) : "–";
    $("pvLod").textContent = `${p.positive} of ${p.nuclei} nuclei positive. The automatic rule would call ${pct(p.lod_frac_positive)} positive.`;
  } catch (e) { toast(e.message, true); }
}
function tuneStep(d) { const T = S.tune; T.idx = Math.max(0, Math.min(T.sample.fields.length - 1, T.idx + d)); setP({ field: T.idx, cutoff: T.cut.toFixed(2) }); }
function openTuneFull() {
  openModal(`<div class="modal-head"><span class="title">Field ${S.tune.idx + 1} · cut-off ${S.tune.cut.toFixed(1)}</span><span class="spacer"></span><button onclick="closeModal()">Close</button></div>
    <img class="full" src="${tuneImgUrl()}" style="max-height:84vh">`);
}
function tuneLock() {
  const T = S.tune, v = T.cut;
  confirmBox(`<h3 style="margin-top:0">Lock ${esc(T.ch)} cut-off = ${v.toFixed(1)}?</h3>
    <p>This fixes the ${T.measure === "p90" ? "90th-percentile" : "mean"} cut-off for <b>${T.scope ? esc(T.scope) : "every n"}</b>, tuned on ${T.sample.fields.length} blinded fields.</p>
    <p><b>Lock it before looking at any statistics</b>: a cut-off chosen after seeing which condition “wins” is no longer an independent choice.</p>
    <p class="small muted">The value is stored in <code>cutoffs.json</code> with the fields it was tuned on. It takes effect at the next run (detections are reused, so the run is quick); the generations are rebuilt automatically.</p>`,
    async () => {
      const r = await post("/api/tuning/lock", { experiment: S.exp, channel: T.ch, cutoff: v, measure: T.measure, n: T.scope, field_ids: T.sample.fields.map(f => f.field_id) });
      T.info = await get("/api/tuning/info", { exp: S.exp }); S.intens = null; S.ov = null;
      render();
      confirmBox(`Locked ${esc(r.lock.channel)} = ${(+r.lock.cutoff).toFixed(1)}. <b>Run the analysis now</b> to apply it?`, () => startJob("run"), "Start the run", "Later");
    }, "Lock the cut-off");
}
function tuneUnlock() {
  const T = S.tune;
  confirmBox(`Remove the locked ${esc(T.ch)} cut-off for ${T.scope ? esc(T.scope) : "every n"}? The next run then uses the automatic rule (or negative-control wells).`, async () => {
    await post("/api/tuning/unlock", { experiment: S.exp, channel: T.ch, measure: T.measure, n: T.scope });
    T.info = await get("/api/tuning/info", { exp: S.exp }); S.intens = null; S.ov = null; render();
  }, "Unlock");
}
