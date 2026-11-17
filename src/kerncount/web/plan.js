"use strict";
// Planner (sampling/gN/planner.json through planner.load / evaluate / curve; chart drawn here in SVG, no
// libraries) and plate maps (kerncount.platemap, in-process).

// ---------------------------------------------------------------- Planner

RENDER.planner = async function (b, e, tok) {
  const P = S.params;
  if (!S.plan || S.plan.gen !== (P.gen || S.plan.gen) || S.plan.exp !== S.exp || (P.cl && S.plan.cell_line !== P.cl)) {
    S.plan = await get("/api/planner", { exp: S.exp, gen: P.gen || "", cl: P.cl || "" });
    S.plan.exp = S.exp;
  }
  if (tok !== S.renderTok) return;
  const pl = S.plan, gens = e.generations;
  const runCard = `<div class="card"><h3>Sampling analysis</h3><p class="small muted" style="margin-top:0">Splits the variation into field / well / n, estimates power for other designs and subsamples your data. The planner below is built from it.</p>
    <div class="row"><label>Generation <select id="sGen">${gens.map(g => `<option ${g === gens[gens.length - 1] ? "selected" : ""}>${esc(g)}</option>`).join("")}</select></label>
      <label>Subsamples per design <input id="sIt" value="1000" size="6"></label>
      <button class="primary" ${S.running || !gens.length ? "disabled" : ""} ${on("click", function (event) { startJob('sampling',{generation:$('sGen').value,iterations:$('sIt').value}) })}>Run sampling analysis</button></div></div>`;
  if (!pl.gen) { b.innerHTML = `<div class="notice info">No planner yet: run the sampling analysis on a generation first.</div>${runCard}`; return; }
  const files = pl.files || [], base = pl.base || `sampling/${pl.gen}/`;
  const plots = files.filter(f => f.path.endsWith(".png")).map(f => ({ ...f, kind: "model" }));
  const tables = files.filter(f => /\.(csv|txt)$/.test(f.path));
  const statics = `<h2>Recommendations for the next experiment</h2>
    ${pl.recommendations ? `<pre class="text">${esc(pl.recommendations)}</pre>` : `<div class="muted">No recommendations in REPORT.txt.</div>`}
    ${plots.length ? `<h2>Planner plots</h2><div class="plots">${plots.map(f => plotCard(f, base)).join("")}</div>` : ""}
    <h3>Tables and report</h3><div class="tbllist">${tables.map(f => `<button class="btn" title="${esc(f.description)}" ${on("click", function (event) { previewFile((base + f.path)) })}>${esc(f.title)}</button>`).join("")}</div>
    ${runCard}`;
  const head = `<div class="selrow"><label class="f">Based on sampling of ${pl.generations.length > 1 ? `<div class="pills">${pl.generations.map(g => `<button class="pill ${g === pl.gen ? "active" : ""}" ${on("click", function (event) { S.plan=null;setP({gen:String(g)}) })}>${g}</button>`).join("")}</div>` : `<span class="tag">${esc(pl.gen)}</span>`}</label>
      ${(pl.cell_lines || []).length > 1 ? `<label class="f">Cell line <div class="pills">${pl.cell_lines.map(c => `<button class="pill ${c.cell_line === pl.cell_line ? "active" : ""}" ${on("click", () => { S.plan = null; setP({ cl: c.cell_line }); })}>${esc(c.cell_line)} <span class="muted">(${c.n_count} n)</span></button>`).join("")}</div></label>` : (pl.cell_line ? `<span class="tag">${esc(pl.cell_line)}</span>` : "")}
      <span class="small muted">${esc(pl.metric || "")} · ${pl.design.n_count} n · ${pl.design.conditions} conditions · ${pl.design.wells} wells × ${pl.design.fields_per_well} fields now${(pl.excluded_ns || []).length ? ` · left out (acceptance enforced): ${esc(pl.excluded_ns.join(", "))}` : ""}</span></div>`;
  if (!pl.ok) {
    b.innerHTML = head + `<div class="notice n1"><b>Not enough data to plan yet${pl.design.n_count < 2 ? ": needs ≥ 2 n" : ""}.</b> The planner needs at least two n's that contain the reference and the treatments, because the n-to-n variation can only be estimated from several independent experiments. (${esc(pl.error || "")})</div>` + statics;
    return;
  }
  const d = pl.design, o = pl.options;
  const refs = d.references, ref = refs.includes(P.ref) ? P.ref : d.main_reference;
  const effs = pl.effects.filter(x => x.reference === ref && x.n_pairs >= 2);
  const v = {
    n: +(P.n || Math.max(d.n_count + 1, 3)), wells: +(P.wells || d.wells), fields: +(P.fields || d.fields_per_well),
    effect: P.effect || String(o.effects.includes(0.7) ? 0.7 : o.effects[0]), rw: P.rw || "", alpha: +(P.alpha || o.alpha),
    power: +(P.power || o.target_power), vary: ["n", "wells", "fields"].includes(P.vary) ? P.vary : "n", show: P.show === "ci" ? "ci" : "power",
    positions: P.pos || "", adjust: P.adj === "0" ? false : P.adj === "1" ? true : o.adjust_for_family !== false,
  };
  const isTreat = v.effect.startsWith("t:");
  const slider = (key, label, min, max, step, val, hint, fmt = x => x) => `<div class="ctl"><div class="lab"><span>${label}</span><span class="val" id="v_${key}">${fmt(val)}</span></div>
    <input type="range" min="${min}" max="${max}" step="${step}" value="${val}" data-k="${key}">${hint ? `<div class="hint">${hint}</div>` : ""}</div>`;
  b.innerHTML = head + `
    <div class="plan">
      <div class="card" id="planCtl">
        <h3>Next experiment</h3>
        ${slider("n", "Number of n (independent experiments)", 2, Math.max(o.max_n, 20), 1, v.n, `now: ${d.n_count}`)}
        ${slider("wells", "Wells per condition", 1, 8, 1, v.wells, `now: ${d.wells}`)}
        ${slider("fields", "Fields (images) per well", 1, Math.max(3 * d.fields_per_well, 24), 1, v.fields, `now: ${d.fields_per_well}`)}
        <div class="ctl"><div class="lab"><span>Effect to detect</span></div>
          <select data-k="effect" style="width:100%">
            <optgroup label="Planning effect (ratio vs reference)">${[0.5, 0.6, 0.7, 0.8, 0.9, ...o.effects].filter((x, i, a) => a.indexOf(x) === i).sort().map(x => `<option value="${x}" ${!isTreat && +v.effect === x ? "selected" : ""}>ratio ${x} (${Math.round((1 - x) * 100)}% ${x < 1 ? "fewer" : "more"})</option>`).join("")}</optgroup>
            <optgroup label="Observed effect of a treatment">${effs.map(x => `<option value="t:${esc(x.treatment)}" ${v.effect === "t:" + x.treatment ? "selected" : ""}>${esc(x.treatment)} (observed ${(+x.ratio).toFixed(2)})</option>`).join("")}</optgroup></select>
          ${isTreat ? `<div class="hint">Observed in few n's: real effects are usually smaller.</div>` : ""}</div>
        <div class="ctl"><div class="lab"><span>Compared with</span></div><select data-k="ref" style="width:100%">${refs.map(r => `<option ${r === ref ? "selected" : ""}>${esc(r)}</option>`).join("")}</select></div>
        <div class="ctl"><div class="lab"><span>Wells per reference</span></div>
          <select data-k="rw" style="width:100%"><option value="">same as treatments</option><option value="auto" ${v.rw === "auto" ? "selected" : ""}>auto (Dunnett √k)</option>${[1, 2, 3, 4, 5, 6, 8].map(x => `<option value="${x}" ${v.rw === String(x) ? "selected" : ""}>${x}</option>`).join("")}</select></div>
        <div class="row"><div class="ctl" style="flex:1"><div class="lab"><span>alpha</span></div><select data-k="alpha">${[0.01, 0.05, 0.1].map(x => `<option ${x === v.alpha ? "selected" : ""}>${x}</option>`).join("")}</select></div>
          <div class="ctl" style="flex:1"><div class="lab"><span>target power</span></div><select data-k="power">${[0.7, 0.8, 0.9, 0.95].map(x => `<option ${x === v.power ? "selected" : ""}>${x}</option>`).join("")}</select></div></div>
        <div class="ctl"><label class="small" title="On: power and smallest detectable effect for a star on the headline plot, i.e. at the multiple-comparison level the statistics use (single-step multivariate t over the ${pl.family_size || ""} comparisons vs the reference). Off: one unadjusted test."><input type="checkbox" data-k="adj" ${v.adjust ? "checked" : ""}> Power for a <b>star</b> (adjusted for the ${pl.family_size || ""} comparisons)</label></div>
        <div class="ctl"><div class="lab"><span>Field positions</span></div><select data-k="pos" style="width:100%"><option value="">as now (${esc(d.positions_default)})</option><option value="standardised" ${v.positions === "standardised" ? "selected" : ""}>standardised (= fixed: same spots in every well)</option><option value="random" ${v.positions === "random" ? "selected" : ""}>random</option></select></div>
      </div>
      <div>
        <div class="metrics4" id="planOut"><div class="mcard"><div class="k">Power</div><div class="v">…</div></div></div>
        <div class="card"><div class="row" style="margin-bottom:6px"><b id="chartTitle">Chart</b><span class="spacer"></span>
          <span class="small muted">vary</span><div class="pills">${["n", "wells", "fields"].map(x => `<button class="pill ${v.vary === x ? "active" : ""}" ${on("click", function (event) { setP({vary:String(x)},false);planUpdate() })}>${x === "n" ? "n" : x}</button>`).join("")}</div>
          <span class="small muted">show</span><div class="pills"><button class="pill ${v.show === "power" ? "active" : ""}" ${on("click", function (event) { setP({show:'power'},false);planUpdate() })}>power</button><button class="pill ${v.show === "ci" ? "active" : ""}" ${on("click", function (event) { setP({show:'ci'},false);planUpdate() })}>CI width</button></div></div>
          <div id="planChart" style="min-height:300px"></div>
          <div class="small muted">Shaded band: pessimistic to optimistic (10–90% of the uncertainty in the variance estimates); line: expected. With few n's plan with the pessimistic end.</div></div>
        <div class="card" id="planNotes"></div>
      </div>
    </div>` + statics;
  $("planCtl").addEventListener("input", ev => {
    const el = ev.target, k = el.dataset.k; if (!k) return;
    if (el.type === "range") $("v_" + k).textContent = el.value;
    setP({ [k]: el.type === "checkbox" ? (el.checked ? "1" : "0") : el.value }, false);
    clearTimeout(S.planTimer); S.planTimer = setTimeout(planUpdate, el.type === "range" ? 40 : 0);
  });
  $("planCtl").addEventListener("change", ev => { if (ev.target.tagName === "SELECT" && ["effect", "ref"].includes(ev.target.dataset.k)) render(); });
  planUpdate();
};

async function planUpdate() {
  const P = S.params, pl = S.plan, d = pl.design, o = pl.options;
  const body = { experiment: S.exp, gen: pl.gen, n: P.n || Math.max(d.n_count + 1, 3), wells: P.wells || d.wells, fields: P.fields || d.fields_per_well,
    effect: (P.effect || String(o.effects.includes(0.7) ? 0.7 : o.effects[0])).replace(/^t:/, ""), reference_wells: P.rw || "",
    alpha: P.alpha || o.alpha, target_power: P.power || o.target_power, reference: P.ref || "", positions: P.pos || "", vary: P.vary || "n",
    cell_line: pl.cell_line || "", adjust_for_family: P.adj === "0" ? false : P.adj === "1" ? true : o.adjust_for_family !== false };
  const seq = (S.planSeq = (S.planSeq || 0) + 1);
  let r;
  try { r = await post("/api/planner/evaluate", body, true); } catch (e) { return; }
  if (seq !== S.planSeq || !$("planOut")) return;
  const x = r.result, tp = +body.target_power;
  const band = (o2, f) => `${f(o2.low)} – ${f(o2.high)}`;
  const pw = v => Math.round(100 * v) + "%", ci = v => "±" + Math.round(v) + "%", rt = v => v.toFixed(2);
  $("planOut").innerHTML = `
    <div class="mcard ${x.power.mid >= tp ? "good" : "bad"}"><div class="k">Power</div><div class="v">${pw(x.power.mid)}</div><div class="b">band ${band(x.power, pw)} · target ${pw(tp)}</div></div>
    <div class="mcard"><div class="k">95% CI of the ratio</div><div class="v">${ci(x.ci_halfwidth_pct.mid)}</div><div class="b">band ${band({ low: x.ci_halfwidth_pct.low, high: x.ci_halfwidth_pct.high }, ci)}</div></div>
    <div class="mcard"><div class="k">Smallest detectable ratio</div><div class="v">${rt(x.mde_ratio.mid)}</div><div class="b">at ${pw(tp)} power · band ${band({ low: x.mde_ratio.low, high: x.mde_ratio.high }, rt)}</div></div>
    <div class="mcard"><div class="k">Effort</div><div class="v">${Math.round(x.cost).toLocaleString()}</div><div class="b">${esc((pl.costs && pl.costs.unit) || "units")} · ${x.wells_per_n} wells per n · ${x.images} images per condition</div></div>`;
  $("planNotes").innerHTML = `<h3>What this means</h3>
    <p style="margin-top:0">With <b>${x.n}</b> n's, <b>${x.wells}</b> well(s) × <b>${x.fields}</b> field(s) per condition and <b>${x.reference_wells}</b> reference well(s), a true ${esc(x.effect_label)} vs ${esc(x.reference)}
      would ${x.adjust_for_family ? `get a <b>star</b> on the headline plot (adjusted p &lt; ${body.alpha}, ${x.family_size} comparisons vs ${esc(x.reference)})` : `be significant in one unadjusted test (p &lt; ${body.alpha}; not the chance of a star)`} in about <b>${pw(x.power.mid)}</b> of experiments (pessimistic ${pw(x.power.low)}).
      ${Math.round(100 * x.measurement_share)}% of the noise in one n's ratio is measurement (wells and fields); the rest is n-to-n biology, which only more n's reduce.</p>
    ${x.notes.length ? `<ul>${x.notes.map(n => `<li>${esc(n)}</li>`).join("")}</ul>` : ""}
    <div class="small muted">Test: ${esc(x.test)} · critical value ${x.t_critical.toFixed(2)} (df ${x.df.toFixed(1)}${x.adjust_for_family ? `; unadjusted ${x.t_critical_unadjusted.toFixed(2)}` : ""}) · treatment class ${esc(x.cls)} · positions ${esc(x.positions)}</div>`;
  drawPlanChart(r.curve, r.vary, P.show === "ci" ? "ci" : "power", +body[r.vary === "n" ? "n" : r.vary], tp);
}

// SVG line + band chart (no libraries)
function drawPlanChart(pts, vary, show, cur, target) {
  const W = 760, H = 300, m = { l: 52, r: 16, t: 12, b: 40 }, iw = W - m.l - m.r, ih = H - m.t - m.b;
  const key = show === "ci" ? "ci_halfwidth_pct" : "power";
  const ys = pts.map(p => p[key]);
  const xmin = pts[0].x, xmax = pts[pts.length - 1].x;
  let ymax = show === "ci" ? Math.min(Math.max(...ys.map(y => y.mid)) * 1.6, Math.max(...ys.map(y => y.high))) : 1;
  let ystep = 0.25;
  if (show === "ci") {  // round tick steps: 10, 20, 25, 50, 100, 200, ... per cent
    ystep = [10, 20, 25, 50, 100, 200, 250, 500, 1000].find(s => ymax / s <= 5) || 2000;
    ymax = Math.max(ystep, Math.ceil(ymax / ystep) * ystep);
  }
  const X = x => m.l + (xmax === xmin ? iw / 2 : (x - xmin) / (xmax - xmin) * iw), Y = y => m.t + ih - Math.min(Math.max(y, 0), ymax) / ymax * ih;
  const fmtY = y => show === "ci" ? "±" + Math.round(y) + "%" : Math.round(y * 100) + "%";
  const yt = Array.from({ length: Math.round(ymax / ystep) + 1 }, (_, i) => i * ystep);
  const step = Math.max(1, Math.ceil((xmax - xmin) / 12));
  const xt = pts.map(p => p.x).filter(x => (x - xmin) % step === 0);
  const area = pts.map(p => `${X(p.x)},${Y(p[key].high)}`).join(" ") + " " + pts.slice().reverse().map(p => `${X(p.x)},${Y(p[key].low)}`).join(" ");
  const line = pts.map(p => `${X(p.x)},${Y(p[key].mid)}`).join(" ");
  const curPt = pts.find(p => p.x === cur);
  const col = show === "ci" ? "var(--k-relative)" : "var(--k-model)";
  const lab = { n: "number of n", wells: "wells per condition", fields: "fields per well" }[vary];
  $("chartTitle").textContent = (show === "ci" ? "95% CI half-width" : "Power") + " vs " + lab;
  $("planChart").innerHTML = `<svg class="chart" viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="${esc($("chartTitle").textContent)}">
    ${yt.map(y => `<line class="gr" x1="${m.l}" x2="${W - m.r}" y1="${Y(y)}" y2="${Y(y)}"/><text x="${m.l - 6}" y="${Y(y) + 4}" text-anchor="end">${fmtY(y)}</text>`).join("")}
    ${xt.map(x => `<text x="${X(x)}" y="${H - m.b + 16}" text-anchor="middle">${x}</text>`).join("")}
    <line class="ax" x1="${m.l}" x2="${W - m.r}" y1="${m.t + ih}" y2="${m.t + ih}"/><line class="ax" x1="${m.l}" x2="${m.l}" y1="${m.t}" y2="${m.t + ih}"/>
    <text x="${m.l + iw / 2}" y="${H - 6}" text-anchor="middle">${lab}</text>
    <polygon points="${area}" fill="${col}" fill-opacity="0.15" stroke="none"/>
    <polyline points="${line}" fill="none" stroke="${col}" stroke-width="2.5"/>
    ${pts.map(p => `<circle cx="${X(p.x)}" cy="${Y(p[key].mid)}" r="3" fill="${col}"><title>${lab} ${p.x}: ${fmtY(p[key].mid)} (${fmtY(p[key].low)} – ${fmtY(p[key].high)})</title></circle>`).join("")}
    ${show === "power" ? `<line x1="${m.l}" x2="${W - m.r}" y1="${Y(target)}" y2="${Y(target)}" stroke="#16a34a" stroke-dasharray="5 4"/><text x="${W - m.r}" y="${Y(target) - 5}" text-anchor="end" style="fill:#16a34a">target ${Math.round(target * 100)}%</text>` : ""}
    ${curPt ? `<line x1="${X(cur)}" x2="${X(cur)}" y1="${m.t}" y2="${m.t + ih}" stroke="#111827" stroke-dasharray="2 3"/><circle cx="${X(cur)}" cy="${Y(curPt[key].mid)}" r="6" fill="#fff" stroke="#111827" stroke-width="2"/>
      <text x="${Math.min(X(cur) + 8, W - 120)}" y="${m.t + 12}" style="font-weight:600;fill:#111827">your design: ${fmtY(curPt[key].mid)}</text>` : ""}
  </svg>`;
}

// ---------------------------------------------------------------- Plate maps

RENDER.plates = async function (b, e, tok) {
  if (!S.plans) S.plans = await get("/api/plans", { exp: S.exp });
  if (tok !== S.renderTok) return;
  const pi = S.plans, df = pi.defaults, P = S.params;
  const sel = pi.plans.find(p => p.name === P.plan) || pi.plans[0];
  b.innerHTML = `
    <p class="muted" style="margin-top:0">Randomised plate maps for the next experiment: one plate per n. Conditions rotate through edge and inner wells, rows and columns across plates, so plate position effects (evaporation at the edge, imaging order) spread over every condition instead of biasing one.</p>
    <div class="cols-side">
      <div class="card"><h3>New plate maps</h3>
        <div class="kv" id="pmForm">
          <label><div class="k">plates (= n's)</div><input id="pm_plates" value="${esc(P.plates || 3)}" size="4"></label>
          <label><div class="k">plate format</div><select id="pm_fmt">${[6, 12, 24, 48, 96].map(x => `<option ${x === +(P.fmt || df.plate_format) ? "selected" : ""}>${x}</option>`).join("")}</select></label>
          <label><div class="k">wells per condition</div><input id="pm_wells" value="${esc(P.wells || 2)}" size="4"></label>
          <label><div class="k">wells per reference</div><select id="pm_rw"><option value="">settings (${esc(df.reference_wells)})</option><option value="auto">auto (Dunnett √k)</option>${[1, 2, 3, 4, 5, 6].map(x => `<option>${x}</option>`).join("")}</select></label>
          <label><div class="k">seed (empty = random)</div><input id="pm_seed" value="${esc(P.seed || "")}" size="10"></label>
          <label><div class="k">negative control</div><input id="pm_nc" value="${esc(df.negative_control || "")}" placeholder="e.g. secondary only" size="14"></label>
          <label><div class="k">negative-control wells</div><select id="pm_ncw"><option>2</option><option>1</option></select></label>
        </div>
        <div style="margin-top:8px"><label class="chk"><input type="checkbox" id="pm_edge" ${df.avoid_edge ? "checked" : ""}> avoid edge wells (inner wells first)</label>
          <label class="chk"><input type="checkbox" id="pm_bal" checked> balance across plates</label></div>
        <p class="small muted">Conditions: ${pi.conditions.map(c => `<span class="tag ${pi.references.includes(c) ? "blue" : ""}">${esc(c)}</span>`).join("") || "none found"}<br>Blue = reference (gets the reference wells). The negative control gets its own wells on every plate and is left out of the statistics.</p>
        <button class="primary" id="pmGo" ${on("click", function (event) { makePlates() })}>Generate plate maps</button> <span class="small muted" id="pmState"></span></div>
      <div class="card"><h3>Earlier plans</h3>
        ${pi.plans.length ? `<table class="grid"><tbody>${pi.plans.map(p => `<tr class="click ${p === sel ? "sel" : ""}" ${on("click", function (event) { setP({plan:(p.name)}) })}><td>${esc(p.name.replace("platemap_", ""))}</td><td class="small muted">${p.plates.length} plate(s)</td></tr>`).join("")}</tbody></table>`
          : `<div class="muted">None yet.</div>`}</div>
    </div>
    ${sel ? `<h2>${esc(sel.name)}</h2>
      <div class="row" style="margin-bottom:10px">${sel.files.filter(f => /\.(pdf|csv)$/.test(f)).map(f => `<a class="btn small" href="${fileUrl(f, 1)}">${esc(baseName(f))}</a>`).join("")}
        <button class="small" ${on("click", function (event) { viewText((sel.path + "/README.txt")) })}>README</button></div>
      <div class="plots">${sel.plates.map(p => `<div class="plot platecard k-design"><div class="ph"><div class="pt">${esc(baseName(p).replace(".png", "").replace("_", " "))}</div></div>
        <img loading="lazy" src="${fileUrl(p)}" ${on("click", function (event) { openImage((p),(baseName(p)),"",(p.replace(/\.png$/, ".pdf"))) })}>
        <div class="pf"><span class="spacer"></span><a class="btn small" href="${fileUrl(p.replace(/\.png$/, ".pdf"), 1)}">PDF (A4)</a><a class="btn small" href="${fileUrl(p.replace(/\.png$/, ".csv"), 1)}">CSV</a></div></div>`).join("")}</div>
      ${sel.edge ? `<h3>Edge-well use per condition (over all plates)</h3><div class="tablewrap" style="max-height:none"><table class="grid"><thead><tr>${sel.edge.columns.map(c => `<th>${esc(c)}</th>`).join("")}</tr></thead>
        <tbody>${sel.edge.rows.map(r => `<tr>${r.map(v => `<td>${esc(fmtCell(v))}</td>`).join("")}</tr>`).join("")}</tbody></table></div>` : ""}
      <details style="margin-top:10px"><summary class="small">README (seed and the command that recreates these maps)</summary><pre class="text">${esc(sel.readme)}</pre></details>` : ""}`;
};
async function makePlates() {
  const v = id => $(id).value.trim();
  $("pmGo").disabled = true; $("pmState").textContent = "drawing the plates…";
  try {
    const r = await post("/api/platemap", { experiment: S.exp, plates: v("pm_plates"), plate_format: v("pm_fmt"), wells: v("pm_wells"),
      reference_wells: v("pm_rw"), seed: v("pm_seed"), negative_control: v("pm_nc"), negative_control_wells: v("pm_ncw"),
      avoid_edge: $("pm_edge").checked, balance: $("pm_bal").checked });
    S.plans = null; setP({ plan: r.name });
  } catch (e) { if ($("pmGo")) { $("pmGo").disabled = false; $("pmState").textContent = ""; } }
}
