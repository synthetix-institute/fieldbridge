/* The page around the instrument: definitions, stepped sequences, the co-discovery plots and the sources. */
(function () {
  'use strict';
  const S = window.FIELDBRIDGE_SITE;
  const $ = id => document.getElementById(id);
  if (!S) {
    document.addEventListener('DOMContentLoaded', () => { $('instrument').innerHTML = '<p class="note">The data of the page did not load. Build it with <code>python3 -B -m fieldbridge demo</code> and open build/site/index.html.</p>'; });
    return;
  }
  const V = window.FieldBridgeViews;
  const esc = s => String(s).replace(/[&<>"']/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[c]));
  const strip = s => String(s).replace(/<[^>]+>/g, '');
  const repo = p => (S.repo || 'https://github.com/synthetix-institute/fieldbridge') + '/blob/main/' + p;
  const I = () => window.FieldBridgeInstrument;
  const seq = {current: null, index: -1};

  // ---------------------------------------------------------------------------------------- the expression
  function definitions() {
    $('definitions').innerHTML = Object.entries(S.slots).map(([k, s]) =>
      `<div data-slot="${k}"><dt><span class="s">${s.symbol}</span>${esc(s.name)}</dt><dd>${esc(s.definition)}</dd></div>`).join('');
    document.querySelectorAll('.formula .sym').forEach(b => {
      const d = S.slots[b.dataset.slot];
      if (d) b.title = `${d.name}: ${d.definition}`;
      b.addEventListener('mouseenter', () => I().light(b.dataset.slot));
      b.addEventListener('mouseleave', () => I().light(null));
      b.addEventListener('click', () => {
        const row = $('slot-' + b.dataset.slot);
        if (!row) return;
        row.scrollIntoView({behavior: matchMedia('(prefers-reduced-motion: reduce)').matches ? 'auto' : 'smooth', block: 'center'});
        row.classList.remove('changed'); void row.offsetWidth; row.classList.add('changed');
      });
    });
  }

  // ---------------------------------------------------------------------------------------- sequences
  function sequences() {
    $('sequence-picker').innerHTML = S.sequences.map(s => `<button type="button" data-seq="${s.id}" aria-pressed="false">${esc(s.title)}</button>`).join('');
    document.querySelectorAll('[data-seq]').forEach(b => b.addEventListener('click', () => begin(b.dataset.seq)));
    $('sequence-next').addEventListener('click', () => {
      const s = seq.current;
      if (s && seq.index === s.steps.length - 1) begin(S.sequences[(S.sequences.indexOf(s) + 1) % S.sequences.length].id);
      else step(seq.index + 1);
    });
    $('sequence-prev').addEventListener('click', () => step(seq.index - 1));
    $('sequence-close').addEventListener('click', leave);
  }
  function begin(id) {
    seq.current = S.sequences.find(s => s.id === id);
    document.querySelectorAll('[data-seq]').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.seq === id)));
    $('sequence-player').hidden = false;
    step(0);
    $('instrument').scrollIntoView({behavior: matchMedia('(prefers-reduced-motion: reduce)').matches ? 'auto' : 'smooth', block: 'start'});
  }
  function apply(st) {
    if (st.kind === 'start') I().start(st.node, st.text);
    else if (st.kind === 'edge') { const e = I().byId[st.edge]; I().applyEdge(e, st.reverse, st.text); }
    else if (st.kind === 'prepare') { I().state.message = {text: st.text}; I().setPreparation(st.preparation); }
  }
  function step(k) {
    const s = seq.current;
    if (!s) return;
    k = Math.max(0, Math.min(s.steps.length - 1, k));
    if (k <= seq.index || seq.index < 0 || k > seq.index + 1) {   // replay from the start without animation delays
      for (let i = 0; i <= k; i++) apply(s.steps[i]);
    } else apply(s.steps[k]);
    seq.index = k;
    const st = s.steps[k];
    $('sequence-text').innerHTML = `<b>${k + 1} / ${s.steps.length}.</b> ${st.text}`;
    $('sequence-progress').innerHTML = s.steps.map((x, i) => {
      const slot = x.kind === 'edge' ? I().byId[x.edge].slot : x.kind === 'prepare' ? 'P' : null;
      return `<span class="${slot ? 'slot' : ''}${i <= k ? ' done' : ''}"${slot ? ` data-slot="${slot}"` : ''} title="${slot ? esc(S.slots[slot].name) : 'start'}"></span>`;
    }).join('');
    $('sequence-prev').disabled = k === 0;
    const lastStep = k === s.steps.length - 1;
    $('sequence-next').textContent = lastStep ? 'Next sequence' : 'Next step';
  }
  function leave() {
    seq.current = null; seq.index = -1;
    $('sequence-player').hidden = true;
    document.querySelectorAll('[data-seq]').forEach(b => b.setAttribute('aria-pressed', 'false'));
  }

  // ---------------------------------------------------------------------------------------- co-discovery
  const FIELD_COLORS = ['--xi', '--p', '--c', '--r', '--omega', '--a', '--good', '--warn'];
  function collapse() {
    const cd = S.codiscovery, grid = $('collapse-grid'), figs = [];
    const fieldsOf = rows => [...new Set(rows.map(r => r.field))];
    if (cd.rotation && cd.rotation.length) figs.push({id: 'rotation', title: 'Bloch rotation', rows: cd.rotation, draw: drawRotation,
      caption: rows => `<b>Bloch rotation.</b> The exact signal of ${rows.length} realizations from ${fieldsOf(rows).length} fields, rescaled by its own angle: (f − cos²θ)/sin²θ against |Ω|t lies on cos(|Ω|t). Largest deviation ${V.fmt(Math.max(...rows.map(r => r.residual)), 2)}.`});
    const T = [['symmetric-write', 'Symmetric write', 'εx − x³', 1.3313, 'π<sup>¼</sup> = 1.331', 'swept-write constant'],
               ['threshold-write', 'One-sided write', 'μ + x²', 1.0187929716, '|a₁′| = 1.0188', 'delay constant'],
               ['phase-locking', 'Phase locking', '−sin φ', 1.0, '1', 'half-width of the locking range in units of K']];
    for (const [key, title, form, expected, expText, constName] of T) {
      const rows = cd[key] || [];
      if (!rows.length) continue;
      figs.push({id: key, title, rows, expected, key, draw: drawCanonical,
        caption: rs => `<b>${title}.</b> ${rs.length} realizations from ${fieldsOf(rs).length} fields; the reduced drift at the end of each derivation in its canonical units, against ${form}. Below: the ${constName}, expected ${expText}.${rs.some(r => r.same_as) ? ' An open circle is a realization that simulates the equations of another: a replicate, not a separate model.' : ''}`});
    }
    const lawRefs = id => { const r = (S.law_reading || {})[id], M = window.FieldBridgeMechanisms; return r && M.reading ? `<span class="refs">Law ${M.reading(r)}.</span>` : ''; };
    grid.innerHTML = figs.map(f => `<figure><canvas id="cd-${f.id}" aria-label="${esc(f.title)}"></canvas>${f.key ? `<canvas id="cd-${f.id}-law" class="strip" aria-label="${esc(f.title)}: law constants" style="aspect-ratio:4/2.1"></canvas>` : ''}<figcaption>${f.caption(f.rows)} ${lawRefs(f.id)}</figcaption></figure>`).join('');
    const redraw = () => figs.forEach(f => f.draw(f));
    redraw();
    new ResizeObserver(redraw).observe(grid);
    // a point of a strip of law constants names its realization and opens it in the instrument
    grid.querySelectorAll('canvas.strip').forEach(c => {
      const near = ev => (c._points || []).find(p => Math.hypot(p.x - ev.offsetX, p.y - ev.offsetY) < 14);
      c.addEventListener('pointermove', ev => { const p = near(ev); c.style.cursor = p ? 'pointer' : 'default'; c.title = p ? p.label : ''; });
      c.addEventListener('click', ev => { const p = near(ev); if (p) open(p.id); });
    });
  }
  function color(i) { return getComputedStyle(document.documentElement).getPropertyValue(FIELD_COLORS[i % FIELD_COLORS.length]).trim(); }
  function drawRotation(f) {
    const canvas = $('cd-rotation'), g = V.fit(canvas), {ctx, w, h, c} = g, xmax = 4 * Math.PI;
    const fr = V.frame(g, {x: 0, y: 0, w, h}, [0, xmax], [-1.15, 1.15], {xLabel: '|Ω| t', yLabel: '(f − cos²θ) / sin²θ', yTicks: [-1, 0, 1], xTicks: [0, Math.PI, 2 * Math.PI, 3 * Math.PI, 4 * Math.PI], xFormat: v => v === 0 ? '0' : `${Math.round(v / Math.PI)}π`});
    const pts = []; for (let k = 0; k <= 200; k++) { const x = xmax * k / 200; pts.push([x, Math.cos(x)]); }
    V.line(ctx, pts, fr.X, fr.Y, c.faint, 5);
    f.rows.forEach((r, i) => {
      const th = r.theta * Math.PI / 180, s2 = Math.sin(th) ** 2, c2 = Math.cos(th) ** 2;
      if (s2 < 0.05) return;
      V.line(ctx, r.times.map((t, k) => [r.rate * t, (r.f[k] - c2) / s2]).filter(p => p[0] <= xmax), fr.X, fr.Y, color(i), 1.2);
    });
  }
  function drawCanonical(f) {
    const canvas = $('cd-' + f.id), g = V.fit(canvas), {ctx, w, h, c} = g;
    const rows = f.rows.filter(r => r.canonical);
    if (f.key === 'phase-locking') {
      const fr = V.frame(g, {x: 0, y: 0, w, h}, [-Math.PI, Math.PI], [-1.3, 1.3], {xLabel: 'canonical phase φ', yLabel: 'phase gained per period / K', yTicks: [-1, 0, 1], xTicks: [-Math.PI, 0, Math.PI], xFormat: v => v === 0 ? '0' : v > 0 ? 'π' : '−π'});
      const ref = []; for (let k = 0; k <= 120; k++) { const x = -Math.PI + 2 * Math.PI * k / 120; ref.push([x, -Math.sin(x)]); }
      V.line(ctx, ref, fr.X, fr.Y, c.faint, 5);
      rows.forEach((r, i) => r.canonical.phi.forEach((x, k) => V.dot(ctx, fr.X(x), fr.Y(r.canonical.g[k]), 2.2, color(i))));
    } else {
      const xs = rows.flatMap(r => r.canonical.x), lim = Math.max(0.5, Math.min(f.key === 'symmetric-write' ? 3 : 1.5, Math.max(...xs.map(Math.abs))));
      const refF = f.key === 'symmetric-write' ? x => -x * x * x : x => x * x;
      const ys = rows.flatMap(r => r.canonical.g.filter((_, k) => Math.abs(r.canonical.x[k]) <= lim));
      const yr = V.autoRange([...ys, refF(lim), refF(-lim)].filter(y => Math.abs(y) < 50), 0.08);
      const fr = V.frame(g, {x: 0, y: 0, w, h}, [-lim, lim], yr, {xLabel: 'x (canonical units)', yLabel: 'reduced drift'});
      const ref = []; for (let k = 0; k <= 100; k++) { const x = -lim + 2 * lim * k / 100; ref.push([x, refF(x)]); }
      V.line(ctx, ref, fr.X, fr.Y, c.faint, 5);
      rows.forEach((r, i) => V.line(ctx, r.canonical.x.map((x, k) => [x, r.canonical.g[k]]).filter(p => Math.abs(p[0]) <= lim), fr.X, fr.Y, color(i), 1.6));
    }
    const lawRows = f.rows.filter(r => r.law && Number.isFinite(r.law.constant));
    const strip = $('cd-' + f.id + '-law');
    if (!strip) return;
    const gs = V.fit(strip);
    if (!lawRows.length) { V.label(gs.ctx, 'law constants not computed for this build', 8, gs.h / 2, gs.c.muted, 'left', '12px'); return; }
    const vals = lawRows.flatMap(r => [r.law.constant - r.law.stderr, r.law.constant + r.law.stderr]).concat([f.expected]);
    const fr = V.frame(gs, {x: 0, y: 0, w: gs.w, h: gs.h}, [-0.5, lawRows.length - 0.5], V.autoRange(vals, 0.15), {xTicks: [], left: 54, bottom: 64, noGrid: true});
    gs.ctx.strokeStyle = gs.c.omega; gs.ctx.setLineDash([5, 4]); gs.ctx.beginPath(); gs.ctx.moveTo(fr.l, fr.Y(f.expected)); gs.ctx.lineTo(fr.r, fr.Y(f.expected)); gs.ctx.stroke(); gs.ctx.setLineDash([]);
    // a realization that simulates the equations of another is a replicate of that model, drawn open
    const sameAs = r => r.same_as && S.nodes[r.same_as] ? ` Same equations as the ${String(S.nodes[r.same_as].name).replace(/<[^>]+>/g, '')}: a replicate of that model.` : '';
    strip._points = lawRows.map((r, i) => ({x: fr.X(i), y: fr.Y(r.law.constant), id: r.id,
      label: `${String(r.name).replace(/<[^>]+>/g, '')}: ${V.pm(r.law.constant, r.law.stderr)}.${sameAs(r)} Select to open it above.`}));
    lawRows.forEach((r, i) => {
      const x = fr.X(i), col = color(f.rows.indexOf(r));
      gs.ctx.strokeStyle = col; gs.ctx.lineWidth = 1.5; gs.ctx.beginPath(); gs.ctx.moveTo(x, fr.Y(r.law.constant - r.law.stderr)); gs.ctx.lineTo(x, fr.Y(r.law.constant + r.law.stderr)); gs.ctx.stroke();
      if (r.same_as) V.dot(gs.ctx, x, fr.Y(r.law.constant), 3.2, gs.c.panel, col);
      else V.dot(gs.ctx, x, fr.Y(r.law.constant), 3, col);
      gs.ctx.save(); gs.ctx.translate(x, fr.b + 4); gs.ctx.rotate(-Math.PI / 4); gs.ctx.fillStyle = gs.c.muted; gs.ctx.font = '10px ' + getComputedStyle(document.body).fontFamily; gs.ctx.textAlign = 'right'; gs.ctx.textBaseline = 'middle';
      gs.ctx.fillText(S.nodes[r.id] && S.nodes[r.id].short || r.id, 0, 0); gs.ctx.restore();
    });
  }

  // ---------------------------------------------------------------------------------------- realizations
  const matrix = {field: '', query: '', selected: null};
  const searchable = r => strip(`${r.name} ${r.short} ${r.field} ${r.source || ''} ${S.classes[r.class] || ''}`).toLowerCase();
  function realizations() {
    const M = window.FieldBridgeMechanisms, all = Object.values(S.nodes);
    const fields = [...new Set(all.filter(r => !r.universal).map(r => r.field))].sort((a, b) => a.localeCompare(b));
    $('field-filter').innerHTML = `<button type="button" data-field="" aria-pressed="true">all fields</button>`
      + fields.map(f => `<button type="button" data-field="${esc(f)}" aria-pressed="false">${esc(f)}</button>`).join('');
    const chip = r => `<button type="button" class="rchip${r.universal ? ' universal' : ''}" data-id="${r.id}" data-field="${esc(r.universal ? '' : r.field)}" title="${esc(strip(r.name))}">${esc(r.short || r.id)}</button>`;
    $('matrix').innerHTML = Object.keys(S.mechanisms).map(k => {
      const m = S.mechanisms[k], rs = all.filter(r => r.class === k).sort((a, b) => strip(a.name).localeCompare(strip(b.name)));
      const groups = {};
      rs.filter(r => !r.universal).forEach(r => (groups[r.field] ||= []).push(r));
      const canon = rs.filter(r => r.universal);
      return `<div class="mrow" data-class="${k}">
        <div class="mhead">${M.glyph(k)}<div><b>${esc(S.classes[k] || k)}</b><span>${m.canonical || ''}</span></div></div>
        <div class="mcells">${canon.length ? `<span class="fgroup canonical"><span class="fname">canonical form</span>${canon.map(chip).join('')}</span>` : ''}${
          Object.keys(groups).sort((a, b) => a.localeCompare(b)).map(f => `<span class="fgroup" data-field="${esc(f)}"><span class="fname">${esc(f)}</span>${groups[f].map(chip).join('')}</span>`).join('')}</div></div>`;
    }).join('');
    document.querySelectorAll('#field-filter [data-field]').forEach(b => b.addEventListener('click', () => { matrix.field = b.dataset.field; filter(); }));
    $('matrix-search').addEventListener('input', ev => { matrix.query = ev.target.value.trim().toLowerCase(); filter(); });
    document.querySelectorAll('#matrix .rchip').forEach(b => b.addEventListener('click', () => select(b.dataset.id)));
    filter();
    sourcesTable();
  }
  function filter() {
    const {field, query} = matrix;
    let n = 0;
    const classes = new Set();
    document.querySelectorAll('#field-filter [data-field]').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.field === field)));
    document.querySelectorAll('#matrix .rchip').forEach(b => {
      const r = S.nodes[b.dataset.id];
      const on = (!field || b.dataset.field === field) && (!query || searchable(r).includes(query));
      b.classList.toggle('dim', !on);
      if (on) { n++; classes.add(r.class); }
    });
    document.querySelectorAll('#matrix .mrow').forEach(row => row.classList.toggle('dim', !classes.has(row.dataset.class)));
    document.querySelectorAll('#matrix .fgroup').forEach(g => g.classList.toggle('dim', !!field && g.dataset.field !== field));
    const total = Object.keys(S.nodes).length;
    $('matrix-summary').textContent = field || query
      ? `${n} of ${total} realizations${field ? ` in ${field}` : ''}${query ? ` matching “${query}”` : ''}, in ${classes.size} mechanism${classes.size === 1 ? '' : 's'}.`
      : `${total} realizations of ${Object.keys(S.mechanisms).length} mechanisms in ${new Set(Object.values(S.nodes).filter(r => !r.universal).map(r => r.field)).size} fields.`;
  }
  function select(id) {
    const r = S.nodes[id], M = window.FieldBridgeMechanisms;
    matrix.selected = id;
    document.querySelectorAll('#matrix .rchip').forEach(b => b.classList.toggle('selected', b.dataset.id === id));
    const box = $('matrix-detail'), row = document.querySelector(`#matrix .mrow[data-class="${r.class}"]`);
    if (row && row.insertAdjacentElement) row.insertAdjacentElement('afterend', box);  // under the row of the realization
    box.hidden = false;
    box.innerHTML = `${M.glyph(r.class, 'big')}<div><p class="card-label">${esc(S.classes[r.class] || r.class)}</p><h3>${r.name}</h3>
      <p class="muted">${r.universal ? 'canonical form, no field' : esc(r.field)}${r.question ? ' · ' + esc(r.question) : ''}</p>
      <p>${r.source ? esc(r.source) + ' · ' : ''}<a href="${repo(r.spec || r.base_spec)}" target="_blank" rel="noopener">${r.spec ? 'specification' : 'derived from ' + esc((r.base_spec || '').split('/').pop())}</a>${r.tutorial ? ` · <a href="${repo(r.tutorial)}" target="_blank" rel="noopener">tutorial</a>` : ''}</p>
      <button type="button" class="primary" id="matrix-open">Open it in the instrument</button></div>`;
    $('matrix-open').addEventListener('click', () => open(id));
  }
  function open(id) {
    leave();
    I().walkTo(id);
    $('instrument').scrollIntoView({behavior: matchMedia('(prefers-reduced-motion: reduce)').matches ? 'auto' : 'smooth', block: 'start'});
  }
  // ---------------------------------------------------------------------------------------- memory in model materials
  function memory() {
    const M = window.FieldBridgeMechanisms, ids = (S.materials || []).filter(id => S.nodes[id]);
    const box = $('memory-groups');
    if (!ids.length) { $('memory').hidden = true; return; }
    const fields = new Set(ids.map(id => S.nodes[id].field)), laws = {};
    ids.forEach(id => { const l = S.nodes[id].facts.loss || 'not calculated'; laws[l] = (laws[l] || 0) + 1; });
    $('memory-summary').innerHTML = `${ids.length} materials from ${fields.size} fields. Retention: `
      + Object.entries(laws).sort((a, b) => b[1] - a[1]).map(([l, n]) => `<span class="law"><b>${n}</b> ${M.lawLinked(l)}</span>`).join('');
    const R = S.retention_laws;
    if (R) $('memory-laws').innerHTML = `<p>A stored state is lost by one of three laws, set by the form of the landscape at the state (<a href="${repo(R.link)}" target="_blank" rel="noopener">Module 4, Section 1.3 ↗</a>):</p><ul>`
      + Object.entries(R.laws).map(([n, text]) => `<li><b>Law ${n}</b> ${text}${R.reading && R.reading[n] ? ` <span class="refs">(${M.reading(R.reading[n], {only: ['sources']})})</span>` : ''}</li>`).join('')
      + `</ul><p>Laws 1 and 2 are the loss for the two signs κ &gt; 0 and κ = 0 of the curvature along the written direction. Where κ &lt; 0 the information decreases only while the state is near the unstable point: once the expansion outruns the noise it stops decreasing, at ½ ln(1 + W), and this is how a state is written. Law 3 needs wells separated by a barrier. A pattern written in a field is lost mode by mode, each mode relaxing at its own rate κ(k) by Law 1; a conserved density, whose rates vanish at long wavelengths, loses it as a power of time (the two cards of loss in a field on the <a href="#instrument">map</a>). A state is written where a rate κ passes through zero as a control is varied, at the write points of the map; the writes have laws of their own, whose constants are the same in every field (<a href="#codiscovery">One mechanism in different fields</a>).</p>`;
    const groups = [];
    ids.forEach(id => { const k = S.nodes[id].class, g = groups.find(x => x.k === k); if (g) g.ids.push(id); else groups.push({k, ids: [id]}); });
    const tile = id => {
      const r = S.nodes[id], f = r.facts, card = r.card && r.card.image;
      const thumb = card ? `<span class="thumb"><img src="${esc(card)}" alt="" loading="lazy"></span>` : `<span class="thumb glyph-only">${M.glyph(r.class, 'big')}</span>`;
      // the retention law of a material without a stable state already says so
      const states = f.states != null && !String(f.loss || '').startsWith('no stable state') ? `${f.states_text || f.states} stable state${f.states === 1 ? '' : 's'}` : '';
      return `<article class="material" data-class="${r.class}">
        <button type="button" class="material-open" data-open="${id}" aria-label="Open ${esc(strip(r.name))} in the instrument">${thumb}
          <span class="m-text"><b>${r.name}</b><span class="m-field">${esc(r.field)}</span>${r.question ? `<span class="m-question">${esc(r.question)}</span>` : ''}</span>
        </button>
        <p class="m-facts">${states ? esc(states) + ' · ' : ''}${M.lawLinked(f.loss || '')}</p>
        ${card ? `<a class="m-card" href="${esc(card)}" target="_blank" rel="noopener">memory card ↗</a>` : ''}</article>`;
    };
    const absent = S.absent || {};
    box.innerHTML = groups.map(g => `<div class="memory-group" data-class="${g.k}">
        <h3>${M.glyph(g.k, 'small')}${esc(S.classes[g.k] || g.k)}${absent[g.k] ? `<span class="absent-tag" title="${esc(absent[g.k][1])}">${esc(absent[g.k][0])}</span>` : ''}<span class="muted">${g.ids.length} material${g.ids.length === 1 ? '' : 's'}</span></h3>
        <div class="memory-grid">${g.ids.map(tile).join('')}</div></div>`).join('');
    box.querySelectorAll('[data-open]').forEach(b => b.addEventListener('click', () => open(b.dataset.open)));
  }
  // ---------------------------------------------------------------------------------------- mechanisms in preparation
  function planned() {
    const M = window.FieldBridgeMechanisms, items = S.planned || [];
    if (!items.length) { $('planned').hidden = true; return; }
    $('planned-count').textContent = ['No', 'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven', 'Eight', 'Nine'][items.length] || String(items.length);
    $('planned-grid').innerHTML = items.map(p => `<article class="planned-card" data-planned="${p.id}">
        <div class="planned-head">${M.glyph(p.id)}<div><p class="card-label">${esc(p.group)}</p><h3>${esc(p.name)}</h3></div></div>
        <p class="planned-law">${p.law}</p>
        <dl><dt>What exists</dt><dd>${p.exists}${p.tutorial ? ` (<a href="${repo(p.tutorial)}" target="_blank" rel="noopener">tutorial ↗</a>)` : ''}</dd>
        <dt>What is missing</dt><dd>${p.missing}</dd>
        ${p.sources && p.sources.length ? `<dt>Original publications</dt><dd>${M.reading(p, {max: 6, only: ['sources', 'doc'], prefix: ''})}</dd>` : ''}</dl></article>`).join('');
  }
  function sourcesTable() {
    const order = Object.keys(S.mechanisms);
    const rows = Object.values(S.nodes).sort((a, b) => order.indexOf(a.class) - order.indexOf(b.class) || strip(a.name).localeCompare(strip(b.name)));
    $('sources-table').innerHTML = '<thead><tr><th>Realization</th><th>Field</th><th>Source and files</th></tr></thead><tbody>'
      + rows.map(r => `<tr><td><button type="button" class="plain" data-go="${r.id}">${r.name}</button><br><span class="muted">${esc(S.classes[r.class] || '')}</span></td><td>${r.universal ? 'canonical form' : esc(r.field)}</td><td>${r.source ? esc(r.source) + '<br>' : ''}<a href="${repo(r.spec || r.base_spec)}" target="_blank" rel="noopener">${r.spec ? 'specification' : 'derived from ' + esc((r.base_spec || '').split('/').pop())}</a>${r.tutorial ? ` · <a href="${repo(r.tutorial)}" target="_blank" rel="noopener">tutorial</a>` : ''}</td></tr>`).join('') + '</tbody>';
    document.querySelectorAll('[data-go]').forEach(b => b.addEventListener('click', () => open(b.dataset.go)));
    $('boundary').textContent = S.boundary;
    const law = S.law || {};
    $('provenance').textContent = `Calculated with FieldBridge (memory implementation ${S.provenance.memory_implementation_sha256.slice(0, 12)}; Python ${S.provenance.versions.python}, numpy ${S.provenance.versions.numpy}, scipy ${S.provenance.versions.scipy}). `
      + (law.computed ? 'Law constants computed for this build.' : law.record_implementation_sha256 ? `Law constants from docs/site/law_constants.json, computed with memory implementation ${law.record_implementation_sha256.slice(0, 12)} for unchanged specifications.` : 'Law constants not computed for this build.');
  }

  document.addEventListener('DOMContentLoaded', () => {
    definitions(); sequences(); realizations(); memory(); planned();
    // ?r=<realization> opens that realization directly, for links from the tutorial
    const want = new URLSearchParams(location.search).get('r');
    I().start(want && S.nodes[want] ? want : S.start);
    collapse();
  });
  window.FieldBridgeSite = {leaveSequence: leave, begin};
})();
