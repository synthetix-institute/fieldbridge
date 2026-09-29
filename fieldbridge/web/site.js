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
        caption: rs => `<b>${title}.</b> ${rs.length} realizations from ${fieldsOf(rs).length} fields; the reduced drift at the end of each derivation in its canonical units, against ${form}. Below: the ${constName}, expected ${expText}.`});
    }
    grid.innerHTML = figs.map(f => `<figure><canvas id="cd-${f.id}" aria-label="${esc(f.title)}"></canvas>${f.key ? `<canvas id="cd-${f.id}-law" class="strip" aria-label="${esc(f.title)}: law constants" style="aspect-ratio:4/2.1"></canvas>` : ''}<figcaption>${f.caption(f.rows)}</figcaption></figure>`).join('');
    const redraw = () => figs.forEach(f => f.draw(f));
    redraw();
    new ResizeObserver(redraw).observe(grid);
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
    lawRows.forEach((r, i) => {
      const x = fr.X(i), col = color(f.rows.indexOf(r));
      gs.ctx.strokeStyle = col; gs.ctx.lineWidth = 1.5; gs.ctx.beginPath(); gs.ctx.moveTo(x, fr.Y(r.law.constant - r.law.stderr)); gs.ctx.lineTo(x, fr.Y(r.law.constant + r.law.stderr)); gs.ctx.stroke();
      V.dot(gs.ctx, x, fr.Y(r.law.constant), 3, col);
      gs.ctx.save(); gs.ctx.translate(x, fr.b + 4); gs.ctx.rotate(-Math.PI / 4); gs.ctx.fillStyle = gs.c.muted; gs.ctx.font = '10px ' + getComputedStyle(document.body).fontFamily; gs.ctx.textAlign = 'right'; gs.ctx.textBaseline = 'middle';
      gs.ctx.fillText(S.nodes[r.id] && S.nodes[r.id].short || r.id, 0, 0); gs.ctx.restore();
    });
  }

  // ---------------------------------------------------------------------------------------- sources
  function sources() {
    const order = ['unitary', 'dissipative', 'field', 'stochastic'];
    const rows = Object.values(S.nodes).sort((a, b) => order.indexOf(a.family) - order.indexOf(b.family) || strip(a.name).localeCompare(strip(b.name)));
    $('sources-table').innerHTML = '<thead><tr><th>Realization</th><th>Field</th><th>Source and files</th></tr></thead><tbody>'
      + rows.map(r => `<tr><td><button type="button" class="plain" data-go="${r.id}">${r.name}</button><br><span class="muted">${esc(S.classes[r.class] || '')}</span></td><td>${esc(r.field)}</td><td>${r.source ? esc(r.source) + '<br>' : ''}<a href="${repo(r.spec || r.base_spec)}" target="_blank" rel="noopener">${r.spec ? 'specification' : 'derived from ' + esc((r.base_spec || '').split('/').pop())}</a>${r.tutorial ? ` · <a href="${repo(r.tutorial)}" target="_blank" rel="noopener">tutorial</a>` : ''}</td></tr>`).join('') + '</tbody>';
    document.querySelectorAll('[data-go]').forEach(b => b.addEventListener('click', () => { leave(); I().walkTo(b.dataset.go); $('instrument').scrollIntoView({block: 'start'}); }));
    $('boundary').textContent = S.boundary;
    const law = S.law || {};
    $('provenance').textContent = `Calculated with FieldBridge (memory implementation ${S.provenance.memory_implementation_sha256.slice(0, 12)}; Python ${S.provenance.versions.python}, numpy ${S.provenance.versions.numpy}, scipy ${S.provenance.versions.scipy}). `
      + (law.computed ? 'Law constants computed for this build.' : law.record_implementation_sha256 ? `Law constants from docs/site/law_constants.json, computed with memory implementation ${law.record_implementation_sha256.slice(0, 12)} for unchanged specifications.` : 'Law constants not computed for this build.');
  }

  document.addEventListener('DOMContentLoaded', () => {
    definitions(); sequences(); sources();
    // ?r=<realization> opens that realization directly, for links from the tutorial
    const want = new URLSearchParams(location.search).get('r');
    I().start(want && S.nodes[want] ? want : S.start);
    collapse();
  });
  window.FieldBridgeSite = {leaveSequence: leave, begin};
})();
