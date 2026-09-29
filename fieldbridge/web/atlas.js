/* The map: every realization, grouped by mechanism, and every single-component change between them. */
(function () {
  'use strict';
  const S = window.FIELDBRIDGE_SITE;
  if (!S) return;
  const NS = 'http://www.w3.org/2000/svg';
  const esc = s => String(s).replace(/[&<>"']/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[c]));
  const strip = s => String(s).replace(/<[^>]+>/g, '');
  let pos = {}, lastWidth = 0;

  // the columns of one mechanism family share a band; a column with many realizations wraps into sub-columns
  const BANDS = [['rotation', 'conserved', 'obstructed'], ['single-state', 'symmetric-write', 'threshold-write', 'subcritical-write', 'field-write'],
                 ['oscillation', 'neutral-cycles'], ['exponential-loss', 'power-loss', 'convention']];
  const BAND_TITLES = ['quantum: closed unitary evolution', 'memory: dissipative, with a bath', 'oscillators', 'fields and conventions'];
  function layout(width) {
    const cols = Object.fromEntries(S.atlas.columns.map(c => [c.class, c])), narrow = width < 760;
    pos = {};
    const labels = [], rowH = 23, maxRows = narrow ? 99 : 8;
    let y = 8;
    BANDS.forEach((band, b) => {
      const present = band.filter(k => cols[k]);
      if (!present.length) return;
      labels.push({x: 8, y: y + 12, text: BAND_TITLES[b], anchor: 'start', band: true});
      y += 26;
      if (narrow) {
        present.forEach(k => {
          const c = cols[k], perRow = Math.max(1, Math.floor((width - 16) / 150)), cellW = (width - 16) / perRow;
          labels.push({x: 8, y: y + 8, text: S.classes_short[k] || k, full: S.classes[k], anchor: 'start'});
          y += 22;
          c.nodes.forEach((id, i) => { pos[id] = {x: 14 + cellW * (i % perRow), y: y + rowH * Math.floor(i / perRow), band: b, col: k}; });
          y += rowH * Math.ceil(c.nodes.length / perRow) + 6;
        });
        y += 8;
        return;
      }
      const subs = present.map(k => Math.ceil(cols[k].nodes.length / maxRows)), total = subs.reduce((a, v) => a + v, 0);
      const unit = (width - 16) / total;
      let x = 8, deepest = 0;
      present.forEach((k, i) => {
        const c = cols[k], rows = Math.ceil(c.nodes.length / subs[i]);
        labels.push({x: x + 6, y: y + 8, text: S.classes_short[k] || k, full: S.classes[k], anchor: 'start'});
        c.nodes.forEach((id, j) => { pos[id] = {x: x + 12 + unit * Math.floor(j / rows), y: y + 30 + rowH * (j % rows), band: b, col: k}; });
        deepest = Math.max(deepest, rows);
        x += unit * subs[i];
      });
      y += 30 + rowH * deepest + 14;
    });
    return {width, height: y, narrow, labels};
  }
  function curve(a, b) {
    const ax = a.x + 4, bx = b.x + 4;
    if (Math.abs(a.y - b.y) < 1) return `M${ax},${a.y} Q${(ax + bx) / 2},${a.y - 16} ${bx},${b.y}`;
    if (Math.abs(ax - bx) < 1) { const bend = 26 + Math.abs(a.y - b.y) * 0.15; return `M${ax},${a.y} C${ax - bend},${a.y} ${bx - bend},${b.y} ${bx},${b.y}`; }
    const my = (a.y + b.y) / 2;
    return a.band === b.band ? `M${ax},${a.y} C${(ax + bx) / 2},${a.y} ${(ax + bx) / 2},${b.y} ${bx},${b.y}` : `M${ax},${a.y} C${ax},${my} ${bx},${my} ${bx},${b.y}`;
  }
  function draw() {
    const svg = document.getElementById('atlas');
    if (!svg) return;
    const width = Math.max(320, svg.parentElement.clientWidth - 16);
    lastWidth = width;
    const L = layout(width);
    svg.setAttribute('viewBox', `0 0 ${L.width} ${L.height}`);
    svg.setAttribute('height', L.height);
    const edges = S.edges.filter(e => pos[e.from] && pos[e.to]).map(e =>
      `<path class="edge" data-edge="${e.id}" data-slot="${e.slot}" d="${curve(pos[e.from], pos[e.to])}"><title>${esc(strip(e.change))}</title></path>`).join('');
    const nodes = Object.entries(pos).map(([id, p]) => {
      const r = S.nodes[id];
      return `<g class="node" data-node="${id}" tabindex="0" role="button" aria-label="${esc(strip(r.name))}: ${esc(S.classes[r.class] || '')}">
        <circle cx="${p.x + 4}" cy="${p.y}" r="5"/><text x="${p.x + 13}" y="${p.y + 4}">${esc(r.short || id)}</text>
        <title>${esc(strip(r.name))} · ${esc(r.field)}</title></g>`;
    }).join('');
    const labels = L.labels.map(l => `<text class="${l.band ? 'band-label' : 'col-label'}" x="${l.x}" y="${l.y}" text-anchor="${l.anchor || 'middle'}"><title>${esc(l.full || l.text)}</title>${esc(l.text)}</text>`).join('');
    svg.innerHTML = `<g class="edges">${edges}</g>${labels}<g class="nodes">${nodes}</g>`;
    svg.querySelectorAll('.node').forEach(g => {
      const go = () => window.FieldBridgeInstrument && (window.FieldBridgeSite && window.FieldBridgeSite.leaveSequence(), window.FieldBridgeInstrument.walkTo(g.dataset.node));
      g.addEventListener('click', go);
      g.addEventListener('keydown', ev => { if (ev.key === 'Enter' || ev.key === ' ') { ev.preventDefault(); go(); } });
    });
    highlight();
  }
  function highlight(extra) {
    const I = window.FieldBridgeInstrument, svg = document.getElementById('atlas');
    if (!I || !svg) return;
    const st = I.state, onPath = new Set(st.path.filter(p => p.edge).map(p => p.edge)), visited = new Set(st.path.map(p => p.node));
    svg.querySelectorAll('.node').forEach(g => { g.classList.toggle('current', g.dataset.node === st.node); g.classList.toggle('visited', visited.has(g.dataset.node) && g.dataset.node !== st.node); });
    const lit = extra && extra.light;
    svg.querySelectorAll('.edge').forEach(p => { p.classList.toggle('on-path', onPath.has(p.dataset.edge)); p.classList.toggle('dim', !!lit && p.dataset.slot !== lit); });
  }
  function legend() {
    const el = document.getElementById('map-legend');
    if (el) el.innerHTML = Object.entries(S.slots).map(([k, s]) => `<span data-slot="${k}"><i></i>${s.symbol} ${esc(s.name)}</span>`).join('')
      + '<span>● current realization · thick lines: the path of changes</span>';
  }
  document.addEventListener('DOMContentLoaded', () => {
    legend(); draw();
    window.FieldBridgeInstrument && window.FieldBridgeInstrument.onChange((st, extra) => highlight(extra));
    let t = null;
    new ResizeObserver(() => { clearTimeout(t); t = setTimeout(() => { const svg = document.getElementById('atlas'); if (svg && Math.abs(svg.parentElement.clientWidth - 16 - lastWidth) > 20) draw(); }, 120); }).observe(document.getElementById('atlas').parentElement);
  });
  window.FieldBridgeAtlas = {draw, highlight};
})();
