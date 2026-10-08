/* Mechanisms: a drawing for each class of mechanism, and the map of mechanisms. A node of the map is a mechanism
 * written without a field; a line joins two mechanisms when a verified change of one component (in its colour) takes a
 * realization of the one to a realization of the other. Selecting a mechanism walks the shortest path of changes to
 * the realization the map opens for it: its canonical form where the page has one. */
(function () {
  'use strict';
  const esc = s => String(s).replace(/[&<>"']/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[c]));
  const strip = s => String(s).replace(/<[^>]+>/g, '');

  // ---------------------------------------------------------------------------------------- drawings, 48 × 32
  const AXES = '<path class="axis" d="M4 4 V28 H44"/>';
  const GLYPHS = {
    'rotation': '<ellipse cx="24" cy="20" rx="18" ry="7"/><path d="M24 20 L34 5"/><circle class="dot" cx="34" cy="5" r="2.4"/><path d="M13.5 26.4 l-3.8 -1.3 l2.6 -2.9"/>',
    'conserved': AXES + '<path d="M4 11 H44"/><circle class="dot" cx="17" cy="11" r="2.4"/><circle class="dot" cx="31" cy="11" r="2.4"/>',
    'obstructed': AXES + '<path d="M5 16 C8 5 11 5 14 14 C16 21 18 26 21 18 C23 12 25 7 28 13 C31 21 33 27 37 18 C40 11 42 13 44 15"/>',
    'single-state': '<path d="M6 5 C12 31 36 31 42 5"/><circle class="dot" cx="24" cy="24.3" r="3.2"/>',
    'symmetric-write': '<path d="M3 5 C8 31 16 31 24 15 C32 31 40 31 45 5"/><circle class="dot" cx="12.4" cy="23.6" r="3.2"/>',
    'threshold-write': '<path d="M3 4 C8 33 16 33 23 17 C28 24 32 22 35 17 C38 12 42 8 45 4"/><circle class="dot" cx="11.6" cy="24.6" r="3.2"/>',
    'subcritical-write': '<path d="M3 4 C6 31 11 31 15 16 C18 24 30 24 33 16 C37 31 42 31 45 4"/><circle class="dot" cx="24" cy="18.6" r="3.2"/><circle class="dot faint" cx="9.2" cy="23.5" r="2.4"/>',
    'field-write': '<path d="M9 5 L15 17 M21 5 L27 17 M33 5 L39 17"/><path d="M6 25 H40"/><path d="M36 21.5 L41 25 L36 28.5"/>',
    'oscillation': '<ellipse cx="24" cy="16" rx="17" ry="10.5"/><path d="M25.5 2.8 l4.4 2.7 l-4.4 2.8"/>',
    'neutral-cycles': '<ellipse cx="24" cy="16" rx="19" ry="12.5"/><ellipse cx="24" cy="16" rx="12" ry="7.8"/><ellipse cx="24" cy="16" rx="5" ry="3.2"/>',
    'exponential-loss': AXES + '<path d="M5 6 C11 22 18 26 44 26.5"/>',
    'power-loss': '<path class="axis" d="M4 28 H44"/><path class="faint" d="M15 27.5 C19.5 27.5 21 5 24 5 C27 5 28.5 27.5 33 27.5"/><path d="M4 27 C13 26 16 17 24 17 C32 17 35 26 44 27"/>',
    'convention': '<path d="M4 25 L9 22 L13 24 L18 17 L22 19 L27 12 L31 14 L36 8 L40 9 L44 5"/><path class="faint" d="M4 25 L9 24 L13 26 L18 21 L22 23 L27 19 L31 21 L36 17 L40 19 L44 16"/>',
    // mechanisms in preparation
    'frustrated-loops': '<path class="faint" d="M11 25 L24 6 L37 25 Z"/><path d="M5 21 L17 29"/><path d="M18 7 L30 5"/><path d="M31 29 L43 21"/>',
    'retention-rewriting': AXES + '<path d="M5 26 C18 25 30 18 44 5"/><path class="faint" d="M5 26 C14 20 24 16 44 14"/>',
    'hopf-onset': '<ellipse class="faint" cx="24" cy="16" rx="17" ry="11"/><path d="M24 16 c2 -1 4 1 2 3 c-3 2 -7 -1 -5 -5 c3 -4 10 -2 10 4 c0 6 -8 9 -13 5"/>',
    'kuramoto': '<ellipse class="faint" cx="24" cy="16" rx="16" ry="12"/><circle class="dot" cx="38" cy="12" r="2.3"/><circle class="dot" cx="39.6" cy="17.5" r="2.3"/><circle class="dot" cx="35" cy="7" r="2.3"/><circle class="dot" cx="37" cy="23" r="2.3"/><circle class="dot faint" cx="9" cy="21" r="2.3"/>',
  };
  // the return to a turning point: an excursion that closes on its turning point, or one that ends beside it
  GLYPHS['return-point'] = AXES + '<path d="M6 26 C22 26 24 6 43 6"/><path class="faint" d="M43 6 C27 6 25 26 6 26"/><path d="M18 20 C25 17 28 14 33 12 C27 16 24 18 18 20"/><circle class="dot" cx="33" cy="12" r="2.2"/>';
  GLYPHS['no-return'] = AXES + '<path d="M6 26 C22 26 24 6 43 6"/><path class="faint" d="M43 6 C27 6 25 26 6 26"/><path d="M18 21 C25 18 28 15 33 12 C28 17 26 20 23 23"/><circle class="dot" cx="33" cy="12" r="2.2"/><circle class="dot faint" cx="23" cy="23" r="2.2"/>';
  // the return to a set point: the output after a step of the input returns to its level, part of the way, or not at all
  GLYPHS['perfect-adaptation'] = AXES + '<path class="faint" d="M6 20 H43"/><path d="M6 20 H12 C14 8 17 6 20 12 C24 20 30 20 43 20"/>';
  GLYPHS['fine-tuned-adaptation'] = AXES + '<path class="faint" d="M6 20 H43"/><path d="M6 20 H12 C14 8 17 6 20 12 C24 20 30 20 43 20"/><circle class="dot" cx="38" cy="20" r="2.2"/>';
  GLYPHS['partial-adaptation'] = AXES + '<path class="faint" d="M6 20 H43"/><path d="M6 20 H12 C14 8 17 6 20 11 C24 16 30 15 43 15"/>';
  GLYPHS['no-adaptation'] = AXES + '<path class="faint" d="M6 20 H43"/><path d="M6 20 H12 C16 10 22 9 43 9"/>';
  // computation: the capacity against the delay (bars), at degree 1 only, odd degrees, all degrees, or none
  GLYPHS['linear-memory'] = AXES + '<path d="M9 30 V8 M15 30 V13 M21 30 V18 M27 30 V22 M33 30 V25 M39 30 V27"/>';
  GLYPHS['odd-capacity'] = AXES + '<path d="M9 30 V9 M15 30 V15 M21 30 V20"/><path class="faint" d="M27 30 V26 M33 30 V28"/><circle class="dot" cx="39" cy="22" r="2.2"/>';
  GLYPHS['nonlinear-capacity'] = AXES + '<path d="M9 30 V10 M15 30 V16 M21 30 V21"/><path class="faint" d="M27 30 V18 M33 30 V22 M39 30 V25"/>';
  GLYPHS['integrating'] = AXES + '<path d="M6 28 C14 26 22 20 30 14 C35 10 39 8 43 7"/><path class="faint" d="M9 30 V29 M15 30 V29 M21 30 V29"/>';
  // the quantum write: a swept parametric oscillator chooses one of two states; in the linear stage the chosen branch
  // grows from the seed, at a few photons the two shallow wells exchange population before the choice is frozen
  GLYPHS['linear-stage-write'] = AXES + '<path class="faint" d="M6 12 C16 12 22 17 43 20"/><path d="M6 12 C16 12 22 7 43 4"/><circle class="dot" cx="10" cy="12" r="1.6"/>';
  GLYPHS['equilibrium-write'] = AXES + '<path d="M6 6 C12 22 18 22 24 9 C30 22 36 22 43 6"/><circle class="dot" cx="15" cy="18" r="2"/><circle class="dot faint" cx="33" cy="18" r="2"/>';
  function glyph(klass, cls = '') {
    return `<svg class="glyph ${cls}" viewBox="0 0 48 32" aria-hidden="true" focusable="false">${GLYPHS[klass] || GLYPHS['single-state']}</svg>`;
  }

  // ---------------------------------------------------------------------------------------- the map
  // Two arrangements of the same cards. Wide: columns by family, the families grouped in modules. Narrow: two columns,
  // placed so that every line joins neighbouring cards and none crosses another. The map shows one module at a time,
  // the module of the selected mechanism unless another is chosen, or every module.
  const WIDE = {
    w: 2374, h: 318, card: {w: 124, h: 92}, col: [76, 240, 414, 578, 742, 916, 1090, 1264, 1428, 1602, 1776, 1950, 2124, 2298], row: [78, 222],
    place: {'rotation': [0, 0], 'conserved': [0, 1], 'obstructed': [1, 0],
            'symmetric-write': [2, 0], 'single-state': [2, 1], 'threshold-write': [3, 0], 'subcritical-write': [3, 1], 'field-write': [4, 1],
            'return-point': [5, 0], 'no-return': [5, 1],
            'oscillation': [6, 0], 'neutral-cycles': [6, 1], 'exponential-loss': [7, 0], 'power-loss': [7, 1], 'convention': [8, 0.5],
            'perfect-adaptation': [9, 0], 'fine-tuned-adaptation': [9, 1], 'partial-adaptation': [10, 0], 'no-adaptation': [10, 1],
            'linear-memory': [11, 0], 'odd-capacity': [11, 1], 'integrating': [12, 0], 'nonlinear-capacity': [12, 1],
            'linear-stage-write': [13, 0], 'equilibrium-write': [13, 1]},
    families: [{label: 'quantum: closed evolution', from: 0, to: 1, module: 'quantum'},
               {label: 'memory: writing a state', from: 2, to: 4, module: 'memory'},
               {label: 'memory: turning points', from: 5, to: 5, module: 'memory'}, {label: 'memory: phase', from: 6, to: 6, module: 'memory'},
               {label: 'memory: retention', from: 7, to: 7, module: 'memory'}, {label: 'noise', from: 8, to: 8, module: 'memory'},
               {label: 'regulation: set point', from: 9, to: 10, module: 'regulation'},
               {label: 'computation: capacity', from: 11, to: 12, module: 'computation'},
               {label: 'quantum: open evolution', from: 13, to: 13, module: 'quantum'}],
  };
  const NARROW = {
    w: 360, h: 1324, card: {w: 150, h: 78}, col: [86, 274], row: [50, 152],
    place: {'conserved': [0, 0], 'rotation': [1, 0], 'obstructed': [0, 1], 'single-state': [1, 1],
            'subcritical-write': [0, 2], 'symmetric-write': [1, 2], 'threshold-write': [0, 3], 'oscillation': [1, 3],
            'field-write': [0, 4], 'neutral-cycles': [1, 4], 'exponential-loss': [0, 5], 'power-loss': [1, 5], 'convention': [0, 6],
            'return-point': [1, 6], 'no-return': [0, 7],
            'perfect-adaptation': [0, 8], 'fine-tuned-adaptation': [1, 8], 'partial-adaptation': [0, 9], 'no-adaptation': [1, 9],
            'linear-memory': [0, 10], 'integrating': [1, 10], 'odd-capacity': [0, 11], 'nonlinear-capacity': [1, 11],
            'linear-stage-write': [0, 12], 'equilibrium-write': [1, 12]},
    families: [],
  };
  const MODULES = [{key: 'memory', label: 'memory'}, {key: 'quantum', label: 'quantum'}, {key: 'regulation', label: 'regulation'},
                   {key: 'computation', label: 'computation'}, {key: 'all', label: 'all mechanisms'}];
  const moduleOf = k => { const f = WIDE.place[k] && WIDE.families.find(f => WIDE.place[k][0] >= f.from && WIDE.place[k][0] <= f.to); return f ? f.module : null; };
  let shown = 'memory';
  let layout = WIDE;
  let view = {place: {}, col: WIDE.col, row: WIDE.row, w: WIDE.w, h: WIDE.h, families: []};
  const visible = k => shown === 'all' || moduleOf(k) === shown;
  /** The cards shown, packed: the wide map drops the columns without a card and keeps the spacing of neighbouring
   *  columns (a wider gap between families); the narrow map drops the rows without a card. */
  function arrange(keys) {
    const wide = layout === WIDE, axis = wide ? 0 : 1;
    const used = [...new Set(keys.map(k => layout.place[k][axis]))].sort((a, b) => a - b);
    const index = new Map(used.map((v, i) => [v, i]));
    const place = {};
    keys.forEach(k => { const [c, r] = layout.place[k]; place[k] = wide ? [index.get(c), r] : [c, index.get(r)]; });
    if (!wide) {
      const pitch = layout.row[1] - layout.row[0];
      return {place, col: layout.col, row: layout.row, w: layout.w, h: used.length ? 2 * layout.row[0] + pitch * (used.length - 1) : layout.h, families: []};
    }
    const col = [];
    used.forEach((c, i) => col.push(i === 0 ? layout.col[0] : col[i - 1] + (c === used[i - 1] + 1 ? layout.col[c] - layout.col[c - 1] : 174)));
    const families = layout.families.map(f => ({label: f.label, cols: used.filter(c => c >= f.from && c <= f.to)}))
      .filter(f => f.cols.length).map(f => ({label: f.label, from: index.get(f.cols[0]), to: index.get(f.cols[f.cols.length - 1])}));
    return {place, col, row: layout.row, w: used.length ? col[col.length - 1] + layout.col[0] : layout.w, h: layout.h, families};
  }
  const at = k => { const [c, r] = view.place[k]; return {x: view.col[c], y: view.row[0] + (view.row[1] - view.row[0]) * r}; };
  const ORDER = Object.keys(WIDE.place);

  /** The single-component changes between classes: {a, b, slots: [...], edges: [...]} with a, b in reading order. */
  function links(S) {
    const out = new Map();
    for (const e of S.edges) {
      const ca = S.nodes[e.from].class, cb = S.nodes[e.to].class;
      if (ca === cb || !WIDE.place[ca] || !WIDE.place[cb]) continue;
      const [a, b] = ORDER.indexOf(ca) < ORDER.indexOf(cb) ? [ca, cb] : [cb, ca], id = a + '|' + b;
      if (!out.has(id)) out.set(id, {a, b, slots: [], edges: []});
      const l = out.get(id);
      if (!l.slots.includes(e.slot)) l.slots.push(e.slot);
      l.edges.push(e.id);
    }
    return [...out.values()];
  }
  /** A line from card to card, bent by `bend` perpendicular to it; the label sits on its middle. */
  function path(a, b, bend, shift = 0) {
    const p = at(a), q = at(b), dx = q.x - p.x, dy = q.y - p.y, n = Math.hypot(dx, dy) || 1;
    const ox = -dy / n, oy = dx / n, mx = (p.x + q.x) / 2, my = (p.y + q.y) / 2;
    const cx = mx + ox * bend, cy = my + oy * bend;
    return {d: `M${p.x + ox * shift},${p.y + oy * shift} Q${cx + ox * shift},${cy + oy * shift} ${q.x + ox * shift},${q.y + oy * shift}`,
            lx: (p.x + 2 * cx + q.x) / 4, ly: (p.y + 2 * cy + q.y) / 4};
  }
  function counts(S, klass) {
    const rs = Object.values(S.nodes).filter(r => r.class === klass);
    const fields = new Set(rs.filter(r => !r.universal).map(r => r.field));
    return {n: rs.length, fields: fields.size};
  }
  function draw() {
    const S = window.FIELDBRIDGE_SITE, svg = document.getElementById('mechanisms');
    if (!S || !svg) return;
    layout = (svg.parentElement && svg.parentElement.clientWidth && svg.parentElement.clientWidth < 640) ? NARROW : WIDE;
    const present = Object.keys(S.mechanisms || {}).filter(k => layout.place[k] && visible(k));
    view = arrange(present);
    const {w: W, h: H} = view, CARD = layout.card;
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    svg.classList.toggle('narrow', layout === NARROW);
    svg.dataset.module = shown;
    // the wide map is drawn at its design size when the frame is wide enough and scales down below that
    svg.style.width = layout === WIDE ? `${W}px` : '';
    svg.style.minWidth = layout === WIDE ? `${Math.min(860, W)}px` : '';
    const fam = view.families.map(f => {
      const x0 = view.col[f.from] - CARD.w / 2 - 14, x1 = view.col[f.to] + CARD.w / 2 + 14;
      return `<g class="family"><rect x="${x0}" y="4" width="${x1 - x0}" height="${H - 8}" rx="14"/><text x="${x0 + 8}" y="${H - 14}">${esc(f.label)}</text></g>`;
    }).join('');
    // a line along a row that would cross a card of that row is bent around it: above the upper row, below the lower
    const crosses = (a, b) => { const [ca, ra] = view.place[a], [cb, rb] = view.place[b];
      return ra === rb && present.some(k => view.place[k][1] === ra && view.place[k][0] > Math.min(ca, cb) && view.place[k][0] < Math.max(ca, cb)); };
    const lines = links(S).filter(l => view.place[l.a] && view.place[l.b]).map(l => {
      const dir = Math.sign(at(l.b).x - at(l.a).x) || 1;
      const around = crosses(l.a, l.b) ? (view.place[l.a][1] === 0 ? -1 : 1) * dir * 118 : 0;
      const n = l.slots.length, mid = path(l.a, l.b, around), pw = 18 * n + 4;
      // one line for each component, side by side; one label for the pair, with the symbol of each component
      const strokes = l.slots.map((slot, i) => `<path data-slot="${slot}" d="${path(l.a, l.b, around, (i - (n - 1) / 2) * 4).d}"/>`).join('');
      const symbols = l.slots.map((slot, i) => `<text data-slot="${slot}" x="${mid.lx + (i - (n - 1) / 2) * 18}" y="${mid.ly + 4.5}" text-anchor="middle">${esc(S.slots[slot].symbol)}</text>`).join('');
      const names = l.slots.map(slot => S.slots[slot].name).join(' or ');
      return `<g class="link" data-a="${l.a}" data-b="${l.b}" data-slots="${l.slots.join(' ')}"${n === 1 ? ` data-slot="${l.slots[0]}"` : ''}>${strokes}`
        + `<rect class="tag" x="${mid.lx - pw / 2}" y="${mid.ly - 10}" width="${pw}" height="20" rx="10"/>${symbols}`
        + `<title>A change of the ${esc(names)}: ${esc(S.classes[l.a])} ↔ ${esc(S.classes[l.b])}</title></g>`;
    }).join('');
    const cards = present.map(k => {
      const p = at(k), c = counts(S, k), m = S.mechanisms[k], absent = (S.absent || {})[k], narrow = layout === NARROW;
      const count = c.fields ? `${c.n} in ${c.fields} field${c.fields === 1 ? '' : 's'}` : `${c.n} realization${c.n === 1 ? '' : 's'}`;
      // an outcome without the mechanism of its group is hatched and labelled: on the wide map in a third line, on the
      // narrow one before the count
      const sub = absent && narrow ? `${absent[0]} · ${count}` : count, lift = absent && !narrow ? 5 : 0;
      return `<g class="mech${absent ? ' absent' : ''}" data-class="${k}" tabindex="0" role="button" aria-label="${esc(S.classes[k])}${absent ? `, ${esc(absent[0])}` : ''}: ${esc(strip(m.text || ''))}">
        <rect x="${p.x - CARD.w / 2}" y="${p.y - CARD.h / 2}" width="${CARD.w}" height="${CARD.h}" rx="12"/>
        <svg x="${p.x - 24}" y="${p.y - CARD.h / 2 + (narrow ? 5 : 9) - lift}" width="48" height="32" viewBox="0 0 48 32" class="glyph">${GLYPHS[k] || ''}</svg>
        <text class="name${(S.classes_short[k] || k).length > 16 ? ' long' : ''}" x="${p.x}" y="${p.y + (narrow ? 13 : 16) - lift}" text-anchor="middle">${esc(S.classes_short[k] || k)}</text>
        <text class="sub" x="${p.x}" y="${p.y + (narrow ? 29 : 33) - lift}" text-anchor="middle">${esc(sub)}</text>
        ${absent && !narrow ? `<text class="absent-label" x="${p.x}" y="${p.y + 40}" text-anchor="middle">${esc(absent[0])}</text>` : ''}
        <title>${esc(S.classes[k])}${absent ? ` (${esc(absent[0])}: ${esc(absent[1])})` : ''}. ${esc(strip(m.text || ''))} Canonical form: ${esc(strip(m.canonical || ''))}. ${c.n} realization${c.n === 1 ? '' : 's'}${c.fields ? ` in ${c.fields} field${c.fields === 1 ? '' : 's'}` : ''}.</title></g>`;
    }).join('');
    const hatch = '<defs><pattern id="absent-hatch" width="7" height="7" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
      + '<rect class="hatch-ground" width="7" height="7"/><path class="hatch-line" d="M0 0 V7"/></pattern></defs>';
    svg.innerHTML = hatch + fam + `<g class="links">${lines}</g><g class="mechs">${cards}</g>`;
    // a family label wider than its box (a family of one column) widens the box, and the map when the box is the last
    let right = W;
    svg.querySelectorAll('.family').forEach(g => {
      const t = g.querySelector('text'), r = g.querySelector('rect');
      if (!t || !r || typeof t.getComputedTextLength !== 'function') return;
      const need = Math.ceil(t.getComputedTextLength()) + 16, x = Number(r.getAttribute('x'));
      if (need > Number(r.getAttribute('width'))) { r.setAttribute('width', need); right = Math.max(right, x + need); }
    });
    if (right > W) { svg.setAttribute('viewBox', `0 0 ${right} ${H}`); if (layout === WIDE) svg.style.width = `${right}px`; }
    svg.querySelectorAll('.mech').forEach(g => {
      const go = () => {
        const I = window.FieldBridgeInstrument, m = S.mechanisms[g.dataset.class];
        if (!I || !m) return;
        window.FieldBridgeSite && window.FieldBridgeSite.leaveSequence();
        I.walkTo(m.node);
      };
      g.addEventListener('click', go);
      g.addEventListener('keydown', ev => { if (ev.key === 'Enter' || ev.key === ' ') { ev.preventDefault(); go(); } });
    });
    chips();
    highlight({keep: true});
  }
  /** The buttons above the map: one for each module with mechanisms on the map, and one for all of them. */
  function chips() {
    const S = window.FIELDBRIDGE_SITE, box = document.getElementById('map-modules');
    if (!S || !box) return;
    const all = Object.keys(S.mechanisms || {}).filter(k => WIDE.place[k]);
    box.innerHTML = MODULES.map(m => ({...m, n: m.key === 'all' ? all.length : all.filter(k => moduleOf(k) === m.key).length}))
      .filter(m => m.n).map(m => `<button type="button" class="module${m.key === shown ? ' on' : ''}" data-module="${m.key}" `
        + `aria-pressed="${m.key === shown}" title="${m.key === 'all' ? 'Every mechanism on the map' : `The mechanisms of the ${esc(m.label)} module`}">`
        + `${esc(m.label)} <small>${m.n}</small></button>`).join('');
    box.querySelectorAll('.module').forEach(b => b.addEventListener('click', () => module(b.dataset.module)));
  }
  /** The module shown on the map; with a key, shows that module ('all' for every mechanism) and returns it. */
  function module(key) {
    if (key !== undefined && key !== shown && MODULES.some(m => m.key === key)) { shown = key; draw(); }
    return shown;
  }
  /** Current mechanism filled; the mechanisms one change away outlined in the colour of that change. When the
   *  selection moves to a mechanism outside the module shown, the map follows it (not after a choice by hand: `keep`). */
  function highlight(extra) {
    const S = window.FIELDBRIDGE_SITE, I = window.FieldBridgeInstrument, svg = document.getElementById('mechanisms');
    if (!S || !I || !svg || !I.state.node) return;
    const here = S.nodes[I.state.node].class, next = new Map();
    if (!(extra && extra.keep) && shown !== 'all' && moduleOf(here) && moduleOf(here) !== shown) { shown = moduleOf(here); draw(); return; }
    for (const e of (I.out[I.state.node] || [])) next.set(S.nodes[e.to].class, e.slot);
    for (const e of (I.into[I.state.node] || [])) next.set(S.nodes[e.from].class, e.slot);
    const seen = new Set(I.state.path.map(p => S.nodes[p.node].class));
    const lit = extra && extra.light;
    svg.querySelectorAll('.mech').forEach(g => {
      const k = g.dataset.class;
      g.classList.toggle('current', k === here);
      g.classList.toggle('next', k !== here && next.has(k));
      g.classList.toggle('seen', k !== here && seen.has(k));
      if (next.has(k) && k !== here) g.setAttribute('data-slot', next.get(k)); else g.removeAttribute('data-slot');
    });
    svg.querySelectorAll('.link').forEach(l => {
      const touches = l.dataset.a === here || l.dataset.b === here;
      l.classList.toggle('near', touches);
      l.classList.toggle('dim', !!lit && !(l.dataset.slots || '').split(' ').includes(lit));
    });
  }
  /** A change under the pointer: its target mechanism is marked on the map, or, when the target lies outside the
   *  module shown, the button of its module. */
  function preview(klass) {
    const svg = document.getElementById('mechanisms'), box = document.getElementById('map-modules');
    if (svg) svg.querySelectorAll('.mech').forEach(g => g.classList.toggle('preview', !!klass && g.dataset.class === klass));
    if (box) box.querySelectorAll('.module').forEach(b => b.classList.toggle('preview', !!klass && !visible(klass) && b.dataset.module === moduleOf(klass)));
  }
  document.addEventListener('DOMContentLoaded', () => {
    draw();
    window.FieldBridgeInstrument && window.FieldBridgeInstrument.onChange((st, extra) => highlight(extra));
    // the arrangement follows the width of the frame
    const svg = document.getElementById('mechanisms');
    if (svg && svg.parentElement && window.ResizeObserver)
      new ResizeObserver(() => { const narrow = svg.parentElement.clientWidth < 640; if (narrow !== (layout === NARROW)) draw(); }).observe(svg.parentElement);
  });
  /** A text that names a law of retention, with every '(Law N ...)' linked to its definition (Module 4, Section 1.3)
   *  and explained on hover. The text is escaped. */
  function lawLinked(text) {
    const S = window.FIELDBRIDGE_SITE, R = S && S.retention_laws;
    const html = esc(text || '');
    if (!R) return html;
    const href = (S.repo || 'https://github.com/synthetix-institute/fieldbridge') + '/blob/main/' + R.link;
    return html.replace(/\(Law ([123])([^)]*)\)/g, (all, n, rest) => R.laws[n]
      ? `(<a class="law-link" href="${href}" target="_blank" rel="noopener" title="Law ${n}: ${esc(strip(R.laws[n].replace(/<sup>/g, '^').replace(/<\/?su[bp]>/g, '')).replace(/&gt;/g, '>').replace(/&lt;/g, '<'))}">Law ${n}</a>${rest})` : all);
  }
  /** Where a mechanism or a law is defined and derived, and its original publications (site_references.py): links to
   *  the tutorial, to the DOI of each publication and to its section of docs/mechanisms.md. `only` selects parts. */
  function reading(r, {max = 5, only = ['defined', 'derived', 'sources', 'doc'], prefix = 'original publications: '} = {}) {
    if (!r) return '';
    const S = window.FIELDBRIDGE_SITE, base = (S && S.repo || 'https://github.com/synthetix-institute/fieldbridge') + '/blob/main/';
    const plain = t => strip(t || '').replace(/&gt;/g, '>').replace(/&lt;/g, '<').replace(/&amp;/g, '&');
    const a = (href, label, title) => `<a href="${esc(href)}" target="_blank" rel="noopener"${title ? ` title="${esc(plain(title))}"` : ''}>${label}</a>`;
    const parts = [];
    if (only.includes('defined') && r.defined) parts.push(`defined in ${a(base + r.defined.link, esc(r.defined.label), r.defined.what)}`);
    if (only.includes('derived') && r.derived && r.derived.length)
      parts.push(`derived in ${r.derived.map(d => a(base + d.link, esc(d.label), d.what)).join(', ')}`);
    const src = r.sources || [];
    if (only.includes('sources') && src.length)
      parts.push(`${prefix}${src.slice(0, max).map(s => s.url ? a(s.url, esc(s.short), s.cite + (s.what ? ': ' + s.what : ''))
        : `<span title="${esc(plain(s.cite + (s.what ? ': ' + s.what : '')))}">${esc(s.short)}</span>`).join('; ')}${src.length > max ? '; …' : ''}`);
    if (only.includes('doc') && r.doc) parts.push(a(base + r.doc, 'all references ↗', 'docs/mechanisms.md'));
    return parts.join(' · ');
  }
  window.FieldBridgeMechanisms = {glyph, draw, highlight, preview, links, lawLinked, reading, module, modules: MODULES.map(m => m.key)};
})();
