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
  };
  function glyph(klass, cls = '') {
    return `<svg class="glyph ${cls}" viewBox="0 0 48 32" aria-hidden="true" focusable="false">${GLYPHS[klass] || GLYPHS['single-state']}</svg>`;
  }

  // ---------------------------------------------------------------------------------------- the map
  // columns by family: closed evolution | relaxation and writing | cycles | fields and noise
  const W = 1256, H = 318, CARD = {w: 136, h: 92};
  const COL = [78, 230, 392, 544, 696, 858, 1020, 1172], ROW = [78, 222];
  const PLACE = {
    'rotation': [0, 0], 'conserved': [0, 1], 'obstructed': [1, 0],
    'symmetric-write': [2, 0], 'single-state': [2, 1], 'threshold-write': [3, 0], 'subcritical-write': [3, 1], 'field-write': [4, 1],
    'oscillation': [5, 0], 'neutral-cycles': [5, 1], 'exponential-loss': [6, 0], 'power-loss': [6, 1], 'convention': [7, 0.5],
  };
  const FAMILIES = [{label: 'closed evolution', from: 0, to: 1}, {label: 'relaxation and writing', from: 2, to: 4},
                    {label: 'cycles', from: 5, to: 5}, {label: 'fields and noise', from: 6, to: 7}];
  const at = k => { const [c, r] = PLACE[k]; return {x: COL[c], y: ROW[0] + (ROW[1] - ROW[0]) * r}; };

  /** The single-component changes between classes: {a, b, slots: Set, edges: [...]} with a, b in reading order. */
  function links(S) {
    const out = new Map(), keys = Object.keys(PLACE);
    for (const e of S.edges) {
      const ca = S.nodes[e.from].class, cb = S.nodes[e.to].class;
      if (ca === cb || !PLACE[ca] || !PLACE[cb]) continue;
      const [a, b] = keys.indexOf(ca) < keys.indexOf(cb) ? [ca, cb] : [cb, ca], id = a + '|' + b;
      if (!out.has(id)) out.set(id, {a, b, slots: new Set(), edges: []});
      out.get(id).slots.add(e.slot); out.get(id).edges.push(e.id);
    }
    return [...out.values()];
  }
  function path(a, b, bend) {
    const p = at(a), q = at(b), mx = (p.x + q.x) / 2, my = (p.y + q.y) / 2;
    const dx = q.x - p.x, dy = q.y - p.y, n = Math.hypot(dx, dy) || 1;
    const cx = mx - dy / n * bend, cy = my + dx / n * bend;
    return {d: `M${p.x},${p.y} Q${cx},${cy} ${q.x},${q.y}`, lx: (p.x + 2 * cx + q.x) / 4, ly: (p.y + 2 * cy + q.y) / 4};
  }
  function counts(S, klass) {
    const rs = Object.values(S.nodes).filter(r => r.class === klass);
    const fields = new Set(rs.filter(r => !r.universal).map(r => r.field));
    return {n: rs.length, fields: fields.size, canonical: rs.some(r => r.universal)};
  }
  function draw() {
    const S = window.FIELDBRIDGE_SITE, svg = document.getElementById('mechanisms');
    if (!S || !svg) return;
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    const present = Object.keys(S.mechanisms || {}).filter(k => PLACE[k]);
    const fam = FAMILIES.map(f => {
      const x0 = COL[f.from] - CARD.w / 2 - 10, x1 = COL[f.to] + CARD.w / 2 + 10;
      return `<g class="family"><rect x="${x0}" y="4" width="${x1 - x0}" height="${H - 8}" rx="14"/><text x="${x0 + 12}" y="${H - 14}">${esc(f.label)}</text></g>`;
    }).join('');
    // a line along a row that would cross a card of that row is bent around it: above the upper row, below the lower
    const crosses = (a, b) => { const [ca, ra] = PLACE[a], [cb, rb] = PLACE[b];
      return ra === rb && present.some(k => PLACE[k][1] === ra && PLACE[k][0] > Math.min(ca, cb) && PLACE[k][0] < Math.max(ca, cb)); };
    const lines = links(S).map(l => {
      const slots = [...l.slots], dir = Math.sign(COL[PLACE[l.b][0]] - COL[PLACE[l.a][0]]) || 1;
      const around = crosses(l.a, l.b) ? (PLACE[l.a][1] === 0 ? -1 : 1) * dir * 118 : 0;
      return slots.map((slot, i) => {
        const bend = around + (i - (slots.length - 1) / 2) * 16, g = path(l.a, l.b, bend);
        return `<g class="link" data-a="${l.a}" data-b="${l.b}" data-slot="${slot}"><path d="${g.d}"/>`
          + `<circle cx="${g.lx}" cy="${g.ly}" r="10"/><text x="${g.lx}" y="${g.ly + 4.5}" text-anchor="middle">${esc(S.slots[slot].symbol)}</text>`
          + `<title>${esc(S.slots[slot].name)}: ${esc(S.classes[l.a])} ↔ ${esc(S.classes[l.b])}</title></g>`;
      }).join('');
    }).join('');
    const cards = present.map(k => {
      const p = at(k), c = counts(S, k), m = S.mechanisms[k];
      const sub = `${c.n} realization${c.n === 1 ? '' : 's'}${c.fields ? ` · ${c.fields} field${c.fields === 1 ? '' : 's'}` : ''}`;
      return `<g class="mech" data-class="${k}" tabindex="0" role="button" aria-label="${esc(S.classes[k])}: ${esc(strip(m.text || ''))}">
        <rect x="${p.x - CARD.w / 2}" y="${p.y - CARD.h / 2}" width="${CARD.w}" height="${CARD.h}" rx="12"/>
        <svg x="${p.x - 24}" y="${p.y - CARD.h / 2 + 9}" width="48" height="32" viewBox="0 0 48 32" class="glyph">${GLYPHS[k] || ''}</svg>
        <text class="name" x="${p.x}" y="${p.y + 16}" text-anchor="middle">${esc(S.classes_short[k] || k)}</text>
        <text class="sub" x="${p.x}" y="${p.y + 33}" text-anchor="middle">${esc(sub)}</text>
        <title>${esc(S.classes[k])}. ${esc(strip(m.text || ''))} Canonical form: ${esc(strip(m.canonical || ''))}.</title></g>`;
    }).join('');
    svg.innerHTML = fam + `<g class="links">${lines}</g><g class="mechs">${cards}</g>`;
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
    highlight();
  }
  /** Current mechanism filled; the mechanisms one change away outlined in the colour of that change. */
  function highlight(extra) {
    const S = window.FIELDBRIDGE_SITE, I = window.FieldBridgeInstrument, svg = document.getElementById('mechanisms');
    if (!S || !I || !svg || !I.state.node) return;
    const here = S.nodes[I.state.node].class, next = new Map();
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
      l.classList.toggle('dim', !!lit && l.dataset.slot !== lit);
    });
  }
  /** A change under the pointer: its target mechanism is marked on the map. */
  function preview(klass) {
    const svg = document.getElementById('mechanisms');
    if (!svg) return;
    svg.querySelectorAll('.mech').forEach(g => g.classList.toggle('preview', !!klass && g.dataset.class === klass));
  }
  document.addEventListener('DOMContentLoaded', () => {
    draw();
    window.FieldBridgeInstrument && window.FieldBridgeInstrument.onChange((st, extra) => highlight(extra));
  });
  window.FieldBridgeMechanisms = {glyph, draw, highlight, preview, links};
})();
