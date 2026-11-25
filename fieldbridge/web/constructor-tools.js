/* Construction tools for the realizations written as equations: a term added to one equation with a check of the
 * carrier, the consequence of an edit shown before it is kept, the history of the construction with export and
 * replay, a plane of two parameters shaded by the number of stable states, and a protocol of segments for one
 * parameter. Everything is computed in the browser with model-physics.js and dissipative.js; a realization built
 * here is not a FieldBridge specification and is marked as such. */
(function () {
  'use strict';
  const S = window.FIELDBRIDGE_SITE, I = window.FieldBridgeInstrument, D = window.FieldBridgeDissipative;
  const P = window.FieldBridgeModels, MML = window.FieldBridgeMath, V = window.FieldBridgeViews;
  if (!S || !I || !D || !P || !MML) return;
  const $ = id => document.getElementById(id);
  const fmt = (v, d = 3) => (V && V.fmt ? V.fmt(v, d) : Number(v).toPrecision(d));
  const esc = s => String(s).replace(/[&<>"']/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[c]));
  const strip = s => String(s).replace(/<[^>]+>/g, '');
  const css = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim() || '#888';
  const RESERVED = new Set(['E', 'N', 'pi', 'e', 'i']);
  const usable = rec => rec && rec.family === 'dissipative' && rec.engine && rec.engine.kind_dynamics === 'equations';

  // ------------------------------------------------------------------------------------------- terms and types
  const num = c => ['number', c], name = n => ['name', n];
  /** Candidate monomials for the equation of variable i, by the kind of carrier. Each has a key (for replay), a
   *  tree without its coefficient and, for an amount, the sign of coefficient that keeps the amount non-negative. */
  function monomials(e, i) {
    const v = e.variables, n = v.length, kind = e.carrier.kind, out = [];
    const add = (key, tree, label) => out.push({key, tree, label});
    add('1', null, 'constant');
    for (let j = 0; j < n; j++) add(`x${j}`, name(v[j]), 'linear');
    for (let j = 0; j < n; j++) add(`x${j}^2`, ['pow', name(v[j]), num(2)], 'quadratic');
    for (let j = 0; j < n; j++) for (let k = j + 1; k < n && out.length < 40; k++) add(`x${j}x${k}`, ['mul', name(v[j]), name(v[k])], 'product');
    for (let j = 0; j < n; j++) add(`x${j}^3`, ['pow', name(v[j]), num(3)], 'cubic');
    for (let j = 0; j < n; j++) add(`sin${j}`, ['sin', name(v[j])], 'periodic');
    for (let j = 0; j < n; j++) add(`cos${j}`, ['cos', name(v[j])], 'periodic');
    if (kind === 'torus') for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) if (j !== k) add(`sin${j}-${k}`, ['sin', ['sub', name(v[j]), name(v[k])]], 'coupling of phases');
    return out;
  }
  const withCoef = (c, tree) => (tree ? ['mul', name(c), tree] : name(c));
  /** The check of a term against the carrier. Returns {ok, sign, reason}: sign is +1 or -1 when only one sign of the
   *  coefficient is allowed. The check evaluates the term; it does not read its label. */
  function check(e, i, mono, params) {
    const v = e.variables, n = v.length, c = e.carrier, random = D.rng(97), scale = c.scale || 1;
    const env = q => ({...params, ...Object.fromEntries(v.map((x, k) => [x, q[k]]))});
    const t = mono.tree || num(1);
    const val = q => P.value(t, env(q));
    if (c.kind === 'torus') {
      const per = c.period || 2 * Math.PI;
      for (let s = 0; s < 12; s++) {
        const q = v.map(() => random() * per);
        for (let j = 0; j < n; j++) {
          const q2 = q.slice(); q2[j] += per;
          if (Math.abs(val(q2) - val(q)) > 1e-9 * Math.max(1, Math.abs(val(q)))) {
            return {ok: false, reason: `${v[j]} is an angle on a circle of period ${fmt(per, 3)}. The term changes when ${v[j]} goes once around the circle, so it is not a function on the carrier. An angle enters only through periodic functions such as sin and cos.`};
          }
        }
      }
      return {ok: true, reason: 'periodic in every angle: a function on the torus.'};
    }
    if (c.kind === 'orthant') {
      let lo = Infinity, hi = -Infinity;
      for (let s = 0; s < 40; s++) {
        const q = v.map(() => 2 * scale * random()); q[i] = 0;
        const x = val(q); if (!Number.isFinite(x)) return {ok: false, reason: 'the term is not finite at some non-negative states.'};
        lo = Math.min(lo, x); hi = Math.max(hi, x);
      }
      if (Math.abs(lo) < 1e-12 && Math.abs(hi) < 1e-12) return {ok: true, reason: `the term vanishes at ${v[i]} = 0, so it cannot make ${v[i]} negative: either sign of the coefficient is allowed.`};
      if (lo >= -1e-12) return {ok: true, sign: 1, reason: `at ${v[i]} = 0 the term is positive. It is a source of ${v[i]}, and its coefficient must be non-negative, or the term would drive ${v[i]} below zero.`};
      if (hi <= 1e-12) return {ok: true, sign: -1, reason: `at ${v[i]} = 0 the term is negative, so its coefficient must be non-positive for ${v[i]} to stay an amount.`};
      return {ok: false, reason: `at ${v[i]} = 0 the term takes both signs, so no coefficient keeps ${v[i]} non-negative. A loss of an amount needs a factor ${v[i]}.`};
    }
    for (let s = 0; s < 20; s++) {
      const q = v.map(() => (2 * random() - 1) * 2 * scale);
      if (!Number.isFinite(val(q))) return {ok: false, reason: 'the term is not finite at some states.'};
    }
    return {ok: true, reason: 'a smooth function of the state on the carrier.'};
  }
  function freshName(e) {
    const taken = new Set([...Object.keys(e.params || {}), ...e.variables]);
    for (let k = 1; ; k++) { const c = 'c' + k; if (!taken.has(c) && !RESERVED.has(c)) return c; }
  }
  /** A realization built from rec by one added term: stored among the nodes so that the constructor treats it like
   *  any other, with the realization it came from as its reference. */
  function localNode(rec, i, mono, value) {
    const e = rec.engine, c = freshName(e), id = `local:${rec.id}|${i}|${mono.key}|${c}`;
    const st = I.state, params = {...(st.node === rec.id ? st.params : e.params), [c]: value};
    const chk = check(e, i, mono, params), mag = Math.max(1, 2 * Math.abs(value));
    const range = chk.sign === 1 ? [0, mag, 0.005] : chk.sign === -1 ? [-mag, 0, 0.005] : [-mag, mag, 0.005];
    const term = {id: `${e.variables[i]}-add-${c}`, variable: i, label: `${e.variables[i]}: added ${mono.label} term`,
                  text: `${c} * ${mono.key}`, tree: withCoef(c, mono.tree)};
    const engine = {...e, params, ranges: {...(e.ranges || {}), [c]: range}, terms: [...e.terms, term],
                    links: [...(e.links || []), ...(mono.tree ? referenced(e, mono.tree).map(j => ({source: j, target: i, term: term.id})) : [])]};
    delete engine.potential; delete engine._frame;
    const base = rec.local ? rec.local.base : rec;
    const added = mono.tree ? `${c}·${monoText(e, mono)}` : c;
    const node = {...rec, id, engine, name: `${rec.name} + ${esc(added)}`, short: 'built here',
                  facts: {states: null, loss: ''}, spec: null, base_spec: base.spec || base.base_spec, tutorial: null, card: null,
                  local: {base, parent: rec.id, add: {variable: i, key: mono.key, coefficient: c, value},
                          text: `Built here from ${strip(rec.name)} by adding ${added} to d${e.variables[i]}/dt. The carrier allows the term: ${chk.reason} These equations are not a FieldBridge specification and have not been verified; the stable states and the motion shown are computed in the browser for the equations as written.`}};
    S.nodes[id] = node;
    return node;
  }
  const referenced = (e, tree) => e.variables.map((v, j) => (JSON.stringify(tree).includes(`["name","${v}"]`) ? j : -1)).filter(j => j >= 0);
  function monoText(e, mono) {
    if (!mono.tree) return '1';
    const t = s => (s[0] === 'name' ? s[1] : s[0] === 'number' ? String(s[1]) : s[0] === 'pow' ? `${t(s[1])}^${t(s[2])}`
      : s[0] === 'mul' ? `${t(s[1])} ${t(s[2])}` : s[0] === 'sub' ? `${t(s[1])} − ${t(s[2])}` : `${s[0]}(${t(s[1])})`);
    return t(mono.tree);
  }
  /** Add the term to the current realization and move to the result, as one step of the construction. */
  function addTerm(i, key, value) {
    const rec = S.nodes[I.state.node]; if (!usable(rec)) return false;
    const mono = monomials(rec.engine, i).find(m => m.key === key); if (!mono) return false;
    const chk = check(rec.engine, i, mono, {...I.state.params, c: value}); if (!chk.ok) return false;
    if (chk.sign && Math.sign(value) && Math.sign(value) !== chk.sign) value = -value;
    const node = localNode(rec, i, mono, value), from = rec.id;
    I.state.path.push({node: node.id, add: {variable: i, key, value}, local: true});
    I.load(node.id, {text: null, from, slot: 'Omega'});
    return true;
  }

  // ------------------------------------------------------------------------------------------- stable states
  function stableStates(e, params, removed, seeds) {
    const pts = D.equilibria(e, params, removed || [], seeds || e.states || [], e.variables.length > 4 ? 6 : 10, 5).filter(p => p.stable);
    if (pts.length) return {points: pts.map(p => p.q), moving: 0};
    const a = D.attractors(e, params, removed || [], {n: 6, T: 50, seed: 5});
    return {points: a.points.map(p => p.q), moving: a.moving || 0, unbounded: a.unbounded || 0};
  }
  const describe = (e, st) => !st.points.length ? (st.unbounded ? 'no stable state: preparations run away' : 'no stable state of rest: the preparations keep moving')
    : `${st.points.length} stable state${st.points.length === 1 ? '' : 's'}: ${st.points.slice(0, 4).map(q => '(' + q.map(x => fmt(x, 3)).join(', ') + ')').join(', ')}${st.points.length > 4 ? ', …' : ''}`;

  // ------------------------------------------------------------------------------------------- history
  const history = [];
  let lastPath = 0, lastNode = null;
  function logStep(kind, text) { history.push({kind, text, node: I.state.node, time: new Date().toISOString()}); renderHistory(); }
  function observe(st, info) {
    if (info && info.light !== undefined) return;
    const n = st.path.length;
    if (n === lastPath && st.node === lastNode) return;
    const p = st.path[n - 1], rec = S.nodes[st.node];
    if (n === 1) logStep('start', `start at ${strip(rec.name)}`);
    else if (n < lastPath) logStep('undo', `back to ${strip(rec.name)}`);
    else if (p.edge) { const e = I.byId[p.edge]; logStep('verified', `${strip(p.reverse ? 'reverse of ' + e.change : e.change)} → ${strip(rec.name)} (verified change of ${S.slots[e.slot].name})`); }
    else if (p.add) logStep('built', `added ${strip(rec.name.split(' + ').pop())} to d${rec.engine.variables[p.add.variable]}/dt (built here, not verified)`);
    lastPath = n; lastNode = st.node;
    onNode(rec);
  }
  function exportConstruction() {
    const st = I.state;
    return {fieldbridge_construction: 1, start: st.path[0].node,
            steps: st.path.slice(1).map(p => (p.edge ? {edge: p.edge, reverse: !!p.reverse} : {add: p.add})),
            params: {...st.params}, protocol: committed ? {param: committed.param, segments: committed.segments, variable: committed.variable} : null,
            log: history.map(h => ({kind: h.kind, text: h.text}))};
  }
  function replay(doc) {
    if (!doc || doc.fieldbridge_construction !== 1 || !S.nodes[doc.start]) throw new Error('not a FieldBridge construction');
    I.start(doc.start);
    for (const s of doc.steps || []) {
      if (s.edge) { const e = I.byId[s.edge]; if (!e || !I.applyEdge(e, s.reverse)) throw new Error('the change ' + s.edge + ' does not leave ' + I.state.node); }
      else if (s.add) { if (!addTerm(s.add.variable, s.add.key, s.add.value)) throw new Error('the term ' + s.add.key + ' cannot be added here'); }
    }
    if (doc.params) { Object.assign(I.state.params, doc.params); I.render(); }
    if (doc.protocol) { draft = {...doc.protocol, segments: doc.protocol.segments.map(s => ({...s}))}; committed = draft; draft = clone(committed); drawProtocol(); }
    logStep('replay', `replayed a construction of ${(doc.steps || []).length} step${(doc.steps || []).length === 1 ? '' : 's'}`);
  }
  function renderHistory() {
    const list = $('ct-history'); if (!list) return;
    list.innerHTML = history.map(h => `<li class="ct-${h.kind}"><span class="ct-tag">${{start: 'start', undo: 'undo', verified: 'verified', built: 'built here', plane: 'values', protocol: 'protocol', replay: 'replay'}[h.kind] || h.kind}</span>${esc(h.text)}</li>`).join('');
    list.scrollTop = list.scrollHeight;
  }

  // ------------------------------------------------------------------------------------------- the panel
  function panel() {
    const host = document.querySelector('#instrument .expression-card'); if (!host || $('ct')) return;
    const el = document.createElement('div');
    el.className = 'ct'; el.id = 'ct';
    el.innerHTML = `
      <div class="ct-tabs" role="tablist">
        <button type="button" role="tab" data-ct="add" aria-selected="true">Add a term</button>
        <button type="button" role="tab" data-ct="plane" aria-selected="false">Two parameters</button>
        <button type="button" role="tab" data-ct="protocol" aria-selected="false">Protocol</button>
        <button type="button" role="tab" data-ct="history" aria-selected="false">History</button>
        <button type="button" role="tab" data-ct="examples" aria-selected="false">Examples</button>
      </div>
      <div class="ct-pane" data-pane="add">
        <p class="ct-note">A term is added to one equation of the drift Ω. The carrier Ξ decides which terms are allowed: an angle enters only through periodic functions, and the equation of an amount may contain only a term that cannot drive that amount below zero. The consequence is computed before the term is kept.</p>
        <div class="ct-row"><label for="ct-eq">equation</label><select id="ct-eq"></select>
          <label for="ct-coef">coefficient</label><input id="ct-coef" type="number" step="0.05" value="0.1"></div>
        <div class="ct-terms" id="ct-terms"></div>
        <div class="ct-preview" id="ct-preview" aria-live="polite"></div>
        <div class="ct-row"><button type="button" class="primary" id="ct-keep" disabled>Add the term</button></div>
      </div>
      <div class="ct-pane" data-pane="plane" hidden>
        <p class="ct-note">Each cell is shaded by the number of stable states found by Newton's method from several starts at those two values, with all other values as set in A. A click sets both values.</p>
        <div class="ct-row"><label for="ct-px">horizontal</label><select id="ct-px"></select><label for="ct-py">vertical</label><select id="ct-py"></select>
          <button type="button" id="ct-plane-run">Compute</button></div>
        <canvas id="ct-plane" class="ct-canvas" role="img" aria-label="Number of stable states over two parameters"></canvas>
        <p class="ct-legend" id="ct-plane-legend"></p>
      </div>
      <div class="ct-pane" data-pane="protocol" hidden>
        <p class="ct-note">A protocol holds one parameter at a sequence of values. The preparations start in each stable state at the values set in A; the drawing shows one variable against time, dashed for the protocol being edited and solid for the one applied, and states where each preparation ends.</p>
        <div class="ct-row"><label for="ct-pp">parameter</label><select id="ct-pp"></select><label for="ct-pv">variable</label><select id="ct-pv"></select></div>
        <table class="ct-seg"><thead><tr><th>segment</th><th>value</th><th>duration</th><th></th></tr></thead><tbody id="ct-seg"></tbody></table>
        <div class="ct-row"><button type="button" id="ct-seg-add">Add a segment</button><button type="button" class="primary" id="ct-seg-apply">Apply</button><button type="button" id="ct-seg-revert">Return to the applied protocol</button></div>
        <canvas id="ct-proto" class="ct-canvas" role="img" aria-label="Response to the protocol"></canvas>
        <p class="ct-legend" id="ct-proto-legend" aria-live="polite"></p>
      </div>
      <div class="ct-pane" data-pane="history" hidden>
        <p class="ct-note">Every step of the construction, with its source: a verified change of FieldBridge, or a term built here. The construction can be saved as a file and replayed.</p>
        <ol class="ct-history" id="ct-history"></ol>
        <div class="ct-row"><button type="button" id="ct-export">Save the construction</button>
          <label class="ct-file">Replay a saved construction<input type="file" id="ct-import" accept="application/json,.json"></label></div>
        <p class="ct-legend" id="ct-import-msg" aria-live="polite"></p>
      </div>
      <div class="ct-pane" data-pane="examples" hidden>
        <ol class="ct-examples">
          <li><b>The cusp from the pitchfork.</b> The pitchfork normal form dx/dt = εx − x³ has two stable states for ε &gt; 0. A constant field h added to the equation breaks the symmetry x → −x. Over the plane of ε and h the region with two states is bounded by two lines of folds that meet at the origin: the cusp. <button type="button" data-example="cusp">Build it</button></li>
          <li><b>Erasing by annealing under a field.</b> With the field h of the first example, the control ε is lowered below zero and raised again. Whichever state the system held, it ends in the state favoured by the sign of h: the stored state is erased and rewritten. <button type="button" data-example="anneal">Build it</button></li>
          <li><b>A term the carrier refuses.</b> On a ring of phases a term proportional to a phase itself is not a function on the circle; the panel refuses it and states why. <button type="button" data-example="refuse">Show it</button></li>
        </ol>
        <p class="ct-legend" id="ct-example-msg" aria-live="polite"></p>
      </div>`;
    host.insertAdjacentElement('afterend', el);
    el.querySelectorAll('[data-ct]').forEach(b => b.addEventListener('click', () => tab(b.dataset.ct)));
    $('ct-eq').addEventListener('change', () => { selected = null; renderTerms(); });
    $('ct-coef').addEventListener('input', () => { if (selected) preview(selected); });
    $('ct-keep').addEventListener('click', () => {
      if (!selected) return;
      const ok = addTerm(Number($('ct-eq').value), selected, Number($('ct-coef').value));
      if (ok) { selected = null; }
    });
    $('ct-plane-run').addEventListener('click', () => plane());
    $('ct-plane').addEventListener('click', planeClick);
    $('ct-plane').addEventListener('mousemove', planeHover);
    $('ct-pp').addEventListener('change', () => { resetProtocol(); });
    $('ct-pv').addEventListener('change', () => { if (draft) draft.variable = Number($('ct-pv').value); drawProtocol(); });
    $('ct-seg-add').addEventListener('click', () => { const last = draft.segments[draft.segments.length - 1]; draft.segments.push({value: last ? last.value : 0, T: 10}); renderSegments(); drawProtocol(); });
    $('ct-seg-apply').addEventListener('click', () => { committed = clone(draft); logStep('protocol', `protocol of ${draft.param}: ${draft.segments.map(s => `${fmt(s.value, 3)} for ${fmt(s.T, 3)}`).join(', ')}`); drawProtocol(); });
    $('ct-seg-revert').addEventListener('click', () => { if (committed) { draft = clone(committed); renderSegments(); drawProtocol(); } });
    $('ct-export').addEventListener('click', () => {
      const blob = new Blob([JSON.stringify(exportConstruction(), null, 2)], {type: 'application/json'});
      const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = 'fieldbridge-construction.json'; a.click();
      setTimeout(() => URL.revokeObjectURL(a.href), 1000);
    });
    $('ct-import').addEventListener('change', ev => {
      const file = ev.target.files[0]; if (!file) return;
      file.text().then(t => { try { replay(JSON.parse(t)); $('ct-import-msg').textContent = 'Replayed.'; } catch (err) { $('ct-import-msg').textContent = 'Not replayed: ' + err.message + '.'; } });
    });
    el.querySelectorAll('[data-example]').forEach(b => b.addEventListener('click', () => example(b.dataset.example)));
    new ResizeObserver(() => { drawPlane(); drawProtocol(); }).observe(el);
  }
  function tab(key) {
    document.querySelectorAll('#ct [data-ct]').forEach(b => b.setAttribute('aria-selected', String(b.dataset.ct === key)));
    document.querySelectorAll('#ct [data-pane]').forEach(p => { p.hidden = p.dataset.pane !== key; });
    if (key === 'plane' && !grid) plane();
    if (key === 'protocol') drawProtocol();
  }
  /** Called after every move of the constructor: the tools follow the current realization. */
  function onNode(rec) {
    const el = $('ct'); if (!el) return;
    el.classList.toggle('ct-off', !usable(rec));
    if (!usable(rec)) { el.querySelector('.ct-tabs').dataset.off = `The tools act on realizations written as equations; ${strip(rec.name)} is computed by its own engine.`; return; }
    const e = rec.engine;
    $('ct-eq').innerHTML = e.variables.map((v, i) => `<option value="${i}">d${esc(v)}/dt</option>`).join('');
    selected = null; renderTerms();
    const keys = Object.keys(e.ranges || {});
    const opts = keys.map(k => `<option value="${esc(k)}">${esc(k)}</option>`).join('');
    const px = $('ct-px'), py = $('ct-py'), keepX = px.value, keepY = py.value;
    px.innerHTML = opts; py.innerHTML = opts;
    px.value = keys.includes(keepX) ? keepX : keys[0] || '';
    py.value = keys.includes(keepY) && keepY !== px.value ? keepY : keys.find(k => k !== px.value) || '';
    grid = null; drawPlane();
    if (!$('ct').querySelector('[data-pane="plane"]').hidden) plane();
    $('ct-pp').innerHTML = opts; $('ct-pv').innerHTML = e.variables.map((v, i) => `<option value="${i}">${esc(v)}</option>`).join('');
    if (!draft || draft.node !== baseOf(rec) || !keys.includes(draft.param)) resetProtocol(); else { $('ct-pp').value = draft.param; $('ct-pv').value = draft.variable; renderSegments(); drawProtocol(); }
  }
  const baseOf = rec => (rec.local ? rec.local.base.id : rec.id);

  // ------------------------------------------------------------------------------------------- add a term
  let selected = null;
  function renderTerms() {
    const rec = S.nodes[I.state.node]; if (!usable(rec)) return;
    const e = rec.engine, i = Number($('ct-eq').value || 0), params = {...I.state.params, c: 1};
    $('ct-terms').innerHTML = monomials(e, i).map(m => {
      const chk = check(e, i, m, params);
      return `<button type="button" class="ct-term${chk.ok ? '' : ' refused'}" data-key="${esc(m.key)}" aria-pressed="${selected === m.key}" title="${esc(chk.reason)}">${MML.math(MML.tree(withCoef('c', m.tree)))}</button>`;
    }).join('');
    $('ct-terms').querySelectorAll('[data-key]').forEach(b => b.addEventListener('click', () => {
      selected = b.dataset.key;
      $('ct-terms').querySelectorAll('[data-key]').forEach(x => x.setAttribute('aria-pressed', String(x.dataset.key === selected)));
      preview(selected);
    }));
    $('ct-preview').innerHTML = '<span class="muted">Choose a term to see its consequence before it is added.</span>';
    $('ct-keep').disabled = true;
  }
  function preview(key) {
    const rec = S.nodes[I.state.node], e = rec.engine, i = Number($('ct-eq').value), mono = monomials(e, i).find(m => m.key === key);
    let value = Number($('ct-coef').value) || 0;
    const chk = check(e, i, mono, {...I.state.params, c: value});
    if (!chk.ok) {
      $('ct-preview').innerHTML = `<p class="ct-refused"><b>Not allowed on this carrier.</b> ${esc(chk.reason)}</p>`;
      $('ct-keep').disabled = true; return;
    }
    let signNote = '';
    if (chk.sign && Math.sign(value) && Math.sign(value) !== chk.sign) { value = -value; signNote = ` The coefficient is taken as ${fmt(value, 3)}.`; }
    const c = freshName(e), params = {...I.state.params, [c]: value};
    const trial = {...e, params, terms: [...e.terms, {id: 'trial', variable: i, tree: withCoef(c, mono.tree)}]};
    const before = stableStates(e, I.state.params, I.state.removed, I.state.seeds), after = stableStates(trial, params, I.state.removed, before.points);
    const eq = `d${e.variables[i]}/dt gains ${MML.math(MML.tree(withCoef(c, mono.tree)))} with ${esc(c)} = ${fmt(value, 3)}`;
    $('ct-preview').innerHTML = `<p>${eq}. <span class="muted">Allowed: ${esc(chk.reason)}${esc(signNote)}</span></p>
      <dl class="ct-compare"><dt>now</dt><dd>${esc(describe(e, before))}</dd><dt>with the term</dt><dd>${esc(describe(trial, after))}</dd></dl>`;
    $('ct-keep').disabled = false;
  }

  // ------------------------------------------------------------------------------------------- two parameters
  const N = 28;
  let grid = null, job = 0, hover = null;
  function plane() {
    const rec = S.nodes[I.state.node]; if (!usable(rec)) return;
    const e = rec.engine, kx = $('ct-px').value, ky = $('ct-py').value;
    if (!kx || !ky || kx === ky) { grid = null; $('ct-plane-legend').textContent = 'Two different parameters with ranges are needed. A term added in the first tab brings its coefficient as a further parameter.'; drawPlane(); return; }
    const [x0, x1] = e.ranges[kx], [y0, y1] = e.ranges[ky], my = ++job;
    grid = {kx, ky, x: [Math.min(x0, I.state.params[kx]), Math.max(x1, I.state.params[kx])], y: [Math.min(y0, I.state.params[ky]), Math.max(y1, I.state.params[ky])], n: new Array(N * N).fill(null), node: rec.id};
    let k = 0, seeds = I.state.seeds || e.states || [];
    const chunk = () => {
      if (my !== job) return;
      const t0 = performance.now();
      while (k < N * N && performance.now() - t0 < 60) {
        const r = Math.floor(k / N), col = r % 2 ? N - 1 - (k % N) : k % N;   // boustrophedon: neighbouring cells pass their states on
        const params = {...I.state.params, [kx]: grid.x[0] + (col + 0.5) * (grid.x[1] - grid.x[0]) / N, [ky]: grid.y[0] + (r + 0.5) * (grid.y[1] - grid.y[0]) / N};
        const pts = D.equilibria(e, params, I.state.removed, seeds, e.variables.length > 3 ? 4 : 6, 3 + k).filter(p => p.stable);
        grid.n[r * N + col] = pts.length; if (pts.length) seeds = pts.map(p => p.q);
        k++;
      }
      drawPlane();
      if (k < N * N) setTimeout(chunk, 0);
      else $('ct-plane-legend').innerHTML = legendPlane();
    };
    $('ct-plane-legend').textContent = 'Computing…';
    setTimeout(chunk, 0);
  }
  const shade = n => (n == null ? null : n === 0 ? 'none' : Math.min(n, 4));
  function legendPlane() {
    const counts = {}; grid.n.forEach(n => { const s = shade(n); counts[s] = (counts[s] || 0) + 1; });
    const parts = [['none', 'no stable state of rest'], [1, 'one'], [2, 'two'], [3, 'three'], [4, 'four or more']].filter(([s]) => counts[s]);
    return parts.map(([s, t]) => `<span class="ct-key ct-key-${s}"></span>${t} (${Math.round(100 * counts[s] / grid.n.length)} %)`).join(' ') + ' · the ring marks the values set now';
  }
  function canvasSize(cv) {
    const r = cv.getBoundingClientRect(), dpr = window.devicePixelRatio || 1, w = Math.max(200, r.width), h = Math.max(160, r.height);
    cv.width = w * dpr; cv.height = h * dpr; const g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0); return {g, w, h};
  }
  const PAD = {l: 46, r: 12, t: 10, b: 32};
  function drawPlane() {
    const cv = $('ct-plane'); if (!cv || cv.offsetParent === null) return;
    const {g, w, h} = canvasSize(cv); g.clearRect(0, 0, w, h);
    if (!grid) return;
    const W = w - PAD.l - PAD.r, H = h - PAD.t - PAD.b, cw = W / N, ch = H / N;
    const ramp = {1: css('--panel-2'), 2: mix(css('--accent'), 0.45), 3: mix(css('--accent'), 0.75), 4: css('--accent')};
    for (let r = 0; r < N; r++) for (let c = 0; c < N; c++) {
      const s = shade(grid.n[r * N + c]); if (s == null) continue;
      const x = PAD.l + c * cw, y = PAD.t + (N - 1 - r) * ch;
      if (s === 'none') { g.fillStyle = css('--grid'); g.fillRect(x, y, cw + 0.5, ch + 0.5); g.strokeStyle = css('--axis'); g.lineWidth = 0.6; g.beginPath(); g.moveTo(x, y + ch); g.lineTo(x + cw, y); g.stroke(); }
      else { g.fillStyle = ramp[s]; g.fillRect(x, y, cw + 0.5, ch + 0.5); }
    }
    axes(g, PAD.l, PAD.t, W, H, grid.x, grid.y, grid.kx, grid.ky);
    const px = PAD.l + (I.state.params[grid.kx] - grid.x[0]) / (grid.x[1] - grid.x[0]) * W, py = PAD.t + H - (I.state.params[grid.ky] - grid.y[0]) / (grid.y[1] - grid.y[0]) * H;
    g.strokeStyle = css('--ink'); g.lineWidth = 2; g.beginPath(); g.arc(px, py, 6, 0, 2 * Math.PI); g.stroke();
    if (hover) { g.fillStyle = css('--ink'); g.font = '12px ' + css('--sans'); g.fillText(hover, PAD.l + 6, PAD.t + 14); }
  }
  function mix(hex, a) {
    const p = h => { const m = h.replace('#', ''); return [0, 2, 4].map(i => parseInt(m.length === 3 ? m[i / 2] + m[i / 2] : m.slice(i, i + 2), 16)); };
    try { const [r, g2, b] = p(hex), [R, G, B] = p(css('--panel')); return `rgb(${Math.round(R + a * (r - R))},${Math.round(G + a * (g2 - G))},${Math.round(B + a * (b - B))})`; } catch (err) { return hex; }
  }
  function axes(g, x, y, W, H, xr, yr, xl, yl) {
    g.strokeStyle = css('--axis'); g.lineWidth = 1; g.strokeRect(x, y, W, H);
    g.fillStyle = css('--muted'); g.font = '12px ' + css('--sans'); g.textAlign = 'center';
    [0, 0.5, 1].forEach(f => { g.fillText(fmt(xr[0] + f * (xr[1] - xr[0]), 3), x + f * W, y + H + 14); });
    g.fillText(xl, x + W / 2, y + H + 28);
    g.textAlign = 'right';
    [0, 1].forEach(f => { g.fillText(fmt(yr[0] + f * (yr[1] - yr[0]), 2), x - 5, y + H - f * H + 4); });
    g.save(); g.translate(14, y + H / 2); g.rotate(-Math.PI / 2); g.textAlign = 'center'; g.fillText(yl, 0, 0); g.restore();
    g.textAlign = 'left';
  }
  function planeAt(ev) {
    const cv = $('ct-plane'), r = cv.getBoundingClientRect(), W = r.width - PAD.l - PAD.r, H = r.height - PAD.t - PAD.b;
    const fx = (ev.clientX - r.left - PAD.l) / W, fy = 1 - (ev.clientY - r.top - PAD.t) / H;
    if (!grid || fx < 0 || fx > 1 || fy < 0 || fy > 1) return null;
    return {vx: grid.x[0] + fx * (grid.x[1] - grid.x[0]), vy: grid.y[0] + fy * (grid.y[1] - grid.y[0]), cell: Math.min(N - 1, Math.floor(fy * N)) * N + Math.min(N - 1, Math.floor(fx * N))};
  }
  function planeHover(ev) {
    const a = planeAt(ev); const n = a ? grid.n[a.cell] : null;
    hover = a ? `${grid.kx} = ${fmt(a.vx, 3)}, ${grid.ky} = ${fmt(a.vy, 3)}: ${n == null ? '…' : n === 0 ? 'no stable state of rest' : n + ' stable state' + (n === 1 ? '' : 's')}` : null;
    drawPlane();
  }
  function planeClick(ev) {
    const a = planeAt(ev); if (!a) return;
    I.state.params[grid.kx] = a.vx; I.state.params[grid.ky] = a.vy;
    logStep('plane', `${grid.kx} = ${fmt(a.vx, 3)}, ${grid.ky} = ${fmt(a.vy, 3)} set on the plane`);
    I.render(); drawPlane(); drawProtocol();
  }

  // ------------------------------------------------------------------------------------------- protocol
  let draft = null, committed = null;
  const clone = p => (p ? {...p, segments: p.segments.map(s => ({...s}))} : null);
  function resetProtocol() {
    const rec = S.nodes[I.state.node]; if (!usable(rec)) return;
    const e = rec.engine, keys = Object.keys(e.ranges || {});
    const param = keys.includes($('ct-pp').value) ? $('ct-pp').value : (e.control && keys.includes(e.control.name) ? e.control.name : keys[0]);
    if (!param) { draft = committed = null; return; }
    const [lo, hi] = e.ranges[param], v = I.state.params[param];
    draft = {node: baseOf(rec), param, variable: 0, segments: [{value: v, T: 10}, {value: Math.abs(v - lo) > Math.abs(v - hi) ? lo : hi, T: 15}, {value: v, T: 25}]};
    committed = null;
    $('ct-pp').value = param; $('ct-pv').value = 0;
    renderSegments(); drawProtocol();
  }
  function renderSegments() {
    if (!draft) return;
    $('ct-seg').innerHTML = draft.segments.map((s, k) => `<tr><td>${k + 1}</td>
      <td><input type="number" step="0.05" value="${s.value}" data-seg="${k}" data-f="value" aria-label="value in segment ${k + 1}"></td>
      <td><input type="number" step="1" min="0.5" value="${s.T}" data-seg="${k}" data-f="T" aria-label="duration of segment ${k + 1}"></td>
      <td>${draft.segments.length > 1 ? `<button type="button" data-del="${k}" aria-label="Remove segment ${k + 1}">×</button>` : ''}</td></tr>`).join('');
    $('ct-seg').querySelectorAll('input').forEach(inp => inp.addEventListener('input', () => {
      const v = Number(inp.value); if (!Number.isFinite(v)) return;
      draft.segments[Number(inp.dataset.seg)][inp.dataset.f] = inp.dataset.f === 'T' ? Math.max(0.5, Math.min(200, v)) : v; drawProtocol();
    }));
    $('ct-seg').querySelectorAll('[data-del]').forEach(b => b.addEventListener('click', () => { draft.segments.splice(Number(b.dataset.del), 1); renderSegments(); drawProtocol(); }));
  }
  /** RK4 through the segments from q0; the carrier is respected as in model-physics.integrate. */
  function runProtocol(e, prot, q0) {
    const dt = Math.min(e.dt || 0.01, 0.02), out = [[0, ...q0]], removed = I.state.removed;
    let q = q0.slice(), t = 0, steps = 0;
    const project = a => (e.carrier.kind === 'orthant' ? a.map(x => Math.max(0, x)) : a);
    for (const s of prot.segments) {
      const f = D.field(e, {...I.state.params, [prot.param]: s.value}, removed), F = a => f(project(a)), end = t + s.T;
      const every = Math.max(1, Math.round(0.1 / dt));
      while (t < end - 1e-9 && steps < 80000) {
        const h = Math.min(dt, end - t), a = F(q), b = F(q.map((x, i) => x + h * a[i] / 2)), c = F(q.map((x, i) => x + h * b[i] / 2)), d = F(q.map((x, i) => x + h * c[i]));
        q = project(q.map((x, i) => x + h * (a[i] + 2 * b[i] + 2 * c[i] + d[i]) / 6)); t += h; steps++;
        if (q.some(x => !Number.isFinite(x) || Math.abs(x) > 1e5)) return {data: out, final: q, status: 'unbounded'};
        if (steps % every === 0) out.push([t, ...q]);
      }
    }
    out.push([t, ...q]);
    return {data: out, final: q, status: 'complete'};
  }
  function drawProtocol() {
    const cv = $('ct-proto'); if (!cv || cv.offsetParent === null || !draft) return;
    const rec = S.nodes[I.state.node]; if (!usable(rec)) return;
    const e = rec.engine, vi = Number($('ct-pv').value || 0);
    const starts = stableStates(e, I.state.params, I.state.removed, I.state.seeds).points.slice(0, 3);
    if (!starts.length) starts.push(e.initial.slice());
    const finalParams = {...I.state.params, [draft.param]: draft.segments[draft.segments.length - 1].value};
    const ends = stableStates(e, finalParams, I.state.removed, starts).points;
    const runs = p => starts.map(q0 => runProtocol(e, p, q0));
    const rd = runs(draft), rc = committed && committed.param === draft.param ? runs(committed) : null;
    const {g, w, h} = canvasSize(cv); g.clearRect(0, 0, w, h);
    const top = 46, W = w - PAD.l - PAD.r, H = h - PAD.t - PAD.b - top;
    const Tt = Math.max(...[draft, rc ? committed : null].filter(Boolean).map(p => p.segments.reduce((a, s) => a + s.T, 0)));
    const vals = [...rd, ...(rc || [])].flatMap(r => r.data.map(row => row[1 + vi])), [ylo, yhi] = V && V.autoRange ? V.autoRange(vals, 0.08) : [Math.min(...vals), Math.max(...vals)];
    const X = t => PAD.l + t / Tt * W, Y = v => PAD.t + top + H - (v - ylo) / ((yhi - ylo) || 1) * H;
    // the protocol: the value of the parameter in each segment, above the response
    const pv = [...draft.segments, ...(committed ? committed.segments : [])].map(s => s.value), plo = Math.min(...pv), phi = Math.max(...pv);
    const PY = v => PAD.t + top - 8 - (phi > plo ? (v - plo) / (phi - plo) : 0.5) * (top - 26);
    const stairs = (p, dash, color) => { g.setLineDash(dash); g.strokeStyle = color; g.lineWidth = 2; g.beginPath(); let t = 0; p.segments.forEach((s, k) => { const y = PY(s.value); if (k) g.lineTo(X(t), y); else g.moveTo(X(0), y); t += s.T; g.lineTo(X(t), y); }); g.stroke(); g.setLineDash([]); };
    if (committed) stairs(committed, [], css('--p'));
    stairs(draft, [5, 4], css('--p'));
    g.fillStyle = css('--muted'); g.font = '12px ' + css('--sans'); g.fillText(`${draft.param}: ${fmt(plo, 3)} to ${fmt(phi, 3)}`, PAD.l + 4, PAD.t + 10);
    axes(g, PAD.l, PAD.t + top, W, H, [0, Tt], [ylo, yhi], 'time', e.variables[vi]);
    const colors = [css('--accent'), css('--r'), css('--c')];
    const curve = (r, k, dash) => { g.setLineDash(dash); g.strokeStyle = colors[k % 3]; g.lineWidth = 2; g.beginPath(); r.data.forEach((row, j) => { const x = X(row[0]), y = Y(row[1 + vi]); if (j) g.lineTo(x, y); else g.moveTo(x, y); }); g.stroke(); g.setLineDash([]); };
    if (rc) rc.forEach((r, k) => curve(r, k, []));
    rd.forEach((r, k) => curve(r, k, [5, 4]));
    const near = q => { let best = -1, bd = Infinity; ends.forEach((s, k) => { const d = P.distance(e, q, s); if (d < bd) { bd = d; best = k; } }); return bd < 0.05 * (e.carrier.scale || 1) + 1e-3 ? best : -1; };
    const label = (r, k) => r.status !== 'complete' ? 'runs away' : near(r.final) >= 0 ? `ends at (${ends[near(r.final)].map(x => fmt(x, 3)).join(', ')}), ${ends.length === 1 ? 'the only stable state' : `one of ${ends.length} stable states`} at the final value` : 'still moving at the end';
    $('ct-proto-legend').innerHTML = (rc ? 'Solid: the applied protocol. ' : '') + 'Dashed: the protocol being edited. ' +
      starts.map((q, k) => `<span style="color:${colors[k % 3]}">●</span> from (${q.map(x => fmt(x, 3)).join(', ')}): ${label(rd[k], k)}${rc ? `; applied: ${label(rc[k], k)}` : ''}`).join('; ') + '.';
  }

  // ------------------------------------------------------------------------------------------- examples
  function example(key) {
    const msg = $('ct-example-msg');
    if (key === 'cusp' || key === 'anneal') {
      I.start('pitchfork');
      addTerm(0, '1', 0.1);
      $('ct-px').value = 'eps'; $('ct-py').value = 'c1';
      if (key === 'cusp') {
        I.state.params.c1 = 0.1; I.render();
        tab('plane'); plane();
        msg.textContent = 'dx/dt = εx − x³ + c1 with c1 = 0.1. On the plane of ε and c1 the shaded wedge with two stable states opens from the origin; its edges are the folds 27 c1² = 4 ε³. Click inside and outside the wedge and watch the response below.';
      } else {
        I.state.params.eps = 1; I.state.params.c1 = 0.05; I.render();
        $('ct-pp').value = 'eps'; resetProtocol();
        draft.segments = [{value: 1, T: 10}, {value: -0.5, T: 15}, {value: 1, T: 30}]; renderSegments();
        tab('protocol'); drawProtocol();
        msg.textContent = 'dx/dt = εx − x³ + c1 with c1 = 0.05. Both stored states are prepared; ε is lowered to −0.5, where only one state exists, and raised again. Both preparations end in the state with the sign of c1. Change the sign of c1 in A and the protocol writes the other state.';
      }
      logStep('replay', key === 'cusp' ? 'example: the cusp from the pitchfork' : 'example: erasing by annealing under a field');
      return;
    }
    const ring = Object.values(S.nodes).find(r => usable(r) && r.engine.carrier.kind === 'torus');
    if (!ring) { msg.textContent = 'No realization on a torus is listed.'; return; }
    I.start(ring.id);
    {
      tab('add'); $('ct-eq').value = 0; renderTerms(); selected = 'x0';
      $('ct-terms').querySelectorAll('[data-key]').forEach(x => x.setAttribute('aria-pressed', String(x.dataset.key === 'x0')));
      preview('x0');
      msg.textContent = `The ${strip(ring.name)} is written on a torus. The linear term c·${ring.engine.variables[0]} is refused, and the first tab states why.`;
    }
  }

  // ------------------------------------------------------------------------------------------- hover preview of verified changes
  let tip = null;
  function equations(rec) {
    const e = rec.engine;
    if (!usable(rec)) return esc(strip(rec.slots && rec.slots.Omega || ''));
    return e.variables.map((v, i) => `<div>${MML.math(`<mfrac><mrow><mi>d</mi>${MML.name(v)}</mrow><mrow><mi>d</mi><mi>t</mi></mrow></mfrac><mo>=</mo>`)} ${
      e.terms.filter(t => t.variable === i).map((t, j) => { const neg = t.tree[0] === 'neg'; return (neg ? ' − ' : j ? ' + ' : '') + MML.math(MML.tree(neg ? t.tree[1] : t.tree)); }).join('') || '0'}</div>`).join('');
  }
  function wireCards() {
    document.querySelectorAll('#slots [data-edge]').forEach(b => {
      if (b.dataset.ctWired) return; b.dataset.ctWired = '1';
      const show = () => {
        const e = I.byId[b.dataset.edge], target = S.nodes[b.dataset.reverse === '1' ? e.from : e.to]; if (!target) return;
        if (!tip) { tip = document.createElement('div'); tip.className = 'ct-tip'; tip.setAttribute('role', 'tooltip'); document.body.appendChild(tip); }
        tip.innerHTML = `<p class="card-label">After this change</p><p><b>${target.name}</b> · ${esc(S.classes[target.class] || target.class)}</p>${equations(target)}<p class="muted">${esc(e.check || '')}</p>`;
        const r = b.getBoundingClientRect(); tip.style.left = Math.max(8, Math.min(window.innerWidth - 380, r.left)) + 'px'; tip.style.top = (window.scrollY + r.bottom + 6) + 'px'; tip.hidden = false;
      };
      const hide = () => { if (tip) tip.hidden = true; };
      b.addEventListener('mouseenter', show); b.addEventListener('focus', show);
      b.addEventListener('mouseleave', hide); b.addEventListener('blur', hide); b.addEventListener('click', hide);
    });
  }

  function init() {
    panel();
    I.onChange((st, info) => { observe(st, info); if (!info) wireCards(); });
    if (I.state.node) { observe(I.state); wireCards(); }
  }
  window.FieldBridgeConstructor = {monomials, check, addTerm, exportConstruction, replay, stableStates, runProtocol, history,
                                   get draft() { return draft; }, get committed() { return committed; }, get grid() { return grid; }};
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
