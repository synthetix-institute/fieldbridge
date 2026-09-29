/* The instrument: the realization expression of the current node, its six components, the verified changes that
 * leave it, and the consequences calculated in the browser for the values the visitor sets. */
(function () {
  'use strict';
  const S = window.FIELDBRIDGE_SITE;
  if (!S) return;
  const U = window.FieldBridgeUnitary, D = window.FieldBridgeDissipative, F = window.FieldBridgeFields;
  const MML = window.FieldBridgeMath, V = window.FieldBridgeViews, M = window.FieldBridgeMechanisms;
  const $ = id => document.getElementById(id);
  const SLOTS = ['Omega', 'Xi', 'C', 'R', 'P', 'A'];
  const reduced = window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;
  const fmt = V.fmt;
  const esc = s => String(s).replace(/[&<>"']/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[c]));
  const repo = p => (S.repo || 'https://github.com/synthetix-institute/fieldbridge') + '/blob/main/' + p;
  const out = {}, into = {};
  S.edges.forEach(e => { (out[e.from] ||= []).push(e); (into[e.to] ||= []).push(e); });
  const byId = Object.fromEntries(S.edges.map(e => [e.id, e]));
  const LETTERS = {
    unitary: {S: ['C', 'sector: a conserved operator restricts the carrier'], A: ['Omega', 'algebra: commutators of H and R until the span closes'],
              K: ['Omega', 'canonical form: [J_a, J_b] = i ε_abc J_c'], L: ['P', 'law: exact evolution against the Rabi law']},
    dissipative: {S: ['Xi', 'symmetry of the drift'], C: ['A', 'continuation of a state along the control'], W: ['P', 'write field toward another state'],
                  R: ['Omega', 'reduction to the critical direction'], U: ['A', 'unfolding to the cusp'], K: ['Xi', 'canonical form of the reduced drift'],
                  L: ['P', 'law, tested by simulating the realization']}};

  const state = {node: null, params: {}, active: [], removed: [], prep: 'top', protocol: 'relax', bias: 0.2, pulse: 1.5,
                 noise: 0, drive: {on: false, amp: 5, nu: 0.5}, path: [], seen: new Set(), view: {az: -0.6, el: 0.32},
                 side: null, playing: !reduced, message: null, sequenceText: null};
  const listeners = [];
  let scene = null, frame = null, last = 0, analysis = null, pending = false;

  // ---------------------------------------------------------------------------------------- navigation
  function load(id, via) {
    const rec = S.nodes[id];
    if (!rec) return;
    state.node = id;
    const e = rec.engine;
    state.params = rec.family === 'field' ? Object.fromEntries(['M', 'T', ...(e.conserved ? [] : ['kappa0', 'Dgrad'])].map(k => [k, e[k]]))
      : {...(e.params || {})};
    state.active = rec.family === 'unitary' ? e.terms.map(() => true) : [];
    state.removed = [];
    state.prep = 'top'; state.protocol = 'relax';
    state.noise = rec.family === 'dissipative' ? (e.noise || 0) : 0;
    state.drive = {on: false, amp: 5, nu: 0.5};
    state.side = null;
    state.message = via || null;
    state.seen.add(rec.class);
    render(via && via.slot);
    listeners.forEach(f => f(state));
  }
  function applyEdge(edge, reverse, text) {
    const to = reverse ? edge.from : edge.to;
    if ((reverse ? edge.to : edge.from) !== state.node) return false;
    state.path.push({edge: edge.id, reverse: !!reverse, node: to});
    load(to, {edge, reverse: !!reverse, slot: edge.slot, text});
    return true;
  }
  function start(id, text) {
    state.path = [{node: id}];
    state.seen = new Set();
    load(id, text ? {text} : null);
  }
  function undo() {
    if (state.path.length < 2) return;
    state.path.pop();
    const prev = state.path[state.path.length - 1];
    load(prev.node, {text: null});
  }
  /** The shortest chain of edges (either direction) from one node to another. */
  function route(from, to) {
    const prev = {[from]: null}, queue = [from];
    while (queue.length) {
      const n = queue.shift();
      if (n === to) break;
      for (const e of out[n] || []) if (!(e.to in prev)) { prev[e.to] = {edge: e, reverse: false, from: n}; queue.push(e.to); }
      for (const e of into[n] || []) if (!(e.from in prev)) { prev[e.from] = {edge: e, reverse: true, from: n}; queue.push(e.from); }
    }
    if (!(to in prev)) return null;
    const steps = [];
    for (let n = to; prev[n]; n = prev[n].from) steps.unshift(prev[n]);
    return steps;
  }
  function walkTo(id) {
    let steps = route(state.node, id);
    if (!steps) { start(id); return; }
    const go = () => {
      const s = steps.shift();
      if (!s) return;
      applyEdge(s.edge, s.reverse);
      if (steps.length) setTimeout(go, reduced ? 0 : 650);
    };
    go();
  }

  // ---------------------------------------------------------------------------------------- the expression
  function render(changedSlot) {
    const rec = S.nodes[state.node], mech = (S.mechanisms || {})[rec.class] || {};
    $('r-name').innerHTML = rec.name;
    $('r-field').textContent = rec.universal ? 'canonical form, no field' : rec.field;
    const src = rec.spec ? `<a href="${repo(rec.spec)}" target="_blank" rel="noopener">specification</a>` :
      `derived from <a href="${repo(rec.base_spec)}" target="_blank" rel="noopener">${esc(rec.base_spec.split('/').pop())}</a>`;
    $('r-source').innerHTML = src + (rec.tutorial ? ` · <a href="${repo(rec.tutorial)}" target="_blank" rel="noopener">tutorial</a>` : '');
    $('r-class').innerHTML = M.glyph(rec.class, 'small') + esc(S.classes_short[rec.class] || rec.class);
    const head = document.querySelector('.mechanism-head');
    if (head && state.shownClass && state.shownClass !== rec.class && !reduced) { head.classList.remove('flash'); void head.offsetWidth; head.classList.add('flash'); }
    state.shownClass = rec.class;
    $('m-glyph').innerHTML = M.glyph(rec.class, 'big');
    $('m-name').textContent = S.classes[rec.class] || rec.class;
    $('m-canonical').innerHTML = mech.canonical ? 'canonical form ' + mech.canonical : '';
    $('m-where').innerHTML = rec.universal ? 'Written without a field.' : `Realized in ${esc(rec.field)}: ${rec.name}.`;
    $('m-bar').innerHTML = M.glyph(rec.class, 'small') + `<b>${esc(S.classes_short[rec.class] || rec.class)}</b><span>${rec.universal ? 'canonical form' : esc(rec.field)}</span>`;
    $('slots').innerHTML = SLOTS.map(slot => `<li class="slot" data-slot="${slot}" id="slot-${slot}">
        <button class="slot-symbol" type="button" data-light="${slot}" aria-label="${S.slots[slot].name}">${S.slots[slot].symbol}</button>
        <span class="slot-name">${S.slots[slot].name}</span>
        <div class="slot-value">${slotValue(rec, slot)}</div>
        <div class="slot-edits">${edits(rec, slot)}</div></li>`).join('');
    if (changedSlot && !reduced) { const row = $('slot-' + changedSlot); if (row) { row.classList.remove('changed'); void row.offsetWidth; row.classList.add('changed'); } }
    wire(rec);
    buildScene(rec);
    recompute();
    renderPath();
  }
  function termChip(label, k, off, slot = 'Omega') {
    return `<button type="button" class="term${off ? ' off' : ''}" data-term="${k}" aria-pressed="${!off}" title="${off ? 'Switch this term on' : 'Switch this term off'}">${label}</button>`;
  }
  function slotValue(rec, slot) {
    const e = rec.engine, fam = rec.family;
    if (fam === 'unitary') {
      if (slot === 'Omega') return `<div class="equation-line"><span class="lhs">${MML.math('<mi>H</mi><mo>=</mo>')}</span>${e.terms.map((t, k) => termChip(MML.math(MML.term(t.tree, t.operator, t.hc)), k, !state.active[k])).join('<span>+</span>')}</div>`;
      if (slot === 'Xi') return V.carrier(rec) + (e.carrier.sector ? `<span class="muted">${esc(rec.slots.Xi)}</span>` : '');
      if (slot === 'C') return `<span>${esc(rec.slots.C)}</span>`;
      if (slot === 'R') return MML.math(e.observable.terms.map((t, i) => {
        const c = String(t.coefficient).trim(), neg = c.startsWith('-');
        return (neg ? '<mo>−</mo>' : i ? '<mo>+</mo>' : '') + MML.term(MML.coefficient(neg ? c.slice(1) : c), t.operator, false);
      }).join(''));
      if (slot === 'P') {
        const preps = [{id: 'top', label: 'top eigenstate of R'}, ...(e.preparations || []).filter(p => p.state)];
        return `<span class="muted">prepared at t = 0; H acts from t = 0</span><div class="choice">${preps.map(p => `<button type="button" data-prep="${p.id}" aria-pressed="${state.prep === p.id}">${p.label}</button>`).join('')}</div>`;
      }
      if (slot === 'A') return params(e.ranges, state.params) || '<span class="muted">numerical coefficients: ' + esc(rec.slots.A) + '</span>';
    }
    if (fam === 'dissipative') {
      if (slot === 'Omega') return drift(rec);
      if (slot === 'Xi') return V.carrier(rec) + `<span class="muted">${esc(rec.slots.Xi)}</span>`;
      if (slot === 'C') return `<span class="muted">${esc(rec.slots.C)}</span>` + params({D: [0, Math.max(0.2, 3 * (e.noise || 0)), 0.001]}, {D: state.noise}, 'noise', {D: 'bath D'});
      if (slot === 'R') return `<span>${esc(rec.slots.R)}</span>`;
      if (slot === 'P') return protocolEditor(rec);
      if (slot === 'A') return params(e.ranges, state.params);
    }
    if (fam === 'field') {
      if (slot === 'A') return params(fieldRanges(e), state.params);
      return `<span>${esc(rec.slots[slot])}</span>`;
    }
    if (fam === 'stochastic') {
      if (slot === 'A') return params(e.ranges, state.params);
      if (slot === 'Omega') return MML.math('<mi>d</mi><mi>X</mi><mo>=</mo><mi>μ</mi><mi>X</mi><mi>d</mi><mi>t</mi><mo>+</mo><mi>σ</mi><mi>X</mi><mi>d</mi><mi>W</mi>');
      return `<span>${esc(rec.slots[slot])}</span>`;
    }
    return esc(rec.slots[slot] || '');
  }
  function fieldRanges(e) {
    const r = {M: [0.2, 3, 0.01], T: [0.1, 2, 0.01]};
    if (!e.conserved) { r.kappa0 = [0.02, 1, 0.01]; r.Dgrad = [0.1, 3, 0.01]; }
    return r;
  }
  function params(ranges, values, group = 'param', labels = {}) {
    const keys = Object.keys(ranges || {});
    if (!keys.length) return '';
    return `<div class="params">${keys.map(k => {
      const [lo, hi, st] = ranges[k], v = values[k];
      return `<label for="${group}-${k}">${labels[k] ? esc(labels[k]) : MML.math(MML.name(k))}</label>
        <input id="${group}-${k}" type="range" min="${Math.min(lo, v)}" max="${Math.max(hi, v)}" step="${st === 1 ? 1 : 'any'}" value="${v}" data-${group}="${k}">
        <output id="${group}-${k}-out">${fmt(v, 4)}</output>`;
    }).join('')}</div>`;
  }
  function drift(rec) {
    const e = rec.engine;
    if (e.kind_dynamics === 'equations') {
      return e.variables.map((v, i) => `<div class="equation-line"><span class="lhs">${MML.math(`<mfrac><mrow><mi>d</mi>${MML.name(v)}</mrow><mrow><mi>d</mi><mi>t</mi></mrow></mfrac><mo>=</mo>`)}</span>${
        e.terms.filter(t => t.variable === i).map((t, j) => {
          const neg = t.tree[0] === 'neg', body = neg ? t.tree[1] : t.tree;
          return `${j ? `<span>${neg ? '−' : '+'}</span>` : neg ? '<span>−</span>' : ''}${termChip(MML.math(MML.tree(body)), t.id, state.removed.includes(t.id))}`;
        }).join('')}</div>`).join('');
    }
    if (e.kind_dynamics === 'gene') return e.terms.map(t => termChip(esc(t.label), t.id, state.removed.includes(t.id))).join(' ');
    return e.terms.map(t => termChip(esc(t.label) + ` <span class="muted">${esc(t.text)}</span>`, t.id, state.removed.includes(t.id))).join(' ');
  }
  const GREEK = {eps: 'ε', alpha: 'α', gamma: 'γ', mu: 'μ', lam: 'λ', psi: 'ψ', delta: 'δ', Delta: 'Δ', sigma: 'σ', omega: 'ω', nu: 'ν', kappa: 'κ', beta: 'β'};
  const symbol = name => GREEK[name] || name;
  function protocolEditor(rec) {
    const e = rec.engine, opts = [['relax', 'relaxation of random preparations']];
    const writeKind = (e.events || [])[0] ? e.events[0].kind : '';
    if (e.control && e.events && e.events.length && !rec.facts.scale_control) opts.push(['sweep', `sweep of ${symbol(e.control.name)} through the write point`]);
    if ((e.states || []).length >= 2 || rec.facts.states >= 2) opts.push(['pulse', 'field pulse toward the other state']);
    if (e.drive) opts.push(['drive', `periodic modulation of ${symbol(e.drive.param)}`]);
    let html = `<div class="choice">${opts.map(([id, label]) => `<button type="button" data-protocol="${id}" aria-pressed="${state.protocol === id}">${label}</button>`).join('')}</div>`;
    if (state.protocol === 'sweep') html += params({bias: [-0.5, 0.5, 0.005]}, {bias: state.bias}, 'proto', {bias: 'bias'}) + `<p class="note">The preparations start in the single state and the control crosses the write point${writeKind ? ' (' + esc(writeKind) + ')' : ''}; the bias chooses the state that is written.</p>`;
    if (state.protocol === 'pulse') html += params({pulse: [0, 3, 0.01]}, {pulse: state.pulse}, 'proto', {pulse: 'field h'}) + `<p class="note">All preparations start in one stored state; the field acts for 6 time units.</p>`;
    if (state.protocol === 'drive') html += params({amp: [0, 20, 0.1], nu: [-3, 3, 0.01]}, {amp: state.drive.amp, nu: state.drive.nu}, 'drive', {amp: 'amplitude × ε', nu: 'detuning ν/K'}) + `<p class="note">The drive frequency is ${e.drive.ratio}(ω₀ − νK); FieldBridge finds locking for |ν| &lt; 1 at weak drive.</p>`;
    return html;
  }
  /** The label of an edge read backwards: "a → b" becomes "b → a" (keeping a prefix such as "h: "), an added term
   *  is removed. */
  function inverse(label) {
    const parts = label.split(' → ');
    if (parts.length === 2) {
      const m = parts[0].match(/^(.*?(?::|=)\s*)(.+)$/);
      return m ? `${m[1]}${parts[1]} → ${m[2]}` : `${parts[1]} → ${parts[0]}`;
    }
    if (/^\+\s/.test(label)) return 'remove ' + label.replace(/^\+\s*/, '');
    return 'undo: ' + label;
  }
  function changeButton(e, reverse) {
    const here = S.nodes[state.node], target = S.nodes[reverse ? e.from : e.to], same = target.class === here.class;
    const where = target.universal ? 'canonical form' : esc(target.field);
    const what = same ? (target.universal ? 'same mechanism, canonical form' : `same mechanism in ${where}`)
      : esc(S.classes_short[target.class] || target.class);
    return `<button type="button" class="change${reverse ? ' back' : ''}${same ? ' same' : ''}" data-edge="${e.id}"${reverse ? ' data-reverse="1"' : ''} data-to-class="${target.class}">
        <span class="change-what">${reverse ? inverse(e.change) : e.change}</span>
        <span class="change-to">${M.glyph(target.class, 'small')}<span><b>${what}</b><small>${esc(target.short || target.id)}${same || target.universal ? '' : ' · ' + where}</small></span></span></button>`;
  }
  function edits(rec, slot) {
    const fwd = (out[state.node] || []).filter(e => e.slot === slot).map(e => changeButton(e, false));
    const back = (into[state.node] || []).filter(e => e.slot === slot).map(e => changeButton(e, true));
    return fwd.concat(back).join('');
  }
  function wire(rec) {
    document.querySelectorAll('#slots [data-edge]').forEach(b => {
      b.addEventListener('click', () => {
        window.FieldBridgeSite && window.FieldBridgeSite.leaveSequence();
        M.preview(null);
        applyEdge(byId[b.dataset.edge], b.dataset.reverse === '1');
      });
      b.addEventListener('mouseenter', () => M.preview(b.dataset.toClass));
      b.addEventListener('focus', () => M.preview(b.dataset.toClass));
      b.addEventListener('mouseleave', () => M.preview(null));
      b.addEventListener('blur', () => M.preview(null));
    });
    document.querySelectorAll('#slots [data-term]').forEach(b => b.addEventListener('click', () => {
      if (rec.family === 'unitary') { const k = Number(b.dataset.term); state.active[k] = !state.active[k]; }
      else { const id = b.dataset.term; state.removed = state.removed.includes(id) ? state.removed.filter(x => x !== id) : [...state.removed, id]; }
      b.classList.toggle('off'); b.setAttribute('aria-pressed', String(!b.classList.contains('off')));
      rebuild();
    }));
    document.querySelectorAll('#slots [data-prep]').forEach(b => b.addEventListener('click', () => { setPreparation(b.dataset.prep); }));
    document.querySelectorAll('#slots [data-protocol]').forEach(b => b.addEventListener('click', () => {
      state.protocol = b.dataset.protocol; state.drive.on = state.protocol === 'drive';
      $('slot-P').querySelector('.slot-value').innerHTML = protocolEditor(rec); wire(rec); rebuild();
    }));
    const onInput = (sel, set) => document.querySelectorAll(sel).forEach(inp => inp.addEventListener('input', () => {
      const k = inp.dataset[Object.keys(inp.dataset)[0]], v = Number(inp.value);
      set(k, v); const o = document.getElementById(inp.id + '-out'); if (o) o.textContent = fmt(v, 4);
      rebuild();
    }));
    onInput('#slots [data-param]', (k, v) => { state.params[k] = v; });
    onInput('#slots [data-noise]', (k, v) => { state.noise = v; });
    onInput('#slots [data-proto]', (k, v) => { state[k] = v; });
    onInput('#slots [data-drive]', (k, v) => { state.drive[k] = v; });
    document.querySelectorAll('#slots [data-light]').forEach(b => b.addEventListener('mouseenter', () => light(b.dataset.light)));
    document.querySelectorAll('#slots [data-light]').forEach(b => b.addEventListener('mouseleave', () => light(null)));
  }
  function light(slot) {
    document.querySelectorAll('.formula .sym').forEach(s => s.classList.toggle('lit', s.dataset.slot === slot));
    listeners.forEach(f => f(state, {light: slot}));
  }
  function setPreparation(id) {
    state.prep = id;
    document.querySelectorAll('#slots [data-prep]').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.prep === id)));
    rebuild();
  }
  function rebuild() {
    if (pending) return;
    pending = true;
    requestAnimationFrame(() => { pending = false; buildScene(S.nodes[state.node]); recompute(); });
  }

  // ---------------------------------------------------------------------------------------- consequences
  function recompute() {
    const rec = S.nodes[state.node], f = rec.facts, e = rec.engine;
    const lead = $('consequence-lead'), text = $('consequence-text'), facts = [];
    const live = (k, v) => facts.push([k, v, true]), known = (k, v) => facts.push([k, v, false]);
    let leadText = S.classes[rec.class] || rec.class, body = '';
    const base = rec.family === 'field' ? Object.fromEntries(Object.keys(state.params).map(k => [k, e[k]])) : (e.params || {});
    const changed = JSON.stringify(state.params) !== JSON.stringify(base) || state.removed.length || state.active.some(x => !x);
    if (rec.family === 'unitary') {
      const a = U.analyze(e, state.params, state.active);
      analysis = a;
      leadText = S.classes[a.klass];
      if (a.reached) body = `The Hamiltonian and the observable generate su(2): a three-vector of expectations rotates at |Ω| = ${fmt(a.rate)} about an axis at θ = ${fmt(a.theta, 3)}° to R, and the measured signal follows cos²θ + sin²θ cos |Ω|t.`;
      else if (a.klass === 'conserved') body = 'R commutes with H. Its closure has dimension 1: the measured value is conserved.';
      else body = `The commutators of H and R span an algebra of dimension ${a.dim}, not a rotating three-vector.` + (a.cause != null ? ` Without the term ${MML.math(MML.term(e.terms[a.cause].tree, e.terms[a.cause].operator, e.terms[a.cause].hc))} the algebra is su(2).` : '');
      live('closure of R', a.closure); live('dimension of the algebra', a.dim);
      if (a.reached) { live('rate |Ω|', fmt(a.rate, 4)); live('angle θ', fmt(a.theta, 4) + '°'); }
      known('carrier', `${f.carrier}${f.sector !== f.hilbert ? `; ${f.sector} of ${f.hilbert} states` : ''}`);
      if (f.rep) known('representation', f.rep);
    } else if (rec.family === 'dissipative') {
      const att = scene && scene.attractors;
      if (att) {
        const k = att.points.length;
        leadText = k >= 2 ? `${k} stable states` : k === 1 ? 'one stable state' : att.moving ? 'no stable state: a cycle' : S.classes[rec.class];
        live('stable states', k + (att.moving ? `; ${att.moving} of ${att.n} preparations keep moving` : ''));
      }
      body = mechanismText(rec);
      if (f.write_point != null) known('write point', `${symbol(f.write_param)} = ${fmt(Math.abs(f.write_point) < 5e-7 ? 0 : f.write_point, 4)}, ${f.write_kind}`);
      for (const [p, name] of [['sym', 'symmetric write'], ['thr', 'one-sided write'], ['lock', 'phase locking']]) {
        if (!f[p + '_status']) continue;
        const reached = String(f[p + '_status']).startsWith('reached');
        let v = reached ? `${f[p + '_status']} (${f[p + '_word']})` : `stops: ${f[p + '_obstruction'] || f[p + '_status']}`;
        if (reached && f[p + '_law'] != null) v += p === 'lock' ? `; half-width ${fmt(f[p + '_law'], 4)} ± ${fmt(f[p + '_law_err'], 2)} K` : `; law constant ${fmt(f[p + '_law'], 5)} ± ${fmt(f[p + '_law_err'], 2)}`;
        if (reached && p === 'lock') v += `; ${f.lock_ratio}:1`;
        known(name, v);
      }
      if (f.operating_param) known('operating point', `${symbol(f.operating_param)} = ${fmt(f.operating_value)} (letter C)`);
      known('retention', f.loss);
      if (state.noise > 0 && state.protocol !== 'drive') live('bath', `D = ${fmt(state.noise, 3)}`);
      if (scene && scene.driveResult) live('drive', scene.driveResult);
    } else if (rec.family === 'field') {
      const big = largeLattice({...e, ...state.params}), t1 = big.L * big.L / 200;
      const times = logspace(Math.log10(20), Math.log10(t1), 24), v = F.snr(big, times), slope = F.exponent(times, v, 20, t1);
      body = e.conserved ? `A conserved density keeps the total of the write; its profile spreads, and the signal decays as a power of time.` : 'Without a conservation law every mode relaxes at least at the rate M κ₀: the write is lost exponentially.';
      live(`slope of ln SNR against ln t, t = 20–${fmt(t1, 3)}, ${big.L}${big.d === 2 ? '²' : ''} sites`, fmt(slope, 3));
      known('predicted', e.conserved ? `t^${f.exponent}` : `exponential, rate ${fmt(f.rate)}`);
      if (f.exponent_fit != null) known('exact spectral sum (large lattice)', fmt(f.exponent_fit, 4));
    } else if (rec.family === 'stochastic') {
      body = `The mean of log X grows at ${f.growth} in the ${f.convention} reading of the noise term.`;
      const mu = state.params.mu, s = state.params.sigma;
      live('growth of ⟨log X⟩', fmt(e.convention === 'ito' ? mu - s * s / 2 : mu, 4));
      known('drift of log X (verified)', f.growth);
    }
    lead.innerHTML = leadText;
    const via = state.message, msg = via && (via.text || (via.edge && !via.reverse ? via.edge.text : null));
    text.innerHTML = (msg ? `<span>${msg}</span> ` : '') + (body && (!msg || changed) ? `<span class="muted">${body}</span>` : '')
      + (changed ? `<p class="note">Recalculated in the browser for the values you set; derivation words and law constants are FieldBridge's for the values listed in the specification.</p>` : '');
    $('facts').innerHTML = facts.map(([k, v, isLive]) => `<dt>${esc(k)}</dt><dd class="${isLive ? 'live' : ''}">${v}</dd>`).join('');
    const word = rec.family === 'unitary' ? f.word : [f.sym_word, f.thr_word, f.lock_word].filter(w => w && w !== '—')[0];
    const letters = LETTERS[rec.family === 'unitary' ? 'unitary' : 'dissipative'];
    $('derivation').innerHTML = word ? `<span>Derivation</span>` + word.split('').map(l => { const [slot, meaning] = letters[l] || ['A', l]; return `<span class="letter" data-slot="${slot}" title="${esc(meaning)}">${l}</span>`; }).join('')
      + `<span>· FieldBridge, at the listed values</span>` : '';
    const card = rec.card;
    $('card-figure').hidden = !card;
    if (card) { $('card-image').src = card.image; $('card-image').alt = 'FieldBridge memory card of ' + rec.name.replace(/<[^>]+>/g, ''); }
  }
  function mechanismText(rec) {
    const f = rec.facts;
    switch (rec.class) {
      case 'symmetric-write': return `Two stable states appear at a supercritical pitchfork of ${symbol(f.write_param)}: sweeping the control through the write point with a weak bias writes the state the bias favours.`;
      case 'threshold-write': return 'A stored state disappears at a fold: a field or a control past the threshold switches the realization to the other state, with a delay set by the Airy law.';
      case 'subcritical-write': return 'The stored state loses stability at a subcritical pitchfork: the realization jumps to a distant branch.';
      case 'field-write': return 'The control only rescales the energy; a state is written by a uniform field and kept by the barriers between states.';
      case 'single-state': return 'Every preparation relaxes to one state: nothing of the preparation is kept.';
      case 'oscillation': return 'The preparations settle on a limit cycle. Its phase is a flat direction; a periodic drive can fix it.';
      case 'neutral-cycles': return 'A conserved quantity fills the plane with closed orbits: no orbit attracts its neighbours, and no drive-independent phase is kept.';
      default: return '';
    }
  }
  const logspace = (a, b, n) => Array.from({length: n}, (_, i) => Math.pow(10, a + (b - a) * i / (n - 1)));

  // ---------------------------------------------------------------------------------------- scenes
  function buildScene(rec) {
    const e = rec.engine;
    const titles = {unitary: 'Expectation of the rotating three-vector (left) and the measured signal (right)'};
    if (rec.family === 'unitary') scene = unitaryScene(rec);
    else if (rec.family === 'dissipative') scene = dissipativeScene(rec);
    else if (rec.family === 'field') scene = fieldScene(rec);
    else scene = stochasticScene(rec);
    $('view-title').textContent = scene.title || titles[rec.family] || '';
    $('view-legend').innerHTML = scene.legend || '';
    scene.t = 0;
    draw();
    loop();
  }
  function unitaryScene(rec) {
    const e = rec.engine, H = U.hamiltonian(e, state.params, state.active), O = U.observable(e);
    const prep = (e.preparations || []).find(p => p.id === state.prep && p.state);
    const psi0 = prep ? {re: Float64Array.from(prep.state.re), im: Float64Array.from(prep.state.im)} : U.topEigenstate(O);
    const a = U.analyze(e, state.params, state.active);
    const period = a.reached ? 2 * Math.PI / a.rate : 6, tMax = 2 * period;
    const times = Array.from({length: 361}, (_, i) => tMax * i / 360);
    const rows = U.signal(e, state.params, state.active, psi0, times);
    let omega = null;
    if (e.frame) {
      const Fm = e._frame || (e._frame = e.frame.map(U.fromJSON));
      const H0 = U.traceless(H);
      omega = Fm.map(Fa => U.inner(H0, Fa) / U.inner(Fa, Fa));
    }
    const law = a.reached && state.prep === 'top' ? U.rabi(a.theta, a.rate) : null;
    return {kind: 'unitary', rows, tMax, law, omega, period, speed: tMax / 9,
            legend: `Blue: exact evolution on the carrier. ${law ? 'Dashed: the Rabi law cos²θ + sin²θ cos |Ω|t. ' : ''}${e.frame ? `Axes: the frame of the rotation${e.frame_from !== rec.id ? ' of ' + esc(S.nodes[e.frame_from].name.replace(/<[^>]+>/g, '')) : ''}, R vertical. Drag the sphere to turn it.` : ''}`,
            title: 'The expectation of the rotating three-vector, and the measured signal'};
  }
  function dissipativeScene(rec) {
    const e = rec.engine, dim = e.variables.length, params = {...state.params};
    const f = D.field(e, params, state.removed);
    // stable states by Newton from the states FieldBridge stored (or from those of the previous values), and a few
    // random starts; only when there is none, a short relaxation tells a cycle from a drift
    const seeds = state.seedsFor === rec.id && state.seeds ? state.seeds : (e.states || []);
    const stable = D.equilibria(e, params, state.removed, seeds, dim > 4 ? 6 : 12, 5).filter(p => p.stable);
    const att = stable.length ? {points: stable.map(p => ({q: p.q, eig: p.eig})), moving: 0, unbounded: 0, n: 0}
      : D.attractors(e, params, state.removed, {n: 6, T: 60, seed: 5});
    state.seeds = att.points.map(p => p.q); state.seedsFor = rec.id;
    const sc = {kind: 'dissipative', attractors: att, dim, t: 0, speed: 1, particles: [], trails: [], side: null, record: 3};
    const torus = e.carrier.kind === 'torus', P = e.carrier.period || 2 * Math.PI;
    const scale = e.carrier.scale || 1;
    // the range of the drawing: the stored states and the scale of the carrier
    const pts = att.points.map(p => p.q);
    // one variable: a range set by the stored states, so that the wells and the barrier between them fill the drawing
    const span = i => {
      const v = pts.map(q => q[i]), reach = v.length ? Math.max(...v.map(Math.abs)) : 0;
      const hi = dim === 1 ? Math.max(0.5 * scale, 1.45 * reach) : Math.max(scale, ...v) * 1.25;
      return e.carrier.kind === 'orthant' ? [0, hi] : torus ? [0, P] : [-hi, hi];
    };
    sc.xr = span(0); sc.yr = dim > 1 ? span(1) : null;
    const random = D.rng(17);
    const n = dim === 1 ? 26 : e.kind_dynamics === 'rotor' ? 1 : 70;
    const sample = () => e.variables.map((_, i) => torus ? random() * P : e.carrier.kind === 'orthant' ? random() * sc.xr[1] * (i ? (sc.yr ? sc.yr[1] / sc.xr[1] : 1) : 1) : (2 * random() - 1) * sc.xr[1]);
    const control = e.control && e.control.name;
    let flow = (q, t) => f(q), pStart = null;
    if (state.protocol === 'sweep' && control && att.points.length >= 1) {
      const [lo, hi] = e.control.range, ends = [lo, hi].map(v => D.attractors(e, {...params, [control]: v}, state.removed, {n: 8, T: 40, seed: 3}).points.length);
      pStart = ends[0] <= ends[1] ? lo : hi;
      const target = att.points.length >= 2 ? att.points[state.bias >= 0 ? 0 : 1].q : null, Tsweep = 40;
      flow = (q, t) => {
        const v = t < Tsweep ? pStart + (params[control] - pStart) * t / Tsweep : params[control];
        const base = D.field(e, {...params, [control]: v}, state.removed)(q);
        if (target && t < Tsweep) { const force = writeForce(e, q, target, Math.abs(state.bias)); return base.map((x, i) => x + force[i]); }
        return base;
      };
      const start = D.attractors(e, {...params, [control]: pStart}, state.removed, {n: 6, T: 40, seed: 2}).points[0];
      sc.initial = () => start ? start.q.map(x => x + 1e-3 * (random() - 0.5)) : sample();
      sc.window = [0, Tsweep]; sc.tMax = 70; sc.target = target;
    } else if (state.protocol === 'pulse' && att.points.length >= 2) {
      const a0 = att.points[0].q, target = att.points[1].q;
      flow = (q, t) => { const base = f(q); if (t >= 4 && t < 10) { const force = writeForce(e, q, target, state.pulse); return base.map((x, i) => x + force[i]); } return base; };
      sc.initial = () => a0.map(x => x + 0.02 * (random() - 0.5)); sc.window = [4, 10]; sc.tMax = 40;
    } else if (state.protocol === 'drive' && e.drive) {
      const d = e.drive, K = state.drive.amp * d.K, eps = state.drive.amp * d.eps, wf = d.ratio * (d.omega - state.drive.nu * K);
      const q0 = att.points.length ? att.points[0].q : e.initial;
      const T = Math.min(600, 50 * d.period);
      const data = D.driven(e, params, state.removed, {param: d.param, p1: d.p1, eps, omega: wf}, q0.slice(), T, Math.min(0.02, e.dt || 0.01), 0.05);
      const phases = D.crossingPhases(data, 0, 2 * Math.PI / wf, d.ratio);
      const tail = phases.slice(-12).map(p => p[1]), spread = tail.length > 3 ? circularSpread(tail) : 1;
      sc.drive = {data, phases, T, locked: spread < 0.03};
      sc.driveResult = sc.drive.locked ? `locked (phase spread ${fmt(spread, 2)} cycle over the last crossings)` : `slipping (phase spread ${fmt(spread, 2)} cycle)`;
      sc.tMax = T;
    }
    if (!sc.initial) sc.initial = sample;
    if (dim > 1 && !torus && e.kind_dynamics !== 'rotor') {
      // the drawing also covers where the dynamics goes: the driven trajectory, or a short run from a few preparations
      const seen = sc.drive ? sc.drive.data.map(r => r.slice(1)) : [];
      if (!seen.length) for (let k = 0; k < 3; k++) {
        const run = window.FieldBridgeModels.integrate(e, {params, removed: state.removed, initial: sample(), duration: 40});
        if (run.status === 'complete') run.data.forEach(r => seen.push(r.slice(1)));
      }
      const widen = (r, vals) => { const [lo, hi] = V.autoRange(vals, 0.08); return [Math.min(r[0], lo), Math.max(r[1], hi)]; };
      if (seen.length) { sc.xr = widen(sc.xr, seen.map(q => q[0])); sc.yr = widen(sc.yr, seen.map(q => q[1])); }
      if (e.carrier.kind === 'orthant') { sc.xr[0] = Math.max(0, sc.xr[0]); sc.yr[0] = Math.max(0, sc.yr[0]); }
    }
    sc.flow = flow;
    sc.particles = Array.from({length: n}, () => sc.initial());
    sc.trails = sc.particles.map(q => [q.slice()]);
    sc.tMax = sc.tMax || (dim === 1 ? 30 : 40);
    sc.dt = Math.min(0.05, (e.dt || 0.01) * 5);
    // side: bifurcation diagram along the control for memory realizations
    if (control && !rec.facts.scale_control && rec.class !== 'oscillation' && rec.class !== 'neutral-cycles') {
      const [lo, hi] = e.control.range, values = Array.from({length: 41}, (_, i) => lo + (hi - lo) * i / 40);
      sc.scan = D.scan(e, params, state.removed, control, values, 10);
      sc.range = [lo, hi];
    }
    if (dim === 1) {
      sc.xs = Array.from({length: 241}, (_, i) => (torus ? -Math.PI : sc.xr[0]) + ((torus ? 2 * Math.PI : sc.xr[1] - sc.xr[0]) * i / 240));
      sc.V = D.potential1D(e, params, state.removed, sc.xs);
    }
    sc.title = e.kind_dynamics === 'rotor' ? 'One configuration of the rotors, and the overlap with the stored states' :
      dim === 1 ? (torus ? 'Directions of the preparations, and the energy along the angle' : 'Preparations in the energy landscape, and the states along the control') :
      state.protocol === 'drive' ? 'Trajectory under the drive, and the phase of the oscillation against the drive' : 'Preparations in the plane of the first two variables, and the states along the control';
    sc.legend = legendFor(rec, sc);
    return sc;
  }
  function circularSpread(xs) { const c = xs.reduce((s, x) => s + Math.cos(2 * Math.PI * x), 0) / xs.length, s = xs.reduce((a, x) => a + Math.sin(2 * Math.PI * x), 0) / xs.length; return 1 - Math.hypot(c, s); }
  function writeForce(e, q, target, h) {
    if (e.carrier.kind === 'torus') { const m = 2 * Math.PI / e.carrier.period; return q.map((x, i) => -h * m * Math.sin(m * (x - target[i]))); }
    const d = target.map((x, i) => x - q[i]), n = Math.hypot(...d) || 1;
    return d.map(x => h * x / n);
  }
  function legendFor(rec, sc) {
    const parts = [];
    if (sc.attractors.points.length) parts.push('Green circles: stable states calculated for the current values.');
    if (state.protocol === 'sweep') parts.push('Shaded: the sweep of the control with the bias.');
    if (state.protocol === 'pulse') parts.push('Shaded: the field pulse.');
    if (state.noise > 0 && state.protocol !== 'drive') parts.push(`Noise D = ${fmt(state.noise, 3)}.`);
    if (sc.scan) parts.push('Right: equilibria along the control; the blue line is its current value, dotted lines the write points FieldBridge located.');
    return parts.join(' ');
  }
  /** The lattice of the exponent check of fields.py: large enough that its size does not end the power law. */
  const largeLattice = e => ({...e, L: e.d === 1 ? 4096 : 256});
  function fieldScene(rec) {
    const e = {...rec.engine, ...state.params}, big = largeLattice(e);
    const times = logspace(-1, e.d === 1 ? 4.5 : 2.8, 90), snr = F.snr(big, times);
    const guide = e.conserved ? times.filter(t => t >= 10).map(t => [t, snr[times.findIndex(x => x >= 100)] * Math.pow(t / 100, rec.facts.exponent_value)]) : null;
    return {kind: 'field', e, times, snr, guide, tMax: 3000, t: 0.1, initial: F.profile(e, 0), speed: 1,
            title: 'The trace of a write along the lattice, and the signal it leaves', legend: `Left: the mean trace on the lattice of the specification (${e.L}${e.d === 2 ? '²' : ''} sites); dashed, at t = 0. Right: signal-to-noise ratio of the best measurement of the write on ${big.L}${big.d === 2 ? '²' : ''} sites, where the size of the lattice plays no role${e.conserved ? `; dashed, t^${rec.facts.exponent}` : ''}.`};
  }
  function stochasticScene(rec) {
    const e = rec.engine, mu = state.params.mu, s = state.params.sigma, random = D.rng(9), T = 10, dt = 0.02;
    const drift = e.convention === 'ito' ? mu - s * s / 2 : mu;
    const gauss = () => { let u = 0, v = 0; while (!u) u = random(); while (!v) v = random(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); };
    const paths = Array.from({length: 7}, () => { let y = 0; const pts = [[0, 0]]; for (let t = dt; t <= T + 1e-9; t += dt) { y += drift * dt + s * Math.sqrt(dt) * gauss(); pts.push([t, y]); } return pts; });
    return {kind: 'stochastic', paths, T, tMax: T, speed: T / 6, drift, mu, s,
            title: 'Paths of log X, and its mean in the two conventions', legend: `Right: mean of log X for Itô (μ − σ²/2) and Stratonovich (μ); the thick line is the convention of this realization.`};
  }

  // ---------------------------------------------------------------------------------------- drawing
  function draw() {
    if (!scene) return;
    const rec = S.nodes[state.node], main = $('view-main'), side = $('view-side');
    if (scene.kind === 'unitary') {
      const k = Math.min(scene.rows.length - 1, Math.round(scene.t / scene.tMax * (scene.rows.length - 1)));
      const row = scene.rows[k], trail = scene.rows.slice(Math.max(0, k - 140), k + 1).filter(r => r.m).map(r => r.m);
      V.sphere(main, {m: row.m, trail, omega: scene.omega, view: state.view, caption: row.m ? '' : 'no rotation frame for this realization'});
      V.signal(side, {tMax: scene.tMax, law: scene.law, cursor: scene.t, curves: [{pts: scene.rows.map(r => [r.t, r.f])}], yLabel: 'measured R, relative to its top value'});
      return;
    }
    if (scene.kind === 'dissipative') {
      const e = rec.engine, dim = scene.dim;
      if (e.kind_dynamics === 'rotor') {
        const q = scene.particles[0];
        V.rotors(main, {positions: e.positions || e.variables.map((_, i) => [i, 0]), angles: q, m: e.graph ? e.graph.m : 2, bonds: (e.links || []).map(l => [l.source, l.target]), caption: `t = ${fmt(scene.t, 3)}`});
        const trail = scene.trails[0];
        const overlaps = scene.attractors.points.slice(0, 3).map((p, i) => ({label: 'state ' + (i + 1), pts: trail.map((qq, k) => [k * scene.dt * scene.record, qq.reduce((s, x, j) => s + Math.cos((2 * Math.PI / e.carrier.period) * (x - p.q[j])), 0) / qq.length])}));
        V.series(side, {tMax: scene.tMax, curves: overlaps, yr: [-1.05, 1.05], yLabel: 'overlap with the stored states', window: scene.window});
        return;
      }
      if (dim === 1) {
        const torus = e.carrier.kind === 'torus', wrapA = x => Math.atan2(Math.sin(x), Math.cos(x));
        if (torus) {
          const psi = e.params.psi != null ? state.params.psi : null, h = state.params.h;
          V.circle(main, {angles: scene.particles.map(q => q[0]), states: scene.attractors.points.map(p => p.q[0]), field: psi != null ? Math.PI / 2 - psi : null, h, fieldLabel: 'field', easyAxis: psi != null, caption: `t = ${fmt(scene.t, 3)}`});
        } else {
          V.landscape(main, {xs: scene.xs, V: scene.V, balls: scene.particles.map(q => q[0]), states: scene.attractors.points.map(p => p.q[0]), xLabel: e.variables[0], yLabel: 'V = −∫ drift'});
        }
        if (scene.scan && state.side !== 'energy') V.bifurcation(side, {scan: scene.scan, range: scene.range, current: state.params[e.control.name], param: symbol(e.control.name), events: e.events, yLabel: torus ? 'angle' : e.variables[0]});
        else if (torus) V.landscape(side, {xs: scene.xs, V: scene.V, balls: scene.particles.map(q => wrapA(q[0])), states: scene.attractors.points.map(p => wrapA(p.q[0])), xLabel: 'angle', yLabel: 'energy'});
        else V.series(side, {tMax: scene.tMax, cursor: scene.t, curves: scene.trails.slice(0, 8).map(tr => ({pts: tr.map((q, k) => [k * scene.dt * scene.record, q[0]])})), window: scene.window});
        return;
      }
      const tr = scene.trails.map(t => ({pts: t.slice(-40).map(q => [q[0], q[1]])}));
      if (scene.drive) {
        const shown = scene.drive.data.filter(r => r[0] <= scene.t).slice(-400);
        V.phasePlane(main, {xr: scene.xr, yr: scene.yr, trails: [{pts: shown.map(r => [r[1], r[2]])}], equilibria: [], xLabel: e.variables[0], yLabel: e.variables[1]});
        V.phaseDrift(side, {phases: scene.drive.phases, tMax: scene.drive.T, cursor: scene.t, locked: scene.drive.locked, caption: scene.drive.locked ? 'locked' : 'slipping'});
        return;
      }
      V.phasePlane(main, {xr: scene.xr, yr: scene.yr, trails: tr, equilibria: scene.attractors.points.map(p => Object.assign([p.q[0], p.q[1]], {stable: true})), xLabel: e.variables[0], yLabel: e.variables[1]});
      if (scene.scan) V.bifurcation(side, {scan: scene.scan, range: scene.range, current: state.params[e.control.name], param: symbol(e.control.name), events: e.events, yLabel: `${e.variables[0]} − mean of the others`});
      else V.series(side, {tMax: scene.tMax, cursor: scene.t, curves: e.variables.map((v, i) => ({label: v, pts: scene.trails[0].map((q, k) => [k * scene.dt * scene.record, q[i]])})), window: scene.window});
      return;
    }
    if (scene.kind === 'field') {
      const t = scene.t, phi = F.profile(scene.e, t);
      V.profile(main, {phi, initial: scene.initial, caption: `t = ${fmt(t, 3)}`});
      V.loglog(side, {t: scene.times, v: scene.snr, guide: scene.guide, cursor: Math.max(scene.times[0], t), note: rec.facts.exponent ? `t^${rec.facts.exponent}` : 'exponential'});
      return;
    }
    if (scene.kind === 'stochastic') {
      const cut = scene.t;
      V.series(main, {tMax: scene.T, cursor: cut, curves: scene.paths.map(p => ({pts: p})), yLabel: 'log X'});
      const ts = Array.from({length: 50}, (_, i) => scene.T * i / 49);
      V.series(side, {tMax: scene.T, curves: [
        {label: 'Itô', pts: ts.map(t => [t, (scene.mu - scene.s * scene.s / 2) * t]), dash: rec.engine.convention === 'ito' ? [] : [5, 4]},
        {label: 'Stratonovich', pts: ts.map(t => [t, scene.mu * t]), dash: rec.engine.convention === 'ito' ? [5, 4] : []}], yLabel: '⟨log X⟩'});
    }
  }
  function advance(dt) {
    if (!scene) return;
    if (scene.kind === 'dissipative' && !scene.drive) {
      const e = S.nodes[state.node].engine, amp = Math.sqrt(2 * state.noise * scene.dt), random = scene.random || (scene.random = D.rng(29));
      const gauss = () => { let u = 0, v = 0; while (!u) u = random(); while (!v) v = random(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); };
      const steps = Math.max(1, Math.round(dt * scene.speed * 6 / scene.dt));
      scene.record = 3;
      for (let s = 0; s < steps && scene.t < scene.tMax; s++) {
        const t = scene.t, h = scene.dt;
        scene.particles = scene.particles.map(q => {
          const a = scene.flow(q, t), b = scene.flow(q.map((x, i) => x + h * a[i] / 2), t + h / 2), c = scene.flow(q.map((x, i) => x + h * b[i] / 2), t + h / 2), d = scene.flow(q.map((x, i) => x + h * c[i]), t + h);
          let r = q.map((x, i) => x + h * (a[i] + 2 * b[i] + 2 * c[i] + d[i]) / 6 + (amp ? amp * gauss() : 0));
          if (e.carrier.kind === 'orthant') r = r.map(x => Math.max(0, x));
          if (e.carrier.kind === 'torus') r = r.map(x => ((x % e.carrier.period) + e.carrier.period) % e.carrier.period);
          return r.map(x => Number.isFinite(x) ? x : 0);
        });
        scene.t += h;
        scene.tick = (scene.tick || 0) + 1;
        if (scene.tick % scene.record === 0) scene.trails.forEach((tr, i) => { tr.push(scene.particles[i].slice()); if (tr.length > 400) tr.shift(); });
      }
      if (scene.t >= scene.tMax && !scene.window && state.noise > 0) scene.tMax += 20;
      return;
    }
    scene.t = Math.min(scene.tMax, scene.t + dt * (scene.kind === 'field' ? scene.t * 0.9 + 0.3 : scene.kind === 'dissipative' ? scene.tMax / 12 : scene.speed));
    if (scene.kind === 'unitary' && scene.t >= scene.tMax) scene.t = 0;
  }
  function loop() {
    cancelAnimationFrame(frame);
    last = 0;
    const tick = time => {
      const dt = last ? Math.min(0.05, (time - last) / 1000) : 0; last = time;
      if (state.playing) advance(dt);
      draw();
      if (state.playing) frame = requestAnimationFrame(tick);
    };
    if (reduced) { advance(scene.tMax || 0); if (scene.kind === 'dissipative') for (let k = 0; k < 40; k++) advance(1); draw(); return; }
    frame = requestAnimationFrame(tick);
  }

  // ---------------------------------------------------------------------------------------- path
  function renderPath() {
    const items = state.path.map((p, i) => {
      const r = S.nodes[p.node], name = M.glyph(r.class, 'small') + `<span>${r.name}</span>`, cur = i === state.path.length - 1;
      const step = p.edge ? `<span class="path-step" data-slot="${byId[p.edge].slot}" title="${esc(byId[p.edge].change.replace(/<[^>]+>/g, ''))}">—${S.slots[byId[p.edge].slot].symbol}${p.reverse ? '↩' : ''}→</span>` : '';
      return `<li>${step}<span class="path-node${cur ? ' current' : ''}" title="${esc(S.classes[r.class] || '')}">${name}</span></li>`;
    });
    $('path-list').innerHTML = items.join('');
    const n = state.path.length - 1, classes = [...state.seen].map(c => S.classes[c]).filter(Boolean);
    $('path-count').textContent = `${n} change${n === 1 ? '' : 's'} · ${classes.length} mechanism${classes.length === 1 ? '' : 's'}: ${classes.join('; ')}`;
    $('path-undo').disabled = n === 0;
    const restart = $('path-restart');
    if (restart) restart.hidden = state.node === S.start && n === 0;
  }

  // ---------------------------------------------------------------------------------------- controls
  function setup() {
    $('path-undo').addEventListener('click', () => { window.FieldBridgeSite && window.FieldBridgeSite.leaveSequence(); undo(); });
    const restart = $('path-restart');
    if (restart) restart.addEventListener('click', () => { window.FieldBridgeSite && window.FieldBridgeSite.leaveSequence(); start(S.start); });
    $('view-play').addEventListener('click', () => {
      state.playing = !state.playing;
      $('view-play').textContent = state.playing ? '❚❚' : '▶';
      $('view-play').setAttribute('aria-label', state.playing ? 'Pause the animation' : 'Play the animation');
      if (state.playing) loop();
    });
    $('view-restart').addEventListener('click', () => { buildScene(S.nodes[state.node]); });
    const main = $('view-main');
    let drag = null;
    main.addEventListener('pointerdown', ev => { if (scene && scene.kind === 'unitary') { drag = [ev.clientX, ev.clientY, state.view.az, state.view.el]; main.setPointerCapture(ev.pointerId); } });
    main.addEventListener('pointermove', ev => { if (!drag) return; state.view.az = drag[2] + (ev.clientX - drag[0]) / 120; state.view.el = Math.max(-1.3, Math.min(1.3, drag[3] + (ev.clientY - drag[1]) / 120)); if (!state.playing) draw(); });
    main.addEventListener('pointerup', () => { drag = null; });
    new ResizeObserver(() => draw()).observe($('views'));
    document.addEventListener('visibilitychange', () => { if (document.hidden) cancelAnimationFrame(frame); else if (state.playing) loop(); });
  }

  window.FieldBridgeInstrument = {state, load, start, applyEdge, undo, route, walkTo, setPreparation, byId, out, into,
                                  onChange: f => listeners.push(f), light};
  document.addEventListener('DOMContentLoaded', setup);
})();
