/* The constructor around the identity of a realization, I_real = ((Ω, Ξ); C, R, P; A). Each component is a socket
 * that shows what it holds. Opening a socket shows its content, the tools that act on it, and the parts that can be
 * attached to it: the same mechanism on another carrier or in another form (detach and attach), or a part that
 * changes the mechanism. A part is attached by a click or by dragging it onto its socket; while it is held over the
 * socket, the mechanism it would give is shown beside the formula. A part fits only its own socket. */
(function () {
  'use strict';
  const S = window.FIELDBRIDGE_SITE, I = window.FieldBridgeInstrument, M = window.FieldBridgeMechanisms;
  if (!S || !I || !M) return;
  const $ = id => document.getElementById(id);
  const esc = s => String(s).replace(/[&<>"']/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[c]));
  const strip = s => String(s == null ? '' : s).replace(/<[^>]+>/g, '');
  const SLOTS = ['Omega', 'Xi', 'C', 'R', 'P', 'A'];
  const TOOL = {Omega: 'add', P: 'protocol', A: 'plane'};
  const GREEK = {eps: 'ε', alpha: 'α', gamma: 'γ', mu: 'μ', lam: 'λ', psi: 'ψ', delta: 'δ', Delta: 'Δ', sigma: 'σ', omega: 'ω', nu: 'ν', kappa: 'κ', beta: 'β'};
  const sym = k => GREEK[k] || k;
  let open = 'Omega', extra = null, card = null, observer = null;

  const clip = (s, n = 38) => { s = strip(s).replace(/\s+/g, ' ').trim(); return s.length > n ? s.slice(0, n - 1) + '…' : s; };
  const parts = slot => {
    const id = I.state.node;
    return [...(I.out[id] || []).filter(e => e.slot === slot).map(e => [e, false]), ...(I.into[id] || []).filter(e => e.slot === slot).map(e => [e, true])];
  };
  function summary(rec, slot) {
    if (slot === 'Omega') {
      const line = document.querySelector('#slot-Omega .equation-line');
      if (line) {
        const n = document.querySelectorAll('#slot-Omega .equation-line').length, c = line.cloneNode(true);
        c.querySelectorAll('button').forEach(b => { const s = document.createElement('span'); s.className = b.classList.contains('off') ? 'cf-off' : ''; s.innerHTML = b.innerHTML; b.replaceWith(s); });
        return c.innerHTML + (n > 1 ? `<span class="cf-more">+ ${n - 1} more equation${n > 2 ? 's' : ''}</span>` : '');
      }
    }
    if (slot === 'A') {
      const p = Object.entries(I.state.params || {});
      if (p.length) return esc(p.slice(0, 3).map(([k, v]) => `${sym(k)} = ${Number(v).toPrecision(3).replace(/\.?0+$/, '')}`).join(', ') + (p.length > 3 ? ', …' : ''));
    }
    if (slot === 'P') {
      const on = document.querySelector('#slot-P [data-protocol][aria-pressed="true"]');
      if (on) return esc(clip(on.textContent));
    }
    return esc(clip(rec.slots && rec.slots[slot]));
  }

  // --------------------------------------------------------------------------------------------- the formula
  function build() {
    card = document.querySelector('#instrument .expression-card');
    if (!card || $('composer')) return;
    const el = document.createElement('div');
    el.className = 'composer'; el.id = 'composer';
    const sock = s => `<button type="button" class="socket" data-socket="${s}" style="--slot: var(--${{Omega: 'omega', Xi: 'xi', C: 'c', R: 'r', P: 'p', A: 'a'}[s]})" aria-pressed="false" aria-controls="slots">
        <span class="socket-sym">${S.slots[s].symbol}</span><span class="socket-name">${esc(S.slots[s].name)}</span>
        <span class="socket-val" id="sv-${s}"></span><span class="socket-n" id="sn-${s}"></span></button>`;
    el.innerHTML = `
      <div class="cf" role="group" aria-label="The realization as six components; open one to change it">
        <div class="cf-row"><span class="cf-lhs">I<sub>real</sub> =</span><span class="cf-p">((</span>${sock('Omega')}<span class="cf-p">,</span>${sock('Xi')}<span class="cf-p">);</span></div>
        <div class="cf-row cf-row2">${sock('C')}<span class="cf-p">,</span>${sock('R')}<span class="cf-p">,</span>${sock('P')}<span class="cf-p">;</span>${sock('A')}<span class="cf-p">)</span></div>
        <p class="cf-note">The core is the operation Ω on the carrier Ξ; C closes the equations, R is measured, P prepares and drives, A are the parameters. Open a component to see it, change it, or attach a part to it.</p>
      </div>
      <p class="drawer-head" id="drawer-head"></p>`;
    card.querySelector('.realization-head').insertAdjacentElement('afterend', el);
    // the resulting mechanism and the construction: at the top of the response beside the formula, or under the
    // formula when the page has one column
    const side = document.createElement('div');
    side.className = 'composer-side'; side.id = 'composer-side';
    side.innerHTML = `      <div class="cid" aria-live="polite">
        <div class="cid-now"><span class="card-label">Resulting mechanism</span><span id="cid-glyph"></span><b id="cid-name"></b><span id="cid-where" class="muted"></span></div>
        <div class="cid-next" id="cid-next" hidden><span class="cid-arrow">→</span><span id="cid-next-glyph"></span><span><b id="cid-next-name"></b><span id="cid-next-where" class="muted"></span></span></div>
      </div>
      <div class="cchain">
        <span class="card-label">Construction</span><ol id="cchain"></ol>
        <span class="cchain-tools">
          <button type="button" id="cc-undo" title="Undo the last change">↶ Undo</button>
          <button type="button" id="cc-save" title="Save the construction as a file">Save</button>
          <button type="button" id="cc-load" title="Replay a saved construction">Replay…</button>
          <button type="button" id="cc-examples" aria-pressed="false">Worked examples</button>
        </span>
      </div>`;
    const place = () => {
      const narrow = window.matchMedia('(max-width: 1080px)').matches, res = $('result');
      if (narrow) el.querySelector('.cf').insertAdjacentElement('afterend', side); else if (res) res.insertAdjacentElement('afterbegin', side);
    };
    place();
    window.matchMedia('(max-width: 1080px)').addEventListener('change', place);
    card.classList.add('composing');
    el.querySelectorAll('[data-socket]').forEach(b => {
      b.addEventListener('click', () => select(b.dataset.socket));
      b.addEventListener('dragover', ev => { const ok = dragged && dragged.slot === b.dataset.socket; if (ok) { ev.preventDefault(); b.classList.add('accept'); } else if (dragged) b.classList.add('refuse'); });
      b.addEventListener('dragleave', () => b.classList.remove('accept', 'refuse'));
      b.addEventListener('drop', ev => { ev.preventDefault(); b.classList.remove('accept', 'refuse'); const d = dragged; dragged = null; if (d && d.slot === b.dataset.socket) attach(d.edge, d.reverse); });
    });
    $('cc-undo').addEventListener('click', () => { window.FieldBridgeSite && window.FieldBridgeSite.leaveSequence(); I.undo(); });
    $('cc-save').addEventListener('click', () => { const b = $('ct-export'); if (b) b.click(); });
    $('cc-load').addEventListener('click', () => { const b = $('ct-import'); if (b) b.click(); });
    $('cc-examples').addEventListener('click', () => { extra = extra === 'examples' ? null : 'examples'; arrange(); });
    // the tools panel of constructor-tools.js moves under the open socket
    const ct = $('ct'); if (ct) { $('slots').insertAdjacentElement('afterend', ct); ct.classList.add('in-drawer'); }
    observer = new MutationObserver(() => requestAnimationFrame(arrange));
    observer.observe($('slots'), {childList: true});
    I.onChange((st, info) => { if (!info) requestAnimationFrame(arrange); });
    arrange();
  }
  function select(slot) {
    extra = null;
    open = open === slot && window.matchMedia('(max-width: 700px)').matches ? null : slot;
    arrange();
    if (open && window.matchMedia('(max-width: 1080px)').matches) $('slot-' + open) && $('slot-' + open).scrollIntoView({block: 'nearest', behavior: 'smooth'});
  }
  function attach(edge, reverse) {
    window.FieldBridgeSite && window.FieldBridgeSite.leaveSequence();
    M.preview(null);
    I.applyEdge(edge, reverse);
  }

  // --------------------------------------------------------------------------------------------- after each render
  function arrange() {
    const rec = S.nodes[I.state.node]; if (!rec || !card) return;
    for (const s of SLOTS) {
      const v = $('sv-' + s), n = $('sn-' + s), k = parts(s).length, b = document.querySelector(`#composer [data-socket="${s}"]`);
      if (v) v.innerHTML = summary(rec, s);
      const built = s === 'Omega' && rec.family === 'dissipative' && rec.engine && rec.engine.kind_dynamics === 'equations';
      if (n) n.textContent = k ? `${k} part${k === 1 ? '' : 's'}` : built ? 'build' : '';
      if (b) { b.setAttribute('aria-pressed', String(open === s && !extra)); b.classList.toggle('has-parts', k > 0); }
    }
    document.querySelectorAll('#slots > li.slot').forEach(li => li.classList.toggle('open', !extra && li.dataset.slot === open));
    // the parts in the open socket: the same mechanism elsewhere (detach and attach), and parts that change it
    const li = open && $('slot-' + open), edits = li && li.querySelector('.slot-edits');
    if (edits && !edits.dataset.grouped) {
      const same = [...edits.querySelectorAll('.change.same')], other = [...edits.querySelectorAll('.change:not(.same)')];
      edits.dataset.grouped = '1'; edits.innerHTML = '';
      const group = (title, list) => { if (!list.length) return; const h = document.createElement('p'); h.className = 'parts-head'; h.textContent = title; edits.append(h, ...list); };
      group(open === 'Xi' ? 'Detach the mechanism and attach it to another carrier' : 'Same mechanism, written another way', same);
      group('Attach a part that changes the mechanism', other);
    }
    document.querySelectorAll('#slots .change').forEach(wireCard);
    const head = $('drawer-head');
    if (head) {
      const k = open ? parts(open).length : 0;
      head.innerHTML = extra === 'examples' ? 'Worked examples, built in the page step by step.'
        : open ? `<b style="color: var(--${{Omega: 'omega', Xi: 'xi', C: 'c', R: 'r', P: 'p', A: 'a'}[open]})">${S.slots[open].symbol} · ${esc(S.slots[open].name)}</b>: ${k ? `${k} part${k === 1 ? '' : 's'} can be attached here. Drag one onto the socket ${S.slots[open].symbol} above, or click it.` : 'no stored part fits this socket for this realization.'}${TOOL[open] && toolUsable(rec) ? ' The tool below builds a new one.' : ''}`
        : 'Open a component of the formula to change it.';
    }
    const ct = $('ct');
    if (ct) {
      const key = extra === 'examples' ? 'examples' : open && TOOL[open];
      const usable = toolUsable(rec);
      ct.hidden = !key || (!usable && key !== 'examples');
      if (key && !ct.hidden) { const t = ct.querySelector(`[data-ct="${key}"]`); if (t && t.getAttribute('aria-selected') !== 'true') t.click(); }
    }
    $('cc-examples').setAttribute('aria-pressed', String(extra === 'examples'));
    $('cc-undo').disabled = I.state.path.length < 2;
    identity(rec);
    chain();
  }
  const toolUsable = rec => rec.family === 'dissipative' && rec.engine && rec.engine.kind_dynamics === 'equations';
  function identity(rec, next) {
    $('cid-glyph').innerHTML = M.glyph(rec.class, 'small');
    $('cid-name').textContent = rec.local ? 'not classified: built here' : (S.classes[rec.class] || rec.class);
    $('cid-where').textContent = ' · ' + strip(rec.name) + (rec.universal ? ', canonical form' : ', ' + rec.field);
  }
  function ghost(target) {
    const box = $('cid-next'); if (!box) return;
    if (!target) { box.hidden = true; return; }
    const here = S.nodes[I.state.node];
    $('cid-next-glyph').innerHTML = M.glyph(target.class, 'small');
    $('cid-next-name').textContent = target.class === here.class ? 'same mechanism' : (S.classes[target.class] || target.class);
    $('cid-next-where').textContent = ' · ' + strip(target.name) + (target.universal ? ', canonical form' : ', ' + target.field);
    box.hidden = false;
  }
  let dragged = null;
  function wireCard(b) {
    if (b.dataset.cmp) return; b.dataset.cmp = '1';
    const edge = I.byId[b.dataset.edge], reverse = b.dataset.reverse === '1', target = edge && S.nodes[reverse ? edge.from : edge.to];
    if (!edge) return;
    b.draggable = true;
    b.addEventListener('dragstart', ev => { dragged = {edge, reverse, slot: edge.slot}; ev.dataTransfer.setData('text/plain', edge.id); ev.dataTransfer.effectAllowed = 'move'; ghost(target); $('composer').classList.add('dragging'); });
    b.addEventListener('dragend', () => { dragged = null; ghost(null); $('composer').classList.remove('dragging'); document.querySelectorAll('#composer .socket').forEach(s => s.classList.remove('accept', 'refuse')); });
    b.addEventListener('mouseenter', () => ghost(target)); b.addEventListener('focus', () => ghost(target));
    b.addEventListener('mouseleave', () => { if (!dragged) ghost(null); }); b.addEventListener('blur', () => ghost(null));
  }
  function chain() {
    const list = $('cchain'); if (!list) return;
    list.innerHTML = I.state.path.map((p, i) => {
      const r = S.nodes[p.node]; if (!r) return '';
      const e = p.edge && I.byId[p.edge], slot = e ? e.slot : p.add ? 'Omega' : null;
      const step = slot ? `<span class="cc-step" style="color: var(--${{Omega: 'omega', Xi: 'xi', C: 'c', R: 'r', P: 'p', A: 'a'}[slot]})" title="${esc(e ? strip(e.change) : 'term added here, not verified')}">${S.slots[slot].symbol}${p.add ? '+' : ''}</span>` : '';
      return `<li>${step}<button type="button" class="cc-node${i === I.state.path.length - 1 ? ' current' : ''}" data-back="${i}" title="${esc(strip(r.name))}">${M.glyph(r.local ? (S.nodes[r.local.parent] || r).class : r.class, 'small')}<span>${esc(clip(r.name, 26))}</span></button></li>`;
    }).join('');
    list.querySelectorAll('[data-back]').forEach(b => b.addEventListener('click', () => {
      const k = Number(b.dataset.back); while (I.state.path.length - 1 > k) I.undo();
    }));
    list.scrollLeft = list.scrollWidth;
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', () => setTimeout(build, 0)); else setTimeout(build, 0);
})();
