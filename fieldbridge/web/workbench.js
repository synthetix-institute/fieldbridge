/* Workbench: three studies on published equations that the visitor changes and runs in the browser.
   1. The drive history of a Josephson bifurcation amplifier selects its state (Siddiqi et al., PRL 93, 207002, 2004).
   2. Two coupled condensates with the same measured imbalance and current reach different outcomes
      (two-mode equations; Raghavan et al., PRA 59, 620, 1999; Albiez et al., PRL 95, 010402, 2005).
   3. The difference appears only when tunnelling and interaction act together.
   Every curve is integrated here; the folds of study 1 and the boundary of study 2 are written in closed form. */
(function () {
  'use strict';
  const $ = (id) => document.getElementById(id);
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim() || '#888';

  function fit(canvas) {
    const box = canvas.getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
    const w = Math.max(10, box.width), h = Math.max(10, box.height);
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) { canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr); }
    const ctx = canvas.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, w, h);
    return { ctx, w, h };
  }
  function axes(g, box, xr, yr, xl, yl) {
    const { ctx } = g, [x0, y0, x1, y1] = box;
    ctx.strokeStyle = css('--axis'); ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x0, y1); ctx.lineTo(x1, y1); ctx.stroke();
    ctx.fillStyle = css('--muted'); ctx.font = '12px ' + css('--sans');
    ctx.textAlign = 'center'; ctx.fillText(xl, (x0 + x1) / 2, y1 + 30);
    ctx.save(); ctx.translate(x0 - 34, (y0 + y1) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(yl, 0, 0); ctx.restore();
    ctx.font = '11px ' + css('--sans'); ctx.fillStyle = css('--faint');
    ctx.textAlign = 'center'; ctx.fillText(fmt(xr[0]), x0, y1 + 14); ctx.fillText(fmt(xr[1]), x1, y1 + 14);
    ctx.textAlign = 'right'; ctx.fillText(fmt(yr[0]), x0 - 5, y1); ctx.fillText(fmt(yr[1]), x0 - 5, y0 + 9);
    const X = (v) => x0 + (v - xr[0]) / (xr[1] - xr[0]) * (x1 - x0), Y = (v) => y1 - (v - yr[0]) / (yr[1] - yr[0]) * (y1 - y0);
    return { X, Y };
  }
  const fmt = (v) => (Math.abs(v) >= 10 || v === 0 || Number.isInteger(v)) ? String(Math.round(v * 10) / 10) : v.toFixed(1);
  function line(ctx, pts, color, width, dash) {
    ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash || []); ctx.beginPath();
    pts.forEach(([x, y], i) => (i ? ctx.lineTo(x, y) : ctx.moveTo(x, y))); ctx.stroke(); ctx.setLineDash([]);
  }
  function rk4(f, y, t0, t1, dt, sample, every) {
    const out = [[t0, ...y]]; let t = t0, k = 0;
    while (t < t1 - 1e-12) {
      const h = Math.min(dt, t1 - t);
      const a = f(t, y), b = f(t + h / 2, y.map((v, i) => v + h / 2 * a[i])), c = f(t + h / 2, y.map((v, i) => v + h / 2 * b[i])), d = f(t + h, y.map((v, i) => v + h * c[i]));
      y = y.map((v, i) => v + h / 6 * (a[i] + 2 * b[i] + 2 * c[i] + d[i])); t += h; k += 1;
      if (sample && k % every === 0) out.push([t, ...y]);
    }
    if (!sample || k % every) out.push([t, ...y]);
    return out;
  }

  /* ---------- study 1: the drive history of the amplifier ---------- */
  const JBA = { U: 0.1, Delta: 2.0, kappa: 1.0, Ffinal: 3.0, seg: 40 };
  const quad = (a, b, c) => [(-b - Math.sqrt(b * b - 4 * a * c)) / (2 * a), (-b + Math.sqrt(b * b - 4 * a * c)) / (2 * a)];
  const NFOLD = quad(3 * JBA.U ** 2, -4 * JBA.U * JBA.Delta, JBA.Delta ** 2 + JBA.kappa ** 2 / 4);
  const F2 = (n) => n * ((JBA.U * n - JBA.Delta) ** 2 + JBA.kappa ** 2 / 4);
  const FUP = Math.sqrt(F2(NFOLD[0])), FDOWN = Math.sqrt(F2(NFOLD[1]));
  function jbaRun(levels) {
    // the drive moves linearly to each level over 8 time units and is held there for the rest of the segment
    const knots = [[0, 0]]; levels.forEach((L, i) => { knots.push([i * JBA.seg + 8, L]); knots.push([(i + 1) * JBA.seg, L]); });
    const drive = (t) => { for (let i = 1; i < knots.length; i++) if (t <= knots[i][0]) { const [ta, fa] = knots[i - 1], [tb, fb] = knots[i]; return fa + (fb - fa) * (t - ta) / Math.max(tb - ta, 1e-9); } return knots[knots.length - 1][1]; };
    const f = (t, [x, y]) => { const r = JBA.U * (x * x + y * y) - JBA.Delta; return [r * y - JBA.kappa / 2 * x, -r * x - drive(t) - JBA.kappa / 2 * y]; };
    const T = levels.length * JBA.seg, s = rk4(f, [0, 0], 0, T, 0.01, true, 10);
    return { t: s.map((q) => q[0]), F: s.map((q) => drive(q[0])), n: s.map((q) => q[1] * q[1] + q[2] * q[2]), T };
  }
  const PRESETS = { direct: [3.0, 3.0, 3.0], write: [4.5, 3.0, 3.0], erase: [4.5, 1.6, 3.0] };
  const st1 = { levels: PRESETS.write.slice(), committed: null, history: [] };
  function outcome(n) { return n > (NFOLD[0] + NFOLD[1]) / 2 ? 'high state' : 'low state'; }
  function draw1() {
    const draft = jbaRun(st1.levels), done = st1.committed && jbaRun(st1.committed);
    // drive against time
    let g = fit($('wb1-drive')), B = [44, 10, g.w - 10, g.h - 42], A = axes(g, B, [0, draft.T], [0, 6], 'time', 'drive F');
    [[FUP, 'low state ends'], [FDOWN, 'high state ends']].forEach(([F, lab]) => { line(g.ctx, [[B[0], A.Y(F)], [B[2], A.Y(F)]], css('--faint'), 1, [3, 4]); g.ctx.fillStyle = css('--faint'); g.ctx.textAlign = 'right'; g.ctx.font = '11px ' + css('--sans'); g.ctx.fillText(lab + ' (' + F.toFixed(2) + ')', B[2], A.Y(F) - 4); });
    if (done) line(g.ctx, done.t.map((t, i) => [A.X(t), A.Y(done.F[i])]), css('--faint'), 2);
    line(g.ctx, draft.t.map((t, i) => [A.X(t), A.Y(draft.F[i])]), css('--p'), 2.4, done ? [6, 4] : null);
    // photon number against time
    g = fit($('wb1-n')); B = [44, 10, g.w - 10, g.h - 42]; A = axes(g, B, [0, draft.T], [0, 30], 'time', 'photons in the resonator');
    if (done) line(g.ctx, done.t.map((t, i) => [A.X(t), A.Y(done.n[i])]), css('--faint'), 2);
    line(g.ctx, draft.t.map((t, i) => [A.X(t), A.Y(draft.n[i])]), css('--xi'), 2.4, done ? [6, 4] : null);
    // steady states against drive, with the path of the run
    g = fit($('wb1-s')); B = [44, 10, g.w - 10, g.h - 42]; A = axes(g, B, [0, 6], [0, 30], 'drive F', 'photons');
    const curve = []; for (let n = 0; n <= 30; n += 0.05) curve.push([Math.sqrt(F2(n)), n]);
    const part = (lo, hi, dash) => line(g.ctx, curve.filter(([, n]) => n >= lo && n <= hi).map(([F, n]) => [A.X(F), A.Y(n)]), css('--rule-strong'), dash ? 1.5 : 3, dash ? [4, 4] : null);
    part(0, NFOLD[0], false); part(NFOLD[0], NFOLD[1], true); part(NFOLD[1], 30, false);
    line(g.ctx, draft.F.map((F, i) => [A.X(F), A.Y(draft.n[i])]), css('--xi'), 1.6);
    const last = draft.n[draft.n.length - 1]; g.ctx.fillStyle = css('--xi'); g.ctx.beginPath(); g.ctx.arc(A.X(JBA.Ffinal), A.Y(last), 5, 0, 2 * Math.PI); g.ctx.fill();
    $('wb1-result').textContent = 'At the final drive F = 3 the resonator holds ' + last.toFixed(2) + ' photons: the ' + outcome(last) + '.' + (done ? ' The previous protocol (grey) ended in the ' + outcome(done.n[done.n.length - 1]) + '.' : '');
    ['wb1-peak', 'wb1-low'].forEach((id, i) => { $(id + '-v').textContent = st1.levels[i].toFixed(2); });
    $('wb1-apply').disabled = st1.committed && st1.committed.join() === st1.levels.join();
    $('wb1-undo').disabled = !st1.history.length;
  }
  function commit1() { if (st1.committed) st1.history.push(st1.committed.slice()); st1.committed = st1.levels.slice(); log('Amplifier: drive raised to ' + st1.levels[0].toFixed(2) + ', lowered to ' + st1.levels[1].toFixed(2) + ', then held at 3.00 — ' + outcome(jbaRun(st1.levels).n.slice(-1)[0]) + '.'); draw1(); }
  function undo1() { if (!st1.history.length) return; st1.committed = st1.history.pop(); st1.levels = st1.committed.slice(); syncSliders(); draw1(); }
  function syncSliders() { $('wb1-peak').value = st1.levels[0]; $('wb1-low').value = st1.levels[1]; }

  /* ---------- studies 2 and 3: two coupled condensates ---------- */
  const st2 = { L: 15, z: 0.4, v: -0.55, J: true, I: true };
  function becRun(z0, v0, cosSign, L, J, I) {
    const r0 = Math.sqrt(1 - z0 * z0), s = Math.max(-1, Math.min(1, -v0 / Math.max(r0, 1e-9))), phi0 = cosSign > 0 ? Math.asin(s) : Math.PI - Math.asin(s);
    const j = J ? 1 : 0, u = I ? L : 0;
    // tunnelling off: the measured current is set by hand and the atoms do not move; the comparison is with the same z and current zero
    const f = (t, [z, p]) => { const r = Math.sqrt(Math.max(1 - z * z, 1e-12)); return [-j * r * Math.sin(p), u * z + j * z * Math.cos(p) / r]; };
    const s2 = rk4(f, [z0, phi0], 0, 30, 0.002, true, 25);
    return { t: s2.map((q) => q[0]), z: s2.map((q) => q[1]) };
  }
  const crosses = (z) => z.some((v) => v * z[0] < 0);
  function openShare(L) { let n = 0, k = 0; for (let i = 0; i < 140; i++) for (let j = 0; j < 140; j++) { const z = -1 + (i + 0.5) / 70, v = -1 + (j + 0.5) / 70; if (z * z + v * v >= 1) continue; k++; if (Math.abs(L * z * z / 2 - 1) < Math.sqrt(1 - z * z - v * v)) n++; } return n / k; }
  function draw2() {
    const L = st2.L, g = fit($('wb2-disc')), R = Math.min(g.w, g.h) / 2 - 22, cx = g.w / 2, cy = g.h / 2 - 6, ctx = g.ctx;
    const img = ctx.createImageData(Math.ceil(2 * R), Math.ceil(2 * R)), cO = hex(css('--r')), cD = hex(css('--xi'));
    for (let i = 0; i < img.width; i++) for (let j = 0; j < img.height; j++) {
      const z = (i - R) / R, v = (R - j) / R; if (z * z + v * v >= 1) continue;
      const open = Math.abs(L * z * z / 2 - 1) < Math.sqrt(1 - z * z - v * v), c = open ? cO : cD, k = 4 * (j * img.width + i);
      img.data[k] = c[0]; img.data[k + 1] = c[1]; img.data[k + 2] = c[2]; img.data[k + 3] = open ? 120 : 40;
    }
    const dpr = window.devicePixelRatio || 1, tmp = document.createElement('canvas'); tmp.width = img.width; tmp.height = img.height; tmp.getContext('2d').putImageData(img, 0, 0);
    ctx.drawImage(tmp, cx - R, cy - R);
    ctx.strokeStyle = css('--rule-strong'); ctx.beginPath(); ctx.arc(cx, cy, R, 0, 2 * Math.PI); ctx.stroke();
    ctx.fillStyle = css('--muted'); ctx.font = '12px ' + css('--sans'); ctx.textAlign = 'center';
    ctx.fillText('imbalance z = (left − right) / total', cx, cy + R + 17);
    ctx.save(); ctx.translate(cx - R - 9, cy); ctx.rotate(-Math.PI / 2); ctx.fillText('current dz/dt', 0, 0); ctx.restore();
    ctx.fillStyle = css('--ink'); ctx.beginPath(); ctx.arc(cx + st2.z * R, cy - st2.v * R, 6, 0, 2 * Math.PI); ctx.fill();
    $('wb2-share').textContent = Math.round(100 * openShare(L)) + '%';
    $('wb2-L-v').textContent = L.toFixed(1);
    $('wb2-state').textContent = 'z = ' + st2.z.toFixed(2) + ', dz/dt = ' + st2.v.toFixed(2);
    drawRuns('wb2-runs', st2.J && true, st2.I && true, 'wb2-result', true);
    drawRuns('wb3-runs', st2.J, st2.I, 'wb3-result', false);
  }
  function drawRuns(id, J, I, resId, full) {
    const a = becRun(st2.z, st2.v, 1, st2.L, full || J, full || I), b = becRun(st2.z, st2.v, -1, st2.L, full || J, full || I);
    const g = fit($(id)), B = [44, 10, g.w - 10, g.h - 42], A = axes(g, B, [0, 30], [-1, 1], 'time (1 / tunnelling rate)', 'imbalance z');
    line(g.ctx, [[B[0], A.Y(0)], [B[2], A.Y(0)]], css('--grid'), 1);
    line(g.ctx, a.t.map((t, i) => [A.X(t), A.Y(a.z[i])]), css('--xi'), 2.2);
    line(g.ctx, b.t.map((t, i) => [A.X(t), A.Y(b.z[i])]), css('--r'), 2.2, [7, 4]);
    const same = a.z.every((v, i) => Math.abs(v - b.z[i]) < 1e-6);
    const word = (z) => (crosses(z) ? 'atoms flow back and forth' : 'atoms stay on one side');
    $(resId).textContent = same ? 'The two preparations give the same curve: ' + word(a.z) + '.' : 'Phase in step (solid): ' + word(a.z) + '. Phase out of step (dashed): ' + word(b.z) + '.';
  }
  function hex(c) { const m = c.match(/^#?([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/i); return m ? [1, 2, 3].map((i) => parseInt(m[i], 16)) : [128, 128, 128]; }

  /* ---------- log of changes, shared ---------- */
  const changes = [];
  function log(text) { changes.push({ time: new Date().toISOString(), change: text }); const li = document.createElement('li'); li.textContent = text; $('wb-log').prepend(li); $('wb-export').disabled = false; }
  function exportLog() {
    const blob = new Blob([JSON.stringify({ schema: 'fieldbridge-workbench-log/1', changes, amplifier: { specification: 'examples/quantum/semiclassical/jba_resonator.json', drive_levels: st1.committed, segment: JBA.seg }, condensates: { Lambda: st2.L, z: st2.z, current: st2.v } }, null, 1)], { type: 'application/json' });
    const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = 'fieldbridge-workbench.json'; a.click(); URL.revokeObjectURL(a.href);
  }

  function init() {
    if (!$('workbench')) return;
    syncSliders(); commit1();
    $('wb1-peak').addEventListener('input', (e) => { st1.levels[0] = +e.target.value; draw1(); });
    $('wb1-low').addEventListener('input', (e) => { st1.levels[1] = +e.target.value; draw1(); });
    $('wb1-apply').addEventListener('click', commit1);
    $('wb1-undo').addEventListener('click', undo1);
    document.querySelectorAll('[data-wb1-preset]').forEach((b) => b.addEventListener('click', () => { st1.levels = PRESETS[b.dataset.wb1Preset].slice(); syncSliders(); draw1(); }));
    $('wb2-L').addEventListener('input', (e) => { st2.L = +e.target.value; draw2(); });
    $('wb2-L').addEventListener('change', () => log('Condensates: interaction set to Λ = ' + st2.L.toFixed(1) + ' (share with an open outcome ' + $('wb2-share').textContent + ').'));
    $('wb2-disc').addEventListener('click', (e) => {
      const box = e.target.getBoundingClientRect(), R = Math.min(box.width, box.height) / 2 - 22, z = (e.clientX - box.left - box.width / 2) / R, v = -(e.clientY - box.top - (box.height / 2 - 6)) / R;
      if (z * z + v * v >= 0.995) return; st2.z = z; st2.v = v; draw2(); log('Condensates: measured state z = ' + z.toFixed(2) + ', dz/dt = ' + v.toFixed(2) + '.');
    });
    $('wb3-J').addEventListener('change', (e) => { st2.J = e.target.checked; draw2(); log('Condensates: tunnelling ' + (st2.J ? 'on' : 'off') + ', interaction ' + (st2.I ? 'on' : 'off') + '.'); });
    $('wb3-I').addEventListener('change', (e) => { st2.I = e.target.checked; draw2(); log('Condensates: tunnelling ' + (st2.J ? 'on' : 'off') + ', interaction ' + (st2.I ? 'on' : 'off') + '.'); });
    $('wb-export').addEventListener('click', exportLog);
    draw2();
    let pending = 0; window.addEventListener('resize', () => { cancelAnimationFrame(pending); pending = requestAnimationFrame(() => { draw1(); draw2(); }); });
    if (window.matchMedia) window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { draw1(); draw2(); });
  }
  window.FieldBridgeWorkbench = { jbaRun, becRun, openShare, FUP, FDOWN };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
