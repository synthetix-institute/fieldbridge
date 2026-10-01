/* Drawings of the page: canvases for the dynamics and the diagrams, and the carrier diagrams of the expression. */
(function (root) {
  'use strict';
  const css = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  function palette() {
    return {ink: css('--ink'), muted: css('--muted'), faint: css('--faint'), rule: css('--rule'), grid: css('--grid'),
            axis: css('--axis'), panel: css('--panel'), accent: css('--accent'), omega: css('--omega'), xi: css('--xi'),
            c: css('--c'), r: css('--r'), p: css('--p'), a: css('--a'), good: css('--good'), warn: css('--warn')};
  }
  const SERIES = ['--xi', '--p', '--c', '--r', '--omega', '--a'];
  function fit(canvas) {
    const box = canvas.getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
    const w = Math.max(1, Math.round(box.width)), h = Math.max(1, Math.round(box.height));
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) { canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr); }
    const ctx = canvas.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, w, h);
    return {ctx, w, h, c: palette()};
  }
  const fmt = (v, n = 3) => { if (!Number.isFinite(v)) return '—'; const a = Math.abs(v); const s = a !== 0 && (a >= 1e4 || a < 1e-3) ? v.toExponential(1) : String(+v.toPrecision(n)); return s.replace('-', '−'); };
  // a value with its uncertainty: the uncertainty to two significant digits, the value to the same decimal place
  const pm = (v, err) => {
    if (!Number.isFinite(v) || !(err > 0) || !Number.isFinite(err)) return fmt(v, 5);
    const d = Math.min(12, Math.max(0, 1 - Math.floor(Math.log10(err))));
    return `${v.toFixed(d)} ± ${err.toFixed(d)}`.replace(/-/g, '−');
  };
  function ticks(lo, hi, n = 4) {
    const span = hi - lo || 1, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0))), f = step0 / mag;
    const step = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * mag, out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9 * span; v += step) out.push(Math.abs(v) < 1e-12 * span ? 0 : v);
    return out;
  }
  /** Tick labels with as many decimals as the spacing of the ticks needs. */
  function stepFormat(t) {
    if (t.length < 2) return v => fmt(v);
    const step = Math.min(...t.slice(1).map((v, i) => Math.abs(v - t[i])));
    const a = Math.max(...t.map(Math.abs));
    if (a >= 1e4 || (a > 0 && a < 1e-3)) return v => fmt(v);
    const d = Math.max(0, Math.min(6, Math.ceil(-Math.log10(step) + 1e-9)));
    return v => v.toFixed(d).replace('-', '−');
  }
  /** Axes in a box; returns the maps from data to pixels. */
  function frame(g, box, xr, yr, o = {}) {
    const {ctx, c} = g, l = box.x + (o.left ?? 40), t = box.y + (o.top ?? 14), r = box.x + box.w - (o.right ?? 10), b = box.y + box.h - (o.bottom ?? 30);
    const X = x => l + (x - xr[0]) / (xr[1] - xr[0]) * (r - l), Y = y => b - (y - yr[0]) / (yr[1] - yr[0]) * (b - t);
    ctx.font = '11px ' + css('--sans'); ctx.lineWidth = 1;
    if (!o.noGrid) {
      ctx.strokeStyle = c.grid;
      for (const v of o.yTicks || ticks(yr[0], yr[1])) { ctx.beginPath(); ctx.moveTo(l, Y(v)); ctx.lineTo(r, Y(v)); ctx.stroke(); }
    }
    ctx.strokeStyle = c.axis; ctx.beginPath(); ctx.moveTo(l, t); ctx.lineTo(l, b); ctx.lineTo(r, b); ctx.stroke();
    ctx.fillStyle = c.muted; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    const yt = o.yTicks || ticks(yr[0], yr[1]), xt = o.xTicks || ticks(xr[0], xr[1]);
    const yf = o.yFormat || stepFormat(yt), xf = o.xFormat || stepFormat(xt);
    for (const v of yt) ctx.fillText(yf(v), l - 5, Y(v));
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (const v of xt) ctx.fillText(xf(v), X(v), b + 5);
    if (o.xLabel) { ctx.textAlign = 'right'; ctx.fillText(o.xLabel, r, b + 17); }
    if (o.yLabel) { ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText(o.yLabel, l + 4, t); }
    return {X, Y, l, r, t, b};
  }
  function line(ctx, pts, X, Y, color, width = 1.8, dash = []) {
    ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash); ctx.beginPath();
    let pen = false;
    for (const [x, y] of pts) { if (!Number.isFinite(y)) { pen = false; continue; } if (pen) ctx.lineTo(X(x), Y(y)); else ctx.moveTo(X(x), Y(y)); pen = true; }
    ctx.stroke(); ctx.setLineDash([]);
  }
  const REF = [3, 3];  // the dash of a reference: the realization before a change
  function ring(ctx, x, y, r, color) { ctx.strokeStyle = color; ctx.lineWidth = 1.3; ctx.setLineDash(REF); ctx.beginPath(); ctx.arc(x, y, r, 0, 2 * Math.PI); ctx.stroke(); ctx.setLineDash([]); }
  function dot(ctx, x, y, r, fill, stroke) { ctx.beginPath(); ctx.arc(x, y, r, 0, 2 * Math.PI); if (fill) { ctx.fillStyle = fill; ctx.fill(); } if (stroke) { ctx.strokeStyle = stroke; ctx.lineWidth = 1.5; ctx.stroke(); } }
  function label(ctx, text, x, y, color, align = 'left', font = '12px') { ctx.font = font + ' ' + css('--sans'); ctx.fillStyle = color; ctx.textAlign = align; ctx.textBaseline = 'middle'; ctx.fillText(text, x, y); }
  function arrow(ctx, x0, y0, x1, y1, color, width = 2) {
    ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = width;
    ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
    const a = Math.atan2(y1 - y0, x1 - x0), s = 5 + width * 1.5;
    if (Math.hypot(x1 - x0, y1 - y0) < s) return;
    ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x1 - s * Math.cos(a - 0.4), y1 - s * Math.sin(a - 0.4)); ctx.lineTo(x1 - s * Math.cos(a + 0.4), y1 - s * Math.sin(a + 0.4)); ctx.closePath(); ctx.fill();
  }

  // ---------------------------------------------------------------------------------------- the Bloch sphere
  /** m: the expectation of the frame operators, in units of the spin; the observable along z, Omega in the x-z plane. */
  function sphere(canvas, s) {
    const g = fit(canvas), {ctx, w, h, c} = g, R = Math.min(w, h) * 0.36, cx = w / 2, cy = h / 2 + 6;
    const az = s.view?.az ?? -0.6, el = s.view?.el ?? 0.32;
    const P = ([x, y, z]) => {
      const x1 = x * Math.cos(az) - y * Math.sin(az), y1 = x * Math.sin(az) + y * Math.cos(az);
      return [cx + R * x1, cy - R * (z * Math.cos(el) - y1 * Math.sin(el)), y1 * Math.cos(el) + z * Math.sin(el)];
    };
    ctx.strokeStyle = c.rule; ctx.lineWidth = 1.2; ctx.beginPath(); ctx.arc(cx, cy, R, 0, 2 * Math.PI); ctx.stroke();
    const circle = f => { const pts = []; for (let k = 0; k <= 96; k++) pts.push(P(f(2 * Math.PI * k / 96))); return pts; };
    for (const pts of [circle(a => [Math.cos(a), Math.sin(a), 0]), circle(a => [Math.cos(a), 0, Math.sin(a)])]) {
      for (let k = 1; k < pts.length; k++) {
        ctx.strokeStyle = c.rule; ctx.setLineDash(pts[k][2] > 0 ? [3, 4] : []); ctx.beginPath();
        ctx.moveTo(pts[k - 1][0], pts[k - 1][1]); ctx.lineTo(pts[k][0], pts[k][1]); ctx.stroke();
      }
    }
    ctx.setLineDash([]);
    const O = P([0, 0, 0]), top = P([0, 0, 1.18]), bottom = P([0, 0, -1.08]);
    ctx.strokeStyle = c.r; ctx.lineWidth = 1.4; ctx.setLineDash([5, 4]); ctx.beginPath(); ctx.moveTo(bottom[0], bottom[1]); ctx.lineTo(top[0], top[1]); ctx.stroke(); ctx.setLineDash([]);
    label(ctx, 'R', top[0] + 6, top[1] - 4, c.r, 'left', '600 13px');
    if (s.omega) {
      const n = Math.hypot(...s.omega) || 1, u = s.omega.map(v => v / n), tip = P(u.map(v => 1.12 * v));
      arrow(ctx, O[0], O[1], tip[0], tip[1], c.omega, 2.2); label(ctx, 'Ω', tip[0] + 7, tip[1], c.omega, 'left', '600 14px');
    }
    if (s.trail && s.trail.length > 1) {
      for (let k = 1; k < s.trail.length; k++) {
        const a = P(s.trail[k - 1]), b = P(s.trail[k]);
        ctx.strokeStyle = c.xi; ctx.globalAlpha = 0.15 + 0.85 * k / s.trail.length; ctx.lineWidth = 2;
        ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
      }
      ctx.globalAlpha = 1;
    }
    if (s.m) {
      const tip = P(s.m); arrow(ctx, O[0], O[1], tip[0], tip[1], c.xi, 2.6); dot(ctx, tip[0], tip[1], 4.2, c.xi);
      const len = Math.hypot(...s.m);
      label(ctx, `|⟨J⟩|/j = ${fmt(len, 3)}`, 10, h - 12, len < 0.995 ? c.warn : c.muted, 'left', '12px');
    }
    if (s.caption) label(ctx, s.caption, 10, 14, c.muted, 'left', '12px');
  }

  // ---------------------------------------------------------------------------------------- plots
  function signal(canvas, s) {
    const g = fit(canvas), {ctx, w, h, c} = g;
    const f = frame(g, {x: 0, y: 0, w, h}, [0, s.tMax], [-1.1, 1.1], {xLabel: 't', yLabel: s.yLabel || 'f(t)', yTicks: [-1, 0, 1]});
    if (s.reference) line(ctx, s.reference, f.X, f.Y, c.faint, 1.6, REF);
    if (s.law) { const pts = []; for (let k = 0; k <= 240; k++) { const t = s.tMax * k / 240; pts.push([t, s.law(t)]); } line(ctx, pts, f.X, f.Y, c.omega, 1.5, [6, 4]); }
    (s.curves || []).forEach((cur, i) => line(ctx, cur.pts.filter(p => p[0] <= (s.cursor ?? Infinity)), f.X, f.Y, cur.color || css(SERIES[i % SERIES.length]), 2));
    if (s.cursor != null) { ctx.strokeStyle = c.faint; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(f.X(s.cursor), f.t); ctx.lineTo(f.X(s.cursor), f.b); ctx.stroke(); }
  }
  function autoRange(values, pad = 0.08, floor = 1e-6) {
    const v = values.filter(Number.isFinite); if (!v.length) return [-1, 1];
    let lo = Math.min(...v), hi = Math.max(...v); if (hi - lo < floor) { lo -= 0.5; hi += 0.5; }
    const d = (hi - lo) * pad; return [lo - d, hi + d];
  }
  function landscape(canvas, s) {
    const g = fit(canvas), {ctx, w, h, c} = g;
    // the range follows the current landscape; a reference that leaves it by more than its height is cut, not shown small
    const lo = Math.min(...s.V), hi = Math.max(...s.V), span = hi - lo || 1;
    const inside = s.reference ? s.reference.filter(v => v >= lo - span && v <= hi + span) : [];
    const yr = autoRange(s.V.concat(inside), 0.12);
    const f = frame(g, {x: 0, y: 0, w, h}, [s.xs[0], s.xs[s.xs.length - 1]], yr, {xLabel: s.xLabel || 'x', yLabel: s.yLabel || 'energy', xFormat: s.xFormat});
    if (s.reference) {
      ctx.save(); ctx.beginPath(); ctx.rect(f.l, f.t, f.r - f.l, f.b - f.t); ctx.clip();
      line(ctx, s.xs.map((x, i) => [x, s.reference[i]]), f.X, f.Y, c.faint, 1.6, REF);
      ctx.restore();
    }
    line(ctx, s.xs.map((x, i) => [x, s.V[i]]), f.X, f.Y, c.ink, 2);
    const Vat = x => { const i = Math.max(0, Math.min(s.xs.length - 1, Math.round((x - s.xs[0]) / (s.xs[1] - s.xs[0])))); return s.V[i]; };
    (s.states || []).forEach(q => dot(ctx, f.X(q), f.Y(Vat(q)), 5, null, c.good));
    (s.balls || []).forEach(x => dot(ctx, f.X(x), f.Y(Vat(x)) - 5, 3.2, c.xi));
  }
  function bifurcation(canvas, s) {
    const g = fit(canvas), {ctx, w, h, c} = g;
    const ys = s.scan.concat(s.reference || []).flatMap(r => r.points.map(p => p.s));
    const f = frame(g, {x: 0, y: 0, w, h}, s.range, autoRange(ys, 0.1), {xLabel: s.param, yLabel: s.yLabel || 'state'});
    for (const row of s.reference || []) for (const p of row.points) if (p.stable) ring(ctx, f.X(row.v), f.Y(p.s), 3.4, c.faint);
    for (const row of s.scan) for (const p of row.points) dot(ctx, f.X(row.v), f.Y(p.s), p.stable ? 2.8 : 2.2, p.stable ? c.ink : null, p.stable ? null : c.faint);
    (s.events || []).forEach(e => { ctx.strokeStyle = c.a; ctx.setLineDash([2, 3]); ctx.beginPath(); ctx.moveTo(f.X(e.value), f.t); ctx.lineTo(f.X(e.value), f.b); ctx.stroke(); ctx.setLineDash([]); });
    if (s.current != null) { ctx.strokeStyle = c.accent; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(f.X(s.current), f.t); ctx.lineTo(f.X(s.current), f.b); ctx.stroke(); }
    label(ctx, 'filled: stable · open: unstable', f.r, f.t + 4, c.muted, 'right', '11px');
  }
  function phasePlane(canvas, s) {
    const g = fit(canvas), {ctx, w, h, c} = g;
    const f = frame(g, {x: 0, y: 0, w, h}, s.xr, s.yr, {xLabel: s.xLabel, yLabel: s.yLabel});
    if (s.arrows) {
      ctx.strokeStyle = c.rule; ctx.lineWidth = 1;
      for (const a of s.arrows) { const x0 = f.X(a[0]), y0 = f.Y(a[1]), n = Math.hypot(a[2], a[3]) || 1, L = 9; ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x0 + L * a[2] / n, y0 - L * a[3] / n); ctx.stroke(); }
    }
    (s.trails || []).forEach((tr, i) => {
      const col = tr.color || c.xi; ctx.strokeStyle = col; ctx.lineWidth = 1.2;
      for (let k = 1; k < tr.pts.length; k++) { ctx.globalAlpha = 0.03 + 0.42 * k / tr.pts.length; ctx.beginPath(); ctx.moveTo(f.X(tr.pts[k - 1][0]), f.Y(tr.pts[k - 1][1])); ctx.lineTo(f.X(tr.pts[k][0]), f.Y(tr.pts[k][1])); ctx.stroke(); }
      ctx.globalAlpha = 1;
      const last = tr.pts[tr.pts.length - 1]; if (last) dot(ctx, f.X(last[0]), f.Y(last[1]), 2.2, col);
    });
    (s.reference || []).forEach(e => ring(ctx, f.X(e[0]), f.Y(e[1]), 8, c.faint));
    (s.equilibria || []).forEach(e => dot(ctx, f.X(e[0]), f.Y(e[1]), 5, e.stable ? c.good : c.panel, e.stable ? null : c.warn));
  }
  function series(canvas, s) {
    const g = fit(canvas), {ctx, w, h, c} = g;
    const vals = s.curves.flatMap(cu => cu.pts.map(p => p[1]));
    const f = frame(g, {x: 0, y: 0, w, h}, [s.t0 ?? 0, s.tMax], s.yr || autoRange(vals), {xLabel: s.xLabel ?? 't', yLabel: s.yLabel || '', xTicks: s.xTicks, xFormat: s.xFormat});
    if (s.window) { ctx.fillStyle = c.p; ctx.globalAlpha = 0.1; ctx.fillRect(f.X(s.window[0]), f.t, f.X(s.window[1]) - f.X(s.window[0]), f.b - f.t); ctx.globalAlpha = 1; }
    s.curves.forEach((cu, i) => line(ctx, cu.pts.filter(p => p[0] <= (s.cursor ?? Infinity)), f.X, f.Y, cu.color || css(SERIES[i % SERIES.length]), 1.8, cu.dash || []));
    let x = f.l + 6; s.curves.forEach((cu, i) => { if (!cu.label) return; label(ctx, cu.label, x, f.t + 6, cu.color || css(SERIES[i % SERIES.length]), 'left', '11px'); x += ctx.measureText(cu.label).width + 14; });
  }
  /** A hysteresis loop of response against drive: the major loop, the loop before a change (dashed), and an excursion
   *  drawn up to the cursor, from its turning point (ringed). Between events the response is constant: steps. */
  function hysteresis(canvas, s) {
    const g = fit(canvas), {ctx, w, h, c} = g;
    const f = frame(g, {x: 0, y: 0, w, h}, s.xr, [-1.08, 1.08], {xLabel: s.xLabel || 'drive H', yLabel: s.yLabel || 'response R'});
    const steps = pts => pts.flatMap((p, i) => (i ? [[p[0], pts[i - 1][1]], p] : [p]))
      .map(([x, y]) => [Math.min(Math.max(x, s.xr[0]), s.xr[1]), y]);
    if (s.reference) for (const b of ['up', 'down']) line(ctx, steps(s.reference[b]), f.X, f.Y, c.faint, 1.4, REF);
    for (const b of ['up', 'down']) line(ctx, steps(s.major[b]), f.X, f.Y, c.rule, 1.6);
    if (s.path && s.path.length) {
      line(ctx, steps(s.path), f.X, f.Y, c.xi, 2.2);
      const last = s.path[s.path.length - 1];
      dot(ctx, f.X(last[0]), f.Y(last[1]), 3.4, c.xi);
    }
    if (s.turn) { ring(ctx, f.X(s.turn[0]), f.Y(s.turn[1]), 6.5, c.ink); label(ctx, 'H₁', f.X(s.turn[0]) + 9, f.Y(s.turn[1]) - 11, c.muted, 'left', '11px'); }
    if (s.caption) label(ctx, s.caption, f.l + 6, f.t + 8, c.muted, 'left', '11px');
  }
  /** The magnetization in the plane of the easy axis (vertical) and the field; each moment an arrow on the circle. */
  function circle(canvas, s) {
    const g = fit(canvas), {ctx, w, h, c} = g, R = Math.min(w, h) * 0.36, cx = w / 2, cy = h / 2 + 4;
    const at = (a, r = 1) => [cx + R * r * Math.cos(a), cy - R * r * Math.sin(a)];
    ctx.strokeStyle = c.rule; ctx.lineWidth = 1.2; ctx.beginPath(); ctx.arc(cx, cy, R, 0, 2 * Math.PI); ctx.stroke();
    if (s.easyAxis !== false) {
      const [x0, y0] = at(-Math.PI / 2, 1.2), [x1, y1] = at(Math.PI / 2, 1.2);
      ctx.strokeStyle = c.faint; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke(); ctx.setLineDash([]);
      label(ctx, 'easy axis', x1 + 6, y1 + 2, c.muted, 'left', '11px');
    }
    if (s.field != null && s.h) { const [x0, y0] = at(s.field + Math.PI, 1.3), [x1, y1] = at(s.field, 1.3), [lx, ly] = at(s.field + 0.14, 1.12); arrow(ctx, x0, y0, x1, y1, c.a, 1.6); label(ctx, s.fieldLabel || 'field', lx, ly - 8, c.a, 'center', '11px'); }
    (s.reference || []).forEach(a => { const [x, y] = at(a); ring(ctx, x, y, 9.5, c.faint); });
    (s.states || []).forEach(a => { const [x, y] = at(a); dot(ctx, x, y, 6, null, c.good); });
    (s.angles || []).forEach(a => { const [x, y] = at(a, 0.92); ctx.globalAlpha = 0.55; arrow(ctx, cx, cy, x, y, c.xi, 1.4); ctx.globalAlpha = 1; });
    if (s.caption) label(ctx, s.caption, 10, 14, c.muted, 'left', '12px');
  }
  function rotors(canvas, s) {
    const g = fit(canvas), {ctx, w, h, c} = g;
    const xs = s.positions.map(p => p[0]), ys = s.positions.map(p => p[1]);
    const [x0, x1] = [Math.min(...xs), Math.max(...xs)], [y0, y1] = [Math.min(...ys), Math.max(...ys)];
    const sc = Math.min((w - 60) / Math.max(x1 - x0, 1e-9), (h - 60) / Math.max(y1 - y0, 1e-9)), L = Math.min(22, sc * 0.35);
    const P = ([x, y]) => [w / 2 + (x - (x0 + x1) / 2) * sc, h / 2 - (y - (y0 + y1) / 2) * sc];
    ctx.strokeStyle = c.rule; ctx.lineWidth = 1;
    (s.bonds || []).forEach(([i, j]) => { const a = P(s.positions[i]), b = P(s.positions[j]); ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke(); });
    s.positions.forEach((p, i) => {
      const [x, y] = P(p), a = s.angles[i], dx = L * Math.cos(a), dy = L * Math.sin(a);
      ctx.strokeStyle = c.xi; ctx.lineWidth = 5; ctx.lineCap = 'round';
      ctx.beginPath(); ctx.moveTo(x - (s.m === 2 ? dx : 0), y + (s.m === 2 ? dy : 0)); ctx.lineTo(x + dx, y - dy); ctx.stroke();
      if (s.m !== 2) dot(ctx, x + dx, y - dy, 3.5, c.p);
      ctx.lineCap = 'butt';
    });
    if (s.caption) label(ctx, s.caption, 10, 14, c.muted, 'left', '12px');
  }
  function loglog(canvas, s) {
    const g = fit(canvas), {ctx, w, h, c} = g;
    const lt = s.t.map(Math.log10), lv = s.v.map(v => Math.log10(Math.max(v, 1e-300)));
    const yr = autoRange(lv.filter(v => v > -12), 0.06);
    const f = frame(g, {x: 0, y: 0, w, h}, [lt[0], lt[lt.length - 1]], yr, {xLabel: 't', yLabel: 'SNR of the write', xFormat: v => '10' + sup(v), yFormat: v => '10' + sup(v), xTicks: intTicks(lt), yTicks: intTicks(yr)});
    if (s.reference) {
      ctx.save(); ctx.beginPath(); ctx.rect(f.l, f.t, f.r - f.l, f.b - f.t); ctx.clip();
      line(ctx, lt.map((x, i) => [x, Math.log10(Math.max(s.reference[i], 1e-300))]), f.X, f.Y, c.faint, 1.6, REF);
      ctx.restore();
    }
    if (s.guide) line(ctx, s.guide.map(([t, v]) => [Math.log10(t), Math.log10(v)]), f.X, f.Y, c.omega, 1.4, [6, 4]);
    line(ctx, lt.map((x, i) => [x, lv[i]]), f.X, f.Y, c.xi, 2.2);
    if (s.cursor != null) { ctx.strokeStyle = c.faint; ctx.beginPath(); ctx.moveTo(f.X(Math.log10(s.cursor)), f.t); ctx.lineTo(f.X(Math.log10(s.cursor)), f.b); ctx.stroke(); }
    if (s.note) label(ctx, s.note, f.r, f.t + 6, c.omega, 'right', '12px');
  }
  const SUP = {'-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹'};
  const sup = v => String(Math.round(v)).split('').map(ch => SUP[ch] || ch).join('');
  function intTicks([lo, hi]) { const out = []; const step = Math.max(1, Math.ceil((hi - lo) / 5)); for (let v = Math.ceil(lo); v <= hi; v += step) out.push(v); return out; }
  function profile(canvas, s) {
    const g = fit(canvas), {ctx, w, h, c} = g, L = s.phi.length, half = Math.floor(L / 2);
    const xs = Array.from({length: L}, (_, i) => i - half), idx = i => ((i % L) + L) % L;
    const vals = [...s.phi, ...(s.initial || [])];
    const f = frame(g, {x: 0, y: 0, w, h}, [-half, L - half - 1], s.yr || autoRange(vals, 0.1), {xLabel: 'site', yLabel: 'trace of the write'});
    if (s.initial) line(ctx, xs.map(x => [x, s.initial[idx(x)]]), f.X, f.Y, c.faint, 1.2, [3, 3]);
    line(ctx, xs.map(x => [x, s.phi[idx(x)]]), f.X, f.Y, c.xi, 2);
    if (s.caption) label(ctx, s.caption, f.r, f.t + 6, c.muted, 'right', '12px');
  }
  function phaseDrift(canvas, s) {
    const g = fit(canvas), {ctx, w, h, c} = g;
    const f = frame(g, {x: 0, y: 0, w, h}, [0, s.tMax], [0, 1], {xLabel: 't', yLabel: 'phase against the drive (cycles)', yTicks: [0, 0.5, 1]});
    for (const [t, p] of s.phases) if (t <= (s.cursor ?? Infinity)) dot(ctx, f.X(t), f.Y(p), 2.3, s.locked ? c.good : c.warn);
    if (s.caption) label(ctx, s.caption, f.r, f.t + 6, s.locked ? c.good : c.warn, 'right', '12px');
  }

  // ---------------------------------------------------------------------------------------- carriers
  const esc = t => String(t).replace(/[&<>"]/g, ch => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[ch]));
  /** A small diagram of the carrier of a realization, for the row of Xi. */
  function carrier(rec) {
    const e = rec.engine, W = 320, H = 70;
    const svg = inner => `<svg class="carrier-diagram" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(rec.facts.carrier || rec.slots.Xi)}">${inner}</svg>`;
    const spin = (x, y, lbl) => `<circle cx="${x}" cy="${y}" r="13" fill="none" stroke="var(--xi)" stroke-width="1.6"/><text x="${x}" y="${y + 4}" text-anchor="middle">${lbl}</text>`;
    if (rec.family === 'unitary') {
      const k = e.carrier.kind;
      if (k === 'qubits') {
        const n = e.carrier.n, xs = Array.from({length: n}, (_, i) => n === 1 ? W / 2 : 40 + (W - 80) * i / (n - 1));
        return svg(xs.slice(1).map((x, i) => `<line x1="${xs[i] + 13}" x2="${x - 13}" y1="30" y2="30" stroke="var(--rule-strong)" stroke-width="2"/>`).join('')
          + xs.map((x, i) => spin(x, 30, i)).join('') + `<text x="${W / 2}" y="64" text-anchor="middle">${n} spin${n > 1 ? 's' : ''}-½${e.carrier.sector ? ' · ' + e.carrier.sector.dimension + ' of ' + e.carrier.sector.full_dimension + ' states' : ''}</text>`);
      }
      if (k === 'spin') {
        const d = e.dim, ys = Array.from({length: d}, (_, i) => 10 + 40 * i / Math.max(d - 1, 1));
        return svg(ys.map((y, i) => `<line x1="${W / 2 - 50}" x2="${W / 2 + 50}" y1="${y}" y2="${y}" stroke="var(--xi)" stroke-width="2"/><text x="${W / 2 + 58}" y="${y + 4}">m = ${fmt((d - 1) / 2 - i)}</text>`).join('')
          + `<text x="${W / 2}" y="66" text-anchor="middle">one spin j = ${(d - 1) % 2 ? (d - 1) + '/2' : (d - 1) / 2}</text>`);
      }
      if (k === 'bosons') {
        const N = e.carrier.sector ? Math.round(e.carrier.sector.value) : e.carrier.max_quanta;
        const well = x => `<path d="M${x - 40},12 Q${x},62 ${x + 40},12" fill="none" stroke="var(--xi)" stroke-width="1.8"/>`;
        return svg(well(W / 2 - 55) + well(W / 2 + 55) + Array.from({length: N}, (_, i) => `<circle cx="${W / 2 - 70 + 10 * i}" cy="36" r="3.5" fill="var(--p)"/>`).join('')
          + `<text x="${W / 2}" y="66" text-anchor="middle">${N} bosons in two wells</text>`);
      }
      if (k === 'fermions') {
        return svg(`<line x1="${W / 2 - 80}" x2="${W / 2 - 20}" y1="26" y2="26" stroke="var(--xi)" stroke-width="2"/><line x1="${W / 2 + 20}" x2="${W / 2 + 80}" y1="26" y2="26" stroke="var(--xi)" stroke-width="2"/><text x="${W / 2 - 50}" y="46" text-anchor="middle">k</text><text x="${W / 2 + 50}" y="46" text-anchor="middle">−k</text><text x="${W / 2}" y="66" text-anchor="middle">a pair level of two fermionic modes</text>`);
      }
      return '';
    }
    if (rec.family === 'dissipative') {
      const vars = e.variables, n = vars.length;
      if (e.kind_dynamics === 'rotor') return svg(`<text x="${W / 2}" y="40" text-anchor="middle">${n} rotors: angles with period ${fmt(e.carrier.period)}</text>`);
      const pos = vars.map((_, i) => n === 1 ? [W / 2, 30] : n === 2 ? [W / 2 - 60 + 120 * i, 30] : [W / 2 + 70 * Math.cos(2 * Math.PI * i / n - Math.PI / 2) * 1.3, 32 + 22 * Math.sin(2 * Math.PI * i / n - Math.PI / 2)]);
      const kind = e.carrier.kind === 'torus' ? 'angle' + (n > 1 ? 's' : '') : e.carrier.kind === 'orthant' ? (n > 1 ? 'non-negative amounts' : 'a non-negative amount') : (n > 1 ? 'real variables' : 'a real variable');
      return svg(pos.map(([x, y], i) => `<circle cx="${x}" cy="${y}" r="14" fill="none" stroke="var(--xi)" stroke-width="1.6"/><text x="${x}" y="${y + 4}" text-anchor="middle">${esc(vars[i].replace(/_/g, ''))}</text>`).join('')
        + `<text x="${W / 2}" y="${n > 2 ? 69 : 64}" text-anchor="middle">${kind}</text>`);
    }
    return '';
  }

  root.FieldBridgeViews = {fit, frame, line, dot, label, arrow, sphere, signal, landscape, bifurcation, phasePlane, series,
                           hysteresis, circle, rotors, loglog, profile, phaseDrift, carrier, autoRange, fmt, pm, palette};
})(window);
