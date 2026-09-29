/* Memory realizations in the browser: states, their stability, writes and drives.
 *
 * Built on model-physics.js (the drift of a specification as an arithmetic tree, integrated by RK4). The same
 * quantities as fieldbridge/memory/analysis.py, for the parameters the visitor sets:
 *   equilibria    zeros of the drift by damped Newton from seeds and random starts, with the eigenvalues of the
 *                 Jacobian (kappa = -eigenvalues: > 0 relaxing, < 0 unstable)
 *   attractors    where many preparations end: stable states, or a cycle when they do not settle
 *   scan          equilibria along the control parameter: the bifurcation diagram
 *   driven        integration with a parameter modulated as p1 + eps cos(omega t) (phase_locking._solve)
 *   langevin      Euler-Maruyama with a seeded generator, for retention under noise
 * Tested against fieldbridge.memory in tests/test_site_engines.py.
 */
(function (root, factory) {
  const models = typeof module === 'object' && module.exports ? require('./model-physics.js') : root.FieldBridgeModels;
  const api = factory(models);
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.FieldBridgeDissipative = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function (M) {
  'use strict';

  function field(model, params, removed) { return M.evaluator(model, params || model.params, removed || []); }
  function wrap(model, q) {
    const c = model.carrier;
    if (c.kind === 'orthant') return q.map(Math.abs);
    if (c.kind === 'torus') return q.map(x => ((x % c.period) + c.period) % c.period);
    return q;
  }
  function diff(model, a, b) {
    const c = model.carrier;
    return b.map((x, i) => { let d = x - a[i]; if (c.kind === 'torus') d = ((d + c.period / 2) % c.period + c.period) % c.period - c.period / 2; return d; });
  }
  const dist = (model, a, b) => Math.hypot(...diff(model, a, b));
  function jacobian(f, q, eps = 1e-6) {
    const n = q.length, J = Array.from({length: n}, () => new Array(n).fill(0));
    for (let j = 0; j < n; j++) {
      const qp = q.slice(), qm = q.slice(); qp[j] += eps; qm[j] -= eps;
      const fp = f(qp), fm = f(qm);
      for (let i = 0; i < n; i++) J[i][j] = (fp[i] - fm[i]) / (2 * eps);
    }
    return J;
  }

  // ---------------------------------------------------------------------------------------- eigenvalues
  function symmetricEigenvalues(A0) {
    const n = A0.length, A = A0.map(r => r.slice());
    for (let sweep = 0; sweep < 80; sweep++) {
      let off = 0; for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += A[p][q] ** 2;
      if (off < 1e-26) break;
      for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-300) continue;
        const tau = (A[q][q] - A[p][p]) / (2 * A[p][q]), t = (tau >= 0 ? 1 : -1) / (Math.abs(tau) + Math.sqrt(1 + tau * tau));
        const c = 1 / Math.sqrt(1 + t * t), s = t * c;
        for (let k = 0; k < n; k++) { const a = A[k][p], b = A[k][q]; A[k][p] = c * a - s * b; A[k][q] = s * a + c * b; }
        for (let k = 0; k < n; k++) { const a = A[p][k], b = A[q][k]; A[p][k] = c * a - s * b; A[q][k] = s * a + c * b; }
      }
    }
    return A.map((r, i) => ({re: r[i], im: 0}));
  }
  // Balancing, reduction to Hessenberg form and the shifted QR algorithm (Numerical Recipes: balanc, elmhes, hqr),
  // on 1-based arrays as in the reference.
  function balance(a, n) {
    const RADIX = 2, sq = RADIX * RADIX;
    let last = false;
    while (!last) {
      last = true;
      for (let i = 1; i <= n; i++) {
        let r = 0, c = 0;
        for (let j = 1; j <= n; j++) if (j !== i) { c += Math.abs(a[j][i]); r += Math.abs(a[i][j]); }
        if (c && r) {
          let g = r / RADIX, f = 1; const s0 = c + r;
          while (c < g) { f *= RADIX; c *= sq; }
          g = r * RADIX;
          while (c > g) { f /= RADIX; c /= sq; }
          if ((c + r) / f < 0.95 * s0) { last = false; g = 1 / f; for (let j = 1; j <= n; j++) a[i][j] *= g; for (let j = 1; j <= n; j++) a[j][i] *= f; }
        }
      }
    }
  }
  function hessenberg(a, n) {
    for (let m = 2; m < n; m++) {
      let x = 0, i = m;
      for (let j = m; j <= n; j++) if (Math.abs(a[j][m - 1]) > Math.abs(x)) { x = a[j][m - 1]; i = j; }
      if (i !== m) {
        for (let j = m - 1; j <= n; j++) [a[i][j], a[m][j]] = [a[m][j], a[i][j]];
        for (let j = 1; j <= n; j++) [a[j][i], a[j][m]] = [a[j][m], a[j][i]];
      }
      if (x) for (i = m + 1; i <= n; i++) {
        let y = a[i][m - 1];
        if (y !== 0) { y /= x; a[i][m - 1] = y; for (let j = m; j <= n; j++) a[i][j] -= y * a[m][j]; for (let j = 1; j <= n; j++) a[j][m] += y * a[j][i]; }
      }
    }
  }
  function hqr(a, n) {
    const wr = new Array(n + 1).fill(0), wi = new Array(n + 1).fill(0);
    let anorm = 0;
    for (let i = 1; i <= n; i++) for (let j = Math.max(i - 1, 1); j <= n; j++) anorm += Math.abs(a[i][j]);
    let nn = n, t = 0;
    const sign = (x, y) => (y >= 0 ? Math.abs(x) : -Math.abs(x));
    while (nn >= 1) {
      let its = 0, l;
      do {
        for (l = nn; l >= 2; l--) {
          let s = Math.abs(a[l - 1][l - 1]) + Math.abs(a[l][l]);
          if (s === 0) s = anorm;
          if (Math.abs(a[l][l - 1]) + s === s) { a[l][l - 1] = 0; break; }
        }
        let x = a[nn][nn];
        if (l === nn) { wr[nn] = x + t; wi[nn--] = 0; }
        else {
          let y = a[nn - 1][nn - 1], w = a[nn][nn - 1] * a[nn - 1][nn];
          if (l === nn - 1) {
            const p = 0.5 * (y - x), q = p * p + w;
            let z = Math.sqrt(Math.abs(q));
            x += t;
            if (q >= 0) { z = p + sign(z, p); wr[nn - 1] = wr[nn] = x + z; if (z) wr[nn] = x - w / z; wi[nn - 1] = wi[nn] = 0; }
            else { wr[nn - 1] = wr[nn] = x + p; wi[nn - 1] = -(wi[nn] = z); }
            nn -= 2;
          } else {
            if (its === 60) throw new Error('eigenvalues: no convergence');
            if (its === 10 || its === 20) {
              t += x;
              for (let i = 1; i <= nn; i++) a[i][i] -= x;
              const s = Math.abs(a[nn][nn - 1]) + Math.abs(a[nn - 1][nn - 2]);
              y = x = 0.75 * s; w = -0.4375 * s * s;
            }
            ++its;
            let m, p, q, r, z, s;
            for (m = nn - 2; m >= l; m--) {
              z = a[m][m]; r = x - z; s = y - z;
              p = (r * s - w) / a[m + 1][m] + a[m][m + 1]; q = a[m + 1][m + 1] - z - r - s; r = a[m + 2][m + 1];
              s = Math.abs(p) + Math.abs(q) + Math.abs(r); p /= s; q /= s; r /= s;
              if (m === l) break;
              const u = Math.abs(a[m][m - 1]) * (Math.abs(q) + Math.abs(r));
              const v = Math.abs(p) * (Math.abs(a[m - 1][m - 1]) + Math.abs(z) + Math.abs(a[m + 1][m + 1]));
              if (u + v === v) break;
            }
            for (let i = m + 2; i <= nn; i++) { a[i][i - 2] = 0; if (i !== m + 2) a[i][i - 3] = 0; }
            for (let k = m; k <= nn - 1; k++) {
              if (k !== m) {
                p = a[k][k - 1]; q = a[k + 1][k - 1]; r = 0;
                if (k !== nn - 1) r = a[k + 2][k - 1];
                if ((x = Math.abs(p) + Math.abs(q) + Math.abs(r)) !== 0) { p /= x; q /= x; r /= x; }
              }
              if ((s = sign(Math.sqrt(p * p + q * q + r * r), p)) !== 0) {
                if (k === m) { if (l !== m) a[k][k - 1] = -a[k][k - 1]; } else a[k][k - 1] = -s * x;
                p += s; x = p / s; y = q / s; z = r / s; q /= p; r /= p;
                for (let j = k; j <= nn; j++) {
                  p = a[k][j] + q * a[k + 1][j];
                  if (k !== nn - 1) { p += r * a[k + 2][j]; a[k + 2][j] -= p * z; }
                  a[k + 1][j] -= p * y; a[k][j] -= p * x;
                }
                const mmin = nn < k + 3 ? nn : k + 3;
                for (let i = l; i <= mmin; i++) {
                  p = x * a[i][k] + y * a[i][k + 1];
                  if (k !== nn - 1) { p += z * a[i][k + 2]; a[i][k + 2] -= p * r; }
                  a[i][k + 1] -= p * q; a[i][k] -= p;
                }
              }
            }
          }
        }
      } while (l < nn - 1);
    }
    return Array.from({length: n}, (_, i) => ({re: wr[i + 1], im: wi[i + 1]}));
  }
  /** Eigenvalues of a real matrix: exact for n <= 2, Jacobi for symmetric matrices (gradient drifts), otherwise the
   * shifted QR algorithm on the balanced Hessenberg form. */
  function eigenvalues(J) {
    const n = J.length;
    if (n === 1) return [{re: J[0][0], im: 0}];
    if (n === 2) {
      const tr = J[0][0] + J[1][1], det = J[0][0] * J[1][1] - J[0][1] * J[1][0], d = tr * tr / 4 - det;
      return d >= 0 ? [{re: tr / 2 - Math.sqrt(d), im: 0}, {re: tr / 2 + Math.sqrt(d), im: 0}]
                    : [{re: tr / 2, im: -Math.sqrt(-d)}, {re: tr / 2, im: Math.sqrt(-d)}];
    }
    let asym = 0, scale = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { asym = Math.max(asym, Math.abs(J[i][j] - J[j][i])); scale = Math.max(scale, Math.abs(J[i][j])); }
    if (asym <= 1e-9 * Math.max(scale, 1)) return symmetricEigenvalues(J);
    const a = Array.from({length: n + 1}, (_, i) => Array.from({length: n + 1}, (_, j) => (i && j ? J[i - 1][j - 1] : 0)));
    balance(a, n); hessenberg(a, n);
    return hqr(a, n).sort((x, y) => x.re - y.re || x.im - y.im);
  }
  function classify(eigs, tol = 1e-7) {
    const scale = Math.max(1, ...eigs.map(e => Math.abs(e.re)));
    return {nUnstable: eigs.filter(e => e.re > tol * scale).length, nFlat: eigs.filter(e => Math.abs(e.re) <= tol * scale).length,
            nStable: eigs.filter(e => e.re < -tol * scale).length, rotating: eigs.some(e => Math.abs(e.im) > 1e-6)};
  }

  // ---------------------------------------------------------------------------------------- equilibria
  function solve(A, b) {   // Gaussian elimination with partial pivoting; least-norm fallback is not needed here
    const n = b.length, M2 = A.map((r, i) => [...r, b[i]]);
    for (let k = 0; k < n; k++) {
      let piv = k; for (let i = k + 1; i < n; i++) if (Math.abs(M2[i][k]) > Math.abs(M2[piv][k])) piv = i;
      if (Math.abs(M2[piv][k]) < 1e-14) return null;
      [M2[k], M2[piv]] = [M2[piv], M2[k]];
      for (let i = k + 1; i < n; i++) { const f = M2[i][k] / M2[k][k]; for (let j = k; j <= n; j++) M2[i][j] -= f * M2[k][j]; }
    }
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) { let s = M2[i][n]; for (let j = i + 1; j < n; j++) s -= M2[i][j] * x[j]; x[i] = s / M2[i][i]; }
    return x;
  }
  const baseLength = model => model.carrier.kind === 'torus' ? model.carrier.period : Math.max(1, model.carrier.scale || 1);
  function newton(model, f, q0, iters = 80) {
    let q = wrap(model, q0.slice()); const cap = 0.5 * baseLength(model);
    for (let it = 0; it < iters; it++) {
      const F = f(q), nf = Math.hypot(...F);
      if (!Number.isFinite(nf)) return {q, ok: false};
      if (nf < 1e-11) return {q, ok: true};
      let step = solve(jacobian(f, q), F);
      if (!step) return {q, ok: false};
      const ns = Math.hypot(...step); if (ns > cap) step = step.map(s => s * cap / ns);
      q = wrap(model, q.map((x, i) => x - step[i]));
    }
    return {q, ok: Math.hypot(...f(q)) < 1e-8};
  }
  function rng(seed) { let s = seed >>> 0 || 1; return () => { s ^= s << 13; s >>>= 0; s ^= s >> 17; s ^= s << 5; s >>>= 0; return s / 4294967296; }; }
  function sample(model, n, random) {
    const c = model.carrier, d = model.variables.length;
    return Array.from({length: n}, () => Array.from({length: d}, () => {
      if (c.kind === 'torus') return random() * c.period;
      if (c.kind === 'orthant') return random() < 0.5 ? random() * (c.scale || 1) : (c.scale || 1) * Math.pow(10, -3 * random());
      return (2 * random() - 1) * (c.scale || 1);
    }));
  }
  function equilibria(model, params, removed, seeds = [], n = 24, seed = 1) {
    const f = field(model, params, removed), found = [], tol = 1e-4 * baseLength(model);
    for (const q0 of [...seeds, ...sample(model, n, rng(seed))]) {
      const {q, ok} = newton(model, f, q0);
      if (!ok || found.some(e => dist(model, e.q, q) < tol)) continue;
      const eig = eigenvalues(jacobian(f, q)), cls = classify(eig);
      found.push({q, eig, ...cls, stable: cls.nUnstable === 0 && cls.nFlat === 0});
    }
    return found;
  }

  // ---------------------------------------------------------------------------------------- attractors
  function relax(model, f, q0, T, dt) {
    // at most 4000 steps: this only finds where a preparation ends, and the page recomputes it at every change
    let q = q0.slice(); const step = Math.max(Math.min(dt, 0.05), T / 4000);
    const rk = (x) => f(model.carrier.kind === 'orthant' ? x.map(v => Math.max(0, v)) : x);
    for (let t = 0; t < T; t += step) {
      const a = rk(q), b = rk(q.map((x, i) => x + step * a[i] / 2)), c = rk(q.map((x, i) => x + step * b[i] / 2)), d = rk(q.map((x, i) => x + step * c[i]));
      q = q.map((x, i) => x + step * (a[i] + 2 * b[i] + 2 * c[i] + d[i]) / 6);
      if (model.carrier.kind === 'orthant') q = q.map(v => Math.max(0, v));
      if (q.some(v => !Number.isFinite(v) || Math.abs(v) > 1e6)) return {q, speed: Infinity};
    }
    return {q: wrap(model, q), speed: Math.hypot(...rk(q))};
  }
  /** Where preparations end: the stable states reached and whether some do not settle (a cycle, or drift). */
  function attractors(model, params, removed, opts = {}) {
    const f = field(model, params, removed), n = opts.n || (model.variables.length > 4 ? 8 : 16), T = opts.T || 80, dt = opts.dt || model.dt || 0.01;
    const ends = sample(model, n, rng(opts.seed || 7)).map(q0 => relax(model, f, q0, T, dt));
    const points = [], tol = 1e-3 * baseLength(model);
    let moving = 0, unbounded = 0;
    for (const e of ends) {
      if (!Number.isFinite(e.speed)) { unbounded++; continue; }
      if (e.speed > 1e-4) { moving++; continue; }
      const {q, ok} = newton(model, f, e.q);
      if (!ok) { moving++; continue; }
      const hit = points.find(p => dist(model, p.q, q) < tol);
      if (hit) hit.count++; else points.push({q, count: 1, eig: eigenvalues(jacobian(f, q))});
    }
    return {points: points.filter(p => classify(p.eig).nUnstable === 0), moving, unbounded, n};
  }
  function label(att) {
    const k = att.points.length;
    if (att.unbounded === att.n) return 'unbounded';
    if (k === 0 && att.moving) return 'no stable state: the preparations keep moving on a cycle';
    if (k === 1) return att.moving ? 'one stable state; some preparations keep moving' : 'one stable state';
    return `${k} stable states` + (att.moving ? '; some preparations keep moving' : '');
  }

  // ---------------------------------------------------------------------------------------- diagrams
  function summary(model, q) {
    if (q.length === 1) return model.carrier.kind === 'torus' ? Math.atan2(Math.sin(q[0] * 2 * Math.PI / model.carrier.period), Math.cos(q[0] * 2 * Math.PI / model.carrier.period)) : q[0];
    return q[0] - q.slice(1).reduce((s, x) => s + x, 0) / (q.length - 1);
  }
  /** Equilibria along a parameter; each value is seeded by the states of the previous one. */
  function scan(model, params, removed, name, values, n = 12) {
    let seeds = [];
    return values.map((v, k) => {
      const eq = equilibria(model, {...params, [name]: v}, removed, seeds, n, 11 + k);
      seeds = eq.map(e => e.q);
      return {v, points: eq.map(e => ({q: e.q, s: summary(model, e.q), stable: e.stable}))};
    });
  }
  /** V(x) = -integral of the drift for a one-variable realization, on a grid (trapezoids). */
  function potential1D(model, params, removed, xs) {
    const f = field(model, params, removed); let V = 0; const out = [0];
    for (let i = 1; i < xs.length; i++) { V -= 0.5 * (f([xs[i - 1]])[0] + f([xs[i]])[0]) * (xs[i] - xs[i - 1]); out.push(V); }
    return out;
  }

  // ---------------------------------------------------------------------------------------- drive, noise
  /** RK4 with a parameter modulated as p1 + eps cos(omega t); records every sample interval. */
  function driven(model, params, removed, drive, q0, T, dt = 0.01, every = 0.05) {
    let q = q0.slice(), t = 0, next = 0; const data = [];
    const at = time => M.evaluator(model, {...params, [drive.param]: drive.p1 + drive.eps * Math.cos(drive.omega * time)}, removed || []);
    while (t < T - 1e-12) {
      const h = Math.min(dt, T - t);
      const a = at(t)(q), b = at(t + h / 2)(q.map((x, i) => x + h * a[i] / 2)), c = at(t + h / 2)(q.map((x, i) => x + h * b[i] / 2)), d = at(t + h)(q.map((x, i) => x + h * c[i]));
      q = q.map((x, i) => x + h * (a[i] + 2 * b[i] + 2 * c[i] + d[i]) / 6); t += h;
      if (model.carrier.kind === 'orthant') q = q.map(v => Math.max(0, v));
      if (t >= next - 1e-12) { data.push([t, ...q]); next += every; }
      if (q.some(v => !Number.isFinite(v))) break;
    }
    return data;
  }
  /** Phase of the oscillation against the drive: at each upward crossing of the mean of a coordinate, the fraction
   * of the drive period (times the ratio) that has elapsed. Constant when locked, drifting when the phase slips. */
  function crossingPhases(data, axis, period, ratio = 1) {
    const vals = data.map(r => r[axis + 1]), mean = vals.reduce((s, v) => s + v, 0) / vals.length, out = [];
    for (let k = 1; k < data.length; k++) {
      if (vals[k - 1] < mean && vals[k] >= mean) {
        const t = data[k - 1][0] + (mean - vals[k - 1]) / (vals[k] - vals[k - 1]) * (data[k][0] - data[k - 1][0]);
        out.push([t, ((t / (period * ratio)) % 1 + 1) % 1]);
      }
    }
    return out;
  }
  function gauss(random) { let u = 0, v = 0; while (u === 0) u = random(); while (v === 0) v = random(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
  function langevin(model, params, removed, D, q0, T, dt = 0.01, every = 0.1, seed = 3) {
    const f = field(model, params, removed), random = rng(seed), amp = Math.sqrt(2 * D * dt), data = [[0, ...q0]];
    let q = q0.slice(), next = every;
    for (let t = dt; t <= T + 1e-12; t += dt) {
      const F = f(q);
      q = q.map((x, i) => x + dt * F[i] + amp * gauss(random));
      if (model.carrier.kind === 'orthant') q = q.map(Math.abs);
      if (t >= next) { data.push([t, ...wrap(model, q)]); next += every; }
    }
    return data;
  }

  return {field, wrap, diff, jacobian, eigenvalues, classify, newton, equilibria, attractors, label, summary, scan,
          potential1D, driven, crossingPhases, langevin, rng};
});
