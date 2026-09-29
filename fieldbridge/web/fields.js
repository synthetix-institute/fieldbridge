/* Fields in the browser: the exact loss of a write in a linear lattice field (fieldbridge/memory/fields.py).
 *
 * On a ring of L^d sites a write of amplitude A at one site (a charge) or a pair of opposite writes at neighbouring
 * sites (a dipole) relaxes mode by mode, at kappa(k) = M lambda(k) for a conserved density (Model B) and
 * kappa(k) = M (kappa0 + D lambda(k)) without conservation (Model A), with lambda(k) = sum 2 (1 - cos k_i). The signal-
 * to-noise ratio of the best measurement of the write is SNR(t) = (1/N) sum_k |dphi_k|^2 exp(-2 kappa t) / sigma_k^2.
 */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.FieldBridgeFields = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  function modes(e) {
    const L = e.L, d = e.d, ks = Array.from({length: L}, (_, i) => 2 * Math.PI * i / L), out = [];
    if (d === 1) for (const kx of ks) out.push({kx, ky: 0, lam: 2 * (1 - Math.cos(kx))});
    else for (const kx of ks) for (const ky of ks) out.push({kx, ky, lam: 2 * (1 - Math.cos(kx)) + 2 * (1 - Math.cos(ky))});
    return out;
  }
  function rates(e, m) {
    const kappa0 = e.kappa0 ?? 0.2, D = e.Dgrad ?? 1;
    return e.conserved ? {rate: e.M * m.lam, variance: e.T} : {rate: e.M * (kappa0 + D * m.lam), variance: e.T / (kappa0 + D * m.lam)};
  }
  /** SNR of the profile: all modes except the conserved total, as fields.spectral_snr(modes="profile"). */
  function snr(e, times) {
    const ms = modes(e), N = ms.length, A2 = e.amplitude * e.amplitude;
    const terms = ms.filter(m => !(e.conserved && m.lam < 1e-12)).map(m => {
      const {rate, variance} = rates(e, m);
      return {rate, w: (e.write === 'charge' ? A2 : A2 * 2 * (1 - Math.cos(m.kx))) / variance};
    });
    return times.map(t => terms.reduce((s, m) => s + m.w * Math.exp(-2 * m.rate * t), 0) / N);
  }
  /** The mean trace of the write along x (at y = 0 in two dimensions): sum_k dphi_k exp(-kappa t) e^{i k x} / N. */
  function profile(e, t) {
    const ms = modes(e), N = ms.length;
    const w = ms.map(m => {
      const {rate} = rates(e, m), decay = Math.exp(-rate * t);
      // charge: A at site 0; dipole: A at site 0 and -A at site 1 along x
      const re = e.write === 'charge' ? e.amplitude : e.amplitude * (1 - Math.cos(m.kx));
      const im = e.write === 'charge' ? 0 : e.amplitude * Math.sin(m.kx);
      return {kx: m.kx, re: re * decay, im: im * decay};
    });
    return Array.from({length: e.L}, (_, x) => w.reduce((s, m) => s + m.re * Math.cos(m.kx * x) - m.im * Math.sin(m.kx * x), 0) / N);
  }
  /** Slope of ln SNR against ln t between t0 and t1. */
  function exponent(times, values, t0, t1) {
    const pts = times.map((t, i) => [Math.log(t), Math.log(values[i])]).filter(([lt, lv], i) => times[i] >= t0 && times[i] <= t1 && Number.isFinite(lv));
    const n = pts.length, mx = pts.reduce((s, p) => s + p[0], 0) / n, my = pts.reduce((s, p) => s + p[1], 0) / n;
    let sxy = 0, sxx = 0; for (const [x, y] of pts) { sxy += (x - mx) * (y - my); sxx += (x - mx) ** 2; }
    return sxy / sxx;
  }
  return {modes, snr, profile, exponent};
});
