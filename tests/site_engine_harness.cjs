// Runs one calculation of the browser engines on a JSON job, for tests/test_site_engines.py.
const fs = require('node:fs');
const path = require('node:path');
const web = path.resolve(__dirname, '..', 'fieldbridge', 'web');
const U = require(path.join(web, 'unitary.js'));
const D = require(path.join(web, 'dissipative.js'));
const F = require(path.join(web, 'fields.js'));
const MML = require(path.join(web, 'mathml.js'));
const job = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const vec = s => s ? {re: Float64Array.from(s.re), im: Float64Array.from(s.im)} : null;
const plain = M => ({n: M.n, re: Array.from(M.re), im: Array.from(M.im)});
const ops = {
  eigh: j => Array.from(U.eigh(U.fromJSON(j.matrix)).values),
  eigenvalues: j => D.eigenvalues(j.matrix),
  unitary_analyze: j => { const a = U.analyze(j.engine, j.params || j.engine.params, j.active); return a.frame ? {...a, frame: a.frame.map(plain)} : a; },
  unitary_axis: j => U.frameRotation(U.hamiltonian(j.engine, j.params || j.engine.params, j.active), j.engine.frame.map(U.fromJSON)),
  unitary_leak: j => U.frameLeak(U.hamiltonian(j.engine, j.params || j.engine.params, j.active), j.engine.frame.map(U.fromJSON)),
  unitary_signal: j => {
    const psi0 = vec(j.psi0) || U.topEigenstate(U.observable(j.engine));
    return U.signal(j.engine, j.params || j.engine.params, j.active, psi0, j.times);
  },
  equilibria: j => D.equilibria(j.model, j.params || j.model.params, [], j.seeds || [], j.n || 24, j.seed || 1)
    .map(e => ({q: e.q, stable: e.stable, eig: e.eig})),
  scan: j => D.scan(j.model, j.params || j.model.params, [], j.name, j.values, j.n || 12),
  attractors: j => D.attractors(j.model, j.params || j.model.params, [], j.opts || {}),
  driven_phase: j => {
    const data = D.driven(j.model, j.params || j.model.params, [], j.drive, j.q0, j.T, j.dt || 0.01, 0.01);
    return D.crossingPhases(data, j.axis || 0, 2 * Math.PI / j.drive.omega, j.ratio || 1);
  },
  fields_snr: j => F.snr(j.engine, j.times),
  fields_profile: j => F.profile(j.engine, j.t),
  mathml_tree: j => MML.tree(j.tree),
  mathml_term: j => MML.term(j.tree, j.operator, j.hc),
  // views.js is written for the browser; its number formats need no document
  value_with_uncertainty: j => {
    global.window = global.window || {};
    require(path.join(web, 'views.js'));
    return j.pairs.map(([v, err]) => window.FieldBridgeViews.pm(v, err));
  },
};
process.stdout.write(JSON.stringify(ops[job.op](job)));
