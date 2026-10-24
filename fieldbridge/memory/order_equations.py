"""The order of writes on a network given by its equations: the material of `order` replaced by drift equations.

The material of `order` is one coordinate per element in a fixed double well, written by force pulses. A network of
physical elements (Kerr parametric oscillators, parametric oscillators of a coherent Ising machine, rf-SQUID loops) is
given instead by its drift equations, a specification of kind ``equations`` in which every node has its own
variables, a readout variable, a control parameter that makes it bistable and a bias parameter that writes it; the
couplings are the parameters that the uncoupled reference sets to zero. Two kinds of write are read:

  ramp    the control of the write's support is taken from its rest value to its working value over `ramp`, with the
          bias of the write's sign on during the ramp and for `hold` after it, then off for `settle` (a pump ramp,
          the anneal of one qubit, a gain ramp: initialization or readout one node after another). A node written
          before is first returned to rest for `reset`; a node never written stays at rest and holds no state.
  pulse   every control is at its working value from the start and the network is prepared in the state of sign
          `initial`; a write applies the bias of its sign to its support for `pulse`, then `settle` without it
          (a flux pulse or a drive pulse on a bistable node).

After the last write the biases are off for `release` and the state is certified as a stable fixed point of the full
equations (Newton on the drift, every eigenvalue of the Jacobian with a negative real part). Its pattern is the sign
of the readout of every node whose control is at its working value; a node at rest reads 0. The graph of the
commutation index is read from the equations: two nodes are joined when the drift of one depends on the variables of
the other.

With `noise` the same protocol is integrated with the noise of the network specification (stochastic Heun, an
increment sqrt(2 D dt) on every variable or on the listed ones). The outcome of an order is then a distribution of
patterns, and the record is the capacity of the channel from the order to the pattern (log2 D where every outcome is
certain), with a permutation test of the hypothesis that the patterns do not depend on the order.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from .spec import SpecError

MAX_WRITES = 6
MODES = ("ramp", "pulse")
DEFAULTS = {"ramp": 20.0, "hold": 10.0, "settle": 5.0, "reset": 0.0, "pulse": 5.0, "release": 50.0}


@dataclass
class Write:
    name: str
    support: np.ndarray          # 0/1 mask over the nodes
    sign: float
    amplitude: float


@dataclass
class Node:
    variables: List[int]         # indices of the node's variables in the state vector
    readout: int
    control: str
    bias: str


def _expand(entry, n: int) -> List[Dict]:
    """The node list of a specification: a list of objects, or one object with {i} in its names and a count."""
    if isinstance(entry, list):
        return entry
    if isinstance(entry, dict) and "count" in entry:
        def sub(v, i):
            return [s.replace("{i}", str(i)) for s in v] if isinstance(v, list) else str(v).replace("{i}", str(i))
        return [{k: sub(v, i) for k, v in entry.items() if k != "count"} for i in range(int(entry["count"]))]
    raise SpecError("nodes must be a list of {variables, readout, control, bias} or one such object with {i} and a count")


class EquationNetwork:
    """A network of nodes given by drift equations, with the writes of an order specification."""

    def __init__(self, name: str, real, nodes: List[Node], coupling: Sequence[str], protocol: Dict, writes: List[Write],
                 couplings: str, noise: Optional[Dict] = None, zero_coupling: bool = False):
        self.name, self.real, self.nodes, self.coupling = name, real, nodes, list(coupling)
        self.protocol, self.writes, self.couplings, self.noise_cfg = protocol, writes, couplings, noise
        self.zero_coupling = zero_coupling
        self.n = len(nodes)
        self.dim = real.carrier.dim
        self.base = dict(real.params)
        if zero_coupling:
            for c in self.coupling:
                self.base[c] = 0.0
        self.edges = self._edges()
        self.signs = [{"attractive": 1.0, "antiferromagnetic": -1.0}.get(couplings, 0.0)] * len(self.edges)
        self._jac = self._jacobian()

    # structure
    def _edges(self):
        sym = self.real.symbolic
        owner = {}
        for k, node in enumerate(self.nodes):
            for v in node.variables:
                owner[v] = k
        found = set()
        for i, expr in enumerate(sym["drift"]):
            for v, s in enumerate(sym["variables"]):
                a, b = owner.get(i), owner.get(v)
                if a is None or b is None or a == b:
                    continue
                if expr.has(s):
                    found.add((min(a, b), max(a, b)))
        return sorted(found)

    def _jacobian(self):
        import sympy as sp
        sym = self.real.symbolic
        J = sp.Matrix(sym["drift"]).jacobian(sp.Matrix(sym["variables"]))
        fn = sp.lambdify(list(sym["variables"]) + list(sym["parameters"]), J, modules="numpy")
        names = [str(p) for p in sym["parameters"]]
        return lambda q, p: np.asarray(fn(*q, *[p[k] for k in names]), float)

    def uncoupled(self) -> "EquationNetwork":
        return EquationNetwork(f"{self.name} (uncoupled reference)", self.real, self.nodes, self.coupling,
                               self.protocol, self.writes, self.couplings, self.noise_cfg, zero_coupling=True)

    # protocol
    def _params(self, control: Dict[int, float], bias: Dict[int, float]) -> Dict[str, float]:
        p = dict(self.base)
        for k, node in enumerate(self.nodes):
            p[node.control] = control.get(k, p[node.control])
            p[node.bias] = bias.get(k, 0.0)
        return p

    def segments(self, names: Sequence[str]):
        """The protocol of one order as (duration, start controls, end controls, biases), and the nodes at work."""
        pr, by = self.protocol, {w.name: w for w in self.writes}
        rest, on = pr["rest"], pr["on"]
        segs = []
        if pr["mode"] == "ramp":
            control = {k: rest for k in range(self.n)}
            written = set()
            for name in names:
                w = by[name]
                S = [k for k in range(self.n) if w.support[k] > 0]
                b = {k: w.sign * w.amplitude for k in S}
                if pr["reset"] > 0 and written & set(S):
                    down = {**control, **{k: rest for k in S}}
                    segs.append((pr["reset"], down, down, {}))
                    control = down
                up = {**control, **{k: on for k in S}}
                segs.append((pr["ramp"], control, up, b))
                control = up
                if pr["hold"] > 0:
                    segs.append((pr["hold"], control, control, b))
                if pr["settle"] > 0:
                    segs.append((pr["settle"], control, control, {}))
                written |= set(S)
        else:
            control = {k: on for k in range(self.n)}
            for name in names:
                w = by[name]
                segs.append((pr["pulse"], control, control, {k: w.sign * w.amplitude for k in range(self.n) if w.support[k] > 0}))
                if pr["settle"] > 0:
                    segs.append((pr["settle"], control, control, {}))
        segs.append((pr["release"], control, control, {}))
        active = [k for k in range(self.n) if abs(control[k] - on) < 1e-12]
        return segs, control, active

    def _at(self, seg, s: float) -> Dict[str, float]:
        _, c0, c1, b = seg
        return self._params({k: c0[k] + s * (c1[k] - c0[k]) for k in c0}, b)

    def initial(self) -> np.ndarray:
        """The rest state (ramp writes) or the state of sign `initial` with every control at work (pulse writes)."""
        pr = self.protocol
        if pr["mode"] == "ramp":
            seg = (pr["release"], {k: pr["rest"] for k in range(self.n)}, {k: pr["rest"] for k in range(self.n)}, {})
            return self._advance(np.zeros(self.dim), seg)
        on = {k: pr["on"] for k in range(self.n)}
        amp = max(w.amplitude for w in self.writes)
        q = self._advance(np.zeros(self.dim), (pr["pulse"] * 4, on, on, {k: pr["initial"] * amp for k in range(self.n)}))
        return self._advance(q, (pr["release"], on, on, {}))

    # deterministic dynamics
    def F(self, q: np.ndarray, p: Dict[str, float]) -> np.ndarray:
        return self.real.drift(np.atleast_2d(q), p)

    def _advance(self, q: np.ndarray, seg) -> np.ndarray:
        from scipy.integrate import solve_ivp
        T = seg[0]
        sol = solve_ivp(lambda t, y: self.F(y, self._at(seg, t / T))[0], (0.0, T), q, method="LSODA", rtol=1e-9,
                        atol=1e-11)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y[:, -1]

    def certify(self, q: np.ndarray, p: Dict[str, float], max_time: float = 4000.0) -> Optional[np.ndarray]:
        """The stable fixed point the state relaxes to, or None if it is not reached within max_time."""
        from scipy.integrate import solve_ivp
        from scipy.optimize import root
        elapsed = 0.0
        while elapsed <= max_time:
            if np.max(np.abs(self.F(q, p)[0])) < 1e-7:
                sol = root(lambda x: self.F(x, p)[0], q, jac=lambda x: self._jac(x, p), tol=1e-13)
                eq = sol.x
                lead = float(np.max(np.linalg.eigvals(self._jac(eq, p)).real))
                if np.max(np.abs(self.F(eq, p)[0])) < 1e-10 and lead < -1e-9 and np.max(np.abs(eq - q)) < 1e-4:
                    return eq
            sol = solve_ivp(lambda t, y: self.F(y, p)[0], (0.0, 200.0), q, method="LSODA", rtol=1e-10, atol=1e-12)
            q = sol.y[:, -1]
            elapsed += 200.0
        return None

    def read(self, q: np.ndarray, active: Sequence[int]) -> str:
        return "".join(("+" if q[nd.readout] > 0 else "-") if k in active else "0" for k, nd in enumerate(self.nodes))

    def pattern(self, names: Sequence[str]) -> Optional[str]:
        """The certified pattern after the writes in the given order and the release."""
        segs, control, active = self.segments(names)
        q = self.initial()
        for seg in segs:
            q = self._advance(q, seg)
        eq = self.certify(q, self._at(segs[-1], 1.0))
        return None if eq is None else self.read(eq, active)

    # stochastic dynamics
    def sample(self, names: Sequence[str], trajectories: int, rng) -> List[str]:
        """Patterns of the given number of noisy trajectories (stochastic Heun with the network's noise)."""
        cfg = self.noise_cfg or {}
        dt = float(cfg.get("dt", self.real.dt or 0.005))
        D = float(cfg.get("D", self.real.noise))
        mask = np.ones(self.dim)
        if cfg.get("variables"):
            mask = np.array([1.0 if v in cfg["variables"] else 0.0 for v in self.real.variables])
        amp = np.sqrt(2.0 * D * dt) * mask
        segs, control, active = self.segments(names)
        q = np.tile(self.initial(), (trajectories, 1))
        for seg in segs:
            steps = max(1, int(np.ceil(seg[0] / dt)))
            h = seg[0] / steps
            for s in range(steps):
                p0, p1 = self._at(seg, s / steps), self._at(seg, (s + 1) / steps)
                dW = amp * rng.standard_normal(q.shape) * np.sqrt(h / dt)
                f0 = self.F(q, p0)
                qp = q + f0 * h + dW
                q = q + 0.5 * (f0 + self.F(qp, p1)) * h + dW
        bad = ~np.all(np.isfinite(q), axis=1)
        return [None if b else self.read(row, active) for row, b in zip(q, bad)]


# ------------------------------------------------------------------------------------------------ specification
def from_spec(spec: Dict, base: Optional[Path] = None) -> EquationNetwork:
    """The equation network of an order specification with material ``equations``; SpecError says what to change."""
    from .spec import load
    o = spec["order"]
    net_spec = o.get("network")
    if isinstance(net_spec, str):
        path = Path(net_spec) if base is None or Path(net_spec).is_absolute() else Path(base) / net_spec
        try:
            import json
            net_spec = json.loads(path.read_text(encoding="utf-8"))
        except OSError as error:
            raise SpecError(f"network: cannot read {path}") from error
    if not isinstance(net_spec, dict) or net_spec.get("kind") != "equations":
        raise SpecError("network must be a specification of kind 'equations' (inline or a path)")
    real = load(net_spec)
    if not hasattr(real, "symbolic"):
        raise SpecError("network must be a specification of kind 'equations'")
    var_index = {v: k for k, v in enumerate(real.variables)}
    nodes = []
    for k, nd in enumerate(_expand(o.get("nodes"), 0)):
        try:
            variables = [var_index[v] for v in nd["variables"]]
            readout = var_index[nd["readout"]]
        except KeyError as error:
            raise SpecError(f"node {k}: {error.args[0]!r} is not a variable of the network") from error
        for key in ("control", "bias"):
            if nd.get(key) not in real.params:
                raise SpecError(f"node {k}: {key} {nd.get(key)!r} must be a parameter of the network")
        if readout not in variables:
            raise SpecError(f"node {k}: the readout must be one of the node's variables")
        nodes.append(Node(variables=variables, readout=readout, control=nd["control"], bias=nd["bias"]))
    n = len(nodes)
    if n < 2:
        raise SpecError("an equation network needs at least two nodes")
    coupling = o.get("coupling", [])
    for c in coupling:
        if c not in real.params:
            raise SpecError(f"coupling {c!r} must be a parameter of the network")
    w = o.get("write")
    if not isinstance(w, dict) or w.get("mode") not in MODES:
        raise SpecError("write must be an object with mode 'ramp' or 'pulse'")
    protocol = {"mode": w["mode"], **{k: float(w.get(k, v)) for k, v in DEFAULTS.items()}}
    for key in ("on",) + (("rest",) if w["mode"] == "ramp" else ()):
        if key not in w:
            raise SpecError(f"write: {key!r} (the control value {'at work' if key == 'on' else 'at rest'}) is required")
    protocol["on"] = float(w["on"])
    protocol["rest"] = float(w.get("rest", w["on"]))
    protocol["initial"] = float(w.get("initial", -1.0))
    amplitude = float(w.get("bias", 0.1))
    writes_spec = o.get("writes")
    if not isinstance(writes_spec, dict) or not 2 <= len(writes_spec) <= MAX_WRITES:
        raise SpecError(f"writes must be an object of two to {MAX_WRITES} named writes")
    writes = []
    for name, ws in writes_spec.items():
        ws = ws or {}
        mask = np.zeros(n)
        for v in ws.get("support", []):
            if not 0 <= int(v) < n:
                raise SpecError(f"write {name}: support nodes are numbered 0 to {n - 1}")
            mask[int(v)] = 1.0
        if not mask.any():
            raise SpecError(f"write {name}: a support of at least one node is required")
        writes.append(Write(name=str(name), support=mask, sign=float(ws.get("sign", 1.0)),
                            amplitude=float(ws.get("amplitude", amplitude))))
    couplings = o.get("couplings", "mixed")
    if couplings not in ("attractive", "antiferromagnetic", "mixed"):
        raise SpecError("couplings must be 'attractive', 'antiferromagnetic' or 'mixed'")
    noise = o.get("noise")
    if noise is not None and not isinstance(noise, dict):
        raise SpecError("noise must be an object {trajectories, seed, dt?, variables?}")
    return EquationNetwork(str(spec.get("name", real.name)), real, nodes, coupling, protocol, writes, couplings, noise)


# ------------------------------------------------------------------------------------------------ noisy record
def capacity(P: np.ndarray, tol: float = 1e-10, max_iter: int = 5000) -> float:
    """Capacity in bits of the channel with rows P[input] = distribution of outputs (Blahut-Arimoto)."""
    P = np.asarray(P, float)
    m = P.shape[0]
    r = np.full(m, 1.0 / m)
    with np.errstate(divide="ignore", invalid="ignore"):
        logP = np.where(P > 0, np.log(P), 0.0)
    for _ in range(max_iter):
        q = r @ P
        with np.errstate(divide="ignore", invalid="ignore"):
            d = np.sum(np.where(P > 0, P * (logP - np.log(np.where(q > 0, q, 1.0))), 0.0), axis=1)
        c = np.exp(d)
        lo, hi = np.log(r @ c), np.log(np.max(c))
        r = r * c / (r @ c)
        if hi - lo < tol:
            break
    return float(lo / np.log(2.0))


def noisy_measure(net: EquationNetwork, trajectories: int, seed: int, permutations: int = 200) -> Dict[str, object]:
    """Pattern counts for every order, the capacity of the order-to-pattern channel and its permutation null."""
    rng = np.random.default_rng(seed)
    names = [w.name for w in net.writes]
    orders = ["".join(p) for p in itertools.permutations(names)]
    samples = {o: net.sample(list(o), trajectories, rng) for o in orders}
    labels = sorted({s for v in samples.values() for s in v if s is not None})
    idx = {s: k for k, s in enumerate(labels)}
    counts = np.zeros((len(orders), max(len(labels), 1)))
    lost = 0
    for a, o in enumerate(orders):
        for s in samples[o]:
            if s is None:
                lost += 1
            else:
                counts[a, idx[s]] += 1
    P = counts / np.maximum(counts.sum(axis=1, keepdims=True), 1)
    C = capacity(P)
    pooled = np.concatenate([np.repeat(np.arange(counts.shape[1]), counts[a].astype(int)) for a in range(len(orders))])
    sizes = counts.sum(axis=1).astype(int)
    null = []
    for _ in range(permutations):
        perm = rng.permutation(pooled)
        Q, start = np.zeros_like(counts), 0
        for a, m in enumerate(sizes):
            Q[a] = np.bincount(perm[start:start + m], minlength=counts.shape[1])
            start += m
        null.append(capacity(Q / np.maximum(Q.sum(axis=1, keepdims=True), 1)))
    null = np.array(null)
    tv = max((0.5 * np.abs(P[a] - P[b]).sum() for a in range(len(orders)) for b in range(a + 1, len(orders))), default=0.0)
    return {"trajectories": trajectories, "seed": seed, "lost": lost,
            "patterns": {o: {labels[k]: int(counts[a, k]) for k in range(len(labels)) if counts[a, k]} for a, o in enumerate(orders)},
            "capacity": C, "null_capacity_mean": float(null.mean()), "null_capacity_95": float(np.quantile(null, 0.95)),
            "p_value": float((1 + np.sum(null >= C)) / (1 + len(null))), "excess": float(C - null.mean()),
            "max_total_variation": float(tv)}
