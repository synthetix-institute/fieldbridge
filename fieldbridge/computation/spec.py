"""Driven bodies from specifications: an input held over intervals and observables measured in each.

Two schemas:

``fieldbridge-computation/1``
    an ``equations`` body as in the memory module (carrier, parameters, drift, read by the same restricted parser:
    nothing is evaluated with eval and no input is passed to sympify) with a block ``computation``;
``fieldbridge-springs/1``
    a network of nonlinear springs and point masses given as data (positions, springs, constants, fixed and input
    nodes), which written-out expressions cannot hold for networks of realistic size; see ``springs``.

The block ``computation``:

    input          a declared parameter, redrawn for every interval and held within it (slot P)
    law            "uniform" on [-A, A] (default), "gaussian" with deviation A, or "binary" (+-A, equal probability)
    amplitude      A
    offset         a constant added to the input (default 0), for example a bias that breaks a symmetry
    hold           the length of an interval
    observables    expressions in the variables and parameters (slot R); for a spring network "lengths"
    virtual_nodes  measurement times per interval, equally spaced, the last at its end (default 1)
    washout, length  intervals discarded and used per stream
    substeps       fourth-order Runge-Kutta steps per interval (default 20)
    solver         "rk4" (default) or "stiff" (backward differentiation per interval, for rates that differ by orders of
                   magnitude)
    form           "flow" (default: the drift is a rate) or "map" (the drift gives the next state, x_t = F(x_{t-1}, u_t))
    initial        a state from which the steady state at the offset is found
    exact_box      the range of each state variable for the exact capacities (default: from a driven run)
    expect         the class that a closed form or the source gives (benchmarks and controls)
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Dict, Tuple, Union

import numpy as np
from scipy.optimize import root

from ..memory import spec as memory_spec

SCHEMA = "fieldbridge-computation/1"
SPRINGS = "fieldbridge-springs/1"
SpecError = memory_spec.SpecError
LAWS = ("uniform", "gaussian", "binary")


def _read(source: Union[str, Path, Dict]) -> Tuple[Dict, str]:
    if isinstance(source, dict):
        return copy.deepcopy(source), "<dict>"
    path = Path(source)
    try:
        return json.loads(path.read_text(encoding="utf-8")), str(path)
    except (OSError, json.JSONDecodeError) as error:
        raise SpecError(f"Cannot read {path}: {error}") from error


class Body:
    """What every driven body provides: the protocol of the block ``computation`` and the drift, observation,
    steady state and simulation built on it."""

    def _protocol(self, comp: Dict) -> None:
        if not isinstance(comp, dict):
            raise SpecError("Missing block 'computation'")
        self.law = comp.get("law", "uniform")
        if self.law not in LAWS:
            raise SpecError(f"computation.law must be one of {', '.join(LAWS)}")
        self.amplitude = memory_spec._number(comp.get("amplitude", 1.0), "computation.amplitude")
        self.offset = memory_spec._number(comp.get("offset", 0.0), "computation.offset")
        self.hold = memory_spec._number(memory_spec._require(comp, "hold"), "computation.hold")
        self.V = int(comp.get("virtual_nodes", 1))
        self.washout = int(comp.get("washout", 200))
        self.length = int(comp.get("length", 3000))
        self.substeps = int(comp.get("substeps", 20))
        self.solver = comp.get("solver", "rk4")
        self.form = comp.get("form", "flow")
        if self.amplitude <= 0 or self.hold <= 0 or self.V < 1 or self.substeps < 1:
            raise SpecError("amplitude, hold, virtual_nodes and substeps must be positive")
        if self.substeps % self.V:
            raise SpecError("computation.substeps must be a multiple of computation.virtual_nodes")
        if self.solver not in ("rk4", "stiff") or self.form not in ("flow", "map"):
            raise SpecError("computation.solver must be rk4 or stiff, computation.form flow or map")
        if self.form == "map" and self.V != 1:
            raise SpecError("a map is measured once per interval")
        self.expect = comp.get("expect")
        self.exact_box = comp.get("exact_box")

    # ------------------------------------------------------------------------------------------- simulation
    def draw(self, rng, shape) -> np.ndarray:
        if self.law == "uniform":
            return rng.uniform(-self.amplitude, self.amplitude, shape)
        if self.law == "gaussian":
            return self.amplitude * rng.standard_normal(shape)
        return self.amplitude * rng.choice([-1.0, 1.0], size=shape)

    def interval(self, X: np.ndarray, u: np.ndarray, record: bool = False):
        """The state after one interval with the input u held; with record, also the states at the V measurement
        times (fourth-order Runge-Kutta, or the map itself)."""
        if self.form == "map":
            X = self.f(X, u)
            return (X, [X]) if record else X
        h = self.hold / self.substeps
        every = self.substeps // self.V
        nodes = []
        for s in range(1, self.substeps + 1):
            k1 = self.f(X, u)
            k2 = self.f(X + h / 2 * k1, u)
            k3 = self.f(X + h / 2 * k2, u)
            k4 = self.f(X + h * k3, u)
            X = X + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            if record and s % every == 0:
                nodes.append(X.copy())
        return (X, nodes) if record else X

    def drive(self, rng, streams: int = 16, length: int | None = None) -> Tuple[np.ndarray, np.ndarray]:
        """Independent streams side by side. Returns the inputs without their offset (streams, intervals) and the
        measured signals (streams, intervals, m V)."""
        T = self.washout + (self.length if length is None else length)
        U = self.draw(rng, (streams, T))
        X = np.repeat(self.steady()[None, :], streams, axis=0)
        Y = np.empty((streams, T, self.m * self.V))
        if self.solver == "stiff" and self.form == "flow":
            return U, self._drive_stiff(U, X, Y)
        for t in range(T):
            u = U[:, t] + self.offset
            X, nodes = self.interval(X, u, record=True)
            for j, Xj in enumerate(nodes):
                Y[:, t, j * self.m:(j + 1) * self.m] = self.observe(Xj, u)
            if not np.all(np.isfinite(X)):
                raise FloatingPointError(f"the state diverged at interval {t}")
        return U, Y

    def _drive_stiff(self, U: np.ndarray, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
        from scipy.integrate import solve_ivp
        from scipy.sparse import block_diag
        S, n = X.shape
        pattern = block_diag([np.ones((n, n))] * S).tocsc()
        nodes = self.hold * (np.arange(self.V) + 1) / self.V
        for t in range(U.shape[1]):
            u = U[:, t] + self.offset
            sol = solve_ivp(lambda _, y: self.f(y.reshape(S, n), u).ravel(), (0.0, self.hold), X.ravel(),
                            method="BDF", t_eval=nodes, jac_sparsity=pattern, rtol=1e-8, atol=1e-12)
            if not sol.success:
                raise FloatingPointError(f"the stiff integration failed at interval {t}: {sol.message}")
            for j in range(self.V):
                Y[:, t, j * self.m:(j + 1) * self.m] = self.observe(sol.y[:, j].reshape(S, n), u)
            X = sol.y[:, -1].reshape(S, n)
        return Y


class Driven(Body):
    """An equations body (schema fieldbridge-computation/1)."""

    def __init__(self, spec: Dict, path: str = "<dict>"):
        if spec.get("schema") != SCHEMA:
            raise SpecError(f"schema must be {SCHEMA!r}")
        comp = spec.get("computation")
        body = {k: v for k, v in spec.items() if k != "computation"}
        body["schema"] = memory_spec.SCHEMA
        self.real = memory_spec.load(body)
        self.spec, self.path = spec, path
        self._protocol(comp)
        self.input = memory_spec._require(comp, "input", str)
        if self.input not in self.real.params:
            raise SpecError("computation.input must be a declared parameter")
        sym = self.real.symbolic
        self.vsyms, self.psyms = sym["variables"], sym["parameters"]
        names = {str(s): s for s in self.vsyms + self.psyms}
        obs = memory_spec._require(comp, "observables", list)
        if not obs:
            raise SpecError("computation.observables must name at least one expression")
        self.obs_exprs = [memory_spec.parse_expression(e, names) for e in obs]
        self.obs_text = [str(e) for e in obs]
        sp = memory_spec._sympy()
        self._obs = sp.lambdify(self.vsyms + self.psyms, self.obs_exprs, modules="numpy")
        self.n, self.m = len(self.vsyms), len(self.obs_exprs)
        self.variables = [str(v) for v in self.vsyms]
        self.initial = comp.get("initial")

    def params(self, u) -> Dict:
        p = dict(self.real.params)
        p[self.input] = u
        return p

    def f(self, X: np.ndarray, u) -> np.ndarray:
        return self.real.drift(np.atleast_2d(X), self.params(u))

    def observe(self, X: np.ndarray, u) -> np.ndarray:
        X = np.atleast_2d(X)
        p = self.params(u)
        with np.errstate(divide="ignore", invalid="ignore"):
            vals = self._obs(*[X[:, k] for k in range(self.n)], *[p[str(s)] for s in self.psyms])
        return np.stack([np.broadcast_to(np.asarray(v, float), X[:, 0].shape) for v in vals], axis=-1)

    def steady(self, u0: float | None = None) -> np.ndarray:
        """The steady state (flow) or the fixed point (map) at the input u0 (default: the offset)."""
        u0 = self.offset if u0 is None else u0
        x0 = np.asarray(self.initial if self.initial is not None else np.zeros(self.n), float)
        if self.form == "map":
            g = lambda x: self.f(x[None, :], u0)[0] - x  # noqa: E731
        else:
            g = lambda x: self.f(x[None, :], u0)[0]  # noqa: E731
        # judged by the residual: started at the steady state itself, the solver makes no progress and reports failure
        small = 1e-10 * (1.0 + float(np.abs(x0).max(initial=0.0)))
        if np.abs(g(x0)).max(initial=0.0) < small:
            return x0
        sol = root(g, x0, tol=1e-12)
        if not sol.success and np.abs(sol.fun).max(initial=0.0) >= small:
            raise SpecError(f"no steady state found from {x0.tolist()} at the input {u0}; give computation.initial")
        return sol.x


def load(source: Union[str, Path, Dict]) -> Body:
    """A driven body from a specification (a path or a dictionary)."""
    spec, path = _read(source)
    for key in ("question", "assumptions", "provenance"):
        if key not in spec:
            raise SpecError(f"Missing field {key!r}")
    if spec.get("schema") == SPRINGS:
        from .springs import SpringNetwork
        body = SpringNetwork(spec, path)
    else:
        body = Driven(spec, path)
    body.name = str(spec.get("name", spec["question"]))
    body.field = str(spec.get("field", ""))
    prov = spec["provenance"]
    body.source = str(prov.get("source", "")) if isinstance(prov, dict) else str(prov)
    body.spec_sha256 = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()
    return body
