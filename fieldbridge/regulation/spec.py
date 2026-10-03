"""Regulated realizations from specifications (schema ``fieldbridge-regulation/1``).

A specification is an ``equations`` body as in the memory module (carrier, parameters, drift, read by the same
restricted parser: nothing is evaluated with eval and no input is passed to sympify) with a block ``regulation``:

    input      a declared parameter that the protocol steps (slot P)
    output     an expression in the variables and parameters: the observable y = h(q, u) (slot R)
    steps      input values; the first is the reference value u0, every other one is a value after a step
    initial    a state from which the steady state at u0 is reached (default: 1 for every variable)
    fixed      parameters held fixed when parameter points are sampled (the input is always held)
    set_point  the set point stated by the source, if any; it is compared with the result, never used
    expect     the class that the source or a closed form gives (benchmarks and controls)
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import sympy as sp

from ..memory import spec as memory_spec

SCHEMA = "fieldbridge-regulation/1"
SpecError = memory_spec.SpecError
parse_expression = memory_spec.parse_expression


@dataclass
class Regulated:
    """A body with an input that is stepped and an output that is measured."""
    name: str
    field: str
    source: str
    path: str
    spec_sha256: str
    variables: List[str]
    params: Dict[str, float]
    input: str
    output_text: str
    steps: List[float]
    initial: np.ndarray
    fixed: tuple
    set_point: Optional[float]
    orthant: bool
    q: list
    p: list
    F: sp.Matrix
    h: sp.Expr
    spec: Dict = field(repr=False, default_factory=dict)

    def __post_init__(self):
        self.pnames = [str(s) for s in self.p]
        self.n = len(self.q)
        self.u_index = self.pnames.index(self.input)
        u = self.p[self.u_index]
        args = list(self.q) + list(self.p)
        self._f = sp.lambdify(args, list(self.F), modules="numpy")
        self._jac = sp.lambdify(args, self.F.jacobian(self.q).tolist(), modules="numpy")
        self._fu = sp.lambdify(args, list(self.F.diff(u)), modules="numpy")
        self._h = sp.lambdify(args, self.h, modules="numpy")
        self._gradh = sp.lambdify(args, [self.h.diff(s) for s in self.q], modules="numpy")
        self._hu = sp.lambdify(args, self.h.diff(u), modules="numpy")
        # structural interaction graph: edge j -> i when the drift of i contains variable j
        self.depends = [[self.F[i].has(self.q[j]) for j in range(self.n)] for i in range(self.n)]
        self.input_enters = [self.F[i].has(u) for i in range(self.n)]
        self.output_reads = [self.h.has(s) for s in self.q]
        self.output_reads_input = self.h.has(u)

    def pvec(self, u: Optional[float] = None, **overrides: float) -> np.ndarray:
        vals = dict(self.params)
        vals.update(overrides)
        if u is not None:
            vals[self.input] = u
        return np.array([vals[k] for k in self.pnames], float)

    @property
    def u0(self) -> float:
        return float(self.steps[0])

    def f(self, q: np.ndarray, p: np.ndarray) -> np.ndarray:
        return np.array(self._f(*q, *p), float)

    def jac(self, q: np.ndarray, p: np.ndarray) -> np.ndarray:
        return np.array(self._jac(*q, *p), float).reshape(self.n, self.n)

    def fu(self, q: np.ndarray, p: np.ndarray) -> np.ndarray:
        return np.array(self._fu(*q, *p), float).reshape(self.n)

    def y(self, q: np.ndarray, p: np.ndarray) -> float:
        return float(self._h(*q, *p))

    def gradh(self, q: np.ndarray, p: np.ndarray) -> np.ndarray:
        return np.array(self._gradh(*q, *p), float).reshape(self.n)

    def hu(self, q: np.ndarray, p: np.ndarray) -> float:
        return float(self._hu(*q, *p))

    def f_batch(self, Q: np.ndarray, P: np.ndarray) -> np.ndarray:
        """Drift at sampled states Q (m x n) with parameter rows P (m x n_p)."""
        with np.errstate(divide="ignore", invalid="ignore"):
            vals = self._f(*Q.T, *P.T)
        return np.stack([np.broadcast_to(np.asarray(v, float), (len(Q),)) for v in vals], axis=1)

    def derivatives_batch(self, Q: np.ndarray, P: np.ndarray):
        """d F/dq (m x n x n), d F/du (m x n), grad h (m x n) and dh/du (m) at sampled states, in vectorized calls."""
        cnt = len(Q)
        with np.errstate(all="ignore"):
            J, Fu, gh, hu = (self._jac(*Q.T, *P.T), self._fu(*Q.T, *P.T), self._gradh(*Q.T, *P.T),
                             self._hu(*Q.T, *P.T))
        b = lambda v: np.broadcast_to(np.asarray(v, float), (cnt,))
        Jm = np.stack([np.stack([b(J[i][j]) for j in range(self.n)], axis=-1) for i in range(self.n)], axis=1)
        return (Jm, np.stack([b(v) for v in Fu], axis=-1), np.stack([b(v) for v in gh], axis=-1), b(hu))

    def y_batch(self, Q: np.ndarray, P: np.ndarray) -> np.ndarray:
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.broadcast_to(np.asarray(self._h(*Q.T, *P.T), float), (len(Q),)).copy()


def input_hash(spec: Dict) -> str:
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()


def load(source: Union[str, Path, Dict]) -> Regulated:
    """Read a specification of schema fieldbridge-regulation/1 (a path or an already parsed dictionary)."""
    path = "" if isinstance(source, dict) else str(source)
    spec = source if isinstance(source, dict) else json.loads(Path(source).read_text(encoding="utf-8"))
    if spec.get("schema") != SCHEMA:
        raise SpecError(f"schema must be {SCHEMA!r}")
    reg = spec.get("regulation")
    if not isinstance(reg, dict):
        raise SpecError("a regulated specification needs a block 'regulation' with input, output and steps")
    if spec.get("kind") != "equations":
        raise SpecError("the regulation module reads specifications of kind 'equations'")
    body = {k: v for k, v in spec.items() if k != "regulation"}
    body["schema"] = memory_spec.SCHEMA          # the body is read by the memory parser
    real = memory_spec.load(body)
    sym = real.symbolic
    q, p = list(sym["variables"]), list(sym["parameters"])
    names = {str(s): s for s in q + p}
    u = reg.get("input")
    if u not in real.params:
        raise SpecError("regulation.input must be a declared parameter")
    out = reg.get("output")
    if not isinstance(out, str):
        raise SpecError("regulation.output must be an expression in the variables and parameters")
    h = parse_expression(out, names)
    steps = [float(s) for s in reg.get("steps", [])]
    if len(steps) < 2 or not all(np.isfinite(steps)):
        raise SpecError("regulation.steps must list the reference input and at least one stepped value")
    variables = list(real.variables)
    init = reg.get("initial", {})
    if set(init) - set(variables):
        raise SpecError("regulation.initial names an undeclared variable")
    initial = np.array([float(init.get(v, 1.0)) for v in variables])
    fixed = tuple(reg.get("fixed", ()))
    if set(fixed) - set(real.params):
        raise SpecError("regulation.fixed names an undeclared parameter")
    sp_val = reg.get("set_point")
    prov = spec.get("provenance", {})
    return Regulated(name=str(spec.get("name", Path(path).stem if path else "")), field=str(spec.get("field", "")),
                     source=str(prov.get("source", "")) if isinstance(prov, dict) else str(prov), path=path,
                     spec_sha256=input_hash(spec), variables=variables, params=dict(real.params), input=u,
                     output_text=out, steps=steps, initial=initial, fixed=tuple(sorted(set(fixed) | {u})),
                     set_point=None if sp_val is None else float(sp_val),
                     orthant=real.carrier.kind == "orthant", q=q, p=p, F=sp.Matrix(sym["drift"]), h=h, spec=spec)
