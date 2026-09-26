"""Realizations from declarative specifications (schema ``fieldbridge-memory/1``).

A specification describes a material as a realization I_real = ((Omega, Xi); C, R, P; A) in one JSON file.
Two kinds are read here:

  equations  a drift expression for every variable of the carrier, in declared variables and parameters,
             optionally a potential (gradient systems)
  network    units of one type (rotor, spin, gene) joined by edges, each with a transport

A third kind, field (a lattice field, conserved or not), is read by fields.from_spec and the command
``fieldbridge memory field``.

Expressions are parsed by a restricted arithmetic parser: numbers, declared names, + - * / **, and the functions
sqrt, exp, log, sin, cos, tanh; the constant pi. Nothing is evaluated with eval, and input is never passed to
sympify. The parsed expressions are compiled to vectorized numpy functions.
"""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Dict, List, Union

import numpy as np

from .carriers import Carrier
from .identity import Realization

SCHEMA = "fieldbridge-memory/1"
FUNCTIONS = ("sqrt", "exp", "log", "sin", "cos", "tanh", "abs")
RESERVED = {"t", "pi", "I", "E", "S", "N", "O", "Q"}


class SpecError(ValueError):
    """The specification is outside the supported contract; the message says what to change."""


def _sympy():
    try:
        import sympy as sp
    except ImportError as error:  # pragma: no cover - dependency message
        raise SpecError("The memory constructor needs the memory extra: pip install -e '.[memory]'") from error
    return sp


def parse_expression(text: str, names: Dict[str, object]):
    """Restricted arithmetic: numbers, declared names, + - * / **, unary minus, FUNCTIONS, pi."""
    sp = _sympy()
    text = str(text)
    if len(text) > 4000:
        raise SpecError("Expression is too long")
    try:
        tree = ast.parse(text, mode="eval")
    except SyntaxError as error:
        raise SpecError(f"Expected an arithmetic expression, not {text!r}") from error
    if sum(1 for _ in ast.walk(tree)) > 400:
        raise SpecError("Expression is too complex")
    funcs = {name: (sp.Abs if name == "abs" else getattr(sp, name)) for name in FUNCTIONS}

    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            return sp.Rational(str(node.value)) if isinstance(node.value, int) else sp.Float(node.value)
        if isinstance(node, ast.Name):
            if node.id == "pi":
                return sp.pi
            if node.id in names:
                return names[node.id]
            raise SpecError(f"Undeclared symbol {node.id!r}: declare it as a variable or a parameter")
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            item = visit(node.operand)
            return -item if isinstance(node.op, ast.USub) else item
        if isinstance(node, ast.BinOp):
            left, right = visit(node.left), visit(node.right)
            ops = {ast.Add: lambda a, b: a + b, ast.Sub: lambda a, b: a - b, ast.Mult: lambda a, b: a * b,
                   ast.Div: lambda a, b: a / b, ast.Pow: lambda a, b: a ** b}
            for op, fn in ops.items():
                if isinstance(node.op, op):
                    return fn(left, right)
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in funcs
                and len(node.args) == 1 and not node.keywords):
            return funcs[node.func.id](visit(node.args[0]))
        raise SpecError(f"Unsupported syntax in {text!r}; allowed: numbers, declared names, + - * / **, "
                        + ", ".join(FUNCTIONS) + ", pi")

    return visit(tree.body)


def _require(spec: Dict, key: str, kind=None):
    if key not in spec:
        raise SpecError(f"Missing field {key!r}")
    if kind is not None and not isinstance(spec[key], kind):
        raise SpecError(f"Field {key!r} must be {kind.__name__ if isinstance(kind, type) else kind}")
    return spec[key]


def _number(value, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SpecError(f"{what} must be a number")
    if not np.isfinite(value):
        raise SpecError(f"{what} must be finite")
    return float(value)


def input_hash(spec: Dict) -> str:
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()


def load(source: Union[str, Path, Dict]) -> Realization:
    spec = json.loads(Path(source).read_text(encoding="utf-8")) if not isinstance(source, dict) else source
    if spec.get("schema") != SCHEMA:
        raise SpecError(f"schema must be {SCHEMA!r}")
    for key in ("question", "assumptions", "provenance"):
        _require(spec, key)
    kind = _require(spec, "kind", str)
    real = _equations(spec) if kind == "equations" else _network(spec) if kind == "network" else None
    if real is None:
        raise SpecError("kind must be 'equations' or 'network' (a 'field' specification is read by "
                        "'fieldbridge memory field')")
    real.spec = spec
    real.spec_sha256 = input_hash(spec)
    return real


# ------------------------------------------------------------------------------------------------ equations
def _carrier(spec: Dict, dim: int) -> Carrier:
    c = _require(spec, "carrier", dict)
    kind = c.get("kind")
    if kind not in ("euclid", "torus", "orthant"):
        raise SpecError("carrier.kind must be euclid (R^n), torus (angles) or orthant (concentrations)")
    period = 2 * np.pi
    if kind == "torus":
        sp = _sympy()
        period = float(parse_expression(c.get("period", "2*pi"), {}).evalf())
        if not period > 0:
            raise SpecError("carrier.period must be positive")
    scale = _number(c.get("scale", 1.0), "carrier.scale")
    return Carrier(str(c.get("name", "state")), dim, kind, period=period, scale=scale)


def _common(spec: Dict, params: Dict[str, float]) -> Dict:
    control = spec.get("control")
    out = {"noise": _number(_require(spec, "noise"), "noise"), "dt": _number(spec.get("dt", 0.01), "dt"),
           "closure": str(spec.get("closure", "")), "observable": str(spec.get("observable", "")),
           "provenance": str(spec["provenance"].get("source", "")) if isinstance(spec["provenance"], dict)
           else str(spec["provenance"]), "control": None, "control_range": None,
           "material": tuple(spec.get("material", ()))}
    if control is not None:
        name = _require(control, "name", str)
        if name not in params:
            raise SpecError(f"control {name!r} must be a declared parameter")
        lo, hi = (_number(v, "control.range") for v in _require(control, "range", list))
        if not lo < hi:
            raise SpecError("control.range must be increasing")
        out["control"], out["control_range"] = name, (lo, hi)
    for m in out["material"]:
        if m not in params:
            raise SpecError(f"material parameter {m!r} must be declared")
    return out


def _equations(spec: Dict) -> Realization:
    sp = _sympy()
    variables = _require(_require(spec, "carrier", dict), "variables", list)
    if not variables or len(set(variables)) != len(variables):
        raise SpecError("carrier.variables must be distinct names")
    params = {k: _number(v, f"parameter {k}") for k, v in _require(spec, "parameters", dict).items()}
    for name in list(variables) + list(params):
        if not str(name).isidentifier() or name in RESERVED or name in FUNCTIONS:
            raise SpecError(f"Invalid or reserved name {name!r}")
    if set(variables) & set(params):
        raise SpecError("A name cannot be both a variable and a parameter")
    vsyms = [sp.Symbol(v, real=True) for v in variables]
    psyms = [sp.Symbol(p, real=True) for p in params]
    names = {**dict(zip(variables, vsyms)), **dict(zip(params, psyms))}
    drift = _require(spec, "drift", dict)
    if set(drift) != set(variables):
        raise SpecError("drift must give one expression for every variable, and no other")
    exprs = [parse_expression(drift[v], names) for v in variables]
    fn = sp.lambdify(vsyms + psyms, exprs, modules="numpy")
    pot_fn = None
    if "potential" in spec:
        pot = parse_expression(spec["potential"], names)
        pot_fn = sp.lambdify(vsyms + psyms, pot, modules="numpy")
        grads = [sp.simplify(-sp.diff(pot, v) - e) for v, e in zip(vsyms, exprs)]
        if any(g != 0 for g in grads):
            raise SpecError("potential must satisfy drift = -grad(potential)")
    pnames = list(params)

    def drift_fn(q, p):
        cols = [q[:, k] for k in range(q.shape[1])]
        with np.errstate(divide="ignore", invalid="ignore"):
            vals = fn(*cols, *[p[k] for k in pnames])
        return np.stack([np.broadcast_to(np.asarray(v, float), cols[0].shape) for v in vals], axis=-1)

    def potential_fn(q, p):
        cols = [q[:, k] for k in range(q.shape[1])]
        return np.broadcast_to(np.asarray(pot_fn(*cols, *[p[k] for k in pnames]), float), cols[0].shape)

    common = _common(spec, params)
    real = Realization(name=str(spec.get("name", spec["question"])), carrier=_carrier(spec, len(variables)),
                       drift=drift_fn, params=params, potential=potential_fn if pot_fn else None,
                       control=common["control"], control_range=common["control_range"], noise=common["noise"],
                       dt=common["dt"], closure=common["closure"], observable=common["observable"],
                       provenance=common["provenance"], material=common["material"], tags=tuple(spec.get("tags", ())))
    real.variables = list(variables)
    real.symbolic = {"variables": vsyms, "parameters": psyms, "drift": exprs}
    return real


# ------------------------------------------------------------------------------------------------ networks
def _network(spec: Dict) -> Realization:
    from . import networks
    net = _require(spec, "network", dict)
    unit = _require(net, "unit", str)
    noise = _number(_require(spec, "noise"), "noise")
    if unit == "rotor":
        pos = np.asarray(_require(net, "positions", list), float)
        edges = []
        for i, j, k in _require(net, "edges", list):
            terms = {k: 1.0} if isinstance(k, str) else k
            if not isinstance(terms, dict) or not terms or any(t not in ("align", "anti", "reflect") for t in terms):
                raise SpecError("a rotor edge is [i, j, kind] or [i, j, {kind: strength}] with kind align, anti or reflect")
            edges.append((int(i), int(j), k if isinstance(k, str) else {t: _number(v, "edge strength")
                                                                         for t, v in terms.items()}))
        real = networks.rotor_network(pos, edges, m=int(net.get("m", 2)), J=_number(net.get("J", 1.0), "J"),
                                      noise=noise)
    elif unit == "spin":
        n = int(_require(net, "n"))
        edges = [(int(i), int(j), int(s)) for i, j, s in _require(net, "edges", list)]
        real = networks.spin_network(n, edges, k=_number(net.get("k", 0.4), "k"), eps=_number(net.get("eps", 1.0), "eps"),
                                     noise=noise)
    elif unit == "gene":
        n = int(_require(net, "n"))
        regs = [(int(i), int(j), int(s)) for i, j, s in _require(net, "regulations", list)]
        real = networks.gene_network(n, regs, alpha=_number(net.get("alpha", 10.0), "alpha"),
                                     hill=_number(net.get("hill", 4.0), "hill"), noise=noise)
    else:
        raise SpecError("network.unit must be rotor, spin or gene")
    for key in ("lam",):
        if key in net:
            real.params[key] = _number(net[key], f"network.{key}")
    control = spec.get("control")
    if control is not None:
        name = _require(control, "name", str)
        if name not in real.params:
            raise SpecError(f"control {name!r} must be a parameter of the network ({sorted(real.params)})")
        lo, hi = (_number(v, "control.range") for v in _require(control, "range", list))
        real.control, real.control_range = name, (lo, hi)
    real.name = str(spec.get("name", real.name))
    real.observable = str(spec.get("observable", real.observable))
    return real
