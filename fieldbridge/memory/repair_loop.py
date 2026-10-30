"""From a recorded law departure to a tested description of smaller dimension.

`learning.learn` tests the exact identity Dh F = G(h). A reduction at a write point closes only order by order, so
that identity fails even where the reduction is accurate. This module adds three things beside it.

The task generator writes a `fieldbridge-learning/1` task for a law of a record: the coordinate u = w.(x - x*) along
the left critical mode, the linear-stage drift the failed law assumes, and development and withheld interventions on a
ladder of radii. It uses the specification, the write point and the sweep only; the interventions are fixed before
any candidate is derived.

The closure modulo order N expands h(x, p) and g(z, p) about the write point (x*, p*) and solves

    Dh F = g(h)   modulo order N + 1 in (x - x*, p - p*),

with p constant in time. The level sets of h are the stable fibres over the centre manifold, so h also applies to
preparations away from that manifold. Monomials of h free of the non-critical modes carry a zero eigenvalue and go to
g; the others are solved order by order, and a vanishing denominator is reported as a resonance. `verify_modulo_order`
recomputes the residual from the parsed drift along rays through the point and does not use the series arithmetic of
the solver. One order and one radius are selected on development interventions, frozen, and tested once on withheld
interventions through `MoveLibrary.reserve_withheld`, the registry of the exact search.

The size rule reports a correction whose coordinates determine the source state as restating the model, whatever its
trajectories do, and the same for one of the source dimension or larger.

Scope: deterministic drift at one critical direction. The specification's noise is set to zero, so a tested reduction
is a statement about the drift near the write point; it does not re-derive the selection probabilities of the failed
law. The series backend accepts Euclidean, torus (in the covering space) and orthant (at interior points) carriers, any
source dimension within the monomial budget and smooth drifts. No result here establishes novelty.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
from itertools import combinations_with_replacement
import hashlib
import json
from math import comb, sqrt
from pathlib import Path
import time
from typing import Dict, List, Optional, Tuple

import mpmath as mp
import numpy as np
import sympy as sp

from . import analysis as an
from . import construct
from . import learning as ln
from . import spec as mspec
from . import completion as cm

ORDER_ENGINE = "polynomial-intertwining-modulo-order/1"
FAILED_LAW = "unexplained_law_departure_requires_convergence_checks"


@dataclass(frozen=True)
class OrderPolicy:
    """The rule of a campaign, fixed before any withheld response is read."""
    max_order: int = 7
    orders: Tuple[int, ...] = ()              # empty: 1 .. max_order
    radii: Tuple[float, ...] = (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)
    per_radius: int = 3
    prediction_tolerance: float = 0.01        # of the amplitude of the coordinate along the trajectory
    development_margin: float = 0.5           # development must pass at this fraction of the tolerance
    refinement_tolerance: float = 1e-5
    minimum_omission: float = 0.1             # the failed-law drift must be at least this wrong
    minimum_relative_improvement: float = 0.9
    identity_tolerance: float = 1e-6
    power_factor: float = 1.25                # the strongest nonlinear term of g is scaled for the power control
    max_monomials: int = 4000
    rate: float = 0.02                        # the sweep rate of `completion.complete`
    denominator_floor: float = 1e-9

    def __post_init__(self):
        if not (1 <= self.max_order <= 11 and self.radii and 1 <= self.per_radius
                and len(self.radii) * self.per_radius <= 24 and 0 < self.prediction_tolerance < 1):
            raise ValueError("Invalid order-closure policy")

    def order_list(self) -> Tuple[int, ...]:
        return tuple(self.orders) or tuple(range(1, self.max_order + 1))


# ------------------------------------------------------------------------------------------- truncated series
Series = Dict[Tuple[int, ...], float]


def _degree(key) -> int:
    return sum(key)


def _add(a: Series, b: Series, scale: float = 1.0) -> Series:
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0.0) + scale * v
    return out


def _mul(a: Series, b: Series, order: int) -> Series:
    out: Series = {}
    by_degree: Dict[int, List] = {}
    for k, v in b.items():
        by_degree.setdefault(_degree(k), []).append((k, v))
    for ka, va in a.items():
        da = _degree(ka)
        for db, items in by_degree.items():
            if da + db > order:
                continue
            for kb, vb in items:
                key = tuple(i + j for i, j in zip(ka, kb))
                out[key] = out.get(key, 0.0) + va * vb
    return out


def _constant(value: float, d: int) -> Series:
    return {(0,) * d: float(value)}


def _taylor(function, a: Series, order: int, d: int) -> Series:
    """f(a) by the Taylor series of f about the constant term of a."""
    zero = (0,) * d
    a0 = a.get(zero, 0.0)
    delta = {k: v for k, v in a.items() if k != zero}
    t = sp.Symbol("t", real=True)
    expression = function(t)
    out = _constant(float(expression.subs(t, a0)), d)
    power = _constant(1.0, d)
    factorial = 1.0
    derivative = expression
    for k in range(1, order + 1):
        derivative = sp.diff(derivative, t)
        factorial *= k
        power = _mul(power, delta, order)
        out = _add(out, power, float(derivative.subs(t, a0)) / factorial)
    return out


def _power(a: Series, n: int, order: int, d: int) -> Series:
    out = _constant(1.0, d)
    for _ in range(n):
        out = _mul(out, a, order)
    return out


def _series_of(expr, values: Dict[sp.Symbol, Series], order: int, d: int) -> Series:
    """The truncated series of a parsed drift expression, with `values` substituted for its symbols."""
    if expr.is_Symbol:
        return values[expr]
    if expr.is_number:
        return _constant(float(expr), d)
    if expr.is_Add:
        out: Series = {}
        for arg in expr.args:
            out = _add(out, _series_of(arg, values, order, d))
        return out
    if expr.is_Mul:
        out = _constant(1.0, d)
        for arg in expr.args:
            out = _mul(out, _series_of(arg, values, order, d), order)
        return out
    if expr.is_Pow:
        base, exponent = expr.args
        if exponent.is_Integer and exponent >= 0:
            return _power(_series_of(base, values, order, d), int(exponent), order, d)
        if not exponent.is_number:
            raise ValueError("A symbolic exponent is outside the series backend")
        value = float(exponent)
        return _taylor(lambda t: t ** sp.Float(value) if value != int(value) else t ** int(value),
                       _series_of(base, values, order, d), order, d)
    if isinstance(expr, sp.Function) and len(expr.args) == 1:
        return _taylor(lambda t: expr.func(t), _series_of(expr.args[0], values, order, d), order, d)
    raise ValueError(f"The series backend does not treat {expr.func}")


# ------------------------------------------------------------------------------------------------- the point
def _carrier_limit(real, state) -> Optional[str]:
    kind = real.carrier.kind
    if kind == "orthant" and np.min(state) <= 0:
        return "the write point lies on the boundary of the orthant, where the chart does not reach both sides"
    return None


def _polish(real, state, over) -> np.ndarray:
    """Newton at fixed parameters, keeping the recorded state where the iteration does not improve it."""
    q = np.asarray(state, float)
    best, best_norm = q, float(np.linalg.norm(real.F(q, **over)))
    for _ in range(12):
        f = real.F(q, **over)
        J = an.jacobian(real, q, over)
        try:
            q = q - np.linalg.lstsq(J, f, rcond=1e-10)[0]
        except np.linalg.LinAlgError:
            break
        norm = float(np.linalg.norm(real.F(q, **over)))
        if norm < best_norm:
            best, best_norm = q, norm
    return best


def _critical_frame(J: np.ndarray):
    """Right and left critical vectors (w.v = 1), a basis S of the kernel of w, and the block on it."""
    ev, R = np.linalg.eig(J)
    i = int(np.argmin(np.abs(ev)))
    if abs(ev[i].imag) > 1e-9 * max(1.0, np.abs(ev).max()):
        raise ValueError("A complex critical pair is outside the one-coordinate closure")
    v = np.real(R[:, i])
    v = v / np.linalg.norm(v)
    evl, L = np.linalg.eig(J.T)
    w = np.real(L[:, int(np.argmin(np.abs(evl - ev[i])))])
    w = w / np.dot(w, v)
    # the complement of v in the kernel of w: the invariant complement of the critical direction
    Q = np.linalg.svd(w.reshape(1, -1))[2][1:].T
    P = np.column_stack([v, Q]) if Q.size else v.reshape(-1, 1)
    Pinv = np.linalg.inv(P)
    block = Pinv @ J @ P
    return float(ev[i].real), v, w, P, Pinv, block[1:, 1:], np.delete(ev, i)


def _monomials(d: int, degree: int) -> List[Tuple[int, ...]]:
    out = []
    for indices in combinations_with_replacement(range(d), degree):
        out.append(tuple(indices.count(i) for i in range(d)))
    return out


# ------------------------------------------------------------------------------------------- the closure
def _solve_closure(real, state, parameters: Dict[str, float], swept: List[str], order: int, policy: OrderPolicy):
    """h and g through `order` at the write point, in the coordinates (u, s, e): u along the critical direction, s the
    other modes after removing their linear response to e, e the swept parameters."""
    over = dict(parameters)
    n = real.carrier.dim
    q = np.asarray(state, float)
    # the symbolic Jacobian: rounding in the critical block would otherwise enter every denominator
    sym = real.symbolic
    x, p = sym["variables"], sym["parameters"]
    point = {**dict(zip(x, q.tolist())), **{s: float({**real.params, **over}[str(s)]) for s in p}}
    J = np.array([[float(sp.diff(f, xi).subs(point)) for xi in x] for f in sym["drift"]])
    lam, v, w, P, Pinv, block, others = _critical_frame(J)
    if others.size and np.min(np.abs(others)) < 1e-6 * max(1.0, np.abs(others).max()):
        raise ValueError("A second eigenvalue is critical at the write point; the closure treats one direction")
    m = len(swept)
    d = n + m
    if sum(comb(d + k - 1, k) for k in range(order + 1)) > policy.max_monomials:
        raise ValueError(f"Order {order} in {d} variables exceeds the budget of {policy.max_monomials} monomials")
    psym = {str(s): s for s in p}
    # linear response of the source to the swept parameters, in the eigenbasis
    B = np.array([[float(sp.diff(f, psym[name]).subs(point)) for name in swept] for f in sym["drift"]]).reshape(n, m)
    Bm = Pinv @ B
    K = -np.linalg.solve(block, Bm[1:, :]) if n > 1 else np.zeros((0, m))
    # substitute x = x* + P (u, s + K e), p = p* + e
    unit = lambda i: {tuple(int(j == i) for j in range(d)): 1.0}
    modal: List[Series] = [unit(0)]
    for i in range(n - 1):
        term = unit(i + 1)
        for j in range(m):
            if K[i, j] != 0.0:
                term = _add(term, unit(n + j), float(K[i, j]))
        modal.append(term)
    values: Dict[sp.Symbol, Series] = {}
    for i, xi in enumerate(x):
        series = _constant(q[i], d)
        for j in range(n):
            if P[i, j] != 0.0:
                series = _add(series, modal[j], float(P[i, j]))
        values[xi] = series
    for s in p:
        name = str(s)
        series = _constant(float({**real.params, **over}[name]), d)
        if name in swept:
            series = _add(series, unit(n + swept.index(name)))
        values[s] = series
    raw = [_series_of(sp.expand(f) if f.is_polynomial(*x) else f, values, order, d) for f in sym["drift"]]
    field_modal: List[Series] = []
    for i in range(n):
        out: Series = {}
        for j in range(n):
            if Pinv[i, j] != 0.0:
                out = _add(out, raw[j], float(Pinv[i, j]))
        field_modal.append(out)

    # the linear field A on (u, s, e): u' = lam u + c.e, s' = block s, e' = 0
    c = [float(Bm[0, j]) for j in range(m)]
    h: Series = {tuple(int(j == 0) for j in range(d)): 1.0}
    g: Series = {}                         # keyed by (power of z, powers of e)
    denominators: List[float] = []
    resonance = None
    for k in range(0, order + 1):
        residual = _order_residual(field_modal, h, g, k, n, m, order, d)
        keys = _monomials(d, k)
        pure, mixed = [], []
        for key in keys:
            (pure if not any(key[1:n]) else mixed).append(key)
        for key in pure:
            value = residual.get(key, 0.0)
            if value != 0.0:
                g[(key[0],) + key[n:]] = g.get((key[0],) + key[n:], 0.0) + value
        if mixed:
            index = {key: i for i, key in enumerate(mixed)}
            L = np.zeros((len(mixed), len(mixed)))
            for key, col in index.items():
                L[col, col] -= lam
                # D(monomial).A, the part along s
                for i in range(n - 1):
                    bi = key[i + 1]
                    if not bi:
                        continue
                    for j in range(n - 1):
                        target = list(key)
                        target[i + 1] -= 1
                        target[j + 1] += 1
                        L[index[tuple(target)], col] += bi * block[i, j]
                # the part along u
                if key[0]:
                    L[col, col] += key[0] * lam
                    for j in range(m):
                        if c[j] != 0.0:
                            target = list(key)
                            target[0] -= 1
                            target[n + j] += 1
                            L[index[tuple(target)], col] += key[0] * c[j]
            # L is triangular in the eigenbasis, so its eigenvalues are the denominators. Its singular values are not:
            # a large bias makes it strongly non-normal, with tiny singular values and no resonance
            smallest = _smallest_denominator(others, lam, k)
            denominators.append(smallest)
            if smallest < policy.denominator_floor * max(1.0, float(np.abs(others).max())):
                if k < 2:
                    raise ValueError(f"resonance at order {k}: a denominator of the homological operator vanishes")
                resonance = k
                break
            rhs = np.array([-residual.get(key, 0.0) for key in mixed])
            solution = np.linalg.solve(L, rhs)
            for key, value in zip(mixed, solution):
                if value != 0.0:
                    h[key] = h.get(key, 0.0) + float(value)
    return {"state": q.tolist(), "parameters": over, "swept": list(swept), "lambda_critical": lam,
            "mode": v.tolist(), "left_mode": w.tolist(), "P_inverse": Pinv.tolist(),
            "K": K.tolist(), "h": h, "g": g, "source_dimension": n, "resonance_at_order": resonance,
            "order": order if resonance is None else resonance - 1,
            "other_spectrum": [[float(z.real), float(z.imag)] for z in others],
            "smallest_denominator": min(denominators) if denominators else None}


def _smallest_denominator(others: np.ndarray, lam: float, k: int) -> float:
    """min |sum beta_i lambda_i + (a - 1) lam| over the monomials of degree k with some non-critical power."""
    smallest = float("inf")
    for b in range(1, k + 1):
        for chosen in combinations_with_replacement(others.tolist(), b):
            total = sum(chosen)
            for a in range(0, k - b + 1):
                smallest = min(smallest, abs(total + (a - 1) * lam))
    return smallest


def _derivative_product(h: Series, fields: List[Series], n: int, order: int) -> Series:
    """Dh . F in the modal coordinates (u, s); the parameters do not evolve."""
    out: Series = {}
    for i in range(n):
        partial: Series = {}
        for key, value in h.items():
            if key[i]:
                lowered = list(key)
                lowered[i] -= 1
                partial[tuple(lowered)] = key[i] * value
        out = _add(out, _mul(partial, fields[i], order))
    return out


def _compose_g(g: Series, h: Series, m: int, n: int, order: int, d: int) -> Series:
    out: Series = {}
    powers = {0: _constant(1.0, d)}
    for a in range(1, max((key[0] for key in g), default=0) + 1):
        powers[a] = _mul(powers[a - 1], h, order)
    for key, value in g.items():
        a, gamma = key[0], key[1:]
        factor = {tuple([0] * n) + tuple(gamma): value}
        out = _add(out, _mul(powers[a], factor, order))
    return out


def _order_residual(fields, h, g, k, n, m, order, d) -> Series:
    """The coefficient at total degree k of Dh F - g(h), with everything known below order k."""
    field_k = [dict(f) for f in fields]
    left = _derivative_product(h, field_k, n, k)
    right = _compose_g(g, h, m, n, k, d)
    out = _add(left, right, -1.0)
    return {key: value for key, value in out.items() if _degree(key) == k}


# ------------------------------------------------------------------------------------------ evaluating a candidate
def truncate(closure: Dict, order: int) -> Dict:
    """The candidate of a lower order from the closure already solved (the same homological equations)."""
    return {**closure, "order": order, "h": {k: v for k, v in closure["h"].items() if sum(k) <= order},
            "g": {k: v for k, v in closure["g"].items() if sum(k) <= order}}


def _arrays(terms: Series):
    keys = list(terms)
    if not keys:
        return np.zeros((0, 1), int), np.zeros(0)
    return np.array(keys, int), np.array([terms[k] for k in keys])


class Evaluator:
    """h(x, p) and g(z, p) from their frozen coefficients; numpy arithmetic, real or complex."""

    def __init__(self, candidate: Dict):
        self.c = candidate
        self.q = np.asarray(candidate["state"], float)
        self.Pinv = np.asarray(candidate["P_inverse"], float)
        self.K = np.asarray(candidate["K"], float)
        self.swept = list(candidate["swept"])
        self.p0 = {k: float(v) for k, v in candidate["parameters"].items()}
        self.n = self.q.size
        self.hk, self.hc = _arrays(candidate["h"])
        self.gk, self.gc = _arrays(candidate["g"])

    def _e(self, params) -> np.ndarray:
        return np.array([params[k] - self.p0[k] for k in self.swept])

    def modal(self, x, params) -> np.ndarray:
        e = self._e(params)
        xi = self.Pinv @ (np.asarray(x) - self.q)
        if self.n > 1 and e.size:
            xi = xi.copy()
            xi[1:] = xi[1:] - self.K @ e
        return np.concatenate([xi, e])

    def h(self, x, params):
        v = self.modal(x, params)
        return np.dot(self.hc, np.prod(v[None, :] ** self.hk, axis=1))

    def gradient_h(self, x, params):
        """Dh in the source coordinates: the chain rule through the frozen linear map."""
        v = self.modal(x, params)
        out = np.zeros(self.n, dtype=v.dtype)
        for i in range(self.n):
            mask = self.hk[:, i] > 0
            lowered = self.hk[mask].copy()
            lowered[:, i] -= 1
            out[i] = np.dot(self.hc[mask] * self.hk[mask, i], np.prod(v[None, :] ** lowered, axis=1))
        return out @ self.Pinv

    def g_coefficients(self, params) -> np.ndarray:
        """The coefficients of g as a polynomial in z at these parameters, lowest first."""
        e = self._e(params)
        top = int(self.gk[:, 0].max()) if self.gk.size else 0
        out = np.zeros(top + 1)
        if self.gk.size:
            weights = np.prod(e[None, :] ** self.gk[:, 1:], axis=1) * self.gc if e.size else self.gc
            np.add.at(out, self.gk[:, 0], weights)
        return out

    def g(self, z, params, *, coefficients=None):
        a = self.g_coefficients(params) if coefficients is None else coefficients
        return np.polyval(a[::-1], z)


def verify_modulo_order(task: Dict, candidate: Dict, policy: OrderPolicy, *, rays: int = 4, seed: int = 20261030):
    """Recompute Dh F and g(h) from the parsed drift along rays through the write point and read the Taylor coefficients
    of each by a discrete Cauchy integral at 60 digits. The series arithmetic of the solver takes no part. Each order
    through N must cancel against the size of its own two sides, or against the rounding of the lower orders."""
    real = mspec.load(task["model"])
    sym = real.symbolic
    names = sym["variables"] + sym["parameters"]
    drift = sp.lambdify(names, sym["drift"], modules="mpmath")
    pnames = [str(s) for s in sym["parameters"]]
    order = candidate["order"]
    q = candidate["state"]
    Pinv = candidate["P_inverse"]
    K = candidate["K"]
    swept = list(candidate["swept"])
    fixed = {**real.params, **candidate["parameters"]}
    n = len(q)
    d = n + len(swept)
    h_terms = list(candidate["h"].items())
    g_terms = list(candidate["g"].items())
    # g(h) has degree order^2 on a polynomial drift of low degree: below that the coefficients alias
    samples, radius = (order + 1) ** 2 + 8, mp.mpf("0.01")
    rng = np.random.default_rng(seed)
    left_size, right_size, residual = [np.zeros(order + 2) for _ in range(3)]
    with mp.workdps(60):
        roots = [mp.expjpi(2 * mp.mpf(j) / samples) for j in range(samples)]
        for _ in range(rays):
            direction = rng.standard_normal(d)
            direction = direction / np.linalg.norm(direction)
            left, right = [], []
            for root in roots:
                s = radius * root
                x = [mp.mpf(q[i]) + s * mp.mpf(direction[i]) for i in range(n)]
                e = [s * mp.mpf(direction[n + j]) for j in range(len(swept))]
                params = {**fixed, **{name: mp.mpf(fixed[name]) + e[j] for j, name in enumerate(swept)}}
                f = drift(*x, *[params[k] for k in pnames])
                f = [f[i] for i in range(n)] if not isinstance(f, (list, tuple)) else list(f)
                xi = [sum(mp.mpf(Pinv[i][k]) * (x[k] - mp.mpf(q[k])) for k in range(n)) for i in range(n)]
                for i in range(1, n):
                    xi[i] -= sum(mp.mpf(K[i - 1][j]) * e[j] for j in range(len(swept)))
                v = xi + e
                top = order + 2
                powers = [[mp.mpf(1)] for _ in v]
                for i, value in enumerate(v):
                    for _power_index in range(top):
                        powers[i].append(powers[i][-1] * value)
                z = mp.mpf(0)
                modal_gradient = [mp.mpf(0) for _ in range(n)]
                for key, c in h_terms:
                    # the product without the factor of variable i gives the partial derivative
                    product = mp.mpf(c)
                    nonzero = [i for i, k in enumerate(key) if k]
                    for i in nonzero:
                        product *= powers[i][key[i]]
                    z += product
                    for i in nonzero:
                        if i < n:
                            rest = mp.mpf(c) * key[i] * powers[i][key[i] - 1]
                            for other in nonzero:
                                if other != i:
                                    rest *= powers[other][key[other]]
                            modal_gradient[i] += rest
                gradient = [sum(modal_gradient[i] * mp.mpf(Pinv[i][k]) for i in range(n)) for k in range(n)]
                left.append(sum(gradient[k] * f[k] for k in range(n)))
                total = mp.mpf(0)
                for key, c in g_terms:
                    term = mp.mpf(c) * z ** key[0]
                    for j, k in enumerate(key[1:]):
                        if k:
                            term *= e[j] ** k
                    total += term
                right.append(total)
            for k in range(order + 2):
                scale = radius ** k * samples
                a = sum(left[j] * mp.expjpi(-2 * mp.mpf(j * k) / samples) for j in range(samples)) / scale
                b = sum(right[j] * mp.expjpi(-2 * mp.mpf(j * k) / samples) for j in range(samples)) / scale
                left_size[k] = max(left_size[k], float(abs(a)))
                right_size[k] = max(right_size[k], float(abs(b)))
                residual[k] = max(residual[k], float(abs(a - b)))
    size = np.maximum(left_size, right_size)
    rounding = 1e-9 * np.maximum.accumulate(np.maximum(size, 1.0))
    limit = np.maximum(policy.identity_tolerance * size, rounding)
    passed = bool(np.all(residual[:order + 1] <= limit[:order + 1]))
    return {"passed": passed, "order": order, "residual_by_order": residual.tolist(),
            "side_size_by_order": size.tolist(), "limit_by_order": limit[:order + 1].tolist(),
            "first_omitted_order": order + 1, "first_omitted_residual": float(residual[order + 1])}


def candidate_from_record(frozen: Dict) -> Dict:
    return {**frozen, "h": {tuple(k): v for k, v in frozen["h"]}, "g": {tuple(k): v for k, v in frozen["g"]}}


def candidate_record(candidate: Dict) -> Dict:
    """A JSON form of a candidate; hashing it freezes it before the withheld test."""
    return {"engine": ORDER_ENGINE, "order": candidate["order"], "state": candidate["state"],
            "parameters": candidate["parameters"], "swept": candidate["swept"],
            "left_mode": candidate["left_mode"], "P_inverse": candidate["P_inverse"], "K": candidate["K"],
            "h": [[list(k), v] for k, v in sorted(candidate["h"].items())],
            "g": [[list(k), v] for k, v in sorted(candidate["g"].items())]}


def description_size(source_dimension: int, retained_dimension: int, rank: Optional[int]) -> Dict:
    """Is a description smaller than its source? Coordinates determining the source state restate the model."""
    if rank is not None and rank >= source_dimension:
        verdict = "restates_model"
    elif retained_dimension >= source_dimension:
        verdict = "not_smaller"
    else:
        verdict = "smaller"
    return {"source_dimension": source_dimension, "retained_dimension": retained_dimension,
            "jacobian_rank": rank, "verdict": verdict}


def assess_exact_report(task: Dict, report: Dict, *, seed: int = 20261030) -> Dict:
    """The size rule on a report of `learning.learn`: the rank of the retained coordinates at generic source points."""
    candidate = report.get("candidate")
    if not candidate:
        return {"verdict": "no_candidate"}
    real = mspec.load(task["model"])
    x, p = real.symbolic["variables"], real.symbolic["parameters"]
    names = {str(s): s for s in x + p}
    obs = [ln.expression(v, names) for v in candidate["coordinates"].values()]
    J = sp.Matrix([[sp.diff(o, xi) for xi in x] for o in obs])
    fn = sp.lambdify(x + p, J, modules="numpy")
    rng = np.random.default_rng(seed)
    rank = 0
    for _ in range(3):
        at = list(rng.normal(size=len(x))) + [real.params[str(s)] for s in p]
        rank = max(rank, int(np.linalg.matrix_rank(np.asarray(fn(*at), float).reshape(len(obs), len(x)), tol=1e-9)))
    return description_size(len(x), len(obs), rank)


def _trajectory(rhs, initial, duration, *, fine=False, reduced=False):
    return ln._trajectory(rhs, np.atleast_1d(np.asarray(initial, float)), duration, fine=fine, reduced=reduced)


# ------------------------------------------------------------------------------------------- generating a task
@dataclass
class WritePoint:
    real: object
    state: np.ndarray
    parameters: Dict[str, float]
    swept: List[str]
    direction: Dict[str, float]
    end_factor: float
    target: str
    failed_law: Dict
    note: str = ""
    extra: Dict = field(default_factory=dict)


def _locate_laws(record: Dict):
    completion = record.get("completion")
    if not isinstance(completion, dict):
        return
    for i, point in enumerate(completion.get("write_points", [])):
        if isinstance(point.get("law"), dict):
            yield f"completion/{i}/law", point, None, point["law"]
        for j, unfolding in enumerate(point.get("unfoldings", [])):
            if isinstance(unfolding.get("law"), dict):
                yield f"completion/{i}/unfolding/{j}", point, unfolding, unfolding["law"]


def _write_point(record: Dict, spec: Dict, target: str, point: Dict, unfolding: Optional[Dict], law: Dict,
                 policy: OrderPolicy) -> WritePoint:
    real = mspec.load({**spec, "noise": 0})
    if not hasattr(real, "symbolic"):
        raise ValueError(f"The series backend treats equations, not kind {spec.get('kind')!r}")
    control = real.control
    wp = point["write_point"]
    state = np.asarray(wp["state"], float)
    if unfolding is None:
        over = {control: float(wp[control])}
        q = _polish(real, state, over)
        nf = an.normal_form(real, q, control, float(wp[control]))
        g = abs(nf["kappa_slope"])
        if g <= 0:
            raise ValueError("The critical rate does not change along the control")
        sign = float(np.sign(-nf["kappa_slope"]))
        return WritePoint(real, q, over, [control], {control: sign},
                          4.0 * sqrt(policy.rate) / g, target, law)
    second = unfolding["second"]
    nf = an.normal_form(real, state, control, float(wp[control]))
    found = construct.cancel_asymmetry(real, nf, second)
    if not found.get("success"):
        raise ValueError("The cusp of this record is no longer found from its write point")
    cusp = {control: float(found["control"]), second: float(found["value"])}
    recorded = unfolding.get("cusp", {})
    drift_of_cusp = max(abs(cusp[k] - float(recorded[k])) for k in cusp) if recorded else None
    q = _polish(real, np.asarray(found["state"], float), cusp)
    path = construct.along_path(real, cusp, found["direction"])
    nfp = construct.sweep_normal_form(path, q)
    g = abs(nfp["kappa_slope"])
    if g <= 0:
        raise ValueError("The critical rate does not change along the sweep")
    sign = float(np.sign(-nfp["kappa_slope"]))
    direction = {k: sign * float(found["direction"][k]) for k in cusp}
    note = "" if drift_of_cusp is None or drift_of_cusp < 1e-3 else (
        f"the cusp found now differs from the record by {drift_of_cusp:.3g}")
    return WritePoint(real, q, cusp, [control, second], direction, 4.0 * sqrt(policy.rate) / g, target, law, note,
                      {"cusp_change_from_record": drift_of_cusp})


def _written_scale(pt: WritePoint, v: np.ndarray, w: np.ndarray, seed: int) -> Dict:
    """The distance along the critical coordinate between the states the sweep ends in."""
    real = pt.real
    end = {**pt.parameters, **{k: pt.parameters[k] + pt.end_factor * pt.direction[k] for k in pt.swept}}
    rng = np.random.default_rng(seed)
    length = float(np.linalg.norm(v)) * max(1e-9, float(an.base_length(real)))
    seeds = [pt.state + sign * f * length * v for sign in (1, -1) for f in (0.05, 0.2, 0.5)]
    stable = cm.stable_states(real, rng, end, n_starts=16, seeds=np.array(seeds))
    coordinates = sorted(float(np.dot(w, s - pt.state)) for s in stable)
    if len(coordinates) >= 2:
        gaps = [abs(c) for c in coordinates]
        above = [c for c in coordinates if c > 0]
        below = [c for c in coordinates if c < 0]
        # the pair on either side of the write point; without one the two nearest
        if above and below:
            scale = min(above) - max(below)
        else:
            order = np.argsort(gaps)
            scale = abs(coordinates[order[0]] - coordinates[order[1]])
        kind = "separation of the written states"
    elif coordinates:
        scale, kind = abs(coordinates[0]), "distance to the one state at the end of the sweep"
    else:
        scale, kind = 0.0, ""
    if not scale > 1e-9:
        scale, kind = 0.25 * float(an.base_length(real)), "a quarter of the base length (no written state found)"
    return {"scale": float(scale), "kind": kind, "stable_critical_coordinates": coordinates,
            "end_parameters": {k: float(end[k]) for k in pt.swept}}


def _interventions(pt: WritePoint, v, w, Q, length, policy: OrderPolicy, seed: int, tag: str):
    real = pt.real
    rng = np.random.default_rng(seed)
    n, rows = pt.state.size, []
    duration = min(100.0, 8.0 / sqrt(policy.rate))
    swept = pt.swept
    transverse = {swept[0]: -pt.direction[swept[1]], swept[1]: pt.direction[swept[0]]} if len(swept) == 2 else None
    for radius in policy.radii:
        for i in range(policy.per_radius):
            sign = 1.0 if i % 3 != 2 else -1.0                # most past the crossing, some before it
            factor = sign * float(rng.uniform(0.4, 1.0)) * radius * pt.end_factor
            params = {k: float(pt.parameters[k] + factor * pt.direction[k]) for k in swept}
            if transverse:
                off = 0.2 * float(rng.uniform(-1, 1)) * radius * pt.end_factor
                params = {k: params[k] + off * transverse[k] for k in swept}
            for name, (lo, hi) in real.bounds.items():
                if name in params:
                    params[name] = float(min(max(params[name], lo if lo is not None else -np.inf),
                                             hi if hi is not None else np.inf))
            along = float(rng.choice([-1.0, 1.0])) * float(rng.uniform(0.3, 1.0)) * radius * length
            displacement = along * v
            if n > 1:
                displacement = displacement + Q @ rng.normal(size=n - 1) * 0.3 * radius * length / np.sqrt(n - 1)
            rows.append({"id": f"{tag}-r{radius:g}-{i}", "initial": (pt.state + displacement).tolist(),
                         "parameters": params, "duration": duration, "radius": radius})
    return rows


def learning_task(record: Dict, target: str, *, root: Path = Path("."), policy: Optional[OrderPolicy] = None) -> Dict:
    """Write the learning task of one law of a record, from its specification, write point and failed law."""
    policy = policy or OrderPolicy()
    spec = json.loads((Path(root) / record["source"]).read_text(encoding="utf-8"))
    for name, point, unfolding, law in _locate_laws(record):
        if name == target:
            break
    else:
        raise KeyError(f"No law {target!r} in this record")
    pt = _write_point(record, spec, target, point, unfolding, law, policy)
    limit = _carrier_limit(pt.real, pt.state)
    if limit:
        raise ValueError(limit)
    return _task_from_point(record, spec, pt, policy)


def _task_from_point(record: Dict, spec: Dict, pt: WritePoint, policy: OrderPolicy) -> Dict:
    real = pt.real
    sym = real.symbolic
    x, p = sym["variables"], sym["parameters"]
    point = {**dict(zip(x, pt.state.tolist())), **{s: float({**real.params, **pt.parameters}[str(s)]) for s in p}}
    psym = {str(s): s for s in p}
    J = np.array([[float(sp.diff(f, xi).subs(point)) for xi in x] for f in sym["drift"]])
    lam, v, w, P, Pinv, block, others = _critical_frame(J)
    Q = P[:, 1:]
    B = np.array([[float(sp.diff(f, psym[name]).subs(point)) for name in pt.swept] for f in sym["drift"]])
    c = w @ B
    scale_info = _written_scale(pt, v, w, seed=int(ln.digest(pt.target)[:8], 16))
    length = scale_info["scale"]
    used = {str(s) for s in x + p}
    name = next(c for c in ("u", "u_critical", "u_critical_mode") if c not in used)
    u_text = "+".join(f"({float(w[i])!r})*({v_name}-({float(pt.state[i])!r}))" for i, v_name in enumerate(map(str, x)))
    terms = [f"({float(lam)!r})*{name}"]
    for swept, value in zip(pt.swept, c):
        terms.append(f"({float(value)!r})*({swept}-({pt.parameters[swept]!r}))")
    slope = _rate_slope(real, pt)
    for swept, value in zip(pt.swept, slope):
        terms.append(f"({float(value)!r})*{name}*({swept}-({pt.parameters[swept]!r}))")
    key = ln.digest({"record": record.get("spec_sha256"), "target": pt.target})
    seed = int(key[:12], 16)
    model = {**spec, "noise": 0}
    task = {"schema": ln.SCHEMA, "model": model, "coordinates": {name: u_text},
            "candidate_drift": {name: "+".join(terms)},
            "development": _interventions(pt, v, w, Q, length, policy, seed, "d"),
            "withheld": _interventions(pt, v, w, Q, length, policy, seed + 1, "w"),
            "write_point": {"state": pt.state.tolist(), "parameters": pt.parameters, "swept": pt.swept,
                            "direction": pt.direction, "sweep_end_factor": pt.end_factor,
                            "scale": scale_info, "note": pt.note, **pt.extra},
            "failed_law": {"target": pt.target, "write_class": pt.failed_law.get("write_class"),
                           "certified": pt.failed_law.get("certified"),
                           "law_kind": pt.failed_law.get("law_kind"),
                           "obstruction": pt.failed_law.get("obstruction_in_nonlinear_stage"),
                           "one_direction": pt.failed_law.get("one_direction")},
            "record": {"name": record.get("name"), "source": record.get("source"),
                       "spec_sha256": record.get("spec_sha256"),
                       "implementation_sha256": record.get("implementation_sha256")},
            "noise_set_to_zero": spec.get("noise")}
    ln._validate_cases(task, real)
    return task


def _rate_slope(real, pt) -> List[float]:
    """d(lambda)/d(e_j) by a central difference of the critical eigenvalue along the branch."""
    out = []
    for name, base in pt.parameters.items():
        if name not in pt.swept:
            continue
        width = 1e-5 * max(1.0, abs(base))
        values = []
        for sign in (1, -1):
            over = {**pt.parameters, name: base + sign * width}
            q = _polish(real, pt.state, over)
            ev = np.linalg.eigvals(an.jacobian(real, q, over))
            values.append(float(ev[np.argmin(np.abs(ev))].real))
        out.append((values[0] - values[1]) / (2 * width))
    return out


# ----------------------------------------------------------------------------------------------- one campaign
def _row(ev: Evaluator, full, fine, rhs, case, policy: OrderPolicy):
    """Prediction error of a reduced drift `rhs` on one intervention, against the same observable of the source."""
    params = case["parameters"]
    duration = case["duration"]
    truth = np.array([ev.h(q, params) for q in fine])
    swing = max(float(np.max(np.abs(truth - truth.mean()))), float(np.max(np.abs(truth))))
    z0 = truth[0]
    try:
        reduced = _trajectory(rhs, [z0], duration, reduced=True)[:, 0]
        refined = _trajectory(rhs, [z0], duration, fine=True, reduced=True)[:, 0]
    except (ArithmeticError, ValueError) as error:
        return {"id": case["id"], "radius": case.get("radius"), "error": float("inf"), "blown_up": str(error)[:120],
                "refinement": float("inf"), "swing": swing, "passed": False}
    # the coarse and fine source trajectories bound the solver error of the comparison
    coarse = np.array([ev.h(q, params) for q in full])
    scale = max(swing, 1e-12)
    error = float(np.max(np.abs(refined - truth))) / scale
    refinement = float(max(np.max(np.abs(coarse - truth)), np.max(np.abs(reduced - refined)))) / scale
    return {"id": case["id"], "radius": case.get("radius"), "error": error, "refinement": refinement,
            "swing": swing, "passed": error <= policy.prediction_tolerance
            and refinement <= policy.refinement_tolerance}


def _source_trajectories(real, cases):
    out = {}
    for case in cases:
        params = {**real.params, **case["parameters"]}
        rhs = lambda t, q, params=params: real.F(q, **params)
        try:
            coarse = _trajectory(rhs, case["initial"], case["duration"])
            fine = _trajectory(rhs, case["initial"], case["duration"], fine=True)
            out[case["id"]] = (coarse, fine, None)
        except (ArithmeticError, ValueError) as error:
            out[case["id"]] = (None, None, str(error)[:120])
    return out


def _rows(candidate: Dict, cases, trajectories, policy: OrderPolicy, *, g_override=None) -> List[Dict]:
    ev = Evaluator(candidate)
    rows = []
    for case in cases:
        coarse, fine, failed = trajectories[case["id"]]
        if failed:
            rows.append({"id": case["id"], "radius": case.get("radius"), "error": float("inf"),
                         "source_failed": failed, "refinement": float("inf"), "passed": False})
            continue
        params = case["parameters"]
        a = ev.g_coefficients(params) if g_override is None else g_override(params)
        rhs = lambda t, y, a=a: np.array([np.polyval(a[::-1], y[0])])
        rows.append(_row(ev, coarse, fine, rhs, case, policy))
    return rows


def _select(candidates, development, trajectories, policy: OrderPolicy):
    """The largest radius at which some order passes every development intervention, then the smallest such order."""
    evaluations = []
    selected = None
    for radius in policy.radii:
        cases = [c for c in development if c["radius"] == radius]
        for order in sorted(candidates):
            rows = _rows(candidates[order], cases, trajectories, policy)
            passed = bool(rows) and all(r["passed"] and r["error"] <= policy.development_margin
                                        * policy.prediction_tolerance for r in rows)
            evaluations.append({"radius": radius, "order": order, "passed": passed,
                                "max_error": max(r["error"] for r in rows)})
            if passed and selected is None:
                selected = (radius, order)
        if selected is not None:
            break
    return selected, evaluations


def _failed_law_coefficients(task: Dict, candidate: Dict):
    """The polynomial in z of the task's proposed drift, at given parameters, for the omission control."""
    real = mspec.load(task["model"])
    p = real.symbolic["parameters"]
    (name, text), = task["candidate_drift"].items()
    z = sp.Symbol(name, real=True)
    names = {str(s): s for s in p}
    names[name] = z
    expr = ln.expression(text, names)
    poly = sp.Poly(expr, z)
    fns = [sp.lambdify(p, c, modules="numpy") for c in poly.all_coeffs()[::-1]]
    order = [str(s) for s in p]
    base = dict(real.params)

    def coefficients(params):
        merged = {**base, **params}
        return np.array([float(f(*[merged[k] for k in order])) for f in fns])
    return coefficients


def _plain_projection(candidate: Dict) -> Dict:
    """The same reduced drift with h replaced by the linear coordinate: no correction of the fibres."""
    unit = tuple(int(i == 0) for i in range(len(next(iter(candidate["h"])))))
    return {**candidate, "h": {unit: 1.0}}


def _strongest_nonlinear(candidate: Dict, radius: float, length: float, end_factor: float):
    best, value = None, 0.0
    for key, c in candidate["g"].items():
        if sum(key) < 2 or key[0] < 2:
            continue
        size = abs(c) * (radius * length) ** key[0] * (radius * end_factor) ** sum(key[1:])
        if size > value:
            best, value = key, size
    return best


def order_recipes(library) -> List[Dict]:
    """Verified order closures of a library. They do not reach `MoveLibrary.recipes`, whose engine they do not carry."""
    with library.connect() as db:
        rows = db.execute("SELECT recipe FROM moves").fetchall()
    return [recipe for recipe in (json.loads(r[0]) for r in rows) if recipe.get("engine") == ORDER_ENGINE]


def close_modulo_order(task: Dict, library=None, policy: Optional[OrderPolicy] = None) -> Dict:
    """Derive a closure, select an order and radius on development, freeze, test once on withheld, and retain it."""
    policy = policy or OrderPolicy()
    started = time.monotonic()
    report: Dict = {"engine": ORDER_ENGINE, "policy": asdict(policy), "status": "unresolved",
                    "input_sha256": ln.digest(task), "novelty_established": False,
                    "validation_scope": "deterministic drift with the noise of the specification set to zero, on the "
                                        "declared interventions; the selection probabilities of the failed law are "
                                        "not re-derived",
                    "failed_law": task.get("failed_law"), "record": task.get("record"),
                    "write_point": task.get("write_point")}
    recipe = None
    try:
        real = mspec.load(task["model"])
        development, withheld = task["development"], task["withheld"]
        wp = task["write_point"]
        length = wp["scale"]["scale"]
        top = max(policy.order_list())
        closure = _solve_closure(real, wp["state"], wp["parameters"], wp["swept"], top, policy)
        report["other_spectrum"] = closure["other_spectrum"]
        report["smallest_denominator"] = closure["smallest_denominator"]
        report["lambda_critical"] = closure["lambda_critical"]
        report["resonance_at_order"] = closure["resonance_at_order"]
        report["g_low_orders"] = {str(k): v for k, v in closure["g"].items() if sum(k) <= 3}
        candidates = {order: truncate(closure, order) for order in policy.order_list() if order <= closure["order"]}
        development_runs = _source_trajectories(real, development)
        selected, evaluations = _select(candidates, development, development_runs, policy)
        report["development_evaluations"] = evaluations
        report["candidate_evaluations"] = len(evaluations)
        if selected is None:
            report["status"] = "no_order_passes_development"
        else:
            radius, order = selected
            candidate = candidates[order]
            frozen = candidate_record(candidate)
            identity = verify_modulo_order(task, candidate, policy)
            report.update(selected_order=order, selected_radius=radius, identity_modulo_order=identity,
                          candidate=frozen, candidate_sha256=ln.digest(frozen),
                          withheld_commitment=ln.digest(withheld),
                          description=description_size(real.carrier.dim, 1, 1))
            if not identity["passed"]:
                report["status"] = "identity_modulo_order_failed"
            else:
                x, p = real.symbolic["variables"], real.symbolic["parameters"]
                model_key = ln.digest({"variables": list(map(str, x)), "parameters": list(map(str, p)),
                                       "drift": list(map(str, real.symbolic["drift"])), "carrier": real.carrier.kind,
                                       "engine": ORDER_ENGINE})
                if library is not None and not library.reserve_withheld(
                        model_key, [{k: r[k] for k in ("initial", "parameters")} for r in withheld]):
                    report["status"] = "withheld_already_used"
                else:
                    report.update(_withheld_test(task, real, candidate, closure, withheld, order, radius, length,
                                                 policy))
                    if report["status"] == "verified_smaller_description":
                        recipe = {"engine": ORDER_ENGINE, "source_dimension": real.carrier.dim, "retained_dimension": 1,
                                  "order": order, "reach": radius,
                                  "kind": (task.get("failed_law") or {}).get("write_class")}
    except (ArithmeticError, np.linalg.LinAlgError) as error:
        report.update(status="numerical_failure", reason=f"{type(error).__name__}: {error}")
    except (ValueError, KeyError, TypeError, mspec.SpecError, sp.PolynomialError) as error:
        report.update(status="unsupported_or_resonant", reason=f"{type(error).__name__}: {error}")
    report["seconds"] = round(time.monotonic() - started, 3)
    report["implementation_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report["learning_sha256"] = hashlib.sha256(Path(ln.__file__).read_bytes()).hexdigest()
    report = _finite(report)
    if library is not None:
        stored = {**report, "status": "verified_extension"} if recipe is not None else report
        kept = library.store(stored, recipe)
        report["attempt_id"] = kept["attempt_id"]
    return report


def _finite(value):
    """A JSON-safe copy that keeps a diverged error as the text "inf"; the libraries refuse non-finite numbers."""
    if isinstance(value, dict):
        return {str(k): _finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite(v) for v in value]
    if isinstance(value, np.ndarray):
        return _finite(value.tolist())
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else str(float(value))
    if isinstance(value, np.integer):
        return int(value)
    return value


def _withheld_test(task, real, candidate, closure, withheld, order, radius, length, policy: OrderPolicy) -> Dict:
    runs = _source_trajectories(real, withheld)
    rows = _rows(candidate, withheld, runs, policy)
    at = [r for r in rows if r["radius"] == radius]
    cases_at = [c for c in withheld if c["radius"] == radius]
    passed = all(r["passed"] for r in at)
    # the failed-law drift on the same observable
    proposed = _failed_law_coefficients(task, candidate)
    omission_rows = _rows(candidate, cases_at, runs, policy, g_override=proposed)
    omission = max(r["error"] for r in omission_rows)
    error = max(r["error"] for r in at)
    improvement = 1 - error / max(omission, 1e-30) if np.isfinite(omission) else 1.0
    # a coefficient of g changed: a test that cannot see this has no power at this radius
    wp = task["write_point"]
    key = _strongest_nonlinear(candidate, radius, length, wp["sweep_end_factor"])
    power = {"term": list(key) if key else None, "factor": policy.power_factor, "rejected": None}
    if key is not None:
        altered = {**candidate, "g": {**candidate["g"], key: candidate["g"][key] * policy.power_factor}}
        power_rows = _rows(altered, cases_at, runs, policy)
        power.update(max_error=max(r["error"] for r in power_rows),
                     rejected=not all(r["passed"] for r in power_rows))
    # the plain coordinate with the same drift: is the correction of the fibres needed?
    plain = _rows(_plain_projection(candidate), cases_at, runs, policy)
    ladder = [{"radius": r, "max_error": max(x["error"] for x in rows if x["radius"] == r)} for r in policy.radii]
    finite = [(l["radius"], l["max_error"]) for l in ladder if l["radius"] <= radius and np.isfinite(l["max_error"])
              and l["max_error"] > 1e-9]
    slope = (float(np.polyfit(np.log([f[0] for f in finite]), np.log([f[1] for f in finite]), 1)[0])
             if len(finite) >= 2 else None)
    out = {"withheld": {"rows": rows, "max_error_at_radius": error, "error_by_radius": ladder,
                        "scaling_slope_of_error_with_radius": slope},
           "omission_error": omission, "relative_improvement": improvement, "power_control": power,
           "plain_coordinate_error": max(r["error"] for r in plain),
           "reaches_written_states": radius >= max(policy.radii)}
    if not passed:
        status = "withheld_failed"
    elif omission < policy.minimum_omission:
        status = "insufficient_excitation"
    elif improvement < policy.minimum_relative_improvement:
        status = "no_improvement_on_the_failed_law_drift"
    elif not power["rejected"]:
        status = "no_test_power"
    else:
        status = "verified_smaller_description"
    return {**out, "status": status}


# ------------------------------------------------------------------------------------------------- the loop
def _record_of(departure: Dict, root: Path):
    """The whole record of a retained outcome, located in its JSONL file by its digest."""
    path = Path(departure["source"])
    path = path if path.is_absolute() else Path(root) / path
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.strip():
                record = json.loads(line)
                if ln.digest(_preserved(record)) == departure["record_sha256"]:
                    return record
    return None


def _preserved(record):
    if isinstance(record, float) and not np.isfinite(record):
        return {"nonfinite_value": str(record)}
    if isinstance(record, dict):
        return {k: _preserved(v) for k, v in record.items()}
    if isinstance(record, list):
        return [_preserved(v) for v in record]
    return record


def close_departures(library, *, root: Path = Path("."), policy: Optional[OrderPolicy] = None,
                     select=None, revisit: Tuple[str, ...] = ()) -> List[Dict]:
    """Every retained outcome gets a row: a law departure is given a task and a campaign (once), and the others keep the
    reason they are not a repair question. `select` restricts the campaign by its write class. `revisit` names earlier
    statuses to run again after a change of the code; the registry refuses any withheld set already consumed, and the
    earlier report stays in the attempts and in the new report."""
    policy = policy or OrderPolicy()
    rows = []
    for key, departure in library.departures():
        if departure.get("kind") != FAILED_LAW:
            rows.append({"departure": key, "target": departure.get("target"), "status": "no_failed_law",
                         "kind": departure.get("kind")})
            continue
        write_class = (departure.get("original_result") or {}).get("write_class")
        if select is not None and write_class not in select:
            rows.append({"departure": key, "target": departure.get("target"), "status": "outside_this_campaign",
                         "write_class": write_class})
            continue
        old = library.campaign(key)
        if old is not None and old.get("status") not in revisit:
            rows.append({"departure": key, "status": "already_processed", "previous_status": old.get("status")})
            continue
        library.claim_campaign(key)
        record = _record_of(departure, root)
        try:
            if record is None:
                raise ValueError("the record of this departure is no longer in its source file")
            task = learning_task(record, departure["target"], root=root, policy=policy)
        except (ValueError, KeyError, ArithmeticError, np.linalg.LinAlgError, mspec.SpecError,
                sp.PolynomialError, TypeError) as error:
            report = {"status": "no_task_generated", "reason": f"{type(error).__name__}: {error}",
                      "target": departure.get("target"), "novelty_established": False}
            report = library.store(report)
        else:
            report = close_modulo_order(task, library, policy)
        if old is not None:
            report = {**report, "earlier_report": {k: old.get(k) for k in ("status", "reason", "attempt_id")}}
        library.store_campaign(key, _finite(report))
        rows.append({"departure": key, "target": departure.get("target"), "record": departure.get("paper_id"),
                     "write_class": write_class, "status": report["status"]})
    return rows


def summary(rows: List[Dict]) -> Dict:
    return dict(Counter(r["status"] for r in rows))
