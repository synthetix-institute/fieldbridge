"""Bodies that grow and divide: what a lineage needs from each kind of specification.

Every body provides
    n, f(Q, L), jac(q, L), steady(q0, L)       the drift at size L, for the reduction at the threshold;
    grow(Q, L0, L1, rng, noise)                the rows of Q integrated while the size grows from L0 to L1, with the
                                               noise of the specification times `noise` (0: noise-free);
    divide(Q, rng)                             the state of the followed daughter of each row;
    diffusion(q, L)                            the diffusion matrix D (dq = F dt + sqrt(2D) dW) at q and size L;
    growth                                     the growth law of the size.

Three kinds:
    equations   an equations body of the memory module with the size as a declared parameter; noise as a number
                (additive on every variable) or one expression per variable for the variance rate (for example the
                chemical Langevin noise of a reaction volume); division keeps or halves each variable;
    field       variables on N cells (key "cells") of a body of length L in one dimension, with local reactions, a mobility times
                the Laplacian (constant or depending on the local variables), faces with zero flux or with an
                extrapolation length; division cuts the body in the middle;
    chain       an inextensible elastic filament of N rigid segments under a dead compressive load, with overdamped
                anisotropic drag and thermal noise; division severs it in the middle.
"""
from __future__ import annotations

from math import exp, log
from typing import Dict, List, Optional

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import root

from ..memory import spec as memory_spec

SpecError = memory_spec.SpecError


class Growth:
    """L(t) = L0 e^{rt} (exponential) or L0 + r t (linear)."""

    def __init__(self, law: str, rate: float):
        if law not in ("exponential", "linear"):
            raise SpecError("lineage.growth.law must be exponential or linear")
        if rate <= 0:
            raise SpecError("lineage.growth.rate must be positive")
        self.law, self.rate = law, float(rate)

    def size(self, L0: float, t):
        return L0 * np.exp(self.rate * t) if self.law == "exponential" else L0 + self.rate * t

    def duration(self, L0: float, L1: float) -> float:
        return log(L1 / L0) / self.rate if self.law == "exponential" else (L1 - L0) / self.rate

    def ramp(self, L: float) -> float:
        """dL/dt at size L."""
        return self.rate * L if self.law == "exponential" else self.rate

    def advance(self, L: float, dt: float) -> float:
        return L * exp(self.rate * dt) if self.law == "exponential" else L + self.rate * dt


def _compile(exprs: Dict[str, str], variables: List[str], params: Dict[str, float], what: str):
    """Expressions in the variables and parameters, compiled for columns of states."""
    sp = memory_spec._sympy()
    vs = [sp.Symbol(v, real=True) for v in variables]
    ps = [sp.Symbol(p, real=True) for p in params]
    names = {**dict(zip(variables, vs)), **dict(zip(params, ps))}
    parsed = [memory_spec.parse_expression(exprs[v], names) for v in variables]
    fn = sp.lambdify(vs + ps, parsed, modules="numpy")
    pnames = list(params)

    def call(cols: List[np.ndarray], p: Dict[str, float]) -> List[np.ndarray]:
        with np.errstate(divide="ignore", invalid="ignore"):
            vals = fn(*cols, *[p[k] for k in pnames])
        return [np.broadcast_to(np.asarray(v, float), cols[0].shape) for v in vals]

    return call, [str(e) for e in parsed]


class _Common:
    """Numerical Jacobian, steady state and the size parameter."""

    def params(self, L: float) -> Dict[str, float]:
        p = dict(self.p)
        p[self.size] = float(L)
        return p

    def jac(self, q, L, h: float = 1e-7) -> np.ndarray:
        q = np.asarray(q, float)
        J = np.empty((self.n, self.n))
        for j in range(self.n):
            e = np.zeros(self.n)
            e[j] = h * max(1.0, abs(q[j]))
            J[:, j] = (self.f(q + e, L)[0] - self.f(q - e, L)[0]) / (2 * e[j])
        return J

    def steady(self, q0, L) -> np.ndarray:
        sol = root(lambda x: self.f(x, L)[0], np.asarray(q0, float), jac=lambda x: self.jac(x, L), tol=1e-13)
        return sol.x


# ------------------------------------------------------------------------------------------------ equations
class EquationsBody(_Common):
    kind = "equations"

    def __init__(self, spec: Dict, size: str, growth: Growth, noise, division: Dict, dt: float):
        mem = {k: v for k, v in spec.items() if k not in ("schema", "body", "lineage")}
        mem.update(schema=memory_spec.SCHEMA, kind="equations", noise=0.0)
        self.real = memory_spec.load(mem)
        self.p, self.size, self.growth, self.dt = dict(self.real.params), size, growth, dt
        if size not in self.p:
            raise SpecError("lineage.size must be a declared parameter")
        self.variables = list(self.real.variables)
        self.n = len(self.variables)
        self.orthant = self.real.carrier.kind == "orthant"
        if isinstance(noise, dict):
            if set(noise) != set(self.variables):
                raise SpecError("lineage.noise must give one expression for every variable, or be a number")
            self._noise, self.noise_text = _compile(noise, self.variables, self.p, "noise")
        else:
            D = memory_spec._number(noise, "lineage.noise")
            self._noise, self.noise_text = (lambda cols, p, D=D: [np.full(cols[0].shape, D)] * len(cols)), [str(D)]
        rules = division.get("variables", {v: "keep" for v in self.variables})
        if set(rules) != set(self.variables) or not set(rules.values()) <= {"keep", "halve"}:
            raise SpecError("lineage.division.variables must give keep or halve for every variable")
        self.halve = np.array([rules[v] == "halve" for v in self.variables])
        self.partition = division.get("partition", "exact")
        self.volume = division.get("volume")
        if self.partition not in ("exact", "binomial"):
            raise SpecError("lineage.division.partition must be exact or binomial")
        if self.partition == "binomial" and self.volume not in self.p:
            raise SpecError("a binomial partition needs lineage.division.volume, a declared parameter")

    def f(self, Q, L):
        return self.real.drift(np.atleast_2d(np.asarray(Q, float)), self.params(L))

    def _D(self, Q: np.ndarray, L: float) -> np.ndarray:
        cols = [Q[:, k] for k in range(self.n)]
        return np.maximum(np.stack(self._noise(cols, self.params(L)), axis=-1), 0.0)

    def diffusion(self, q, L) -> np.ndarray:
        return np.diag(self._D(np.atleast_2d(np.asarray(q, float)), L)[0])

    def grow(self, Q, L0, L1, rng, noise: float = 1.0, dt: Optional[float] = None):
        """Euler-Maruyama (Ito) with the noise of the specification times `noise`."""
        Q = np.array(Q, float)
        T = self.growth.duration(L0, L1)
        steps = max(1, int(np.ceil(T / (dt or self.dt))))
        h = T / steps
        L = L0
        for _ in range(steps):
            drift = self.f(Q, L)
            if noise > 0:
                Q = Q + h * drift + np.sqrt(2.0 * noise * h * self._D(Q, L)) * rng.standard_normal(Q.shape)
            else:
                Q = Q + h * drift
            if self.orthant:
                np.maximum(Q, 0.0, out=Q)
            L = self.growth.advance(L, h)
        return Q

    def divide(self, Q, rng, noise: float = 1.0):
        Q = np.array(Q, float)
        Q[:, self.halve] *= 0.5
        if self.partition == "binomial" and noise > 0:
            # binomial partition of the molecules: variance N/4 per species, x/(4 Omega) in concentrations
            omega = self.p[self.volume]
            Q[:, self.halve] += np.sqrt(np.maximum(Q[:, self.halve], 0.0) / (2.0 * omega)) \
                * rng.standard_normal(Q[:, self.halve].shape) * np.sqrt(noise)
        return np.maximum(Q, 0.0) if self.orthant else Q


# ------------------------------------------------------------------------------------------------ field
class FieldBody(_Common):
    """Variables u_v(x) on N cells of x in [0, 1] (z = L x); the state is variable-major: (u_1 cells, u_2 cells, ...)."""
    kind = "field"

    def __init__(self, spec: Dict, size: str, growth: Growth, noise: float, dt: float):
        fld = memory_spec._require(spec, "cells", dict)
        self.N = int(memory_spec._require(fld, "number", int))
        if self.N < 8 or self.N % 2:
            raise SpecError("cells.number must be an even number of at least 8")
        self.variables = memory_spec._require(fld, "variables", list)
        self.p = {k: memory_spec._number(v, f"parameter {k}")
                  for k, v in memory_spec._require(spec, "parameters", dict).items()}
        for name in list(self.variables) + list(self.p):
            if not str(name).isidentifier() or name in memory_spec.RESERVED or name in memory_spec.FUNCTIONS:
                raise SpecError(f"Invalid or reserved name {name!r}")
        if size not in self.p:
            raise SpecError("lineage.size must be a declared parameter")
        self.size, self.growth, self.dt = size, growth, dt
        self.k = len(self.variables)
        self.n = self.k * self.N
        reaction = memory_spec._require(fld, "reaction", dict)
        if set(reaction) != set(self.variables):
            raise SpecError("cells.reaction must give one expression for every variable")
        self._reaction, self.reaction_text = _compile(reaction, self.variables, self.p, "reaction")
        mob = fld.get("mobility", {v: "1" for v in self.variables})
        if set(mob) != set(self.variables):
            raise SpecError("cells.mobility must give one expression for every variable")
        self._mobility, self.mobility_text = _compile({v: str(mob[v]) for v in self.variables}, self.variables,
                                                      self.p, "mobility")
        ref = fld.get("reference", {})
        self.u0 = np.array([constant(ref.get(v, 0.0), f"cells.reference.{v}") for v in self.variables])
        faces = fld.get("faces", "zero_flux")
        if faces == "zero_flux":
            self.d = None
        else:
            d = memory_spec._require(faces, "extrapolation")
            self.d = self.p[d] if isinstance(d, str) else memory_spec._number(d, "cells.faces.extrapolation")
            if self.d <= 0:
                raise SpecError("cells.faces.extrapolation must be positive")
        self.mode = spec.get("lineage", {}).get("growth", {}).get("mode", "uniform")
        if self.mode not in ("uniform", "faces"):
            raise SpecError("lineage.growth.mode must be uniform or faces")
        self.eps = memory_spec._number(noise, "lineage.noise")
        self.x = (np.arange(self.N) + 0.5) / self.N

    # -- the Laplacian (in z) of the deviation from the face value, with the faces
    def _ghost(self, L: float) -> float:
        if self.d is None:
            return 1.0
        hz = L / self.N
        return (1 - hz / (2 * self.d)) / (1 + hz / (2 * self.d))

    def _lap(self, E: np.ndarray, L: float) -> np.ndarray:
        g, hz = self._ghost(L), L / self.N
        left = np.concatenate([g * E[:, :1], E[:, :-1]], 1)
        right = np.concatenate([E[:, 1:], g * E[:, -1:]], 1)
        return (left - 2 * E + right) / hz ** 2

    def _split(self, Q: np.ndarray) -> List[np.ndarray]:
        return [Q[:, i * self.N:(i + 1) * self.N] for i in range(self.k)]

    def _terms(self, Q: np.ndarray, L: float):
        U = self._split(np.atleast_2d(np.asarray(Q, float)))
        p = self.params(L)
        R = self._reaction(U, p)
        M = self._mobility(U, p)
        laps = [self._lap(U[i] - self.u0[i], L) for i in range(self.k)]
        adv = [self._advection(U[i]) for i in range(self.k)] if self.mode == "faces" else [0.0] * self.k
        return U, R, M, laps, adv

    def _advection(self, u: np.ndarray) -> np.ndarray:
        """r (x - 1/2) u_x: deposition on both faces, in the scaled coordinate (central differences)."""
        left = np.concatenate([u[:, :1], u[:, :-1]], 1)
        right = np.concatenate([u[:, 1:], u[:, -1:]], 1)
        return self.growth.rate * (self.x - 0.5) * (right - left) * (self.N / 2.0)

    def f(self, Q, L):
        U, R, M, laps, adv = self._terms(Q, L)
        return np.concatenate([R[i] + M[i] * laps[i] + adv[i] for i in range(self.k)], 1)

    def diffusion(self, q, L) -> np.ndarray:
        return np.full(self.n, self.eps * self.N / L)

    def _banded(self, L: float, c: float) -> np.ndarray:
        g, hz, N = self._ghost(L), L / self.N, self.N
        kk = c / hz ** 2
        ab = np.zeros((3, N))
        ab[0, 1:] = -kk
        ab[1, :] = 1 + 2 * kk
        ab[1, 0] -= kk * g
        ab[1, -1] -= kk * g
        ab[2, :-1] = -kk
        return ab

    def grow(self, Q, L0, L1, rng, noise: float = 1.0, dt: Optional[float] = None):
        """A linearly implicit step: c_v times the Laplacian implicit, the rest explicit; c_v is the largest mobility
        of variable v over the rows at the start of the step, which keeps the step stable for any dt."""
        Q = np.array(Q, float)
        T = self.growth.duration(L0, L1)
        steps = max(1, int(np.ceil(T / (dt or self.dt))))
        h = T / steps
        L = L0
        for _ in range(steps):
            U, R, M, laps, adv = self._terms(Q, L)
            out = []
            for i in range(self.k):
                c = max(float(np.max(M[i])), 1e-12)
                u = U[i] + h * (R[i] + (M[i] - c) * laps[i] + adv[i])
                if noise > 0 and self.eps > 0:
                    u = u + np.sqrt(2.0 * noise * self.eps * self.N / L * h) * rng.standard_normal(u.shape)
                u = self.u0[i] + solve_banded((1, 1), self._banded(L, h * c), (u - self.u0[i]).T).T
                out.append(u)
            Q = np.concatenate(out, 1)
            L = self.growth.advance(L, h)
        return Q

    def divide(self, Q, rng, noise: float = 1.0):
        """Cut in the middle and follow one half (left or right at random), refined to N cells; the cut becomes a face
        of the same kind."""
        h = self.N // 2
        right = rng.random(len(Q)) < 0.5
        parts = []
        for u in self._split(np.asarray(Q, float)):
            parts.append(np.repeat(np.where(right[:, None], u[:, h:], u[:, :h]), 2, axis=1))
        return np.concatenate(parts, 1)


def constant(value, what: str) -> float:
    """A number, or a constant expression of numbers and pi (for example "pi/2")."""
    if isinstance(value, str):
        return float(memory_spec.parse_expression(value, {}))
    return memory_spec._number(value, what)


# ------------------------------------------------------------------------------------------------ chain
class ChainBody(_Common):
    """Bond angles theta_1..theta_N of an inextensible chain (segment length h = L/N), the first end at height 0 and the
    last held at height 0 by a constraint, both ends free to slide along the axis of the load and free to turn. Energy
    U = (kappa/h) sum (1 - cos(theta_{j+1} - theta_j)) + F (x_N - x_0). Friction metric M = sum_i J_i^T Z_i J_i with
    Z_i = h w_i [zeta_perp (1 - t t) + zeta_par t t]; x_0 is a free coordinate and is solved for, not kept."""
    kind = "chain"
    RELAX = 2.0e3        # relaxation rate of the constraint error (rounding only)

    def __init__(self, spec: Dict, size: str, growth: Growth, noise: float, dt: float):
        ch = memory_spec._require(spec, "chain", dict)
        self.N = int(memory_spec._require(ch, "segments", int))
        if self.N < 4 or self.N % 2:
            raise SpecError("chain.segments must be an even number of at least 4")
        num = lambda key: memory_spec._number(memory_spec._require(ch, key), f"chain.{key}")  # noqa: E731
        self.kappa, self.F, self.zperp, self.zpar = num("bending"), num("load"), num("drag_perp"), num("drag_par")
        if min(self.kappa, self.F, self.zperp, self.zpar) <= 0:
            raise SpecError("chain bending, load and drags must be positive")
        self.p = {k: memory_spec._number(v, f"parameter {k}")
                  for k, v in memory_spec._require(spec, "parameters", dict).items()}
        if size not in self.p:
            raise SpecError("lineage.size must be a declared parameter")
        self.size, self.growth, self.dt = size, growth, dt
        self.kT = memory_spec._number(noise, "lineage.noise")
        self.n = self.N
        w = np.ones(self.N + 1)
        w[0] = w[-1] = 0.5
        self.wb = w
        self.tri = np.tril(np.ones((self.N + 1, self.N)), -1)

    def _metric(self, T: np.ndarray, L: float):
        R, N = T.shape
        h = L / N
        s, c = np.sin(T), np.cos(T)
        Jx = np.zeros((R, N + 1, N + 1))
        Jy = np.zeros((R, N + 1, N + 1))
        Jx[:, :, 0] = 1.0
        Jx[:, :, 1:] = -h * s[:, None, :] * self.tri[None]
        Jy[:, :, 1:] = h * c[:, None, :] * self.tri[None]
        ang = np.empty((R, N + 1))
        ang[:, 0], ang[:, -1] = T[:, 0], T[:, -1]
        ang[:, 1:-1] = np.arctan2(s[:, :-1] + s[:, 1:], c[:, :-1] + c[:, 1:])
        tx, ty = np.cos(ang), np.sin(ang)
        Pt = tx[:, :, None] * Jx + ty[:, :, None] * Jy
        W = h * self.wb
        tr = lambda X: np.swapaxes(X * W[None, :, None], 1, 2)     # noqa: E731
        M = self.zperp * (tr(Jx) @ Jx + tr(Jy) @ Jy) + (self.zpar - self.zperp) * (tr(Pt) @ Pt)
        return M, Jx, Jy, tx, ty, W

    def _grad_U(self, T: np.ndarray, L: float) -> np.ndarray:
        R, N = T.shape
        h = L / N
        d = np.sin(np.diff(T, axis=1))
        g = np.zeros((R, N + 1))
        g[:, 1:] -= self.F * h * np.sin(T)
        g[:, 1:-1] -= self.kappa / h * d
        g[:, 2:] += self.kappa / h * d
        return g

    def _hess_bend(self, T: np.ndarray, L: float) -> np.ndarray:
        R, N = T.shape
        cd = self.kappa / (L / N) * np.cos(np.diff(T, axis=1))
        H = np.zeros((R, N + 1, N + 1))
        j = np.arange(N - 1)
        H[:, 1 + j, 1 + j] += cd
        H[:, 2 + j, 2 + j] += cd
        H[:, 1 + j, 2 + j] -= cd
        H[:, 2 + j, 1 + j] -= cd
        return H

    def random_force(self, T: np.ndarray, L: float, rng, dt: float, kT: float) -> np.ndarray:
        """sum_i J_i^T Z_i^(1/2) sqrt(2 kT/dt) xi_i (covariance 2 kT M / dt)."""
        R, N = T.shape
        _, Jx, Jy, tx, ty, W = self._metric(T, L)
        xi = rng.standard_normal((R, N + 1, 2))
        xt = tx * xi[:, :, 0] + ty * xi[:, :, 1]
        xn = -ty * xi[:, :, 0] + tx * xi[:, :, 1]
        sq = np.sqrt(W)
        fx = sq * (np.sqrt(self.zpar) * xt * tx - np.sqrt(self.zperp) * xn * ty)
        fy = sq * (np.sqrt(self.zpar) * xt * ty + np.sqrt(self.zperp) * xn * tx)
        return np.sqrt(2 * kT / dt) * (np.einsum("rin,ri->rn", Jx, fx) + np.einsum("rin,ri->rn", Jy, fy))

    def velocity(self, T: np.ndarray, L: float, force=None, dt: float = 0.0, implicit: bool = False) -> np.ndarray:
        R, N = T.shape
        h = L / N
        M = self._metric(T, L)[0]
        if implicit:
            M = M + dt * self._hess_bend(T, L)
        gg = np.zeros((R, N + 1))
        gg[:, 1:] = h * np.cos(T)
        gval = h * np.sin(T).sum(1)
        F = -self._grad_U(T, L) if force is None else force - self._grad_U(T, L)
        sol = np.linalg.solve(M, np.stack([F, gg], axis=-1))
        v, A = sol[..., 0], sol[..., 1]
        mu = -(np.einsum("rn,rn->r", gg, v) + self.RELAX * gval) / np.einsum("rn,rn->r", gg, A)
        return (v + mu[:, None] * A)[:, 1:]

    def f(self, Q, L):
        return self.velocity(np.atleast_2d(np.asarray(Q, float)), L)

    def steady(self, q0, L):
        return np.zeros(self.n)                       # the straight filament on the axis

    def diffusion(self, q, L) -> np.ndarray:
        """kT (M^-1 - A A^T / (grad g . A)) on the angles."""
        T = np.atleast_2d(np.asarray(q, float))
        M = self._metric(T, L)[0][0]
        gg = np.r_[0.0, L / self.N * np.cos(T[0])]
        Mi = np.linalg.inv(M)
        A = Mi @ gg
        return self.kT * (Mi - np.outer(A, A) / (gg @ A))[1:, 1:]

    def grow(self, Q, L0, L1, rng, noise: float = 1.0, dt: Optional[float] = None):
        """Steps linearly implicit in the bending energy; the midpoint scheme of Grassia, Hinch and Nitsche (J. Fluid
        Mech. 282, 373, 1995): the random force drawn at the start is used again at the midpoint, which supplies the
        drift kT div(M^-1) of the configuration-dependent metric (also without noise, for the same drift)."""
        T = np.array(Q, float)
        dur = self.growth.duration(L0, L1)
        steps = max(1, int(np.ceil(dur / (dt or self.dt))))
        h = dur / steps
        L = L0
        kT = noise * self.kT
        for _ in range(steps):
            f = self.random_force(T, L, rng, h, kT) if kT > 0 else None
            mid = T + 0.5 * h * self.velocity(T, L, f, h, implicit=True)
            T = T + h * self.velocity(mid, L, f, h, implicit=True)
            L = self.growth.advance(L, h)
        return T

    def divide(self, Q, rng, noise: float = 1.0):
        """Sever in the middle, follow one half (refined to N bonds) and turn it rigidly so that both ends lie on the axis
        of the load: sum sin(theta - a) = 0 with a = atan2(sum sin theta, sum cos theta)."""
        h = self.N // 2
        right = rng.random(len(Q)) < 0.5
        half = np.repeat(np.where(right[:, None], Q[:, h:], Q[:, :h]), 2, axis=1)
        a = np.arctan2(np.sin(half).sum(1), np.cos(half).sum(1))
        return half - a[:, None]
