"""A collective: the mean-field equations of a population of N units, with a noise that falls as 1/N.

    dy = F(y; c, h) dt + sqrt(2 D(y)) dW,   D = (1/2N) sum_j nu_j nu_j^T a_j(y)   (transitions)
                                            D = diag(v(y)) / (2N)                   (a variance rate per variable)

y are the fractions (or concentrations, or collective amplitudes), c the control, h the bias, N the size (a declared
parameter: the number of units, molecules per source volume, or the thermal stability factor of a macrospin).

At the symmetric state (h = 0) the critical eigenvalue crosses zero at c*: with the critical vectors v, w (w.v = 1),
a = d lambda/dc, the bias on the critical mode h_s = w . dF/dh, the cubic coefficient b and D_s = w^T D w, a sweep
dc/dt = r gives P = Phi(pi^(1/4) h_s h / (D_s^(1/2) (a r)^(1/4))) (the law of the memory module's symmetric write), with
D_s proportional to 1/N; Lambda = a r / (|b| D_s) must be large. The law along the actual sweep (`law_history`) takes
lambda, h_s and D_s on the unbiased mean trajectory integrated through the sweep: when the crossing window
(r/a)^(1/2) is not narrow, or when the mean state lags the control (slow symmetric modes), the closed form is its slow
limit only.
"""
from __future__ import annotations

from math import pi, sqrt
from typing import Dict, List, Optional

import numpy as np
from scipy.optimize import brentq, root
from scipy.stats import norm

from ..memory import spec as memory_spec
from . import passage

SpecError = memory_spec.SpecError


def compile_expressions(exprs: List[str], variables: List[str], params: Dict[str, float]):
    """Expressions of the restricted parser in the variables and parameters, compiled for columns of states."""
    sp = memory_spec._sympy()
    vs = [sp.Symbol(v, real=True) for v in variables]
    ps = [sp.Symbol(p, real=True) for p in params]
    names = {**dict(zip(variables, vs)), **dict(zip(params, ps))}
    parsed = [memory_spec.parse_expression(str(e), names) for e in exprs]
    fn = sp.lambdify(vs + ps, parsed, modules="numpy")
    pnames = list(params)

    def call(Y: np.ndarray, p: Dict[str, float]) -> np.ndarray:
        cols = [Y[:, k] for k in range(Y.shape[1])]
        with np.errstate(divide="ignore", invalid="ignore"):
            vals = fn(*cols, *[p[k] for k in pnames])
        return np.stack([np.broadcast_to(np.asarray(v, float), cols[0].shape) for v in vals], axis=1)

    return call, [str(e) for e in parsed]


class Collective:
    kind = "collective"

    def __init__(self, spec: Dict):
        pop = memory_spec._require(spec, "population", dict)
        self.size = memory_spec._require(pop, "size", str)
        control = memory_spec._require(pop, "control", dict)
        self.control = memory_spec._require(control, "name", str)
        self.range = [memory_spec._number(v, "population.control.range")
                      for v in memory_spec._require(control, "range", list)]
        if len(self.range) != 2 or not self.range[0] < self.range[1]:
            raise SpecError("population.control.range must be an increasing pair")
        bias = pop.get("bias")
        self.bias = memory_spec._require(bias, "name", str) if isinstance(bias, dict) else None
        noise = memory_spec._require(pop, "noise", dict)
        mem = {k: v for k, v in spec.items() if k not in ("schema", "kind", "population")}
        mem.update(schema=memory_spec.SCHEMA, kind="equations", noise=0.0)
        trans = noise.get("transitions")
        variables = memory_spec._require(memory_spec._require(spec, "carrier", dict), "variables", list)
        if "drift" not in spec:
            if not trans:
                raise SpecError("a collective needs a drift or a list of transitions")
            mem["drift"] = {v: " + ".join(f"({t['change'][v]})*({t['rate']})" for t in trans if v in t["change"]) or "0"
                            for v in variables}
        self.real = memory_spec.load(mem)
        self.variables = list(self.real.variables)
        self.n = len(self.variables)
        self.p = dict(self.real.params)
        for name in (self.size, self.control, self.bias):
            if name is not None and name not in self.p:
                raise SpecError(f"{name!r} must be a declared parameter")
        if trans:
            if not isinstance(trans, list) or not all(isinstance(t, dict) and "change" in t and "rate" in t
                                                        for t in trans):
                raise SpecError("population.noise.transitions: a list of {change, rate}")
            nus = []
            for t in trans:
                if not set(t["change"]) <= set(self.variables):
                    raise SpecError("a transition changes an undeclared variable")
                nus.append([memory_spec._number(t["change"].get(v, 0), "transition change") for v in self.variables])
            self.nus = np.array(nus)
            self._rates, self.noise_text = compile_expressions([t["rate"] for t in trans], self.variables, self.p)
            self.noise_kind = "transitions"
        else:
            var = memory_spec._require(noise, "variance", dict)
            if set(var) != set(self.variables):
                raise SpecError("population.noise.variance gives one expression for every variable")
            self._rates, self.noise_text = compile_expressions([var[v] for v in self.variables], self.variables,
                                                               self.p)
            self.nus = np.eye(self.n)
            self.noise_kind = "variance"
        guess = pop.get("symmetric_state")
        self.guess = (np.array([memory_spec._number(v, "symmetric_state") for v in guess]) if guess is not None
                      else np.zeros(self.n))
        if len(self.guess) != self.n:
            raise SpecError("population.symmetric_state gives one value per variable")
        self.orthant = self.real.carrier.kind == "orthant"

    # ------------------------------------------------------------------ equations
    def params(self, c: float, h: float = 0.0, N: Optional[float] = None) -> Dict[str, float]:
        p = dict(self.p)
        p[self.control] = float(c)
        if self.bias is not None:
            p[self.bias] = float(h)
        if N is not None:
            p[self.size] = float(N)
        return p

    def F(self, Y, c: float, h: float = 0.0) -> np.ndarray:
        return self.real.drift(np.atleast_2d(np.asarray(Y, float)), self.params(c, h))

    def jac(self, y, c: float, h: float = 0.0, eps: float = 1e-7) -> np.ndarray:
        y = np.asarray(y, float)
        J = np.empty((self.n, self.n))
        for j in range(self.n):
            e = np.zeros(self.n)
            e[j] = eps * max(1.0, abs(y[j]))
            J[:, j] = (self.F(y + e, c, h)[0] - self.F(y - e, c, h)[0]) / (2 * e[j])
        return J

    def steady(self, y0, c: float, h: float = 0.0) -> np.ndarray:
        sol = root(lambda x: self.F(x, c, h)[0], np.asarray(y0, float), jac=lambda x: self.jac(x, c, h), tol=1e-13)
        return sol.x

    def rates(self, Y, c: float, h: float = 0.0, N: Optional[float] = None) -> np.ndarray:
        """Propensities of the transitions (or the variance rates) per unit size, rows of Y."""
        return np.maximum(self._rates(np.atleast_2d(np.asarray(Y, float)), self.params(c, h, N)), 0.0)

    def D(self, y, c: float, N: float, h: float = 0.0) -> np.ndarray:
        a = self.rates(y, c, h, N)[0]
        return (self.nus.T * a) @ self.nus / (2.0 * N)

    # ------------------------------------------------------------------ reduction at the threshold
    def reduce(self, s_max: float = 1e-3) -> Dict:
        """c*, a, b, the critical vectors and h_s of the symmetric state (h = 0) in the control range."""
        lo, hi = self.range

        def lam(c):
            ev = np.linalg.eigvals(self.jac(self.steady(self.guess, c), c))
            return float(ev.real.max())
        try:
            cs = brentq(lam, lo, hi, xtol=1e-13 * max(1.0, abs(hi)))
        except ValueError as error:
            raise SpecError("the symmetric state does not lose stability inside population.control.range") from error
        q = self.steady(self.guess, cs)
        J = self.jac(q, cs)
        ev, V = np.linalg.eig(J)
        v = np.real(V[:, int(np.argmin(np.abs(ev)))])
        v /= np.linalg.norm(v)
        evl, W = np.linalg.eig(J.T)
        w = np.real(W[:, int(np.argmin(np.abs(evl)))])
        w /= w @ v
        dc = 1e-6 * max(1.0, abs(cs))
        a = (lam(cs + dc) - lam(cs - dc)) / (2 * dc)
        hs = float(w @ (self.F(q, cs, 1e-6)[0] - self.F(q, cs, -1e-6)[0]) / 2e-6) if self.bias else 0.0
        if hs < 0:                                  # orient the critical coordinate along the bias
            v, w, hs = -v, -w, -hs
        P = np.eye(self.n) - np.outer(v, w)
        ss = np.linspace(-s_max, s_max, 9)
        ss = ss[ss != 0]
        fs = []
        if self.n > 1:
            U = np.linalg.svd(P)[0][:, : self.n - 1]
        for s in ss:
            if self.n > 1:
                sol = root(lambda z: U.T @ (P @ self.F(q + s * v + U @ z, cs)[0]), np.zeros(self.n - 1), tol=1e-14)
                qq = q + s * v + U @ sol.x
            else:
                qq = q + s * v
            fs.append(float(w @ self.F(qq, cs)[0]))
        coef = np.linalg.lstsq(np.vstack([ss, ss ** 3, ss ** 5]).T, np.array(fs), rcond=None)[0]
        return {"c_star": float(cs), "y_sym": q, "a": float(a), "b": float(-coef[1]), "v": v, "w": w, "h_s": hs}

    def D_s(self, red: Dict, N: float) -> float:
        return float(red["w"] @ self.D(red["y_sym"], red["c_star"], N) @ red["w"])

    # ------------------------------------------------------------------ laws of the swept write
    def law(self, red: Dict, h: float, N: float, r: float) -> Dict[str, float]:
        Ds = self.D_s(red, N)
        z = pi ** 0.25 * red["h_s"] * h / (sqrt(Ds) * (red["a"] * r) ** 0.25)
        return {"P": float(norm.cdf(z)), "probit": float(z), "D_s": Ds,
                "Lambda": red["a"] * r / (abs(red["b"]) * Ds) if red["b"] != 0 else float("inf"),
                "window": sqrt(r / red["a"])}

    def law_history(self, red: Dict, h: float, N: float, r: float, c0: float, c1: float, points: int = 2000,
                    substeps: int = 10) -> Dict[str, float]:
        """The law along the unbiased mean trajectory of the sweep from the steady state at c0: probit = m/sigma with
        the bias and the noise on the critical mode accumulated through the passage, the start in the biased
        equilibrium (mean h_s/|lambda| and variance D_s/|lambda| at c0). The probit is linear in h and in N^(1/2)."""
        m, s2 = self.history_factors(red, N, r, c0, c1, points, substeps)
        z = h * m / sqrt(s2)
        return {"P": float(norm.cdf(z)), "probit": float(z), "m": m, "s2": s2}

    def history_factors(self, red: Dict, N: float, r: float, c0: float, c1: float, points: int = 2000,
                        substeps: int = 10):
        """(m, s2) of the law along the passage at the size N: probit = h m / s2^(1/2)."""
        w = red["w"]
        cs = np.linspace(c0, c1, points)
        ts = (cs - c0) / r
        y = self.steady(red["y_sym"], c0)
        lam, hs, Ds = np.empty(points), np.empty(points), np.empty(points)
        f = lambda y, c: self.F(y, c)[0]                                              # noqa: E731
        for k, c in enumerate(cs):
            if k > 0:
                dt = (ts[k] - ts[k - 1]) / substeps
                cc = cs[k - 1]
                for _ in range(substeps):
                    k1 = f(y, cc)
                    k2 = f(y + dt / 2 * k1, cc + r * dt / 2)
                    k3 = f(y + dt / 2 * k2, cc + r * dt / 2)
                    k4 = f(y + dt * k3, cc + r * dt)
                    y = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
                    cc += r * dt
            ev, V = np.linalg.eig(self.jac(y, c))
            lam[k] = ev[int(np.argmax(np.abs(np.real(V).T @ w)))].real
            hs[k] = float(w @ (self.F(y, c, 1e-6)[0] - self.F(y, c, -1e-6)[0]) / 2e-6)
            Ds[k] = float(w @ self.D(y, c, N) @ w)
        Lam = np.concatenate([[0.0], np.cumsum(0.5 * (lam[1:] + lam[:-1]) * np.diff(ts))])
        trap = getattr(np, "trapezoid", None) or np.trapz
        L0 = Lam.min()
        m = hs[0] / abs(lam[0]) * np.exp(L0) + trap(hs * np.exp(-(Lam - L0)), ts)
        s2 = Ds[0] / abs(lam[0]) * np.exp(2 * L0) + trap(2 * Ds * np.exp(-2 * (Lam - L0)), ts)
        return float(m), float(s2)

    # ------------------------------------------------------------------ stochastic collective
    def _step(self, Y, c, h, N, dt, rng):
        a = self.rates(Y, c, h, N)
        if self.noise_kind == "transitions":
            # chemical Langevin step: each transition with its own noise
            k = a * dt + np.sqrt(a * dt / N) * rng.standard_normal(a.shape)
            Y = Y + (self.F(Y, c, h) * dt) + (k - a * dt) @ self.nus
        else:
            Y = Y + self.F(Y, c, h) * dt + np.sqrt(a * dt / N) * rng.standard_normal(Y.shape)
        return np.maximum(Y, 0.0) if self.orthant else Y

    def sweep(self, red: Dict, h: float, N: float, r: float, c0: float, c1: float, replicas: int, seed: int = 1,
              dt: float = 1e-2) -> Dict[str, float]:
        """Replicas of the stochastic collective, relaxed at c0 (ten relaxation times), swept to c1 at the rate r; the
        fraction whose critical coordinate ends with the sign of the bias."""
        rng = np.random.default_rng(seed)
        y0 = self.steady(red["y_sym"], c0, h)
        lam0 = abs(float(np.linalg.eigvals(self.jac(y0, c0, h)).real.max()))
        Y = np.repeat(y0[None, :], replicas, axis=0)
        for _ in range(int(np.ceil(10.0 / (lam0 * dt)))):
            Y = self._step(Y, c0, h, N, dt, rng)
        steps = int(np.ceil((c1 - c0) / (r * dt)))
        h_dt = (c1 - c0) / (r * steps)
        c = c0
        for _ in range(steps):
            Y = self._step(Y, c, h, N, h_dt, rng)
            c += r * h_dt
        s = (Y - self.steady(red["y_sym"], c1)) @ red["w"]
        P = float(np.mean(np.sign(s) == np.sign(h)) + 0.5 * np.mean(s == 0))
        return {"P": P, "stderr": sqrt(max(P * (1 - P), 1e-12) / replicas), "replicas": replicas}

    # ------------------------------------------------------------------ passage after a step
    def step_exponent(self, red: Dict, c0: float, c1: float, T: float, points: int = 4000) -> Dict:
        """The unbiased mean trajectory after a step c0 -> c1 (from the steady state at c0), the dimension d of the
        leading eigenspace at the end, its left basis, and Lambda(t) from the rate of that eigenspace along the
        trajectory."""
        f = lambda y, c=c1: self.F(y, c)[0]                                          # noqa: E731
        t = np.linspace(0.0, T, points)
        y = self.steady(red["y_sym"], c0)
        ys = [y]
        dt = t[1] - t[0]
        for _ in range(points - 1):
            for _s in range(5):
                h = dt / 5
                k1 = f(y)
                k2 = f(y + h / 2 * k1)
                k3 = f(y + h / 2 * k2)
                k4 = f(y + h * k3)
                y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            ys.append(y)
        ys = np.array(ys)
        Jend = self.jac(ys[-1], c1)
        d = passage.leading_dimension(Jend)
        ev, V = np.linalg.eig(Jend)
        order = np.argsort(-ev.real)[:d]
        basis = np.real(np.column_stack([V[:, order[0]].real, V[:, order[0]].imag])) if (
            d == 2 and abs(ev[order[0]].imag) > 1e-12) else np.real(V[:, order])
        q, _ = np.linalg.qr(basis)
        lam = np.empty(points)
        for k in range(points):
            evk, Vk = np.linalg.eig(self.jac(ys[k], c1))
            lam[k] = evk[int(np.argmax(np.linalg.norm(q.T @ np.real(Vk), axis=0)))].real
        return {"t": t, "mean": ys, "d": d, "basis": q, "exponent": passage.Exponent(t, lam),
                "rate_end": float(lam[-1])}

    def step_times(self, red: Dict, N: float, c0: float, c1: float, threshold: float, T: float, replicas: int,
                   step: Dict, seed: int = 1, dt: float = 1e-2, h: float = 0.0) -> np.ndarray:
        """First times at which the projection of y - mean(t) on the leading eigenspace exceeds `threshold`, after the
        step, for replicas relaxed at c0 (ten relaxation times)."""
        rng = np.random.default_rng(seed)
        y0 = self.steady(red["y_sym"], c0, h)
        lam0 = abs(float(np.linalg.eigvals(self.jac(y0, c0, h)).real.max()))
        Y = np.repeat(y0[None, :], replicas, axis=0)
        for _ in range(int(np.ceil(10.0 / (lam0 * dt)))):
            Y = self._step(Y, c0, h, N, dt, rng)
        times = np.full(replicas, np.inf)
        tm, mean, q = step["t"], step["mean"], step["basis"]
        for k in range(int(np.ceil(T / dt))):
            Y = self._step(Y, c1, h, N, dt, rng)
            tt = (k + 1) * dt
            m = np.array([np.interp(tt, tm, mean[:, j]) for j in range(self.n)])
            hit = (np.linalg.norm((Y - m) @ q, axis=1) > threshold) & ~np.isfinite(times)
            times[hit] = tt
            if np.isfinite(times).all():
                break
        return times
