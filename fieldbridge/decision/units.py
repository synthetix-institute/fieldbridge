"""Units with a diverse parameter, coupled through their mean: dx_i = f(x_i; a_i, X, h, c) dt + sqrt(2 D) dW_i.

The unit drift f is an expression in the unit variable, the diverse parameter a, the mean X of the population, the
bias h, the control c and further parameters. The parameter a has a Gaussian density (centre, width).

- A quantile sample (a_i = centre + width Phi^-1((i - 1/2)/M)) is symmetric when f is odd under (x, a, X, h) -> -(x, a,
  X, h): its collective pitchfork at c* is reduced with M representative units, and a sweep follows the swept law with
  D_s proportional to D/N.
- A random sample is not symmetric. At the reference profile x*(a) of the infinite population every unit feels the
  drift df/dX (X_sample - X_inf), so the sample carries a frozen bias of standard deviation
  s_q = |sum w df/dX / sum w df/dh| (Var_a x* / N)^(1/2) in units of h. When s_q exceeds the thermal spread
  sigma_th = (D_s (a r / pi)^(1/2))^(1/2) / h_s, the outcome is set by the sample: P = Phi(h / (sigma_th^2 + s_q^2)^(1/2))
  tends to Phi(h / s_q), whatever the sweep rate and D.
"""
from __future__ import annotations

from math import pi, sqrt
from typing import Dict

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

from ..memory import spec as memory_spec
from .collective import compile_expressions

SpecError = memory_spec.SpecError


class Units:
    kind = "units"

    def __init__(self, spec: Dict):
        pop = memory_spec._require(spec, "population", dict)
        unit = memory_spec._require(pop, "unit", dict)
        self.x = memory_spec._require(unit, "variable", str)
        self.mean = memory_spec._require(unit, "mean", str)
        het = memory_spec._require(pop, "heterogeneity", dict)
        self.a = memory_spec._require(het, "parameter", str)
        if het.get("density", "gaussian") != "gaussian":
            raise SpecError("units: only a gaussian heterogeneity is supported")
        self.centre = memory_spec._number(het.get("centre", 0.0), "heterogeneity.centre")
        self.width = memory_spec._number(memory_spec._require(het, "width"), "heterogeneity.width")
        control = memory_spec._require(pop, "control", dict)
        self.control = memory_spec._require(control, "name", str)
        self.range = [memory_spec._number(v, "population.control.range")
                      for v in memory_spec._require(control, "range", list)]
        self.bias = memory_spec._require(memory_spec._require(pop, "bias", dict), "name", str)
        self.noise = memory_spec._number(memory_spec._require(pop, "noise"), "population.noise")
        self.p = {k: memory_spec._number(v, f"parameter {k}")
                  for k, v in memory_spec._require(spec, "parameters", dict).items()}
        for name in (self.control, self.bias):
            if name not in self.p:
                raise SpecError(f"{name!r} must be a declared parameter")
        for name in (self.x, self.mean, self.a):
            if name in self.p or not name.isidentifier() or name in memory_spec.RESERVED:
                raise SpecError(f"invalid unit name {name!r}")
        expr = memory_spec._require(unit, "drift", str)
        names = [self.x, self.a, self.mean]
        self._f, self.drift_text = compile_expressions([expr], names, self.p)
        self.M = int(pop.get("representatives", 400))

    # ------------------------------------------------------------------ the unit
    def f(self, x, a, X, c, h=0.0):
        p = dict(self.p)
        p[self.control], p[self.bias] = float(c), float(h)
        x, a = np.broadcast_arrays(np.asarray(x, float), np.asarray(a, float))
        Xc = np.broadcast_to(np.asarray(X, float), x.shape)
        cols = np.stack([x.ravel(), a.ravel(), Xc.ravel()], axis=1)
        return self._f(cols, p)[:, 0].reshape(x.shape)

    def quantiles(self, M: int) -> np.ndarray:
        return self.centre + self.width * norm.ppf((np.arange(M) + 0.5) / M)

    def profile(self, a: np.ndarray, X: float, c: float, h: float = 0.0) -> np.ndarray:
        """x*(a): the stationary unit in the field X, by elementwise Newton steps (bounded to 1) from the sign of a, and
        bracketing on [-span, span] for any unit that has not converged."""
        a = np.asarray(a, float)
        x = 0.5 * np.sign(a - self.centre)
        e = 1e-7
        for _ in range(200):
            fx = self.f(x, a, X, c, h)
            d = (self.f(x + e, a, X, c, h) - self.f(x - e, a, X, c, h)) / (2 * e)
            step = np.clip(np.where(d != 0, fx / d, 0.0), -1.0, 1.0)
            x = x - step
            if np.max(np.abs(step)) < 1e-14:
                break
        bad = np.nonzero(np.abs(self.f(x, a, X, c, h)) > 1e-10)[0]
        span = 10.0 * (1.0 + abs(self.centre) + self.width)
        for i in bad:
            x[i] = brentq(lambda u: float(self.f(np.array([u]), np.array([a[i]]), X, c, h)[0]), -span, span,
                          xtol=1e-15)
        return x

    def collective_state(self, a: np.ndarray, c: float, h: float = 0.0) -> np.ndarray:
        if h == 0:
            return self.profile(a, 0.0, c)
        X = brentq(lambda X: self.profile(a, X, c, h).mean() - X, -10, 10, xtol=1e-14)
        return self.profile(a, X, c, h)

    def jac(self, x: np.ndarray, a: np.ndarray, c: float, eps: float = 1e-7) -> np.ndarray:
        """J = diag(df/dx) + (1/M) df/dX 1^T for the M representatives."""
        X = x.mean()
        dfx = (self.f(x + eps, a, X, c) - self.f(x - eps, a, X, c)) / (2 * eps)
        dfX = (self.f(x, a, X + eps, c) - self.f(x, a, X - eps, c)) / (2 * eps)
        return np.diag(dfx) + np.outer(dfX, np.ones_like(x)) / len(x)

    # ------------------------------------------------------------------ reduction
    def reduce(self) -> Dict:
        a = self.quantiles(self.M)
        lam = lambda c: float(np.linalg.eigvals(self.jac(self.profile(a, 0.0, c), a, c)).real.max())   # noqa: E731
        try:
            cs = brentq(lam, *self.range, xtol=1e-12)
        except ValueError as error:
            raise SpecError("the symmetric state does not lose stability inside population.control.range") from error
        x = self.profile(a, 0.0, cs)
        J = self.jac(x, a, cs)
        ev, V = np.linalg.eig(J)
        v = np.real(V[:, int(np.argmin(np.abs(ev)))])
        evl, W = np.linalg.eig(J.T)
        w = np.real(W[:, int(np.argmin(np.abs(evl)))])
        w /= w @ v
        dc = 1e-6 * max(1.0, abs(cs))
        da = (lam(cs + dc) - lam(cs - dc)) / (2 * dc)
        eps = 1e-6
        dfh = (self.f(x, a, x.mean(), cs, eps) - self.f(x, a, x.mean(), cs, -eps)) / (2 * eps)
        dfX = (self.f(x, a, x.mean() + eps, cs) - self.f(x, a, x.mean() - eps, cs)) / (2 * eps)
        hs = float(w @ dfh)
        if hs < 0:
            v, w, hs = -v, -w, -hs
        # frozen bias of a random sample, in units of h: (sum w df/dX / sum w df/dh) (Var x* / N)^(1/2)
        grid = self.centre + self.width * np.linspace(-9, 9, 4001)
        xs = self.profile(grid, 0.0, cs)
        trap = getattr(np, "trapezoid", None) or np.trapz
        dens = norm.pdf(grid, self.centre, self.width)
        mean_x = float(trap(dens * xs, grid))
        var_x = float(trap(dens * (xs - mean_x) ** 2, grid))
        ratio = float(abs((w @ dfX) / (w @ dfh)))
        return {"c_star": float(cs), "a": float(da), "h_s": hs, "v": v, "w": w, "x_sym": x,
                "frozen_per_sqrtN": ratio * sqrt(var_x), "var_x": var_x}

    def D_s(self, red: Dict, N: float) -> float:
        """w^T (D M/N) w for the M representatives (each stands for N/M units)."""
        return float(self.noise * self.M / N * (red["w"] @ red["w"]))

    def law(self, red: Dict, h: float, N: float, r: float, sample: str = "quantiles") -> Dict[str, float]:
        Ds = self.D_s(red, N)
        sig_th = sqrt(Ds) * (red["a"] * r) ** 0.25 / pi ** 0.25 / red["h_s"]
        s_q = red["frozen_per_sqrtN"] / sqrt(N) if sample == "random" else 0.0
        z = h / sqrt(sig_th ** 2 + s_q ** 2)
        return {"P": float(norm.cdf(z)), "probit": float(z), "sigma_thermal": sig_th, "s_frozen": s_q, "D_s": Ds}

    # ------------------------------------------------------------------ simulation
    def sweep(self, red: Dict, h: float, N: int, r: float, c0: float, c1: float, replicas: int,
              sample: str = "quantiles", seed: int = 1, dt: float = 1e-2, burn: float = 50.0) -> Dict[str, float]:
        """Euler-Maruyama for replicas x N units, relaxed at c0 for `burn`, swept to c1; the fraction ending with the
        sign of the bias (the sign of X at c1)."""
        rng = np.random.default_rng(seed)
        A = (self.quantiles(N)[None, :].repeat(replicas, 0) if sample == "quantiles"
             else self.centre + self.width * rng.standard_normal((replicas, N)))
        x = np.sign(A - self.centre) * np.abs(A - self.centre) ** (1 / 3)
        amp = sqrt(2 * self.noise * dt)
        c = c0
        n_burn = int(burn / dt)
        steps = int(np.ceil((c1 - c0) / (r * dt)))
        for k in range(n_burn + steps):
            X = x.mean(axis=1, keepdims=True)
            x = x + dt * self.f(x, A, X, c, h) + amp * rng.standard_normal(x.shape)
            if k >= n_burn:
                c += r * dt
        s = np.sign(x.mean(axis=1))
        P = float(np.mean(s == np.sign(h)) + 0.5 * np.mean(s == 0))
        return {"P": P, "stderr": sqrt(max(P * (1 - P), 1e-12) / replicas), "replicas": replicas, "N": N,
                "sample": sample}
