"""Limit-cycle units with a spread parameter, coupled all to all through an observable:

    dx_i/dt = F(x_i; p_i) + K (c_bar - c(x_i)) e_k,   c_bar = (1/N) sum_j c(x_j).

Reduction (letters S, C, R, H, K): the cycle, its Floquet rate kappa and phase response Z of the unit at the centre of
the parameter density; the coupling function H(phi) = (1/2 pi) int Z_k(psi) [c(x(psi + phi)) - c(x(psi))] d psi and its
first harmonics a1, b1; the frequency map omega(p) and the density g(omega) = rho(p(omega)) / |d omega/dp| (a change of
variables, not a kernel estimate); K_c and Omega from the dispersion relation (onset.py) and the growth rate mu of
synchrony at K = factor K_c.

Check (letter L): full units started on the cycles of the unit forced by the self-consistent incoherent mean field,
each parameter group spread evenly in time (symmetric rings), with a coherent seed 1e-7; the growth rate of the
amplitude of the mean field at the reference frequency, period by period, fitted after the amplitude has passed its
minimum and grown a hundredfold above it. The ratio to mu tends to 1 as the spread (and K/kappa) goes to zero; the
excess at finite spread is the amplitude correction of order K/kappa.
"""
from __future__ import annotations

import copy
from math import pi
from typing import Callable, Dict

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.stats import norm

from ..memory import phase_locking as pl
from ..memory import spec as memory_spec
from . import onset

SpecError = memory_spec.SpecError


def observable(real, expr: str) -> Callable[[np.ndarray, Dict], np.ndarray]:
    """c(x) compiled from the restricted parser in the unit's variables and parameters."""
    sp = memory_spec._sympy()
    sym = real.symbolic
    names = {str(s): s for s in sym["variables"] + sym["parameters"]}
    e = memory_spec.parse_expression(expr, names)
    fn = sp.lambdify(sym["variables"] + sym["parameters"], e, modules="numpy")
    pn = [str(s) for s in sym["parameters"]]

    def c(X: np.ndarray, p: Dict) -> np.ndarray:
        X = np.atleast_2d(X)
        with np.errstate(all="ignore"):
            v = fn(*[X[:, i] for i in range(X.shape[1])], *[p[k] for k in pn])
        return np.broadcast_to(np.asarray(v, float), X[:, 0].shape)

    return c


class Oscillators:
    kind = "oscillators"

    def __init__(self, spec: Dict):
        pop = memory_spec._require(spec, "population", dict)
        het = memory_spec._require(pop, "heterogeneity", dict)
        self.param = memory_spec._require(het, "parameter", str)
        if het.get("density", "gaussian") != "gaussian":
            raise SpecError("oscillators: only a gaussian heterogeneity is supported")
        self.centre = memory_spec._number(memory_spec._require(het, "centre"), "heterogeneity.centre")
        self.spread = memory_spec._number(memory_spec._require(het, "spread"), "heterogeneity.spread")
        if not 0 < self.spread < 0.2:
            raise SpecError("heterogeneity.spread is the frequency spread sigma_omega/omega, in (0, 0.2)")
        cpl = memory_spec._require(pop, "coupling", dict)
        self.acting_on = memory_spec._require(cpl, "acting_on", str)
        self.through = memory_spec._require(cpl, "through", str)
        self.factors = [memory_spec._number(f, "coupling.factors") for f in cpl.get("factors", [1.5])]
        self.unit_spec = {k: v for k, v in spec.items() if k not in ("schema", "kind", "population")}
        self.unit_spec.update(schema=memory_spec.SCHEMA, kind="equations", noise=0.0)
        self.real = memory_spec.load(self.unit_spec)
        if self.param not in self.real.params:
            raise SpecError("heterogeneity.parameter must be a declared parameter")
        if self.acting_on not in self.real.variables:
            raise SpecError("coupling.acting_on must be a variable of the unit")
        self.k = self.real.variables.index(self.acting_on)
        self.c = observable(self.real, self.through)

    # ------------------------------------------------------------------ reduction
    def cycle(self, real, value: float, seed: int = 0) -> Dict:
        cyc = pl.find_cycle(real, {self.param: float(value)}, np.random.default_rng(seed))
        if cyc is None:
            raise SpecError(f"no limit cycle at {self.param} = {value}")
        cyc = pl.refine_cycle(real, {self.param: float(value)}, cyc)
        cyc["kappa"] = pl.floquet(real, {self.param: float(value)}, cyc)["kappa"]
        return cyc

    def coupling_function(self, n_grid: int = 512) -> Dict:
        cyc = self.cycle(self.real, self.centre)
        over = {self.param: self.centre}
        pr = pl.phase_response(self.real, over, cyc, self.param, n_grid=n_grid)
        T = cyc["period"]
        X = cyc["orbit"](T * np.arange(n_grid) / n_grid).T
        cx = self.c(X, {**self.real.params, **over})
        Zk = pr["Z"][:, self.k]
        H = np.array([np.mean(Zk * (np.roll(cx, -n) - cx)) for n in range(n_grid)])
        phi = 2 * pi * np.arange(n_grid) / n_grid
        return {"omega": 2 * pi / T, "period": T, "kappa": float(cyc["kappa"]), "H": H, "phi": phi,
                "a1": float(2 * np.mean(H * np.cos(phi))), "b1": float(2 * np.mean(H * np.sin(phi)))}

    def frequency(self, value: float) -> float:
        return float(2 * pi / self.cycle(self.real, value)["period"])

    def reduce(self, grid_points: int = 25) -> Dict:
        cf = self.coupling_function()
        dp = 1e-3 * max(1.0, abs(self.centre))
        slope = (self.frequency(self.centre + dp) - self.frequency(self.centre - dp)) / (2 * dp)
        if slope == 0:
            raise SpecError("the frequency does not depend on the spread parameter")
        s = self.spread * cf["omega"] / abs(slope)
        grid = self.centre + s * np.linspace(-6, 6, grid_points)
        w = np.array([self.frequency(v) for v in grid])
        order = np.argsort(w)
        if np.any(np.diff(w[order]) <= 0):
            raise SpecError("omega(p) is not monotone over the spread")
        inv = CubicSpline(w[order], grid[order])
        dinv = inv.derivative()
        lo, hi = float(w.min()), float(w.max())

        # g tabulated on a fine grid (linear interpolation; relative error about 1e-6) for the quadratures
        xs = np.linspace(lo, hi, 8001)
        gs = norm.pdf(inv(xs), self.centre, s) * np.abs(dinv(xs))
        gs[0] = gs[-1] = 0.0

        def g(x):
            return float(np.interp(x, xs, gs, left=0.0, right=0.0))

        scale = self.spread * cf["omega"]
        Kc = onset.critical_coupling(g, cf["a1"], cf["b1"], centre=cf["omega"], scale=scale, support=(lo, hi),
                                     grid=61)
        keep = slice(None, None, max(1, len(cf["H"]) // 128))
        out = {**{k: cf[k] for k in ("omega", "period", "kappa", "a1", "b1")}, "domega_dp": float(slope),
               "H_phi": [float(x) for x in cf["phi"][keep]], "H": [float(x) for x in cf["H"][keep]],
               "g_omega": [float(x) for x in xs[::40]], "g": [float(x) for x in gs[::40]],
               "param_spread": float(s), "density": g, "support": (lo, hi), "scale": scale,
               "K_c": None if Kc is None else Kc["K_c"], "Omega_c": None if Kc is None else Kc["Omega"], "rows": []}
        if Kc is not None:
            for f in self.factors:
                K = f * Kc["K_c"]
                gr = onset.growth_rate(g, cf["a1"], cf["b1"], K, centre=cf["omega"], scale=scale, support=(lo, hi))
                out["rows"].append({"factor": f, "K": K, "mu": gr["mu"], "Omega": gr["Omega"],
                                    "K_over_kappa": K / cf["kappa"]})
        return out

    # ------------------------------------------------------------------ full units
    def _forced(self, K: float, cbar: float):
        s = copy.deepcopy(self.unit_spec)
        s["drift"][self.acting_on] = f"({s['drift'][self.acting_on]}) + Kpop*(cpop - ({self.through}))"
        s["parameters"] = {**s["parameters"], "Kpop": float(K), "cpop": float(cbar)}
        return memory_spec.load(s)

    def simulate(self, red: Dict, factor: float, n_p: int = 16, n_theta: int = 8, steps: int = 400,
                 iterations: int = 3) -> Dict:
        """The measured growth rate of synchrony of n_p x n_theta full units against the reduction's mu."""
        row = next(r for r in red["rows"] if r["factor"] == factor)
        K = row["K"]
        s = red["param_spread"]
        groups = self.centre + s * norm.ppf((np.arange(n_p) + 0.5) / n_p)
        vals = np.repeat(groups, n_theta)
        P = {**self.real.params, self.param: vals}
        ref = self.cycle(self.real, self.centre)
        T = ref["period"]
        orbit = ref["orbit"](T * np.arange(512) / 512).T
        cbar = float(np.mean(self.c(orbit, {**self.real.params, self.param: self.centre})))
        base = 2 * pi * np.arange(n_theta) / n_theta
        th = base + 2e-7 * np.sin(base)
        X = np.empty((len(vals), len(self.real.variables)))
        for _ in range(iterations):
            unit = self._forced(K, cbar)
            means = []
            for gi, v in enumerate(groups):
                cyc = self.cycle(unit, v)
                Tg = cyc["period"]
                X[gi * n_theta:(gi + 1) * n_theta] = cyc["orbit"](Tg * th / (2 * pi)).T
                grid = cyc["orbit"](Tg * np.arange(256) / 256).T
                means.append(float(np.mean(self.c(grid, {**self.real.params, self.param: float(v)}))))
            cbar = float(np.mean(means))

        def drift(X):
            F = self.real.drift(X, P)
            cx = self.c(X, P)
            F[:, self.k] += K * (cx.mean() - cx)
            return F
        dt = T / steps
        n_steps = int(round((18.0 / row["mu"] + 3 * T) / dt))
        cb = np.empty(n_steps + 1)
        for k in range(n_steps + 1):
            cb[k] = float(np.mean(self.c(X, P)))
            if k == n_steps:
                break
            k1 = drift(X)
            k2 = drift(X + 0.5 * dt * k1)
            k3 = drift(X + 0.5 * dt * k2)
            k4 = drift(X + dt * k3)
            X = X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        span = float(np.ptp(self.c(orbit, {**self.real.params, self.param: self.centre})))
        ph = np.exp(-2j * pi * np.arange(steps) / steps)
        nb = len(cb) // steps
        amp = np.array([abs(2.0 / steps * np.sum(cb[i * steps:(i + 1) * steps] * ph)) for i in range(nb)]) / span
        tc = (np.arange(nb) + 0.5) * steps * dt
        above = np.nonzero(amp > 5e-2)[0]
        end = above[0] if len(above) else nb
        imin = int(np.argmin(amp[:end])) if end > 0 else 0
        lo = max(1e-6, 100.0 * float(amp[imin])) if end > 0 else 1e-6
        sel = np.nonzero((np.arange(nb) >= imin) & (np.arange(nb) < end) & (amp > lo) & (amp < 5e-2))[0]
        mu = float(np.polyfit(tc[sel], np.log(amp[sel]), 1)[0]) if len(sel) >= 5 else float("nan")
        return {"factor": factor, "K": K, "mu_reduction": row["mu"], "mu_measured": mu,
                "ratio": mu / row["mu"] if np.isfinite(mu) else float("nan"), "points": int(len(sel)),
                "units": int(len(vals)), "K_over_kappa": row["K_over_kappa"]}
