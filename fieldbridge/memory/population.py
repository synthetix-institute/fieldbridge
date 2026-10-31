"""The time to leave a state of a population whose noise comes from its own births and deaths.

A specification may declare the rates behind its one-variable drift, F(x) = b(x) - d(x), per unit of the population
scale N. The population then moves one individual at a time, and the mean time to leave a state grows as exp(N S). Three
values of S are reported for each way out of a state:

  rates          S = int ln(d/b) dx against the drift: the exponent of the birth-death chain
  state noise    S = int 2 (d - b)/(d + b) dx: a diffusion with the intensity (b + d)/2N of each state
  constant noise S = (barrier of F)/(N D), with D = b/N at the state being left

The second agrees with the first only when births and deaths nearly balance along the path (Doering, Sargsyan and
Sander, 2005), and the third can err in either direction. The mean time of the chain itself is computed exactly at a
given N, with a reflecting state at n = 1 for an upward path and at the largest n for a downward one. Where two states
coexist along the control, the value at which they are equally occupied is found with the rates and with a constant
noise; for a constant noise it does not depend on the intensity.
"""
from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from .identity import Realization


def _functions(real: Realization, over: Optional[Dict[str, float]]):
    rates = getattr(real, "rates", None)
    if rates is None:
        raise ValueError(f"{real.name} declares no birth and death rates (the field 'rates' of a specification)")
    p = {**real.params, **(over or {})}
    b = lambda x: rates["birth"](x, p)
    d = lambda x: rates["death"](x, p)
    return b, d, float(rates["upper"](p))


def _roots(b, d, upper: float) -> List[Dict[str, object]]:
    """The states of the drift in (0, upper], with the state x = 0 left out; stable where F falls with x."""
    grid = np.unique(np.concatenate((np.geomspace(1e-9 * upper, upper, 3000), np.linspace(1e-9 * upper, upper, 3000))))
    F = b(grid) - d(grid)
    found = []
    for i in np.nonzero(F[:-1] * F[1:] < 0)[0]:
        found.append({"x": float(brentq(lambda x: float(b(x) - d(x)), grid[i], grid[i + 1], xtol=1e-12)),
                      "stable": bool(F[i] > 0)})
    return found


def _actions(b, d, start: float, stop: float) -> Dict[str, float]:
    # the integrals have a logarithmic end point at x = 0, which quad handles
    def integral(f):
        return float(quad(f, start, stop, limit=400)[0])
    S = integral(lambda x: np.log(d(x) / b(x)))
    grid = np.linspace(start, stop, 4001)[1:-1]
    return {"rates": S, "state_noise": integral(lambda x: 2 * (d(x) - b(x)) / (d(x) + b(x))),
            "constant_noise": integral(lambda x: d(x) - b(x)) / float(b(start)),
            "largest_log_ratio": float(np.max(np.abs(np.log(d(grid) / b(grid)))))}


def log_exit_time(b, d, upper: float, N: int, start: float, stop: float) -> float:
    """Natural logarithm of the exact mean time of the chain from the state nearest start to the state nearest stop."""
    n0, n1 = int(round(start * N)), int(round(stop * N))
    top = int(np.ceil(upper * N))
    n = np.arange(1, top + 1)
    log_b, log_d = np.log(N * b(n / N)), np.log(N * d(n / N))
    step, total = -np.inf, -np.inf
    if n1 > n0:    # T(n -> n+1) = (1 + d_n T(n-1 -> n)) / b_n, with no way down from n = 1
        for m in range(1, n1):
            step = np.logaddexp(0.0, log_d[m - 1] + step) - log_b[m - 1] if m > 1 else -log_b[0]
            if m >= n0:
                total = np.logaddexp(total, step)
    else:          # T(n -> n-1) = (1 + b_n T(n+1 -> n)) / d_n, with no way up from the largest n
        for m in range(top, n1, -1):
            step = np.logaddexp(0.0, log_b[m - 1] + step) - log_d[m - 1] if m < top else -log_d[top - 1]
            if m <= n0:
                total = np.logaddexp(total, step)
    return float(total)


def equal_occupation(real: Realization, over: Optional[Dict[str, float]] = None) -> Optional[Dict[str, float]]:
    """The value of the control at which two coexisting states are equally occupied, with the rates and with a constant
    noise, or None where the control has no such value."""
    if not real.control:
        return None
    lo, hi = real.control_range

    def difference(value):
        b, d, upper = _functions(real, {**(over or {}), real.control: value})
        stable = [r["x"] for r in _roots(b, d, upper) if r["stable"]]
        if len(stable) != 2:
            return None
        return (float(quad(lambda x: np.log(b(x) / d(x)), stable[0], stable[1], limit=400)[0]),
                float(quad(lambda x: b(x) - d(x), stable[0], stable[1], limit=400)[0]))

    grid = np.linspace(lo, hi, 171)
    values = [difference(v) for v in grid]
    out = {}
    for j, name in enumerate(("rates", "constant_noise")):
        for i in range(len(grid) - 1):
            if values[i] is None or values[i + 1] is None or values[i][j] * values[i + 1][j] >= 0:
                continue
            f = lambda v: difference(v)[j]
            out[name] = float(brentq(f, grid[i], grid[i + 1], xtol=1e-6))
            break
    return out if len(out) == 2 else None


def escape(real: Realization, N: Optional[int] = None, over: Optional[Dict[str, float]] = None) -> Dict[str, object]:
    """Every way out of every stable state, with its three exponents and, at a given N, the exact mean time."""
    b, d, upper = _functions(real, over)
    roots = _roots(b, d, upper)
    xs = [r["x"] for r in roots]
    paths = []
    for i, r in enumerate(roots):
        if not r["stable"]:
            continue
        below = xs[i - 1] if i > 0 else 0.0
        targets = [below] + ([xs[i + 1]] if i + 1 < len(xs) else [])
        for stop in targets:
            row = {"from": r["x"], "to": stop, "extinction": stop == 0.0, **_actions(b, d, r["x"], stop)}
            if N:
                row["N"] = int(N)
                row["ln_time_exact"] = log_exit_time(b, d, upper, int(N), r["x"], stop)
            paths.append(row)
    return {"states": roots, "paths": paths, "equal_occupation": equal_occupation(real, over),
            "scale": getattr(real, "rates")["scale"]}


def markdown(report: Dict[str, object]) -> str:
    rows = ["N is the " + report["scale"] + ".", "",
            "| Path | S from the rates | S with the noise of each state | S with constant noise |"
            + (" ln(exact time) |" if report["paths"] and "N" in report["paths"][0] else ""),
            "| --- | --- | --- | --- |" + (" --- |" if report["paths"] and "N" in report["paths"][0] else "")]
    for p in report["paths"]:
        end = "extinction" if p["extinction"] else f"{p['to']:.4g}"
        rows.append(f"| {p['from']:.4g} to {end} | {p['rates']:.4g} | {p['state_noise']:.4g} | {p['constant_noise']:.4g} |"
                    + (f" {p['ln_time_exact']:.4g} at N = {p['N']} |" if "N" in p else ""))
    rows += ["", "The mean time to take a path grows as exp(N S), with S from the rates. The noise of each state gives the "
             "same S only while births and deaths nearly balance along the path; constant noise is not quantitative."]
    eq = report["equal_occupation"]
    if eq:
        rows += ["", f"The two states are equally occupied at {eq['rates']:.4g} of the control with the rates, and at "
                 f"{eq['constant_noise']:.4g} with a constant noise of any intensity."]
    return "\n".join(rows)
