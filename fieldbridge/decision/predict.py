"""Predictions for a population, from its equations, before any stochastic run.

Classes:
    synchronizes        the incoherent state of the oscillators loses stability at a finite K_c (b1 > 0 or a shifted
                        root of the dispersion relation); the growth rate of synchrony above it is mu(K)
    no-onset            no root of the dispersion relation gives a positive K: no coupling synchronizes (b1 <= 0,
                        for example the neutral in-phase state of an overdamped Josephson array)
    follows-the-bias    a sweep through the collective pitchfork selects the favoured state with
                        P = Phi(pi^(1/4) h_s h / (D_s^(1/2) (a r)^(1/4))), D_s proportional to 1/N, or its version along
                        the actual passage
    set-by-the-sample   diverse units in a random sample: the frozen bias s_q of the sample exceeds the thermal spread
                        at every rate of the protocol, so P = Phi(h / (sigma_th^2 + s_q^2)^(1/2)) hardly depends on the
                        sweep rate
    reflection-seed     after a step one real eigenvalue leads (d = 1): the 10-90% window of passage times is
                        ln 13.09 in the growth exponent
    rotation-seed       a complex or degenerate pair leads (d = 2): the window is ln 4.675
"""
from __future__ import annotations

from math import log, pi, sqrt
from typing import Dict

from scipy.stats import norm

from . import passage

CLASSES = ("synchronizes", "no-onset", "follows-the-bias", "set-by-the-sample", "reflection-seed", "rotation-seed")


def _plain_red(red: Dict) -> Dict:
    return {k: v for k, v in red.items() if k not in ("density", "support", "v", "w", "y_sym", "x_sym")}


def predict(pop) -> Dict:
    b = pop.body
    if pop.target == "synchronization":
        red = b.reduce()
        out = {"reduction": _plain_red(red), "class": "synchronizes" if red["K_c"] is not None else "no-onset"}
        return {**out, "_red": red}
    red = b.reduce()
    if pop.target == "collective-write":
        rows, factors = [], {}
        for N in pop.sizes:
            for r in pop.rates:
                if pop.kind == "units":
                    base = b.law(red, 1.0, N, r, sample="quantiles")
                    for z in pop.z:
                        h = z * base["sigma_thermal"]
                        lw = b.law(red, h, N, r, sample=pop.sample)
                        rows.append({"N": N, "r": r, "z": z, "h": h, "P": lw["P"], "sigma_thermal": lw["sigma_thermal"],
                                     "s_frozen": lw["s_frozen"]})
                else:
                    Ds = b.D_s(red, N)
                    if r not in factors:          # (m, s2) at the first size; s2 scales as 1/N
                        factors[r] = (b.history_factors(red, pop.sizes[0], r, pop.sweep_from, pop.sweep_to, points=800),
                                      pop.sizes[0])
                    (m, s2), N0 = factors[r]
                    for z in pop.z:
                        h = z * sqrt(Ds) * (red["a"] * r) ** 0.25 / (pi ** 0.25 * red["h_s"])
                        lin = b.law(red, h, N, r)
                        zh = h * m / sqrt(s2 * N0 / N)
                        rows.append({"N": N, "r": r, "z": z, "h": h, "P_linear": lin["P"], "P": float(norm.cdf(zh)),
                                     "Lambda": lin["Lambda"], "window": lin["window"], "D_s": lin["D_s"]})
        cls = "follows-the-bias"
        if pop.kind == "units" and pop.sample == "random" and all(r["s_frozen"] > r["sigma_thermal"] for r in rows):
            cls = "set-by-the-sample"
        return {"reduction": _plain_red(red), "rows": rows, "class": cls, "_red": red}
    # seeded passage
    st = b.step_exponent(red, pop.step_from, pop.step_to, pop.duration)
    d = st["d"]
    out = {"reduction": _plain_red(red), "d": d, "rate_end": st["rate_end"], "window_law": log(passage.q_ratio(d)),
           "window_time_at_end": log(passage.q_ratio(d)) / st["rate_end"] if st["rate_end"] > 0 else None,
           "class": "reflection-seed" if d == 1 else "rotation-seed"}
    return {**out, "_red": red, "_step": st}
