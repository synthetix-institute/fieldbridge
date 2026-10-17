"""The write law of a symmetric switch across the noise ratio.

A switch swept slowly through a pitchfork chooses the state favoured by a weak bias h with probability

    P = Phi(c h / (D^{1/2} r^{1/4})),

with D the noise and r the rate at which the growth rate of the critical mode rises. In the scaled variables
tau = r^{1/2} t, xi = x |a3|^{1/2} / r^{1/4} and eta = h |a3|^{1/2} / r^{3/4} the normal form reads

    d xi = (tau xi -+ xi^3 + eta) d tau + (2 Gamma)^{1/2} dW,    Gamma = D |a3| / r,

so that the constant c is a function of the noise ratio Gamma alone. Gamma compares the noise with the barrier that
grows behind the choice. For Gamma << 1 the decision is made in the linear stage at the bifurcation and c = pi^{1/4}
(Kondepudi and Nelson). For Gamma >> 1 the state hops between the wells after the bifurcation until the growing barrier
freezes it; the outcome is a Boltzmann ratio and c approaches the activated limit
sqrt(pi) [ln(1.6 Gamma)]^{1/4} Gamma^{-1/4} (Kramers rates at freeze-out). A first-order switch (positive cubic term)
has the same function within 4%: it escapes over the shrinking barrier before the bifurcation, as far before it as the
continuous switch freezes after it. A first-order switch follows the law only while the choice is made inside its
bistable range (the range of the control in which its far states exist); beyond it the far states decide.

The tabulated values are simulations of the normal form with the bias toward each side in turn (40 000 sweeps per side
for Gamma <= 3 and 20 000 above; standard errors 0.012-0.016); for Gamma <= 1 they agree with pi^{1/4} within their
errors, and the table uses pi^{1/4} there.
"""
from __future__ import annotations

from math import exp, log, pi, sqrt
from statistics import NormalDist
from typing import Dict, Optional, Sequence

import numpy as np

from . import analysis as an
from .identity import Realization

KONDEPUDI = pi ** 0.25
GAMMAS = (1.0, 3.0, 10.0, 30.0, 100.0)
CONTINUOUS = (KONDEPUDI, 1.268, 1.144, 0.949, 0.784)
FIRST_ORDER = (KONDEPUDI, 1.250, 1.110, 0.962, 0.791)
# share of the sweeps of a first-order switch decided before the bifurcation, on the normal form
EARLY_GAMMAS = (0.01, 0.03, 0.1, 0.3, 1.0, 3.0)
EARLY = (0.0, 0.0, 0.031, 0.314, 0.829, 0.997)
LINEAR, CROSSOVER, ACTIVATED = "linear stage", "crossover", "activated"


def activated_constant(gamma: float) -> float:
    """Kramers freeze-out: the difference of the occupations relaxes at 2k = (2^{1/2} tau / pi) e^{-tau^2 / 4 Gamma}
    toward eta tau^{1/2} / Gamma, which gives c = pi^{1/2} Gamma^{-1/4} <[ln(A/K)]^{1/4}> over K ~ Exp(1) on [0, A],
    A = 2^{3/2} Gamma / pi; approximately pi^{1/2} [ln(1.6 Gamma)]^{1/4} Gamma^{-1/4}. An asymptotic estimate for
    Gamma >> 1 (about 10% above the simulated function at Gamma = 30-100)."""
    from scipy.integrate import quad
    A = 2.0 * sqrt(2.0) * gamma / pi
    if A <= 1.0:
        return float("nan")
    val, _ = quad(lambda K: max(log(A / K), 0.0) ** 0.25 * exp(-K), 0.0, A, limit=200)
    return sqrt(pi) * gamma ** -0.25 * val


def write_constant(gamma: float, first_order: bool = False) -> float:
    """c(Gamma): pi^{1/4} up to Gamma = 1, interpolated in ln Gamma through the simulated table up to Gamma = 100, and
    continued above it by the activated limit scaled to the last tabulated value."""
    if gamma <= GAMMAS[0]:
        return KONDEPUDI
    table = FIRST_ORDER if first_order else CONTINUOUS
    if gamma <= GAMMAS[-1]:
        return float(np.interp(log(gamma), [log(g) for g in GAMMAS], table))
    return table[-1] * activated_constant(gamma) / activated_constant(GAMMAS[-1])


def early_share(gamma: float) -> float:
    """Share of the sweeps of a first-order switch decided before its bifurcation (normal form)."""
    if gamma >= EARLY_GAMMAS[-1]:
        return 1.0
    if gamma <= EARLY_GAMMAS[0]:
        return 0.0
    return float(np.interp(log(gamma), [log(g) for g in EARLY_GAMMAS], EARLY))


def regime(gamma: float) -> str:
    """linear stage (c = pi^{1/4} within the simulated errors) up to Gamma = 1, activated (c below 0.75 pi^{1/4})
    from Gamma = 30, crossover between."""
    return LINEAR if gamma <= 1.0 else ACTIVATED if gamma >= 30.0 else CROSSOVER


def probability(z_kondepudi: float, gamma: float, first_order: bool = False) -> float:
    """Accuracy of a write whose Kondepudi-Nelson argument is z = pi^{1/4} h / (D^{1/2} r^{1/4}), at noise ratio
    Gamma."""
    return NormalDist().cdf(z_kondepudi * write_constant(gamma, first_order) / KONDEPUDI)


# ------------------------------------------------------------------------------------------------ the normal form
def normal_form_sweeps(gamma: float, eta: float, n: int, rng, first_order: bool = False) -> Dict[str, float]:
    """Sweep the scaled normal form through the bifurcation; returns the share of sweeps that end on the side of a
    positive bias eta, and for a first-order switch the share decided before the bifurcation (escape past
    1.5 tau_max^{1/2} + 3, beyond the unstable branch at every tau)."""
    s = 1.0 if first_order else -1.0
    span = max(5.0, 7.0 * sqrt(gamma))
    L = 1.5 * sqrt(span) + 3.0
    dt = min(0.004, 0.25 / span)     # resolves the curvature of the wells; an overshoot past L keeps its sign
    amp = sqrt(2.0 * gamma * dt)
    xi = np.zeros(n)
    alive = np.ones(n, bool)
    out = np.zeros(n)
    when = np.full(n, np.nan)
    tau = -span
    while tau < span:
        xi = np.where(alive, xi + dt * (tau * xi + s * xi ** 3 + eta) + amp * rng.standard_normal(n), xi)
        if first_order:
            gone = alive & (np.abs(xi) > L)
            out[gone] = np.sign(xi[gone])
            when[gone] = tau
            alive &= ~gone
            if not alive.any():
                break
        tau += dt
    out[alive] = np.sign(xi[alive])
    return {"toward_bias": float(np.mean(out > 0)),
            "decided_before": float(np.mean(when < 0)) if first_order else 0.0}


def normal_form_constant(gamma: float, n: int, rng, first_order: bool = False) -> Dict[str, float]:
    """c(Gamma) measured on the normal form with the bias toward each side in turn (n sweeps per side)."""
    eta = 0.5 * sqrt(gamma) / write_constant(gamma, first_order)        # accuracies near 0.69
    plus = normal_form_sweeps(gamma, eta, n, rng, first_order)
    minus = normal_form_sweeps(gamma, -eta, n, rng, first_order)
    P = 0.5 * (plus["toward_bias"] + (1.0 - minus["toward_bias"]))
    return {**_constant(P, 2 * n, sqrt(gamma) / eta), "gamma": gamma, "first_order": first_order,
            "decided_before": 0.5 * (plus["decided_before"] + minus["decided_before"])}


def _constant(P: float, n: int, scale: float) -> Dict[str, float]:
    se = sqrt(P * (1.0 - P) / n)
    zP = NormalDist().inv_cdf(min(max(P, 1e-9), 1 - 1e-9))
    return {"measured": P, "stderr": se, "constant": zP * scale,
            "constant_stderr": se / NormalDist().pdf(zP) * scale}


# ------------------------------------------------------------------------------------------------ a realization
def symmetric_point(real: Realization, rng=None) -> Optional[Dict]:
    """The first critical point along the control at which one direction becomes critical and a symmetry removes the
    even terms: a continuous (supercritical) or first-order (subcritical) pitchfork. None if there is none."""
    rng = rng if rng is not None else np.random.default_rng(1)
    for ev in an.locate_writes(real, rng, n_starts=24 if real.carrier.dim <= 4 else 12):
        nf = an.normal_form(real, ev["q"], real.control, ev["v"])
        kind = nf.get("kind", "")
        if (kind.startswith("supercritical pitchfork") or kind.startswith("subcritical pitchfork")) \
                and nf.get("critical_dimension", 1) == 1 and not nf.get("hopf") and abs(nf.get("a3", 0)) > 0:
            return nf
    return None


def material_noise_ratio(real: Realization, nf: Dict, rate: float) -> float:
    """Gamma = D_s |a3| / r of the material's own noise for a sweep at `rate` (growth rate per unit time)."""
    w = np.asarray(nf["left_mode"], float)
    return float(real.noise * float(w @ w) * abs(nf["a3"]) / rate)


def sweeps(real: Realization, nf: Dict, rng, gammas: Sequence[float] = (0.05, 1.0, 10.0), q: float = 5e-4,
           n_traj: int = 400) -> Dict[str, object]:
    """Sweep the full model through its symmetric critical point at several noise ratios and compare the constant of
    the write with c(Gamma). The product D_s |a3| = q is fixed, so that the choice is made near the bifurcation
    (within |eps| of order (q ln Gamma)^{1/2}), and the rate is r = q / Gamma. The bias is a constant force along the
    critical mode, toward each side in turn with n_traj sweeps each; the constant is computed from the mean of the two
    one-sided accuracies, and half their difference is the offset of the sweep. For a first-order switch the share of
    sweeps that pass half the distance to the far state before the bifurcation is recorded. The cost grows as
    Gamma / q^{1/2}: Gamma = 10 takes about 10^5-10^6 steps."""
    first_order = abs(nf["a3"]) > 0 and nf["a3"] > 0
    vc = float(nf["value"])
    qc, vec = np.asarray(nf["state"], float), np.asarray(nf["mode"], float)
    w = np.asarray(nf["left_mode"], float)
    left = w / float(w @ vec)
    g, a3 = abs(nf["kappa_slope"]), abs(nf["a3"])
    D_s = q / a3
    D = D_s / float(w @ w)
    dt = min(an.safe_dt(real), 0.02)
    sign_control = np.sign(-nf["kappa_slope"])
    far = 0.5 * nf["half"] if real.carrier.kind == "torus" else 1e3 * an.base_length(real)
    rows = []
    for G in gammas:
        r = q / G
        E = max(5.0 * sqrt(r), 5.2 * sqrt(q * max(log(1.6 * G + 1.0), 0.5)))
        dv = sign_control * r / g
        c0 = vc - dv * E / r
        steps = int(round(2 * E / r / dt))
        h_s = 0.5 * sqrt(D_s) * r ** 0.25 / write_constant(G, first_order)
        s_far = np.inf
        if first_order:
            # the far state, reached without noise from a small displacement at the end of the sweep; a sweep is
            # decided when its amplitude passes half the distance to it
            qd = qc + 1e-2 * vec * max(1.0, float(np.linalg.norm(qc)))
            for _ in range(int(round(400.0 / dt))):
                qd = real.carrier.wrap(qd + dt * np.asarray(real.F(qd[None], **{real.control: c0 + dv * steps * dt}),
                                                            float)[0])
            s_far = 0.5 * abs(float(real.carrier.diff(qc[None], qd[None])[0] @ left))
        by_sign = {}
        for sg in (1, -1):
            target = qc + sg * far * vec
            if real.carrier.kind == "torus":
                target = real.carrier.wrap(target)
            h_per = abs(float(w @ real.write_force(qc, target, 1.0)))
            x = np.repeat(qc[None], n_traj, axis=0)
            first = np.full(n_traj, np.nan)
            amp = sqrt(2.0 * D * dt)
            for k in range(steps):
                c = c0 + dv * k * dt
                f = np.asarray(real.F(x, **{real.control: c}), float) + real.write_force(x, target, h_s / h_per)
                x = real.carrier.wrap(x + dt * f + amp * rng.standard_normal(x.shape))
                if first_order:
                    s = real.carrier.diff(np.repeat(qc[None], n_traj, axis=0), x) @ left
                    new = np.isnan(first) & (np.abs(s) > s_far)
                    first[new] = (k + 1) * dt * r - E
            s = real.carrier.diff(np.repeat(qc[None], n_traj, axis=0), x) @ left
            by_sign[sg] = {"toward_bias": float(np.mean(np.sign(s) == sg)),
                           "decided_before": float(np.mean(first < 0)) if first_order else 0.0}
        P = 0.5 * (by_sign[1]["toward_bias"] + by_sign[-1]["toward_bias"])
        row = _constant(P, 2 * n_traj, sqrt(D_s) * r ** 0.25 / h_s)
        row.update(gamma=float(G), rate=float(r), D_along_mode=float(D_s), steps=steps,
                   offset=0.5 * (by_sign[1]["toward_bias"] - by_sign[-1]["toward_bias"]),
                   by_sign={str(k): v for k, v in by_sign.items()},
                   expected=write_constant(G, first_order), regime=regime(G),
                   decided_before=0.5 * (by_sign[1]["decided_before"] + by_sign[-1]["decided_before"]))
        if first_order:
            row["decided_before_normal_form"] = early_share(G)
        rows.append(row)
    return {"name": real.name, "first_order": first_order, "control": real.control, "critical_value": vc,
            "a3": float(nf["a3"]), "q": q, "trajectories_per_sign": n_traj, "rows": rows,
            "material": {"noise": float(real.noise),
                         "noise_ratio_at_rate_0.02": material_noise_ratio(real, nf, 0.02),
                         "crossover_rate": material_noise_ratio(real, nf, 1.0) / 10.0}}


def markdown(rep: Dict[str, object]) -> str:
    kind = "first-order (subcritical)" if rep["first_order"] else "continuous (supercritical)"
    lines = [f"# Write law across the noise ratio: {rep['name']}", "",
             f"A {kind} symmetric switch at {rep['control']} = {rep['critical_value']:.6g} (cubic coefficient "
             f"{rep['a3']:+.4g}). D_s|a3| = {rep['q']:g}; {rep['trajectories_per_sign']} sweeps per side of the bias.",
             "", "| Gamma | regime | accuracy | constant c | expected c(Gamma) | offset |"
             + (" decided before |" if rep["first_order"] else ""),
             "|---|---|---|---|---|---|" + ("---|" if rep["first_order"] else "")]
    for r in rep["rows"]:
        line = (f"| {r['gamma']:g} | {r['regime']} | {r['measured']:.3f} | {r['constant']:.3f} +- "
                f"{r['constant_stderr']:.3f} | {r['expected']:.3f} | {r['offset']:+.3f} |")
        if rep["first_order"]:
            line += f" {r['decided_before']:.2f} (normal form {r['decided_before_normal_form']:.2f}) |"
        lines.append(line)
    m = rep["material"]
    lines += ["", f"At its own noise ({m['noise']:g}) a sweep at rate 0.02 has Gamma = "
                  f"{m['noise_ratio_at_rate_0.02']:.3g} ({regime(m['noise_ratio_at_rate_0.02'])}); sweeps slower than "
                  f"{m['crossover_rate']:.3g} (Gamma > 10) lose part of the fourth-root gain."]
    return "\n".join(lines) + "\n"


def law_markdown() -> str:
    lines = ["# The write law across the noise ratio", "",
             "P = Phi(c h / (D^{1/2} r^{1/4})), Gamma = D|a3|/r; c = pi^{1/4} = 1.3313 for Gamma << 1.", "",
             "| Gamma | c, continuous | c, first order | activated limit | regime |", "|---|---|---|---|---|"]
    for G in (0.1, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0):
        act = activated_constant(G) if G >= 10 else float("nan")
        lines.append(f"| {G:g} | {write_constant(G):.3f} | {write_constant(G, True):.3f} | "
                     f"{'' if act != act else f'{act:.3f}'} | {regime(G)} |")
    return "\n".join(lines) + "\n"
