"""The card of a regulated realization: what a step of the input does to the output, and why.

Classes:
  perfect-adaptation       the steady output does not depend on the input (G = 0) at every sampled parameter point
                           with a stable steady state; reported with the integrator that enforces it, or with the
                           note that none was found in the searched coordinates
  fine-tuned-adaptation    G = 0 at the stated parameters but not at the sampled ones
  partial-adaptation       after the step the output returns part of the way (|final| < |peak|); the remaining
                           fraction and the clamps that remove the return are reported
  no-adaptation            the output moves to its new value without returning (|final| = |peak|)
Absent:
  no set point             no stable steady state is reached at the reference input (unstable loop: a steady
                           state exists but is unstable)
  no response              the input does not reach the output
A vanishing gain is judged against the open-loop gain of the clamps (gains.reference_gain), not by the relative
sensitivity S = (u/y) G, which is undefined when the set point is zero; S is reported where it is defined.
"""
from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

from .gains import attenuation, input_reaches_output, reference_gain, static_gain
from .integrator import conservation_laws, find_integrator
from .spec import Regulated
from .step import step_response
from .steady import steady_state

GAIN_ZERO = 1e-8       # |G| / G_reference below this is a steady output independent of the input
RETURN_MIN = 0.02      # a return of less than 2 % of the peak deviation counts as no return


def _fd_gain(m: Regulated, p: np.ndarray, q: np.ndarray, L) -> Optional[float]:
    u0 = float(p[m.u_index])
    du = 1e-4 * max(abs(u0), 1e-3)
    ys = []
    for sgn in (-1, 1):
        pp = p.copy()
        pp[m.u_index] = u0 + sgn * du
        s = steady_state(m, pp, q, L)
        if not (s["converged"] and s["stable"]):
            return None
        ys.append(s["y"])
    return (ys[1] - ys[0]) / (2 * du)


def sample_parameters(m: Regulated, count: int, rng, factor: float = 2.0) -> List[np.ndarray]:
    """Parameter points with every parameter not held fixed multiplied by a factor in [1/factor, factor]."""
    free = [i for i, k in enumerate(m.pnames) if k not in m.fixed]
    base = m.pvec(m.u0)
    out = []
    for _ in range(count):
        p = base.copy()
        p[free] *= np.exp(rng.uniform(-np.log(factor), np.log(factor), size=len(free)))
        out.append(p)
    return out


def card(m: Regulated, samples: int = 32, seed: int = 0, steps: bool = True, first_step_only: bool = False,
         trace: bool = False) -> Dict[str, object]:
    rng = np.random.default_rng(seed)
    p = m.pvec(m.u0)
    out: Dict[str, object] = {"name": m.name, "field": m.field, "source": m.source, "spec": m.path,
                              "spec_sha256": m.spec_sha256, "input": m.input, "output": m.output_text,
                              "u0": m.u0, "steps": m.steps, "seed": seed}
    if not input_reaches_output(m):
        out["class"] = "no response"
        return out
    L = conservation_laws(m, p, m.initial[None, :], rng)
    out["conservation_laws"] = int(len(L))
    st = steady_state(m, p, m.initial, L)
    out["steady_state"] = {v: float(x) for v, x in zip(m.variables, st["q"])}
    out["y0_steady"] = st["y"]
    out["rates"] = [complex(z).real for z in st["rates"]]
    out["rates_imag"] = [complex(z).imag for z in st["rates"]]
    if not (st["converged"] and st["stable"]):
        out["class"] = "unstable loop" if st["equilibrium"] else "no set point"
        out["oscillation_of_output"] = st["oscillation"]
        return out
    q = st["q"]
    g = static_gain(m, p, q, L)
    g_ref = reference_gain(m, p, q, L)
    # S = (u/y) G is defined only where y is not zero on the scale of the open-loop response, G_ref |u|
    defined = abs(st["y"]) > 1e-9 * g_ref * abs(m.u0) and m.u0 != 0
    out.update(G=g["G"], S=g["G"] * m.u0 / st["y"] if defined else None, G_reference=g_ref,
               gain_ratio=abs(g["G"]) / g_ref if g_ref > 0 else None,
               bordered=g["bordered"], G_finite_difference=_fd_gain(m, p, q, L))
    if g_ref == 0:
        out["class"] = "no response"
        return out
    centers = np.vstack([q, m.initial])
    integ = find_integrator(m, p, centers, rng, y0_steady=st["y"])
    out["integrator"] = integ
    out["attenuation"] = attenuation(m, p, q, L)

    # ---- robustness over parameter points
    rows = []
    for pp in sample_parameters(m, samples, rng):
        s = steady_state(m, pp, q, L)
        row = {"stable": bool(s["converged"] and s["stable"])}
        if row["stable"]:
            gs = static_gain(m, pp, s["q"], L)
            ref = reference_gain(m, pp, s["q"], L)
            row["S"] = gs["S"]
            row["gain_ratio"] = abs(gs["G"]) / ref if ref > 0 else None
            # the stage that found the integrator at the stated parameters (stage 1 when none was found there)
            st_nom = integ.get("stage") or 1
            row["integrator"] = bool(find_integrator(m, pp, np.vstack([s["q"], q]), rng, y0_steady=s["y"],
                                                     count=200, stage1=st_nom == 1, stage2=st_nom == 2,
                                                     stage3=st_nom == 3)["found"])
        rows.append(row)
    stable_rows = [r for r in rows if r["stable"]]
    out["robustness"] = {"samples": samples, "stable": len(stable_rows),
                         "gain_zero": sum(1 for r in stable_rows if r["gain_ratio"] is not None
                                          and r["gain_ratio"] < GAIN_ZERO),
                         "integrator_found": sum(1 for r in stable_rows if r["integrator"]),
                         "gain_ratio_median": float(np.median([r["gain_ratio"] for r in stable_rows
                                                               if r["gain_ratio"] is not None]))
                         if stable_rows else None}

    # ---- steps
    stepped = m.steps[1:2] if first_step_only else m.steps[1:]
    resp = [step_response(m, p, q, u1, integ, L, trace=trace) for u1 in stepped] if steps else []
    out["step_responses"] = resp

    # ---- class
    S0 = out["gain_ratio"] is not None and out["gain_ratio"] < GAIN_ZERO
    rob = out["robustness"]
    if S0 and rob["stable"] and rob["gain_zero"] == rob["stable"]:
        out["class"] = "perfect-adaptation"
        out["note"] = (None if integ["found"] else
                       "no integrator in the searched coordinates (linear, logarithmic, gain linear in the state)")
    elif S0:
        out["class"] = "fine-tuned-adaptation"
    elif not steps:
        out["class"] = "partial or no adaptation (no step response computed)"
    else:
        ratios = [r["final_over_peak"] for r in resp if r.get("final_over_peak") is not None]
        returned = bool(ratios) and max(abs(r) for r in ratios) < 1 - RETURN_MIN
        out["class"] = "partial-adaptation" if returned else "no-adaptation"
        if returned:
            out["remaining_fraction_of_peak"] = [float(r) for r in ratios]
    return out


def summary(res: Dict[str, object]) -> str:
    """One line: class, S, integrator and the calibration of the first step."""
    parts = [f"{res['name']}: {res['class']}"]
    if res.get("gain_ratio") is not None:
        parts.append(f"G/G_open = {res['gain_ratio']:.3g}")
    if res.get("S") is not None:
        parts.append(f"S = {res['S']:.3g}")
    integ = res.get("integrator") or {}
    if integ.get("found"):
        w = ", ".join(f"{c:+.4g} {v}" for v, c in integ["coefficients"]["w"].items())
        v = ", ".join(f"{c:+.4g} ln {x}" for x, c in integ["coefficients"]["v"].items())
        parts.append(f"phi = {' '.join(filter(None, [w, v]))}; stage {integ['stage']}; y0 = {integ['set_point']:.6g}")
    rob = res.get("robustness")
    if rob:
        parts.append(f"G = 0 at {rob['gain_zero']}/{rob['stable']} stable samples")
    resp = res.get("step_responses") or []
    cal = [r.get("calibration_ratio") for r in resp if r.get("calibration_ratio") is not None]
    if cal:
        parts.append(f"calibration |ratio - 1| <= {max(abs(c - 1) for c in cal):.1e}")
    return "; ".join(parts)
