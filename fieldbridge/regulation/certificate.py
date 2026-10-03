"""Certificates: what a card establishes, recorded so that it can be checked without repeating the search.

The search for an integrator explores sampled states; checking a recorded integrator does not. A certificate records
the class, the steady state, the integrator (stage, coefficients of phi, gain, set point) and the gain ratio.
``check`` takes a report and the specification it was computed from and verifies, in seconds and without a search:

  1. the specification is the one the report read (hash of the input);
  2. the recorded steady state solves F = 0 to Newton accuracy at the reference input, and it is stable;
  3. the recorded integrator satisfies its identity at fresh samples of states and inputs, and symbolically where the
     drift and the output are rational functions (the identity is then exact for every state, not sampled);
  4. the gain ratio G/G_open at the steady state agrees with the class (zero for perfect adaptation).

A check that passes shows that the report follows from the specification. It does not show that the specification
describes a material.
"""
from __future__ import annotations

from typing import Dict, List

import numpy as np
import sympy as sp

from .card import GAIN_ZERO
from .gains import reference_gain, static_gain
from .integrator import CHECK_TOL, conservation_laws, rate, sample_inputs, sample_states
from .spec import Regulated, input_hash
from .steady import rates, refine


def certificate(res: Dict) -> Dict[str, object]:
    """The part of a card that a check can verify."""
    integ = res.get("integrator") or {}
    keep = ("found", "stage", "gain", "coefficients", "k_I", "gain_coefficients", "set_point")
    return {"class": res["class"], "input": res["input"], "output": res["output"], "u0": res["u0"],
            "steady_state": res.get("steady_state"), "gain_ratio": res.get("gain_ratio"),
            "integrator": {k: integ[k] for k in keep if k in integ} if integ.get("found") else {"found": False}}


def _rationalize(x: float):
    return sp.nsimplify(x, rational=True, tolerance=1e-10 * max(1.0, abs(x)))


def _symbolic_identity(m: Regulated, integ: Dict) -> str:
    """'exact' if the identity simplifies to zero with rational coefficients, 'not exact' if it does not, and
    'not attempted' for stage 3 or for drifts that are not rational functions of the state."""
    if integ.get("stage") not in (1, 2):
        return "not attempted"
    params = {s: _rationalize(m.params[str(s)]) for s in m.p}
    F = [f.subs(params) for f in m.F]
    h = m.h.subs(params)
    if not all(e.is_rational_function(*m.q) for e in F + [h]):
        return "not attempted"
    index = {str(s): i for i, s in enumerate(m.q)}
    lhs = sum(_rationalize(c) * F[index[v]] for v, c in integ["coefficients"]["w"].items())
    lhs += sum(_rationalize(c) * F[index[v]] / m.q[index[v]] for v, c in integ["coefficients"]["v"].items())
    if integ["stage"] == 1:
        gain = _rationalize(integ["k_I"])
    else:
        g = integ["gain_coefficients"]
        gain = _rationalize(g.get("1", 0.0)) + sum(_rationalize(c) * m.q[index[v]] for v, c in g.items() if v != "1")
    expr = sp.cancel(sp.together(lhs - gain * (h - _rationalize(integ["set_point"]))))
    return "exact" if expr == 0 else "not exact"


def check(report: Dict, m: Regulated, seed: int = 1) -> Dict[str, object]:
    cert = report["certificate"]
    rng = np.random.default_rng(seed)
    out: List[Dict[str, object]] = []
    out.append({"check": "the specification is the one the report read",
                "passed": report.get("input_sha256") == input_hash(m.spec)})
    if cert["class"] in ("no response", "no set point", "unstable loop"):
        out.append({"check": "class without a steady state is not checked further", "passed": True})
        return {"passed": all(c["passed"] for c in out), "checks": out}
    p = m.pvec(m.u0)
    q = np.array([cert["steady_state"][v] for v in m.variables], float)
    L = conservation_laws(m, p, q[None, :], rng)
    fine = refine(m, p, q, L, L @ q if len(L) else None)
    lam = rates(m, p, fine["q"], L)
    moved = float(np.linalg.norm(fine["q"] - q) / max(np.linalg.norm(q), 1e-300))
    out.append({"check": "the recorded steady state solves F = 0 and is stable", "moved": moved,
                "residual": fine["residual"], "largest_rate": float(lam.real.max()) if lam.size else None,
                "passed": moved < 1e-6 and (lam.size == 0 or lam.real.max() < 0)})
    integ = cert["integrator"]
    if integ.get("found"):
        Q = sample_states(m, q[None, :], 300, rng, spread=1.5)
        U = sample_inputs(m, 300, rng, spread=1.5)
        worst = 0.0
        if integ["stage"] in (1, 2):
            for qq, uu in zip(Q, U):
                pp = p.copy()
                pp[m.u_index] = uu
                with np.errstate(all="ignore"):
                    F = m.f(qq, pp)
                    terms = [c * F[m.variables.index(v)] for v, c in integ["coefficients"]["w"].items()]
                    terms += [c * F[m.variables.index(v)] / qq[m.variables.index(v)]
                              for v, c in integ["coefficients"]["v"].items()]
                    rhs = rate(integ, m, qq, pp, m.y(qq, pp))
                if not (np.all(np.isfinite(terms)) and np.isfinite(rhs)):
                    continue
                worst = max(worst, abs(sum(terms) - rhs) / max(sum(abs(t) for t in terms) + abs(rhs), 1e-300))
            numeric = worst < CHECK_TOL
        else:
            # stage 3: the rate w.F is a function of the output alone, of one sign on each side of the set point
            w = np.array([integ["coefficients"]["w"].get(v, 0.0) for v in m.variables])
            P = np.repeat(p[None, :], len(Q), axis=0)
            P[:, m.u_index] = U
            Z = m.f_batch(Q, P) @ w
            Y = m.y_batch(Q, P)
            ok = np.isfinite(Z) & np.isfinite(Y) & (np.abs(Y - integ["set_point"]) > 1e-6)
            side = np.sign(Z[ok]) * np.sign(Y[ok] - integ["set_point"])
            numeric = bool(np.all(side > 0))
            worst = float(np.mean(side <= 0))
        out.append({"check": "the recorded integrator satisfies its identity at fresh samples", "stage": integ["stage"],
                    "worst_relative_residual": worst, "passed": bool(numeric)})
        sym = _symbolic_identity(m, integ)
        out.append({"check": "the identity simplifies to zero with rational coefficients", "result": sym,
                    "passed": sym != "not exact"})
    G = static_gain(m, p, fine["q"], L)["G"]
    ref = reference_gain(m, p, fine["q"], L)
    ratio = abs(G) / ref if ref > 0 else None
    if cert["class"] == "perfect-adaptation" or cert["class"] == "fine-tuned-adaptation":
        agrees = ratio is not None and ratio < GAIN_ZERO
    else:
        agrees = ratio is not None and ratio >= GAIN_ZERO and (
            cert["gain_ratio"] is None or abs(ratio - cert["gain_ratio"]) <= 1e-6 * max(cert["gain_ratio"], 1e-12))
    out.append({"check": "the gain ratio at the steady state agrees with the class", "gain_ratio": ratio,
                "passed": bool(agrees)})
    return {"passed": all(c["passed"] for c in out), "checks": out}
