"""Predictions for a body that grows and divides, from its equations, before any lineage with noise.

- The threshold: the size L_c at which the symmetric state loses stability, with a, b and the critical vectors
  (reduce.reduce_at_threshold); b > 0 is a supercritical pitchfork, the case of the law.
- The dip: whether division puts the daughter below L_c (L_div / 2 < L_c).
- The linear gain over a generation, ln G = int lambda dt along the symmetric state (lineage.ln_gain): the order
  survives a generation without noise only if ln G > 0.
- The noise-free lineage: the order at successive divisions, and how far division moves the daughter away from the
  symmetric state of its size in directions other than the order.
- Lambda = a ramp / (b D_s) at L_c: the law of inheritance needs Lambda >> 1.

Classes (by order of the checks):
    kept-above-threshold          the daughter is born above L_c: no dip; the order is kept behind its barrier
    lost-in-the-dip               ln G < 0: the dip removes more order than the regrowth restores; lost without noise
    threshold-moved               ln G > 0 along the symmetric state, but the noise-free lineage loses the order: division
                                  moves the daughter away from the symmetric state whose threshold it would cross
                                  (an order made of a redistributed conserved amount: the halves differ in that amount)
    inherited-through-threshold   the daughter crosses the threshold again and keeps the parent's sign with
                                  P = Phi(phi_c / sigma_c)
"""
from __future__ import annotations

from typing import Dict

import numpy as np

from .lineage import law_for_daughter, ln_gain, noise_free_lineage, noise_s
from .reduce import reduce_at_threshold
from .spec import SpecError

CLASSES = ("inherited-through-threshold", "kept-above-threshold", "lost-in-the-dip", "threshold-moved")


def reduction(lin) -> Dict:
    lo, hi = lin.search
    body = lin.body
    try:
        q0 = body.steady(lin.q_sym0, lo)
        red = reduce_at_threshold(body, q0, lo, hi)
    except ValueError as error:
        raise SpecError(str(error)) from error
    return red


def predict(lin, red: Dict = None, generations: int = 3) -> Dict:
    red = red or reduction(lin)
    Lc = red["L_c"]
    L_div = lin.L_div(Lc)
    out = {"L_c": Lc, "a": red["a"], "b": red["b"], "supercritical": red["b"] > 0, "L_div": L_div,
           "dip": L_div / 2 < Lc, "ramp": lin.growth.ramp(Lc), "D_s": noise_s(lin, red, Lc)}
    out["Lambda"] = (red["a"] * out["ramp"] / (red["b"] * out["D_s"])
                     if red["b"] > 0 and out["D_s"] > 0 else float("inf"))
    if not out["dip"]:
        out["lnG"] = None
        out["class"] = "kept-above-threshold"
    else:
        out["lnG"] = ln_gain(lin, red, L_div)
    path = noise_free_lineage(lin, red, L_div, generations=generations)
    at_div = [abs(path["order"][k]) for k in range(len(path["order"])) if path["L"][k] == L_div]
    out["order_at_divisions"] = at_div
    out["lineage_ratio"] = at_div[-1] / at_div[0] if at_div[0] > 0 else 0.0
    # the daughter of a parent that left the crossing with the noise amplitude: displacement from the symmetric
    # state of its size outside the critical direction, and the law for it
    born = path["births"][0]
    q_sym0 = lin.body.steady(red["q_sym"], L_div / 2)
    d = born - q_sym0
    d -= (red["w"] @ d) * red["v"]
    out["displacement"] = float(np.linalg.norm(d) / max(np.linalg.norm(q_sym0), 1.0))
    if out["dip"]:
        if out["lnG"] <= 0:
            out["class"] = "lost-in-the-dip"
        elif out["lineage_ratio"] < 1e-3:
            out["class"] = "threshold-moved"
        else:
            out["class"] = "inherited-through-threshold"
    if out["class"] == "inherited-through-threshold":
        out["law"] = law_for_daughter(lin, red, born, L_div)
    out["path"] = {k: path[k] for k in ("t", "L", "order")}
    return out
