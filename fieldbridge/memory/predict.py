"""Predictions without calculation: what a realization's structure implies about its memory.

The structure is read from the drift without integrating it: the signs of the interaction graph (off-diagonal
Jacobian entries on sampled states), whether the couplings are reciprocal (symmetric Jacobian), the symmetry group
(permutations of the variables, with sign flips on R^n and reflections on angles, that commute with the drift),
continuous symmetries (directions along which the drift does not change), and what the control parameter does
(multiplies the drift, adds a bias, or reshapes it). For declared networks the transports of the edges are used
directly (networks.predict).

Each prediction states the rule it applies and the structural input it uses, and gives a control that would
falsify it. The rules are classical in their fields (loop signs: Thomas 1981, Soule 2003; monotone systems:
Hirsch 1985, Smith 1995; symmetry and bifurcation: Golubitsky and Stewart 2002; Goldstone directions); what is
predicted is existence and class, not numbers.
"""
from __future__ import annotations

import itertools
from typing import Dict, List, Optional, Tuple

import numpy as np

from . import analysis as an
from .identity import Realization


def _samples(real: Realization, rng, n: int) -> np.ndarray:
    return real.carrier.sample(rng, n)


def interaction_signs(real: Realization, rng, n: int = 24, tol: float = 1e-9) -> np.ndarray:
    """Sign of dF_i/dq_j (i != j) where it is the same on every sampled state; 0 if absent, 2 if it changes."""
    qs = _samples(real, rng, n)
    J = np.stack([an.jacobian(real, q) for q in qs])
    scale = max(1e-12, float(np.abs(J).max()))
    pos = np.all(J > tol * scale, axis=0)
    neg = np.all(J < -tol * scale, axis=0)
    zero = np.all(np.abs(J) <= tol * scale, axis=0)
    S = np.where(pos, 1, np.where(neg, -1, np.where(zero, 0, 2)))
    np.fill_diagonal(S, 0)
    return S


def reciprocal(real: Realization, rng, n: int = 16, tol: float = 1e-6) -> bool:
    """A symmetric Jacobian on every sampled state: the flow is a gradient flow (reciprocal couplings)."""
    for q in _samples(real, rng, n):
        J = an.jacobian(real, q)
        if np.max(np.abs(J - J.T)) > tol * max(1.0, float(np.abs(J).max())):
            return False
    return True


def symmetries(real: Realization, rng, n: int = 12, tol: float = 1e-7, max_dim: int = 6) -> List[Dict[str, object]]:
    """Linear symmetries g (permutation with signs) with F(g q) = g F(q) on sampled states.

    On R^n a variable may change sign; on angles a sign change is the reflection theta -> -theta; on concentrations
    only permutations are admissible. Above max_dim only cyclic shifts and the reversal are tried."""
    d = real.carrier.dim
    qs = _samples(real, rng, n)
    F = real.F(qs)
    scale = max(1e-12, float(np.abs(F).max()))
    perms = (list(itertools.permutations(range(d))) if d <= max_dim else
             [tuple(np.roll(np.arange(d), k)) for k in range(d)] + [tuple(np.arange(d)[::-1])])
    signs = [np.ones(d)] + ([-np.ones(d)] if real.carrier.kind in ("euclid", "torus") else [])
    if real.carrier.kind == "euclid" and d <= 4:
        signs = [np.array(s, float) for s in itertools.product((1, -1), repeat=d)]
    found = []
    for perm in perms:
        for sgn in signs:
            if all(p == i for i, p in enumerate(perm)) and np.all(sgn > 0):
                continue
            gq = sgn[None, :] * qs[:, perm]
            if real.carrier.kind == "torus":
                gq = real.carrier.wrap(gq)
            gF = sgn[None, :] * F[:, perm]
            if np.max(np.abs(real.F(gq) - gF)) <= tol * scale:
                order = _order(perm, sgn)
                found.append({"permutation": list(perm), "signs": sgn.astype(int).tolist(), "order": order})
    if real.carrier.kind == "torus":
        # a uniform shift by half a period: every orientation turned by the same angle
        shifted = real.carrier.wrap(qs + 0.5 * real.carrier.period)
        if np.max(np.abs(real.F(shifted) - F)) <= tol * scale:
            found.append({"shift": 0.5 * real.carrier.period, "order": 2,
                          "meaning": "every angle turned by half a period; stored states come in partner pairs"})
    return found


def _order(perm, sgn) -> int:
    d = len(perm)
    P = np.eye(d)[list(perm)] * np.asarray(sgn)[:, None]
    M = P.copy()
    for k in range(1, 25):
        if np.allclose(M, np.eye(d)):
            return k
        M = M @ P
    return 0


def continuous_symmetries(real: Realization, rng, n: int = 8, tol: float = 1e-7) -> List[List[int]]:
    """Constant directions v (entries 0 or +-1) with F(q + a v) = F(q) for all sampled q and a."""
    d = real.carrier.dim
    if d > 12:
        return []
    qs = _samples(real, rng, n)
    F = real.F(qs)
    scale = max(1e-12, float(np.abs(F).max()))
    shifts = (0.37, 1.3)
    out = []
    cands = ([np.array((1,) + s, float) for s in itertools.product((1, -1), repeat=d - 1)]
             if real.carrier.kind == "torus" else [np.ones(d)] + list(np.eye(d)))
    for v in cands:
        if real.carrier.kind == "orthant":
            continue
        ok = all(np.max(np.abs(real.F(real.carrier.wrap(qs + a * v[None])) - F)) <= tol * scale for a in shifts)
        if ok:
            out.append(v.astype(int).tolist())
    return out


def linear(real: Realization, rng, n: int = 8) -> bool:
    """F(a q1 + b q2) = a F(q1) + b F(q2) on sampled states (a linear, homogeneous drift)."""
    if real.carrier.kind != "euclid":
        return False
    q1, q2 = _samples(real, rng, n), _samples(real, rng, n)
    F = real.F(0.7 * q1 - 1.3 * q2)
    G = 0.7 * real.F(q1) - 1.3 * real.F(q2)
    return bool(np.max(np.abs(F - G)) <= 1e-9 * max(1.0, float(np.abs(G).max())))


def control_role(real: Realization, rng, n: int = 8) -> Optional[str]:
    if not real.control:
        return None
    if an.is_scale_control(real, rng):
        return "scale"
    qs = _samples(real, rng, n)
    v = float(real.params[real.control]) or 1.0
    d1 = real.F(qs, **{real.control: 1.5 * v}) - real.F(qs)
    d2 = real.F(qs + 0.1 * real.carrier.scale, **{real.control: 1.5 * v}) - real.F(qs + 0.1 * real.carrier.scale)
    if np.max(np.abs(d1 - d1[0][None])) < 1e-9 * max(1.0, float(np.abs(d1).max())) and np.allclose(d1, d2):
        return "bias"
    return "shape"


def _cycles(S: np.ndarray, max_len: int = 8) -> List[Tuple[List[int], int]]:
    n = S.shape[0]
    cycles = []
    for start in range(n):
        stack = [(start, [start], 1)]
        while stack:
            node, path, sign = stack.pop()
            for nxt in range(n):
                s = S[nxt, node]  # dF_nxt / dq_node: node acts on nxt
                if s == 0 or s == 2:
                    continue
                if nxt == start:
                    cycles.append((list(path), sign * int(s)))
                elif nxt > start and nxt not in path and len(path) < max_len:
                    stack.append((nxt, path + [nxt], sign * int(s)))
    return cycles


def _sign_representation(group: List[Dict[str, object]]) -> bool:
    """True if some symmetry can act as -1 on a mode: an element of even order, or a sign-reversing one."""
    return any(g["order"] % 2 == 0 for g in group if g["order"] > 0)


def predict(real: Realization, rng=None) -> Dict[str, object]:
    rng = rng if rng is not None else np.random.default_rng(0)
    declared = None
    if getattr(real, "network", None) is not None:
        from . import networks
        declared = networks.predict(real)
    if linear(real, rng):
        return {"source": "structure", "linear": True, "predictions": {
            "holding": {"prediction": "a single state: memory only as a kernel of unobserved coordinates (Law 1)",
                        "rule": "a linear drift has one steady state and no write point",
                        "uses": "the drift is linear in the state", "falsifier": "two stable states"}}}
    S = interaction_signs(real, rng)
    recip = reciprocal(real, rng)
    group = symmetries(real, rng)
    flat = continuous_symmetries(real, rng)
    role = control_role(real, rng)
    cyc = _cycles(S) if real.carrier.dim > 1 else []
    mixed = bool(np.any(S == 2))
    positive = [c for c in cyc if c[1] > 0]
    negative = [c for c in cyc if c[1] < 0]
    preds: Dict[str, Dict[str, object]] = {}
    # storage and oscillation
    if recip:
        preds["oscillation"] = {"prediction": "impossible", "rule": "a gradient flow has no cycles",
                                "uses": "symmetric Jacobian (reciprocal couplings)",
                                "falsifier": "a sustained oscillation in simulation"}
    elif real.carrier.dim > 1 and not mixed:
        preds["oscillation"] = {"prediction": ("possible (a negative loop)" if negative else
                                               "not expected (no negative loop)"),
                                "rule": "a negative loop is necessary for sustained oscillation (Thomas)",
                                "uses": "signs of the interaction loops", "falsifier": "an oscillation without a negative loop"}
    if real.carrier.dim > 1 and not mixed and not recip:
        preds["multistability"] = {"prediction": ("possible (a positive loop)" if positive else
                                                  "impossible: at most one steady state"),
                                   "rule": "a positive loop is necessary for multistability (Thomas; Soule 2003)",
                                   "uses": "signs of the interaction loops",
                                   "falsifier": "two stable states without a positive loop"}
    # write kind from the symmetry group
    if role == "scale":
        preds["write"] = {"prediction": "the control changes no state, only the depth of the landscape against the "
                                        "noise; a state is written by a field, which removes the stored state at a "
                                        "threshold field",
                          "rule": "a control that multiplies the whole drift leaves every fixed point in place",
                          "uses": "the drift is proportional to the control",
                          "falsifier": "a fixed point that moves or vanishes as the control changes"}
    elif not group:
        hopf = (not recip) and bool(negative)
        preds["write"] = {"prediction": ("one-sided writes at folds; a symmetric write needs one tuned parameter (a cusp)"
                                         + ("; with a negative loop the state can instead start to oscillate (Hopf)"
                                            if hopf else "")),
                          "rule": "without symmetry, a state generically appears or vanishes at a fold (codimension one); "
                                  "a pitchfork is codimension two"
                                  + ("; a non-reciprocal negative loop allows a Hopf crossing" if hopf else ""),
                          "uses": "no symmetry of the drift was found",
                          "falsifier": "a pitchfork write point without tuning a second parameter"}
    elif _sign_representation(group):
        preds["write"] = {"prediction": "a symmetric state that loses stability along a mode the symmetry reverses "
                                        "does so at a pitchfork, where a weak bias decides the written state; no second "
                                        "parameter needs tuning",
                          "rule": "a symmetry that reverses a mode makes the pitchfork a codimension-one loss of "
                                  "stability of a symmetric state; whether such a loss occurs in the range of the "
                                  "control is not decided by the structure",
                          "uses": f"{len(group)} symmetry element(s), some of even order",
                          "falsifier": "a fold at a symmetric state along a mode the symmetry reverses"}
    else:
        preds["write"] = {"prediction": "no symmetric pitchfork; a symmetric state loses stability to a rotating pair "
                                        "(Hopf) or at folds",
                          "rule": "a cyclic symmetry of odd order has no representation that reverses a single mode",
                          "uses": f"{len(group)} symmetry element(s), all of odd order",
                          "falsifier": "a pitchfork of the symmetric state"}
    # loss law and the relation between retention and writing times
    osc_possible = (not recip) and bool(negative)
    if osc_possible and not positive:
        preds["holding"] = {"prediction": "no two stored states; if the loop oscillates, its phase is a flat direction "
                                          "that retains the timing of a pulse (Law 2; the ratio of retention time to "
                                          "writing time grows linearly with the write strength)",
                            "rule": "a negative loop allows a limit cycle, and time-translation symmetry leaves its phase "
                                    "without a restoring force",
                            "uses": "signs of the interaction loops", "falsifier": "a restoring force on the phase"}
    elif flat:
        preds["holding"] = {"prediction": "a family of states along a flat direction: Law 2 (information ~ 1/t); the ratio "
                                          "of retention time to writing time grows linearly with the work of the write",
                            "rule": "a continuous symmetry leaves a direction with no restoring force",
                            "uses": f"invariant direction(s) {flat}", "falsifier": "a restoring force along that direction"}
    else:
        preds["holding"] = {"prediction": "discrete states: activated loss (Law 3) if two or more coexist, relaxation "
                                          "(Law 1) otherwise; behind a barrier the ratio of retention time to writing "
                                          "time grows exponentially with the work of the write",
                            "rule": "no continuous symmetry: every state is isolated",
                            "uses": "no invariant direction found", "falsifier": "a flat direction at a stored state"}
    # role of the control in writing
    esc = {"scale": "raising the control while the field acts changes every barrier together and is not limited "
                    "in this way",
           "shape": "sweeping the control through the write point writes at the instability and is not limited in "
                    "this way",
           "bias": "the control biases the choice between states, so that writing by the control alone is one-sided"}
    preds["lock"] = {"prediction": "a change of mobility alone rescales retention and rewriting times by the same "
                                   "factor" + (f"; {esc[role]}" if role else ""),
                     "rule": "retention and rewriting are governed by the same landscape and noise",
                     "uses": f"control role: {role or 'none'}",
                     "falsifier": "a ratio of retention to rewriting time changed by the mobility alone"}
    if declared is not None:
        preds["network"] = {"prediction": declared, "rule": "holonomy of the declared bond transports",
                            "uses": declared.get("uses", "declared edges"),
                            "falsifier": "a satisfied ground state where a flux is nonzero, or the reverse"}
    return {"source": "structure", "reciprocal": recip, "interaction_signs": S.tolist(), "mixed_signs": mixed,
            "loops": [{"nodes": c[0], "sign": c[1]} for c in cyc], "symmetries": group,
            "continuous_symmetries": flat, "control_role": role,
            "period": real.carrier.period if real.carrier.kind == "torus" else None, "predictions": preds}


PITCHFORKS = ("supercritical pitchfork", "subcritical pitchfork")


def _reversed_at(pred: Dict[str, object], event: Dict[str, object], tol: float = 1e-3) -> bool:
    """True if a symmetry of the drift fixes the state of a write point and reverses its unstable mode."""
    q = np.asarray(event.get("state", []), float)
    v = np.asarray(event.get("mode", []), float)
    if q.size == 0 or v.size != q.size:
        return False
    period = pred.get("period")
    for g in pred.get("symmetries", []):
        if "permutation" not in g:
            continue
        perm, sgn = list(g["permutation"]), np.asarray(g["signs"], float)
        dq = sgn * q[perm] - q
        if period:
            dq = (dq + 0.5 * period) % period - 0.5 * period
        if (np.linalg.norm(dq) <= tol * max(1.0, float(np.linalg.norm(q)))
                and np.linalg.norm(sgn * v[perm] + v) <= 1e-2 * np.linalg.norm(v)):
            return True
    return False


def compare(pred: Dict[str, object], card: Dict[str, object]) -> Dict[str, object]:
    """Check the structural predictions against a calculated card (discovery.evaluate)."""
    out = {}
    P = pred.get("predictions", {})
    ev = [e for e in card["construct"]["events"] if e["source"] == "control"]
    kinds = {e["kind"].split(":")[0].split(" (")[0] for e in ev}
    if "write" in P and P["write"]["prediction"].startswith("the control changes no state"):
        out["write"] = {"predicted": "no write point along the control", "calculated": sorted(kinds) or ["none"],
                        "consistent": not ev}
    elif "write" in P and ev:
        text = P["write"]["prediction"]
        if text.startswith("a symmetric state that loses stability"):
            # the prediction is conditional: it concerns write points at a symmetric state along a reversed mode
            tested = [e for e in ev if not e.get("hopf") and _reversed_at(pred, e)]
            ok = all(e["kind"].startswith(PITCHFORKS) for e in tested)
            out["write"] = {"predicted": text, "calculated": sorted(kinds), "tested_write_points": len(tested),
                            "consistent": bool(ok)}
        else:
            ok = not any(k.startswith(PITCHFORKS) for k in kinds)
            out["write"] = {"predicted": text, "calculated": sorted(kinds), "consistent": bool(ok)}
    st = card["states"]
    if "oscillation" in P:
        osc = st["count"] == 0 and st["unconverged"] > 0
        text = P["oscillation"]["prediction"]
        ok = (not osc) if (text.startswith("impossible") or text.startswith("not expected")) else True
        out["oscillation"] = {"predicted": text, "calculated": "oscillates" if osc else "settles", "consistent": bool(ok)}
    if "multistability" in P:
        text = P["multistability"]["prediction"]
        ok = st["count"] <= 1 if text.startswith("impossible") else True
        out["multistability"] = {"predicted": text, "calculated": f"{st['count']} stable state(s)", "consistent": bool(ok)}
    if "holding" in P and not pred.get("linear"):
        flat_pred = P["holding"]["prediction"].startswith(("a family", "no two stored states"))
        law = card["loss_law"]
        flat_calc = "Law 2" in law
        phase_pred = P["holding"]["prediction"].startswith("no two stored states")
        ok = (flat_pred == flat_calc) or (phase_pred and ("limit cycle" in law or "Law 1" in law))
        out["holding"] = {"predicted": "Law 2" if flat_pred else "Law 1 or 3", "calculated": law, "consistent": bool(ok)}
    out["all_consistent"] = all(v["consistent"] for v in out.values() if isinstance(v, dict))
    return out
