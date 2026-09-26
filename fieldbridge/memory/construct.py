"""The construction step, run automatically on a target realization.

  1  carrier      the target's state space Xi_t (given by the realization).
  2  generator    where the target's drift supplies kappa < 0: write points located by continuation along its
                  control parameter and along a bounded write field; at each, the drift along the critical
                  direction reduced onto the slow manifold, ds/dt = a0 + a1 s + a2 s^2 + a3 s^3.
  3  obstruction  Delta, the difference between that drift and the pitchfork generator eps s - s^3 + h:
                  a quadratic term (a one-sided write), a bias carried by the control, a positive cubic (no
                  saturation near the write point), or a rotating crossing (Hopf). -Delta is what the material
                  must add. When the realization names further material parameters, the point at which the
                  quadratic term vanishes (a cusp) is solved for, together with the direction of a sweep through
                  it that carries no bias; the normal form along that sweep is then recomputed.

For a write point that realizes the generator, the swept-write law of the source,
P = Phi(pi^{1/4} h_s / (sqrt(D_s) r^{1/4})) (Kondepudi and Nelson), is transferred with h_s and D_s projected on
the critical direction and r the rate at which the sweep drives a1, and checked by simulating the target.
"""
from __future__ import annotations

from dataclasses import replace
from math import erf, sqrt
from typing import Dict, List, Optional

import numpy as np

from . import analysis as an
from .identity import Realization

TOL = 0.1


def obstruction(nf: Dict) -> Dict[str, object]:
    """Difference between the target's reduced drift and the pitchfork generator, and what cancels it."""
    if nf.get("hopf"):
        return {"hopf": True, "asymmetric": False, "biased": False, "subcritical": False, "realizes_generator": False,
                "terms": {"rotation frequency": nf["omega"]},
                "required": ["The critical eigenvalues cross as a complex pair: the unwritten state starts to oscillate "
                             "and no state is chosen. A memory needs the loop of couplings to compose to a map with a "
                             "fixed point (a positive loop), or reciprocal couplings."]}
    asym, biased = nf["asymmetry"] > TOL, nf["bias"] > TOL
    sub = (not asym) and nf["a3"] > 0
    terms, req = {}, []
    if asym:
        terms["quadratic a2"] = nf["a2"]
        req.append("A quadratic term makes the write one-sided (a threshold): the material must supply the mirror "
                   "term -a2 s^2, from a symmetric partner of the critical mode or a second parameter tuned to the cusp.")
    if biased:
        terms["bias per unit control"] = nf["bias_along_mode"]
        req.append("Changing the control also pushes the critical mode, so a sweep writes a fixed state: the protocol "
                   "must add the opposite field, or the material must be symmetric under reversal of the mode.")
    if sub:
        terms["cubic a3 > 0"] = nf["a3"]
        req.append("The cubic term does not saturate the instability, so the state jumps to a distant branch: a "
                   "saturating term (a finite resource, Hill saturation, an anisotropy) is required for a local write.")
    if not req:
        req.append("No obstruction: along the critical direction the target realizes the pitchfork generator "
                   "eps s - |a3| s^3 + h, with eps proportional to the distance of the control from its critical value.")
    return {"hopf": False, "asymmetric": asym, "biased": biased, "subcritical": sub,
            "realizes_generator": not (asym or biased or sub), "terms": terms, "required": req}


def along_path(real: Realization, origin: Dict[str, float], direction: Dict[str, float]) -> Realization:
    """The realization with its parameters moved together, p = origin + t direction; t is the new control."""
    base = real

    def drift(q, p):
        moved = {**p, **{k: origin[k] + p["t"] * direction[k] for k in origin}}
        return base.drift(q, moved)

    r = replace(real, name=real.name + " along a combined sweep", drift=drift, potential=None,
                params={**real.params, **origin, "t": 0.0}, control="t", control_range=(-1.0, 1.0))
    return r


def _null_vectors(J):
    if not np.all(np.isfinite(J)):
        raise np.linalg.LinAlgError("non-finite Jacobian")
    U, S, Vt = np.linalg.svd(J)
    vec, w = Vt[-1], U[:, -1]
    if float(w @ vec) < 0:
        w = -w
    return vec, w


def _fold_equations(real, param, second, h):
    """Residuals of a write point (F = 0, det J = 0) and its quadratic coefficient a2 along the null direction."""
    def a2_of(q, over, vec_ref=None):
        J = an.jacobian(real, q, over)
        vec, w = _null_vectors(J)
        if vec_ref is not None and float(vec @ vec_ref) < 0:
            vec, w = -vec, -w
        F = real.F(q, **over)
        d2 = (real.F(q + h * vec, **over) + real.F(q - h * vec, **over) - 2 * F) / h ** 2
        return 0.5 * float(w @ d2) / max(abs(float(w @ vec)), 1e-12), vec

    def fold(y, pv):
        n = y.size - 1
        q, over = y[:n], {param: y[n], second: pv}
        res = np.concatenate([real.F(q, **over), [np.linalg.det(an.jacobian(real, q, over))]])
        return res if np.all(np.isfinite(res)) else np.full(n + 1, 1e3)

    return fold, a2_of


def _follow_fold(real, param, second, q0, v0, p0, h, max_steps: int = 400):
    """Follow the write point (a fold) as the second parameter changes, in both directions, until its quadratic
    coefficient a2 changes sign or nearly vanishes. |a2| need not fall at the first step, so each direction is
    followed over its whole range before the directions are compared. Returns the point nearest the cusp."""
    from scipy.optimize import root
    fold, a2_of = _fold_equations(real, param, second, h)
    y = np.concatenate([np.asarray(q0, float), [v0]])
    n = y.size - 1
    try:
        a2_0, vec0 = a2_of(y[:n], {param: v0, second: p0})
    except np.linalg.LinAlgError:
        return None
    scale = max(abs(p0), 1.0)
    best = None
    for direction in (+1.0, -1.0):
        dp = 0.01 * scale * direction
        yy, pv, a2p, vref = y.copy(), p0, a2_0, vec0
        path = [(pv, float(yy[n]), a2p)]
        for _ in range(max_steps):
            if abs(dp) < 1e-7 * scale or abs(pv - p0) > 4 * scale:
                break
            trial = pv + dp
            sol = root(fold, yy, args=(trial,), method="hybr")
            if not (sol.success and np.all(np.isfinite(sol.x)) and np.max(np.abs(sol.fun)) < 1e-7
                    and np.linalg.norm(sol.x - yy) < 0.2 * an.base_length(real) * max(1.0, np.sqrt(n))):
                dp /= 2
                continue
            try:
                a2n, vn = a2_of(sol.x[:n], {param: sol.x[n], second: trial}, vref)
            except np.linalg.LinAlgError:
                dp /= 2
                continue
            if real.carrier.kind == "orthant" and np.min(sol.x[:n]) < 1e-3 * real.carrier.scale:
                break  # the fold reached the boundary of the orthant, where the drift need not be smooth
            path.append((trial, float(sol.x[n]), a2n))
            if np.sign(a2n) != np.sign(a2p):
                # a zero crossing passes through small values; a jump through a singularity does not
                if max(abs(a2n), abs(a2p)) < 0.5 * abs(a2_0):
                    cand = {"y": sol.x if abs(a2n) < abs(a2p) else yy, "p": trial if abs(a2n) < abs(a2p) else pv,
                            "path": list(path), "a2": min(abs(a2n), abs(a2p))}
                    if best is None or cand["a2"] < abs(best["a2"]):
                        best = cand
                break
            if best is None or abs(a2n) < abs(best["a2"]):
                best = {"y": sol.x, "p": trial, "path": list(path), "a2": a2n}
            yy, pv, a2p, vref = sol.x, trial, a2n, vn
            dp *= 1.3
    if best is not None and abs(best["a2"]) < 0.2 * abs(a2_0):
        return best
    return None


def cancel_asymmetry(real: Realization, nf: Dict, second: str) -> Dict[str, object]:
    """Cancel a one-sided write by a second material parameter. The write point (F = 0, det J = 0) is followed
    as the second parameter changes until its quadratic coefficient a2 vanishes: the cusp, where the write becomes
    symmetric. The cusp is then refined from F = 0, det J = 0, a2 = 0, and the direction in (control, second)
    that keeps the bias zero and makes kappa negative is computed; the normal form along it is recomputed."""
    from scipy.optimize import least_squares, root
    param = nf["param"]
    q0 = np.asarray(nf["state"], float)
    n = q0.size
    h = 1e-3 * an.base_length(real)
    p0 = float(real.params[second])

    def eqs(x):
        q, over = x[:n], {param: x[n], second: x[n + 1]}
        try:
            F = real.F(q, **over)
            J = an.jacobian(real, q, over)
            vec, w = _null_vectors(J)
        except np.linalg.LinAlgError:
            return np.full(n + 2, 1e3)
        d2 = (real.F(q + h * vec, **over) + real.F(q - h * vec, **over) - 2 * F) / h ** 2
        a2 = 0.5 * float(w @ d2) / max(abs(float(w @ vec)), 1e-12)
        res = np.concatenate([F, [np.linalg.det(J)], [a2]])
        return res if np.all(np.isfinite(res)) else np.full(n + 2, 1e3)

    out: Dict[str, object] = {"second": second, "success": False}
    with np.errstate(all="ignore"):
        followed = _follow_fold(real, param, second, q0, nf["value"], p0, h)
        starts = []
        if followed is not None:
            path = followed["path"]
            out["fold_path"] = [{"second": a, "control": b, "a2": c} for a, b, c in path]
            y_end = np.asarray(followed["y"], float)
            if len(path) >= 3:
                # continuation stalls where the fold equations turn singular, next to the cusp: extrapolate a2 -> 0
                P, V, A = (np.array([r[i] for r in path[-3:]]) for i in range(3))
                if np.ptp(A) > 0:
                    k = np.polyfit(A, P, 1)
                    p_star = float(np.polyval(k, 0.0))
                    v_star = float(np.polyval(np.polyfit(P, V, 1), p_star))
                    starts.append(np.concatenate([y_end[:n], [v_star, p_star]]))
            starts.append(np.concatenate([y_end, [followed["p"]]]))
            for a, b, _ in path[-4:-1][::-1]:
                starts.append(np.concatenate([y_end[:n], [b, a]]))
        starts.append(np.concatenate([q0, [nf["value"], p0]]))
        lb = np.concatenate([np.full(n, 1e-9 if real.carrier.kind == "orthant" else -np.inf), [-np.inf, -np.inf]])
        x = None
        for x0 in starts:
            for method in ("hybr", "lsq"):
                try:
                    if method == "hybr":
                        sol = root(eqs, x0, method="hybr", options={"xtol": 1e-12})
                    else:
                        sol = least_squares(eqs, np.maximum(x0, lb + 1e-9), bounds=(lb, np.inf), xtol=1e-15,
                                            ftol=1e-15, gtol=1e-15)
                except (np.linalg.LinAlgError, ValueError):
                    continue
                ok = bool(np.all(np.isfinite(sol.x)) and np.max(np.abs(sol.fun)) < 1e-7
                          and (real.carrier.kind != "orthant" or np.all(sol.x[:n] > 1e-3 * real.carrier.scale)))
                if ok:
                    x = sol.x
                    break
            if x is not None:
                break
    if x is None:
        return out
    q, v, pv = x[:n], float(x[n]), float(x[n + 1])
    over = {param: v, second: pv}
    vec, w = _null_vectors(an.jacobian(real, q, over))
    w = w / float(w @ vec)
    e = 1e-6 * max(1.0, abs(v), abs(pv))
    dF = {k: (real.F(q, **{**over, k: over[k] + e}) - real.F(q, **{**over, k: over[k] - e})) / (2 * e) for k in over}
    dJ = {k: (an.jacobian(real, q, {**over, k: over[k] + e}) - an.jacobian(real, q, {**over, k: over[k] - e})) / (2 * e)
          for k in over}
    g = np.array([float(w @ dF[param]), float(w @ dF[second])])
    grow = np.array([float(w @ dJ[param] @ vec), float(w @ dJ[second] @ vec)])
    d = np.array([-g[1], g[0]]) if np.linalg.norm(g) > 1e-9 * max(1.0, np.linalg.norm(grow)) else grow.copy()
    if float(d @ grow) < 0:
        d = -d
    d = d / max(np.linalg.norm(d), 1e-300)
    path = along_path(real, {param: v, second: pv}, {param: float(d[0]), second: float(d[1])})
    nf_path = an.normal_form(path, q, "t", 0.0)
    out.update(success=True, state=q.tolist(), control=v, value=pv,
               direction={param: float(d[0]), second: float(d[1])},
               normal_form_along_sweep={k: nf_path[k] for k in ("kind", "a2", "a3", "asymmetry", "bias", "kappa_slope")
                                        if k in nf_path},
               obstruction_along_sweep=obstruction(nf_path))
    return out


def _phi(x: float) -> float:
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))


def swept_write_check(real: Realization, nf: Dict, rng, target_p: float = 0.8, rate: float = 0.02,
                      gamma: float = 0.05, n_traj: int = 400) -> Dict[str, float]:
    """Transfer the swept-write law of the pitchfork to the target and test it. The control is swept through its
    critical value so that the linear coefficient of the critical mode grows at `rate`; a bias of fixed size
    along the critical direction is chosen so that the transferred law, P = Phi(pi^{1/4} h_s / (sqrt(D_s) r^{1/4})),
    predicts target_p. It is tested at a noise set by Gamma = D_s |a3| / r (the noise relative to the growth of
    the barrier behind the choice) and at the material's own noise."""
    from statistics import NormalDist
    vc = nf["value"]
    qc, vec = np.asarray(nf["state"], float), np.asarray(nf["mode"], float)
    w = np.asarray(nf["left_mode"], float)
    g, a3 = abs(nf["kappa_slope"]), abs(nf["a3"])
    if g <= 0 or a3 <= 0:
        return {}
    # the law assumes a bias of fixed size along the critical direction: on R^n or R_+^n the write points at a
    # distant target along that direction, so that its direction does not turn while the state moves
    far = 0.5 * nf["half"] if real.carrier.kind == "torus" else 1e3 * an.base_length(real)
    target = qc + far * vec if real.carrier.kind != "torus" else real.carrier.wrap(qc + far * vec)
    h_per = float(w @ real.write_force(qc, target, 1.0))
    if abs(h_per) < 1e-12:
        return {}
    rv = rate / g
    t_half = 4.0 / sqrt(rate)
    dv = np.sign(-nf["kappa_slope"]) * rv
    path = lambda t: vc - dv * t_half + dv * t
    z = NormalDist().inv_cdf(target_p)

    def run(D):
        D_s = D * float(w @ w)
        h_s = z * sqrt(D_s) * rate ** 0.25 / np.pi ** 0.25
        q = np.repeat(qc[None], n_traj, axis=0)
        q, _, _ = an.integrate(real, q, 2 * t_half, D, rng, over={}, write=(target, h_s / h_per), control_path=path)
        m = float(np.mean((real.carrier.diff(np.repeat(qc[None], n_traj, axis=0), q) @ w) > 0))
        return {"D": float(D), "D_along_mode": D_s, "gamma": D_s * a3 / rate, "h": float(h_s / h_per),
                "h_along_mode": float(h_s), "measured": m, "stderr": float(sqrt(m * (1 - m) / n_traj))}

    D_law = gamma * rate / (a3 * float(w @ w))
    inside = run(D_law)
    own = run(real.noise)
    return {"rate": rate, "sweep_rate_control": float(rv), "predicted": float(target_p),
            "measured": inside["measured"], "stderr": inside["stderr"], "gamma": inside["gamma"],
            "h_along_mode": inside["h_along_mode"], "D_along_mode": inside["D_along_mode"],
            "at_material_noise": own}


def construct(real: Realization, rng, states: Optional[List[np.ndarray]] = None, n_scan: int = 7,
              max_dim: int = 12, check_write: bool = True) -> Dict[str, object]:
    out: Dict[str, object] = {"events": []}
    scale = an.is_scale_control(real, rng)
    if scale:
        out["control_role"] = (f"{real.control} multiplies the whole drift: it changes no state, only the depth of "
                               "the landscape against the noise (an inverse temperature)")
    if real.control and real.carrier.dim <= max_dim and not scale:
        values = np.linspace(*real.control_range, n_scan)
        for ev in an.locate_writes(real, rng, real.control, values, n_starts=24 if real.carrier.dim <= 4 else 12):
            nf = an.normal_form(real, ev["q"], real.control, ev["v"])
            nf.update(event=ev["type"], source="control")
            nf["obstruction"] = obstruction(nf)
            ob = nf["obstruction"]
            if (ob["asymmetric"] or ob["biased"]) and real.material:
                for second in real.material:
                    try:
                        with np.errstate(all="ignore"):
                            c = cancel_asymmetry(real, nf, second)
                    except np.linalg.LinAlgError:
                        c = {"second": second, "success": False}
                    if c.get("success"):
                        nf["cancellation"] = c
                        break
            if ob["realizes_generator"] and check_write:
                nf["write_law"] = swept_write_check(real, nf, rng)
            out["events"].append(nf)
    if states is not None and len(states) >= 2 and real.carrier.dim <= max_dim:
        thr = an.restoring_threshold(real, states[0], states[1])
        rw = an.with_write_field(real, states[1], (0.0, 3.0 * thr))
        branch, ev = an.track(rw, "h", states[0], 0.0, 3.0 * thr)
        wf: Dict[str, object] = {"threshold": thr, "branch": [
            {"h": b["v"], "kappa_min": b["kappa_min"],
             "overlap": float(real.overlap(np.asarray(b["q"]), states[1]))} for b in branch]}
        if ev is not None:
            nf = an.normal_form(rw, ev["q"], "h", ev["v"])
            nf.update(event=ev["type"], source="write field")
            nf["obstruction"] = obstruction(nf)
            wf["h_c"] = float(ev["v"])
            out["events"].append(nf)
            fine, _ = an.track(rw, "h", states[0], 0.0, 0.999 * float(ev["v"]), n_init=24)
            wf["branch"] = [{"h": b["v"], "kappa_min": b["kappa_min"],
                             "overlap": float(real.overlap(np.asarray(b["q"]), states[1]))} for b in fine]
        out["write_field"] = wf
    return out
