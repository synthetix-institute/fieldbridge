"""The law at an obstruction: complete the reduction that the obstruction names and derive the write law of the
completed model.

The derivation of a write stops at an obstruction when the drift at the write point does not reduce to the pitchfork
eps s - s^3 + h on one critical line (construct.obstruction). The obstruction says what the one-coordinate reduction
discarded:

  degenerate critical mode   k directions become critical together: keep all k (the critical space);
  vanishing cubic            the cubic coefficient is zero: keep the leading term of higher order;
  subcritical                the cubic does not saturate: keep the distant branch to which the state jumps;
  one-sided write            a quadratic term or a bias (a fold, a transcritical or an imperfect pitchfork): keep a
                             second parameter and sweep through the point where the quadratic term vanishes (a cusp);
  rotating crossing (Hopf)   keep the phase (the phase-locking derivation, not this module).

The completed model is not reduced to a canonical form whose law is known in advance; its write law is derived from
the sweep itself, in the same way for every completion:

  linear stage   the exact Gaussian moments of the drift linearized about the unwritten state along the sweep, with
                 the bias and the noise, up to the time at which the spread reaches a small fraction of the distance
                 to the written states (the end of the linear stage);
  outcome map    the deterministic flow of the full drift, with the bias, from samples of that Gaussian to the end of
                 the sweep; each final state is assigned to the nearest stable state there.

The probabilities of the stable states are the law. It is tested against a simulation of the full stochastic
dynamics over the whole sweep (the certificate). A deviation beyond the statistical error is itself an obstruction:
the choice is not complete at the end of the linear stage (noise acts during the nonlinear stage, or the nonlinear
stage brings another mode into play), and it is reported as such.

For one critical direction the law is reported in the units of the swept pitchfork: with the bias toward either
state, z(+-h) = Phi^-1(P(state of the bias)); the response constant c = (z(+h) + z(-h)) / 2 * sqrt(D_s) r^(1/4) / h_s
(pi^(1/4) for the pitchfork) and the offset z0 = (z(+h) - z(-h)) / 2 (zero where a symmetry reverses the critical
mode). For k critical directions the law is the table of probabilities for biases in several directions of the
critical space.
"""
from __future__ import annotations

from math import pi, sqrt
from statistics import NormalDist
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import analysis as an
from . import construct
from .identity import Realization

ND = NormalDist()
GROWTH_STOP = 5.0           # the linear stage ends once the critical rate integrated since the crossing reaches this
                            # (the fluctuations then present are amplified e^5-fold, later noise no longer counts) ...
GROWTH_MAP = 1.0            # the outcome map starts once the critical deviation grows (one e-fold since the
                            # crossing), close to where the choice is made, before terms that grow along the sweep act
LINEAR_FRACTION = 0.2       # the largest size of the mapped deviation (3 sigma) as a fraction of the distance from
                            # the unwritten state to the nearest written state
DEVIATION_SIGMA = 4.0       # law and certificate differ when a probability differs by more than this many errors ...
DEVIATION_FLOOR = 0.02      # ... and by more than this absolute amount (the law has its own sampling and time-step
                            # error; near-certain outcomes have binomial errors far below it)
CURVATURES = (-0.15, -0.1, -0.05, 0.0, 0.05)  # curvatures of the sweep through a cusp tried by
                            # zero_offset_curvature (on the cusp of a driven Kerr resonator its offset
                            # is close to linear in the curvature over this range)


# ------------------------------------------------------------------------------------------------ the completion
def completion_kind(nf: Dict, ob: Dict) -> str:
    """What the obstruction says must be kept."""
    if nf.get("hopf"):
        return "phase"
    if ob.get("degenerate") or int(nf.get("critical_dimension", 1)) > 1:
        return "critical space"
    if ob.get("degenerate_cubic"):
        return "leading order"
    if ob.get("asymmetric") or ob.get("biased"):
        return "second parameter"
    if ob.get("subcritical"):
        return "distant branch"
    return "none"


def critical_space(real: Realization, nf: Dict, k: int) -> Tuple[np.ndarray, np.ndarray]:
    """Real bases R (n x k) of the critical eigenspace at the write point and L of its dual (L^T R = I)."""
    q = np.asarray(nf["state"], float)
    over = {nf["param"]: nf["value"]}
    J = an.jacobian(real, q, over)
    def basis(M):
        ev, V = np.linalg.eig(M)
        idx = np.argsort(np.abs(ev.real))[:k]
        B = np.concatenate([np.real(V[:, idx]), np.imag(V[:, idx])], axis=1)
        U, s, _ = np.linalg.svd(B, full_matrices=False)
        return U[:, :k]
    R, L = basis(J), basis(J.T)
    return R, L @ np.linalg.inv(R.T @ L)         # dual basis: L^T R = I


def _relaxed(real: Realization, q: np.ndarray, over: Dict) -> np.ndarray:
    """The unwritten state at the start of the sweep (relaxed from the write point)."""
    qq, _ = an.relax(real, q, over)
    return np.asarray(qq, float)[0]


def stable_states(real: Realization, rng, over: Dict, n_starts: int = 40, seeds=None) -> List[np.ndarray]:
    """Stable fixed points at `over`, from random starts and from the given seeds (states that a sweep reaches: the
    written states can lie outside the box from which random starts are drawn)."""
    out = []
    for q, kap in an.fixed_points(real, rng, n_starts, over, seeds=seeds):
        if an.n_unstable(kap) == 0:
            out.append(np.asarray(q, float))
    return out


def swept_seeds(real: Realization, sweep: "Sweep", q0: np.ndarray, R: np.ndarray, t_end: float,
                dt: Optional[float] = None) -> np.ndarray:
    """End states of deterministic sweeps started from the unwritten state displaced along each critical direction,
    by several amplitudes and both signs, at the crossing and half-way after it: the states a write can reach. (A
    sweep that itself favours one state carries small displacements to that state; larger ones reach the other.)"""
    dt = dt or 0.5 * an.safe_dt(real)
    b = an.base_length(real)
    out = []
    for t0 in (sweep.t_cross, 0.5 * (sweep.t_cross + t_end)):
        if t0 >= t_end:
            continue
        q_ref = q0.copy(); t = 0.0
        while t < t0 - 1e-12:
            h = min(dt, t0 - t)
            q_ref = _rk4_step(lambda x, s_: real.F(x, **sweep.over(s_)), q_ref, t, h)
            t += h
        q = np.array([q_ref + sg * a * b * R[:, j] for j in range(R.shape[1]) for sg in (1, -1)
                      for a in (0.01, 0.1, 0.5, 2.0)])
        t = t0
        while t < t_end - 1e-12:
            h = min(dt, t_end - t)
            q = real.carrier.wrap(_rk4_step(lambda x, s_: real.F(x, **sweep.over(s_)), q, t, h))
            t += h
        out.append(q[np.all(np.isfinite(q), axis=1)])
    return np.concatenate(out) if out else np.zeros((0, q0.size))


# ------------------------------------------------------------------------------------------------ the sweep
class Sweep:
    """The control swept through the write point so that the critical rate grows at `rate` (as in
    construct.swept_write_check), from 4/sqrt(rate) before the crossing to as long after it."""

    def __init__(self, real: Realization, nf: Dict, rate: float):
        g = abs(nf["kappa_slope"])
        if g <= 0:
            raise ValueError("the critical rate does not change along the control")
        self.real, self.param, self.vc, self.rate = real, nf["param"], float(nf["value"]), rate
        self.dv = float(np.sign(-nf["kappa_slope"]) * rate / g)
        self.t_cross = 4.0 / sqrt(rate)
        self.T = 2 * self.t_cross

    def value(self, t: float) -> float:
        return self.vc - self.dv * self.t_cross + self.dv * t

    def over(self, t: float) -> Dict[str, float]:
        return {self.param: self.value(t)}


def _rk4_step(f, q, t, h):
    k1 = f(q, t); k2 = f(q + 0.5 * h * k1, t + 0.5 * h); k3 = f(q + 0.5 * h * k2, t + 0.5 * h); k4 = f(q + h * k3, t + h)
    return q + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


# ------------------------------------------------------------------------------------------------ the law
def linear_stage_law(real: Realization, sweep: Sweep, q0: np.ndarray, R: np.ndarray, L: np.ndarray,
                     targets: Sequence[np.ndarray], strengths: Sequence[float], D: float, stable: List[np.ndarray],
                     n_samples: int, rng, dt: Optional[float] = None) -> List[Dict]:
    """For each bias (target, strength): the Gaussian of the linear stage and the outcome map.

    The moments of the deviation from the unwritten state (no bias, no noise) are integrated through the crossing
    until the critical rate integrated since the crossing reaches GROWTH_STOP: the fluctuations present then decide
    the write, and noise added later is amplified less by that factor. Their distribution in the critical space gives
    the statistics of the decision (in one direction, the ratio of mean to spread; in k directions, also the angles).
    In the linear stage the critical deviation only grows in size, so the samples are mapped back to the size the
    distribution had where it was still linear (3 sigma = LINEAR_FRACTION of the distance to the nearest written
    state), and from there the deterministic flow of the full drift, with the bias, to the end of the sweep
    assigns each to a stable state. Where the spread reaches that size before the crossing, the noise is too strong
    for the deviation to stay linear through the crossing (recorded as linear_size_before_crossing); the certificate
    then decides whether the law still holds."""
    dt = dt or 0.5 * an.safe_dt(real)
    n = q0.size
    ell = min(float(np.linalg.norm(s - q0)) for s in stable) if stable else an.base_length(real)
    results = []
    for target, h in zip(targets, strengths):
        # reference: the unwritten branch itself, followed by Newton through the crossing (after it, an unstable
        # fixed point). Where the branch moves along the sweep, the state lags behind it; the lag, -dq*/dt, enters
        # the mean as a deterministic forcing (a sweep that itself favours one state). Where Newton fails, the
        # reference is carried by the drift.
        q_ref = q0.copy()
        m = np.zeros(n); C = np.zeros((n, n))
        t, t_star, grown, ended_by = 0.0, None, 0.0, "end of the sweep"
        hist_t, hist_spread, hist_q, hist_g = [], [], [], []
        steps = int(np.ceil(sweep.T / dt))
        for i in range(steps):
            over = sweep.over(t)
            J = an.jacobian(real, q_ref, over)
            if t > sweep.t_cross:
                grown += dt * max(float(np.max(np.linalg.eigvals(L.T @ J @ R).real)), 0.0)
            b = real.write_force(q_ref, target, h)
            q_new, ok = an.newton(real, q_ref, sweep.over(t + dt))
            q_new = np.asarray(q_new, float)
            if not ok or not np.all(np.isfinite(q_new)) or np.linalg.norm(q_new - q_ref) > 0.1 * an.base_length(real):
                q_new = _rk4_step(lambda q, s: real.F(q, **sweep.over(s)), q_ref, t, dt)
                lag = np.zeros(n)
            else:
                lag = -(q_new - q_ref) / dt
            m = m + dt * (J @ m + b + lag)
            C = C + dt * (J @ C + C @ J.T + 2 * D * np.eye(n))
            q_ref = q_new
            t += dt
            Cc = L.T @ C @ L
            spread = 3 * sqrt(max(float(np.max(np.linalg.eigvalsh(0.5 * (Cc + Cc.T)))), 0.0))
            hist_t.append(t); hist_spread.append(spread + float(np.linalg.norm(L.T @ m))); hist_q.append(q_ref.copy())
            hist_g.append(grown if t > sweep.t_cross else 0.0)
            if t > sweep.t_cross and grown >= GROWTH_STOP:
                t_star, ended_by = t, "growth"
                break
        if t_star is None:
            t_star = t
        mc = L.T @ m
        Cc = 0.5 * (L.T @ C @ L + (L.T @ C @ L).T)
        w, V = np.linalg.eigh(Cc)
        z = mc[None] + rng.standard_normal((n_samples, mc.size)) @ (V * np.sqrt(np.clip(w, 0, None))).T
        size = LINEAR_FRACTION * ell
        reach = next((hist_t[i] for i, sp in enumerate(hist_spread) if sp >= size), t_star)
        # the decision statistics are those at t_star (they contain the bias and the noise of the whole linear stage);
        # they are placed where the deviation has just started to grow (GROWTH_MAP e-folds after the crossing), at its
        # size there capped at the size of the linear stage, and carried by the drift alone: the map decides only to
        # which state each direction leads under the full nonlinear drift
        i_map = next((i for i, g in enumerate(hist_g) if g >= GROWTH_MAP), len(hist_t) - 1)
        t_map = hist_t[i_map]
        size_map = min(hist_spread[i_map], size)
        scale = size_map / max(hist_spread[-1], 1e-300)
        q = hist_q[i_map][None] + (scale * z) @ R.T
        def fmap(q, t):
            return real.F(q, **sweep.over(t))
        tt = t_map
        while tt < sweep.T - 1e-12:
            h_step = min(dt, sweep.T - tt)
            q = real.carrier.wrap(_rk4_step(fmap, q, tt, h_step))
            tt += h_step
        results.append({"t_star": t_star, "t_cross": sweep.t_cross, "t_map": t_map, "ended_by": ended_by,
                        "growth": grown, "linear_size_before_crossing": bool(reach < sweep.t_cross),
                        "mean_critical": mc.tolist(), "cov_critical": Cc.tolist(), "final": q})
    return results


def reduced_law(real: Realization, sweep: Sweep, q0: np.ndarray, targets, strengths, D: float,
                stable: List[np.ndarray], n_traj: int, rng, n_grid: int = 49) -> List[List[float]]:
    """The one-coordinate law along the moving branch, for a write point whose completion keeps a second parameter.

    The unwritten branch q*(t) is followed through the crossing; at n_grid times its reduction to the slow direction
    (analysis.normal_form: s' = a1 s + a2 s^2 + a3 s^3, left vector w) and the projections of the lag
    beta = -w.dq*/dt, of the bias and of the noise D_s = D w.w give
        ds/dt = a1(t) s + a2(t) s^2 + a3(t) s^3 + beta(t) + h_s(t) + sqrt(2 D_s(t)) xi,
    which is solved from s = 0. The final s is assigned to the stable state whose coordinate along w is nearest. The
    lag of the state behind the moving branch, which a linear stage overestimates near the crossing, is limited here by
    the quadratic and cubic terms."""
    ts_all = np.linspace(0.0, sweep.T, n_grid)
    qs, coef, real_from = [], [], None
    q = q0.copy()
    for i, t in enumerate(ts_all):
        qn, ok = an.newton(real, q, sweep.over(t))
        q = np.asarray(qn, float) if ok else q
        nf = an.normal_form(real, q, sweep.param, sweep.value(t))
        if nf.get("hopf") or "mode" not in nf:
            # far from the crossing the slowest relaxation can be a complex pair (a spiral): the one-coordinate
            # reduction exists only from the point where the slow direction becomes a real line
            if real_from is not None and t > sweep.t_cross:
                raise ValueError("the slow direction along the sweep is not a real line after the crossing")
            qs, coef, real_from = [], [], None
            continue
        if real_from is None:
            real_from = i
        qs.append(q.copy())
        coef.append((nf["a1"], nf["a2"], nf["a3"], np.asarray(nf["mode"]), np.asarray(nf["left_mode"])))
    if real_from is None or ts_all[real_from] > sweep.t_cross:
        raise ValueError("no real slow direction before the crossing")
    ts = ts_all[real_from:]
    qs = np.array(qs)
    dq = np.gradient(qs, ts, axis=0)
    a1 = np.array([c[0] for c in coef]); a2 = np.array([c[1] for c in coef]); a3 = np.array([c[2] for c in coef])
    W = np.array([c[4] for c in coef])
    # keep the orientation of the slow direction continuous along the sweep
    V = np.array([c[3] for c in coef])
    for i in range(1, len(ts)):
        if V[i] @ V[i - 1] < 0:
            V[i] = -V[i]; W[i] = -W[i]; a2[i] = -a2[i]
    beta = -np.einsum("ij,ij->i", W, dq)
    Ds = D * np.einsum("ij,ij->i", W, W)
    s_states = [float(W[-1] @ real.carrier.diff(qs[-1], st)) for st in stable]
    dt = 0.25 * an.safe_dt(real)
    out = []
    for target, h in zip(targets, strengths):
        hs = np.array([float(W[i] @ real.write_force(qs[i], target, h)) for i in range(len(ts))])
        # start where the reduction begins, from the stationary distribution of the relaxing coordinate there
        lam0 = min(float(a1[0]), -1e-9)
        x = (hs[0] + beta[0]) / (-lam0) + np.sqrt(Ds[0] / (-lam0)) * rng.standard_normal(n_traj)
        t = float(ts[0])
        while t < sweep.T - 1e-12:
            st = min(dt, sweep.T - t)
            c1, c2, c3, b, hh, d = (np.interp(t, ts, y) for y in (a1, a2, a3, beta, hs, Ds))
            x = x + st * (c1 * x + c2 * x * x + c3 * x ** 3 + b + hh) + np.sqrt(2 * d * st) * rng.standard_normal(n_traj)
            x = np.clip(x, -1e3, 1e3)
            t += st
        lab = np.argmin(np.abs(x[:, None] - np.array(s_states)[None, :]), axis=1)
        out.append([float(np.mean(lab == j)) for j in range(len(stable))])
    return out


def certificate(real: Realization, sweep: Sweep, q0: np.ndarray, targets, strengths, D: float, n_traj: int,
                rng) -> List[np.ndarray]:
    """Full stochastic dynamics over the sweep. Where more than 1% of the trajectories diverge (strong noise against a
    stiff drift, as a Kerr frequency shift at large amplitude), the run is repeated at a quarter of the time step;
    trajectories that still diverge are returned as non-finite and counted separately."""
    out = []
    dt0 = an.safe_dt(real)
    for target, h in zip(targets, strengths):
        for dt in (dt0, dt0 / 4):
            q = np.repeat(q0[None], n_traj, axis=0)
            q, _, _ = an.integrate(real, q, sweep.T, D, rng, over={}, write=(target, h),
                                   control_path=lambda s: sweep.value(s), dt=dt)
            if np.mean(~np.all(np.isfinite(q), axis=1)) <= 0.01:
                break
        out.append(q)
    return out


def assign(real: Realization, finals: np.ndarray, stable: List[np.ndarray]) -> np.ndarray:
    """Index of the nearest stable state; -1 for a final state that is not finite (a diverged trajectory)."""
    S = np.array(stable)
    ok = np.all(np.isfinite(finals), axis=1)
    lab = np.full(finals.shape[0], -1, dtype=int)
    if np.any(ok):
        d = np.linalg.norm(real.carrier.diff(finals[ok][:, None, :], S[None, :, :]), axis=-1)
        lab[ok] = np.argmin(d, axis=1)
    return lab


def _probs(labels: np.ndarray, m: int) -> List[float]:
    """Probabilities of the stable states among the trajectories that stay finite."""
    ok = labels >= 0
    n = max(int(np.sum(ok)), 1)
    return [float(np.sum(labels[ok] == j)) / n for j in range(m)]


def hold_check(real: Realization, sweep: "Sweep", finals: List[np.ndarray], targets, strengths, D: float, rng,
               stable: List[np.ndarray], hold: float) -> Dict[str, object]:
    """Continue the certificate's final states at the end value of the control, with the same noise and bias, for
    `hold` time units. The fraction of trajectories that change state says whether the written states retain the choice
    made during the sweep (a write, and the deviation from the law belongs to the nonlinear stage) or exchange it with
    the noise (an equilibrium write: the full dynamics gives the steady-state selection under the bias, as the few-photon
    oscillator does on the open quantum carrier)."""
    over = sweep.over(sweep.T)
    dt0 = an.safe_dt(real)
    switched, after = [], []
    for qf, target, h in zip(finals, targets, strengths):
        ok = np.all(np.isfinite(qf), axis=1)
        if not np.any(ok):
            switched.append(0.0)
            after.append([float("nan")] * len(stable))
            continue
        q = qf[ok]
        lab0 = assign(real, q, stable)
        q1, _, _ = an.integrate(real, q, hold, D, rng, over=over, write=(target, h), dt=dt0)
        lab1 = assign(real, q1, stable)
        fin = (lab1 >= 0) & (lab0 >= 0)
        switched.append(float(np.mean((lab0 != lab1)[fin])) if np.any(fin) else 0.0)
        after.append(_probs(lab1, len(stable)))
    return {"hold_time": float(hold), "switching_fraction": float(np.mean(switched)) if switched else 0.0,
            "switching_by_row": switched, "after_hold": after}


SWITCHING_EQUILIBRIUM = 0.1   # fraction of the certificate's trajectories that change state during a hold as long as the sweep


# ------------------------------------------------------------------------------------------------ one write point
def law_at(real: Realization, nf: Dict, rng, rate: float = 0.02, noise: Optional[float] = None, floors=(1.0,),
           n_samples: int = 4000, n_traj: int = 2000, certify: bool = True, gamma: float = 0.05,
           reduce_one: bool = False) -> Dict[str, object]:
    """The write law at one write point of `real` (whose normal form is nf), for biases along the critical space.

    The noise is `noise` if given; otherwise it is set, as in construct.swept_write_check, by
    Gamma = D_s |a3| / rate (the noise against the growth of the barrier behind the choice) where the cubic coefficient
    is resolved, and is the material's own noise where it is not (a vanishing cubic)."""
    k = max(1, int(nf.get("critical_dimension", 1)))
    R, L = critical_space(real, nf, k)
    sweep = Sweep(real, nf, rate)
    qc = np.asarray(nf["state"], float)
    q0 = _relaxed(real, qc, sweep.over(0.0))
    stable = stable_states(real, rng, sweep.over(sweep.T), seeds=swept_seeds(real, sweep, q0, R, sweep.T))
    if len(stable) < 2:
        # the sweep leaves the range in which a choice exists (a cusp sweep that crosses the far fold): it ends at the
        # latest time after the crossing at which two or more states are stable
        for t_end in np.linspace(sweep.T, sweep.t_cross, 25)[1:-1]:
            st = stable_states(real, rng, sweep.over(float(t_end)),
                               seeds=swept_seeds(real, sweep, q0, R, float(t_end)))
            if len(st) >= 2:
                sweep.T, stable = float(t_end), st
                break
    if len(stable) < 2:
        return {"status": "no choice at the end of the sweep", "stable_states": [s.tolist() for s in stable]}
    # a marginal direction at every written state (a Goldstone mode of a continuous symmetry, as the phase of a
    # non-degenerate parametric oscillator): the written states form a continuous family, and the write chooses a
    # point on it rather than one of separate states
    marg = []
    for st in stable:
        ev = np.linalg.eigvals(an.jacobian(real, st, sweep.over(sweep.T)))
        marg.append(float(np.min(np.abs(ev.real))) <= 1e-6 * max(float(np.max(np.abs(ev))), 1e-12))
    if all(marg):
        return {"status": "continuous family of written states (a marginal direction at every one): the write sets a "
                          "phase on the family, not a choice between states; its law is that of a phase (Adler)",
                "stable_states": [s.tolist() for s in stable[:8]], "completion_needed": "phase"}
    dfac = float(np.mean(np.diag(L.T @ L)))
    a3 = abs(float(nf.get("a3", 0.0) or 0.0))
    if noise is not None:
        D, noise_rule = float(noise), "given"
    elif a3 > 0 and "degenerate cubic" not in nf.get("kind", "") and "line of states" not in nf.get("kind", ""):
        D, noise_rule = gamma * rate / (a3 * dfac), f"Gamma = {gamma}"
    else:
        D, noise_rule = float(real.noise), "material noise"
    D_s = D * dfac
    floor = sqrt(D_s) * rate ** 0.25 / pi ** 0.25          # the bias at which the pitchfork write follows it (z = 1)
    far = 1e3 * an.base_length(real)
    # directions in the critical space: +-each axis for k = 1; 2k directions on the unit circle of the first two axes
    # and +-each further axis otherwise
    dirs: List[np.ndarray] = []
    if k == 1:
        dirs = [np.array([1.0]), np.array([-1.0])]
    else:
        for j in range(8):
            e = np.zeros(k); e[0], e[1] = np.cos(j * pi / 4), np.sin(j * pi / 4)
            dirs.append(e)
        for a in range(2, k):
            for sg in (1, -1):
                e = np.zeros(k); e[a] = sg; dirs.append(e)
    targets, strengths, keys = [], [], []
    for c in floors:
        for e in dirs:
            u = R @ e
            u = u / np.linalg.norm(u)
            tgt = qc + far * u
            per = float(np.linalg.norm(L.T @ real.write_force(q0, tgt, 1.0)))
            targets.append(tgt); strengths.append(c * floor / max(per, 1e-12)); keys.append((c, e))
    lin = linear_stage_law(real, sweep, q0, R, L, targets, strengths, D, stable, n_samples, rng)
    reduced = None
    if k == 1 and reduce_one:
        try:
            reduced = reduced_law(real, sweep, q0, targets, strengths, D, stable, n_samples, rng)
        except (ValueError, np.linalg.LinAlgError):
            reduced = None
    rows = []
    for (c, e), r in zip(keys, lin):
        rows.append({"bias_floors": c, "direction": e.tolist(),
                     "linear_stage_law": _probs(assign(real, r["final"], stable), len(stable)), "t_star": r["t_star"],
                     "t_cross": r["t_cross"], "t_map": r["t_map"], "linear_stage_ended_by": r["ended_by"],
                     "growth": r["growth"], "linear_size_before_crossing": r["linear_size_before_crossing"],
                     "mean_critical": r["mean_critical"]})
    for i, row in enumerate(rows):
        if reduced is not None:
            row["reduced_law"] = reduced[i]
        # the law reported: the reduced one-coordinate law where the completion keeps a second parameter (a sweep
        # that itself favours one state), the linear-stage law otherwise
        row["law"] = row["reduced_law"] if reduced is not None else row["linear_stage_law"]
    sims: List[np.ndarray] = []
    if certify:
        sims = certificate(real, sweep, q0, targets, strengths, D, n_traj, rng)
        for row, qf in zip(rows, sims):
            lab = assign(real, qf, stable)
            row["certificate"] = _probs(lab, len(stable))
            row["certificate_diverged"] = float(np.mean(lab < 0))
            err = [sqrt(max(p * (1 - p), 1.0 / n_traj) / n_traj) for p in row["certificate"]]
            row["deviation_sigma"] = float(max(abs(a - b) / e for a, b, e in zip(row["law"], row["certificate"], err)))
            row["deviation"] = float(max(abs(a - b) for a, b in zip(row["law"], row["certificate"])))
    # the states written by +-axis biases, in critical coordinates (for k = 1: which stable state is +mode)
    coords = [(L.T @ real.carrier.diff(qc, s)).tolist() for s in stable]          # diff(a, b) = b - a
    out = {"status": "derived", "law_kind": "reduced (one coordinate along the moving branch)" if reduced is not None
           else "linear stage and outcome map", "critical_dimension": k, "rate": rate, "sweep_end": sweep.T, "noise": D,
           "noise_rule": noise_rule,
           "noise_along_critical": D_s,
           "bias_floor": floor, "stable_states": [s.tolist() for s in stable], "stable_critical_coordinates": coords,
           "rows": rows}
    if k == 1:
        out.update(_one_direction(rows, coords, floor))
    if certify:
        bad = [r for r in rows if r["deviation_sigma"] > DEVIATION_SIGMA and r["deviation"] > DEVIATION_FLOOR]
        out["certified"] = not bad
        if bad:
            # does the end of the sweep retain the choice? continue the certificate at the end value of the control
            hold = hold_check(real, sweep, sims, targets, strengths, D, rng, stable, sweep.T)
            out["hold"] = hold
            dev = f"{max(r['deviation'] for r in bad):.3f} ({max(r['deviation_sigma'] for r in bad):.1f} errors)"
            if hold["switching_fraction"] > SWITCHING_EQUILIBRIUM:
                out["write_class"] = "equilibrium write"
                out["obstruction_in_nonlinear_stage"] = (
                    f"the law differs from the full dynamics by up to {dev}, and the written states at the end of the "
                    f"sweep do not retain the choice: {100 * hold['switching_fraction']:.0f} % of the trajectories change "
                    f"state during a hold of {sweep.T:.3g} time units at the end value of the control. The full dynamics "
                    "gives the steady-state selection under the bias, an equilibrium write, not the linear-stage law; "
                    "a retained write needs deeper states at the end of the sweep (a longer sweep, or less noise)")
            else:
                out["write_class"] = "nonlinear stage"
                out["obstruction_in_nonlinear_stage"] = (
                    f"the law differs from the full dynamics by up to {dev}: the choice is not complete at the end of "
                    "the linear stage (noise during the nonlinear stage, or a further mode brought into play), and the "
                    f"written states retain it ({100 * hold['switching_fraction']:.0f} % change state during a hold of "
                    f"{sweep.T:.3g} time units)")
        else:
            out["write_class"] = "linear-stage write"
    return out


def _one_direction(rows: List[Dict], coords: List[List[float]], floor: float) -> Dict[str, object]:
    """Response constant and offset of a write along one critical direction, from the law and the certificate."""
    plus = int(np.argmax([c[0] for c in coords])); minus = int(np.argmin([c[0] for c in coords]))
    out = {}
    for src in ("law", "certificate"):
        if src not in rows[0]:
            continue
        for c in sorted({r["bias_floors"] for r in rows}):
            rp = next(r for r in rows if r["bias_floors"] == c and r["direction"][0] > 0)
            rm = next(r for r in rows if r["bias_floors"] == c and r["direction"][0] < 0)
            pp = min(max(rp[src][plus], 1e-4), 1 - 1e-4)
            pm = min(max(rm[src][minus], 1e-4), 1 - 1e-4)
            zp, zm = ND.inv_cdf(pp), ND.inv_cdf(pm)
            out.setdefault(src, []).append({"bias_floors": c, "P_toward_plus": pp, "P_toward_minus": pm,
                                            "response_constant": 0.5 * (zp + zm) / c * pi ** 0.25,
                                            "offset": 0.5 * (zp - zm)})
    return {"one_direction": out}


# ------------------------------------------------------------------------------------------------ the realization
def zero_offset_curvature(real: Realization, cusp: Dict, rng, rate: float = 0.02, gamma: float = 0.05,
                          curvatures: Sequence[float] = CURVATURES, n_samples: int = 2000, n_traj: int = 2000
                          ) -> Dict[str, object]:
    """The curvature of the sweep through a cusp at which its intrinsic offset vanishes.

    A straight sweep through a cusp is tangent to the line on which the two written states are mirror images and leaves
    it at second order; the noise, acting on the quadratic term that grows away from the cusp, adds a bias of its own.
    Both appear as an offset of the write: without bias it does not choose either state with probability one half. The
    sweep is curved, p = cusp + t d + c2 t^2 n with n the direction normal to d in the plane of the two parameters,
    and the offset is read from the full stochastic dynamics (the certificate of law_at, a bias of one noise floor
    toward either state) at each curvature. A straight line fitted to the offset against c2 gives the curvature c2*
    at which it vanishes, and a last run at c2* checks it. `cusp` is a successful result of construct.cancel_asymmetry
    together with the name of its second parameter under "second"."""
    second = cusp["second"]
    origin = {real.control: float(cusp["control"]), second: float(cusp["value"])}
    q = np.asarray(cusp["state"], float)

    def offset_at(c2: float) -> Dict[str, object]:
        path = construct.along_curved_path(real, origin, cusp["direction"], c2)
        nfp = construct.sweep_normal_form(path, q)
        row: Dict[str, object] = {"curvature": float(c2), "kind_along_sweep": nfp["kind"]}
        if "pitchfork" not in str(nfp["kind"]).lower():
            row["status"] = "not a pitchfork along the curved sweep"
            return row
        L = law_at(path, nfp, rng, rate, None, (1.0,), n_samples, n_traj, True, gamma, reduce_one=False)
        cert = ((L.get("one_direction") or {}).get("certificate") or [None])[0]
        if not cert:
            row["status"] = str(L.get("status", "no certificate"))
            return row
        p_plus, p_minus = cert["P_toward_plus"], cert["P_toward_minus"]
        err = 0.5 * sum(sqrt(max(x * (1 - x), 1.0 / n_traj) / n_traj) / max(ND.pdf(ND.inv_cdf(x)), 1e-6)
                        for x in (p_plus, p_minus))
        row.update(status="measured", offset=float(cert["offset"]), offset_error=float(err),
                   response_constant=float(cert["response_constant"]), noise_ratio=float(gamma))
        return row

    rows = [offset_at(c) for c in curvatures]
    ok = [r for r in rows if r.get("status") == "measured"]
    out: Dict[str, object] = {"second": second, "rate": rate, "noise_ratio": gamma, "rows": rows,
                              "curvature_units": "the sweep parameter t of p = cusp + t d + c2 t^2 n, with d the "
                                                 "direction of cancel_asymmetry and n normal to it"}
    if len(ok) < 2:
        out["status"] = "fewer than two curvatures with a measured offset"
        return out
    x = np.array([r["curvature"] for r in ok]); y = np.array([r["offset"] for r in ok])
    w = 1.0 / np.maximum(np.array([r["offset_error"] for r in ok]), 1e-3)
    slope, intercept = np.polyfit(x, y, 1, w=w)
    out.update(offset_slope=float(slope), offset_straight=float(intercept))
    if abs(slope) < 1e-9:
        out["status"] = "the offset does not depend on the curvature"
        return out
    c_star = -intercept / slope
    span = float(np.max(x) - np.min(x))
    if not (np.min(x) - span <= c_star <= np.max(x) + span):
        out.update(status="the zero lies far outside the curvatures tried", zero_offset_curvature=float(c_star))
        return out
    check = offset_at(float(c_star))
    out.update(status="found", zero_offset_curvature=float(c_star), check=check)
    if check.get("status") == "measured":
        out["cancelled"] = bool(abs(check["offset"]) <= max(3 * check["offset_error"], 0.1))
    return out


def complete(real: Realization, rng, rate: float = 0.02, noise: Optional[float] = None, floors=(1.0,),
             n_samples: int = 4000, n_traj: int = 2000, certify: bool = True, gamma: float = 0.05,
             curvature: bool = False) -> Dict[str, object]:
    """Every write point of `real` along its control: the obstruction, what its completion keeps, and the law of the
    completed model (with its certificate)."""
    if not real.control:
        return {"write_points": [], "note": "no control parameter"}
    events = an.locate_writes(real, rng, n_starts=24 if real.carrier.dim <= 4 else 12)
    points = []
    for ev in events:
        nf = an.normal_form(real, ev["q"], real.control, ev["v"])
        ob = construct.obstruction(nf)
        kind = completion_kind(nf, ob)
        entry = {"write_point": {real.control: float(ev["v"]), "state": list(map(float, ev["q"]))},
                 "kind": nf["kind"], "obstruction": ob.get("terms", {}), "completion": kind}
        if nf.get("leading_term"):
            entry["leading_term"] = nf["leading_term"]
        try:
            if kind in ("none", "critical space", "leading order", "distant branch"):
                entry["law"] = law_at(real, nf, rng, rate, noise, floors, n_samples, n_traj, certify, gamma)
            elif kind == "second parameter":
                entry["unfoldings"] = []
                for second in real.material:
                    c = construct.cancel_asymmetry(real, nf, second)
                    u = {"second": second, "success": bool(c.get("success"))}
                    if c.get("success"):
                        path = construct.along_path(real, {real.control: c["control"], second: c["value"]},
                                                    c["direction"])
                        # the sweep is read on a window narrowed until the quartic term of the drift no longer appears
                        # as a quadratic one (construct.sweep_normal_form): on the default window the cusp of a driven
                        # Kerr oscillator reads as transcritical
                        nfp = construct.sweep_normal_form(path, np.asarray(c["state"]))
                        u.update(cusp={real.control: c["control"], second: c["value"]}, direction=c["direction"],
                                 kind_along_sweep=nfp["kind"], asymmetry_at_fold=float(nf.get("asymmetry", 0.0)),
                                 asymmetry_along_sweep=float(nfp.get("asymmetry", 0.0)),
                                 window_halvings=int(nfp.get("window_halvings", 0)))
                        if "pitchfork" not in str(nfp["kind"]).lower():
                            # the path through the cusp does not cross a symmetric point: the second parameter has
                            # not cancelled the asymmetry along the sweep (a transcritical point, two states exchanging
                            # stability, or a fold), and the linear-stage law does not apply
                            u["law"] = {"status": f"no law: the sweep through the cusp found along {second} is "
                                                  f"{str(nfp['kind']).split(':')[0]}, not a pitchfork; the asymmetry is "
                                                  "not cancelled on the path"}
                        else:
                            # the one-coordinate law needs a saturating cubic along the sweep; a subcritical sweep
                            # through the cusp keeps the linear-stage law
                            u["law"] = law_at(path, nfp, rng, rate, noise, floors, n_samples, n_traj, certify, gamma,
                                              reduce_one=float(nfp.get("a3", 0.0)) < 0)
                            if curvature and noise is None and u["law"].get("status") == "derived":
                                # the curvature of the sweep that cancels its intrinsic offset
                                u["zero_offset_curvature"] = zero_offset_curvature(
                                    real, {**c, "second": second}, rng, rate=rate, gamma=gamma,
                                    n_samples=min(n_samples, 2000), n_traj=n_traj)
                    entry["unfoldings"].append(u)
            else:
                entry["law"] = {"status": "phase: derived by the phase-locking target (memory codiscover --target "
                                          "phase-locking)"}
        except (ValueError, np.linalg.LinAlgError) as e:
            entry["law"] = {"status": f"not derived: {e}"}
        points.append(entry)
    return {"write_points": points}
