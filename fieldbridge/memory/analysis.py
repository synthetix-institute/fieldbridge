"""Measure how a realization writes, retains, is observed and rewritten, and where its writes happen.

  stored_states      stable states reached by relaxation from random starts, with basin counts
  kappa_spectrum     kappa = -eig(Jacobian): > 0 relaxing, ~ 0 flat (symmetry), < 0 unstable (writes)
  loss_law           the law by which a stored state is lost, read from the spectra and the number of states
  neb_barrier        barrier between two stored states along a minimum-energy path, relaxed to a stated tolerance;
                     for states half a period apart on a torus, the lowest over the straight paths of lowest
                     energy, with the paths relaxed and their barriers (gradient systems)
  basin_balls        a ball around each stored state that lies in its basin (checked by the noise-free flow)
  basin_of           the stored state whose basin contains a configuration (noise-free flow into those balls)
  hold_time          Monte Carlo time for a stored state to lose its identity at noise D: the trajectories that
                     have left its basin
  write_test         bounded write toward a target at fixed control or while sweeping the control
                     from the end of its range where the barriers are low; accuracy after release (the fraction
                     in the basin of the target) and the rewriting time
  fixed_points       zeros of the drift (Newton from random starts), stable and unstable, with kappa spectra
  track              continuation of a stable state along a parameter until it loses stability or vanishes
  locate_writes      parameter values at which stored states lose stability or appear (write points)
  series_reduction   drift on the centre manifold of one critical direction, closed order by order at the write point
  normal_form        drift along the critical direction at a write point, reduced onto the slow manifold, and the
                     number of directions that become critical there together
  with_write_field   the realization plus a bounded write toward a target, with the field WRITE_FIELD as parameter
  memory_kernel      for a linear realization with an observed block, K(t) = B e^{Dt} C (Mori-Zwanzig)
"""
from __future__ import annotations

from dataclasses import replace
from typing import Dict, List, Optional, Sequence, Tuple

from math import sqrt
import numpy as np

from .identity import Realization
from .series import NotAnalytic, Series, coefficients as series_coefficients

P_STAR = 0.9
# the parameter of the added write field; reserved in specifications (spec.RESERVED), so that it cannot replace a
# parameter of the material such as an applied field h
WRITE_FIELD = "h_write"
# resolution of the continuation in units of the scanned range: locate_writes takes write points closer than this for
# one, and normal_form counts every eigenvalue that crosses zero within it of a write point as critical there
WRITE_RESOLUTION = 2e-3
SERIES_ORDER = 7  # order of the reduction at a write point; a Kerr parametric oscillator at zero detuning needs 5


def orientation_index(v: np.ndarray, rtol: float = 1e-4) -> int:
    """The component that fixes the sign of a direction: the first one within rtol of the largest in magnitude. The
    sign of an eigenvector is arbitrary and differs between LAPACK builds; making this component positive gives every
    machine the same direction. A symmetric direction such as (a, -a) has components equal up to rounding, so the
    plain largest one would be chosen by rounding: the first within the tolerance is chosen instead."""
    a = np.abs(np.asarray(v, float))
    return int(np.flatnonzero(a >= (1.0 - rtol) * a.max())[0])


def base_length(real: Realization) -> float:
    return real.carrier.period / (2 * np.pi) if real.carrier.kind == "torus" else real.carrier.scale


def _merge_tol(real: Realization) -> float:
    return 0.05 * base_length(real) * np.sqrt(real.carrier.dim)


def _canon(real: Realization, q: np.ndarray) -> np.ndarray:
    return real.canonical(q) if real.canonical is not None else q


# ------------------------------------------------------------------------------------------------ dynamics
def _explicit_steps(real: Realization) -> Tuple[float, float]:
    """The step of a fourth-order Runge-Kutta integration and the step of an explicit Euler integration, from the
    eigenvalues of the Jacobian at sampled states and at the ends of the control range (see rate_dt and safe_dt)."""
    rng = np.random.default_rng(12345)
    rate = 0.0
    focus = np.inf
    overs = [{}]
    if real.control and real.control_range:
        overs += [{real.control: real.control_range[0]}, {real.control: real.control_range[1]}]
    for over in overs:
        for q in real.carrier.sample(rng, 6):
            J = jacobian(real, q, over)
            if np.all(np.isfinite(J)):
                ev = np.linalg.eigvals(J)
                rate = max(rate, float(np.max(np.abs(ev))))
                stable = ev[ev.real < -1e-3]
                if stable.size:
                    focus = min(focus, float(np.min(np.abs(stable.real) / np.abs(stable) ** 2)))
    step = float(real.dt if rate <= 0 else min(real.dt, 0.3 / rate))
    return step, min(step, focus)


def rate_dt(real: Realization) -> float:
    """The declared time step, reduced to 0.3 / max |eigenvalue| over sampled states and the control range: the step of
    a fourth-order Runge-Kutta integration (correction 35). Runge-Kutta multiplies an eigenvalue l by R(l dt), with
    |R(i y)|^2 = 1 - y^6 / 72 + y^8 / 576 on the imaginary axis, so at |l| dt = 0.3 it keeps a stable focus stable
    however slowly it decays, and needs none of the bound that safe_dt adds for Euler. Under correction 34 the attractor
    step took safe_dt, and a fixed number of its steps covered a time that fell with that bound, by a factor of 7e-5 in
    one corpus model."""
    cached = getattr(real, "_rate_dt", None)
    if cached is not None:
        return cached
    step, euler = _explicit_steps(real)
    real._rate_dt = step
    if getattr(real, "_safe_dt", None) is None:
        real._safe_dt = euler
    return step


def safe_dt(real: Realization) -> float:
    """The declared time step, reduced where the fastest local rate would make an explicit Euler step unstable
    (dt <= 0.3 / max |eigenvalue| over sampled states and the control range), and where a stable focus would not stay
    stable under it (correction 34): Euler multiplies an eigenvalue l = a + i b by |1 + l dt|, which stays below 1 only
    for dt < 2 |a| / |l|^2. A focus that turns much faster than it decays (b >> |a|) needs a step far below 0.3 / |l|,
    and at that step the explicit prefix of the relaxation made its starts grow by |1 + i b dt| each step until they
    were no longer finite. The bound dt <= |a| / |l|^2, half the limit, is taken over the eigenvalues with a < -1e-3 (the
    card's resolution of rates); for a real eigenvalue it is 1 / |l|, never below 0.3 / |l|, so it binds only on
    spirals. A Runge-Kutta integration takes rate_dt instead (correction 35)."""
    cached = getattr(real, "_safe_dt", None)
    if cached is not None:
        return cached
    step, euler = _explicit_steps(real)
    real._safe_dt = euler
    if getattr(real, "_rate_dt", None) is None:
        real._rate_dt = step
    return euler


def refine_steps(real: Realization, factor: float) -> None:
    """Divide both explicit steps of a realization by factor: a run at a finer numerical resolution (correction 35;
    before it the engine divided safe_dt alone)."""
    rate, euler = rate_dt(real), safe_dt(real)
    real._rate_dt, real._safe_dt = rate / factor, euler / factor


# explicit steps a relaxation or a noiseless trajectory may take at the step of the fastest rate (about the engine's time
# limit of 900 s at 100 microseconds a step); a longer one is integrated by an implicit, adaptive method instead, after an
# explicit prefix within which most stiff relaxations that settle at all have settled (a fast rate settles in a few
# hundred steps)
EXPLICIT_BUDGET = 10_000_000
EXPLICIT_PREFIX = 100_000
# calls of the drift, those of the settling event included, that one implicit run of a relaxation may make (correction
# 30): the stiff relaxations of the development data need at most 484 calls of the drift a run, while a run that creeps
# along a singular line of the drift (an inspiral at |e| -> 1, a < 0) made 1.5e5 calls in 20 s without ending
IMPLICIT_CALLS = 20_000
# explicit steps a relaxation may take before it goes to the implicit method (correction 31). The card relaxes 48 starts
# at once, and more than once; at 1e-4 s a step of 48 starts (2 to 10 variables on the cluster), 3e6 steps take about
# 300 s, a third of the engine's time limit, where EXPLICIT_BUDGET took about 1000 s. A budget in steps, not in measured
# time, keeps the records independent of the load of the machine. Noisy and noiseless trajectories keep EXPLICIT_BUDGET
RELAX_BUDGET = 3_000_000


class _Spent(Exception):
    """An implicit run of a relaxation that has made IMPLICIT_CALLS calls of the drift without settling."""


def _implicit_runs(real: Realization, q: np.ndarray, over: Dict, T: float, tol: Optional[float] = None,
                   t_eval: Optional[np.ndarray] = None, unfinished: Optional[List] = None):
    """Each state of q integrated over (0, T] by backward differentiation (BDF) with the central-difference Jacobian:
    the step is set by accuracy, not by the fastest rate, so a stiff drift (rates of 1e5 to 1e15 against slow ones) is
    integrated in a number of steps that does not grow with its stiffness. With tol, a run stops where its speed falls
    below tol / 2 (the relaxation's criterion); with t_eval, the states at those times are returned as well. The drift
    is evaluated at the wrapped state (an orthant reflects at zero). A run whose drift or Jacobian stops being finite
    ends there, with its state NaN, which every caller reads as not converged. A run of a relaxation (with tol) that has
    made IMPLICIT_CALLS calls of the drift without settling ends there in the same way: it has not settled (correction
    30). The noiseless trajectories of integrate (with t_eval) have no such budget. Given a list as unfinished, the runs
    of a relaxation that end without settling, at the end of their time or of their budget, are appended to it as
    (index, last state) (correction 34)."""
    from scipy.integrate import solve_ivp
    calls = [0]
    last = [None]

    def fun(t, y):
        if tol is not None:
            calls[0] += 1
            if calls[0] > IMPLICIT_CALLS:
                raise _Spent()
            last[0] = y
        return real.F(real.carrier.wrap(y)[None], **over)[0]

    def jac(t, y):
        return jacobian(real, real.carrier.wrap(y), over)

    events = None
    if tol is not None:
        def settled(t, y):
            return float(np.linalg.norm(fun(t, y))) - 0.5 * tol
        settled.terminal, settled.direction = True, -1.0
        events = settled
    atol = 1e-9 * base_length(real)
    ends = np.full_like(q, np.nan)
    paths = None if t_eval is None else np.full((len(t_eval),) + q.shape, np.nan)
    for i, y0 in enumerate(q):
        calls[0] = 0
        last[0] = None
        try:
            with np.errstate(all="ignore"):
                sol = solve_ivp(fun, (0.0, T), y0, method="BDF", jac=jac, rtol=1e-6, atol=atol, events=events,
                                t_eval=t_eval)
        except (ValueError, ZeroDivisionError, OverflowError, np.linalg.LinAlgError):
            continue  # the drift or its Jacobian is not finite along the run (an escape)
        except _Spent:
            if unfinished is not None and last[0] is not None and np.all(np.isfinite(last[0])):
                unfinished.append((i, real.carrier.wrap(np.array(last[0], float))))
            continue  # the budget of the run is spent: it has not settled
        if sol.status < 0 or not np.all(np.isfinite(sol.y[:, -1])):
            continue
        ends[i] = sol.y[:, -1]
        if unfinished is not None and tol is not None and not len(sol.t_events[0]):
            unfinished.append((i, real.carrier.wrap(sol.y[:, -1])))
        if paths is not None:
            paths[:sol.y.shape[1], i] = sol.y.T
    return real.carrier.wrap(ends), paths


def _zero_of_unfinished_run(real: Realization, start: np.ndarray, end: np.ndarray,
                            over: Dict) -> Tuple[Optional[np.ndarray], bool]:
    """The stable zero of the drift near the end of an implicit run of a relaxation that ended without reaching the
    speed criterion, and whether the run has settled at it; (None, False) where there is none (corrections 34 and 35).

    Newton from the end of the run must converge to a zero that has no growing direction, that attracts in at least one
    direction (a rate with Re below -1e-3, the resolution of rates), and that is not a center: a pair of rates that turn
    without decay (|Re| at most the resolution, Im above it), or a zero all of whose rates are flat, as the origin of
    x'' = -x^3, would let the numerical damping of BDF on a neutral orbit pass for an approach. Where the search capped
    at half the base length does not converge, Newton runs again without the cap, its steps halved until the drift
    decreases (correction 35). Such a zero is a stable state by its own rates, whether or not the run reaches it.

    The run has settled at it when it moved towards it (its end is closer to it than the start of the implicit run) and
    the zero does not turn faster than it decays: no stable rate has |Im| > Re (correction 35). On a zero that turns
    faster than it decays, BDF damps the turn as well as the decay. On the fast Bautin point of set 30 it carried inwards
    the runs that start beyond the unstable cycle, which move outwards to the stable cycle, and they were read as settled
    at the origin. There the endpoints of a run cannot tell an approach from that damping, and the run is settled only
    by the speed criterion."""
    with np.errstate(all="ignore"):
        z, ok = newton(real, end, over)
        if not ok:
            z, ok = newton(real, end, over, descent=True)
        if not ok or not np.all(np.isfinite(z)):
            return None, False
        k = kappa_spectrum(real, z, over)
    if not np.all(np.isfinite(k)):
        return None, False
    c = classify(k)
    if c["unstable"] > 0 or c["stable"] == 0 or np.any((np.abs(k.real) <= 1e-3) & (np.abs(k.imag) > 1e-3)):
        return None, False
    stable = k[k.real > 1e-3]
    turns = bool(np.any(np.abs(stable.imag) > stable.real))
    closer = real.carrier.distance(_canon(real, end), _canon(real, z)) < real.carrier.distance(_canon(real, start),
                                                                                                _canon(real, z))
    return z, bool(closer and not turns)


def relax(real: Realization, q0: np.ndarray, over: Optional[Dict] = None, max_time: float = 300.0,
          tol: float = 1e-5, dt: Optional[float] = None, zeros: Optional[List] = None):
    """Relax every start by explicit Euler for max_time at the step dt (by default the step of the fastest rate,
    safe_dt), stopping when every speed is below tol. When the whole time would take more than RELAX_BUDGET steps,
    Euler runs as before for at most EXPLICIT_PREFIX steps; if the starts have not settled by then, each is integrated
    over the time that remains by an implicit, adaptive method with the same speed criterion. A relaxation that settles
    within the prefix, or needs no more than the budget, is unchanged; the states and their stability do not depend on
    the integrator. The budget of a relaxation is RELAX_BUDGET (correction 31). Given a list as zeros, the stable zeros
    near the ends of implicit runs that have not settled at them are appended to it as (index, zero) (correction 35)."""
    over, dt = over or {}, dt or safe_dt(real)
    q = real.carrier.wrap(np.atleast_2d(np.array(q0, float)))
    n = int(max_time / dt)
    implicit = n > RELAX_BUDGET
    for k in range(min(n, EXPLICIT_PREFIX) if implicit else n):
        f = real.F(q, **over)
        q = real.carrier.wrap(q + dt * f)
        if k % 50 == 0 and np.all(np.linalg.norm(f, axis=-1) < tol):
            break
    else:
        if implicit:
            # a run that ends without settling, at the end of its time or of its budget, while it approaches a zero
            # without a growing direction has settled at that zero (correction 34, second part): a stable focus that
            # turns much faster than it decays meets the speed criterion only at a distance of order tol / (its rate of
            # turning), and the implicit runs of an epidemic at the scale of a population ended with the infection
            # still decaying. Near a zero that turns faster than it decays the run is not settled, and the zero is
            # kept in zeros (correction 35)
            q_mid, unfinished = q.copy(), []
            q, _ = _implicit_runs(real, q, over, max_time - EXPLICIT_PREFIX * dt, tol=tol, unfinished=unfinished)
            for i, end in unfinished:
                with np.errstate(all="ignore"):
                    if float(np.linalg.norm(real.F(end[None], **over)[0])) < 1e-4:
                        continue  # slow enough to count as settled by stored_states, as before
                z, settles = _zero_of_unfinished_run(real, q_mid[i], end, over)
                if settles:
                    q[i] = z
                elif z is not None and zeros is not None:
                    zeros.append((i, z))
            with np.errstate(all="ignore"):
                return q, np.linalg.norm(real.F(q, **over), axis=-1)
    return q, np.linalg.norm(real.F(q, **over), axis=-1)


def integrate(real: Realization, q: np.ndarray, T: float, D: float, rng, over: Optional[Dict] = None,
              write=None, control_path=None, record=None, n_record: int = 40, dt: Optional[float] = None):
    """Euler-Maruyama; write = (target, h); control_path(t) sets real.control; record(q) -> value."""
    free = dt is None
    implicit = free and D == 0.0 and write is None and control_path is None
    over, dt = dict(over or {}), dt or safe_dt(real)
    steps = max(1, int(round(T / dt)))
    every = max(1, steps // n_record)
    if implicit and steps > EXPLICIT_BUDGET:
        # a noiseless trajectory longer than the budget: the same samples, integrated by the implicit method of relax
        k = np.arange(1, steps // every + 1)
        if every * int(k[-1]) < 2 ** 62:
            times = k * every * dt
        else:
            # beyond a 64-bit integer the product of the counts would be an array of Python integers, and the grid an
            # array of objects, which the implicit runs cannot test for finite values: at a step of 2e-36, where the
            # fastest rate over the control range is 1e35, a trajectory of 200 time units has 1e38 steps (correction 33).
            # The counts are then multiplied in floating point
            times = k * float(every) * dt
        t_end = steps * dt
        grid = times if times[-1] >= t_end else np.append(times, t_end)
        q, paths = _implicit_runs(real, real.carrier.wrap(np.atleast_2d(np.array(q, float))), over, float(grid[-1]),
                                  t_eval=grid)
        out_v = [record(real.carrier.wrap(x)) for x in paths[:len(times)]] if record is not None else []
        return q, (times if record is not None else np.array([])), np.array(out_v)
    stiff = free and not implicit and steps > EXPLICIT_BUDGET and real.dt > dt
    if stiff:
        # a noisy, written or driven trajectory longer than the budget: the drift-implicit step at the declared step
        dt = float(real.dt)
        steps = max(1, int(round(T / dt)))
        every = max(1, steps // n_record)
    amp = np.sqrt(2.0 * D * dt)
    out_t, out_v = [], []
    for k in range(steps):
        if control_path is not None:
            over[real.control] = control_path(k * dt)
        f = real.F(q, **over)
        if write is not None:
            f = f + real.write_force(q, write[0], write[1])
        if stiff:
            q = real.carrier.wrap(q + _drift_implicit(real, q, over, dt, dt * f + amp * rng.standard_normal(q.shape)))
        else:
            q = real.carrier.wrap(q + dt * f + amp * rng.standard_normal(q.shape))
        if record is not None and (k + 1) % every == 0:
            out_t.append((k + 1) * dt)
            out_v.append(record(q))
    return q, np.array(out_t), np.array(out_v)


def _drift_implicit(real: Realization, q: np.ndarray, over: Dict, dt: float, kick: np.ndarray) -> np.ndarray:
    """The linearly implicit (drift-implicit) Euler-Maruyama increment (I - dt J(q))^-1 (dt f + sqrt(2 D dt) xi) for
    every row of q, with the Jacobian of the drift at q (correction 25 of the discovery engine). A rate lambda with
    lambda dt >> 1 relaxes within one step instead of making it unstable: the fast variable is slaved to the slow ones,
    and its fluctuations, of variance D / lambda, are averaged out. Rates with lambda dt << 1 are integrated as by
    Euler-Maruyama. A row whose Jacobian is not finite becomes NaN, as an escape."""
    J = _jacobians(real, q, over)
    n = q.shape[1]
    A = np.eye(n)[None] - dt * J
    out = np.full_like(kick, np.nan)
    ok = np.all(np.isfinite(A), axis=(1, 2)) & np.all(np.isfinite(kick), axis=1)
    if ok.any():
        out[ok] = np.linalg.solve(A[ok], kick[ok][..., None])[..., 0]
    return out


def _jacobians(real: Realization, Q: np.ndarray, over: Dict, eps: float = 1e-6) -> np.ndarray:
    """jacobian() at every row of Q in one evaluation of the drift, with the same arithmetic."""
    m, n = Q.shape
    E = eps * np.eye(n)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        Fp = real.F((Q[:, None, :] + E[None]).reshape(m * n, n), **over).reshape(m, n, n)
        Fm = real.F((Q[:, None, :] - E[None]).reshape(m * n, n), **over).reshape(m, n, n)
    J = ((Fp - Fm) / (2 * eps)).transpose(0, 2, 1)
    bad = ~np.isfinite(J)
    if bad.any():
        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            F0 = real.F(Q, **over)
            fwd = ((Fp - F0[:, None, :]) / eps).transpose(0, 2, 1)
            bwd = ((F0[:, None, :] - Fm) / eps).transpose(0, 2, 1)
        J = np.where(bad & np.isfinite(fwd), fwd, J)
        J = np.where(~np.isfinite(J) & np.isfinite(bwd), bwd, J)
    return J


# ------------------------------------------------------------------------------------------------ spectra
def jacobian(real: Realization, q: np.ndarray, over: Optional[Dict] = None, eps: float = 1e-6) -> np.ndarray:
    """Central differences. Where one side of the stencil leaves the set on which the drift is finite (a state on the
    edge of an orthant under a non-integer power, x = 0 under sqrt(x)) and the drift at q and on the other side is
    finite, that entry is the one-sided difference from the side inside the set; where every value is finite the
    arithmetic is unchanged. At an edge where the one-sided derivative is infinite (sqrt(x) at 0) the difference is a
    large finite rate of order eps^(-1/2)."""
    over, q = over or {}, np.asarray(q, float)
    n = q.size
    E = eps * np.eye(n)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        Fp = real.F(q[None] + E, **over)
        Fm = real.F(q[None] - E, **over)
    J = ((Fp - Fm) / (2 * eps)).T
    bad = ~np.isfinite(J)
    if bad.any():
        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            F0 = real.F(q[None], **over)[0]
            fwd = ((Fp - F0[None]) / eps).T
            bwd = ((F0[None] - Fm) / eps).T
        J = np.where(bad & np.isfinite(fwd), fwd, J)
        J = np.where(~np.isfinite(J) & np.isfinite(bwd), bwd, J)
    return J


def kappa_spectrum(real: Realization, q: np.ndarray, over: Optional[Dict] = None) -> np.ndarray:
    """kappa = -eigenvalues of the Jacobian, sorted by real part (complex for non-normal, rotating flows)."""
    k = -np.linalg.eigvals(jacobian(real, q, over))
    return k[np.argsort(k.real)]


def n_unstable(kappa: np.ndarray, scale: Optional[float] = None) -> int:
    scale = scale if scale is not None else max(1.0, float(np.abs(kappa.real).max()))
    return int(np.sum(kappa.real < -1e-7 * scale))


def classify(kappa: np.ndarray, tol: float = 1e-3) -> Dict[str, object]:
    re = kappa.real
    return {"unstable": int(np.sum(re < -tol)), "flat": int(np.sum(np.abs(re) <= tol)),
            "stable": int(np.sum(re > tol)), "kappa_min": float(re.min()), "kappa_max": float(re.max()),
            "rotating": bool(np.any(np.abs(kappa.imag) > tol))}


# ------------------------------------------------------------------------------------------------ states
def stored_states(real: Realization, rng, n_starts: int = 48, over: Optional[Dict] = None):
    """Relax from random starts. Returns states, basin counts, spectra, the number of starts that did not
    settle, and whether the states form a continuum (flat directions beyond the declared symmetries).

    A stable zero near the end of an implicit run that has not settled at it is a state as well (correction 35): its
    own rates make it one, while the run, which turns about it faster than it decays, stays among the starts that did
    not settle. Its basin count is the number of such runs that end near it."""
    q0 = real.carrier.sample(rng, n_starts)
    zeros: List = []
    q, speed = relax(real, q0, over, zeros=zeros)
    conv = speed < 1e-4
    for i in np.flatnonzero(conv):
        conv[i] = _settled(real, q[i], over)
    states: List[np.ndarray] = []
    counts: List[int] = []
    tol = _merge_tol(real)
    for qi in list(q[conv]) + [z for i, z in zeros if not conv[i]]:
        ci = _canon(real, qi)
        hit = next((k for k, s in enumerate(states) if real.carrier.distance(ci, _canon(real, s)) < tol), None)
        if hit is None:
            states.append(qi)
            counts.append(1)
        else:
            counts[hit] += 1
    states, counts = _one_state_per_zero(real, states, counts, over)
    order = np.argsort(counts, kind="stable")[::-1]
    states = [states[i] for i in order]
    counts = [counts[i] for i in order]
    spectra = [kappa_spectrum(real, s, over) for s in states]
    keep = [i for i, k in enumerate(spectra) if classify(k)["unstable"] == 0]
    states, counts, spectra = [states[i] for i in keep], [counts[i] for i in keep], [spectra[i] for i in keep]
    n_null = 0 if real.null_modes is None else len(np.atleast_2d(real.null_modes))
    flat = all(classify(k)["flat"] > n_null for k in spectra)
    # many flat states: starts on a circle of states also merge, so that fewer states than half of the starts remain
    continuum = bool(len(states) >= max(6, int(conv.sum()) // 3) and flat)
    if not continuum and flat and len(states) >= 2:
        continuum = all(_on_a_set_of_equilibria(real, s, over) for s in states)
    return states, counts, spectra, int(np.sum(~conv)), continuum


def _settled(real: Realization, q: np.ndarray, over: Optional[Dict]) -> bool:
    """Whether a start whose speed fell below the tolerance has settled near a zero of the drift. A small speed also
    marks slow motion far from any zero: along the slow manifold of a relaxation oscillation (x' = mu (y - x^3/3 + x),
    y' = -x / mu, whose speed there is |x| / mu) or a slow drift that no zero ends. Such starts were read as stored
    states, and their noisy hold time ran at the step of the fast rate. Where one least-squares step of Newton, J d = F,
    is shorter than the merge tolerance and accounts for at least half the drift, the start is near a zero and nothing
    changes. Otherwise Newton runs from the start: it has settled if Newton converges to a zero within the merge
    tolerance, or to a zero without a growing direction, which the start approaches slowly. A start that Newton takes to
    no zero, or to a zero with a growing direction far from it, has not settled.

    Near a zero the same question is asked (correction 26). The slow branch of a relaxation oscillation can pass within
    the merge tolerance of its unstable steady state, as the left branch of the Oregonator passes 0.04 from it. A start
    there relaxes without a growing direction, at the rates of the attracting branch, while the zero it is near grows:
    it is passing the zero, not settled at it, and it was read as a stable state with a flat direction. A start near a
    zero with a growing direction whose own rates have none has not settled. A start whose own rates show the growth, as
    one that stopped near a saddle, is settled as before; it is not kept as a state, since its rates grow."""
    with np.errstate(all="ignore"):
        f = real.F(q[None], **(over or {}))[0]
        J = jacobian(real, q, over)
    tol = _merge_tol(real)
    if np.all(np.isfinite(f)) and np.all(np.isfinite(J)):
        d = np.linalg.lstsq(J, f, rcond=None)[0]
        near = np.linalg.norm(d) < tol and np.linalg.norm(J @ d - f) <= 0.5 * np.linalg.norm(f)
    else:
        return True  # nothing to tell from: as before
    with np.errstate(all="ignore"):
        z, ok = newton(real, q, over)
    if not ok:
        return near
    if not near and real.carrier.distance(_canon(real, z), _canon(real, q)) >= tol:
        return classify(kappa_spectrum(real, z, over))["unstable"] == 0
    if classify(-np.linalg.eigvals(J))["unstable"] > 0:
        return True
    return classify(kappa_spectrum(real, z, over))["unstable"] == 0


def _one_state_per_zero(real: Realization, states: List[np.ndarray], counts: List[int], over: Optional[Dict]):
    """Relaxed states that approach one zero of the drift are one state. Where the restoring force vanishes at linear
    order (a degenerate state, as x' = -x^3), the approach is a power of time, and the relaxation stops where the speed
    falls below its tolerance, at distances from the zero that differ with the start and from one side to the other.
    The starts then cluster into two or more states with no barrier between them. Newton from each relaxed state finds
    the zero it approaches. States whose zeros coincide, at a zero with no growing direction, are counted as one state,
    placed at that zero, with the starts of all of them. A state whose zero is its own, as every state of a landscape
    whose states are not degenerate, is kept as it was relaxed."""
    if len(states) < 2:
        return states, counts
    tol = _merge_tol(real)
    zeros = []
    for s in states:
        z, ok = newton(real, s, over)
        zeros.append(z if ok else None)
    groups: List[List[int]] = []
    for i, z in enumerate(zeros):
        g = None if z is None else next(
            (g for g in groups if zeros[g[0]] is not None
             and real.carrier.distance(_canon(real, z), _canon(real, zeros[g[0]])) < tol), None)
        if g is None:
            groups.append([i])
        else:
            g.append(i)
    out_states: List[np.ndarray] = []
    out_counts: List[int] = []
    for g in groups:
        if len(g) > 1 and classify(kappa_spectrum(real, zeros[g[0]], over))["unstable"] == 0:
            out_states.append(zeros[g[0]])
            out_counts.append(int(sum(counts[i] for i in g)))
        else:
            out_states += [states[i] for i in g]
            out_counts += [counts[i] for i in g]
    return out_states, out_counts


def _on_a_set_of_equilibria(real: Realization, s: np.ndarray, over: Optional[Dict]) -> bool:
    """Whether a state lies on a connected set of equilibria (a line or a surface of them, as the levels of a conserved
    total or a variable that does not move). On such a set the relaxed starts merge into fewer clusters than the count
    of the continuum requires, and were read as separate states with no barrier between them. The zero that the state
    approaches (Newton) is moved by half the merge tolerance along one of its flat directions beyond the declared
    symmetries, and corrected within the plane normal to that direction: a set of equilibria crosses the plane, so the
    correction ends on a zero at that distance; an isolated zero, also a degenerate one (x' = -x^3), does not."""
    over = over or {}
    z, ok = newton(real, s, over)
    if not ok:
        return False
    ev, vec = np.linalg.eig(jacobian(real, z, over))
    basis = None
    if real.null_modes is not None:
        basis = np.linalg.qr(np.atleast_2d(np.asarray(real.null_modes, float)).T)[0]
    delta = 0.5 * _merge_tol(real)
    for k in np.flatnonzero(np.abs(ev.real) <= 1e-3):
        for v in (vec[:, k].real, vec[:, k].imag):
            if basis is not None:
                v = v - basis @ (basis.T @ v)
            if np.linalg.norm(v) < 1e-6:
                continue
            v = v / np.linalg.norm(v)
            for sign in (1.0, -1.0):
                p = real.carrier.wrap(z + sign * delta * v)
                q, ok = _zero_in_plane(real, p, v, over)
                if ok and real.carrier.distance(q, z) >= 0.5 * delta and real.carrier.distance(q, p) <= 0.5 * delta:
                    return True
    return False


def _zero_in_plane(real: Realization, p: np.ndarray, v: np.ndarray, over: Dict, iters: int = 80, tol: float = 1e-11):
    """Newton for a zero of the drift on the plane through p normal to the unit vector v (minimal-norm steps of the
    drift's equations together with the condition that the step has no component along v)."""
    q = np.array(p, float)
    cap = 0.5 * base_length(real)
    for _ in range(iters):
        f = real.F(q, **over)
        if not np.all(np.isfinite(f)):
            return q, False  # outside the set where the drift is finite (see newton)
        if np.linalg.norm(f) < tol:
            return q, True
        A = np.vstack([jacobian(real, q, over), v[None]])
        if not np.all(np.isfinite(A)):
            return q, False
        step = np.linalg.lstsq(A, np.concatenate([f, [0.0]]), rcond=None)[0]
        n = np.linalg.norm(step)
        if n > cap:
            step *= cap / n
        q = real.carrier.wrap(q - step)
    return q, bool(np.linalg.norm(real.F(q, **over)) < 1e-8)


NO_STATE = "no stable state: the trajectories do not settle"
LIMIT_CYCLE = "no stable state: a limit cycle, whose phase diffuses (Law 2 along the cycle)"


def loss_law(n_states: int, spectra: Sequence[np.ndarray], n_unconverged: int, n_starts: int,
             continuum: bool = False, n_null: int = 0) -> str:
    if n_states == 0 and n_unconverged > 0:
        # a start that does not settle oscillates, escapes or blows up in finite time; that it moves on a limit cycle
        # is read only where the search for a cycle finds one (discovery._cycle_isolation, correction 29)
        return NO_STATE
    if continuum:
        return "a continuum of neutral states: diffusion along the flat directions (Law 2)"
    flat = max((classify(k)["flat"] - n_null for k in spectra), default=0)
    if n_states >= 2:
        s = "activation between stored states (Law 3)"
        return s + (f", with {flat} flat direction(s) diffusing (Law 2)" if flat > 0 else "")
    if flat > 0:
        return f"diffusion along {flat} flat direction(s) (Law 2)"
    return "relaxation to a single state (Law 1)" + (
        "; part of the trajectories do not settle (oscillation)" if n_unconverged > n_starts // 4 else "")


def newton(real: Realization, q: np.ndarray, over: Optional[Dict] = None, iters: int = 80, tol: float = 1e-11,
           descent: bool = False):
    """Damped Newton for a zero of the drift. The search ends as not converged at an iterate where the drift or its
    central-difference Jacobian is not finite: a start or a step outside the set where the drift is defined (a
    non-integer power or the square root of a negative value, a pole, a control value that divides the drift). Least
    squares on such a matrix raised "SVD did not converge", and every caller takes the eigenvalues of the Jacobian at a
    converged zero. Since correction 23 a stencil with one side outside that set takes the one-sided difference, so a
    zero on the edge of the set converges, with the rate on its inner side.

    Each step is capped at half the base length of the carrier, so 80 steps reach no zero more than 40 base lengths away.
    With descent the step is not capped: it is halved until the drift decreases, at most 40 times, and the search ends as
    not converged where no halving decreases it (correction 35). A pendulum that turns at 2e4 on a carrier of scale 1
    ends its relaxation runs at velocities of order 1e4, from which the capped search reached no zero."""
    over = over or {}
    q = real.carrier.wrap(np.array(q, float))
    cap = 0.5 * base_length(real)
    for _ in range(iters):
        f = real.F(q, **over)
        if not np.all(np.isfinite(f)):
            return q, False
        J = jacobian(real, q, over)
        if not np.all(np.isfinite(J)):
            return q, False
        if np.linalg.norm(f) < tol:
            return q, True
        step = np.linalg.lstsq(J, f, rcond=None)[0]
        if descent:
            size, t = float(np.linalg.norm(f)), 1.0
            for _ in range(40):
                trial = real.carrier.wrap(q - t * step)
                with np.errstate(all="ignore"):
                    ft = real.F(trial, **over)
                if np.all(np.isfinite(ft)) and float(np.linalg.norm(ft)) < size:
                    break
                t *= 0.5
            else:
                return q, False
            q = trial
            continue
        n = np.linalg.norm(step)
        if n > cap:
            step *= cap / n
        q = real.carrier.wrap(q - step)
    f = real.F(q, **over)
    return q, bool(np.linalg.norm(f) < 1e-8 and np.all(np.isfinite(jacobian(real, q, over))))


def fixed_points(real: Realization, rng, n_starts: int = 40, over: Optional[Dict] = None, seeds=None):
    """Zeros of the drift by damped Newton from random starts (and from given seeds, such as stored states),
    stable and unstable alike."""
    starts = real.carrier.sample(rng, n_starts)
    if seeds is not None and len(seeds):
        starts = np.concatenate([np.atleast_2d(np.asarray(seeds, float)), starts])
    return zeros_from(real, starts, over)


def zeros_from(real: Realization, starts: np.ndarray, over: Optional[Dict] = None):
    """Zeros of the drift by damped Newton from the given starts; zeros closer than the merge tolerance of stored
    states are one."""
    over = over or {}
    found: List[np.ndarray] = []
    tol = _merge_tol(real)
    for q in starts:
        q, ok = newton(real, q, over)
        if ok and not any(real.carrier.distance(_canon(real, q), _canon(real, s)) < tol for s in found):
            found.append(q)
    return [(q, kappa_spectrum(real, q, over)) for q in found]


# ------------------------------------------------------------------------------------------------ barrier
NEB_TOL = 1e-5     # a band is relaxed when the largest force on an image is below this (units of the drift)
MAX_SENSES = 12    # with more coordinates at half a period, not every combination of senses is examined


def _band_force(real: Realization, path: np.ndarray, over: Dict, k_spring: float, climbing: bool):
    """Energies of the images and forces on the interior ones: the drift perpendicular to the band and a spring
    force along it, with the tangent taken toward the neighbour of higher energy (Henkelman & Jonsson 2000). If
    climbing, the highest interior image instead feels the drift with its component along the band reversed, which
    moves it to the saddle (Henkelman, Uberuaga & Jonsson 2000). Also returns the largest drift on an image."""
    c = real.carrier
    E, F = real.V(path, **over), real.F(path, **over)[1:-1]
    fwd, back = c.diff(path[1:-1], path[2:]), c.diff(path[:-2], path[1:-1])
    up, down = E[2:] - E[1:-1], E[:-2] - E[1:-1]
    big = np.maximum(np.abs(up), np.abs(down))[:, None]
    small = np.minimum(np.abs(up), np.abs(down))[:, None]
    tau = np.where(((up > 0) & (down < 0))[:, None], fwd, np.where(((up < 0) & (down > 0))[:, None], back, np.where(
        (up > down)[:, None], fwd * big + back * small, fwd * small + back * big)))
    tau = tau / np.maximum(np.linalg.norm(tau, axis=1, keepdims=True), 1e-300)
    along = np.sum(F * tau, axis=1, keepdims=True) * tau
    force = F - along + k_spring * (np.linalg.norm(fwd, axis=1) - np.linalg.norm(back, axis=1))[:, None] * tau
    if climbing:
        top = int(np.argmax(E[1:-1]))
        force[top] = F[top] - 2.0 * along[top]
    return E, force, float(np.max(np.linalg.norm(F, axis=1)))


def _relax_band(real: Realization, path: np.ndarray, over: Dict, k_spring: float, tol: float, iters: int):
    """Relax a band in two stages. The images first follow the force on the band with the time step of the
    dynamics (safe_dt), which decides the path that the band approaches. Once the largest force on an image has
    fallen to 2% of the largest drift on the band (or after half of the iterations), the highest image climbs and
    the relaxation is completed by damped dynamics with an adaptive step (FIRE; Bitzek et al. 2006), until that
    force is below tol. Damped dynamics applied from the straight path can carry a band over a ridge to another
    path: between the two most populated states of the dipole patch it gave the barrier 6.62, and the descent
    4.19. Returns the band, its energies, whether the tolerance was reached, the iterations used and the largest
    remaining force."""
    c = real.carrier
    step = safe_dt(real)
    dt_max = 1.0 / np.sqrt(0.3 / step + 4.0 * k_spring)  # from the largest rates of the drift and of the springs
    max_move = 0.5 * float(np.mean(np.linalg.norm(c.diff(path[:-1], path[1:]), axis=1)))
    v = np.zeros_like(path[1:-1])
    dt, mix, n_down, damped = 0.1 * dt_max, 0.1, 0, False
    for it in range(iters):
        E, force, drift = _band_force(real, path, over, k_spring, damped)
        res = float(np.max(np.linalg.norm(force, axis=1)))
        if damped and res < tol:
            return path, E, True, it, res
        if not damped and (res < max(tol, 2e-2 * drift) or it >= iters // 2):
            damped = True
            continue
        if not damped:
            move = step * force
        else:
            if np.sum(force * v) >= 0.0:
                v = (1.0 - mix) * v + mix * np.linalg.norm(v) * force / np.linalg.norm(force)
                n_down += 1
                if n_down > 5:
                    dt, mix = min(1.1 * dt, dt_max), 0.99 * mix
            else:
                v, dt, mix, n_down = np.zeros_like(v), 0.5 * dt, 0.1, 0
            v = v + dt * force
            move = dt * v
            move *= min(1.0, max_move / max(float(np.max(np.linalg.norm(move, axis=1))), 1e-300))
        path = np.concatenate([path[:1], c.wrap(path[1:-1] + move), path[-1:]])
    E, force, _ = _band_force(real, path, over, k_spring, True)
    return path, E, False, iters, float(np.max(np.linalg.norm(force, axis=1)))


def _half_period(real: Realization, a: np.ndarray, b: np.ndarray, over: Dict) -> np.ndarray:
    """Coordinates of a torus in which two states differ by half a period, so that the move from one state to the
    other is equally short in both senses. The comparison allows for the accuracy to which the states are known:
    ten times the Newton correction of each state, and not less than 1e-9 of the period."""
    c = real.carrier
    if c.kind != "torus":
        return np.zeros(0, int)
    slack = 1e-9 * c.period + 10.0 * sum(
        np.abs(np.linalg.lstsq(jacobian(real, q, over), real.F(q, **over), rcond=None)[0]) for q in (a, b))
    return np.flatnonzero(0.5 * c.period - np.abs(c.diff(a, b)) <= slack)


def _straight_paths(real: Realization, a: np.ndarray, b: np.ndarray, over: Dict, images: int, limit: int):
    """Displacements of the shortest straight paths from a to b. With k coordinates at half a period (_half_period)
    there are 2^k of them, one for every combination of senses. They are ranked by the highest energy met on the
    path, and the `limit` lowest are returned, lowest first, each with its senses and that energy above a. For
    k <= MAX_SENSES every combination is examined; beyond, those met when one sense at a time is reversed as long
    as the highest energy decreases, starting from all senses positive and from all negative. Returns the
    coordinates at half a period, the number of combinations examined and the selected paths."""
    c = real.carrier
    d = c.diff(a, b)
    half = _half_period(real, a, b, over)
    k = len(half)
    plus = np.where(d[half] > 0, d[half], d[half] + c.period)
    minus, t = plus - c.period, np.linspace(0, 1, images)[:, None]
    e0 = float(real.V(a, **over))

    def peaks(bits):  # bits: one row per path, 1 where the coordinate moves in the negative sense
        out = np.empty(len(bits))
        for i in range(0, len(bits), 256):
            disp = np.repeat(d[None], len(bits[i:i + 256]), axis=0)
            disp[:, half] = np.where(bits[i:i + 256] > 0, minus, plus)
            pts = c.wrap(a[None, None] + t[None] * disp[:, None])
            out[i:i + 256] = real.V(pts.reshape(-1, a.size), **over).reshape(len(disp), images).max(axis=1) - e0
        return out

    if k <= MAX_SENSES:
        bits = (np.arange(2 ** k)[:, None] >> np.arange(k)[None]) & 1
    else:
        seen = set()
        for cur in (np.zeros(k, int), np.ones(k, int)):
            best = float(peaks(cur[None])[0])
            while True:
                near = (cur[None] + np.eye(k, dtype=int)) % 2
                values = peaks(near)
                seen.update(map(tuple, near.tolist()), [tuple(cur.tolist())])
                if not values.min() < best:
                    break
                cur, best = near[int(np.argmin(values))], float(values.min())
        bits = np.array(sorted(seen))
    top = peaks(bits)
    # paths of equal highest energy (images of each other under a symmetry) are ordered by their senses
    order = np.argsort(top, kind="stable")
    tier = np.concatenate([[0], np.cumsum(np.diff(top[order]) > 1e-9 * max(1.0, float(np.ptp(top))))])
    chosen = []
    for i in sorted(range(len(order)), key=lambda j: (tier[j], bits[order[j]].tolist()))[:limit]:
        disp = d.copy()
        disp[half] = np.where(bits[order[i]] > 0, minus, plus)
        chosen.append((disp, "".join("-" if s else "+" for s in bits[order[i]]), float(top[order[i]])))
    return half, len(bits), chosen


def neb_barrier(real: Realization, a: np.ndarray, b: np.ndarray, over: Optional[Dict] = None,
                images: int = 24, iters: int = 20000, k_spring: float = 1.0, tol: float = NEB_TOL,
                n_paths: int = 8) -> Dict[str, object]:
    """Barrier between two stable states of a gradient system along a minimum-energy path.

    A band of images is relaxed from the straight path between the states (nudged elastic band with a climbing
    image, _relax_band) until the largest force on an image is below tol, for at most iters iterations. The record
    states whether the tolerance was reached, the iterations and the remaining force, and the number of unstable
    directions at the highest image (one at the saddle of a minimum-energy path).

    On a torus the shortest straight path is not unique if the states differ by half a period in k coordinates:
    each of them reaches its value in the other state in either sense, and the sense read from the rounded
    difference is arbitrary. The n_paths straight paths of lowest maximum energy are then relaxed
    (_straight_paths), and the barrier is the lowest among the relaxed bands; `senses` gives the senses of these
    coordinates in the reported band, `paths` lists every band, `distinct_barriers` the different barriers found,
    and `method` states in words what was done. The search does not establish the lowest barrier over all paths
    between the two states."""
    over = over or {}
    a, b = np.asarray(a, float), np.asarray(b, float)
    half, n_examined, starts = _straight_paths(real, a, b, over, images, n_paths)
    t = np.linspace(0, 1, images)[:, None]
    rows = []
    for disp, senses, peak in starts:
        path, E, ok, n_it, res = _relax_band(real, real.carrier.wrap(a[None] + t * disp[None]), over, k_spring,
                                             tol, iters)
        top = int(np.argmax(E))
        rows.append({"senses": senses, "straight_barrier": peak, "barrier": float(E[top] - E[0]),
                     "barrier_back": float(E[top] - E[-1]), "saddle_index": top, "converged": ok,
                     "iterations": n_it, "residual": res, "energies": E.tolist(),
                     "unstable_directions": n_unstable(kappa_spectrum(real, path[top], over))})
    pool = [r for r in rows if r["converged"]] or rows
    low = min(r["barrier"] for r in pool)
    same = 1e-6 * max(1.0, abs(low))
    best = next(r for r in pool if r["barrier"] <= low + same)
    distinct: List[float] = []
    for value in sorted(r["barrier"] for r in pool):
        if not distinct or value - distinct[-1] > same:
            distinct.append(value)
    out = {key: best[key] for key in ("barrier", "barrier_back", "energies", "saddle_index", "converged",
                                      "iterations", "residual")}
    out.update(tolerance=tol, unstable_directions=best["unstable_directions"], images=images,
               half_period_coordinates=half.tolist(), senses=best["senses"], straight_paths=2 ** len(half),
               straight_paths_examined=n_examined,
               paths=[{key: value for key, value in r.items() if key != "energies"} for r in rows],
               distinct_barriers=distinct)
    out["method"] = _barrier_method(out)
    return out


def _barrier_method(r: Dict[str, object]) -> str:
    """The calculation of a barrier in words, for the card."""
    n, k, m = len(r["paths"]), len(r["half_period_coordinates"]), len(r["distinct_barriers"])
    band = (f"a nudged elastic band of {r['images']} images (descent along the force on the band, then damped "
            f"dynamics with a climbing image, until the largest force on an image is below {r['tolerance']:g})")
    which = "the band" if n == 1 else "the reported band"
    if k == 0:
        text = f"The straight path between the two states was relaxed as {band}."
    else:
        text = (f"The two states differ by half a period in {k} coordinate{'s' if k > 1 else ''}, each of which "
                f"reaches its value in the other state in either sense: {r['straight_paths']} straight paths are "
                f"equally short"
                + ("" if r["straight_paths_examined"] == r["straight_paths"] else
                   f", of which {r['straight_paths_examined']} were examined") + ". ")
        if n == 1:
            text += (f"The one of lowest maximum energy was relaxed as {band}. The minimum-energy path is not "
                     f"unique, and bands in other senses were not relaxed.")
        else:
            done = sum(p["converged"] for p in r["paths"])
            text += (("Each was" if n == r["straight_paths"] else f"The {n} of lowest maximum energy were each")
                     + f" relaxed as {band}; {done} of {n} reached the tolerance. The barrier is the lowest among "
                     + ("these bands" if done in (0, n) else "the bands that reached it")
                     + f", which have {m} distinct barrier{'s' if m > 1 else ''}. The minimum-energy path is not "
                       f"unique: the bands differ in the senses of these coordinates.")
        text += " A lower barrier on a path that was not relaxed is not excluded."
    if not r["converged"]:
        text += (f" The tolerance was not reached by {which} in {r['iterations']} iterations (remaining force "
                 f"{r['residual']:.2g}), so the barrier is not converged.")
    elif r["unstable_directions"] == 1:
        text += f" The highest image of {which} is a saddle with one unstable direction."
    elif r["unstable_directions"] > 1:
        text += (f" The highest image of {which} has {r['unstable_directions']} unstable directions: the band lies "
                 f"on a ridge, and the barrier is an upper bound for the paths beside it.")
    else:
        text += f" The energy along {which} is highest at one of the two states."
    return text


def restoring_threshold(real: Realization, a: np.ndarray, b: np.ndarray, n: int = 64) -> float:
    """Largest force opposing a straight move from a to b: a bounded write must exceed it to remove the
    barrier on that path; a weaker write needs thermal activation."""
    t = np.linspace(0, 1, n)[:, None]
    d = real.carrier.diff(a, b)
    path = real.carrier.wrap(a[None] + t * d[None])
    u = d / max(np.linalg.norm(d), 1e-12)
    return float(max(1e-6, np.max(-(real.F(path) @ u))))


# ------------------------------------------------------------------------------------------------ hold / write
def basin_balls(real: Realization, states: Sequence[np.ndarray], over: Optional[Dict] = None, n_random: int = 16,
                halvings: int = 8) -> np.ndarray:
    """Radius of a ball around each stored state that lies in the basin of the state. The radius starts at 0.3 of the
    distance to the nearest other state, so that the balls are disjoint, and is halved until every test point on the
    boundary (along the coordinate axes, toward the other states and along random directions) relaxes to the state
    without noise. A state for which no radius passes gets 0, and its configurations are identified when they settle."""
    over = over or {}
    rng = np.random.default_rng(12345)  # fixed test directions; the generator of the simulation is not advanced
    canon = [_canon(real, np.asarray(s, float)) for s in states]
    radii = np.zeros(len(states))
    for i, s in enumerate(states):
        s = np.asarray(s, float)
        others = [j for j in range(len(states)) if j != i]
        d = min((float(real.carrier.distance(canon[i], canon[j])) for j in others), default=base_length(real))
        u = np.concatenate([np.eye(s.size), -np.eye(s.size), rng.standard_normal((n_random, s.size))]
                           + [real.carrier.diff(s, states[j])[None] for j in others])
        u /= np.maximum(np.linalg.norm(u, axis=1, keepdims=True), 1e-300)
        r = 0.3 * d
        for _ in range(halvings + 1):
            end, _ = relax(real, s[None] + r * u, over)
            if np.all(real.carrier.distance(_canon(real, end), canon[i][None]) < 0.5 * r):
                radii[i] = r
                break
            r /= 2
    return radii


def basin_of(real: Realization, q: np.ndarray, states: Sequence[np.ndarray], radii: np.ndarray,
             over: Optional[Dict] = None, t_max: float = 300.0, tol: float = 1e-5) -> np.ndarray:
    """Index of the stored state whose basin contains each configuration. A configuration inside a ball of basin_balls
    belongs to its state; any other is carried by the noise-free flow until it enters a ball, or until it settles at a
    stable state (then the stored state within the merging tolerance of stored_states). -1 for a configuration that
    settles at a stable state that is not stored, or that has done neither within t_max. The nearest stored state is
    the same only where the separatrix is the bisector of the states, as for two states related by a reflection."""
    over, dt = over or {}, safe_dt(real)
    canon = [_canon(real, np.asarray(s, float)) for s in states]
    radii = np.asarray(radii, float)
    tol_state = _merge_tol(real)

    def ball(x):
        inside = np.stack([real.carrier.distance(_canon(real, x), c[None]) for c in canon], axis=-1) < radii
        return np.where(inside.any(axis=-1), np.argmax(inside, axis=-1), -1)

    x = real.carrier.wrap(np.atleast_2d(np.asarray(q, float)))
    label = ball(x)
    todo = np.flatnonzero(label < 0)
    x = x[todo]
    # beyond the explicit step budget the flow takes the drift-implicit step at the declared step, as the trajectories
    # of integrate do (correction 25 of the discovery engine): at a slaved rate of 1e7 the explicit step is 3e-8, and
    # t_max is 1e10 such steps; a rate far above 1 / dt relaxes within one implicit step, and the slow flow is as before
    stiff = t_max / dt > EXPLICIT_BUDGET and real.dt > dt
    if stiff:
        dt = float(real.dt)
    for k in range(1, int(t_max / dt) + 1):
        if not todo.size:
            break
        f = real.F(x, **over)
        x = real.carrier.wrap(x + (_drift_implicit(real, x, over, dt, dt * f) if stiff else dt * f))
        if k % 10:
            continue
        hit = ball(x)
        done = (hit >= 0) | ~np.all(np.isfinite(x), axis=-1)
        # a configuration that has stopped outside every ball has settled at a stable state, unless it stopped near a
        # saddle, from which the flow carries it on
        for j in np.flatnonzero(~done & (np.linalg.norm(f, axis=-1) < tol)):
            if n_unstable(kappa_spectrum(real, x[j], over)) == 0:
                d = [float(real.carrier.distance(_canon(real, x[j]), c)) for c in canon]
                done[j], hit[j] = True, int(np.argmin(d)) if min(d) < tol_state else -1
        label[todo[done]] = hit[done]
        todo, x = todo[~done], x[~done]
    return label


def _fraction_in_basin(real: Realization, snaps: np.ndarray, states, radii, over: Optional[Dict], idx: int):
    """Fraction of the recorded configurations (records x trajectories x dim) in the basin of states[idx] at each
    record. The records are identified together, after the simulation, in one integration of the noise-free flow."""
    snaps = np.asarray(snaps, float)
    label = basin_of(real, snaps.reshape(-1, snaps.shape[-1]), states, radii, over)
    return np.mean(label.reshape(snaps.shape[:-1]) == idx, axis=-1)


def hold_time(real: Realization, idx: int, states, D: float, rng, n_traj: int = 200, t_max: float = 400.0,
              over: Optional[Dict] = None, radii: Optional[np.ndarray] = None) -> Dict[str, object]:
    """Fraction of the trajectories started in states[idx] that remain in its basin (basin_of), and the time at which
    it falls to P_STAR."""
    radii = basin_balls(real, states, over) if radii is None else radii
    q = np.repeat(states[idx][None], n_traj, axis=0)
    _, t, snaps = integrate(real, q, t_max, D, rng, over, record=np.copy, n_record=200)
    frac = _fraction_in_basin(real, snaps, states, radii, over, idx)
    below = np.flatnonzero(frac <= P_STAR)
    return {"time": float(t[below[0]]) if below.size else float(t_max), "censored": bool(below.size == 0),
            "t": t.tolist(), "fraction": frac.tolist()}


def write_test(real: Realization, i_old: int, i_new: int, states, h: float, T_w: float, T_rel: float, D: float,
               rng, n_traj: int = 200, mode: str = "fixed", sweep_from: Optional[float] = None,
               over: Optional[Dict] = None, radii: Optional[np.ndarray] = None) -> Dict[str, object]:
    """Fraction of the trajectories in the basin of states[i_new] (basin_of, without the field and at the parameters
    of the stored states) during a bounded write from states[i_old] and after its release."""
    over = dict(over or {})
    radii = basin_balls(real, states, over) if radii is None else radii
    q = np.repeat(states[i_old][None], n_traj, axis=0)
    target = states[i_new]
    path = None
    if mode == "sweep":
        v1 = over.get(real.control, real.params[real.control])
        path = lambda t: sweep_from + (v1 - sweep_from) * min(1.0, t / T_w)
    q, t1, s1 = integrate(real, q, T_w, D, rng, over, write=(target, h), control_path=path, record=np.copy)
    q, t2, s2 = integrate(real, q, T_rel, D, rng, over, record=np.copy)
    frac = _fraction_in_basin(real, np.concatenate([s1, s2]), states, radii, over, i_new)
    f1, f2 = frac[:len(t1)], frac[len(t1):]
    reach = np.flatnonzero(f1 >= P_STAR)
    return {"mode": mode, "h": h, "T_w": T_w, "accuracy": float(f2[-1]) if f2.size else float("nan"),
            "rewrite_time": float(t1[reach[0]]) if reach.size else float("nan"),
            "t": np.concatenate([t1, T_w + t2]).tolist(), "fraction": np.concatenate([f1, f2]).tolist(),
            "mean_overlap": float(np.mean(real.overlap(q, target)))}


# ------------------------------------------------------------------------------------------------ write points
def _dq_dv(real: Realization, q: np.ndarray, param: str, v: float) -> np.ndarray:
    h = 1e-6 * max(1.0, abs(v))
    J = jacobian(real, q, {param: v})
    dF = (real.F(q, **{param: v + h}) - real.F(q, **{param: v - h})) / (2 * h)
    if not (np.all(np.isfinite(J)) and np.all(np.isfinite(dF))):
        return np.full(q.size, np.nan)  # no predictor: track halves its step, as for a corrector that fails
    return -np.linalg.lstsq(J, dF, rcond=None)[0]


def track(real: Realization, param: str, q0: np.ndarray, v0: float, v1: float, n_init: int = 16,
          max_halvings: int = 18):
    """Follow a state from v0 toward v1 by predictor-corrector continuation. Stops where it loses
    stability (an eigenvalue crosses zero) or vanishes (no nearby solution beyond this value)."""
    q = np.asarray(q0, float)
    v = float(v0)
    k = kappa_spectrum(real, q, {param: v})
    scale = max(1.0, float(np.abs(k.real).max()))
    nu = n_unstable(k, scale)
    branch = [{"v": v, "q": q.tolist(), "kappa_min": float(k.real.min())}]
    full = float(v1 - v0)
    step = full / n_init
    min_step = abs(full) / 2 ** max_halvings
    max_move = 0.5 * _merge_tol(real)
    while abs(v1 - v) > 1e-12 * max(1.0, abs(v1)):
        st = step if abs(step) <= abs(v1 - v) else v1 - v
        dq = _dq_dv(real, q, param, v) * st
        accepted = False
        if np.linalg.norm(dq) <= max_move:
            qn, ok = newton(real, real.carrier.wrap(q + dq), {param: v + st})
            accepted = ok and real.carrier.distance(qn, q) <= 2.0 * np.linalg.norm(dq) + 0.05 * _merge_tol(real)
        if not accepted:
            if abs(st) / 2 < min_step:
                return branch, {"type": "vanishes", "v": v, "q": q.tolist()}
            step = st / 2
            continue
        kn = kappa_spectrum(real, qn, {param: v + st})
        if n_unstable(kn, scale) != nu:
            return branch, _bisect_loss(real, param, q, v, v + st, nu, scale)
        v, q = v + st, qn
        branch.append({"v": v, "q": q.tolist(), "kappa_min": float(kn.real.min())})
        step = np.sign(full) * min(abs(st) * 1.5, abs(full) / n_init)
    return branch, None


def _bisect_loss(real, param, q, va, vb, nu, scale, iters: int = 40):
    qa = q
    for _ in range(iters):
        vm = 0.5 * (va + vb)
        qm, ok = newton(real, qa, {param: vm})
        if ok and real.carrier.distance(qm, qa) < 0.5 * _merge_tol(real) \
                and n_unstable(kappa_spectrum(real, qm, {param: vm}), scale) == nu:
            va, qa = vm, qm
        else:
            vb = vm
    return {"type": "loses stability", "v": 0.5 * (va + vb), "q": qa.tolist()}


def locate_writes(real: Realization, rng, param: Optional[str] = None, values: Optional[Sequence[float]] = None,
                  n_starts: int = 24) -> List[Dict[str, object]]:
    """Every value of the parameter at which a stable state loses stability or appears (tracked forward
    and backward between neighbouring scan values).

    The stable states at a scan value are found by Newton from n_starts random starts, drawn once for each interval
    and direction in the order of the scan, and each is followed to the neighbouring scan value. The starts cover the
    carrier scale, and Newton from them need not reach a state beyond it: on the driven Kerr resonator of bifurcation
    readout (scale 1) the high-amplitude branch (n = 19 photons at its fold) is reached from none of them where it
    coexists with the low-amplitude branch, and above the low-amplitude fold, where it is the only state, Newton
    converges from none. The scan therefore continues every branch it holds. A state followed to a scan value without
    an event is followed on into the next interval, and a state found from the starts of one interval is also
    followed into the other interval at its scan value. At a scan value where no stable state is held after that and
    Newton converged from none of the starts, the starts are relaxed along the flow (relax), and the states they
    settle in are followed in both directions. Where Newton reaches an unstable state, as the focus inside a limit
    cycle, the starts are not relaxed: on a limit cycle the relaxation does not settle and runs to its end. Neither
    step draws random numbers, so the generator is left in the state in which the random starts leave it. The write
    points the two steps add are appended to those found from the random starts, which are kept where both give the
    same write point."""
    param = param or real.control
    if values is None:
        values = np.linspace(*real.control_range, 7)
    values = np.asarray(values, float)
    width = float(abs(values[-1] - values[0]))
    lo, hi = float(np.min(values)), float(np.max(values))
    last, tol = len(values) - 1, _merge_tol(real)
    # held[j, d]: the states at values[j] whose branch has been followed toward values[j + d]
    held: Dict[tuple, List[np.ndarray]] = {(j, d): [] for j in range(last + 1) for d in (1, -1)}

    def stable(q: np.ndarray, k: np.ndarray) -> bool:
        # a state on the boundary of an orthant (a species or a tube absent) sits where the drift
        # need not be smooth; its continuation is not attempted
        on_edge = real.carrier.kind == "orthant" and float(np.min(np.abs(q))) < 1e-6 * real.carrier.scale
        return n_unstable(k) == 0 and not on_edge

    def follow(j: int, d: int, q: np.ndarray):
        """The event on the branch of q from values[j] toward values[j + d], or None and the state reached there."""
        v0, v1 = values[j], values[j + d]
        # each interval is tracked slightly past its end, so that an event exactly at a scan value (where the
        # state found there is already marginal) is still located
        v_end = float(np.clip(v1 + np.sign(v1 - v0) * 1e-2 * abs(v1 - v0), lo, hi))
        held[j, d].append(_canon(real, q))
        branch, ev = track(real, param, q, v0, v_end)
        return ev, (None if ev is not None else np.asarray(branch[-1]["q"], float))

    def reach(j: int, d: int, end: np.ndarray, queue: list) -> None:
        """A branch followed in direction d without an event, at the scan value values[j] it has reached."""
        q, ok = newton(real, end, {param: values[j]})
        if ok:
            held[j, -d].append(_canon(real, q))
            if 0 <= j + d <= last:
                queue.append((j, d, q))

    def run(queue: list) -> None:
        while queue:
            j, d, q = queue.pop(0)
            c = _canon(real, q)
            if any(real.carrier.distance(c, s) < tol for s in held[j, d]) \
                    or not stable(q, kappa_spectrum(real, q, {param: values[j]})):
                continue
            ev, end = follow(j, d, q)
            if ev is not None:
                added.append(ev)
            else:
                reach(j + d, d, end, queue)

    events: List[Dict[str, object]] = []
    added: List[Dict[str, object]] = []
    starts_at: Dict[int, List[np.ndarray]] = {j: [] for j in range(last + 1)}
    zero_at = [False] * (last + 1)
    reached, found = [], []
    for i in range(last):
        for j, d in ((i, 1), (i + 1, -1)):
            starts_at[j].append(real.carrier.sample(rng, n_starts))
            zeros = zeros_from(real, starts_at[j][-1], {param: values[j]})
            zero_at[j] = zero_at[j] or bool(zeros)
            for q, k in zeros:
                if stable(q, k):
                    ev, end = follow(j, d, q)
                    found.append((j, d, q))
                    if ev is not None:
                        events.append(ev)
                    else:
                        reached.append((j + d, d, end))
    # the continuation of the branches held: a state reached at values[j] moving in direction d has its branch toward
    # values[j - d] followed already and goes on toward values[j + d]; a state found from the starts toward
    # values[j + d] goes toward values[j - d] too
    queue: list = []
    for j, d, end in reached:
        reach(j, d, end, queue)
    queue += [(j, -d, q) for j, d, q in found if 0 <= j - d <= last]
    run(queue)
    for j in range(last + 1):
        if held[j, 1] or held[j, -1] or zero_at[j] or not starts_at[j]:
            continue
        over = {param: values[j]}
        with np.errstate(over="ignore", invalid="ignore"):
            q, speed = relax(real, np.concatenate(starts_at[j]), over)
        for qi in q[speed < 1e-4]:
            qn, ok = newton(real, qi, over)
            if ok:
                queue += [(j, d, qn) for d in (1, -1) if 0 <= j + d <= last]
        run(queue)
    out: List[Dict[str, object]] = []
    for group in (events, added):
        for ev in sorted(group, key=lambda e: e["v"]):
            if not any(abs(ev["v"] - e["v"]) < WRITE_RESOLUTION * width
                       and real.carrier.distance(_canon(real, np.array(ev["q"])), _canon(real, np.array(e["q"])))
                       < 4 * tol for e in out):
                out.append(ev)
    return sorted(out, key=lambda e: e["v"])


def _slow_manifold_drift(real, q, vec, w, s, over, points: bool = False):
    """Drift along vec on the slow manifold: for each s, the transverse coordinates are solved so that the
    drift has no component off vec (adiabatic elimination of the relaxing directions). Returns the drift and
    the largest remaining transverse drift relative to the largest drift along the grid, and with points=True also
    the states on the manifold at the amplitudes s."""
    n = q.size
    if n == 1:
        f = np.array([float(w[0] * real.F(real.carrier.wrap(q + si * vec), **over)[0]) for si in s])
        return (f, 0.0, np.array([real.carrier.wrap(q + si * vec) for si in s])) if points else (f, 0.0)
    from scipy.linalg import null_space
    C = null_space(vec[None, :])
    f = np.empty(len(s))
    X = np.empty((len(s), n))
    z_side = {1: np.zeros(n - 1), -1: np.zeros(n - 1)}
    trans, total = 0.0, 0.0
    cap = 0.5 * float(np.max(np.abs(s)))
    for j in np.argsort(np.abs(s), kind="stable"):
        side = 1 if s[j] >= 0 else -1
        z = z_side[side].copy()
        for _ in range(60):
            x = real.carrier.wrap(q + s[j] * vec + C @ z)
            g = C.T @ real.F(x, **over)
            if np.linalg.norm(g) < 1e-13:
                break
            Jc = C.T @ jacobian(real, x, over) @ C
            if not (np.all(np.isfinite(g)) and np.all(np.isfinite(Jc))):
                break  # outside the set where the drift is finite: the window of the normal form is halved
            dz = np.linalg.lstsq(Jc, g, rcond=None)[0]
            nz = np.linalg.norm(dz)
            z = z - (dz * cap / nz if nz > cap else dz)
        z_side[side] = z
        x = real.carrier.wrap(q + s[j] * vec + C @ z)
        X[j] = x
        Fx = real.F(x, **over)
        f[j] = float(w @ Fx)
        trans = max(trans, float(np.linalg.norm(C.T @ Fx)))
        total = max(total, float(np.linalg.norm(Fx)))
    return (f, trans / max(total, 1e-300), X) if points else (f, trans / max(total, 1e-300))


def leading_order(real, q, vec, w, over, half: float, rate: float, scale: float,
                  levels: int = 8) -> Optional[Dict[str, float]]:
    """The leading nonlinear term c s^n of the drift on the slow manifold where its quadratic and cubic coefficients
    are not resolved at the write point and the series of series_reduction does not exist. The drift less its linear part, rate s, is split into its odd and even parts
    at the amplitudes half 2^-k; the larger part is resolved where it exceeds 1e-6 of the largest rate per unit s (the
    floor of the degenerate kind of normal_form), and n is the slope of its logarithm between the two smallest
    resolved amplitudes, rounded to an integer within 0.25. None where no two amplitudes resolve it (the drift
    vanishes along the critical direction, a line of states) or the slope is not near an integer."""
    amps = half * 2.0 ** -np.arange(levels)
    f, _ = _slow_manifold_drift(real, q, vec, w, np.concatenate([amps, -amps]), over)
    if not np.all(np.isfinite(f)):
        return None
    fp, fm = f[:levels], f[levels:]
    odd, even = 0.5 * (fp - fm) - rate * amps, 0.5 * (fp + fm)
    part = odd if np.max(np.abs(odd)) >= np.max(np.abs(even)) else even
    resolved = np.abs(part) > 1e-6 * scale * amps
    # the two smallest resolved amplitudes, adjacent so that their ratio is 2
    k = next((j for j in range(levels - 1, 0, -1) if resolved[j] and resolved[j - 1]), None)
    if k is None:
        return None
    slope = float(np.log2(abs(part[k - 1]) / abs(part[k])))
    n = int(round(slope))
    if n < 2 or abs(slope - n) > 0.25:
        return None
    return {"order": n, "coefficient": float(part[k] / amps[k] ** n), "slope": slope, "amplitude": float(amps[k])}


def _series_drift(real, q0, direction, over, size: int) -> np.ndarray:
    """The coefficients in s of the drift at q0 + direction(s); direction holds one row of coefficients per coordinate,
    constant term excluded. One row per coordinate of the result."""
    y = np.empty((1, q0.size), dtype=object)
    for i in range(q0.size):
        c = np.zeros(size)
        c[0] = q0[i]
        c[1:] = direction[i]
        y[0, i] = Series(c)
    values = real.drift(y, {**real.params, **over})
    return np.array([series_coefficients(v, size) for v in np.asarray(values, dtype=object).reshape(-1)])


def series_reduction(real: Realization, q: np.ndarray, vec: np.ndarray, w: np.ndarray, over: Dict,
                     order: int = SERIES_ORDER) -> Dict[str, object]:
    """The drift on the centre manifold of one critical direction, as a power series in s = w.(x - q), closed order by
    order with no window.

    With x(s) = q + s v + h(s), w.h = 0, h of second order, the invariance of the manifold reads
    (v + h'(s)) ds/dt = F(x(s)), so that ds/dt = w.F(x(s)) and, at order s^k, with a1 the linear coefficient,
        (J - k a1) h_k = [h'(s) ds/dt]_k (from h_2 .. h_(k-1)) - [F]_k (at h_k = 0), w.h_k = 0,
    where the left side is solved as a bordered system (Govaerts, Numerical Methods for Bifurcations of Dynamical
    Equilibria, 2000), which needs no rank decision for the slow directions of a stiff drift. A direction of a continuous
    symmetry (Realization.null_modes) remains in the kernel of the bordered matrix, and the least-squares solution has
    no component along it. The Jacobian comes from the same arithmetic, and vec and w are refined against it to rounding,
    so that a term that vanishes is not given a size by the finite differences of analysis.jacobian.

    Returns {"covered": True, "order": order, "coefficients": [a0 .. a_order], "rounding": per order, the size of a
    coefficient that is rounding of the terms it is a sum of}, or {"covered": False, "reason": ...} where the drift has no power series at q (abs at a kink, sqrt or
    log at zero) or is not written in operations that a Series supports. The coefficient a2 is invariant among the
    coordinates with the same tangent; a3 is invariant only where a2 vanishes, and at a fold it depends on the coordinate."""
    q, v, w = np.asarray(q, float), np.array(vec, float), np.array(w, float)
    n, size = q.size, order + 1
    try:
        J = np.zeros((n, n))
        for i in range(n):
            e = np.zeros((n, 1))
            e[i] = 1.0
            J[:, i] = _series_drift(real, q, e, over, 2)[:, 1]
        if not np.all(np.isfinite(J)):
            raise NotAnalytic("the Jacobian at the write point is not finite")
        lam = float(w @ J @ v) / float(w @ v)

        def bordered(M, a, b, rhs):
            # solves M x + a mu = rhs, b.x = 0
            B = np.zeros((n + 1, n + 1))
            B[:n, :n], B[:n, n], B[n, :n] = M, a, b
            return np.linalg.lstsq(B, np.append(rhs, 0.0), rcond=1e-13)[0][:n]

        for _ in range(3):
            v = v + bordered(J - lam * np.eye(n), v, v, -(J - lam * np.eye(n)) @ v)
            v = v / np.linalg.norm(v)
            w = w + bordered((J - lam * np.eye(n)).T, w, w, -(J - lam * np.eye(n)).T @ w)
            w = w / float(w @ v)
            lam = float(w @ J @ v)
        h = np.zeros((n, order))              # h[:, j] is the coefficient of s^(j + 1) of the tangent part and h
        h[:, 0] = v
        for k in range(2, order + 1):
            f = _series_drift(real, q, h, over, size)
            ds = w @ f
            # h'(s) ds/dt, without the tangent, which the projection removes
            dh = np.zeros(n)
            for j in range(2, k):
                dh += j * h[:, j - 1] * ds[k - j + 1]
            h[:, k - 1] = bordered(J - k * lam * np.eye(n), v, w, dh - f[:, k])
        f = _series_drift(real, q, h, over, size)
    except NotAnalytic as error:
        return {"covered": False, "reason": str(error)}
    except (TypeError, AttributeError, ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        return {"covered": False, "reason": f"the drift is not written in operations of a power series ({error})"}
    a = w @ f
    if not np.all(np.isfinite(a)):
        return {"covered": False, "reason": "a coefficient of the drift at the write point is not finite"}
    return {"covered": True, "order": order, "coefficients": [float(x) for x in a],
            "rounding": [float(1e-13 * x) for x in np.max(np.abs(f), axis=0)]}


def series_reach(reduced: Dict[str, object], half: float, rate: float) -> float:
    """The amplitude, at most `half`, up to which the series describes the drift: where a term of order n >= 4 gives a
    rate a_n s^(n - 1) as large as `rate`, that of the slowest relaxing direction (the largest rate where there is none),
    the reduction to one coordinate is no longer accurate, nor is a term of the series a measure of the drift."""
    a = reduced["coefficients"]
    return min([half] + [(rate / abs(a[n])) ** (1.0 / (n - 1)) for n in range(4, len(a)) if a[n] != 0.0])


def leading_series_term(reduced: Dict[str, object], half: float, scale: float) -> Optional[Dict[str, float]]:
    """The leading term c s^n of a reduced drift whose quadratic and cubic terms are not resolved at the write point: the
    term of lowest order from the quartic on that, at the amplitude `half`, is within a tenth of the largest of them,
    among those above rounding (1e-13 of the largest rate, or of the terms of the drift it is a sum of). The tenth is
    the tolerance of the asymmetry: a write point located a distance d along the mode from the zero of c s^n has the
    lower terms of c (s + d)^n, which are not a term of the material. None where no term up to the order of the series
    is resolved: the drift vanishes along the critical direction to that order."""
    a, rounding = reduced["coefficients"], reduced["rounding"]
    resolved = {n: abs(a[n]) * half ** n for n in range(4, len(a))
                if abs(a[n]) > rounding[n] + 1e-13 * scale * half ** (1 - n)}
    if not resolved:
        return None
    top = max(resolved.values())
    n = min(m for m, t in resolved.items() if t >= 0.1 * top)
    return {"order": n, "coefficient": float(a[n])}


def slow_spectrum(ev: np.ndarray, cand: np.ndarray, lam: complex, slope: float, resolution: float, scale: float,
                  J: np.ndarray, R: np.ndarray) -> Dict[str, object]:
    """What the record needs beyond the count of critical directions: the rates of the other directions at the write
    point (`spectrum`, real and imaginary parts sorted by the size of the real part, the eight slowest), the rate of
    the slowest direction that is not critical (`second_rate`, the splitting of the critical rates), and the rank of
    the Jacobian restricted to the critical eigenspace (`critical_block_rank`: zero for a semisimple degeneracy, as
    the equivariant one of a ring, k - 1 for a nilpotent one, Bogdanov-Takens). A direction whose rate is below the
    square root of the sweep rate of a write check relaxes less than once while the choice is made, so the write then
    takes place in more directions than the continuation counts as critical: completion.law_at counts them."""
    crit = critical_directions(ev, cand, lam, slope, resolution, scale)
    others = [j for j in cand if j not in set(crit.tolist())]
    order = sorted(cand.tolist(), key=lambda j: abs(ev[j].real))
    out: Dict[str, object] = {
        "spectrum": [[float(ev[j].real), float(ev[j].imag)] for j in order[:8]],
        "second_rate": float(min(abs(ev[j].real) for j in others)) if others else None,
    }
    k = int(crit.size)
    if k >= 2:
        B = np.real(R[:, crit])
        Bp = np.linalg.pinv(B)
        block = Bp @ J @ B
        out["critical_block_rank"] = int(np.linalg.matrix_rank(block, tol=1e-6 * max(scale, 1e-300)))
    return out


SLOW_COUPLING_WINDOW = 3.0   # the passage window, in units of sqrt(rate) / |kappa_slope| on each side of the crossing
SLOW_COUPLING_TOL = 1e-6     # of max(|J|, sqrt(rate)): the accuracy of the finite-difference Jacobian


def _unwritten_branch(real: Realization, q: np.ndarray, over: Dict, steps: int = 30) -> Optional[np.ndarray]:
    """The unwritten branch at `over`, by Newton from q (unstable after the crossing, which Newton does not mind)."""
    q = np.array(q, float)
    for _ in range(steps):
        f = np.asarray(real.F(q[None], **over))[0]
        if not np.all(np.isfinite(f)):
            return None
        try:
            dq = np.linalg.lstsq(jacobian(real, q, over), -f, rcond=None)[0]
        except np.linalg.LinAlgError:
            return None
        q = q + dq
        if np.linalg.norm(dq) < 1e-12 * max(1.0, np.linalg.norm(q)):
            break
    return q


def slow_coupling(real: Realization, nf: Dict, rate: float, factor: float = 1.0) -> Dict[str, object]:
    """The slow directions at a write point (rate below factor x sqrt(rate)) and, for each, its largest linear coupling
    to the critical space over the passage window, in the eigenbasis at the write point, along the unwritten branch.

    A slow direction that is not linearly coupled to the critical space anywhere in the window does not take part in
    the choice: the critical coordinate then obeys a closed one-dimensional linear equation through the window, and
    the slow variable enters only at second order, as an outside modulation of the rate or the push, which the law of
    one direction carries. A coupled slow direction is part of the choice, and so is one whose own rate reaches zero
    along the branch inside the window: it is then a second critical direction of the same passage, not a modulator."""
    k_static = max(1, int(nf.get("critical_dimension", 1) or 1))
    q_c = np.asarray(nf["state"], float)
    param, v_c = nf["param"], float(nf["value"])
    J0 = jacobian(real, q_c, {param: v_c})
    ev, V = np.linalg.eig(J0)
    order = np.argsort(np.abs(ev.real))
    crit = list(order[:k_static])
    slow = [int(j) for j in order[k_static:] if abs(ev[j].real) <= factor * sqrt(rate)]
    out: Dict[str, object] = {"slow": [], "coupled": 0}
    if not slow:
        return out
    W = np.linalg.inv(V)
    g = abs(float(nf.get("kappa_slope", 0.0) or 0.0))
    width = SLOW_COUPLING_WINDOW * sqrt(rate) / g if g > 0 else 0.0
    scale = max(float(np.max(np.abs(J0))), sqrt(rate))
    worst = {j: 0.0 for j in slow}
    top = {j: -np.inf for j in slow}
    for v in np.linspace(v_c - width, v_c + width, 13):
        over = {param: float(v)}
        qb = _unwritten_branch(real, q_c, over)
        if qb is None:
            continue
        M = W @ jacobian(real, qb, over) @ V
        for j in slow:
            worst[j] = max(worst[j], float(np.max(np.abs(M[crit, j]))), float(np.max(np.abs(M[j, crit]))))
            top[j] = max(top[j], float(M[j, j].real))
    for j in slow:
        c, crosses = worst[j] / scale, bool(top[j] >= -SLOW_COUPLING_TOL * scale)
        out["slow"].append({"rate": float(ev[j].real), "coupling": c, "largest_rate_in_window": float(top[j]),
                            "crosses": crosses, "coupled": bool(c > SLOW_COUPLING_TOL or crosses)})
    out["coupled"] = sum(1 for x in out["slow"] if x["coupled"])
    return out


def dynamical_critical_dimension(nf: Dict, rate: float, factor: float = 1.0, real: Optional[Realization] = None) -> int:
    """The number of directions that take part in a write swept at `rate`: the critical ones and those whose rate at
    the write point is below factor x sqrt(rate), which relax less than once during the passage (the choice is made
    over a time of order 1/sqrt(rate) around the crossing).

    With the realization, a slow direction counts only if it is linearly coupled to the critical space somewhere in
    the passage window (slow_coupling); an independent slow modulator does not make the write two-dimensional. Without
    it, every slow direction counts."""
    k = int(nf.get("critical_dimension", 1) or 1)
    spec = nf.get("spectrum") or []
    if not spec:
        return k
    slow = sum(1 for re, _ in spec if abs(re) <= factor * sqrt(rate))
    if real is None or slow <= k:
        return max(k, int(slow))
    return max(k, max(1, k) + int(slow_coupling(real, nf, rate, factor)["coupled"]))


def critical_directions(ev: np.ndarray, cand: np.ndarray, lam: complex, slope: float, resolution: float,
                        scale: float) -> np.ndarray:
    """Indices of the eigenvalues ev (among cand, the directions outside the declared null modes) whose real part is
    zero at a write point to the accuracy of the continuation. The continuation stops where the critical eigenvalue
    lam is still |Re lam| from zero, and it does not separate write points closer than its resolution in the
    parameter, over which lam moves by |slope| x resolution; below 1e-7 of the largest rate the Jacobian itself is
    not accurate. An eigenvalue within that distance of zero crosses at the same write point."""
    tol = abs(lam.real) + max(abs(slope) * resolution, 1e-7 * scale)
    return cand[np.abs(ev.real[cand]) <= tol]


def _relaxes_on_a_side(real: Realization, q: np.ndarray, param: str, v: float, width: float, scale: float) -> bool:
    """Whether the state continues to one side of v along the control as a state with no growing direction and at
    least one relaxing direction, within a tenth of the base length."""
    lo, hi = real.control_range if real.control_range is not None else (-np.inf, np.inf)
    for w in (v - 0.01 * width, v + 0.01 * width):
        if not lo <= w <= hi:
            continue
        over = {param: w}
        qn, ok = newton(real, q, over)
        if not ok or real.carrier.distance(q, qn) > 0.1 * base_length(real):
            continue
        re = np.linalg.eigvals(jacobian(real, qn, over)).real
        if re.max() <= 1e-7 * scale and re.min() < -1e-7 * scale:
            return True
    return False


def _on_one_family(real: Realization, a: np.ndarray, b: np.ndarray, over: Dict, nodes: int = 24) -> bool:
    """Whether two states at the same control are joined by a connected set of equilibria: the nodes of the straight
    path between them, each corrected by Newton (minimal norm, so a node already on the set stays), stay equilibria and
    keep their spacing. Mirror images have an unstable or isolated state between them, and the nodes jump."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = b - a
    if real.carrier.kind == "torus":
        d = (d + real.carrier.period / 2) % real.carrier.period - real.carrier.period / 2
    step = float(np.linalg.norm(d)) / (nodes - 1)
    last = a
    for t in np.linspace(0.0, 1.0, nodes)[1:]:
        q, ok = newton(real, a + t * d, over)
        if not ok or real.carrier.distance(last, q) > 3 * step + 1e-9:
            return False
        last = q
    return True


def distinct_write_points(real: Realization, events: List[Dict]) -> Tuple[List[Dict], List[Dict]]:
    """The events of locate_writes that are distinct places where a state loses stability, and those set aside.

    A state with no relaxing direction that does not continue along the control as a state with one is set aside: at
    the ends of the control range of a dissipationless free top the drift vanishes on a whole plane, nothing there
    is stable, and nothing loses stability. Random starts sample a continuum of equilibria at one control value many
    times; the samples count as one event, which keeps the number of them as `samples`."""
    param = real.control
    width = abs(real.control_range[1] - real.control_range[0]) if real.control_range is not None else 1.0
    kept: List[Dict] = []
    aside: List[Dict] = []
    for ev in events:
        q, v = np.asarray(ev["q"], float), float(ev["v"])
        over = {param: v}
        spectrum = np.linalg.eigvals(jacobian(real, q, over))
        re = spectrum.real
        scale = max(1e-12, float(np.abs(spectrum).max()))
        if re.min() >= -1e-7 * scale and not _relaxes_on_a_side(real, q, param, v, width, scale):
            aside.append({**ev, "reason": "no relaxing direction on either side of the control value"})
            continue
        for e in kept:
            if abs(float(e["v"]) - v) < WRITE_RESOLUTION * width and _on_one_family(real, np.asarray(e["q"], float), q, over):
                e["samples"] = int(e.get("samples", 1)) + 1
                break
        else:
            kept.append(dict(ev))
    return kept, aside


def _without_neutral(real: Realization, q: np.ndarray, param: str, v: float, ev: np.ndarray, R: np.ndarray,
                     cand: np.ndarray, over: Dict) -> np.ndarray:
    """The candidates without the directions that are neutral on both sides of v along the control while another
    candidate crosses there. A conserved quantity (the spin length of an undamped atom) or a line of states keeps a zero
    eigenvalue along the whole branch; counted with the eigenvalues that cross, it reads as a degenerate critical mode.
    Where nothing else is near zero, or the state does not continue to both sides, the candidates are returned as they
    are."""
    width = (abs(real.control_range[1] - real.control_range[0])
             if real.control_range is not None and param == real.control else max(1.0, abs(v)))
    neutral: Optional[np.ndarray] = None
    moved = np.inf
    for w in (v - 5e-3 * width, v + 5e-3 * width):
        side = {**over, param: w}
        qn, ok = newton(real, q, side)
        if not ok or real.carrier.distance(q, qn) > 0.1 * base_length(real):
            return cand
        ev_side, R_side = np.linalg.eig(jacobian(real, qn, side))
        # zero against the rates the declared step resolves, at most 0.3 / dt (correction 27 of the discovery engine): at
        # a slaved rate of 1e6 the bound was 0.1, and the eigenvalue that crosses read as neutral on both sides
        zero = np.abs(ev_side.real) <= 1e-7 * max(1e-12, min(float(np.abs(ev_side).max()), 0.3 / real.dt))
        if not zero.all():
            moved = min(moved, float(np.abs(ev_side.real[~zero]).min()))
        if neutral is None or zero.sum() < neutral.shape[1]:
            neutral = R_side[:, zero]
    m = neutral.shape[1]
    if m == 0 or not np.isfinite(moved):
        return cand
    # the crossing eigenvalue is within what it moves over this step of the control
    near = cand[np.abs(ev.real[cand]) <= moved]
    if m >= near.size:
        return cand
    basis = np.linalg.qr(neutral)[0]
    overlap = [float(np.linalg.norm(basis.conj().T @ R[:, j]) / np.linalg.norm(R[:, j])) for j in near]
    return np.setdiff1d(cand, near[np.argsort(overlap, kind="stable")[::-1][:m]])


def normal_form(real: Realization, q: np.ndarray, param: str, v: float, extra: Optional[Dict] = None,
                npts: int = 41, half: Optional[float] = None) -> Dict[str, object]:
    """Reduce the drift at a write point to one coordinate s along the critical direction:
    ds/dt = a0 + a1 s + a2 s^2 + a3 s^3 (on the slow manifold), with the unfolding in the parameter
    (bias = the part of dF/dparam along the critical direction, kappa slope = d(a1)/dparam).

    The reduction holds when one real eigenvalue crosses zero (a complex pair at a Hopf point). critical_dimension
    counts the eigenvalues that cross together (critical_directions); when a symmetry makes the critical mode
    degenerate, as the twisted modes of a ring of an odd number of cells, it exceeds one, the critical eigenspace is
    not a line, and the kind is a degenerate critical mode. The coefficients along one direction of that eigenspace
    are kept as a diagnostic.

    Where one direction is critical, the coefficients a0, a1, a2, a3 and the terms up to SERIES_ORDER are those of the
    series closed order by order at the write point (series_reduction, in the coordinate s = w.(x - q)), and the kind is
    read from them, with the amplitude of the comparison (asymmetry, bias) set by the default window and series_reach
    and not by a fit. Where the series does not exist or is not covered (several critical directions, abs at a kink), the
    coefficients are those of the cubic fitted to the drift on the slow manifold, and `series` says why.

    That drift is sampled in either case, for the record of its shape (s, f and the coefficients `fit`), over s in
    [-half, half], by default a quarter of the base length (an eighth of the period on a torus), and the window is halved
    further where the cubic does not fit."""
    q = np.asarray(q, float)
    over = {**(extra or {}), param: v}
    n = q.size
    J = jacobian(real, q, over)
    ev, R = np.linalg.eig(J)
    cand = np.arange(n)
    if real.null_modes is not None:
        Nm = np.atleast_2d(np.asarray(real.null_modes, float))
        Nm = Nm / np.linalg.norm(Nm, axis=1, keepdims=True)
        overl = np.array([np.linalg.norm(Nm @ (R[:, i] / np.linalg.norm(R[:, i]))) for i in range(n)])
        if np.any(overl < 0.9):
            cand = np.flatnonzero(overl < 0.9)
    cand = _without_neutral(real, q, param, v, ev, R, cand, over)
    i = int(cand[np.argmin(np.abs(ev.real[cand]))])
    lam = ev[i]
    scale = max(1e-12, float(np.abs(ev).max()))
    width = (abs(real.control_range[1] - real.control_range[0])
             if real.control_range is not None and param == real.control else max(1.0, abs(v)))
    dv = 1e-5 * width
    # the derivatives in the control are central differences over v -+ dv. Where the drift is not defined on one side, as
    # a negative variable under a power whose exponent is the control (correction 28), they are taken on the other side;
    # where it is defined on neither, they are not defined, and the point is named by its drift along the critical
    # direction alone. Where every value is finite nothing changes
    Fp, Fm = real.F(q, **{**over, param: v + dv}), real.F(q, **{**over, param: v - dv})
    dF = (Fp - Fm) / (2 * dv)
    if not np.all(np.isfinite(dF)):
        F0 = real.F(q, **over)
        dF = ((Fp - F0) / dv if np.all(np.isfinite(Fp)) else (F0 - Fm) / dv if np.all(np.isfinite(Fm))
              else np.zeros_like(F0))
    out: Dict[str, object] = {"param": param, "value": float(v), "state": q.tolist(),
                              "kappa_critical": [float(-lam.real), float(-lam.imag)]}

    def critical(slope: float) -> Dict[str, object]:
        crit = critical_directions(ev, cand, lam, slope, WRITE_RESOLUTION * width, scale)
        return {"critical_dimension": int(crit.size),
                "critical_kappa": [[float(-ev[j].real), float(-ev[j].imag)] for j in crit]}

    # over the resolution of the write point the critical eigenvalue moves by |slope| x resolution; a frequency below
    # that cannot be told from two real eigenvalues crossing together (where a degenerate pair is not diagonal in the
    # coordinates, rounding gives it imaginary parts of order 1e-13, and when the whole spectrum is critical, 1e-6 of
    # the largest rate is no bound)
    Jn, side = jacobian(real, q, {**over, param: v + dv}), 1.0
    if not np.all(np.isfinite(Jn)):
        Jn, side = jacobian(real, q, {**over, param: v - dv}), -1.0
    if np.all(np.isfinite(Jn)):
        ev_next = np.linalg.eigvals(Jn)
        slope_re = float(side * (ev_next[np.argmin(np.abs(ev_next - lam))].real - lam.real) / dv)
    else:
        slope_re = 0.0
    if abs(lam.imag) > max(1e-6 * scale, abs(slope_re) * WRITE_RESOLUTION * width):
        # a complex pair: two directions are critical at a simple Hopf point
        out.update(kind="Hopf: the unwritten state starts to oscillate", hopf=True, omega=float(abs(lam.imag)),
                   asymmetry=0.0, bias=0.0, a3=0.0, **critical(0.0))
        return out
    no_dual = False
    if lam.imag != 0.0:
        # a pair below the resolution crosses as two real directions. Its eigenvectors from eig are complex with an
        # arbitrary phase, and their real parts alone depend on that phase (with NumPy 2 the left and right real parts
        # came out orthogonal). The plane of the pair is real, spanned by the real and imaginary parts of its right
        # eigenvector, and the plane of its left eigenvector is dual to it: the reduction is taken along one direction
        # of the first, with the left vector of the second that is one on it and zero on the other direction.
        r, ell = R[:, i], np.linalg.inv(R)[i]
        P = sorted([np.real(r), np.imag(r)], key=np.linalg.norm, reverse=True)
        vec = P[0] / np.linalg.norm(P[0])
        vec2 = P[1] - np.dot(P[1], vec) * vec
        Lp = [np.real(ell), np.imag(ell)]
        M = np.array([[np.dot(Lp[0], vec), np.dot(Lp[1], vec)], [np.dot(Lp[0], vec2), np.dot(Lp[1], vec2)]])
        a, b = np.linalg.solve(M, [1.0, 0.0])
        w = a * Lp[0] + b * Lp[1]
    else:
        vec = np.real(R[:, i])
        vec /= np.linalg.norm(vec)
        evl, L = np.linalg.eig(J.T)
        w = np.real(L[:, int(np.argmin(np.abs(evl - lam)))])
        if abs(np.dot(w, vec)) <= 1e-8 * np.linalg.norm(w):
            # a degenerate eigenvalue: eig chooses the bases of the right and the left eigenspace independently, and
            # the left vector of the same eigenvalue can be orthogonal to the right one (at the write point of the
            # two-component model with the symmetry of a square, where the Jacobian is eps times the identity to
            # rounding, it was for some seeds; the coefficients became NaN and the point read as a subcritical
            # pitchfork with no critical direction). The row of R^-1 is the left vector dual to vec: one on it and
            # zero on the other eigenvectors.
            # Where the eigenvalue is defective (a Jordan block, as x' = y, y' = k x at k = 0; a free rigid body at
            # B = A) there are fewer eigenvectors than repeated eigenvalues, R is singular, and no left vector pairs
            # with the right one: the critical eigenspace is not a line, and there are no coefficients along one.
            same = np.flatnonzero(np.abs(ev - lam) <= 1e-6 * max(1.0, float(np.abs(J).max())))
            rank = int(np.linalg.matrix_rank(J - lam.real * np.eye(n), tol=1e-9 * max(1.0, float(np.abs(J).max()))))
            no_dual = n - rank < same.size
            if not no_dual:
                try:
                    w = np.real(np.linalg.inv(R)[i])
                except np.linalg.LinAlgError:
                    no_dual = True
        with np.errstate(all="ignore"):
            w = np.full(n, np.nan) if no_dual else w / np.dot(w, vec)
        no_dual = no_dual or not np.all(np.isfinite(w))
    # the sign of an eigenvector is arbitrary and differs between LAPACK builds: orient the critical direction (see
    # orientation_index). +mode then names the same state on every machine, and the random streams of a test
    # assigned to the bias toward +mode drive the same state.
    if vec[orientation_index(vec)] < 0:
        vec, w = -vec, -w
    if half is None:
        half = 0.25 * base_length(real) if real.carrier.kind != "torus" else real.carrier.period / 8
    if real.carrier.kind == "orthant":
        nz = np.abs(vec) > 1e-9
        room = np.abs(q)[nz] / np.abs(vec)[nz]
        if room.size and np.min(room) > 1e-9:
            half = min(half, 0.5 * float(np.min(room)))
    half0 = half
    for _ in range(8):
        s = np.linspace(-1.0, 1.0, npts) * half
        f, transverse = _slow_manifold_drift(real, q, vec, w, s, over)
        coef = np.polyfit(s, f, 3)
        misfit = float(np.max(np.abs(np.polyval(coef, s) - f))) / max(float(np.max(np.abs(f))), 1e-300)
        if np.all(np.isfinite(f)) and transverse < 1e-6 and misfit < 0.05:
            break
        half /= 2
    a3, a2, a1, a0 = (float(c) for c in coef)
    fit = {"a0": a0, "a1": a1, "a2": a2, "a3": a3}
    Jp = jacobian(real, q, {**over, param: v + dv})
    Jm = jacobian(real, q, {**over, param: v - dv})
    fp, fm = bool(np.all(np.isfinite(Jp))), bool(np.all(np.isfinite(Jm)))
    dJ = ((Jp - Jm) / (2 * dv) if fp and fm else (Jp - J) / dv if fp else (J - Jm) / dv if fm else np.zeros_like(J))
    growth_partial = float(w @ dJ @ vec)
    # along the branch the state moves too, dq/dv = -J^+ dF/dv (no component along the critical direction)
    dq = -np.linalg.lstsq(J, dF, rcond=None)[0]
    dq = dq - vec * float(w @ dq)
    eq = 1e-5 * base_length(real)
    nq = float(np.linalg.norm(dq))
    if nq > 0:
        DJ = (jacobian(real, q + eq * dq / nq, over) - jacobian(real, q - eq * dq / nq, over)) / (2 * eq) * nq
        growth = growth_partial + float(w @ DJ @ vec)
    else:
        growth = growth_partial
    bias_along = float(w @ dF)
    quad, cub = abs(a2) * half ** 2, abs(a3) * half ** 3
    asym = quad / (quad + cub) if quad + cub > 1e-12 * scale * half else 0.0
    # at a fold the control pushes along the critical direction and the eigenvalue goes as the square root of the
    # distance from the point: the growth along the branch diverges, is no rate of crossing, and took in every eigenvalue
    # of a cavity magnon and of a phase qubit. There the rate is the part that does not follow the state.
    fold = asym > 0.1 and abs(bias_along) > max(0.1 * abs(growth_partial) * half, 1e-9 * scale)
    crossing = growth_partial if fold else growth
    if no_dual:
        same = np.flatnonzero(np.abs(ev - lam) <= 1e-6 * max(1.0, float(np.abs(J).max())))
        crit = {"critical_dimension": int(same.size),
                "critical_kappa": [[float(-ev[j].real), float(-ev[j].imag)] for j in same]}
        rank = int(np.linalg.matrix_rank(J - lam.real * np.eye(n), tol=1e-9 * max(1.0, float(np.abs(J).max()))))
        out["defective"] = bool(n - rank < same.size)  # fewer eigenvectors than repeated eigenvalues
    else:
        crit = critical(crossing)
    k = crit["critical_dimension"]
    crit.update(slow_spectrum(ev, cand, lam, crossing, WRITE_RESOLUTION * width, scale, J, R))
    series = (series_reduction(real, q, vec, w, over) if k == 1 else
              {"covered": False, "reason": f"{k} critical directions: the series is closed along one"})
    measure, window = half, half
    if series["covered"]:
        a0, a1, a2, a3 = series["coefficients"][:4]
        out.update({f"a{j}": c for j, c in enumerate(series["coefficients"]) if j > 3})
        measure, window = series_reach(series, half0, crit["second_rate"] or scale), half0
    bias = abs(bias_along) / (abs(bias_along) + abs(growth) * measure + 1e-300)
    quad, cub = abs(a2) * measure ** 2, abs(a3) * measure ** 3
    asym = quad / (quad + cub) if quad + cub > 1e-12 * scale * measure else 0.0
    fold = asym > 0.1 and abs(bias_along) > max(0.1 * abs(growth_partial) * measure, 1e-9 * scale)
    if k > 1:
        kind = (f"degenerate critical mode: {k} directions become critical together; the critical eigenspace is "
                f"{k}-dimensional, not a line")
    elif abs(a2) * window ** 2 + abs(a3) * window ** 3 <= 1e-6 * min(scale, 0.3 / real.dt) * window:
        # neither the quadratic nor the cubic term is resolved against what the rates of the model make over the window.
        # Those are the rates the declared step resolves, at most 0.3 / dt, as in safe_dt (correction 27 of the discovery
        # engine): a variable slaved at a faster rate is eliminated by the reduction to the slow manifold, and at the rate
        # 1e6 it put the bound at the order of the reduced drift itself, so that every fold and pitchfork with a fast
        # variable was read as degenerate. Where neither term is resolved, the leading term is of higher order (a Kerr
        # parametric oscillator at zero detuning saturates at fifth order)
        lead = (leading_series_term(series, measure, scale) if series["covered"] else
                leading_order(real, q, vec, w, over, half0, float(lam.real), scale))
        if lead is None:
            kind = "degenerate: the drift vanishes along the critical direction (a line of states)"
        else:
            out["leading_term"] = lead
            kind = (f"degenerate cubic: the quadratic and cubic terms vanish at the write point; the leading term along "
                    f"the critical direction is of order {lead['order']}")
    elif asym > 0.1 and (fold or bias > 0.1):
        kind = "saddle-node (fold): a state appears or vanishes; written by a pulse past the threshold"
    elif asym > 0.1:
        kind = "transcritical: two states exchange stability; the write is one-sided"
    elif bias > 0.1:
        kind = "imperfect pitchfork: the parameter itself biases the choice"
    elif a3 < 0:
        kind = "supercritical pitchfork: symmetric write; a small bias chooses the state"
    else:
        kind = "subcritical pitchfork: symmetric, but the state jumps to a distant branch"
    out.update(hopf=False, kind=kind, mode=vec.tolist(), left_mode=w.tolist(), a0=a0, a1=a1, a2=a2, a3=a3,
               bias_along_mode=bias_along, bias=bias, kappa_slope=-growth, kappa_slope_partial=-growth_partial,
               asymmetry=asym,
               transverse_residual=transverse, cubic_misfit=misfit, s=s.tolist(), f=f.tolist(), half=half, fit=fit,
               series=series, **crit)
    return out


def with_write_field(real: Realization, target: np.ndarray, h_range) -> Realization:
    """The realization with a bounded write toward target added to its drift; the field WRITE_FIELD is the control."""
    target = np.asarray(target, float)
    base = real

    def drift(q, p):
        return base.drift(q, p) + base.write_force(q, target, p[WRITE_FIELD])

    r = replace(real, name=real.name + " + write field", drift=drift, potential=None,
                params={**real.params, WRITE_FIELD: 0.0}, control=WRITE_FIELD, control_range=tuple(h_range))
    for attr in ("graph", "linear", "loop"):
        if hasattr(real, attr):
            setattr(r, attr, getattr(real, attr))
    return r


def branch_summary(real: Realization, q: np.ndarray) -> float:
    """One number per fixed point for a bifurcation diagram: q_1 in one dimension, q_1 - mean(q_rest)
    otherwise (separates mirror states such as u-high/v-high), and the global order |mean e^{i m q}| on a torus."""
    q = np.asarray(q, float)
    if real.carrier.kind == "torus":
        return float(np.abs(np.mean(np.exp(1j * real.m * q))))
    return float(q[0] if q.size == 1 else q[0] - np.mean(q[1:]))


def bifurcation_scan(real: Realization, rng, values: Sequence[float], n_starts: int = 30, seeds=None):
    rows = []
    for v in values:
        fps = fixed_points(real, rng, n_starts, {real.control: v}, seeds=seeds)
        rows.append({"value": float(v), "points": [
            {"q": q.tolist(), "kappa_min": float(k.real.min()), "stable": n_unstable(k) == 0,
             "summary": branch_summary(real, q), "mean": float(np.mean(q))} for q, k in fps]})
    return rows


def is_scale_control(real: Realization, rng=None, n: int = 8) -> bool:
    """True if the drift is proportional to the control (a coupling scale, or an inverse temperature when the
    noise is fixed): then the control changes no state, only the depth of the landscape relative to the noise."""
    if not real.control:
        return False
    rng = rng if rng is not None else np.random.default_rng(0)
    q = real.carrier.sample(rng, n)
    v = float(real.params[real.control]) or 1.0
    F1 = real.F(q, **{real.control: v})
    F2 = real.F(q, **{real.control: 2.0 * v})
    return bool(np.max(np.abs(F2 - 2.0 * F1)) <= 1e-9 * max(1.0, float(np.max(np.abs(F1)))))


def memory_kernel(real: Realization, t: np.ndarray) -> Optional[np.ndarray]:
    """For a linear realization with blocks [[A, B], [C, D]] (observed, unobserved): K(t) = B e^{Dt} C."""
    lin = getattr(real, "linear", None)
    if lin is None:
        return None
    from scipy.linalg import expm
    return np.array([(lin["B"] @ expm(lin["D"] * ti) @ lin["C"]).item() for ti in t])
