"""Measure how a realization writes, retains, is observed and rewritten, and where its writes happen.

  stored_states      stable states reached by relaxation from random starts, with basin counts
  kappa_spectrum     kappa = -eig(Jacobian): > 0 relaxing, ~ 0 flat (symmetry), < 0 unstable (writes)
  loss_law           the law by which a stored state is lost, read from the spectra and the number of states
  neb_barrier        barrier between two stored states along the minimum-energy path (gradient systems)
  hold_time          Monte Carlo time for a stored state to lose its identity at noise D
  write_test         bounded write toward a target at fixed control or while sweeping the control
                     from the end of its range where the barriers are low; accuracy after release and the
                     rewriting time
  fixed_points       zeros of the drift (Newton from random starts), stable and unstable, with kappa spectra
  track              continuation of a stable state along a parameter until it loses stability or vanishes
  locate_writes      parameter values at which stored states lose stability or appear (write points)
  normal_form        drift along the critical direction at a write point, reduced onto the slow manifold
  with_write_field   the realization plus a bounded write toward a target, with the field h as parameter
  memory_kernel      for a linear realization with an observed block, K(t) = B e^{Dt} C (Mori-Zwanzig)
"""
from __future__ import annotations

from dataclasses import replace
from typing import Dict, List, Optional, Sequence

import numpy as np

from .identity import Realization

P_STAR = 0.9


def base_length(real: Realization) -> float:
    return real.carrier.period / (2 * np.pi) if real.carrier.kind == "torus" else real.carrier.scale


def _merge_tol(real: Realization) -> float:
    return 0.05 * base_length(real) * np.sqrt(real.carrier.dim)


def _canon(real: Realization, q: np.ndarray) -> np.ndarray:
    return real.canonical(q) if real.canonical is not None else q


# ------------------------------------------------------------------------------------------------ dynamics
def safe_dt(real: Realization) -> float:
    """The declared time step, reduced where the fastest local rate would make an explicit Euler step unstable
    (dt <= 0.3 / max |eigenvalue| over sampled states and the control range)."""
    cached = getattr(real, "_safe_dt", None)
    if cached is not None:
        return cached
    rng = np.random.default_rng(12345)
    rate = 0.0
    overs = [{}]
    if real.control and real.control_range:
        overs += [{real.control: real.control_range[0]}, {real.control: real.control_range[1]}]
    for over in overs:
        for q in real.carrier.sample(rng, 6):
            J = jacobian(real, q, over)
            if np.all(np.isfinite(J)):
                rate = max(rate, float(np.max(np.abs(np.linalg.eigvals(J)))))
    dt = float(real.dt if rate <= 0 else min(real.dt, 0.3 / rate))
    real._safe_dt = dt
    return dt


def relax(real: Realization, q0: np.ndarray, over: Optional[Dict] = None, max_time: float = 300.0,
          tol: float = 1e-5, dt: Optional[float] = None):
    over, dt = over or {}, dt or safe_dt(real)
    q = real.carrier.wrap(np.atleast_2d(np.array(q0, float)))
    for k in range(int(max_time / dt)):
        f = real.F(q, **over)
        q = real.carrier.wrap(q + dt * f)
        if k % 50 == 0 and np.all(np.linalg.norm(f, axis=-1) < tol):
            break
    return q, np.linalg.norm(real.F(q, **over), axis=-1)


def integrate(real: Realization, q: np.ndarray, T: float, D: float, rng, over: Optional[Dict] = None,
              write=None, control_path=None, record=None, n_record: int = 40, dt: Optional[float] = None):
    """Euler-Maruyama; write = (target, h); control_path(t) sets real.control; record(q) -> value."""
    over, dt = dict(over or {}), dt or safe_dt(real)
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
        q = real.carrier.wrap(q + dt * f + amp * rng.standard_normal(q.shape))
        if record is not None and (k + 1) % every == 0:
            out_t.append((k + 1) * dt)
            out_v.append(record(q))
    return q, np.array(out_t), np.array(out_v)


# ------------------------------------------------------------------------------------------------ spectra
def jacobian(real: Realization, q: np.ndarray, over: Optional[Dict] = None, eps: float = 1e-6) -> np.ndarray:
    over, q = over or {}, np.asarray(q, float)
    n = q.size
    E = eps * np.eye(n)
    Fp = real.F(q[None] + E, **over)
    Fm = real.F(q[None] - E, **over)
    return ((Fp - Fm) / (2 * eps)).T


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
    settle, and whether the states form a continuum (flat directions beyond the declared symmetries)."""
    q0 = real.carrier.sample(rng, n_starts)
    q, speed = relax(real, q0, over)
    conv = speed < 1e-4
    states: List[np.ndarray] = []
    counts: List[int] = []
    tol = _merge_tol(real)
    for qi in q[conv]:
        ci = _canon(real, qi)
        hit = next((k for k, s in enumerate(states) if real.carrier.distance(ci, _canon(real, s)) < tol), None)
        if hit is None:
            states.append(qi)
            counts.append(1)
        else:
            counts[hit] += 1
    order = np.argsort(counts, kind="stable")[::-1]
    states = [states[i] for i in order]
    counts = [counts[i] for i in order]
    spectra = [kappa_spectrum(real, s, over) for s in states]
    keep = [i for i, k in enumerate(spectra) if classify(k)["unstable"] == 0]
    states, counts, spectra = [states[i] for i in keep], [counts[i] for i in keep], [spectra[i] for i in keep]
    n_null = 0 if real.null_modes is None else len(np.atleast_2d(real.null_modes))
    continuum = bool(len(states) >= max(6, int(conv.sum()) // 2)
                     and all(classify(k)["flat"] > n_null for k in spectra))
    return states, counts, spectra, int(np.sum(~conv)), continuum


def loss_law(n_states: int, spectra: Sequence[np.ndarray], n_unconverged: int, n_starts: int,
             continuum: bool = False, n_null: int = 0) -> str:
    if n_states == 0 and n_unconverged > 0:
        return "no stable state: a limit cycle, whose phase diffuses (Law 2 along the cycle)"
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


def newton(real: Realization, q: np.ndarray, over: Optional[Dict] = None, iters: int = 80, tol: float = 1e-11):
    over = over or {}
    q = real.carrier.wrap(np.array(q, float))
    cap = 0.5 * base_length(real)
    for _ in range(iters):
        f = real.F(q, **over)
        if np.linalg.norm(f) < tol:
            return q, True
        step = np.linalg.lstsq(jacobian(real, q, over), f, rcond=None)[0]
        n = np.linalg.norm(step)
        if n > cap:
            step *= cap / n
        q = real.carrier.wrap(q - step)
    return q, bool(np.linalg.norm(real.F(q, **over)) < 1e-8)


def fixed_points(real: Realization, rng, n_starts: int = 40, over: Optional[Dict] = None, seeds=None):
    """Zeros of the drift by damped Newton from random starts (and from given seeds, such as stored states),
    stable and unstable alike."""
    over = over or {}
    found: List[np.ndarray] = []
    tol = _merge_tol(real)
    starts = real.carrier.sample(rng, n_starts)
    if seeds is not None and len(seeds):
        starts = np.concatenate([np.atleast_2d(np.asarray(seeds, float)), starts])
    for q in starts:
        q, ok = newton(real, q, over)
        if ok and not any(real.carrier.distance(_canon(real, q), _canon(real, s)) < tol for s in found):
            found.append(q)
    return [(q, kappa_spectrum(real, q, over)) for q in found]


# ------------------------------------------------------------------------------------------------ barrier
def neb_barrier(real: Realization, a: np.ndarray, b: np.ndarray, over: Optional[Dict] = None,
                images: int = 24, iters: int = 4000, k_spring: float = 1.0) -> Dict[str, object]:
    over = over or {}
    step = 0.5 * safe_dt(real)
    t = np.linspace(0, 1, images)[:, None]
    path = real.carrier.wrap(a[None] + t * real.carrier.diff(a, b)[None])
    for it in range(iters):
        E = real.V(path, **over)
        F = real.F(path, **over)
        climb = int(np.argmax(E[1:-1])) + 1
        new = path.copy()
        for i in range(1, images - 1):
            tau = real.carrier.diff(path[i - 1], path[i + 1])
            tau /= max(np.linalg.norm(tau), 1e-12)
            f_par = np.dot(F[i], tau) * tau
            if i == climb and it > iters // 4:
                force = F[i] - 2.0 * f_par
            else:
                spring = k_spring * (np.linalg.norm(real.carrier.diff(path[i], path[i + 1]))
                                     - np.linalg.norm(real.carrier.diff(path[i - 1], path[i]))) * tau
                force = F[i] - f_par + spring
            new[i] = path[i] + step * force
        path = real.carrier.wrap(new)
    E = real.V(path, **over)
    return {"barrier": float(E.max() - E[0]), "barrier_back": float(E.max() - E[-1]), "energies": E.tolist(),
            "saddle_index": int(np.argmax(E))}


def restoring_threshold(real: Realization, a: np.ndarray, b: np.ndarray, n: int = 64) -> float:
    """Largest force opposing a straight move from a to b: a bounded write must exceed it to remove the
    barrier on that path; a weaker write needs thermal activation."""
    t = np.linspace(0, 1, n)[:, None]
    d = real.carrier.diff(a, b)
    path = real.carrier.wrap(a[None] + t * d[None])
    u = d / max(np.linalg.norm(d), 1e-12)
    return float(max(1e-6, np.max(-(real.F(path) @ u))))


# ------------------------------------------------------------------------------------------------ hold / write
def nearest(real: Realization, q: np.ndarray, states: Sequence[np.ndarray]) -> np.ndarray:
    d = np.stack([real.carrier.distance(q, s[None]) for s in states], axis=-1)
    return np.argmin(d, axis=-1)


def hold_time(real: Realization, idx: int, states, D: float, rng, n_traj: int = 200, t_max: float = 400.0,
              over: Optional[Dict] = None) -> Dict[str, object]:
    q = np.repeat(states[idx][None], n_traj, axis=0)
    _, t, frac = integrate(real, q, t_max, D, rng, over, record=lambda x: np.mean(nearest(real, x, states) == idx),
                           n_record=200)
    below = np.flatnonzero(frac <= P_STAR)
    return {"time": float(t[below[0]]) if below.size else float(t_max), "censored": bool(below.size == 0),
            "t": t.tolist(), "fraction": frac.tolist()}


def write_test(real: Realization, i_old: int, i_new: int, states, h: float, T_w: float, T_rel: float, D: float,
               rng, n_traj: int = 200, mode: str = "fixed", sweep_from: Optional[float] = None,
               over: Optional[Dict] = None) -> Dict[str, object]:
    over = dict(over or {})
    q = np.repeat(states[i_old][None], n_traj, axis=0)
    target = states[i_new]
    rec = lambda x: np.mean(nearest(real, x, states) == i_new)
    path = None
    if mode == "sweep":
        v1 = over.get(real.control, real.params[real.control])
        path = lambda t: sweep_from + (v1 - sweep_from) * min(1.0, t / T_w)
    q, t1, f1 = integrate(real, q, T_w, D, rng, over, write=(target, h), control_path=path, record=rec)
    q, t2, f2 = integrate(real, q, T_rel, D, rng, over, record=rec)
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
    and backward between neighbouring scan values)."""
    param = param or real.control
    if values is None:
        values = np.linspace(*real.control_range, 7)
    values = np.asarray(values, float)
    width = float(abs(values[-1] - values[0]))
    events: List[Dict[str, object]] = []
    for a, b in zip(values[:-1], values[1:]):
        for v0, v1 in ((a, b), (b, a)):
            for q, k in fixed_points(real, rng, n_starts, {param: v0}):
                # a state on the boundary of an orthant (a species or a tube absent) sits where the drift
                # need not be smooth; its continuation is not attempted
                on_edge = real.carrier.kind == "orthant" and float(np.min(np.abs(q))) < 1e-6 * real.carrier.scale
                if n_unstable(k) == 0 and not on_edge:
                    _, ev = track(real, param, q, v0, v1)
                    if ev is not None:
                        events.append(ev)
    out: List[Dict[str, object]] = []
    for ev in sorted(events, key=lambda e: e["v"]):
        if not any(abs(ev["v"] - e["v"]) < 2e-3 * width
                   and real.carrier.distance(_canon(real, np.array(ev["q"])), _canon(real, np.array(e["q"])))
                   < 4 * _merge_tol(real) for e in out):
            out.append(ev)
    return out


def _slow_manifold_drift(real, q, vec, w, s, over):
    """Drift along vec on the slow manifold: for each s, the transverse coordinates are solved so that the
    drift has no component off vec (adiabatic elimination of the relaxing directions). Returns the drift and
    the largest remaining transverse drift relative to the largest drift along the grid."""
    n = q.size
    if n == 1:
        return np.array([float(w[0] * real.F(real.carrier.wrap(q + si * vec), **over)[0]) for si in s]), 0.0
    from scipy.linalg import null_space
    C = null_space(vec[None, :])
    f = np.empty(len(s))
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
            dz = np.linalg.lstsq(C.T @ jacobian(real, x, over) @ C, g, rcond=None)[0]
            nz = np.linalg.norm(dz)
            z = z - (dz * cap / nz if nz > cap else dz)
        z_side[side] = z
        x = real.carrier.wrap(q + s[j] * vec + C @ z)
        Fx = real.F(x, **over)
        f[j] = float(w @ Fx)
        trans = max(trans, float(np.linalg.norm(C.T @ Fx)))
        total = max(total, float(np.linalg.norm(Fx)))
    return f, trans / max(total, 1e-300)


def normal_form(real: Realization, q: np.ndarray, param: str, v: float, extra: Optional[Dict] = None,
                npts: int = 41) -> Dict[str, object]:
    """Reduce the drift at a write point to one coordinate s along the critical direction:
    ds/dt = a0 + a1 s + a2 s^2 + a3 s^3 (on the slow manifold), with the unfolding in the parameter
    (bias = the part of dF/dparam along the critical direction, kappa slope = d(a1)/dparam)."""
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
    i = int(cand[np.argmin(np.abs(ev.real[cand]))])
    lam = ev[i]
    scale = max(1e-12, float(np.abs(ev).max()))
    width = (abs(real.control_range[1] - real.control_range[0])
             if real.control_range is not None and param == real.control else max(1.0, abs(v)))
    dv = 1e-5 * width
    dF = (real.F(q, **{**over, param: v + dv}) - real.F(q, **{**over, param: v - dv})) / (2 * dv)
    out: Dict[str, object] = {"param": param, "value": float(v), "state": q.tolist(),
                              "kappa_critical": [float(-lam.real), float(-lam.imag)]}
    if abs(lam.imag) > 1e-6 * scale:
        out.update(kind="Hopf: the unwritten state starts to oscillate", hopf=True, omega=float(abs(lam.imag)),
                   asymmetry=0.0, bias=0.0, a3=0.0)
        return out
    vec = np.real(R[:, i])
    vec /= np.linalg.norm(vec)
    evl, L = np.linalg.eig(J.T)
    w = np.real(L[:, int(np.argmin(np.abs(evl - lam)))])
    w = w / np.dot(w, vec)
    half = 0.25 * base_length(real) if real.carrier.kind != "torus" else real.carrier.period / 8
    if real.carrier.kind == "orthant":
        nz = np.abs(vec) > 1e-9
        room = np.abs(q)[nz] / np.abs(vec)[nz]
        if room.size and np.min(room) > 1e-9:
            half = min(half, 0.5 * float(np.min(room)))
    for _ in range(8):
        s = np.linspace(-1.0, 1.0, npts) * half
        f, transverse = _slow_manifold_drift(real, q, vec, w, s, over)
        coef = np.polyfit(s, f, 3)
        misfit = float(np.max(np.abs(np.polyval(coef, s) - f))) / max(float(np.max(np.abs(f))), 1e-300)
        if np.all(np.isfinite(f)) and transverse < 1e-6 and misfit < 0.05:
            break
        half /= 2
    a3, a2, a1, a0 = (float(c) for c in coef)
    Jp = jacobian(real, q, {**over, param: v + dv})
    Jm = jacobian(real, q, {**over, param: v - dv})
    growth_partial = float(w @ ((Jp - Jm) / (2 * dv)) @ vec)
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
    bias = abs(bias_along) / (abs(bias_along) + abs(growth) * half + 1e-300)
    quad, cub = abs(a2) * half ** 2, abs(a3) * half ** 3
    asym = quad / (quad + cub) if quad + cub > 1e-12 * scale * half else 0.0
    if cub + quad <= 1e-6 * scale * half:
        kind = "degenerate: the drift vanishes along the critical direction (a line of states)"
    elif bias > 0.1 and asym > 0.1:
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
               transverse_residual=transverse, cubic_misfit=misfit, s=s.tolist(), f=f.tolist(), half=half)
    return out


def with_write_field(real: Realization, target: np.ndarray, h_range) -> Realization:
    """The realization with a bounded write toward target added to its drift; the field h is the control."""
    target = np.asarray(target, float)
    base = real

    def drift(q, p):
        return base.drift(q, p) + base.write_force(q, target, p["h"])

    r = replace(real, name=real.name + " + write field", drift=drift, potential=None,
                params={**real.params, "h": 0.0}, control="h", control_range=tuple(h_range))
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
