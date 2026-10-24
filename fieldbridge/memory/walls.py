"""Memory in a chain of bistable elements: a written pattern is lost by the motion of domain walls.

The field module (fieldbridge.memory.fields) has one stable state per site: every mode relaxes at its own rate and
a write is lost mode by mode (Laws 1 and 2). Here each site has two stable states, and nearest neighbours are
coupled. A written pattern of signs is then a stable state of the whole chain, and it is lost by noise-activated
events (Law 3) whose rate depends on the pattern around the site that switches:

  m = 0   both neighbours aligned with the site: a switch creates a pair of walls (nucleation);
  m = 1   one neighbour aligned, one opposite: the site sits at a wall, and its switch moves the wall by one site
          (a hop);
  m = 2   both neighbours opposite: the site is a domain of length one, and its switch removes two walls
          (annihilation).

The chain is reduced to a kinetic chain of signs whose rates are those of the classes of the local pattern (the
five sites centred on the switching site, or only m), measured on the full stochastic chain by competing risks:
from a pattern, every site is exposed to its own rate until the first switch anywhere, so the number of first
switches in a class divided by the exposure of that class (the time to the first switch times the number of sites
of the class, summed over trajectories) estimates its rate. The kinetic chain then gives the loss of any written
pattern, and it is compared with the full chain from the same patterns.

Two carriers:

  KPOChain      Kerr parametric oscillators coupled by photon hopping, H = sum_j [(K/2) a_j+2 a_j2 - (p/2)(a_j+2 +
                a_j2) + Delta a_j+ a_j - E (a_j + a_j+)] - J sum_j (a_j+ a_j+1 + h.c.), loss kappa = 1 per mode, in the
                classical image of the Lindbladian: the mean field of beta = sqrt(K) a,
                    dbeta_j/dt = -i (|beta_j|^2 beta_j - p beta_j* + Delta beta_j - J (beta_j-1 + beta_j+1) - E)
                                 - beta_j / 2,
                with the vacuum noise D = K / 8 per quadrature. The deterministic dynamics does not depend on K, which
                sets the noise only (about sqrt(p) / K photons per mode). `ordered` adds the Weyl ordering shift of
                the Kerr term (Delta -> Delta - K), which makes the image the truncated Wigner equation. The field E
                (in units of sqrt(K)) is the longitudinal field of the cat qubit: it lowers the energy of x > 0. The
                flow is Hamiltonian plus a uniform contraction at rate 1/2, so the eigenvalues of every equilibrium
                come in pairs symmetric about -1/2 and there is no potential: rates are measured, not computed from a
                barrier.
  RFSquidChain  overdamped rf-SQUID loops near half a flux quantum, theta = phi - pi, U = sum_j [(theta_j -
                theta_x)^2 / 2 + beta_L cos theta_j] - J sum_j (theta_j - theta_x)(theta_j+1 - theta_x), dtheta =
                -dU/dtheta dt + sqrt(2 T) dW (the loop current is proportional to theta - theta_x). This chain is a gradient system; its rates follow from the barriers of U (Langer), which
                checks the measurement of rates. Inductive coupling of flux loops is antiferromagnetic (J < 0); on a
                bipartite chain without a field it maps onto J > 0 by reversing every other loop.

Limits. A site is read by the sign of its coordinate with hysteresis: a switch is counted when the coordinate
crosses to -threshold times the well position of the other well. The reduction holds when switches are rare
compared with the relaxation of the chain (activated regime); its certificate is the comparison with the full chain.
Patterns that are not stable states of the chain (a short domain at strong coupling) are reported as lost by the
deterministic dynamics, before any noise acts.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy import optimize


# ------------------------------------------------------------------------------------------------ carriers
def _bonds(N: int, ring: bool):
    """Ordered pairs (j, k) of neighbours, each bond in both directions."""
    pairs = [(j, j + 1) for j in range(N - 1)] + ([(N - 1, 0)] if ring and N > 2 else [])
    return pairs + [(k, j) for j, k in pairs]


@dataclass
class KPOChain:
    N: int = 12
    p: float = 4.0          # two-photon drive, units of kappa
    Delta: float = 0.0      # detuning of the mode from half the pump frequency
    J: float = 0.25         # photon hopping between neighbours
    E: float = 0.0          # longitudinal field, units of sqrt(K)
    K: float = 0.5          # Kerr constant: noise D = K / 8 per quadrature of beta
    ring: bool = True
    ordered: bool = False   # Weyl ordering shift Delta -> Delta - K (truncated Wigner)
    complex_state: bool = field(default=True, init=False, repr=False)

    @property
    def D(self) -> float:
        return self.K / 8.0

    @property
    def delta(self) -> float:
        return self.Delta - (self.K if self.ordered else 0.0)

    def well(self) -> float:
        """Amplitude of the stable state of an uncoupled element without field, |beta|^2 = sqrt(p^2 - 1/4) - Delta."""
        r2 = np.sqrt(max(self.p ** 2 - 0.25, 0.0)) - self.delta
        if r2 <= 0:
            raise ValueError("the element is not bistable: sqrt(p^2 - 1/4) <= Delta")
        return float(np.sqrt(r2))

    def neighbours(self, b: np.ndarray) -> np.ndarray:
        nb = np.roll(b, 1, -1) + np.roll(b, -1, -1)
        if not self.ring:
            nb[..., 0] = b[..., 1]
            nb[..., -1] = b[..., -2]
        return nb

    def drift(self, b: np.ndarray) -> np.ndarray:
        return (-1j * (np.abs(b) ** 2 * b - self.p * np.conj(b) + self.delta * b - self.J * self.neighbours(b) - self.E)
                - 0.5 * b)

    def to_real(self, b: np.ndarray) -> np.ndarray:
        return np.concatenate([b.real, b.imag])

    def from_real(self, x: np.ndarray) -> np.ndarray:
        return x[:self.N] + 1j * x[self.N:]

    def coordinate(self, b: np.ndarray) -> np.ndarray:
        """The coordinate whose sign is the bit: the quadrature x."""
        return b.real

    def start(self, signs: np.ndarray) -> np.ndarray:
        return self.well() * np.asarray(signs, float) + 0j

    def jacobian(self, b: np.ndarray) -> np.ndarray:
        """Real Jacobian in (x_1..x_N, y_1..y_N) from df = A dbeta + B dbeta*."""
        N = self.N
        A = -1j * (2 * np.abs(b) ** 2 + self.delta) - 0.5
        B = -1j * (b ** 2 - self.p)
        M = np.zeros((2 * N, 2 * N))
        i = np.arange(N)
        M[i, i] = (A + B).real
        M[N + i, i] = (A + B).imag
        M[i, N + i] = -(A - B).imag
        M[N + i, N + i] = (A - B).real
        # hopping: A_nb = i J, B_nb = 0 -> du/dy_k = -J, dv/dx_k = J
        for j, k in _bonds(N, self.ring):
            M[j, N + k] += -self.J
            M[N + j, k] += self.J
        return M

    def noise(self, rng, shape) -> np.ndarray:
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    def default_dt(self) -> float:
        return 0.005 * 4.0 / max(self.p + 2 * abs(self.J), 4.0)


@dataclass
class RFSquidChain:
    N: int = 12
    beta_L: float = 2.0     # screening parameter (double well for beta_L > 1)
    theta_x: float = 0.0    # flux bias from half a flux quantum (the field)
    J: float = 0.1          # coupling, J > 0 aligns neighbours
    T: float = 0.1          # noise temperature (units of the inductive energy)
    ring: bool = True
    complex_state: bool = field(default=False, init=False, repr=False)

    @property
    def D(self) -> float:
        return self.T

    def well(self) -> float:
        """theta_0 > 0 with theta_0 = beta_L sin theta_0 (uncoupled, unbiased)."""
        if self.beta_L <= 1:
            raise ValueError("the loop is not bistable: beta_L <= 1")
        return float(optimize.brentq(lambda t: t - self.beta_L * np.sin(t), 1e-6, np.pi))

    def neighbours(self, th: np.ndarray) -> np.ndarray:
        nb = np.roll(th, 1, -1) + np.roll(th, -1, -1)
        if not self.ring:
            nb[..., 0] = th[..., 1]
            nb[..., -1] = th[..., -2]
        return nb

    def energy(self, th: np.ndarray) -> np.ndarray:
        u = th - self.theta_x
        bond = u * np.roll(u, -1, -1)
        if not self.ring:
            bond = bond[..., :-1]
        return np.sum(u ** 2 / 2 + self.beta_L * np.cos(th), -1) - self.J * np.sum(bond, -1)

    def drift(self, th: np.ndarray) -> np.ndarray:
        u = th - self.theta_x
        return -u + self.beta_L * np.sin(th) + self.J * self.neighbours(u)

    def to_real(self, th: np.ndarray) -> np.ndarray:
        return np.asarray(th, float)

    def from_real(self, x: np.ndarray) -> np.ndarray:
        return np.asarray(x, float)

    def coordinate(self, th: np.ndarray) -> np.ndarray:
        return th

    def start(self, signs: np.ndarray) -> np.ndarray:
        return self.well() * np.asarray(signs, float)

    def jacobian(self, th: np.ndarray) -> np.ndarray:
        N = self.N
        M = np.diag(-1.0 + self.beta_L * np.cos(th))
        for j, k in _bonds(N, self.ring):
            M[j, k] += self.J
        return M

    def noise(self, rng, shape) -> np.ndarray:
        return rng.standard_normal(shape)

    def default_dt(self) -> float:
        return 0.01


# ------------------------------------------------------------------------------------------------ deterministic
def fixed_point(chain, signs: Sequence[int], guess: Optional[np.ndarray] = None) -> Optional[Dict[str, object]]:
    """The stable equilibrium of the chain that carries the pattern `signs` (Newton from the uncoupled wells), or
    None if Newton does not converge to an equilibrium with these signs or the equilibrium is unstable."""
    s = np.asarray(signs, int)
    x0 = chain.to_real(chain.start(s) if guess is None else guess)
    f = lambda x: chain.to_real(chain.drift(chain.from_real(x)))
    sol = optimize.root(f, x0, jac=lambda x: chain.jacobian(chain.from_real(x)), method="hybr", tol=1e-13)
    z = chain.from_real(sol.x)
    if not sol.success or np.abs(f(sol.x)).max() > 1e-9:
        return None
    if not np.array_equal(np.sign(chain.coordinate(z)).astype(int), s):
        return None
    ev = np.linalg.eigvals(chain.jacobian(z))
    if ev.real.max() >= 0:
        return None
    i = int(np.argmax(ev.real))
    return {"state": z, "slowest": complex(ev[i]), "eigenvalues": ev}


def loss_threshold(chain, signs: Sequence[int], name: str, lo: float, hi: float, tol: float = 1e-4) -> Dict:
    """The value of parameter `name` between lo (pattern stable) and hi at which the pattern stops being a stable
    equilibrium, by continuation from lo and bisection of the last step."""
    s = np.asarray(signs, int)
    fp = fixed_point(replace(chain, **{name: lo}), s)
    if fp is None:
        raise ValueError(f"the pattern is not stable at {name} = {lo}")
    grid = np.linspace(lo, hi, 61)
    guess, last, keep = fp["state"], lo, fp
    for v in grid[1:]:
        fp = fixed_point(replace(chain, **{name: float(v)}), s, guess)
        if fp is None:
            a, b = last, float(v)
            while b - a > tol:
                c = 0.5 * (a + b)
                f2 = fixed_point(replace(chain, **{name: c}), s, guess)
                if f2 is None:
                    b = c
                else:
                    a, guess, keep = c, f2["state"], f2
            # the slowest eigenvalue of the last stable state says how it is lost: real (a fold: the state meets a
            # saddle) or complex (an oscillatory instability)
            return {"parameter": name, "value": 0.5 * (a + b), "error": 0.5 * (b - a), "reached": True,
                    "last_stable": a, "slowest_before": [keep["slowest"].real, keep["slowest"].imag],
                    "state_before": keep["state"]}
        guess, last, keep = fp["state"], float(v), fp
    return {"parameter": name, "value": None, "beyond": hi, "reached": False}


def wall_velocity(chain, E: float, t_max: float = 400.0, dt: Optional[float] = None,
                  threshold: float = 0.6) -> Dict[str, object]:
    """Deterministic motion of a single wall in the field E (open chain, + domain on the first quarter, started at
    the stable state of the pattern at zero field): the time at which each site switches (sign with hysteresis, as in
    first_switches), and the velocity, sites per unit time, from the switches of the sites between N/4 + 3 and
    7N/8. `ordered` is False when sites switch out of order (the minority domain is lost as a whole rather than by
    the motion of the wall); `pinned` when fewer than four sites switch within t_max."""
    field_name = "E" if hasattr(chain, "E") else "theta_x"
    ch = replace(chain, ring=False, **{field_name: 0.0})
    s = np.where(np.arange(ch.N) < ch.N // 4, 1, -1)
    fp = fixed_point(ch, s)
    if fp is None:
        raise ValueError("the single-wall pattern is not a stable state at zero field")
    ch = replace(ch, **{field_name: E})
    dt = dt or ch.default_dt()
    z = fp["state"].copy()
    bits = s.copy()
    lim = threshold * ch.well()
    when = np.full(ch.N, np.nan)
    t = 0.0
    lo, hi = ch.N // 4 + 3, (7 * ch.N) // 8
    while t < t_max and np.isnan(when[hi]):
        a0 = ch.drift(z)
        z = z + 0.5 * dt * (a0 + ch.drift(z + dt * a0))
        t += dt
        flip = ch.coordinate(z) * bits < -lim
        if flip.any():
            bits = np.where(flip, -bits, bits)
            when[flip & np.isnan(when)] = t
    switched = np.flatnonzero(~np.isnan(when))
    if switched.size < 4:
        return {"E": E, "velocity": 0.0, "pinned": True, "switched": int(switched.size), "times": when.tolist()}
    sel = np.arange(lo, hi + 1)
    sel = sel[~np.isnan(when[sel])]
    ordered = bool(np.all(np.diff(when[ch.N // 4:][~np.isnan(when[ch.N // 4:])]) >= -1e-9))
    v = float(np.polyfit(when[sel], sel, 1)[0]) if sel.size >= 3 else float("nan")
    return {"E": E, "velocity": v, "pinned": False, "ordered": ordered, "switched": int(switched.size),
            "times": when.tolist()}


def langer_rate(chain, signs: Sequence[int], site: int) -> Dict[str, object]:
    """Rate at which `site` of the pattern switches in a gradient chain (RFSquidChain), from the barrier of U:
    Gamma = |lambda_u| / (2 pi) sqrt(det H_min / |det H_saddle|) exp(-dU / T) (Langer, overdamped). The saddle is
    found by Newton from the stable state with the switching site moved to the top of its barrier, and must have one
    unstable direction."""
    if not hasattr(chain, "energy"):
        raise TypeError("the Langer rate needs a gradient chain")
    s = np.asarray(signs, int)
    fp = fixed_point(chain, s)
    if fp is None:
        raise ValueError("the pattern is not a stable state of the chain")
    guess = np.array(fp["state"], float)
    guess[site] = 0.0
    sol = optimize.root(chain.drift, guess, jac=chain.jacobian, method="hybr", tol=1e-13)
    Hs = -chain.jacobian(sol.x)
    ev_s = np.linalg.eigvalsh(Hs)
    if not sol.success or (ev_s < 0).sum() != 1:
        raise ValueError("no saddle with one unstable direction near the switch of this site")
    Ha = -chain.jacobian(fp["state"])
    ev_a = np.linalg.eigvalsh(Ha)
    dU = float(chain.energy(sol.x) - chain.energy(fp["state"]))
    pref = abs(ev_s[0]) / (2 * np.pi) * np.sqrt(np.prod(ev_a) / abs(np.prod(ev_s)))
    return {"barrier": dU, "prefactor": float(pref), "rate": float(pref * np.exp(-dU / chain.T)),
            "saddle": sol.x}


# ------------------------------------------------------------------------------------------------ classes
_DIGIT = {-1: 0, 0: 1, 1: 2}
_CHAR = {-1: "-", 0: ".", 1: "+"}


def _code(window: Sequence[int]) -> int:
    return int(sum(_DIGIT[int(s)] * 3 ** k for k, s in enumerate(window)))


def class_name(window: Sequence[int], symmetric: bool) -> str:
    """Canonical name of a five-site window centred on the switching site: under reflection, and (without a field)
    under reversal of every sign, so that the centre reads '+'. '.' is a site beyond the end of an open chain."""
    w = np.asarray(window, int)
    if symmetric and w[2] < 0:
        w = -w
    a = "".join(_CHAR[int(v)] for v in w)
    return min(a, a[::-1])


def nn_name(window: Sequence[int], symmetric: bool) -> str:
    """Nearest-neighbour class: m = number of opposite neighbours, the centre sign with a field, 'e' at an end."""
    w = np.asarray(window, int)
    c = w[2]
    m = int(w[1] == -c) + int(w[3] == -c)
    end = "e" if (w[1] == 0 or w[3] == 0) else ""
    return f"m{m}{end}" + ("" if symmetric else ("+" if c > 0 else "-"))


def windows(signs: Sequence[int], ring: bool) -> np.ndarray:
    s = np.asarray(signs, int)
    N = s.size
    if ring:
        return np.stack([np.roll(s, -k) for k in (-2, -1, 0, 1, 2)], 1)
    pad = np.concatenate([[0, 0], s, [0, 0]])
    return np.stack([pad[j:j + 5] for j in range(N)])


_ALL_WINDOWS = [tuple(int(v) for v in w) for w in np.array(np.meshgrid(*[[-1, 0, 1]] * 5, indexing="ij")).reshape(5, -1).T]


def is_symmetric(chain) -> bool:
    return (getattr(chain, "E", 0.0) == 0.0) and (getattr(chain, "theta_x", 0.0) == 0.0)


# ------------------------------------------------------------------------------------------------ stochastic chain
def _step(chain, z: np.ndarray, dt: float, rng) -> np.ndarray:
    """Stochastic Heun step with additive noise of D per coordinate."""
    xi = np.sqrt(2 * chain.D * dt) * chain.noise(rng, z.shape)
    a0 = chain.drift(z)
    g = z + dt * a0 + xi
    return z + 0.5 * dt * (a0 + chain.drift(g)) + xi


def stationary_sample(chain, fp: Dict[str, object], n: int, rng) -> np.ndarray:
    """n states drawn from the linearized stationary fluctuations about a stable equilibrium (covariance S from
    A S + S A^T + 2 D I = 0), so that a run starts in the quasi-stationary state of the basin instead of at its
    bottom, which would add the relaxation time to every first-passage time."""
    from scipy.linalg import solve_continuous_lyapunov
    A = chain.jacobian(fp["state"])
    S = solve_continuous_lyapunov(A, -2 * chain.D * np.eye(A.shape[0]))
    L = np.linalg.cholesky(0.5 * (S + S.T))
    x = chain.to_real(fp["state"])[None] + rng.standard_normal((n, A.shape[0])) @ L.T
    return np.stack([chain.from_real(v) for v in x])


def first_switches(chain, signs: Sequence[int], n_traj: int, t_max: float, rng, dt: Optional[float] = None,
                   threshold: float = 0.6, start: Optional[np.ndarray] = None) -> Dict[str, np.ndarray]:
    """Time to the first switch anywhere in the chain and the site that switched (-1 if none by t_max), for n_traj
    chains started in the linearized stationary fluctuations about the stable state of the pattern."""
    s = np.asarray(signs, int)
    if start is None:
        fp = fixed_point(chain, s)
        if fp is None:
            raise ValueError("the pattern is not a stable state of the chain")
        z = stationary_sample(chain, fp, n_traj, rng)
    else:
        z = np.repeat(np.asarray(start)[None], n_traj, 0)
    dt = dt or chain.default_dt()
    lim = threshold * chain.well()
    T = np.full(n_traj, float(t_max))
    site = np.full(n_traj, -1)
    alive = np.arange(n_traj)
    t = 0.0
    while t < t_max - 1e-12 and alive.size:
        z[alive] = _step(chain, z[alive], dt, rng)
        t += dt
        crossed = chain.coordinate(z[alive]) * s[None] < -lim
        hit = crossed.any(1)
        if hit.any():
            k = alive[hit]
            T[k] = t
            site[k] = np.argmax(crossed[hit], 1)
            alive = alive[~hit]
    return {"time": T, "site": site}


def flip_rates(chain, patterns: Sequence[Sequence[int]], n_traj: int, t_max: float, rng, dt: Optional[float] = None,
               threshold: float = 0.6) -> Dict[str, object]:
    """Rates of the classes of local pattern by competing risks over the given patterns (see the module text).
    Returns, for the five-site classes and the nearest-neighbour classes, events, exposure, rate and its Poisson
    error. Patterns that are not stable states are listed under `lost`."""
    sym = is_symmetric(chain)
    ev5, ex5, evn, exn = {}, {}, {}, {}
    lost, used = [], 0
    for s in patterns:
        s = np.asarray(s, int)
        try:
            r = first_switches(chain, s, n_traj, t_max, rng, dt, threshold)
        except ValueError:
            lost.append("".join("+" if v > 0 else "-" for v in s))
            continue
        used += 1
        W = windows(s, chain.ring)
        names5 = [class_name(w, sym) for w in W]
        namesn = [nn_name(w, sym) for w in W]
        exposure = float(r["time"].sum())
        for n5, nn in zip(names5, namesn):
            ex5[n5] = ex5.get(n5, 0.0) + exposure
            exn[nn] = exn.get(nn, 0.0) + exposure
        for k in r["site"][r["site"] >= 0]:
            ev5[names5[k]] = ev5.get(names5[k], 0) + 1
            evn[namesn[k]] = evn.get(namesn[k], 0) + 1

    def table(ev, ex):
        return {c: {"events": int(ev.get(c, 0)), "exposure": ex[c], "rate": ev.get(c, 0) / ex[c],
                    "error": np.sqrt(max(ev.get(c, 0), 1)) / ex[c]} for c in sorted(ex)}

    return {"five_site": table(ev5, ex5), "nearest": table(evn, exn), "patterns_used": used, "lost": lost,
            "symmetric": sym, "n_traj": n_traj, "t_max": t_max, "threshold": threshold}


def evolve_pattern(chain, signs: Sequence[int], times: Sequence[float], n_traj: int, rng, dt: Optional[float] = None,
                   threshold: float = 0.6) -> Dict[str, object]:
    """Full stochastic chain from the stable state of a pattern: the bits (sign with hysteresis) at the given times,
    summarized as the overlap with the written pattern, the fraction of chains that read the whole pattern, the
    fraction with no switch yet, and the number of walls."""
    s = np.asarray(signs, int)
    fp = fixed_point(chain, s)
    if fp is None:
        raise ValueError("the pattern is not a stable state of the chain")
    dt = dt or chain.default_dt()
    z = stationary_sample(chain, fp, n_traj, rng)
    bits = np.repeat(s[None], n_traj, 0)
    clean = np.ones(n_traj, bool)
    lim = threshold * chain.well()
    out = {"times": [], "overlap": [], "word": [], "no_switch": [], "walls": []}
    t = 0.0
    for target in np.asarray(times, float):
        while t < target - 1e-12:
            z = _step(chain, z, dt, rng)
            t += dt
            x = chain.coordinate(z)
            flip = x * bits < -lim
            if flip.any():
                bits = np.where(flip, -bits, bits)
                clean &= ~flip.any(1)
        _record(out, target, bits, s, clean, chain.ring)
    return out


def _record(out, t, bits, s, clean, ring):
    out["times"].append(float(t))
    out["overlap"].append(float(np.mean(bits * s[None])))
    out["word"].append(float(np.mean(np.all(bits == s[None], 1))))
    out["no_switch"].append(float(np.mean(clean)))
    d = bits != np.roll(bits, -1, 1)
    if not ring:
        d = d[:, :-1]
    out["walls"].append(float(np.mean(d.sum(1))))


# ------------------------------------------------------------------------------------------------ kinetic chain
class KineticChain:
    """Chain of signs switching at the rates of their local class: five-site classes where measured, otherwise the
    nearest-neighbour class."""

    def __init__(self, N: int, ring: bool, rates: Dict[str, object], symmetric: bool, use_five: bool = True,
                 min_events: int = 20):
        self.N, self.ring, self.sym = N, ring, symmetric
        five = {c: v["rate"] for c, v in rates["five_site"].items() if v["events"] >= min_events} if use_five else {}
        nn = {c: v["rate"] for c, v in rates["nearest"].items()}
        self.table = np.full(3 ** 5, np.nan)
        self.missing = set()
        for w in _ALL_WINDOWS:
            if w[2] == 0:
                continue
            c5, cn = class_name(w, symmetric), nn_name(w, symmetric)
            if c5 in five:
                self.table[_code(w)] = five[c5]
            elif cn in nn and rates["nearest"][cn]["events"] > 0:
                self.table[_code(w)] = nn[cn]
            else:
                self.missing.add(cn)
        self.weights = 3 ** np.arange(5)

    def site_rates(self, s: np.ndarray) -> np.ndarray:
        """s: (n, N) signs."""
        if self.ring:
            W = np.stack([np.roll(s, -k, 1) for k in (-2, -1, 0, 1, 2)], -1)
        else:
            pad = np.pad(s, ((0, 0), (2, 2)))
            W = np.stack([pad[:, k:k + self.N] for k in range(5)], -1)
        code = ((W + 1) * self.weights).sum(-1)
        r = self.table[code]
        if np.isnan(r).any():
            raise ValueError(f"no rate for some local classes (nearest-neighbour classes {sorted(self.missing)})")
        return r

    def evolve(self, signs: Sequence[int], times: Sequence[float], n: int, rng) -> Dict[str, object]:
        s0 = np.asarray(signs, int)
        times = np.asarray(times, float)
        s = np.repeat(s0[None], n, 0)
        t = np.zeros(n)
        clean = np.ones(n, bool)
        out = {"times": [], "overlap": [], "word": [], "no_switch": [], "walls": []}
        snap = np.zeros((times.size, n, self.N), int)
        snap_clean = np.zeros((times.size, n), bool)
        nxt = np.zeros(n, int)  # index of the next recording time per chain
        active = np.arange(n)
        while active.size:
            r = self.site_rates(s[active])
            R = r.sum(1)
            tau = rng.exponential(1.0 / np.maximum(R, 1e-300))
            t_new = t[active] + tau
            # record every time passed before this switch
            for i in range(times.size):
                rec = active[(nxt[active] == i) & (t_new >= times[i])]
                if rec.size:
                    snap[i, rec] = s[rec]
                    snap_clean[i, rec] = clean[rec]
                    nxt[rec] = i + 1
            done = nxt[active] >= times.size
            go = ~done
            idx = active[go]
            if idx.size:
                u = rng.random(idx.size) * R[go]
                site = np.argmax(np.cumsum(r[go], 1) >= u[:, None], 1)
                s[idx, site] *= -1
                clean[idx] = False
                t[idx] = t_new[go]
            active = idx
        for i, tt in enumerate(times):
            _record(out, tt, snap[i], s0, snap_clean[i], self.ring)
        return out


# ------------------------------------------------------------------------------------------------ retention laws
def word_rate(kc: KineticChain, signs: Sequence[int]) -> float:
    """Rate of the first error of a written word: the sum of the rates of its sites (exact for the kinetic chain)."""
    return float(kc.site_rates(np.asarray(signs, int)[None])[0].sum())


def glauber_overlap(t: np.ndarray, hop: float, nucleation: float = 0.0) -> np.ndarray:
    """Overlap of a random pattern in Glauber's chain whose site at a wall switches at `hop` and an aligned site at
    `nucleation`: q(t) = exp(-2 hop t) I_0(2 (hop - nucleation) t) (Glauber 1963)."""
    from scipy.special import ive
    t = np.asarray(t, float)
    a, g = 2 * hop, 2 * (hop - nucleation)
    return np.exp(-(a - g) * t) * ive(0, g * t)


def ring_reversal_rate(N: int, nucleation: float, hop: float, annihilation: float) -> float:
    """Rate at which an aligned word on a ring of N sites reverses through one nucleus (nearest-neighbour kinetic
    chain, one nucleus at a time): N Gamma_0 h, with h = p / (N (1 - p) + 4 p - 2) the probability that a one-site
    nucleus grows to the whole ring, p = 2 Gamma_1 / (Gamma_2 + 2 Gamma_1) the probability that it escapes
    annihilation at each visit of length one (and symmetrically at length N - 1). The overlap of the word decays at
    twice this rate; for large N it tends to 4 Gamma_0 Gamma_1 / Gamma_2."""
    p = 2 * hop / (annihilation + 2 * hop)
    return float(N * nucleation * p / (N * (1 - p) + 4 * p - 2))


def mean_pattern_field(hop: float, nucleation: float, L: int):
    """The linear field of fieldbridge.memory.fields whose relaxation is that of the mean pattern of Glauber's chain:
    d<s_j>/dt = -2 Gamma_0 <s_j> + (Gamma_1 - Gamma_0) lap <s_j>, i.e. a non-conserved field with M = 1,
    kappa0 = 2 Gamma_0 (the mass, from nucleation) and Dgrad = Gamma_1 - Gamma_0 (the mobility, from the hops of
    walls). Exact for Glauber's rates (Gamma_2 = 2 Gamma_1 - Gamma_0); for other rates the kinetic chain is the
    reference and this field is its Glauber approximation."""
    from .fields import FieldModel
    return FieldModel(d=1, L=L, conserved=False, M=1.0, kappa0=2 * nucleation, Dgrad=hop - nucleation, T=1.0)


def mean_pattern_overlap(times: Sequence[float], hop: float, nucleation: float, L: int = 4096) -> np.ndarray:
    """Overlap of a random written pattern with its mean evolution under mean_pattern_field: (1/L) sum_k
    exp(-kappa(k) t), kappa(k) = 2 Gamma_0 + (Gamma_1 - Gamma_0) 2 (1 - cos k)."""
    from .fields import _lambda
    m = mean_pattern_field(hop, nucleation, L)
    lam, _ = _lambda(1, L)
    rate = m.M * (m.kappa0 + m.Dgrad * lam)
    return np.array([float(np.mean(np.exp(-rate * t))) for t in np.asarray(times, float)])


def from_spec(spec: Dict):
    """A chain from a specification of kind 'bistable-chain': {"carrier": "kpo" | "rf-squid", "chain": {...}}."""
    from .spec import SCHEMA, SpecError
    if spec.get("schema") != SCHEMA or spec.get("kind") != "bistable-chain":
        raise SpecError(f"a bistable-chain specification needs schema {SCHEMA!r} and kind 'bistable-chain'")
    for key in ("question", "assumptions", "provenance"):
        if key not in spec:
            raise SpecError(f"Missing required field {key!r}")
    cls = {"kpo": KPOChain, "rf-squid": RFSquidChain}.get(spec.get("carrier"))
    if cls is None:
        raise SpecError("carrier must be 'kpo' or 'rf-squid'")
    entries = spec.get("chain", {})
    allowed = {k for k, f in cls.__dataclass_fields__.items() if f.init}
    unknown = set(entries) - allowed
    if unknown:
        raise SpecError(f"unknown chain entries {sorted(unknown)}; allowed: {sorted(allowed)}")
    ch = cls(**entries)
    if ch.N < 3:
        raise SpecError("chain.N must be at least 3")
    ch.well()
    return ch
