"""Two identical copies that exchange one variable: what holds them apart.

Two copies of a material, du/dt = F(u), exchange the variable u_j at the rate k:

    du_a/dt = F(u_a) - k E_j (u_a - u_b),        du_b/dt = F(u_b) - k E_j (u_b - u_a).

At an equal state u_a = u_b = p the shared state m = (u_a + u_b)/2 relaxes with the Jacobian J of one copy and the
half difference d = (u_a - u_b)/2 with J - 2 k E_j. The exchange can make the difference grow, at a rate lam along v,
while one copy alone is stable (the two-cell instability of Turing). To second order

    dm/dt = J (m - p) + H(d, d)/2,         dd/dt = (J - 2 k E_j) d + H(m - p, d) + T(d, d, d)/6,

so a difference of either sign moves the shared state the same way, m - p = -J^{-1} H(v, v) A^2 / 2 for d = A v, and
the moved shared state acts back on the difference. The amplitude obeys dA/dt = lam A + g A^3 with

    g = g_own + g_shared,   g_own = w . T(v, v, v) / 6,   g_shared = -w . H(J^{-1} H(v, v), v) / 2,

w the left vector (w . v = 1). With g < 0 the copies settle apart at |u_a - u_b| = 2 (lam/|g|)^(1/2). Holding the
shared state at p (external fluxes on both copies) removes g_shared: the held difference is 2 (lam/|g_own|)^(1/2), or
none near p when g_own > 0.

`balance` computes these numbers from the equations of one copy; `integrate` follows the pair, free or with the
shared state held, to its end. The expansion is local: the free prediction needs lam < |slowest rate of J| (the shared
state follows the difference), the held one lam < |next rate of J - 2 k E_j|. A distant unequal state, a run-away,
an oscillating difference and the lifetime of an unequal state under noise are outside it.
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np
from scipy.integrate import solve_ivp


def _f(real, mu):
    if real.control is None or mu is None:
        return lambda q: np.asarray(real.F(np.asarray(q, float)), float)
    return lambda q: np.asarray(real.F(np.asarray(q, float), **{real.control: float(mu)}), float)


def _jac(f, q, h):
    n = len(q); J = np.empty((n, n))
    for i in range(n):
        e = np.zeros(n); e[i] = h
        J[:, i] = (f(q + e) - f(q - e)) / (2 * h)
    return J


def equal_state(real, mu=None, starts: int = 12, seed: int = 0) -> Optional[np.ndarray]:
    """A stable state of one copy (relaxation from random starts, then Newton)."""
    f = _f(real, mu); n = real.carrier.dim; s = float(real.carrier.scale or 1.0)
    rng = np.random.default_rng(seed)
    q0 = rng.uniform(-2, 2, (starts, n)) * s
    if real.carrier.kind == "orthant":
        q0 = np.abs(q0)
    for q in np.vstack([np.zeros((1, n)), q0]):
        sol = solve_ivp(lambda t, y: f(y), (0, 200.0), q, method="LSODA", rtol=1e-8, atol=1e-11 * s)
        p = sol.y[:, -1]
        for _ in range(30):
            F = f(p)
            if not np.all(np.isfinite(F)):
                break
            p = p - np.linalg.lstsq(_jac(f, p, 1e-6 * s), F, rcond=None)[0]
        if np.all(np.isfinite(p)) and np.abs(f(p)).max() < 1e-9 * s and np.linalg.eigvals(_jac(f, p, 1e-6 * s)).real.max() < 0:
            return p
    return None


def balance(real, var: str, rate: float, mu=None, p=None) -> Dict[str, object]:
    """The growth of the difference and the two parts of its cubic coefficient at the equal state p."""
    names = list(real.spec["carrier"]["variables"]) if getattr(real, "spec", None) else None
    j = names.index(var) if names else int(var)
    f = _f(real, mu); s = float(real.carrier.scale or 1.0)
    p = equal_state(real, mu) if p is None else np.asarray(p, float)
    if p is None:
        return {"status": "no stable state of one copy found"}
    n = len(p); J = _jac(f, p, 1e-6 * s)
    L = J.copy(); L[j, j] -= 2 * rate
    ev, V = np.linalg.eig(L); o = np.argsort(-ev.real); lam = ev[o[0]]
    out: Dict[str, object] = {"exchanged": var, "rate": rate, "control": mu, "equal_state": p.tolist(),
                              "growth_rate": float(lam.real),
                              "slowest_rate_of_one_copy": float(np.linalg.eigvals(J).real.max()),
                              "next_rate_of_the_difference": float(ev[o[1]].real) if n > 1 else None}
    if lam.real <= 0:
        return {**out, "status": "the difference does not grow"}
    if abs(lam.imag) > 1e-9 * abs(lam.real):
        return {**out, "status": "the difference grows as an oscillation", "frequency": float(abs(lam.imag))}
    v = V[:, o[0]].real; v /= np.linalg.norm(v)
    el, W = np.linalg.eig(L.T); w = W[:, int(np.argmin(abs(el - lam)))].real; w /= w @ v
    h = 3e-3 * s

    def H(a, b):
        return (f(p + h * (a + b)) - f(p + h * (a - b)) - f(p - h * (a - b)) + f(p - h * (a + b))) / (4 * h * h)

    T3 = (f(p + 2 * h * v) - 2 * f(p + h * v) + 2 * f(p - h * v) - f(p - 2 * h * v)) / (2 * h ** 3)
    x2 = -np.linalg.solve(J, H(v, v)) / 2
    g_own = float(w @ T3) / 6; g_shared = float(w @ H(x2, v)); g = g_own + g_shared
    lam = float(lam.real)
    out.update({"status": "the difference grows", "direction": v.tolist(), "g": g, "g_own": g_own, "g_shared": g_shared,
                "shift_of_the_shared_state_per_A2": x2.tolist(),
                "free_valid": bool(lam < -out["slowest_rate_of_one_copy"]),
                "held_valid": bool(out["next_rate_of_the_difference"] is not None and lam < -out["next_rate_of_the_difference"])})
    if g < 0:
        A = (lam / -g) ** 0.5
        out.update({"free_difference": float(2 * A), "free_shift_of_the_shared_state": (x2 * A * A).tolist()})
    if g_own < 0:
        out["held_difference"] = float(2 * (lam / -g_own) ** 0.5)
    return out


def integrate(real, var: str, rate: float, mu=None, p=None, hold: bool = False, T: float = 3000.0) -> Dict[str, object]:
    """The pair from the equal state displaced along the growing direction, free or with the shared state held."""
    b = balance(real, var, rate, mu, p)
    if "direction" not in b:
        return {"status": b["status"]}
    f = _f(real, mu); p = np.asarray(b["equal_state"]); n = len(p); v = np.asarray(b["direction"])
    j = list(real.spec["carrier"]["variables"]).index(var); s = float(real.carrier.scale or 1.0)

    def free(t, y):
        ua, ub = y[:n], y[n:]
        fa, fb = f(ua), f(ub)
        fa[j] -= rate * (ua[j] - ub[j]); fb[j] -= rate * (ub[j] - ua[j])
        return np.concatenate([fa, fb])

    def held(t, d):
        r = (f(p + d) - f(p - d)) / 2
        r[j] -= 2 * rate * d[j]
        return r
    far = lambda t, y: np.abs(y).max() - 1e3 * s   # noqa: E731
    far.terminal = True
    d0 = 1e-3 * s * v
    if hold:
        sol = solve_ivp(held, (0, T), d0, method="LSODA", rtol=1e-9, atol=1e-12 * s, events=far)
        if sol.status == 1:
            return {"endpoint": "runs away"}
        d = sol.y[:, -1]
        rest = np.abs(held(0, d)).max() < 1e-7 * s
        return {"endpoint": "rest" if rest else "not at rest", "difference": (2 * d).tolist(), "size": float(2 * np.linalg.norm(d))}
    sol = solve_ivp(free, (0, T), np.concatenate([p + d0, p - d0]), method="LSODA", rtol=1e-9, atol=1e-12 * s, events=far)
    if sol.status == 1:
        return {"endpoint": "runs away"}
    y = sol.y[:, -1]; rest = np.abs(free(0, y)).max() < 1e-7 * s
    return {"endpoint": "rest" if rest else "not at rest", "difference": (y[:n] - y[n:]).tolist(),
            "size": float(np.linalg.norm(y[:n] - y[n:])), "shift_of_the_shared_state": ((y[:n] + y[n:]) / 2 - p).tolist()}


def report(real, var: str, rate: float, mu=None) -> Dict[str, object]:
    b = balance(real, var, rate, mu)
    out = {"balance": b}
    if "direction" in b:
        out["free"] = integrate(real, var, rate, mu, b["equal_state"])
        out["held"] = integrate(real, var, rate, mu, b["equal_state"], hold=True)
    return out


def markdown(name: str, r: Dict[str, object]) -> str:
    b = r["balance"]
    lines = [f"# Two copies of {name} exchanging {b.get('exchanged')} at rate {b.get('rate')}", "",
             f"- growth of the difference: {b.get('growth_rate'):.4g} ({b.get('status')})"]
    if "g" in b:
        lines += [f"- cubic coefficient g = {b['g']:.4g}: the copy itself {b['g_own']:.4g}, the shared state {b['g_shared']:.4g}",
                  f"- predicted |u_a - u_b|: free {b.get('free_difference', 'none near the equal state')}, "
                  f"shared state held {b.get('held_difference', 'none near the equal state')}",
                  f"- free pair: {r['free'].get('endpoint')} at {r['free'].get('size')}",
                  f"- shared state held: {r['held'].get('endpoint')} at {r['held'].get('size')}"]
    return "\n".join(lines) + "\n"
