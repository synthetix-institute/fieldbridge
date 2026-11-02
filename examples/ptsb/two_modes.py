"""More than two products: two transverse modes that lose stability along one path.

    V = 8/3 x^3 - 4 x^2 + (1/2) y^2 [1 - 4 x] + (1/2) z^2 [1 - 4 (x - D)] + x (y^4 + z^4 + 2 c y^2 z^2)
        + eps_y y + eps_z z.

The mode y loses stability at x = 1/4 and the mode z at x = 1/4 + D. The quartic coupling c decides the products:

* c < 1: the minima lie off the axes, one in each quadrant (sign y, sign z): four products, labelled by two signs.
* c > 1: the minima lie on the axes, (+-a, 0) and (0, +-b): the mode that grows first suppresses the other. Four
  products, or three when a strong tilt eps_z removes one of the two z products.

The law. In the linear stage y and z are independent Gaussians along the same reference path x_ref(t), each obeying
the equation of ptsb.linear_law with its own stiffness and force. Then
  * for c < 1 the sign of each mode is P = Phi(mean / spread), and the quadrant shares are products of the two;
  * for c > 1 the family is decided by which mode reaches the saturation amplitude first, P(y family) =
    P(|y| > |z|), read when the deviation of the leading mode reaches `A_READ`, and the sign within a family by the
    Gaussian conditioned on that family. For equal means of zero, P(y family) = (2/pi) arctan(sigma_y / sigma_z).
"""
from __future__ import annotations

from math import sqrt
from typing import Dict, Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from ptsb import Phi

A_READ = 0.4          # deviation at which the leading mode is taken to have saturated (products at about 0.85)


class TwoModes:
    def __init__(self, c: float = 0.5, D: float = 0.0, eps_y: float = 0.0, eps_z: float = 0.0):
        self.c, self.D, self.eps_y, self.eps_z = float(c), float(D), float(eps_y), float(eps_z)
        self._geo = None

    def V(self, x, y, z):
        y2, z2 = y * y, z * z
        return 8 / 3 * x ** 3 - 4 * x * x + 0.5 * y2 * (1 - 4 * x) + 0.5 * z2 * (1 - 4 * (x - self.D)) \
            + x * (y2 * y2 + z2 * z2 + 2 * self.c * y2 * z2) + self.eps_y * y + self.eps_z * z

    def grad(self, x, y, z):
        y2, z2 = y * y, z * z
        return (8 * x * x - 8 * x - 2 * y2 - 2 * z2 + y2 * y2 + z2 * z2 + 2 * self.c * y2 * z2,
                y * (1 - 4 * x) + x * (4 * y2 * y + 4 * self.c * y * z2) + self.eps_y,
                z * (1 - 4 * (x - self.D)) + x * (4 * z2 * z + 4 * self.c * z * y2) + self.eps_z)

    def hessian(self, x, y, z):
        y2, z2, c = y * y, z * z, self.c
        Vxy = -4 * y + 4 * y2 * y + 4 * c * y * z2
        Vxz = -4 * z + 4 * z2 * z + 4 * c * z * y2
        return np.array([[16 * x - 8, Vxy, Vxz],
                         [Vxy, (1 - 4 * x) + x * (12 * y2 + 4 * c * z2), 8 * c * x * y * z],
                         [Vxz, 8 * c * x * y * z, (1 - 4 * (x - self.D)) + x * (12 * z2 + 4 * c * y2)]])

    def stationary(self, guess) -> Optional[np.ndarray]:
        sol = root(lambda q: np.array(self.grad(*q)), np.asarray(guess, float), jac=lambda q: self.hessian(*q),
                   tol=1e-12)
        return sol.x if sol.success and np.linalg.norm(self.grad(*sol.x)) < 1e-9 else None

    def ts1(self) -> np.ndarray:
        q = self.stationary([0.0, -self.eps_y, -self.eps_z / (1 + 4 * self.D)])
        return q

    def geometry(self) -> Dict[str, object]:
        """Minima and index-1 saddles beyond the VRIs, by multistart; the commit level lies 30 % of the way from the
        highest saddle between products down to the highest minimum."""
        if self._geo is not None:
            return self._geo
        found = []
        vals = [0.0, 0.1, 0.25, 0.4, 0.6, 0.75, 0.9]
        for x0 in (0.95, 1.1):
            for a in vals + [-v for v in vals[1:]]:
                for b in vals + [-v for v in vals[1:]]:
                    q = self.stationary([x0, a, b])
                    if q is None or q[0] < 0.6 or q[0] > 1.5:
                        continue
                    if any(np.linalg.norm(q - f) < 1e-5 for f in found):
                        continue
                    found.append(q)
        minima, saddles = [], []
        for q in found:
            ev = np.linalg.eigvalsh(self.hessian(*q))
            if (ev > 0).all():
                minima.append(q)
            elif ev[0] < 0 < ev[1]:
                saddles.append(q)
        vm = max(self.V(*q) for q in minima)
        vs = max(self.V(*q) for q in saddles)
        self._geo = {"minima": minima, "saddles": saddles, "level": vs - 0.3 * (vs - vm),
                     "xc": 0.5 * (0.25 + max(self.D, 0.0) + 1.0), "V_saddle": vs, "V_min": vm}
        return self._geo


def sample(s: TwoModes, n: int, rng, kT: float):
    """Thermal flux on the dividing plane through TS1: x momentum Rayleigh, y and z thermal in their own wells."""
    q = s.ts1()
    H = s.hessian(*q)
    y = q[1] + rng.normal(0, sqrt(kT / H[1, 1]), n)
    z = q[2] + rng.normal(0, sqrt(kT / H[2, 2]), n)
    return np.full(n, q[0]), y, z, rng.rayleigh(sqrt(kT), n), rng.normal(0, sqrt(kT), n), rng.normal(0, sqrt(kT), n)


def run(s: TwoModes, init, gamma: float, kT: float, dt: float = 0.005, t_max: float = 40.0, check_every: int = 4,
        backend: str = "numpy", seed: int = 0, device: str = "cuda") -> np.ndarray:
    """BAOAB. Returns (code, sy, sz): code 0 returned behind TS1 or escaped, 1 reached a product, 2 undecided; sy, sz
    the signs of y - y_c and z - z_c on arrival, and |y| and |z| through `first` (1 if |y| > |z|)."""
    g = s.geometry()
    level, xc = g["level"], g["xc"]
    x, y, z, px, py, pz = init
    c1 = float(np.exp(-gamma * dt))
    c2 = float(np.sqrt((1 - c1 * c1) * kT))
    if backend == "torch":
        import torch
        dev = torch.device(device)
        T = lambda a: torch.as_tensor(np.asarray(a, float), dtype=torch.float64, device=dev)
        x, y, z, px, py, pz = (T(v) for v in (x, y, z, px, py, pz))
        gen = torch.Generator(device=dev).manual_seed(seed)
        normal = lambda: torch.randn(x.shape, generator=gen, dtype=torch.float64, device=dev)
        full = lambda v: torch.full(x.shape, v, dtype=torch.int64, device=dev)
        where = torch.where
        to_np = lambda v: v.cpu().numpy()
    else:
        rng = np.random.default_rng(seed)
        normal = lambda: rng.standard_normal(x.shape)
        full = lambda v: np.full(x.shape, v, np.int64)
        where = np.where
        to_np = lambda v: v
    code, y_out, z_out = full(2), None, None
    y_out, z_out = y * 0, z * 0
    alive = full(1) > 0
    gx, gy, gz = s.grad(x, y, z)
    for k in range(int(round(t_max / dt))):
        px, py, pz = px - 0.5 * dt * gx, py - 0.5 * dt * gy, pz - 0.5 * dt * gz
        x, y, z = x + 0.5 * dt * px, y + 0.5 * dt * py, z + 0.5 * dt * pz
        px, py, pz = c1 * px, c1 * py, c1 * pz
        if c2 > 0:
            px, py, pz = px + c2 * normal() * alive, py + c2 * normal() * alive, pz + c2 * normal() * alive
        x, y, z = x + 0.5 * dt * px, y + 0.5 * dt * py, z + 0.5 * dt * pz
        gx, gy, gz = s.grad(x, y, z)
        px, py, pz = px - 0.5 * dt * gx, py - 0.5 * dt * gy, pz - 0.5 * dt * gz
        if (k + 1) % check_every == 0:
            inside = alive & (x > xc) & (s.V(x, y, z) < level)
            back = alive & ((x < -0.4) | (abs(x) > 5) | (abs(y) > 5) | (abs(z) > 5))
            code = where(back, full(0), code)
            code = where(inside, full(1), code)
            y_out = where(inside, y, y_out)
            z_out = where(inside, z, z_out)
            alive = alive & ~back & ~inside
            if not bool(alive.any()):
                break
            m = alive.double() if backend == "torch" else alive.astype(float)
            px, py, pz, gx, gy, gz = px * m, py * m, pz * m, gx * m, gy * m, gz * m
    return to_np(code), to_np(y_out), to_np(z_out)


def shares(s: TwoModes, code, y, z) -> Dict[str, float]:
    """Product shares among trajectories that reached a product, assigned to the nearest minimum of the surface."""
    M = np.array(s.geometry()["minima"])[:, 1:]
    done = code == 1
    Y = np.c_[y[done], z[done]]
    near = np.argmin(((Y[:, None, :] - M[None, :, :]) ** 2).sum(2), axis=1)
    out = {"n": int(code.size), "done": int(done.sum()), "reactant": int((code == 0).sum()),
           "undecided": int((code == 2).sum()), "minima": M.tolist(),
           "share": (np.bincount(near, minlength=len(M)) / max(done.sum(), 1)).tolist(),
           "P_y_plus": float((y[done] > 0).mean()) if done.any() else float("nan"),
           "P_z_plus": float((z[done] > 0).mean()) if done.any() else float("nan"),
           "P_both_plus": float(((y[done] > 0) & (z[done] > 0)).mean()) if done.any() else float("nan"),
           "P_y_family": float((np.abs(y[done]) > np.abs(z[done])).mean()) if done.any() else float("nan")}
    return out


# ------------------------------------------------------------------------------------------------------------ law
def law(s: TwoModes, speed: float, gamma: float, kT: float, t_end: float = 14.0, n_eval: int = 2801) -> Dict[str, float]:
    """Mean and spread of y and z along one reference path. Returns the sign laws (read at the largest amplification
    of each mode in the first excursion beyond its VRI) and the family law (read when the leading mode reaches
    A_READ, or at the first wall contact if it does not)."""
    q = s.ts1()
    H = s.hessian(*q)
    om_y, om_z = H[1, 1], H[2, 2]

    def rhs(t, Z):
        x, v, my, py_, mz, pz_, a11, a12, a21, a22, b11, b12, b21, b22, syy, syp, spp, szz, szp, sqq = Z
        ky = 1 - 4 * x
        kz = 1 - 4 * (x - s.D)
        Vx = 8 * x * x - 8 * x
        return [v, -Vx - gamma * v, py_, -ky * my - gamma * py_ - s.eps_y, pz_, -kz * mz - gamma * pz_ - s.eps_z,
                a21, a22, -ky * a11 - gamma * a21, -ky * a12 - gamma * a22,
                b21, b22, -kz * b11 - gamma * b21, -kz * b12 - gamma * b22,
                2 * syp, spp - ky * syy - gamma * syp, -2 * ky * syp - 2 * gamma * spp + 2 * gamma * kT,
                2 * szp, sqq - kz * szz - gamma * szp, -2 * kz * szp - 2 * gamma * sqq + 2 * gamma * kT]

    def behind(t, Z):
        return Z[0] + 0.4
    behind.terminal = True

    Z0 = [q[0], speed, q[1], 0, q[2], 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0]
    t = np.linspace(0, t_end, n_eval)
    sol = solve_ivp(rhs, (0, t_end), Z0, t_eval=t, rtol=1e-8, atol=1e-11, method="LSODA", events=(behind,))
    x, v, my, _, mz, _, a11, a12, _, _, b11, b12, _, _, syy, _, _, szz, _, _ = sol.y
    vy = a11 ** 2 * kT / om_y + a12 ** 2 * kT + syy
    vz = b11 ** 2 * kT / om_z + b12 ** 2 * kT + szz
    out = {}

    def sign_read(mean, var, vri):
        beyond = x > vri
        if not beyond.any():
            return float("nan")
        k0 = int(np.argmax(beyond))
        back = np.where(~beyond[k0:])[0]
        k1 = k0 + (int(back[0]) if back.size else len(x) - k0)
        kb = k0 + int(np.argmax(var[k0:k1]))
        return float(Phi(mean[kb] / sqrt(var[kb])))

    out["P_y_plus"] = sign_read(my, vy, 0.25)
    out["P_z_plus"] = sign_read(mz, vz, 0.25 + s.D)
    dev = np.maximum(np.sqrt(my ** 2 + vy), np.sqrt(mz ** 2 + vz))
    turn = np.where((v[:-1] > 0) & (v[1:] <= 0))[0]
    k_wall = int(turn[0]) + 1 if turn.size else len(x) - 1
    for A in (A_READ, 0.2, 0.6, 0.8):
        k = np.where(dev >= A)[0]
        kr = int(k[0]) if k.size else k_wall
        f = family(my[kr], sqrt(vy[kr]), mz[kr], sqrt(vz[kr]))
        if A == A_READ:
            out.update(f)
            out.update(x_read=float(x[kr]), s_y=float(sqrt(vy[kr])), s_z=float(sqrt(vz[kr])))
        out[f"P_y_family_{A:g}"] = f["P_y_family"]
    return out


def family(my, sy, mz, sz, grid: int = 2001) -> Dict[str, float]:
    """P(|y| > |z|) and the signs within each family for independent Gaussians."""
    span = 6 * max(sy, sz) + abs(my) + abs(mz)
    u = np.linspace(-span, span, grid)
    du = u[1] - u[0]
    fy = np.exp(-0.5 * ((u - my) / sy) ** 2) / (sy * sqrt(2 * np.pi))
    fz = np.exp(-0.5 * ((u - mz) / sz) ** 2) / (sz * sqrt(2 * np.pi))
    Fz_in = Phi((np.abs(u) - mz) / sz) - Phi((-np.abs(u) - mz) / sz)    # P(|z| < |u|)
    Fy_in = Phi((np.abs(u) - my) / sy) - Phi((-np.abs(u) - my) / sy)    # P(|y| < |u|)
    yplus = float((fy * Fz_in * (u > 0)).sum() * du)
    yminus = float((fy * Fz_in * (u < 0)).sum() * du)
    zplus = float((fz * Fy_in * (u > 0)).sum() * du)
    zminus = float((fz * Fy_in * (u < 0)).sum() * du)
    tot = yplus + yminus + zplus + zminus
    return {"P_y_family": (yplus + yminus) / tot, "y+": yplus / tot, "y-": yminus / tot, "z+": zplus / tot,
            "z-": zminus / tot}


def law_thermal(s: TwoModes, gamma: float, kT: float, n_speeds: int = 16, **kw) -> Dict[str, float]:
    u = (np.arange(n_speeds) + 0.5) / n_speeds
    v = np.sqrt(-2.0 * kT * np.log(1.0 - u))
    L = [law(s, float(vi), gamma, kT, **kw) for vi in v]
    return {k: float(np.nanmean([l[k] for l in L])) for k in L[0]}


# --------------------------------------------------------------------- the family law with the nonlinear separatrix
def _winner(s: TwoModes, x0, v0, y0, z0, ry, rz, gamma, t_end=30.0):
    """Noise-free motion from a point of the linear stage to a product: True if the y family wins."""
    g = s.geometry()
    level, xc = g["level"], g["xc"]

    def rhs(t, q):
        gx, gy, gz = s.grad(q[0], q[1], q[2])
        return [q[3], q[4], q[5], -gx - gamma * q[3], -gy - gamma * q[4], -gz - gamma * q[5]]

    def arrived(t, q):
        return s.V(q[0], q[1], q[2]) - level if q[0] > xc else 1.0
    arrived.terminal = True

    def escaped(t, q):
        return 5.0 - max(abs(q[0]), abs(q[1]), abs(q[2]))
    escaped.terminal = True
    sol = solve_ivp(rhs, (0, t_end), [x0, y0, z0, v0, ry * y0, rz * z0], method="LSODA", rtol=1e-6, atol=1e-8,
                    events=(arrived, escaped))
    return abs(sol.y[1, -1]) > abs(sol.y[2, -1])


def separatrix_law(s: TwoModes, speed: float, gamma: float, kT: float, A0: float = 0.2, t_end: float = 14.0,
                   grid: int = 401) -> Dict[str, float]:
    """P(y family) with the boundary between the families found on the surface. The linear stage gives the Gaussian
    (y, z) at the point where the deviation reaches A0 beyond the first VRI, with the velocity along the growing
    solution. From that radius a bisection on the angle theta = atan(|z| / |y|), separately for z > 0 and z < 0, finds
    the direction that separates the families in the noise-free motion. Nothing is fitted; noise after the read is
    neglected."""
    q = s.ts1()
    H = s.hessian(*q)
    om_y, om_z = H[1, 1], H[2, 2]

    def rhs(t, Z):
        x, v, my, py_, mz, pz_, syy, syp, spp, szz, szp, sqq = Z
        ky = 1 - 4 * x
        kz = 1 - 4 * (x - s.D)
        return [v, -(8 * x * x - 8 * x) - gamma * v, py_, -ky * my - gamma * py_ - s.eps_y,
                pz_, -kz * mz - gamma * pz_ - s.eps_z,
                2 * syp, spp - ky * syy - gamma * syp, -2 * ky * syp - 2 * gamma * spp + 2 * gamma * kT,
                2 * szp, sqq - kz * szz - gamma * szp, -2 * kz * szp - 2 * gamma * sqq + 2 * gamma * kT]

    def behind(t, Z):
        return Z[0] + 0.4
    behind.terminal = True

    Z0 = [q[0], speed, q[1], 0, q[2], 0, kT / om_y, 0, kT, kT / om_z, 0, kT]
    t = np.linspace(0, t_end, 2801)
    sol = solve_ivp(rhs, (0, t_end), Z0, t_eval=t, rtol=1e-8, atol=1e-11, method="LSODA", events=(behind,))
    x, v, my, _, mz, _, vy, cy, _, vz, cz, _ = sol.y
    dev = np.maximum(np.sqrt(my ** 2 + vy), np.sqrt(mz ** 2 + vz))
    ok = np.where((dev >= A0) & (x >= 0.25))[0]
    if not ok.size:
        return {"P_y_family": float("nan")}
    k = int(ok[0])
    ry, rz = cy[k] / vy[k], cz[k] / vz[k]

    def boundary(sign):
        lo, hi = 0.0, 0.5 * np.pi
        for _ in range(11):
            th = 0.5 * (lo + hi)
            if _winner(s, x[k], v[k], A0 * np.cos(th), sign * A0 * np.sin(th), ry, rz, gamma):
                lo = th
            else:
                hi = th
        return 0.5 * (lo + hi)

    th_up, th_dn = boundary(+1.0), boundary(-1.0)
    sy, sz = sqrt(vy[k]), sqrt(vz[k])
    Y = np.linspace(my[k] - 6 * sy, my[k] + 6 * sy, grid)
    W = np.linspace(mz[k] - 6 * sz, mz[k] + 6 * sz, grid)
    YY, WW = np.meshgrid(Y, W, indexing="ij")
    rho = np.exp(-0.5 * ((YY - my[k]) / sy) ** 2 - 0.5 * ((WW - mz[k]) / sz) ** 2)
    ang = np.arctan2(np.abs(WW), np.abs(YY))
    in_y = np.where(WW > 0, ang < th_up, ang < th_dn)
    return {"P_y_family": float((rho * in_y).sum() / rho.sum()), "theta_up": float(th_up), "theta_dn": float(th_dn),
            "x_read": float(x[k])}


def separatrix_thermal(s: TwoModes, gamma: float, kT: float, n_speeds: int = 12, **kw) -> Dict[str, float]:
    u = (np.arange(n_speeds) + 0.5) / n_speeds
    v = np.sqrt(-2.0 * kT * np.log(1.0 - u))
    L = [separatrix_law(s, float(vi), gamma, kT, **kw) for vi in v]
    return {k: float(np.nanmean([l.get(k, np.nan) for l in L])) for k in ("P_y_family", "theta_up", "theta_dn")}
