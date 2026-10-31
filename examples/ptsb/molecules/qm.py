"""DFT along the reaction path of a post-transition-state bifurcation: Hessian at TS1, IRC, transverse direction, curvature scans.

    python qm.py hess 31 && python qm.py irc 31 && python qm.py vrt 31 && python qm.py scan 31    (one GPU, 4 threads)
    python qm.py traj 31                                                                          (noise-free trajectory from TS1)
    python qm.py hessframe 31 --k 10,20,30                                                        (Hessians at trajectory frames, times in fs)
    python qm.py trajy 31 --k=-2,-1,1,2                                                           (trajectories displaced along the transverse direction)
    python qm.py hessat 31 --k 3,9,13                                                             (Hessians at IRC points)

Level of theory as in Yamaki et al. for entries 31 and 32: UB3LYP-D3(BJ)/def2-SVP, SMD (1,4-dioxane), a (99,590) grid, with density
fitting, which is the one difference from their calculation. Internally: coordinates in bohr, mass-weighted q = sqrt(m) x in
sqrt(amu) bohr, energies in hartree. Everything is written to results/ as it is computed, so that a job can be restarted.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from yamaki_geometries import ELEMENTS, coords  # noqa: E402

BOHR = 0.529177210903
MASS = {"H": 1.00782503, "C": 12.0, "F": 18.99840316, "N": 14.0030740, "O": 15.9949146}
# forming bond (radical carbon, central carbon of the allene) that fixes the product side of the imaginary mode
BOND = {31: (10, 2), 32: (6, 2)}
OUT = HERE / "results"


class Calc:
    """Energy and gradient with the density of the previous point as the starting guess."""

    def __init__(self, entry):
        self.entry = entry
        self.elements = ELEMENTS[entry]
        self.dm = None
        self.calls = 0

    def __call__(self, x_bohr):
        from gpu4pyscf.dft import uks
        from pyscf import gto
        atom = [(e, tuple(r)) for e, r in zip(self.elements, np.asarray(x_bohr).reshape(-1, 3))]
        mol = gto.M(atom=atom, unit="Bohr", basis="def2-svp", spin=1, charge=0, verbose=0)
        mf = uks.UKS(mol, xc="b3lyp").density_fit()
        mf.disp = "d3bj"
        mf.grids.atom_grid = (99, 590)
        mf.conv_tol = 1e-10
        ms = mf.SMD()
        ms.with_solvent.solvent = "1,4-dioxane"
        ms.with_solvent.lebedev_order = 29
        e = ms.kernel(dm0=self.dm)
        self.dm = ms.make_rdm1()
        if not ms.converged:
            ms = ms.newton() if hasattr(ms, "newton") else ms
            e = ms.kernel(dm0=self.dm)
            self.dm = ms.make_rdm1()
        g = ms.nuc_grad_method().kernel()
        self.calls += 1
        return float(e), np.asarray(g).ravel()


def masses(entry):
    return np.repeat([MASS[e] for e in ELEMENTS[entry]], 3)


def tr_basis(q, m):
    """Orthonormal mass-weighted translations and rotations at q."""
    sm = np.sqrt(m)
    x = (q / sm).reshape(-1, 3)
    xc = x - (m[::3, None] * x).sum(0) / m[::3].sum()
    B = []
    for a in range(3):
        v = np.zeros((len(x), 3))
        v[:, a] = 1.0
        B.append((v.ravel() * sm))
    for a in np.eye(3):
        B.append((np.cross(a, xc).ravel() * sm))
    Q, _ = np.linalg.qr(np.array(B).T)
    return Q


def project(q, m):
    Q = tr_basis(q, m)
    return np.eye(q.size) - Q @ Q.T


def numeric_hessian(calc, x, h=0.01, ckpt=None):
    """Central differences of the gradient in Cartesian coordinates, written row by row to ckpt."""
    n = x.size
    H = np.full((n, n), np.nan)
    if ckpt is not None and Path(ckpt).exists():
        H = np.load(ckpt)
    for i in range(n):
        if np.isfinite(H[i]).all():
            continue
        d = np.zeros(n)
        d[i] = h
        gp = calc(x + d)[1]
        gm = calc(x - d)[1]
        H[i] = (gp - gm) / (2 * h)
        if ckpt is not None:
            np.save(ckpt, H)
    return H


def mw_modes(H, x, m):
    """Normal modes of a Cartesian Hessian with translations and rotations projected out: frequencies in cm-1 (negative = imaginary),
    columns of L in mass-weighted coordinates."""
    H = 0.5 * (H + H.T)
    sm = np.sqrt(m)
    Hm = H / np.outer(sm, sm)
    q = x * sm
    P = project(q, m)
    lam, L = np.linalg.eigh(P @ Hm @ P)
    # hartree/(bohr^2 amu) -> cm-1
    conv = 5140.4871
    keep = np.argsort(np.abs(lam))[6:]
    keep = keep[np.argsort(lam[keep])]
    return np.sign(lam[keep]) * np.sqrt(np.abs(lam[keep])) * conv, L[:, keep], lam[keep]


def part_hess(entry):
    calc = Calc(entry)
    x = coords(entry, "TS1").ravel() / BOHR
    m = masses(entry)
    e0, g0 = calc(x)
    t0 = time.time()
    H = numeric_hessian(calc, x, ckpt=OUT / f"hess_ts_{entry}.npy")
    w, L, lam = mw_modes(H, x, m)
    print(f"entry {entry}: E {e0:.6f} |g|max {np.abs(g0).max():.5f}; {calc.calls} gradients, {time.time() - t0:.0f} s; "
          f"asymmetry {np.abs(H - H.T).max():.2e}")
    print("frequencies", np.round(w, 1).tolist())
    np.savez(OUT / f"ts_{entry}.npz", x=x, m=m, E=e0, g=g0, H=H, w=w, L=L, lam=lam)


def lqa_step(g, lam, V, ds):
    """Local quadratic approximation (Page and McIver): the step of arc length ds along the steepest-descent line of the quadratic
    model with mass-weighted gradient g and Hessian V diag(lam) V^T."""
    gi = V.T @ g

    def dq(t):
        small = np.abs(lam * t) < 1e-8
        f = np.where(small, -t, (np.exp(-lam * t) - 1.0) / np.where(small, 1.0, lam))
        return V @ (f * gi)

    ts = np.geomspace(1e-4, 1e4, 1200)
    speed = np.array([np.linalg.norm(V @ (np.exp(-lam * t) * gi)) for t in ts])
    arc = np.r_[0.0, np.cumsum(0.5 * (speed[1:] + speed[:-1]) * np.diff(ts))]
    if arc[-1] < ds:
        return dq(ts[-1])
    return dq(float(np.interp(ds, arc, ts)))


def bofill(H, dq, dg):
    r = dg - H @ dq
    dd, rd, rr = dq @ dq, r @ dq, r @ r
    if dd < 1e-14 or rr < 1e-20:
        return H
    phi = rd * rd / (rr * dd)
    ms = np.outer(r, r) / rd if abs(rd) > 1e-14 else 0 * H
    psb = (np.outer(r, dq) + np.outer(dq, r)) / dd - rd * np.outer(dq, dq) / dd ** 2
    return H + phi * ms + (1 - phi) * psb


def part_irc(entry, ds=0.1, steps=500):
    """IRC from TS1 to the product side, in mass-weighted coordinates (arc step ds in sqrt(amu) bohr)."""
    ts = np.load(OUT / f"ts_{entry}.npz")
    x0, m, L, H0 = ts["x"], ts["m"], ts["L"], ts["H"]
    sm = np.sqrt(m)
    u = L[:, 0]
    a, b = BOND[entry]
    r0 = x0.reshape(-1, 3)
    xp = (x0 + 0.05 * u / sm).reshape(-1, 3)
    if np.linalg.norm(xp[a] - xp[b]) > np.linalg.norm(r0[a] - r0[b]):
        u = -u
    calc = Calc(entry)
    q = x0 * sm + 0.05 * u
    Hm = 0.5 * (H0 + H0.T) / np.outer(sm, sm)
    Q, E, G = [x0 * sm], [float(ts["E"])], [ts["g"] / sm]
    s = [0.0]
    gprev = None
    qprev = None
    t0 = time.time()
    for k in range(steps):
        e, gc = calc(q / sm)
        g = gc / sm
        P = project(q, m)
        g = P @ g
        if gprev is not None:
            Hm = bofill(Hm, q - qprev, g - gprev)
        Q.append(q.copy())
        E.append(e)
        G.append(g.copy())
        s.append(s[-1] + np.linalg.norm(q - Q[-2]))
        np.savez(OUT / f"irc_{entry}.npz", q=np.array(Q), E=np.array(E), g=np.array(G), s=np.array(s), m=m, Hm=Hm)
        if k % 10 == 0:
            d = np.linalg.norm((q / sm).reshape(-1, 3)[a] - (q / sm).reshape(-1, 3)[b]) * BOHR
            print(f"step {k} s {s[-1]:.2f} dE {(e - E[0]) * 627.509:.2f} kcal/mol |g| {np.linalg.norm(g):.5f} R {d:.3f} "
                  f"{time.time() - t0:.0f} s", flush=True)
        if np.linalg.norm(g) < 3e-4 or (k > 20 and E[-1] > E[-2] + 1e-6):
            break
        lam, V = np.linalg.eigh(P @ (0.5 * (Hm + Hm.T)) @ P)
        qprev, gprev = q.copy(), g.copy()
        q = q + lqa_step(g, lam, V, ds)


def transverse_curvature(Hm, tangent, q, m):
    """Lowest eigenvalue and eigenvector of the Hessian with the path tangent, translations and rotations projected out."""
    t = tangent / np.linalg.norm(tangent)
    P = project(q, m)
    P = P - np.outer(P @ t, P @ t)
    lam, V = np.linalg.eigh(P @ (0.5 * (Hm + Hm.T)) @ P)
    idx = np.argsort(np.abs(lam))[7:]
    idx = idx[np.argsort(lam[idx])]
    return lam[idx[0]], V[:, idx[0]], lam[idx], V[:, idx]


def part_vrt(entry):
    """The point of the IRC with the softest transverse curvature, from the updated Hessians it is not available, so an exact Hessian is
    computed at the geometry of the IRC nearest the VRT of Yamaki et al. Its lowest transverse mode is the coordinate y."""
    irc = np.load(OUT / f"irc_{entry}.npz")
    q, m, Es, s = irc["q"], irc["m"], irc["E"], irc["s"]
    sm = np.sqrt(m)
    target = align(coords(entry, "VRT").ravel() / BOHR, q[0] / sm, m) * sm
    k = int(np.argmin(np.linalg.norm(q - target, axis=1)))
    print(f"nearest IRC point to the published VRT: index {k}, s {s[k]:.3f}, distance {np.linalg.norm(q[k] - target):.3f} sqrt(amu) bohr",
          flush=True)
    calc = Calc(entry)
    x = q[k] / sm
    calc(x)
    H = numeric_hessian(calc, x, ckpt=OUT / f"hess_vrt_{entry}.npy")
    Hm = 0.5 * (H + H.T) / np.outer(sm, sm)
    t = q[min(k + 1, len(q) - 1)] - q[max(k - 1, 0)]
    lam0, y, lam, V = transverse_curvature(Hm, t, q[k], m)
    conv = 5140.4871
    print("lowest transverse curvatures (cm-1, negative imaginary)", np.round(np.sign(lam[:5]) * np.sqrt(np.abs(lam[:5])) * conv, 1).tolist())
    np.savez(OUT / f"vrt_{entry}.npz", k=k, y=y, lam=lam, V=V, H=H)


def align(x, xref, m):
    """Mass-weighted Kabsch rotation of x onto xref (both Cartesian vectors), translation removed."""
    a, b = x.reshape(-1, 3), xref.reshape(-1, 3)
    w = m[::3, None]
    ca, cb = (w * a).sum(0) / w.sum(), (w * b).sum(0) / w.sum()
    A, B = a - ca, b - cb
    U, _, Vt = np.linalg.svd((w * A).T @ B)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return (A @ R + cb).ravel()


def hvp_rows(calc, q, y, m, eps=0.02):
    """H y at q from a central difference of the mass-weighted gradient along y, and the energy at q."""
    sm = np.sqrt(m)
    e0, _ = calc(q / sm)
    gp = calc((q + eps * y) / sm)[1] / sm
    gm = calc((q - eps * y) / sm)[1] / sm
    return e0, (gp - gm) / (2 * eps)


def part_scan(entry, every=4):
    """Curvature along y on the relaxed path, and the row H y, every few IRC points."""
    irc = np.load(OUT / f"irc_{entry}.npz")
    v = np.load(OUT / f"vrt_{entry}.npz")
    q, m, s = irc["q"], irc["m"], irc["s"]
    y = v["y"]
    calc = Calc(entry)
    idx = list(range(0, len(q), every))
    rows = []
    for j, k in enumerate(idx):
        e, Hy = hvp_rows(calc, q[k], y, m)
        rows.append(Hy)
        print(f"{k} s {s[k]:.2f} y'Hy {y @ Hy:.5f}", flush=True)
        np.savez(OUT / f"scan_{entry}.npz", idx=np.array(idx[:j + 1]), Hy=np.array(rows), y=y)


def part_traj(entry, dt_fs=0.5, t_end_fs=220.0, kT=9.44e-4):
    """Velocity Verlet from TS1 with the mean flux speed along the imaginary mode (Rayleigh mean, 298 K) and every other mode at rest:
    the motion of each mode with no thermal noise, to compare with its position on the IRC."""
    ts = np.load(OUT / f"ts_{entry}.npz")
    x0, m, L = ts["x"], ts["m"], ts["L"]
    sm = np.sqrt(m)
    u = L[:, 0]
    a, b = BOND[entry]
    xp = (x0 + 0.05 * u / sm).reshape(-1, 3)
    r0 = x0.reshape(-1, 3)
    if np.linalg.norm(xp[a] - xp[b]) > np.linalg.norm(r0[a] - r0[b]):
        u = -u
    dt = dt_fs / 1.03275                          # sqrt(amu) bohr / sqrt(hartree) time unit in fs
    q = x0 * sm
    v = np.sqrt(np.pi * kT / 2) * u
    calc = Calc(entry)
    e, gc = calc(q / sm)
    g = gc / sm
    Q, V_, E = [q.copy()], [v.copy()], [e]
    t0 = time.time()
    for k in range(int(round(t_end_fs / dt_fs))):
        v = v - 0.5 * dt * g
        q = q + dt * v
        e, gc = calc(q / sm)
        g = gc / sm
        v = v - 0.5 * dt * g
        Q.append(q.copy())
        V_.append(v.copy())
        E.append(e)
        if k % 20 == 0:
            d = np.linalg.norm((q / sm).reshape(-1, 3)[a] - (q / sm).reshape(-1, 3)[b]) * BOHR
            np.savez(OUT / f"traj_{entry}.npz", q=np.array(Q), v=np.array(V_), E=np.array(E), dt_fs=dt_fs)
            print(f"t {k * dt_fs:.1f} fs R {d:.3f} dE {(e - E[0]) * 627.509:.2f} K {0.5 * (v @ v) * 627.509:.2f} kcal/mol "
                  f"total {(e + 0.5 * v @ v - E[0] - 0.5 * V_[0] @ V_[0]) * 627.509:.3f} {time.time() - t0:.0f} s", flush=True)
    np.savez(OUT / f"traj_{entry}.npz", q=np.array(Q), v=np.array(V_), E=np.array(E), dt_fs=dt_fs)


def part_hessat(entry, ks):
    """Exact Hessians at chosen IRC points: the true lowest curvature across the path there, with no direction assumed."""
    irc = np.load(OUT / f"irc_{entry}.npz")
    q, m, s = irc["q"], irc["m"], irc["s"]
    sm = np.sqrt(m)
    calc = Calc(entry)
    for k in ks:
        x = q[k] / sm
        calc(x)
        H = numeric_hessian(calc, x, ckpt=OUT / f"hess_{entry}_{k}.npy")
        Hm = 0.5 * (H + H.T) / np.outer(sm, sm)
        t = q[min(k + 1, len(q) - 1)] - q[max(k - 1, 0)]
        lam0, y, lam, V = transverse_curvature(Hm, t, q[k], m)
        print(f"k {k} s {s[k]:.2f} lowest (cm-1) {np.round(np.sign(lam[:4]) * np.sqrt(np.abs(lam[:4])) * 5140.4871, 1).tolist()}", flush=True)
        np.savez(OUT / f"hessat_{entry}_{k}.npz", k=k, s=s[k], lam=lam, V=V[:, :3], asym=np.abs(H - H.T).max())


def part_hessframe(entry, times):
    """Exact Hessians at frames of the trajectory: the curvatures the molecule meets where it actually is, across its velocity."""
    tr = np.load(OUT / f"traj_{entry}.npz")
    m = np.load(OUT / f"irc_{entry}.npz")["m"]
    sm = np.sqrt(m)
    a, b = BOND[entry]
    calc = Calc(entry)
    for t in times:
        i = int(round(t / float(tr["dt_fs"])))
        q, v = tr["q"][i], tr["v"][i]
        x = q / sm
        e, gc = calc(x)
        H = numeric_hessian(calc, x, ckpt=OUT / f"hess_{entry}_t{t:g}.npy")
        Hm = 0.5 * (H + H.T) / np.outer(sm, sm)
        lam0, y, lam, V = transverse_curvature(Hm, v, q, m)
        r = x.reshape(-1, 3)
        R = np.linalg.norm(r[a] - r[b]) * BOHR
        print(f"t {t:g} fs R {R:.3f} lowest (cm-1) {np.round(np.sign(lam[:4]) * np.sqrt(np.abs(lam[:4])) * 5140.4871, 1).tolist()}", flush=True)
        np.savez(OUT / f"hessframe_{entry}_{t:g}.npz", t=t, R=R, E=e, g=gc / sm, lam=lam, V=V[:, :3], asym=np.abs(H - H.T).max())


def part_trajy(entry, displacements, t_end_fs=130.0, dt_fs=0.5, kT=9.44e-4, vy=0.0):
    """Trajectories from TS1 displaced along the slow transverse direction, in units of its thermal width sqrt(kT) / w. The displacement is made in
    internal coordinates, as a turn is not a straight line, and the initial velocity is the same flux along the imaginary mode, also in internal coordinates.
    vy adds a velocity along the same direction in units of its thermal width sqrt(kT). A file that already holds a shorter run is continued."""
    from internal import Curvilinear
    from yamaki_geometries import atoms
    ts = np.load(OUT / f"ts_{entry}.npz")
    irc = np.load(OUT / f"irc_{entry}.npz")
    y = np.load(OUT / f"vrt_{entry}.npz")["y"]
    m, L, q0, w = ts["m"], ts["L"], irc["q"][0], np.abs(ts["w"]) / 5140.4871
    sm = np.sqrt(m)
    c = L[:, 1:3].T @ y
    c = c / np.linalg.norm(c)
    sigma = np.sqrt(kT * ((c / w[1:3]) ** 2).sum())
    C = Curvilinear([a[0] for a in atoms(entry, "TS1")], [coords(entry, n) for n in ("TS1", "PM", "P1", "P2")], BOND[entry], q0, m, L)
    a, b = BOND[entry]
    calc = Calc(entry)
    e0 = calc(q0 / sm)[0]
    dt = dt_fs / 1.03275
    for d in displacements:
        target = np.zeros(L.shape[1])
        target[1:3] = c * d * sigma
        q = q0.copy()
        for it in range(60):
            res = target - C(q)
            if np.linalg.norm(res) < 1e-5:
                break
            q = q + L @ res
        # velocity: the flux along the imaginary mode, in internal coordinates at the displaced geometry
        B = np.zeros((len(C.xi0), len(q)))
        for i in range(len(q)):
            h = np.zeros(len(q))
            h[i] = 1e-3
            B[:, i] = C.wrap(C.xi(q + h) - C.xi(q - h)) / 2e-3
        u = np.linalg.pinv(B, rcond=1e-6) @ C.BL[:, 0]
        P = project(q, m)
        u = P @ u
        u = u / np.linalg.norm(u)
        r0 = (q / sm).reshape(-1, 3)
        rp = ((q + 0.05 * u) / sm).reshape(-1, 3)
        if np.linalg.norm(rp[a] - rp[b]) > np.linalg.norm(r0[a] - r0[b]):
            u = -u
        v = np.sqrt(np.pi * kT / 2) * u + vy * np.sqrt(kT) * (L[:, 1:3] @ c)
        name = OUT / (f"trajy_{entry}_{d:+.1f}.npz" if vy == 0 else f"trajy_{entry}_{d:+.1f}_v{vy:+.1f}.npz")
        start = 0
        if name.exists():
            o = np.load(name)
            q, v, start = o["q"][-1], o["v"][-1], len(o["q"]) - 1
        e, gc = calc(q / sm)
        g = gc / sm
        print(f"d {d:+.2f} vy {vy:+.1f} from frame {start}, sigma {sigma:.2f}: {it} iterations, residual {np.linalg.norm(res):.1e}, energy above TS1 {(e - e0) * 627.509:.2f} kcal/mol, harmonic {0.5 * kT * d ** 2 * 627.509:.2f}", flush=True)
        Q, V_, E = ([q.copy()], [v.copy()], [e]) if start == 0 else (list(o["q"]), list(o["v"]), list(o["E"]))
        t0 = time.time()
        for k in range(start, int(round(t_end_fs / dt_fs))):
            v = v - 0.5 * dt * g
            q = q + dt * v
            e, gc = calc(q / sm)
            g = gc / sm
            v = v - 0.5 * dt * g
            Q.append(q.copy())
            V_.append(v.copy())
            E.append(e)
            if k % 20 == 19 or k == int(round(t_end_fs / dt_fs)) - 1:
                np.savez(name, d=d, vy=vy, sigma=sigma, q=np.array(Q), v=np.array(V_), E=np.array(E), dt_fs=dt_fs, e0=e0)
                xi = C.xi(q)
                print(f"d {d:+.1f} t {(k + 1) * dt_fs:.0f} R {np.linalg.norm((q / sm).reshape(-1, 3)[a] - (q / sm).reshape(-1, 3)[b]) * BOHR:.3f} total {(e + 0.5 * v @ v - E[0] - 0.5 * V_[0] @ V_[0]) * 627.509:.3f} {time.time() - t0:.0f} s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("part", choices=["hess", "irc", "vrt", "scan", "traj", "hessat", "hessframe", "trajy"])
    ap.add_argument("entry", type=int)
    ap.add_argument("--k", default="")
    ap.add_argument("--vy", type=float, default=0.0)
    ap.add_argument("--tend", type=float, default=130.0)
    args = ap.parse_args()
    if args.part == "trajy":
        OUT.mkdir(exist_ok=True)
        part_trajy(args.entry, [float(v) for v in args.k.split(",")], t_end_fs=args.tend, vy=args.vy)
        return
    if args.part == "hessframe":
        OUT.mkdir(exist_ok=True)
        part_hessframe(args.entry, [float(v) for v in args.k.split(",")])
        return
    if args.part == "hessat":
        OUT.mkdir(exist_ok=True)
        part_hessat(args.entry, [int(v) for v in args.k.split(",")])
        return
    OUT.mkdir(exist_ok=True)
    {"hess": part_hess, "irc": part_irc, "vrt": part_vrt, "scan": part_scan, "traj": part_traj}[args.part](args.entry)


if __name__ == "__main__":
    main()
