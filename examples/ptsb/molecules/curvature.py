"""The curvature across the motion, on the IRC and on the trajectory, from exact Hessians; and the frequency of each mode along the path.

    python curvature.py 31     # reads results/{ts,irc,vrt,traj,modes}_31.*, hess_31_*.npy, hessat_31_*.npz, hessframe_31_*.npz; writes results/curvature_31.json

Across the motion means with the path tangent (on the IRC) or the velocity (on the trajectory) projected out, and translations and rotations.
A slow mode that stays behind leaves the molecule on a surface that is not the relaxed one. If that matters, the curvature met on the trajectory
differs from the curvature on the IRC at the same forming bond.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from internal import Curvilinear  # noqa: E402
from modes import TORSION, W_CM, TIME_FS, bond_length  # noqa: E402
from qm import BOND, OUT, transverse_curvature  # noqa: E402
from yamaki_geometries import atoms, coords  # noqa: E402


def signed_cm(lam):
    return float(np.sign(lam) * np.sqrt(abs(lam)) * W_CM)


def torsion_direction(C, q, tors):
    """Unit vector in mass-weighted coordinates along which the torsion changes fastest at q."""
    n = C.nb + C.na + tors
    g = np.zeros(len(q))
    for i in range(len(q)):
        e = np.zeros(len(q))
        e[i] = 1e-3
        g[i] = C.wrap(C.xi(q + e) - C.xi(q - e))[n] / 2e-3
    return g / np.linalg.norm(g), float(np.linalg.norm(g))


def response(a, lam, t):
    """Q'' = -lam(t) (Q - a(t)) from rest, by velocity Verlet on the frames. lam may change sign."""
    Q, P = np.zeros_like(a), 0.0
    for i in range(len(t) - 1):
        h = t[i + 1] - t[i]
        P += -0.5 * h * lam[i] * (Q[i] - a[i])
        Q[i + 1] = Q[i] + h * P
        P += -0.5 * h * lam[i + 1] * (Q[i + 1] - a[i + 1])
    return Q


def main(entry):
    ts = np.load(OUT / f"ts_{entry}.npz")
    irc = np.load(OUT / f"irc_{entry}.npz")
    tr = np.load(OUT / f"traj_{entry}.npz")
    mo = np.load(OUT / f"modes_{entry}.npz")
    vrt = np.load(OUT / f"vrt_{entry}.npz")
    m, L, wcm = ts["m"], ts["L"], ts["w"]
    sm = np.sqrt(m)
    q, s, E = irc["q"], irc["s"], irc["E"]
    R = mo["R"]
    elements = [a[0] for a in atoms(entry, "TS1")]
    C = Curvilinear(elements, [coords(entry, n) for n in ("TS1", "PM", "P1", "P2")], BOND[entry], q[0], m, L)
    tors = C.g[2].index(TORSION[entry])
    y = vrt["y"]

    # the IRC
    path, Hs = [], {}
    pts = [(0, OUT / f"hess_ts_{entry}.npy"), (int(vrt["k"]), OUT / f"hess_vrt_{entry}.npy")]
    pts += [(int(p.stem.split("_")[-1]), p) for p in OUT.glob(f"hess_{entry}_[0-9]*.npy")]
    for k, f in sorted(set(pts)):
        H = np.load(f)
        Hm = 0.5 * (H + H.T) / np.outer(sm, sm)
        Hs[k] = Hm
        t = q[min(k + 1, len(q) - 1)] - q[max(k - 1, 0)]
        lam0, v0, lam, V = transverse_curvature(Hm, t, q[k], m)
        d, _ = torsion_direction(C, q[k], tors)
        path.append(dict(k=k, s=float(s[k]), R=float(R[k]), E_kcal=float((E[k] - E[0]) * 627.509), lowest_cm=[signed_cm(v) for v in lam[:3]],
                         overlap_y=float((v0 @ y) ** 2), overlap_torsion=float((v0 @ d) ** 2), overlap_s_q=float(abs(v0 @ L[:, 1])),
                         mode_shares=[float(v) for v in ((L.T @ v0) ** 2)[:12]]))

    # where the lowest curvature across the path passes zero, by interpolation in the eigenvalue
    lamk = np.array([np.sign(p["lowest_cm"][0]) * (p["lowest_cm"][0] / W_CM) ** 2 for p in path])
    sk, Rk = np.array([p["s"] for p in path]), np.array([p["R"] for p in path])
    cross = int(np.argmax(lamk < 0))
    fz = lamk[cross - 1] / (lamk[cross - 1] - lamk[cross])
    s_zero, R_zero = sk[cross - 1] + fz * (sk[cross] - sk[cross - 1]), Rk[cross - 1] + fz * (Rk[cross] - Rk[cross - 1])

    # the trajectory
    frames = []
    dt_fs = float(tr["dt_fs"])
    for f in sorted(OUT.glob(f"hessframe_{entry}_*.npz"), key=lambda p: float(p.stem.split("_")[-1])):
        d = np.load(f)
        i = int(round(float(d["t"]) / dt_fs))
        Hc = np.load(OUT / f"hess_{entry}_t{float(d['t']):g}.npy")
        Hm = 0.5 * (Hc + Hc.T) / np.outer(sm, sm)
        v = tr["v"][i]
        lam0, v0, lam, V = transverse_curvature(Hm, v, tr["q"][i], m)
        dt, _ = torsion_direction(C, tr["q"][i], tors)
        g = d["g"]
        vhat = v / np.linalg.norm(v)
        gt = g - vhat * (vhat @ g)
        Rf = float(d["R"])
        # the same bond length on the IRC before the bond settles
        kR = int(mo["kR"])
        lk = np.interp(Rf, Rk[::-1], lamk[::-1]) if Rk.min() <= Rf <= Rk.max() else np.nan
        frames.append(dict(t=float(d["t"]), R=Rf, lowest_cm=[signed_cm(x) for x in lam[:3]], overlap_torsion=float((v0 @ dt) ** 2),
                           overlap_y=float((v0 @ y) ** 2), side_force=float(np.linalg.norm(gt)), force_on_lowest=float(abs(v0 @ g)),
                           irc_lowest_at_same_R_cm=None if np.isnan(lk) else signed_cm(lk), mode_shares=[float(x) for x in ((L.T @ v0) ** 2)[:12]],
                           torsion_deg=float(mo["xi_t"][i]), V_kcal=float((tr["E"][i] - tr["E"][0]) * 627.509)))

    # the frequency of each TS1 mode along the path, against the frozen frequency, in the linear response of each mode
    keys = sorted(Hs)
    kR, end = int(mo["kR"]), int(mo["end"])
    near = mo["near"]
    tt = mo["tt"][: end + 1] / TIME_FS
    a = mo["a"]
    Q = mo["Q"]
    rows = []
    for j in range(1, L.shape[1]):
        lam_k = np.array([L[:, j] @ Hs[k] @ L[:, j] for k in keys])
        lam_t = np.interp(near, keys, lam_k)
        aj = a[near, j]
        den = float(aj @ aj)
        F = float(Q[: end + 1, j] @ aj / den)
        Ql = response(aj, lam_t, tt)
        k_end = int(np.argmin(np.abs(np.array(keys) - kR)))
        rows.append(dict(j=j, w=float(wcm[j]), a_R=float(a[kR, j]), F=F, F_lin_path=float(Ql @ aj / den), w_at_bond_cm=signed_cm(lam_k[k_end]),
                         w_along_cm=[signed_cm(v) for v in lam_k]))
    res = dict(entry=entry, keys=keys, s_zero=float(s_zero), R_zero=float(R_zero), path=path, frames=frames, modes=rows)
    (OUT / f"curvature_{entry}.json").write_text(json.dumps(res, indent=1))

    print(f"entry {entry}: lowest curvature across the IRC passes zero at s {s_zero:.2f}, R {R_zero:.3f}")
    print("IRC  k    s     R    E kcal   lowest cm-1            ov_y ov_tors")
    for p in path:
        print(f"{p['k']:4d} {p['s']:5.2f} {p['R']:5.3f} {p['E_kcal']:7.2f}  {np.round(p['lowest_cm'], 1)}  {p['overlap_y']:.2f} {p['overlap_torsion']:.2f}")
    print("traj t fs   R    lowest cm-1   IRC at same R  ov_tors ov_y  torsion  force on lowest")
    for f in frames:
        print(f"{f['t']:5.0f} {f['R']:5.3f} {np.round(f['lowest_cm'], 1)} {f['irc_lowest_at_same_R_cm']} {f['overlap_torsion']:.2f} {f['overlap_y']:.2f} {f['torsion_deg']:.0f} {f['force_on_lowest']:.4f}")
    print("mode cm-1 a_R F F_lin_path w_at_bond")
    moving = [r for r in rows if abs(r["a_R"]) >= 0.05]
    for r in moving:
        print(f"{r['j']:3d} {r['w']:7.0f} {r['a_R']:7.3f} {r['F']:6.2f} {r['F_lin_path']:6.2f} {r['w_at_bond_cm']:7.0f}")
    print("mean abs F - F_lin_path, modes moving", float(np.mean([abs(r["F"] - r["F_lin_path"]) for r in moving])))


if __name__ == "__main__":
    main(int(sys.argv[1]))
