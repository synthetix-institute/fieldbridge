"""Which modes of a molecule follow the reaction path, which stay behind and which lag, from DFT.

    python modes.py 31         # reads results/{ts,irc,vrt,traj}_31.npz; writes results/modes_31.npz and results/modes_31.json

The modes are those of TS1, read in internal coordinates (internal.py), because in straight-line coordinates a turning group
counts as a stretch of its bonds. Progress is the length R of the forming bond, the counterpart of x in the model of ../anharm.py:
the trajectory at time t is compared with the IRC point of the same R, up to the first time R reaches its value on the IRC after
the bond has formed (the passage, duration T).

For each TS1 mode j (frequency w_j):
  a_j   displacement of the mode on the IRC at the same R
  Q_j   displacement on a noise-free DFT trajectory from TS1
  F_j   the following fraction over the passage, the regression slope of Q_j on a_j
  F_lin the same from the linear response of a harmonic mode forced by a_j(t), Q'' = -w^2 (Q - a)
F_j near 1 is a mode that follows, near 0 a mode that stays behind. These are the two limits of the model.
g_j is the size of the coupling in the units of that model: x = 1 at the pseudo-minimum, frequency 1 at the transverse mode of TS1.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from internal import Curvilinear  # noqa: E402
from qm import BOHR, BOND, OUT, align  # noqa: E402
from yamaki_geometries import RATIO, atoms, coords  # noqa: E402

KCAL = 627.509
TIME_FS = 1.03275                                 # sqrt(amu) bohr / sqrt(hartree) in fs
NEG = 0.05                                        # sqrt(amu) bohr: a mode moving less than this on the IRC is not counted
W_CM = 5140.4871
KT = 9.44e-4
ROWS = ("stays behind", "between", "follows", "past the path", "transverse", "uninvolved")
TORSION = {31: (6, 0, 2, 3), 32: (1, 0, 2, 3)}    # the turn of the CHR end of the allyl group, which separates P1 from P2


def linear_response(a, t, w):
    """Harmonic mode of frequency w (per unit time, same units as t) forced by a(t), at rest at the start, exact on each interval."""
    Q, P = np.zeros_like(a), np.zeros_like(a)
    for i in range(len(t) - 1):
        h = t[i + 1] - t[i]
        am = 0.5 * (a[i] + a[i + 1])
        c, s = np.cos(w * h), np.sin(w * h)
        y = Q[i] - am
        Q[i + 1] = am + y * c + P[i] / w * s
        P[i + 1] = -w * y * s + P[i] * c
    return Q


def bond_length(q, m, pair):
    x = (q / np.sqrt(m)).reshape(len(q), -1, 3)
    return np.linalg.norm(x[:, pair[0]] - x[:, pair[1]], axis=1) * BOHR


def main(entry):
    ts = np.load(OUT / f"ts_{entry}.npz")
    irc = np.load(OUT / f"irc_{entry}.npz")
    vrt = np.load(OUT / f"vrt_{entry}.npz")
    tr = np.load(OUT / f"traj_{entry}.npz")
    m, L, wcm = ts["m"], ts["L"], ts["w"]
    sm = np.sqrt(m)
    w = np.abs(wcm) / W_CM                        # per sqrt(amu) bohr / sqrt(hartree)
    q, E, s = irc["q"], irc["E"], irc["s"]
    q0 = q[0]
    elements = [a[0] for a in atoms(entry, "TS1")]
    C = Curvilinear(elements, [coords(entry, n) for n in ("TS1", "PM", "P1", "P2")], BOND[entry], q0, m, L)
    a = np.array([C(p) for p in q])               # (points, modes); mode 0 is the imaginary mode
    a_cart = (q - q0) @ L
    over = (L.T @ vrt["y"]) ** 2
    jy = int(np.argmax(over[1:]) + 1)

    pm = align(coords(entry, "PM").ravel() / BOHR, q0 / sm, m) * sm
    vr = align(coords(entry, "VRT").ravel() / BOHR, q0 / sm, m) * sm
    kpm = int(np.argmin(np.linalg.norm(q - pm, axis=1)))
    kv = int(np.argmin(np.linalg.norm(q - vr, axis=1)))

    # progress by the forming bond
    R = bond_length(q, m, BOND[entry])
    R_stop = R[kpm] + 0.01
    kR = int(np.argmax(R <= R_stop))
    Rt = bond_length(tr["q"], m, BOND[entry])
    end = int(np.argmax(Rt <= R_stop))
    near = np.array([int(np.argmin(np.abs(R[:kR + 1] - r))) for r in Rt[:end + 1]])

    dt = float(tr["dt_fs"]) / TIME_FS
    tt = np.arange(len(Rt)) * dt
    Q = np.array([C(f) for f in tr["q"]])
    P = np.gradient(Q, dt, axis=0)
    A = a[near]
    Tp = tt[end]

    rows = []
    for j in range(L.shape[1]):
        aj, Qj = A[:, j], Q[:end + 1, j]
        den = float(aj @ aj)
        Ql = linear_response(aj, tt[:end + 1], w[j])
        half = np.nonzero(np.abs(Q[:, j]) >= 0.5 * abs(a[kR, j]))[0] if abs(a[kR, j]) >= NEG else []
        rows.append(dict(j=j, w=float(wcm[j]), period_fs=float(2 * np.pi / w[j] * TIME_FS), omega_T=float(w[j] * Tp),
                         a_passage=float(a[kR, j]), a_pm=float(a[kpm, j]), a_pm_cartesian=float(a_cart[kpm, j]), Q_end=float(Q[end, j]),
                         F=float(Qj @ aj / den), F_lin=float(Ql @ aj / den), residual=float(np.sqrt(np.mean((Qj - Qj @ aj / den * aj) ** 2))),
                         half_way_fs=float(tt[half[0]] * TIME_FS) if len(half) else None,
                         overlap_y=float(over[j]), thermal_width=float(np.sqrt(KT) / w[j]),
                         energy_held_kcal=float((0.5 * P[end, j] ** 2 + 0.5 * w[j] ** 2 * (Q[end, j] - A[end, j]) ** 2) * KCAL),
                         relaxation_energy_kcal=float(0.5 * w[j] ** 2 * a[kpm, j] ** 2 * KCAL)))
    modes = rows[1:]
    def row_of(r):
        if r["j"] == jy:
            return "transverse"
        if abs(r["a_passage"]) < NEG:
            return "uninvolved"
        if r["F"] < 0.3:
            return "stays behind"
        if r["F"] <= 0.7:
            return "between"
        return "follows" if r["F"] <= 1.3 else "past the path"
    for r in modes:
        r["row"] = row_of(r)
    wy = w[jy]
    X = abs(a[kpm, 0])
    for r in modes:
        r["g"] = float(2 * abs(r["a_pm"]) / X * (w[r["j"]] / wy) ** 1.5)
    moving = [r for r in modes if r["row"] not in ("uninvolved", "transverse")]

    # the inertia that following modes would add, sum (d a_j)^2 / sum (d a_0)^2 over the passage
    da = np.diff(a[:kR + 1], axis=0)
    inertia = (da ** 2).sum(0) / (da[:, 0] ** 2).sum()

    # how far the trajectory is above the path at the same R, against the harmonic estimate from the modes
    dV = (tr["E"][:end + 1] - np.interp(Rt[:end + 1], R[:kR + 1][::-1], E[:kR + 1][::-1])) * KCAL
    harm = 0.5 * ((w[1:] * (Q[:end + 1, 1:] - A[:, 1:])) ** 2).sum(1) * KCAL
    g_m = [r["g"] for r in modes]
    band = [r["g"] for r in moving if abs(r["w"]) / abs(wcm[jy]) <= 4]

    tors = C.g[2].index(TORSION[entry])
    xi_t = np.degrees(np.array([C.xi(f)[C.nb + C.na + tors] for f in tr["q"]]))
    xi_p = np.degrees(np.array([C.xi(f)[C.nb + C.na + tors] for f in q]))
    summary = {
        "entry": entry, "ratio_published": RATIO[entry], "imaginary_cm": float(wcm[0]), "transverse_mode": jy,
        "transverse_cm": float(wcm[jy]), "transverse_overlap": float(over[jy]),
        "s_pm": float(s[kpm]), "s_bond_formed": float(s[kR]), "s_vrt": float(s[kv]),
        "drop_ts_to_pm_kcal": float((E[0] - E[kpm]) * KCAL), "drop_ts_to_vrt_kcal": float((E[0] - E[kv]) * KCAL),
        "passage_fs": float(Tp * TIME_FS), "R_stop": float(R_stop), "R_min_first_fs": float(np.argmin(Rt[:end + 40]) * dt * TIME_FS), "R_min": float(Rt[:end + 40].min()),
        "omega_y_T": float(wy * Tp), "w_imag_over_w_y": float(w[0] / wy), "X": float(X),
        "energy_drop_over_wy2X2": float((E[0] - E[kpm]) / (wy ** 2 * X ** 2)),
        "thermal_width_y_over_X": float(np.sqrt(KT) / wy / X),
        "inertia_all": float(inertia[1:].sum()), "inertia_following_rows": float(sum(inertia[r["j"]] for r in modes if r["row"] in ("follows", "past the path"))),
        "dV_above_path_at_end_kcal": float(dV[end]), "dV_above_path_max_kcal": float(dV.max()),
        "harmonic_estimate_at_end_kcal": float(harm[end]), "kinetic_at_end_kcal": float(0.5 * tr["v"][end] @ tr["v"][end] * KCAL),
        "energy_held_sum_kcal": float(sum(r["energy_held_kcal"] for r in modes)),
        "energy_held_by_row_kcal": {k: float(sum(r["energy_held_kcal"] for r in modes if r["row"] == k)) for k in ROWS},
        "relaxation_energy_by_row_kcal": {k: float(sum(r["relaxation_energy_kcal"] for r in modes if r["row"] == k)) for k in ROWS},
        "count_by_row": {k: sum(r["row"] == k for r in modes) for k in ROWS},
        "displacement_share_by_row": {k: float(sum(r["a_pm"] ** 2 for r in modes if r["row"] == k) / sum(r["a_pm"] ** 2 for r in modes)) for k in ROWS},
        "g_median_moving": float(np.median([r["g"] for r in moving])), "g_max": float(max(g_m)),
        "g_band_up_to_4_wy": band,
        "F_vs_F_lin_mean_abs": float(np.mean([abs(r["F"] - r["F_lin"]) for r in moving])),
        "F_vs_F_lin_low_modes_mean_abs": float(np.mean([abs(r["F"] - r["F_lin"]) for r in moving if abs(r["w"]) < 600])),
        "torsion_ts_deg": float(xi_p[0]), "torsion_pm_deg": float(xi_p[kpm]),
        "torsion_traj_deg": {str(t): float(xi_t[int(t / float(tr["dt_fs"]))]) for t in (10, 20, 30, 40, 50, 60, 70, 80, 100, 150, 200)},
    }
    np.savez(OUT / f"modes_{entry}.npz", near=near, end=end, kv=kv, kpm=kpm, kR=kR, a=a, Q=Q, tt=tt * TIME_FS, wcm=wcm, R=R, Rt=Rt, s=s,
             xi_t=xi_t, xi_p=xi_p, dV=dV, harm=harm, F=np.array([r["F"] for r in rows]), Fl=np.array([r["F_lin"] for r in rows]))
    (OUT / f"modes_{entry}.json").write_text(json.dumps({"summary": summary, "modes": modes}, indent=1))
    print(json.dumps(summary, indent=0))
    print(f"{'j':>3} {'cm-1':>7} {'T fs':>6} {'wT':>6} {'a_R':>7} {'a_pm':>7} {'a_cart':>7} {'Q_end':>7} {'F':>6} {'F_lin':>6} {'g':>6} {'E':>7} {'Erel':>6} {'t1/2':>5} row")
    for r in sorted(modes, key=lambda r: r["w"]):
        hw = r["half_way_fs"]
        print(f"{r['j']:3d} {r['w']:7.0f} {r['period_fs']:6.0f} {r['omega_T']:6.1f} {r['a_passage']:7.3f} {r['a_pm']:7.3f} {r['a_pm_cartesian']:7.3f} {r['Q_end']:7.3f} {r['F']:6.2f} "
              f"{r['F_lin']:6.2f} {r['g']:6.2f} {r['energy_held_kcal']:7.2f} {r['relaxation_energy_kcal']:6.2f} {'' if hw is None else f'{hw:5.0f}':>5} {r['row']}")


if __name__ == "__main__":
    main(int(sys.argv[1]))
