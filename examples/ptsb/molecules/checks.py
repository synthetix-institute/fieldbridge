"""What the transverse direction is, how much of the following is geometry, and how far the instability grows.

    python checks.py 31        # reads results/{ts,irc,vrt,traj,modes,curvature}_31.*; writes results/checks_31.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from internal import Curvilinear  # noqa: E402
from modes import NEG, TIME_FS, W_CM, linear_response  # noqa: E402
from qm import BOND, OUT, align, project  # noqa: E402
from yamaki_geometries import atoms, coords  # noqa: E402

GROUPS = {31: dict(end=(0, 2, [1, 6, 7, 8, 9]), ch2=(3, 2, [4, 5]), radical=(10, [11, 12, 13])),
          32: dict(end=(0, 2, [1, 8, 9, 10, 11]), ch2=(3, 2, [4, 5]), radical=(6, [7, 12, 13]))}


def rigid_rotation(q, m, atom_set, centre, axis):
    sm = np.sqrt(m)
    x = (q / sm).reshape(-1, 3)
    v = np.zeros_like(x)
    for i in atom_set:
        v[i] = np.cross(axis, x[i] - centre)
    u = (v.ravel()) * sm
    u = project(q, m) @ u
    return u / np.linalg.norm(u)


def reference_motions(q, m, entry):
    """Named motions at q: the turn of each end of the allene about its bond to the centre, and the tilt and sideways shift of the radical."""
    g = GROUPS[entry]
    x = (q / np.sqrt(m)).reshape(-1, 3)
    out = {}
    a, b, s = g["end"]
    out["turn of the CHR end"] = rigid_rotation(q, m, s, x[a], (x[a] - x[b]) / np.linalg.norm(x[a] - x[b]))
    a, b, s = g["ch2"]
    out["turn of the CH2 end"] = rigid_rotation(q, m, s, x[a], (x[a] - x[b]) / np.linalg.norm(x[a] - x[b]))
    c, h = g["radical"]
    rad = [c] + h
    bond = x[2] - x[c]
    bond = bond / np.linalg.norm(bond)
    ref = np.array([1.0, 0, 0]) if abs(bond[0]) < 0.8 else np.array([0, 1.0, 0])
    e1 = np.cross(bond, ref)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(bond, e1)
    for name, ax in (("tilt of the radical, axis 1", e1), ("tilt of the radical, axis 2", e2), ("spin of the radical", bond)):
        out[name] = rigid_rotation(q, m, rad, x[c], ax)
    for name, ax in (("shift of the radical, axis 1", e1), ("shift of the radical, axis 2", e2)):
        v = np.zeros_like(x)
        v[rad] = ax
        u = project(q, m) @ (v.ravel() * np.sqrt(m))
        out[name] = u / np.linalg.norm(u)
    return out


def turn_branch(lam, V, u):
    """Among the lowest directions across the motion, the one most like the turn u: its signed frequency and its overlap."""
    ov = (u @ V[:, :8]) ** 2
    i = int(np.argmax(ov))
    return float(np.sign(lam[i]) * np.sqrt(abs(lam[i])) * W_CM), float(ov[i])


def main(entry):
    ts = np.load(OUT / f"ts_{entry}.npz")
    irc = np.load(OUT / f"irc_{entry}.npz")
    tr = np.load(OUT / f"traj_{entry}.npz")
    mo = np.load(OUT / f"modes_{entry}.npz")
    vrt = np.load(OUT / f"vrt_{entry}.npz")
    cu = json.loads((OUT / f"curvature_{entry}.json").read_text())
    m, L, wcm = ts["m"], ts["L"], ts["w"]
    sm = np.sqrt(m)
    q = irc["q"]
    q0 = q[0]
    C = Curvilinear([a[0] for a in atoms(entry, "TS1")], [coords(entry, n) for n in ("TS1", "PM", "P1", "P2")], BOND[entry], q0, m, L)
    end, near, kR = int(mo["end"]), mo["near"], int(mo["kR"])
    a, Q = mo["a"], mo["Q"]

    # what y is
    k = int(vrt["k"])
    refs = reference_motions(q[k], m, entry)
    shares = {n: float((vrt["y"] @ u) ** 2) for n, u in refs.items()}
    refs0 = reference_motions(q0, m, entry)
    mode_names = {}
    for j in range(1, 13):
        mode_names[f"{wcm[j]:.0f}"] = {n: round(float((L[:, j] @ u) ** 2), 2) for n, u in refs0.items() if (L[:, j] @ u) ** 2 > 0.1}

    # how much of the following is geometry: a straight line along the imaginary mode, read in the same coordinates
    frames = np.array([align(f / sm, q0 / sm, m) * sm for f in tr["q"][:end + 1]])
    rect = (frames - q0) @ L
    Qkin = np.array([C(q0 + L[:, 0] * x) for x in rect[:, 0]])
    a_rect = (q - q0) @ L
    rows = []
    for j in range(1, L.shape[1]):
        aj = a[near, j]
        den = float(aj @ aj)
        ar = a_rect[near, j]
        Ql = linear_response(ar, mo["tt"][: end + 1] / TIME_FS, abs(wcm[j]) / W_CM)
        rows.append(dict(j=j, w=float(wcm[j]), a_R=float(a[kR, j]), a_R_cart=float(a_rect[kR, j]), F=float(Q[: end + 1, j] @ aj / den), F_kin=float(Qkin[:, j] @ aj / den),
                         F_cart=float(rect[:, j] @ ar / (ar @ ar)), F_cart_lin=float(Ql @ ar / (ar @ ar)), Q_end=float(Q[end, j]), Q_kin_end=float(Qkin[-1, j]),
                         Q_cart_end=float(rect[-1, j])))

    # energy above the path at the same bond length
    dV = mo["dV"]
    dV_at = {str(t): float(dV[int(t * 2)]) for t in (10, 15, 20, 25) if t * 2 <= end}
    dV_at["end"] = float(dV[-1])
    harm_at = {str(t): float(mo["harm"][int(t * 2)]) for t in (10, 15, 20, 25) if t * 2 <= end}

    # growth of the instability, the integral of sqrt(-curvature) dt over the frames with lowest curvature negative
    fr = cu["frames"]
    tf = np.array([f["t"] for f in fr]) / TIME_FS
    lam = np.array([np.sign(f["lowest_cm"][0]) * (f["lowest_cm"][0] / W_CM) ** 2 for f in fr])
    tg = np.linspace(tf[0], tf[-1], 2001)
    lg = np.interp(tg, tf, lam)
    rate = np.sqrt(np.clip(-lg, 0, None))
    growth = {}
    for t_fs in (30, 40, 50, 60):
        mask = tg <= t_fs / TIME_FS
        growth[str(t_fs)] = float(np.sum(0.5 * (rate[mask][1:] + rate[mask][:-1]) * np.diff(tg[mask])))
    # what the lowest direction is, along the path and on the trajectory
    lowest = []
    from qm import transverse_curvature
    for p in cu["path"]:
        kk = p["k"]
        f = OUT / ("hess_ts_%d.npy" % entry if kk == 0 else "hess_vrt_%d.npy" % entry if kk == k else f"hess_{entry}_{kk}.npy")
        H = np.load(f)
        Hm = 0.5 * (H + H.T) / np.outer(sm, sm)
        t = q[min(kk + 1, len(q) - 1)] - q[max(kk - 1, 0)]
        _, v0, lam_, V = transverse_curvature(Hm, t, q[kk], m)
        ref = reference_motions(q[kk], m, entry)
        tb = turn_branch(lam_, V, ref["turn of the CHR end"])
        lowest.append(dict(where="path", s=p["s"], R=p["R"], cm=p["lowest_cm"][:2], turn_cm=tb[0], turn_overlap=tb[1], first={n: round(float((v0 @ u) ** 2), 2) for n, u in ref.items()},
                           second={n: round(float((V[:, 1] @ u) ** 2), 2) for n, u in ref.items()}))
    for fr_ in cu["frames"]:
        i = int(round(fr_["t"] / float(tr["dt_fs"])))
        H = np.load(OUT / f"hess_{entry}_t{fr_['t']:g}.npy")
        Hm = 0.5 * (H + H.T) / np.outer(sm, sm)
        _, v0, lam_, V = transverse_curvature(Hm, tr["v"][i], tr["q"][i], m)
        ref = reference_motions(tr["q"][i], m, entry)
        tb = turn_branch(lam_, V, ref["turn of the CHR end"])
        gf = np.load(OUT / f"hessframe_{entry}_{fr_['t']:g}.npz")["g"]
        lowest.append(dict(where="trajectory", t=fr_["t"], R=fr_["R"], cm=fr_["lowest_cm"][:2], turn_cm=tb[0], turn_overlap=tb[1],
                           turn_force=float(-(gf @ ref["turn of the CHR end"])), first={n: round(float((v0 @ u) ** 2), 2) for n, u in ref.items()},
                           second={n: round(float((V[:, 1] @ u) ** 2), 2) for n, u in ref.items()}))
    mov = [r for r in rows if abs(r["a_R_cart"]) >= NEG]
    tf_ = np.array([l["t"] for l in lowest if l["where"] == "trajectory"]) / TIME_FS
    lam_t = np.array([np.sign(l["turn_cm"]) * (l["turn_cm"] / W_CM) ** 2 for l in lowest if l["where"] == "trajectory"])
    tg_ = np.linspace(tf_[0], tf_[-1], 2001)
    rate_ = np.sqrt(np.clip(-np.interp(tg_, tf_, lam_t), 0, None))
    growth_turn = {str(t_fs): float(np.sum(0.5 * (rate_[tg_ <= t_fs / TIME_FS][1:] + rate_[tg_ <= t_fs / TIME_FS][:-1]) * np.diff(tg_[tg_ <= t_fs / TIME_FS]))) for t_fs in (30, 40)}
    row_of = lambda r: "stays behind" if r["F_cart"] < 0.3 else ("between" if r["F_cart"] <= 0.7 else "follows")
    for r in rows:
        r["row"] = row_of(r) if abs(r["a_R_cart"]) >= NEG else "uninvolved"
    rows_by = {k: [round(r["w"]) for r in rows if r["row"] == k] for k in ("stays behind", "between", "follows")}
    res = dict(entry=entry, growth_turn=growth_turn, rows_by=rows_by, lowest=lowest, F_cart_vs_lin_mean_abs=float(np.mean([abs(r["F_cart"] - r["F_cart_lin"]) for r in mov])),
               F_cart_vs_lin_corr=float(np.corrcoef([r["F_cart"] for r in mov], [r["F_cart_lin"] for r in mov])[0, 1]), y_shares=shares, mode_names=mode_names, modes=rows, dV_above_path_kcal=dV_at, harmonic_kcal=harm_at, growth_exponent=growth)
    (OUT / f"checks_{entry}.json").write_text(json.dumps(res, indent=1))
    print("y at the VRT:", {n: round(v, 2) for n, v in shares.items() if v > 0.03})
    print("modes:", json.dumps(mode_names))
    print("dV", dV_at, "harm", harm_at, "growth", growth)
    print("mean abs F_cart - F_cart_lin", res["F_cart_vs_lin_mean_abs"], "corr", res["F_cart_vs_lin_corr"], "modes", len(mov))
    print(" cm-1 a_cart F_cart F_lin | a_int F_int F_kin")
    for r in mov:
        print(f"{r['w']:6.0f} {r['a_R_cart']:6.3f} {r['F_cart']:6.2f} {r['F_cart_lin']:6.2f} | {r['a_R']:6.3f} {r['F']:6.2f} {r['F_kin']:6.2f}")
    print("rows", rows_by, "growth of the turn", growth_turn)
    for l in lowest:
        print(l["where"], round(l.get("s", l.get("t")), 2), round(l["R"], 2), np.round(l["cm"], 0), "turn", round(l["turn_cm"]), round(l["turn_overlap"], 2), "force", round(l.get("turn_force", 0), 4), {n[:22]: v for n, v in l["first"].items() if v >= 0.2}, "")


if __name__ == "__main__":
    main(int(sys.argv[1]))
