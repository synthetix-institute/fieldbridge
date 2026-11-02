"""The linear law for the turn of the molecule, with its curvature and side force taken from exact Hessians at frames of the undisplaced trajectory.

    python law_molecule.py     # reads results/{ts,irc,vrt,traj}_E.npz, hessframe_E_*.npz, hess_E_t*.npy, displaced.json; writes results/law_molecule.json

    y'' = -k(t) y + f(t),   y(0) = y0, y'(0) = p0

k(t) is the curvature of the lowest direction across the velocity that is like the turn, f(t) the force along that direction, which is orthogonal to the
velocity and so free of the steep gradient along the path. The solution is y(t) = m(t) + G(t) y0 + H(t) p0. The start that separates the products is
y0* = -m / G, and the velocity counts against the position by H sigma_p / (G sigma_q). The frames past bond formation are beyond the linear stage.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from checks import reference_motions  # noqa: E402
from modes import KT, TIME_FS, W_CM  # noqa: E402
from qm import OUT, transverse_curvature  # noqa: E402
from yamaki_geometries import RATIO  # noqa: E402


def main():
    out = {}
    disp = json.loads((OUT / "displaced.json").read_text())
    for e in (31, 32):
        ts, tr, vrt = np.load(OUT / f"ts_{e}.npz"), np.load(OUT / f"traj_{e}.npz"), np.load(OUT / f"vrt_{e}.npz")
        m, L, w = ts["m"], ts["L"], np.abs(ts["w"]) / W_CM
        sm = np.sqrt(m)
        c = L[:, 1:3].T @ vrt["y"]
        c = c / np.linalg.norm(c)
        yhat = L[:, 1:3] @ c
        sigma_q = np.sqrt(KT * ((c / w[1:3]) ** 2).sum())
        sign0 = np.sign(reference_motions(tr["q"][0], m, e)["turn of the CHR end"] @ yhat)
        ts_, ks, fs, ov = [0.0], [float(((c * w[1:3]) ** 2).sum())], [0.0], [1.0]
        for t in (10, 15, 20, 25, 30, 40):
            i = int(round(t / float(tr["dt_fs"])))
            d = np.load(OUT / f"hessframe_{e}_{t}.npz")
            H = np.load(OUT / f"hess_{e}_t{t:g}.npy")
            Hm = 0.5 * (H + H.T) / np.outer(sm, sm)
            _, _, lam, V = transverse_curvature(Hm, tr["v"][i], tr["q"][i], m)
            u = reference_motions(tr["q"][i], m, e)["turn of the CHR end"]
            o = u @ V[:, :8]
            j = int(np.argmax(o ** 2))
            v = V[:, j] * np.sign(o[j]) * sign0            # positive along +y of the displaced starts
            ts_.append(t / TIME_FS)
            ks.append(float(lam[j]))
            fs.append(float(-(d["g"] @ v)))
            ov.append(float(o[j] ** 2))
        ts_, ks, fs = np.array(ts_), np.array(ks), np.array(fs)
        grid = np.linspace(0, ts_[-1], 4001)
        h = grid[1] - grid[0]
        kg, fg = np.interp(grid, ts_, ks), np.interp(grid, ts_, fs)
        def integrate(y, p, force):
            ys = [y]
            for i in range(len(grid) - 1):
                p += 0.5 * h * (-kg[i] * y + (fg[i] if force else 0.0))
                y += h * p
                p += 0.5 * h * (-kg[i + 1] * y + (fg[i + 1] if force else 0.0))
                ys.append(y)
            return np.array(ys)
        M, G, Hg = integrate(0.0, 0.0, True), integrate(1.0, 0.0, False), integrate(0.0, 1.0, False)
        curve = dict(t_fs=(grid[::40] * TIME_FS).tolist(), m_over_width=(M[::40] / sigma_q).tolist())
        res = {}
        for t in (30.0, 40.0):
            i = int(np.argmin(np.abs(grid - t / TIME_FS)))
            star = -M[i] / G[i]
            ratio = abs(Hg[i]) * np.sqrt(KT) / (abs(G[i]) * sigma_q)
            res[str(int(t))] = dict(m=float(M[i]), G=float(G[i]), H=float(Hg[i]), start_boundary_widths=float(star / sigma_q), velocity_over_position=float(ratio),
                                    z=float(abs(star) / np.sqrt(sigma_q ** 2 + (ratio * sigma_q) ** 2)), share=float(norm.cdf(abs(star) / np.sqrt(sigma_q ** 2 + (ratio * sigma_q) ** 2))))
        # measured: the two starts that bracket the boundary, and the torsion response at 60 fs
        rows = disp[str(e)]
        pos = sorted([r for r in rows if r["vy"] == 0], key=lambda r: r["d"])
        lo, hi = [(a["d"], b["d"]) for a, b in zip(pos, pos[1:]) if np.sign(a["side_P1"]) != np.sign(b["side_P1"])][0]
        tors = lambda r, t: r["torsion_deg"][str(t)]
        by = {(r["d"], r["vy"]): r for r in rows}
        resp = {}
        for t in (40, 60):
            dq = np.unwrap(np.radians([tors(by[(-1.0, 0.0)], t), tors(by[(1.0, 0.0)], t)]))
            dp = np.unwrap(np.radians([tors(by[(0.0, -1.0)], t), tors(by[(0.0, 1.0)], t)]))
            resp[str(t)] = float(abs((dp[1] - dp[0]) / (dq[1] - dq[0])))
        n1, n2 = RATIO[e]
        pub = n1 / (n1 + n2)
        out[str(e)] = dict(sigma_q=float(sigma_q), curve=curve, frames=dict(t_fs=(ts_ * TIME_FS).tolist(), k_cm=[float(np.sign(v) * np.sqrt(abs(v)) * W_CM) for v in ks], f=fs.tolist(), overlap=ov),
                           law=res, measured_boundary_widths=[lo, hi], measured_velocity_over_position=resp,
                           share_from_boundary=[float(norm.cdf(min(abs(lo), abs(hi)))), float(norm.cdf(max(abs(lo), abs(hi))))],
                           published=dict(counts=[n1, n2], share=pub, stderr=float(np.sqrt(pub * (1 - pub) / (n1 + n2))), z=float(norm.ppf(pub))))
        print(json.dumps({k: v for k, v in out[str(e)].items() if k != "curve"}, indent=0))
    (OUT / "law_molecule.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
