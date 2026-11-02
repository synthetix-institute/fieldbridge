"""Which product a trajectory from a displaced start reaches, and how the turn moves.

    python displaced.py        # reads results/traj_E.npz, trajy_E_*.npz; writes results/displaced.json

The turn of the CHR end separates the products: the torsion is near 180 or 0 degrees in P1 and P2, and near +-90 at TS1. The reading is the cosine of the
torsion averaged over the last 40 fs, because the torsion swings widely inside a well.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from internal import Curvilinear  # noqa: E402
from modes import TORSION  # noqa: E402
from qm import BOND, OUT, align  # noqa: E402
from yamaki_geometries import atoms, coords  # noqa: E402

P1_COS = {31: -1.0, 32: 1.0}                      # cos(torsion) in P1
PHI = {31: (6, 0, 2, 3), 32: (8, 0, 2, 3)}        # the same torsion in both, through the carbon on the turning end: 180 in Z (P1), 0 in E (P2)


def main():
    out = {}
    for e in (31, 32):
        ts, irc = np.load(OUT / f"ts_{e}.npz"), np.load(OUT / f"irc_{e}.npz")
        C = Curvilinear([a[0] for a in atoms(e, "TS1")], [coords(e, n) for n in ("TS1", "PM", "P1", "P2")], BOND[e], irc["q"][0], ts["m"], ts["L"])
        n = C.nb + C.na + C.g[2].index(TORSION[e])
        nphi = C.nb + C.na + C.g[2].index(PHI[e])
        rows = []
        for f in [OUT / f"traj_{e}.npz"] + sorted(OUT.glob(f"trajy_{e}_*.npz")):
            d = np.load(f)
            th = np.array([C.xi(q)[n] for q in d["q"][::4]])
            phi = np.array([C.xi(q)[nphi] for q in d["q"][::4]])
            t = np.arange(len(th)) * 4 * float(d["dt_fs"])
            last = t >= t[-1] - 40
            side = float(np.mean(np.cos(th[last]))) * P1_COS[e]
            first = t[np.argmax(np.cos(th) * P1_COS[e] > 0.5)] if (np.cos(th) * P1_COS[e] > 0.5).any() else None
            other = t[np.argmax(np.cos(th) * P1_COS[e] < -0.5)] if (np.cos(th) * P1_COS[e] < -0.5).any() else None
            rows.append(dict(d=float(d["d"]) if "d" in d else 0.0, vy=float(d["vy"]) if "vy" in d else 0.0, t_end=float(t[-1]), side_P1=side,
                             reaches_P1_fs=None if first is None else float(first), reaches_P2_fs=None if other is None else float(other),
                             start_excess_kcal=float((d["E"][0] - d["e0"]) * 627.509) if "e0" in d else 0.0, phi_start=float(np.degrees(phi[0])), t=t.tolist(), phi=np.degrees(phi).tolist(),
                             torsion_deg={str(int(tt)): float(np.degrees(x)) for tt, x in zip(t, th) if int(tt) % 20 == 0 and tt == int(tt)}))
        out[str(e)] = sorted(rows, key=lambda r: (r["vy"] != 0, r["d"], r["vy"]))
        # the undisplaced trajectory along the fixed transverse direction
        tr, vrt = np.load(OUT / f"traj_{e}.npz"), np.load(OUT / f"vrt_{e}.npz")
        m, L, q0 = ts["m"], ts["L"], irc["q"][0]
        c = L[:, 1:3].T @ vrt["y"]
        c = c / np.linalg.norm(c)
        frames = np.array([align(f / np.sqrt(m), q0 / np.sqrt(m), m) * np.sqrt(m) for f in tr["q"][::4]])
        out[f"y_{e}"] = dict(t=(np.arange(len(frames)) * 4 * float(tr["dt_fs"])).tolist(), y=((frames - q0) @ (L[:, 1:3] @ c)).tolist())
        print("entry", e)
        for r in out[str(e)]:
            print(f"  d {r['d']:+.1f} vy {r['vy']:+.1f} end {r['t_end']:.0f} fs, mean cos towards P1 {r['side_P1']:+.2f}, P1 at {r['reaches_P1_fs']}, P2 at {r['reaches_P2_fs']}, start {r['start_excess_kcal']:.2f} kcal/mol",
                  [round(v) for v in list(r["torsion_deg"].values())[:11]])
    (OUT / "displaced.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
