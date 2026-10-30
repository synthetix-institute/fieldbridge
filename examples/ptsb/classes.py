"""Static connections from TS1 to the product the IRC does not reach, and where that product's valley is born.

For each surface: the IRC from TS1 (which product it reaches); a double-ended string between TS1 and the other
product, classified as downhill (energy never rises after TS1) or small barrier (an interior maximum after a dip),
with the height of that barrier above the preceding minimum of the string and relative to TS1; and the points along
the path at which a transverse minimum is born (`folds`): on a symmetric surface the pitchfork at the VRI, on an
asymmetric one the fold at which the minor valley appears, with its energy relative to TS1.

    python classes.py --out results/classes.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ptsb import Collins, Minimal, folds, irc, string_path, ts_frame  # noqa: E402


def classify(s, images=81):
    p = s.points()
    path = irc(s)
    end = path[-1]
    major = min(("P+", "P-"), key=lambda k: np.inf if p[k] is None else np.linalg.norm(p[k] - end))
    near = np.linalg.norm(p[major] - end) < 1e-2
    other = "P-" if major == "P+" else "P+"
    q0 = p["TS1"]
    out = {"surface": repr(s), "IRC_end": end.tolist(), "IRC_product": major if near else "TS2 (ridge)",
           "V_TS1": float(s.V(*q0)), "V_TS2": float(s.V(*p["TS2"]))}
    targets = [other] if near else ["P+", "P-"]
    out["connections"] = []
    for tgt in targets:
        if p[tgt] is None:
            out["connections"].append({"to": tgt, "exists": False})
            continue
        # start from the IRC up to the VRI region then straight to the target, so that the string starts on the
        # reactive side of TS1 rather than across the reactant region
        k = int(np.argmin(np.abs(path[:, 0] - s.vri())))
        start = np.r_[np.linspace(q0, path[k], images // 3, endpoint=False),
                      np.linspace(path[k], p[tgt], images - images // 3)]
        P, E, it = string_path(s, q0, p[tgt], start=start, iters=60000, step=2e-3)
        E = np.asarray(E) - out["V_TS1"]
        rises = np.where(np.diff(E) > 1e-7)[0]
        if rises.size == 0:
            cls, barrier, above = "downhill", 0.0, None
        else:
            i_top = rises[0] + 1 + int(np.argmax(E[rises[0] + 1:]))
            i_dip = int(np.argmin(E[:i_top + 1]))
            cls = "small barrier"
            barrier, above = float(E[i_top] - E[i_dip]), float(E[i_top])
            # the interior maximum of a converged string is a saddle; name it
            near_ts2 = np.linalg.norm(P[i_top] - p["TS2"]) < 0.05
            out_saddle = "TS2" if near_ts2 else f"({P[i_top][0]:.3f}, {P[i_top][1]:.3f})"
        out["connections"].append({"to": tgt, "exists": True, "class": cls, "barrier_above_dip": barrier,
                                   "top_relative_to_TS1": above,
                                   "saddle": (out_saddle if cls == "small barrier" else None),
                                   "iterations": it, "path": P[::4].tolist(), "energy": E[::4].tolist()})
    fl = folds(s, 0.0, p["TS2"][0] + 0.3)
    out["births"] = [{"x": f[0], "y": f[1], "V_minus_TS1": f[2] - out["V_TS1"]} for f in fl]
    out["x_VRI"] = s.vri()
    out["TS1"] = q0.tolist()
    out["omega_TS1"] = ts_frame(s)[3]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    surfaces = [Minimal(e, 0.0) for e in (0.0, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.24)] + \
               [Minimal(0.0, c) for c in (0.1, 0.2, 0.4)] + [Collins(x) for x in (2.0, 2.05, 2.1)]
    rows = []
    for s in surfaces:
        r = classify(s)
        rows.append(r)
        con = "; ".join(f"{c['to']}: {c.get('class')} (barrier {c.get('barrier_above_dip', 0):.4f}, "
                        f"top {c.get('top_relative_to_TS1')}, at {c.get('saddle')})" for c in r["connections"])
        births = ", ".join(f"x={b['x']:.3f} y={b['y']:+.3f} V-V_TS1={b['V_minus_TS1']:+.3f}" for b in r["births"])
        print(f"{r['surface']}: IRC -> {r['IRC_product']}; {con}; births {births}; VRI {r['x_VRI']:.3f}", flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
