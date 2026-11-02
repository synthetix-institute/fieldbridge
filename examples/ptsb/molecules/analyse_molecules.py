"""Figures for the two radical additions to methylallene.

    python analyse_molecules.py    # reads results/{modes,checks,curvature,displaced,law_molecule}*.json; writes ../figures/molecule_modes.png and molecule_turn.png
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from analyse import CAT, FIG, INK, INK2  # noqa: E402

OUT = HERE / "results"
NAME = {31: "CH3 radical", 32: "CH2F radical"}


def J(name):
    return json.loads((OUT / name).read_text())


def main():
    M = {e: J(f"modes_{e}.json") for e in (31, 32)}
    K = {e: J(f"checks_{e}.json") for e in (31, 32)}
    D, L = J("displaced.json"), J("law_molecule.json")

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))
    ax = axes[0]
    for i, e in enumerate((31, 32)):
        p = [l for l in K[e]["lowest"] if l["where"] == "path" and l["turn_overlap"] >= 0.7]
        f = [l for l in K[e]["lowest"] if l["where"] == "trajectory" and l["turn_overlap"] >= 0.7 and l["t"] <= 40]
        ax.plot([l["R"] for l in p], [l["turn_cm"] for l in p], "-", lw=1.6, color=CAT[i], label=f"{NAME[e]}, path")
        ax.plot([l["R"] for l in f], [l["turn_cm"] for l in f], "o", ms=6, color=CAT[i], mec="white", mew=0.8, label=f"{NAME[e]}, trajectory")
    ax.axhline(0, color=INK2, lw=0.8)
    ax.invert_xaxis()
    ax.set_xlabel("forming C-C bond (angstrom)")
    ax.set_ylabel("curvature of the turn (cm-1, negative imaginary)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.45), fontsize=7)
    ax = axes[1]
    for i, e in enumerate((31, 32)):
        r = [x for x in K[e]["modes"] if x["row"] != "uninvolved"]
        o = np.argsort([x["w"] for x in r])
        w = np.array([r[j]["w"] for j in o])
        ax.plot(w, [r[j]["F_cart_lin"] for j in o], "x", ms=4.5, color=CAT[i], mew=1.1, label=f"{NAME[e]}, harmonic mode")
        ax.plot(w, [r[j]["F_cart"] for j in o], "o", ms=5, color=CAT[i], mec="white", mew=0.8, label=f"{NAME[e]}, trajectory")
    ax.set_xscale("log")
    ax.set_xlabel("frequency of the mode at TS1 (cm-1)")
    ax.set_ylabel("share of the path displacement reached")
    ax.legend(loc="upper left", fontsize=7)
    ax = axes[2]
    for i, e in enumerate((31, 32)):
        r = sorted(M[e]["modes"], key=lambda x: x["w"])
        r = [x for x in r if x["g"] >= 0.03 and abs(x["a_pm"]) >= 0.03]
        ax.plot([x["w"] for x in r], [x["g"] for x in r], "o", ms=5, color=CAT[i], mec="white", mew=0.8, label=NAME[e])
    ax.axhline(1, color=INK2, lw=0.8, ls="--")
    ax.axvspan(50, 240, color=INK2, alpha=0.08, lw=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("frequency of the mode at TS1 (cm-1)")
    ax.set_ylabel("coupling g in the units of the model")
    ax.legend(loc="upper left", fontsize=8)
    fig.suptitle("Left: curvature of the turn while it stays one direction. Centre: how far each mode gets by bond formation. "
                 "Right: dashed g = 1, shaded the frequencies of the model", color=INK, fontsize=8.5)
    fig.subplots_adjust(wspace=0.32)
    fig.savefig(FIG / "molecule_modes.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))
    for ax, e in zip(axes[:2], (31, 32)):
        for r in D[str(e)]:
            if r["vy"] != 0:
                continue
            phi = np.abs(np.array(r["phi"]))
            z = r["side_P1"] > 0
            ax.plot(r["t"], phi, "-", lw=1.5 if r["d"] != 0 else 2.2, color=CAT[0] if z else CAT[1])
        ax.axhline(87, color=INK2, lw=0.8, ls="--")
        ax.set_xlim(0, 200)
        ax.set_ylim(0, 180)
        ax.set_yticks([0, 45, 90, 135, 180])
        ax.set_xlabel("time from TS1 (fs)")
        ax.set_title(NAME[e], color=INK, fontsize=9)
    axes[0].set_ylabel("torsion of the turning end (degrees): 180 Z, 0 E")
    axes[0].plot([], [], "-", color=CAT[0], label="ends in Z")
    axes[0].plot([], [], "-", color=CAT[1], label="ends in E")
    axes[0].legend(loc="lower right", fontsize=8)
    ax = axes[2]
    for i, e in enumerate((31, 32)):
        s = L[str(e)]["sigma_q"]
        y = D[f"y_{e}"]
        ax.plot(y["t"], np.array(y["y"]) / s, "-", lw=1.6, color=CAT[i], label=f"{NAME[e]}, trajectory")
        c = L[str(e)]["curve"]
        ax.plot(c["t_fs"], c["m_over_width"], ":", lw=2.0, color=CAT[i], label=f"{NAME[e]}, linear law")
    ax.axvline(30, color=INK2, lw=0.8, ls="--")
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xlim(0, 100)
    ax.set_xlabel("time from TS1 (fs); dashed: bond formed")
    ax.set_ylabel("displacement along the turn (thermal widths)")
    ax.legend(loc="lower left", fontsize=7)
    fig.suptitle("Left and centre: trajectories started at different torsions in the thermal range of TS1; dashed, the perpendicular geometry. "
                 "Right: the undisplaced trajectory", color=INK, fontsize=8.5)
    fig.subplots_adjust(wspace=0.32, top=0.86)
    fig.savefig(FIG / "molecule_turn.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
