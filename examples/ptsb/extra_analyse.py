"""Figures and summary numbers for the bath modes and the two transverse modes.

    python extra_analyse.py      # reads results/bath_*.json, two_modes_*.json; writes figures/bath_modes.png,
                                 # figures/two_modes.png and results/extra_summary.json
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from analyse import CAT, FIG, HERE, INK, INK2, SEQ  # noqa: E402,F401


def J(name):
    return json.loads((HERE / "results" / name).read_text())


def main():
    B = {p["id"]: p for p in J("bath_traj.json")["points"]}
    BL = J("bath_law.json")["laws"]
    TM = J("two_modes_traj.json")["points"]
    TL = J("two_modes_law.json")["laws"]
    TS = J("two_modes_sep.json")["laws"]
    S = {}

    # ---- bath modes
    Ms = (1, 2, 4, 8, 16, 32)
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
    for ax, kT in zip(axes[:2], (0.01, 0.05)):
        for i, g in enumerate((0.3, 1.0, 3.0)):
            P = [1 - B[f"Bk_g{g}_M{M}_T{kT}"]["outcome"]["P_plus"] for M in Ms]
            se = [B[f"Bk_g{g}_M{M}_T{kT}"]["outcome"]["stderr"] for M in Ms]
            ax.errorbar(Ms, P, yerr=se, fmt="o-", ms=4.5, lw=1.4, color=CAT[i], mec="white", mew=0.8,
                        label=f"gamma = {g:g}")
            ax.axhline(1 - B[f"Bk_g{g}_M0_T{kT}"]["outcome"]["P_plus"], color=CAT[i], lw=1, ls="--")
        ax.set_xscale("log", base=2)
        ax.set_xticks(Ms)
        ax.set_xticklabels([str(m) for m in Ms])
        ax.set_xlabel("modes per coordinate (band up to 4)")
        ax.set_title(f"curved valley, curvature 0.1, kT = {kT:g}", color=INK, fontsize=9)
    axes[0].set_ylabel("share of the product outside the bend")
    axes[0].legend(loc="upper right", fontsize=8)
    fig.subplots_adjust(wspace=0.3)
    ax = axes[2]
    Mf = (1, 4, 16, 32)
    for i, g in enumerate((0.3, 1.0, 3.0)):
        P = [B[f"Bf_g{g}_M{M}_v3.0"]["outcome"]["P_plus"] for M in Mf]
        ax.plot(Mf, P, "o-", ms=4.5, lw=1.4, color=CAT[i], mec="white", mew=0.8)
        ax.axhline(B[f"Bf_g{g}_M0_v3.0"]["outcome"]["P_plus"], color=CAT[i], lw=1, ls="--")
        ax.plot(Mf, [BL[f"Bf_g{g}_M{M}_v3.0"]["law"]["P_law"] for M in Mf], ":", color=CAT[i], lw=1.4)
    ax.set_xscale("log", base=2)
    ax.set_xticks(Mf)
    ax.set_xticklabels([str(m) for m in Mf])
    ax.set_xlabel("number of modes")
    ax.set_ylabel("share of the minor product")
    ax.set_title("straight valley, tilt 0.04, speed 3, kT = 0.005", color=INK, fontsize=9)
    fig.suptitle("Harmonic modes in place of friction. Points: trajectories with modes; dashed: Langevin trajectories "
                 "at the same gamma; dotted: exact linear law with the modes", color=INK, fontsize=9.5)
    fig.savefig(FIG / "bath_modes.png")
    plt.close(fig)

    S["bath_curved_mean_abs_vs_langevin"] = {
        str(M): float(np.mean([abs(B[f"Bk_g{g}_M{M}_T{kT}"]["outcome"]["P_plus"]
                                   - B[f"Bk_g{g}_M0_T{kT}"]["outcome"]["P_plus"])
                               for g in (0.3, 1.0, 3.0) for kT in (0.01, 0.05)])) for M in Ms}
    S["bath_fast_mean_abs_vs_langevin"] = {
        str(M): float(np.mean([abs(B[f"Bf_g{g}_M{M}_v{v}"]["outcome"]["P_plus"]
                                   - B[f"Bf_g{g}_M0_v{v}"]["outcome"]["P_plus"])
                               for g in (0.3, 1.0, 3.0) for v in (1.0, 3.0)])) for M in Mf}
    S["bath_thermal_tilt_mean_abs_vs_law"] = {
        str(M): float(np.mean([abs(B[f"B_g{g}_M{M}_T{kT}"]["outcome"]["P_plus"]
                                   - BL[f"B_g{g}_M{M}_T{kT}"]["law"]["P_law"])
                               for g in (0.3, 1.0, 3.0) for kT in (0.01, 0.05)])) for M in Ms}
    S["bath_band"] = {f"g{g}_T{kT}": {"Langevin": 1 - B[f"Bk_g{g}_M0_T{kT}"]["outcome"]["P_plus"],
                                      "band1": 1 - B[f"Bw_g{g}_W1.0_T{kT}"]["outcome"]["P_plus"],
                                      "band4": 1 - B[f"Bk_g{g}_M32_T{kT}"]["outcome"]["P_plus"],
                                      "band16": 1 - B[f"Bw_g{g}_W16.0_T{kT}"]["outcome"]["P_plus"]}
                      for g in (0.3, 1.0, 3.0) for kT in (0.01, 0.05)}

    # ---- two modes
    T = [p for p in TM if p["family"] == "T"]
    W = [p for p in TM if p["family"] == "W"]
    S["quadrants"] = {}
    for g in (0.0, 1.0, 3.0):
        rows = [p for p in T if p["gamma"] == g]
        S["quadrants"][str(g)] = {
            "sign_mean_abs": float(np.mean([abs(p["result"][k] - TL[p["id"]][k]) for p in rows
                                            for k in ("P_y_plus", "P_z_plus")])),
            "sign_max_abs": float(np.max([abs(p["result"][k] - TL[p["id"]][k]) for p in rows
                                          for k in ("P_y_plus", "P_z_plus")])),
            "joint_vs_product_of_marginals": float(np.mean([abs(p["result"]["P_both_plus"]
                                                               - p["result"]["P_y_plus"] * p["result"]["P_z_plus"])
                                                           for p in rows]))}
    S["axes"] = {}
    for g in (0.0, 1.0, 3.0):
        for ez in (0.0, 0.05, 0.2):
            rows = [p for p in W if p["gamma"] == g and p["surface"]["eps_z"] == ez]
            S["axes"][f"g{g}_ez{ez}"] = {
                "simple_mean_abs": float(np.mean([abs(p["result"]["P_y_family"] - TL[p["id"]]["P_y_family"])
                                                  for p in rows])),
                "boundary_mean_abs": float(np.nanmean([abs(p["result"]["P_y_family"] - TS[p["id"]]["P_y_family"])
                                                       for p in rows]))}

    Ds = (0.0, 0.02, 0.05, 0.1, 0.2)
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    for ax, g in zip(axes, (0.0, 1.0, 3.0)):
        for i, kT in enumerate((0.01, 0.05)):
            rows = sorted([p for p in W if p["gamma"] == g and p["kT"] == kT and p["surface"]["eps_z"] == 0.0],
                          key=lambda p: p["surface"]["D"])
            ax.plot(Ds, [p["result"]["P_y_family"] for p in rows], "o", ms=5, color=CAT[i], mec="white", mew=0.8,
                    label=f"kT = {kT:g}")
            ax.plot(Ds, [TL[p["id"]]["P_y_family"] for p in rows], ":", lw=1.4, color=CAT[i])
            ax.plot(Ds, [TS[p["id"]]["P_y_family"] for p in rows], "-", lw=1.8, color=CAT[i])
        ax.axhline(0.5, color=INK2, lw=0.8, ls=":")
        ax.set_xlabel("delay D of the second instability")
        ax.set_title(f"friction gamma = {g:g}", color=INK, fontsize=9)
    axes[0].set_ylabel("share of the two products of the earlier mode")
    axes[0].legend(loc="lower right", fontsize=8)
    fig.suptitle("Two modes with products on the axes. Points: trajectories; dotted: comparison |y| > |z| in the linear "
                 "stage; solid: with the boundary found on the surface", color=INK, fontsize=9.5)
    fig.savefig(FIG / "two_modes.png")
    plt.close(fig)

    (HERE / "results/extra_summary.json").write_text(json.dumps(S, indent=1))
    print(json.dumps(S, indent=0)[:3500])


if __name__ == "__main__":
    main()
