"""Figures and summary numbers for the cubic coupling to further modes.

    python anharm_analyse.py      # reads results/anharm_traj.json, results/anharm_law.json; writes figures/anharm_*.png and
                                  # results/anharm_summary.json
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from analyse import CAT, FIG, HERE, INK, INK2  # noqa: E402
from anharm import frequencies  # noqa: E402

SEEDS = (1, 2, 3)


def share(x):
    return x["outcome"]["P_plus"]


def returns(x):
    return x["outcome"]["reactant"] / x["outcome"]["n"]


def friction(x):
    """Energy lost by (x, y) against the Langevin balance. Includes the reversible energy of following."""
    if x.get("dE_surface") is None:
        return float("nan")
    return -x["dE_surface"] / (x["int_p2"] - 2 * x["kT"] * x["t_mean"])


def boltzmann_entropic(g, M=16, y=0.88):
    """Product share of the entropic family if the products were populated in proportion to their free energies."""
    w = frequencies(M)
    r = np.exp(-0.5 * np.log((w ** 2 + g * w * y) / (w ** 2 - g * w * y)).sum())
    return float(r / (1 + r))


def main():
    B = {p["id"]: p for p in json.loads((HERE / "results/anharm_traj.json").read_text())["points"]}
    L = json.loads((HERE / "results/anharm_law.json").read_text())["laws"]
    S = {}

    def F(cond, st, g, fam="F"):
        return [B[f"{fam}_{cond}_{st}_g{g}_s{s}"] for s in SEEDS if f"{fam}_{cond}_{st}_g{g}_s{s}" in B
                and "outcome" in B[f"{fam}_{cond}_{st}_g{g}_s{s}"]]

    conds = ["curved_T0.01", "curved_T0.05", "fast_T0.005", "tilt_T0.01", "tilt_T0.05"]

    # references
    S["references_friction_from_energy"] = {k: friction(B[k]) for k in
                                            ("R_curved_T0.01_g0.3", "R_curved_T0.01_g1.0", "R_curved_T0.01_g3.0",
                                             "R_tilt_T0.05_g0.3", "R_tilt_T0.05_g1.0", "R_tilt_T0.05_g3.0")}
    S["bare"] = {c: {"share": share(B[f"R_{c}"]), "returns": returns(B[f"R_{c}"])} for c in conds}

    # quadratic forces
    quad = [p for p in B.values() if p["family"] == "F" and p["set"] in ("xQQ", "yQQ", "quad") and "outcome" in p]
    S["quadratic_largest_shift_of_share"] = max(abs(share(p) - share(B[f"R_{p['name']}_T{p['kT']}"])) for p in quad)
    S["quadratic_share_stderr"] = float(np.mean([p["outcome"]["stderr"] for p in quad]))
    S["quadratic_unbound_tables"] = sum(1 for p in B.values() if p.get("unbound"))
    S["quadratic_friction_vs_equivalent"] = {
        f"{st}_g{g}_T{kT}": {"cubic": float(np.mean([friction(x) for x in F(f"curved_T{kT}", st, g)])),
                             "equivalent": float(np.mean([friction(x) for x in F(f"curved_T{kT}", st, g, "E")]))}
        for st, g in (("xQQ", 1.0), ("quad", 0.5), ("yQQ", 1.0)) for kT in (0.01, 0.05)}
    S["K0_at_kT0.01"] = {st: float(np.mean([p["scale"][f"K0_{st}"] for p in B.values() if p["family"] == "F"
                                            and p.get("set") == st and p["g"] == g and p["kT"] == 0.01
                                            and "scale" in p]))
                         for st, g in (("xQQ", 0.5), ("yQQ", 0.5), ("xxQ", 1.0), ("xyQ", 1.0), ("yyQ", 1.0))}

    # displacement forces
    S["displacement_friction"] = {f"{st}_g{g}": {f"T{kT}": float(np.mean([friction(x) for x in F(f"curved_T{kT}", st, g)]))
                                                 for kT in (0.01, 0.05)}
                                  for st in ("xxQ", "xyQ", "yyQ") for g in (0.3, 1.0, 3.0)}
    S["table_scatter_all_g1"] = {c: [share(x) for x in F(c, "all", 1.0)] for c in conds}
    S["table_scatter_M"] = {f"curved_T0.01_M{M}": [share(B[f"M_curved_T0.01_M{M}_s{s}"]) for s in SEEDS] for M in (4, 8, 32)}
    gs = [0, 0.3, 1.0, 3.0]
    dev = {}
    for kT in (0.01, 0.05):
        Ps = [share(B[f"R_curved_T{kT}"])] + [share(B[f"R_curved_T{kT}_g{g}"]) for g in gs[1:]]
        rows = [x for st in ("xxQ", "xyQ", "yyQ", "disp", "all") for g in (0.3, 1.0, 3.0)
                for x in F(f"curved_T{kT}", st, g) if np.isfinite(share(x))]
        d = [share(x) - np.interp(max(friction(x), 0), gs, Ps) for x in rows]
        dev[f"T{kT}"] = {"mean_abs": float(np.mean(np.abs(d))), "max_abs": float(np.max(np.abs(d)))}
    S["share_minus_langevin_at_the_friction_read_from_energy"] = dev
    S["occupations"] = {f"{n}_{o}": [share(B[f"Q_{n}_{o}_s{s}"]) for s in SEEDS]
                        for n in ("curved", "tilt") for o in ("classical", "bose", "wigner")}
    S["mode_mode_coupling"] = {f"g{g}": {"without": [share(x) for x in F("curved_T0.01", "all", g)],
                                         "with": [share(x) for x in F("curved_T0.01", "allQ3", g)],
                                         "returns_with": [returns(x) for x in F("curved_T0.01", "allQ3", g)],
                                         "modes_escaped": [x["bad"] for x in F("curved_T0.01", "allQ3", g)]}
                               for g in (0.3, 1.0, 3.0)}

    # the limits
    lim = {}
    for kT in (0.01, 0.05):
        for st in ("xxQ", "xyQ", "disp"):
            lo = [(share(B[f"S_T{kT}_{st}_low_s{s}"]), share(B[f"Z_T{kT}_{st}_s{s}"])) for s in SEEDS]
            hi = [(share(B[f"S_T{kT}_{st}_high_s{s}"]), share(B[f"ZG_T{kT}_{st}_s{s}"])) for s in SEEDS]
            lim[f"{st}_T{kT}"] = {"slow_full_vs_frozen": lo, "fast_full_vs_following": hi,
                                  "lam8_vs_following": [(share(B[f"S_T{kT}_{st}_lam8.0_s{s}"]),
                                                         share(B[f"G_curved_T{kT}_{st}_g1.0_s{s}"])) for s in SEEDS]}
    S["limits"] = lim
    S["limits_mean_abs"] = {k: float(np.mean([abs(a - b) for st in ("xxQ", "xyQ", "disp") for kT in (0.01, 0.05)
                                              for a, b in lim[f"{st}_T{kT}"][k]]))
                            for k in ("slow_full_vs_frozen", "fast_full_vs_following", "lam8_vs_following")}
    S["limits_mean_abs_xyQ"] = {k: float(np.mean([abs(a - b) for kT in (0.01, 0.05) for a, b in lim[f"xyQ_T{kT}"][k]]))
                                for k in ("slow_full_vs_frozen", "fast_full_vs_following", "lam8_vs_following")}
    S["one_coordinate"] = {f"T{kT}": {lab: {ax: [share(B[f"{pre}_T{kT}_{ax}_g{g}"]) for g in (0.3, 1.0, 3.0)] for ax in axs}
                                      for pre, axs, lab in (("D", "xy", "curved_cartesian"), ("Dp", "sn", "curved_path"),
                                                            ("Dt", "xy", "straight"))} for kT in (0.01, 0.05)}
    S["path_frame_xyQ_g1"] = {f"T{kT}": {"path": [share(B[f"P_T{kT}_xyQ_g1.0_s{s}"]) for s in SEEDS],
                                         "cartesian": [share(x) for x in F(f"curved_T{kT}", "xyQ", 1.0)]}
                              for kT in (0.01, 0.05)}

    # the law
    S["law_vs_trajectories"] = {}
    for c in ("tilt_T0.01", "tilt_T0.05", "fast_T0.005"):
        S["law_vs_trajectories"][c] = {"bare": share(B[f"R_{c}"])}
        for st in ("xxQ", "xyQ", "yyQ", "quad", "disp", "all"):
            rows = [x for g in (0.1, 0.3, 0.5, 1.0, 3.0) for x in F(c, st, g) if x["id"] in L and "P_law" in L[x["id"]]
                    and np.isfinite(L[x["id"]]["P_law"]) and np.isfinite(share(x))]
            S["law_vs_trajectories"][c][st] = float(np.mean([abs(share(x) - L[x["id"]]["P_law"]) for x in rows]))
    S["entropic"] = {f"g{g}_T{kT}": {"trajectories": share(B[f"A_g{g}_T{kT}"]), "se": B[f"A_g{g}_T{kT}"]["outcome"]["stderr"],
                                     "law": L[f"A_g{g}_T{kT}"]["P_law"],
                                     "law_mean_force": L[f"A_g{g}_T{kT}"]["P_law_mean_force"]}
                     for g in (0.1, 0.3) for kT in (0.005, 0.01, 0.02, 0.05, 0.1)}
    S["entropic_mean_abs"] = float(np.mean([abs(v["trajectories"] - v["law"]) for v in S["entropic"].values()]))
    S["entropic_boltzmann"] = {f"g{g}": boltzmann_entropic(g) for g in (0.1, 0.3)}

    # ---- figure 1: strength of each kind of force
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.5))
    Gd = (0.1, 0.3, 1.0, 3.0)
    Gq = (0.1, 0.3, 0.5, 1.0)
    series = [("xxQ", "x²Q: mode displaced along the path", Gd, CAT[0]), ("xyQ", "xyQ: displaced across the path", Gd, CAT[1]),
              ("yyQ", "y²Q", Gd, CAT[2]), ("quad", "xQQ, yQQ: frequency changes", Gq, CAT[3])]
    ax = axes[0]
    for st, lab, G_, col in series:
        rows = [(g, x) for g in G_ for x in F("curved_T0.01", st, g)]
        ax.plot([g for g, _ in rows], [share(x) for _, x in rows], "o", ms=4, color=col, mec="white", mew=0.6, alpha=0.9)
        ax.plot([g for g in G_ if F("curved_T0.01", st, g)],
                [np.mean([share(x) for x in F("curved_T0.01", st, g)]) for g in G_ if F("curved_T0.01", st, g)],
                "-", lw=1.5, color=col, label=lab)
    ax.axhline(share(B["R_curved_T0.01"]), color=INK2, lw=1, ls="--")
    ax.text(0.5, share(B["R_curved_T0.01"]) - 0.035, "no modes", color=INK2, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("cubic strength g")
    ax.set_ylabel("share of the product inside the bend")
    ax.set_title("curved valley, kT = 0.01", color=INK, fontsize=9)
    ax = axes[1]
    for st, lab, G_, col in series:
        for kT, sty, mk in ((0.01, "-", "o"), (0.05, "--", "s")):
            gg = [g for g in G_ if F(f"curved_T{kT}", st, g)]
            yy = [max(np.mean([friction(x) for x in F(f"curved_T{kT}", st, g)]), 1e-5) for g in gg]
            ax.plot(gg, yy, sty, marker=mk, ms=4, lw=1.3, color=col, mec="white", mew=0.5)
    ax.plot([], [], "k-", marker="o", ms=4, lw=1.2, label="kT = 0.01")
    ax.plot([], [], "k--", marker="s", ms=4, lw=1.2, label="kT = 0.05")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(3e-4, 3)
    ax.set_xlabel("cubic strength g")
    ax.set_ylabel("friction read from the energy lost")
    ax.set_title("only the frequency terms scale with temperature", color=INK, fontsize=9)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    ax = axes[2]
    for st, lab, G_, col in series:
        gg = [g for g in G_ if F("curved_T0.01", st, g)]
        ax.plot(gg, [100 * np.mean([returns(x) for x in F("curved_T0.01", st, g)]) for g in gg], "o-", ms=4, lw=1.4, color=col,
                mec="white", mew=0.6)
    ax.axhline(100 * returns(B["R_curved_T0.01"]), color=INK2, lw=1, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("cubic strength g")
    ax.set_ylabel("trajectories returning over TS1 (%)")
    ax.set_title("returns: a small energy loss stops them", color=INK, fontsize=9)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=4, fontsize=8, frameon=False)
    fig.subplots_adjust(wspace=0.3, bottom=0.25, top=0.9, left=0.06, right=0.985)
    fig.savefig(FIG / "anharm_forces.png")
    plt.close(fig)

    # ---- figure 2: the limits and the coordinates
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.6))
    ax = axes[0]
    marks = {"xxQ": "o", "xyQ": "s", "disp": "^"}
    for st in ("xxQ", "xyQ", "disp"):
        for kT in (0.01, 0.05):
            lo, hi = np.array(lim[f"{st}_T{kT}"]["slow_full_vs_frozen"]), np.array(lim[f"{st}_T{kT}"]["fast_full_vs_following"])
            ax.plot(lo[:, 1], lo[:, 0], marks[st], ms=5.5, color=CAT[0], mec="white", mew=0.6)
            ax.plot(hi[:, 1], hi[:, 0], marks[st], ms=5.5, color=CAT[1], mec="white", mew=0.6)
    ax.plot([0, 0.6], [0, 0.6], color=INK2, lw=0.9, ls=":")
    ax.plot([], [], "o", color=CAT[0], label="modes below 2: against frozen modes")
    ax.plot([], [], "o", color=CAT[1], label="modes above 2: against exact following")
    ax.plot([], [], "o", color=INK2, label="x²Q")
    ax.plot([], [], "s", color=INK2, label="xyQ")
    ax.plot([], [], "^", color=INK2, label="both, and y²Q")
    ax.set_xlabel("share in the limit (no mode dynamics)")
    ax.set_ylabel("share with the modes moving")
    ax.set_title("each band against its limit, g = 1", color=INK, fontsize=9)
    ax.legend(loc="upper left", fontsize=7.5, frameon=False)
    ax = axes[1]
    lams = (1.0, 2.0, 4.0, 8.0)
    for i, st in enumerate(("xxQ", "xyQ")):
        for s in SEEDS:
            yy = [share(B[f"F_curved_T0.01_{st}_g1.0_s{s}"])] + [share(B[f"S_T0.01_{st}_lam{lam}_s{s}"]) for lam in lams[1:]]
            ax.plot(lams, yy, "o-", ms=4, lw=1.2, color=CAT[i], mec="white", mew=0.6, alpha=0.9,
                    label=("x²Q" if st == "xxQ" else "xyQ") if s == 1 else None)
            ax.plot([11], [share(B[f"G_curved_T0.01_{st}_g1.0_s{s}"])], "D", ms=4.5, color=CAT[i], mec="white", mew=0.6)
    ax.axhline(share(B["R_curved_T0.01"]), color=INK2, lw=1, ls="--")
    ax.set_xscale("log", base=2)
    ax.set_xticks([1, 2, 4, 8, 11])
    ax.set_xticklabels(["1", "2", "4", "8", "exact"])
    ax.set_xlabel("frequencies multiplied by (displacements fixed)")
    ax.set_ylabel("share inside the bend")
    ax.set_title("faster modes converge to exact following", color=INK, fontsize=9)
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ax = axes[2]
    Gs = (0.3, 1.0, 3.0)
    for lab, pre, a_, col, sty in (("friction on s alone", "Dp", "s", CAT[0], "-"), ("friction on n alone", "Dp", "n", CAT[1], "-"),
                                   ("friction on x alone", "D", "x", CAT[0], ":"), ("friction on y alone", "D", "y", CAT[1], ":")):
        ax.plot(Gs, [share(B[f"{pre}_T0.01_{a_}_g{g}"]) for g in Gs], sty, marker="o", ms=4.5, lw=1.5, color=col, mec="white",
                mew=0.6, label=lab)
    ax.axhline(share(B["R_curved_T0.01"]), color=INK2, lw=1, ls="--")
    ax.axhline(0.5, color=INK2, lw=0.7, ls=":")
    ax.set_xscale("log")
    ax.set_xticks(Gs)
    ax.set_xticklabels(["0.3", "1", "3"])
    ax.minorticks_off()
    ax.set_xlabel("friction gamma on one coordinate")
    ax.set_ylabel("share inside the bend")
    ax.set_title("friction along and across the path, and along fixed axes", color=INK, fontsize=9)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    fig.subplots_adjust(wspace=0.3, bottom=0.15, top=0.91, left=0.06, right=0.985)
    fig.savefig(FIG / "anharm_limits.png")
    plt.close(fig)

    # ---- figure 3: modes that stiffen on one side
    kTs = (0.005, 0.01, 0.02, 0.05, 0.1)
    fig, ax = plt.subplots(figsize=(5.2, 3.5))
    for i, g in enumerate((0.1, 0.3)):
        ax.errorbar(kTs, [share(B[f"A_g{g}_T{kT}"]) for kT in kTs], yerr=[B[f"A_g{g}_T{kT}"]["outcome"]["stderr"] for kT in kTs],
                    fmt="o", ms=5, color=CAT[i], mec="white", mew=0.7, label=f"trajectories, g = {g:g}")
        ax.plot(kTs, [L[f"A_g{g}_T{kT}"]["P_law"] for kT in kTs], "-", lw=1.5, color=CAT[i], label="law" if i == 0 else None)
        ax.axhline(boltzmann_entropic(g), color=CAT[i], lw=1, ls="--", label="free energies of the products" if i == 0 else None)
    ax.axhline(0.5, color=INK2, lw=0.7, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("kT")
    ax.set_ylabel("share of the side where the modes stiffen")
    ax.set_title("no energy difference: modes stiffening on one side", color=INK, fontsize=9)
    ax.legend(loc="center left", bbox_to_anchor=(0.01, 0.33), fontsize=7.5, frameon=False)
    fig.subplots_adjust(left=0.13, right=0.98, top=0.92, bottom=0.14)
    fig.savefig(FIG / "anharm_entropic.png")
    plt.close(fig)

    (HERE / "results/anharm_summary.json").write_text(json.dumps(S, indent=1))
    print(json.dumps({k: v for k, v in S.items() if k not in ("limits", "entropic")}, indent=0))


if __name__ == "__main__":
    main()
