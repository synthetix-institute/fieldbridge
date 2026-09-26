"""Figures for the memory modules of the tutorial (Modules 1-9).

Run the tutorial commands first; they write their results under build/tut (see the modules). Then

    python3 -B docs/tutorial/figures/memory/make_figures.py [--build build/tut]

copies the card, loop and co-discovery figures and draws the phase, regime and field figures from the saved JSON results.
No calculation is repeated here: every curve is read from a result file.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
COPIES = {"m1_toggle/card.png": "m1_toggle_card.png", "m3_loops/loops.png": "m3_loops.png",
          "m4_colloid/card.png": "m4_colloid_card.png", "m9_codiscover/codiscover.png": "m9_codiscovery.png",
          "m9_threshold/codiscover.png": "m9_threshold.png", "m9_phase/codiscover.png": "m9_phase.png"}
COLORS = {"contracting": "#1f77b4", "flat": "#2ca02c", "expanding": "#d62728"}


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def label(ax, text: str) -> None:
    ax.text(-0.14, 1.04, text, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def phase_figure(build: Path, out: Path) -> None:
    p = load(build / "m6_phase" / "phase.json")
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    ax = axes[0]
    ax.plot(p["response_curve"]["phase"], p["response_curve"]["shift"], "o-", ms=4, color="#1f77b4")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("phase at which the kick arrives (rad)")
    ax.set_ylabel("phase shift (rad)")
    ax.set_title("Phase response to a kick", fontsize=10)
    label(ax, "a")
    ax = axes[1]
    h = p["hold"]
    t = np.asarray(h["time"])
    ax.plot(t, h["separation"], "o-", ms=3, color="#1f77b4", label="separation of the written phase")
    ax.plot(t, np.sqrt(h["variance"]), "s-", ms=3, color="#ff7f0e", label="spread of the phase")
    ax.plot(t, np.sqrt(2 * h["phase_diffusion"] * t), "--", color="#ff7f0e", lw=1,
            label=r"$(2D_\phi t)^{1/2}$")
    ax.set_xlabel("time after the write")
    ax.set_ylabel("rad")
    ax.set_title("Field-free retention of a phase shift", fontsize=10)
    ax.set_ylim(0, 2.4)
    ax.legend(fontsize=8, frameon=False, loc="upper center", ncol=1)
    label(ax, "b")
    ax = axes[2]
    w = p["lock"]["writes"]
    hs = np.array([x["h"] for x in w])
    ratio = np.array([x["ratio"] for x in w])
    ax.loglog(hs, ratio, "o", color="#2ca02c", label="calculated")
    ax.loglog(hs, ratio[0] * hs / hs[0], "--", color="0.4", lw=1, label="slope 1")
    ax.set_xlabel("relative change of the promoter strength")
    ax.set_ylabel("retention time / writing time")
    ax.set_title(f"Log-log slope {p['lock']['log_log_slope']:.2f}", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    label(ax, "c")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def regimes_figure(build: Path, out: Path) -> None:
    r = load(build / "m8_regimes" / "regimes.json")
    t = np.asarray(r["times"])
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    names = {"contracting": r"$\kappa > 0$ (relaxation)", "flat": r"$\kappa = 0$ (diffusion)",
             "expanding": r"$\kappa < 0$ (instability)"}
    ax = axes[0]
    for key, reg in r["regimes"].items():
        ax.plot(t, reg["gaussian_information"], color=COLORS[key], label=names[key])
    ax.axhline(r["plateau_gaussian"], color=COLORS["expanding"], ls=":", lw=1)
    ax.text(t[0], r["plateau_gaussian"] * 1.07, r"$\frac{1}{2}\ln(1+W)$", color=COLORS["expanding"], fontsize=9)
    ax.set_xscale("log")
    ax.set_ylim(0, 2.1)
    ax.set_xlabel(r"time after the write ($1/|\kappa|$)")
    ax.set_ylabel("information (nats)")
    ax.set_title("Continuous observable: Eq. (1)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    label(ax, "a")
    ax = axes[1]
    for key, reg in r["regimes"].items():
        if reg.get("binary_exact_linear"):
            ax.plot(t, reg["binary_exact_linear"], color=COLORS[key], lw=1.2)
        ax.plot(reg["simulation"]["times"], reg["simulation"]["information"], "o", ms=3.5, color=COLORS[key],
                mfc="none", label=names[key])
    ax.axhline(np.log(2), color="0.6", lw=0.8)
    ax.text(t[0], np.log(2) * 1.03, r"$\ln 2$", color="0.4", fontsize=9)
    ax.axhline(r["binary_plateau"]["information"], color=COLORS["expanding"], ls=":", lw=1)
    ax.text(t[-1], r["binary_plateau"]["information"] * 1.04, r"$\Phi(W^{1/2})$", color=COLORS["expanding"],
            fontsize=9, ha="right")
    ax.set_xscale("log")
    ax.set_ylim(0, 0.8)
    ax.set_xlabel(r"time after the write ($1/|\kappa|$)")
    ax.set_ylabel("information (nats)")
    ax.set_title("Sign of the observable: exact (lines), simulated (circles)", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    label(ax, "b")
    ax = axes[2]
    for key, reg in r["regimes"].items():
        ax.plot(reg["simulation"]["times"], reg["simulation"]["accuracy"], "o-", ms=3.5, color=COLORS[key],
                label=names[key])
    ax.axhline(0.5, color="0.6", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel(r"time after the write ($1/|\kappa|$)")
    ax.set_ylabel("fraction identified correctly")
    ax.set_title(f"Simulation, {r['input']['n']} trajectories per regime", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    label(ax, "c")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def fields_figure(build: Path, out: Path) -> None:
    runs = [("m8_conserved_1d_charge", r"$d=1$, adds material", "#1f77b4"),
            ("m8_conserved_2d_charge", r"$d=2$, adds material", "#9467bd"),
            ("m8_conserved_1d_dipole", r"$d=1$, moves material", "#2ca02c"),
            ("m8_conserved_2d_dipole", r"$d=2$, moves material", "#d62728")]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
    ax = axes[0]
    for name, lab, col in runs:
        f = load(build / name / "field.json")
        ex = f["exponent_check"]
        tt, snr = np.asarray(ex["times"]), np.asarray(ex["snr"])
        ax.loglog(tt, snr / snr[0], "o", ms=3.5, color=col, label=f"{lab}: {ex['exponent']:.2f}")
        ax.loglog(tt, (tt / tt[0]) ** f["prediction"]["exponent"], "-", lw=0.8, color=col)
    ax.set_xlabel("time after the write")
    ax.set_ylabel("signal-to-noise ratio (relative)")
    ax.set_title("Linear field, large lattice: exponents", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False, loc="lower left")
    label(ax, "a")
    ax = axes[1]
    for name, lab, col in runs:
        f = load(build / name / "field.json")
        ex, sim, ren = f["exact_same_size"], f["simulation"], f["renormalized"]
        ax.loglog(ex["times"], ex["snr_profile"], ":", color=col, lw=1)
        ax.loglog(sim["times"], sim["snr_profile"], "o", ms=3.5, color=col, mfc="none")
        ax.loglog(ex["times"], ren["snr_profile"], "-", color=col, lw=1.1)
    ax.set_xlabel("time after the write")
    ax.set_ylabel("signal-to-noise ratio")
    ax.set_title("Nonlinear field (circles), linear (dotted),\nrenormalized linear (solid)", fontsize=10)
    label(ax, "b")
    ax = axes[2]
    f = load(build / "m8_nonconserved_1d" / "field.json")
    ex, sim = f["exact_same_size"], f["simulation"]
    keep = np.asarray(sim["snr_profile"]) > 0
    ax.semilogy(np.asarray(ex["times"]), ex["snr_profile"], ":", color="#1f77b4", label="linear field")
    ax.semilogy(np.asarray(sim["times"])[keep], np.asarray(sim["snr_profile"])[keep], "o", ms=3.5, color="#d62728",
                mfc="none", label="nonlinear field")
    t0, t1 = f["window"]
    tw = np.linspace(t0, t1, 20)
    i0 = int(np.argmin(np.abs(np.asarray(sim["times"]) - t0)))
    s0 = sim["snr_profile"][i0]
    ax.semilogy(tw, s0 * np.exp(-f["hartree"]["rate"] * (tw - sim["times"][i0])), "-", color="#d62728", lw=1,
                label=f"self-consistent mass: rate {f['hartree']['rate']:.2f}")
    ax.set_ylim(1e-4, 10)
    ax.set_xlabel("time after the write")
    ax.set_ylabel("signal-to-noise ratio")
    ax.set_title("Non-conserved field, $d=1$", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False)
    label(ax, "c")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--build", default="build/tut", help="directory with the tutorial results")
    args = ap.parse_args()
    build = Path(args.build)
    for src, dst in COPIES.items():
        shutil.copyfile(build / src, HERE / dst)
    phase_figure(build, HERE / "m6_phase.png")
    regimes_figure(build, HERE / "m8_regimes.png")
    fields_figure(build, HERE / "m8_fields.png")
    print("figures written to", HERE)


if __name__ == "__main__":
    main()
