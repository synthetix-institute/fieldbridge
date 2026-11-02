#!/usr/bin/env python3
"""Figure of results.json: finite populations of agents against their mean-field equations.

Usage: python3 -B figure.py [results.json] [agents_vs_equations.png]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
EQUATION = "#222222"
SIZES = ["#9ecae1", "#4292c6", "#08519c", "#08306b"]   # darker for a larger population
SECOND = "#c8553d"


def by_size(ax, x, series, marker="o", ms=4):
    for colour, (N, y) in zip(SIZES, series.items()):
        ax.plot(x, y, marker, color=colour, ms=ms, label=f"N = {int(N):,}".replace(",", " "), zorder=3)


def panel_kuramoto(ax, k):
    K = np.linspace(0.3, 3.05, 400)
    ax.plot(K, np.sqrt(np.clip(1 - 2 * k["delta"] / K, 0, None)), color=EQUATION, lw=1.5, label="equation")
    by_size(ax, k["K"], k["r"])
    ax.set(xlabel="coupling K", ylabel="order parameter r", title="Oscillators: continuous transition")
    ax.legend(frameon=False, fontsize=8)


def panel_cycle(ax, c):
    mu = np.array(c["mu"])
    ok = np.array(c["meanfield_rms"]) > 1e-3
    ax.plot(mu[ok], np.array(c["meanfield_rms"])[ok], color=EQUATION, lw=1.5, label="equation (limit cycle)")
    by_size(ax, mu, c["rms"])
    above = [i for i, v in enumerate(c["linear_noise_variance_times_N"]) if v is not None]
    for j, (colour, N) in enumerate(zip(SIZES, c["N"])):
        y = np.sqrt(np.array([c["linear_noise_variance_times_N"][i] for i in above]) / N)
        ax.plot(mu[above], y, "--", color=colour, lw=1, label="linear noise" if j == len(c["N"]) - 1 else None)
    ax.axvline(c["mu_h"], color="#999999", lw=0.8, zorder=0)
    ax.set_yscale("log")
    ax.set(xlabel="mutation rate μ", ylabel="rms deviation of a frequency", title="Cyclic dominance: Hopf bifurcation")
    ax.legend(frameon=False, fontsize=8)


def panel_error(ax, e):
    q = np.linspace(0.4, 0.9, 200)
    ax.plot(q, np.clip((e["f0"] * q - 1) / (e["f0"] - 1), 0, None), color=EQUATION, lw=1.5, label="equation")
    by_size(ax, e["q"], e["x"])
    ax.set(xlabel="copying fidelity q", ylabel="master fraction x", title="Error threshold: transcritical bifurcation")
    ax.legend(frameon=False, fontsize=8)


def panel_hysteresis(ax, b):
    r = b["r"]
    f = lambda u, k: r * u * (1 - u / k) - u * u / (1 + u * u)
    k = np.linspace(3, 20, 600)
    u = np.linspace(0.05, 20, 4000)
    # the states of the equation are the sign changes of f along u at each capacity
    for kv, style in (("stable", "-"), ("unstable", ":")):
        xs, ys = [], []
        for kk in k:
            v = f(u, kk)
            for i in np.nonzero(v[:-1] * v[1:] < 0)[0]:
                root = brentq(f, u[i], u[i + 1], args=(kk,))
                stable = v[i] > 0
                if stable == (kv == "stable"):
                    xs.append(kk)
                    ys.append(root)
        ax.plot(xs, ys, ".", color=EQUATION if kv == "stable" else "#999999", ms=1.5,
                label="stable states" if kv == "stable" else "unstable state")
    run = next(s for s in b["sweeps"] if s["k"] is not None)
    ax.plot(run["k"], run["up"], color=SIZES[1], lw=0.8, label="agents, k rising")
    ax.plot(run["k"][::-1], run["down"], color=SECOND, lw=0.8, label="agents, k falling")
    ax.set(xlabel="capacity k", ylabel="population density u", yscale="log", ylim=(0.4, 25),
           title=f"Budworm: hysteresis at N = {run['N']:,}".replace(",", " "))
    ax.legend(frameon=False, fontsize=8, markerscale=4)


def constant_noise_line(N, ln_tau, slope):
    """A line with the slope of the constant-noise action, starting from the exact value at the smallest N."""
    return ln_tau[0] + slope * (N - N[0])


def panel_lifetime(ax, e):
    rows = e["lifetime"]
    N = np.array([x["N"] for x in rows], float)
    tau = np.array([x["ln_tau_exact"] for x in rows])
    ax.plot(N, tau, color=EQUATION, lw=1.5, label="exact chain")
    sim = [x for x in rows if "ln_tau_sim" in x]
    ax.errorbar([x["N"] for x in sim], [x["ln_tau_sim"] for x in sim], [x["ln_tau_sim_se"] for x in sim], fmt="o", ms=4,
                color=SIZES[2], label="agents", zorder=3)
    ax.plot(N, constant_noise_line(N, tau, e["constant_noise_action"]), "--", color=SECOND, lw=1.2,
            label="constant noise")
    ax.set(xlabel="population size N", ylabel="ln (time to lose the master)", title="Error threshold: lifetime at q = 0.8")
    ax.legend(frameon=False, fontsize=8)


def panel_switching(ax, b):
    rows = b["switching"]
    N = np.array([x["N"] for x in rows], float)
    for key, label, slope, colour in (("down", "outbreak to low", b["constant_noise_action_down"], SIZES[2]),
                                      ("up", "low to outbreak", b["constant_noise_action_up"], SIZES[1])):
        tau = np.array([x[f"ln_tau_{key}"] for x in rows])
        ax.plot(N, tau, color=EQUATION, lw=1.5, label="exact chain" if key == "down" else None)
        ax.plot(N, constant_noise_line(N, tau, slope), "--", color=SECOND, lw=1.2,
                label="constant noise" if key == "down" else None)
        sim = [x for x in rows if f"ln_tau_{key}_sim" in x]
        ax.errorbar([x["N"] for x in sim], [x[f"ln_tau_{key}_sim"] for x in sim], [x[f"ln_tau_{key}_sim_se"] for x in sim],
                    fmt="o", ms=3.5, color=colour, label=f"agents, {label}", zorder=3)
    ax.set(xlabel="population size N", ylabel="ln (time to switch)", title="Budworm: switching at k = 8")
    ax.legend(frameon=False, fontsize=8)


def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "results.json"
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "agents_vs_equations.png"
    res = json.loads(src.read_text(encoding="utf-8"))
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.4))
    panel_kuramoto(axes[0, 0], res["kuramoto"])
    panel_cycle(axes[0, 1], res["cycle"])
    panel_error(axes[0, 2], res["error"])
    panel_hysteresis(axes[1, 0], res["budworm"])
    panel_lifetime(axes[1, 1], res["error"])
    panel_switching(axes[1, 2], res["budworm"])
    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)
        ax.title.set_fontsize(11)
    fig.tight_layout()
    fig.savefig(out, dpi=170)
    print("wrote", out)


if __name__ == "__main__":
    main()
