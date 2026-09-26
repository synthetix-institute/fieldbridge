"""Figures and an HTML gallery for evaluated realizations and composed loops.

card_figure(card, path)   four panels: (a) stable states on the carrier; (b) stable and unstable states against
                          the control, with the write point (or the stable state under a growing uniform field);
                          (c) writing protocols and field-free retention; (d) the drift along the unstable direction
                          at the write point, compared with the pitchfork normal form
loops_figure(rows, path)  prediction from the holonomy against measurement, for genes, spins and rotors
codiscovery_figure(report, path)  derivations of one target in realizations from different fields, and the invariants
                          (symmetric write, threshold write or phase locking)
html_report(cards, path)  one self-contained page with every card, its identity slots and its classification
"""
from __future__ import annotations

import base64
import html
import io
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

INK, INK2, GRID = "#1f2430", "#6b7280", "#e5e7eb"
BLUE, ORANGE, GREEN, PURPLE, RED, CYAN = "#2563eb", "#ea580c", "#059669", "#7c3aed", "#dc2626", "#0891b2"
STATE_COLORS = [BLUE, ORANGE, GREEN, PURPLE, RED, CYAN]
KIND_COLOR = {"supercritical": GREEN, "subcritical": PURPLE, "saddle": ORANGE, "transcritical": RED,
              "imperfect": CYAN, "Hopf": RED, "degenerate": INK2}


def _kind_key(kind: str) -> str:
    for k in KIND_COLOR:
        if kind.startswith(k):
            return k
    return "degenerate"


def _style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7, colors=INK2)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for lab in (ax.xaxis.label, ax.yaxis.label):
        lab.set_fontsize(7.5)
        lab.set_color(INK)


def _title(ax, letter: str, text: str):
    ax.set_title(f"{letter}   {text}", loc="left", fontsize=8.5, fontweight="bold", color=INK)


def _note(ax, text: str, y: float = 0.5, size: float = 7.5):
    ax.text(0.5, y, text, transform=ax.transAxes, ha="center", va="center", fontsize=size, color=INK2, wrap=True)


def _sticks(ax, pos, theta, m, color, lw, length):
    for (x, y), th in zip(pos, theta):
        dx, dy = length * np.cos(th), length * np.sin(th)
        if m == 2:
            ax.plot([x - dx, x + dx], [y - dy, y + dy], color=color, lw=lw, solid_capstyle="round")
        else:
            ax.annotate("", xy=(x + dx, y + dy), xytext=(x - dx, y - dy),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=7))


# ------------------------------------------------------------------------------------------------ card panels
def panel_states(ax, card):
    st, car, geo = card["states"], card["carrier"], card.get("geometry")
    qs = [np.asarray(q) for q in st["q"]]
    if not qs and card.get("trajectory"):
        tr = card["trajectory"]
        Q = np.asarray(tr["q"])
        for i in range(Q.shape[1]):
            ax.plot(tr["t"], Q[:, i], lw=1.1, color=STATE_COLORS[i % 6], label=f"q_{i + 1}")
        ax.set_xlabel("time (noise-free)"); ax.set_ylabel("state")
        ax.legend(fontsize=6.3, loc="upper right", frameon=False, ncol=min(4, Q.shape[1]))
        ax.text(0.02, 0.02, "no stable state: sustained oscillation", transform=ax.transAxes, fontsize=7, color=INK2)
        _style(ax)
        return
    m = int(round(2 * np.pi / car["period"])) if car["kind"] == "torus" else 0
    if geo and car["kind"] == "torus" and qs:
        pos = np.asarray(geo["pos"])
        for i, j in zip(geo["src"], geo["tgt"]):
            ax.plot(pos[[i, j], 0], pos[[i, j], 1], color=GRID, lw=1.0, zorder=0)
        bond = float(np.median([np.linalg.norm(pos[i] - pos[j]) for i, j in zip(geo["src"], geo["tgt"])]))
        for k, q in enumerate(qs[:2]):
            _sticks(ax, pos, q, m, STATE_COLORS[k], 2.4 if k == 0 else 1.3, 0.3 * bond)
        ax.plot([], [], color=BLUE, lw=2.4, label="stable state 1")
        if len(qs) > 1:
            ax.plot([], [], color=ORANGE, lw=1.3, label="stable state 2")
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.spines[:].set_visible(False)
        ax.legend(fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, 0.02), ncol=2, frameon=False)
        return
    if car["dim"] == 1 and "landscape" in card:
        ls = card["landscape"]
        x = np.asarray(ls["x"])
        if "V" in ls:
            y = np.asarray(ls["V"])
            ax.plot(x, y, color=INK, lw=1.2)
            ax.set_ylabel("potential V")
            val = lambda z: float(np.interp(z, x, y))
        else:
            y = np.asarray(ls["F"])
            ax.plot(x, y, color=INK, lw=1.2)
            ax.axhline(0, color=INK2, lw=0.7, ls="--")
            ax.set_ylabel("drift F")
            val = lambda z: 0.0
        for k, q in enumerate(qs):
            ax.plot(q[0], val(q[0]), "o", color=STATE_COLORS[k % 6], ms=6)
        ax.set_xlabel("carrier coordinate")
        _style(ax)
        return
    if car["dim"] == 2 and "landscape" in card:
        ls = card["landscape"]
        g = np.asarray(ls["grid"])
        U, V = np.asarray(ls["U"]), np.asarray(ls["V"])
        speed = np.hypot(U, V)
        ax.streamplot(g, g, U, V, color=np.log10(speed + 1e-9), cmap="Greys", density=0.9, linewidth=0.6, arrowsize=0.6)
        ax.contour(g, g, U, levels=[0], colors=[BLUE], linewidths=1.0)
        ax.contour(g, g, V, levels=[0], colors=[ORANGE], linewidths=1.0)
        for k, q in enumerate(qs):
            ax.plot(q[0], q[1], "o", color=INK, ms=6, zorder=5)
        ax.plot([], [], color=BLUE, lw=1, label="nullcline 1")
        ax.plot([], [], color=ORANGE, lw=1, label="nullcline 2")
        ax.legend(fontsize=6.3, loc="upper right", frameon=False)
        ax.set_xlabel("coordinate 1"); ax.set_ylabel("coordinate 2")
        _style(ax)
        return
    if qs:
        Q = np.asarray(qs[:12])
        cmap = "viridis" if car["kind"] == "orthant" else "RdBu_r"
        lim = float(np.max(np.abs(Q))) or 1.0
        im = ax.imshow(Q, aspect="auto", cmap=cmap, vmin=0 if car["kind"] == "orthant" else -lim, vmax=lim)
        ax.set_xlabel("coordinate"); ax.set_ylabel("stable state")
        ax.set_yticks(range(len(Q))); ax.set_yticklabels([str(i + 1) for i in range(len(Q))], fontsize=6)
        plt.colorbar(im, ax=ax, fraction=0.05, pad=0.02).ax.tick_params(labelsize=6)
        ax.tick_params(labelsize=7)
        return
    ax.axis("off")
    _note(ax, card["loss_law"])


def panel_kappa(ax, card):
    st = card["states"]
    if not st["kappa"]:
        ax.axis("off")
        _note(ax, "no stable state:\nevery fixed point has kappa < 0 along some direction")
        return
    for k, spec in enumerate(st["kappa"][:6]):
        re = sorted(x[0] for x in spec)
        ax.plot(np.arange(1, len(re) + 1), re, marker="o", ms=3, lw=0.8, color=STATE_COLORS[k],
                label=f"state {k + 1}" if k < 4 else None)
    for e in card["construct"]["events"]:
        kc = e["kappa_critical"][0]
        ax.plot(0, kc, marker="*", ms=9, color=KIND_COLOR[_kind_key(e["kind"])], zorder=6)
    ax.axhline(0, color=INK2, lw=0.8, ls="--")
    vals = [abs(x[0]) for spec in st["kappa"] for x in spec if abs(x[0]) > 1e-6]
    if vals and max(vals) / min(vals) > 30:
        ax.set_yscale("symlog", linthresh=max(min(vals), 1e-3))
    top = max(vals + [1.0])
    ax.set_ylim(min(-0.08 * top, ax.get_ylim()[0]), 1.15 * top if ax.get_yscale() == "linear" else ax.get_ylim()[1])
    ax.set_xlim(-0.5, max(len(sp) for sp in st["kappa"]) + 0.5)
    ax.set_xlabel("direction (sorted); star at 0: kappa of the critical direction at a write point")
    ax.set_ylabel("kappa = -Re eigenvalue")
    ax.legend(fontsize=6.3, loc="best", frameon=False)
    _style(ax)


def panel_bifurcation(ax, card):
    scan = card.get("bifurcation")
    evs = card["construct"]["events"]
    if scan:
        allv = [p["summary"] for row in scan for p in row["points"]]
        key = "summary" if (allv and np.ptp(allv) > 1e-6) or card["carrier"]["dim"] == 1 else "mean"
        for row in scan:
            for p in row["points"]:
                ax.plot(row["value"], p[key], "o", ms=3.2, color=INK if p["stable"] else "white",
                        markeredgecolor=INK if p["stable"] else INK2, markeredgewidth=0.8)
        ymax = ax.get_ylim()[1]
        for e in evs:
            if e["source"] != "control":
                continue
            c = KIND_COLOR[_kind_key(e["kind"])]
            ax.axvline(e["value"], color=c, lw=1.0, ls=":")
            ax.text(e["value"], ymax, " write point: " + e["kind"].split(":")[0].split(" (")[0], color=c, fontsize=6.3,
                    rotation=90, va="top")
        car = card["carrier"]
        ax.set_ylabel("mean state coordinate" if key == "mean" else
                      ("net orientational order |<exp(i m q)>|" if car["kind"] == "torus" else
                       ("state coordinate q" if car["dim"] == 1 else "q_1 - mean of the other coordinates")))
        ax.set_xlabel(f"control parameter {card['control']}  (filled: stable, open: unstable)")
        _style(ax)
        return
    wf = card["construct"].get("write_field")
    if wf and wf["branch"]:
        h = [b["h"] for b in wf["branch"]]
        ax.plot(h, [b["overlap"] for b in wf["branch"]], color=INK, lw=1.2, label="stable state 1: overlap with state 2")
        ax2 = ax.twinx()
        ax2.plot(h, [b["kappa_min"] for b in wf["branch"]], color=BLUE, lw=1.0, ls="--")
        ax2.set_ylabel("smallest relaxation rate kappa", color=BLUE, fontsize=7.5)
        ax2.tick_params(labelsize=7, colors=BLUE)
        if "h_c" in wf:
            ax.axvline(wf["h_c"], color=ORANGE, ls=":", lw=1.0)
            ax.annotate(f"state 1 ceases to exist\nat h = {wf['h_c']:.3g}", xy=(wf["h_c"], ax.get_ylim()[1]),
                        xytext=(-4, -4), textcoords="offset points", ha="right", va="top", color=ORANGE,
                        fontsize=6.8)
        ax.set_xlabel("uniform field h toward state 2")
        ax.set_ylabel("overlap with state 2")
        ax.legend(fontsize=6.3, loc="lower right", frameon=False)
        _style(ax)
        return
    ax.axis("off")
    _note(ax, "no control parameter\nand fewer than two stable states")


def panel_normal_form(ax, card):
    evs = card["construct"]["events"]
    real_evs = [e for e in evs if not e.get("hopf") and "s" in e]
    hopf = [e for e in evs if e.get("hopf")]
    if not real_evs:
        ax.axis("off")
        if hopf:
            e = hopf[0]
            _note(ax, f"Hopf bifurcation at {e['param']} = {e['value']:.4g}\na complex pair of eigenvalues crosses,"
                      f" frequency {e['omega']:.3g}\nno stable state is selected; the phase of the cycle\ncan retain the "
                      f"timing of a pulse", size=7)
        else:
            _note(ax, "no bifurcation along the control")
        return
    e = sorted(real_evs, key=lambda x: x["source"] != "control")[0]
    s, f = np.asarray(e["s"]), np.asarray(e["f"])
    fit = e["a0"] + e["a1"] * s + e["a2"] * s ** 2 + e["a3"] * s ** 3
    ax.plot(s, f, "o", ms=2.3, color=INK, label="drift along the unstable direction")
    ax.plot(s, fit, color=BLUE, lw=1.1, label="cubic fit")
    ax.plot(s, e["a1"] * s + e["a3"] * s ** 3, color=ORANGE, lw=1.0, ls="--", label="symmetric part a1 s + a3 s^3")
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xlabel(f"coordinate s along the unstable direction ({e['param']} = {e['value']:.4g})")
    ax.set_ylabel("ds/dt")
    lines = [e["kind"].split(":")[0], f"a2 = {e['a2']:.3g},  a3 = {e['a3']:.3g},  bias {e['bias']:.2f}"]
    c = e.get("cancellation")
    if c:
        nf = c["normal_form_along_sweep"]
        d = ", ".join(f"d{k} = {v:+.2f}" for k, v in c["direction"].items())
        lines += [f"cusp: {c['second']} = {c['value']:.4g}, {e['param']} = {c['control']:.4g}",
                  f"along the sweep {d}: {nf['kind'].split(':')[0]}"]
    law = e.get("write_law")
    if law:
        lines.append(f"swept write: accuracy predicted {law['predicted']:.2f}, measured {law['measured']:.2f} "
                     f"+/- {law['stderr']:.2f} (Gamma = {law['gamma']:.2g})")
    ax.text(0.02, 0.98, "\n".join(lines), transform=ax.transAxes, va="top", fontsize=6.3, color=INK,
            bbox=dict(facecolor="white", edgecolor=GRID, boxstyle="round,pad=0.3", alpha=0.9))
    ax.legend(fontsize=6.0, loc="lower right", frameon=False)
    _style(ax)


def _wrap(text: str, width: int) -> List[str]:
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width:
            lines.append(cur)
            cur = w
        else:
            cur = (cur + " " + w).strip()
    return lines + [cur]


def panel_writes(ax, card):
    from .discovery import sweep_label
    writes = card.get("writes")
    if not writes:
        ax.axis("off")
        _note(ax, "fewer than two discrete stable states:\nno writing test")
        return
    cols = [ORANGE, BLUE, GREEN]
    for w, c in zip(writes, cols):
        lab = f"field {w['h_over_threshold']:.1f} x threshold, {sweep_label(card, w)}: accuracy {w['accuracy']:.2f}"
        ax.plot(w["t"], w["fraction"], color=c, lw=1.2, label=lab)
    ax.axvline(writes[0]["T_w"], color=INK2, ls=":", lw=0.8)
    ax.text(writes[0]["T_w"], 1.0, " release", fontsize=6.3, color=INK2, va="top")
    hold = card.get("hold")
    if hold:
        ax.plot(hold["t"], hold["fraction"], color=INK2, lw=1.0, ls="--",
                label=f"no field: fraction retained in state 1 (retention time {'> ' if hold['censored'] else ''}"
                      f"{hold['time']:.3g})")
    ax.set_ylim(-0.02, 1.05)
    tmax = max([max(w["t"]) for w in writes] + ([max(hold["t"])] if hold else []))
    ax.set_xscale("symlog", linthresh=max(1e-3, writes[0]["T_w"] / 50))
    ax.set_xlim(0, tmax)
    ax.set_xlabel("time"); ax.set_ylabel("fraction in the target state")
    ax.legend(fontsize=6.0, loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False, ncol=1)
    _style(ax)


def panel_tests(ax, card):
    if card.get("transfers"):
        tr = card["transfers"]
        names = ["\n".join(_wrap(t["note"] or t["target"], 26)) for t in tr]
        y = np.arange(len(tr))
        ax.barh(y - 0.18, [max(t["relative_defect"], 1e-16) for t in tr], height=0.34, color=INK2, label="raw")
        ax.barh(y + 0.18, [max(t["relative_defect_rescaled"], 1e-16) for t in tr], height=0.34, color=BLUE,
                label="after time rescaling")
        for yi, t in zip(y, tr):
            ax.text(1.5e-16, yi + 0.42, f"c = {t['time_rescaling']:.3g}" + ("  exact" if t["exact"] else ""), fontsize=6.3,
                    color=GREEN if t["exact"] else INK)
        ax.set_xscale("log"); ax.set_xlim(1e-16, 2)
        ax.set_yticks(y); ax.set_yticklabels(names, fontsize=6.3)
        ax.set_xlabel("relative intertwining defect |F_t(a(q)) - Da F_s(q)|")
        ax.legend(fontsize=6.0, loc="lower right", frameon=False)
        _style(ax)
        return
    if card.get("kernel"):
        k = card["kernel"]
        t, K = np.asarray(k["t"]), np.asarray(k["K"])
        ax.plot(t, K, color=BLUE, lw=1.4, label="K(t) = B exp(D t) C")
        b, c, dz = k["blocks"]["B"][0][0], k["blocks"]["C"][0][0], k["blocks"]["D"][0][0]
        ax.plot(t, b * c * np.exp(dz * t), color=ORANGE, lw=1.0, ls="--", label=f"b c exp(-{-dz:g} t)")
        ax.set_xlabel("time"); ax.set_ylabel("memory kernel of the observable")
        ax.legend(fontsize=6.5, frameon=False)
        _style(ax)
        return
    evs = [e for e in card["construct"]["events"] if e.get("write_law")]
    if evs:
        law = evs[0]["write_law"]
        own = law["at_material_noise"]
        ax.bar([0, 1, 2], [law["predicted"], law["measured"], own["measured"]], yerr=[0, law["stderr"], own["stderr"]],
               color=[ORANGE, BLUE, INK2], width=0.55, capsize=3)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["pitchfork law,\ntransferred", f"target,\nGamma = {law['gamma']:.2g}",
                            f"target at its noise,\nGamma = {own['gamma']:.2g}"], fontsize=6.8)
        ax.set_ylim(0.4, 1.0)
        ax.set_ylabel("accuracy of a swept write")
        ax.text(0.5, 0.96, "Gamma = D_s |a3| / r: noise against the growth of the barrier", transform=ax.transAxes,
                ha="center", fontsize=6.3, color=INK2)
        _style(ax)
        return
    writes = card.get("writes")
    if writes and any(w.get("lock_ratio") for w in writes):
        lab = [f"{w['mode']}\n{w['h_over_threshold']:.1f}x" for w in writes]
        val = [w["lock_ratio"] or np.nan for w in writes]
        ax.bar(range(len(val)), val, color=[ORANGE, BLUE, GREEN][:len(val)], width=0.55)
        ax.set_xticks(range(len(val))); ax.set_xticklabels(lab, fontsize=6.8)
        ax.set_yscale("log"); ax.set_ylabel("retention time / rewriting time")
        from matplotlib.ticker import FuncFormatter, NullFormatter
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
        for i, (w, v) in enumerate(zip(writes, val)):
            ax.text(i, v if np.isfinite(v) else ax.get_ylim()[0], f" acc {w['accuracy']:.2f}" + ("" if np.isfinite(v) else "\n not rewritten"),
                    ha="center", va="bottom", fontsize=6.3, color=INK)
        _style(ax)
        return
    ax.axis("off")
    _note(ax, "no transfer, kernel or write test\nfor this realization")


def card_figure(card: Dict, path: Path):
    """Four panels read in order: which states exist, where a new state can be written, how writing and retention
    behave, and how the local dynamics at the write point compares with the pitchfork normal form. The relaxation
    spectra, transfer defects and memory kernels are reported in card.json and card.md."""
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 8.2), constrained_layout=True)
    (a, b), (c, d) = axes
    panel_states(a, card); _title(a, "a", "Stable states on the carrier")
    panel_bifurcation(b, card)
    by_field = not card.get("bifurcation") and card["construct"].get("write_field", {}).get("branch")
    _title(b, "b", "Stored state under a uniform field" if by_field else "States against the control parameter")
    panel_writes(c, card); _title(c, "c", "Writing protocols and field-free retention")
    panel_normal_form(d, card); _title(d, "d", "Drift along the unstable direction at the write point")
    fig.suptitle(card["name"], fontsize=10.5, color=INK, x=0.01, ha="left", fontweight="bold")
    fig.savefig(path, dpi=120)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ loops
def loops_figure(rows: List[Dict], path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.4), constrained_layout=True)
    (a, b), (c, d) = axes
    _grid_panel(a, [r for r in rows if r["group"] == "sign" and not r["reciprocal"]], "gene")
    _title(a, "a", "Gene loops: number of stable states")
    _grid_panel(b, [r for r in rows if r["group"] == "sign" and r["reciprocal"]], "spin")
    _title(b, "b", "Loops of bistable spins: unsatisfied bonds")
    rot = [r for r in rows if r["group"] == "circle"]
    phi = np.linspace(-np.pi, np.pi, 300)
    markers = {3: "^", 4: "s", 5: "D", 6: "o"}
    for N in (3, 4, 5, 6):
        c.plot(np.abs(phi), 1 - np.cos(phi / N), color=STATE_COLORS[N - 3], lw=0.9, ls="--")
        sel = [r for r in rot if r["N"] == N]
        refl = [r for r in sel if r["holonomy"]["holonomy"] == "reflection"]
        rotn = [r for r in sel if r["holonomy"]["holonomy"] == "rotation"]
        c.plot([abs(r["holonomy"]["flux"]) for r in rotn], [r["measured"] for r in rotn], markers[N], ms=5,
               color=STATE_COLORS[N - 3], label=f"N = {N}")
        c.plot([0.0] * len(refl), [r["measured"] for r in refl], markers[N], ms=5, markerfacecolor="white",
               color=STATE_COLORS[N - 3])
    c.set_xlabel("loop mismatch |Phi| (from the bond directions)")
    c.set_ylabel("excess energy per bond of the ground state")
    c.text(0.03, 0.95, "dashed: 1 - cos(Phi/N), calculated from the bond directions\nopen symbols at 0: odd number of reflections",
           transform=c.transAxes, fontsize=6.8, va="top", color=INK2)
    c.legend(fontsize=6.5, loc="center left", frameon=False)
    _style(c)
    _title(c, "c", "Rotor loops: excess energy against the loop mismatch")
    _example_loops(d, rot)
    _title(d, "d", "Ground states of three rotor loops")
    ok = sum(r["consistent"] for r in rows)
    fig.suptitle(f"Loops of genes, spins and rotors: structural predictions and calculated states "
                 f"({ok} of {len(rows)} loops consistent)", fontsize=10, color=INK, x=0.01, ha="left", fontweight="bold")
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ co-discovery
SLOTS_BY_TARGET = {"symmetric-write": ["S", "C", "R1", "U", "R2", "K", "L"],
                   "threshold-write": ["S", "C", "R1", "W", "R2", "K", "L"],
                   "phase-locking": ["S", "C", "R1", "K", "L"]}
SLOT_HEADS = {"S": "symmetry", "C": "continuation", "R1": "reduction", "U": "unfolding\nto the cusp",
              "W": "write field", "R2": "reduction", "K": "canonical\nform", "L": "law"}
SLOT_HEADS_BY_TARGET = {"phase-locking": {"C": "into the\noscillation", "R1": "reduction\nto the phase",
                                          "K": "averaging:\nAdler form"}}
LETTER_COLOR = {"S": BLUE, "C": INK2, "R": GREEN, "U": ORANGE, "W": ORANGE, "K": PURPLE, "L": CYAN}
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "*"]


def _tile_note(st: Dict, key: str) -> str:
    L = st["letter"]
    if L == "S":
        m = str(st.get("meaning", ""))
        return f"cyclic {st['order']}" if m.startswith("cyclic") else m.split(" ")[0]
    if L == "C":
        return str(st.get("control", "")) + (f" = {st['value']:.3g}" if "value" in st else "")
    if L == "R" and key == "phase-locking":
        return "neutral" if "neutral" in str(st.get("kind", "")) else "phase"
    if L == "R":
        kind = str(st.get("kind", ""))
        return ("fold" if "fold" in kind else "Hopf" if "Hopf" in kind else "transcrit." if "transcritical" in kind
                else "subcrit." if "subcritical" in kind else "pitchfork" if "pitchfork" in kind else kind[:9])
    if L == "U":
        return str(st.get("second", ""))
    if L == "W":
        return "" if st.get("to_state") is None else f"to state {st['to_state']}"
    if L == "K":
        if key == "phase-locking":
            return f"{st.get('ratio', 1)}:1"
        if key == "threshold-write":
            return "mu + x^2"
        return "-x^3" if st.get("cubic_sign", -1) < 0 else "+x^3"
    if L == "L":
        c = st.get("law_constant")
        return "" if c is None else (f"{c:.4f}" if key == "phase-locking" else f"{c:.4g}")
    return ""


def _derivation_panel(ax, rows: List[Dict], class_color: Dict[str, str], key: str = "symmetric-write"):
    from matplotlib.patches import FancyBboxPatch
    order = SLOTS_BY_TARGET[key]
    heads = {**SLOT_HEADS, **SLOT_HEADS_BY_TARGET.get(key, {})}
    xpos = {slot: i for i, slot in enumerate(order)}
    for slot, x in xpos.items():
        ax.text(x, 0.75, f"{slot[0]}\n{heads[slot]}", ha="center", va="bottom", fontsize=6.4, color=INK2,
                linespacing=1.0)
    for i, r in enumerate(rows):
        y = -i
        xs, after_uw = [], False
        for st in r["steps"]:
            L = st["letter"]
            if L == "R":
                slot = "R2" if after_uw else "R1"
            else:
                slot = L
                after_uw = after_uw or L in ("U", "W")
            x = xpos[slot]
            xs.append(x)
            col = LETTER_COLOR[L]
            ax.add_patch(FancyBboxPatch((x - 0.36, y - 0.3), 0.72, 0.6, boxstyle="round,pad=0.02,rounding_size=0.08",
                                        facecolor=col, alpha=0.14, edgecolor="none"))
            ax.add_patch(FancyBboxPatch((x - 0.36, y - 0.3), 0.72, 0.6, boxstyle="round,pad=0.02,rounding_size=0.08",
                                        facecolor="none", edgecolor=col, lw=0.8))
            ax.text(x, y + 0.08, L, ha="center", va="center", fontsize=8, fontweight="bold", color=col)
            ax.text(x, y - 0.17, _tile_note(st, key), ha="center", va="center", fontsize=5.2, color=INK)
        if len(xs) > 1:
            ax.plot(xs, [y] * len(xs), color=GRID, lw=1.2, zorder=0)
        name_col = class_color.get(r.get("class"), INK2) if r["status"] != "obstructed" else INK2
        ax.text(-0.6, y + 0.1, r["name"], ha="right", va="center", fontsize=7, color=INK)
        ax.text(-0.6, y - 0.2, r["field"], ha="right", va="center", fontsize=6, color=name_col)
        end = (max(xs) if xs else -0.4) + 0.55
        x_text = len(order) - 0.25
        if r["status"] == "obstructed":
            ax.text(end, y, "x", ha="center", va="center", fontsize=8, fontweight="bold", color=RED)
            ax.text(x_text, y, r["obstruction_short"], ha="left", va="center", fontsize=6.4, color=RED)
        else:
            ax.text(x_text, y, r["status"], ha="left", va="center", fontsize=6.4, color=GREEN)
    ax.set_xlim(-4.6, len(order) + 3.6)
    ax.set_ylim(-len(rows) + 0.4, 1.55)
    ax.axis("off")


def _law_panel_pitchfork(c, reached, class_color, s):
    from .codiscovery import LAW_CONSTANT
    with_law = [r for r in reached if r.get("law_constant")]
    if not with_law:
        c.axis("off")
        _note(c, "the law was not simulated (--no-law)")
        return
    for j, r in enumerate(with_law):
        k = r["law_constant"]
        c.errorbar(k["constant"], -j, xerr=k["stderr"], fmt="o", ms=4, color=class_color[r["class"]], capsize=2, lw=1.0)
    c.set_yticks([-j for j in range(len(with_law))])
    c.set_yticklabels([r["name"] for r in with_law], fontsize=6.5)
    c.axvline(LAW_CONSTANT, color=INK, lw=1.0, ls="--")
    m, e = s["law_constant_mean"], s["law_constant_stderr"]
    c.axvspan(m - e, m + e, color=GRID, alpha=0.8, lw=0)
    c.text(LAW_CONSTANT, 0.6, " pi^(1/4)", fontsize=6.5, color=INK, va="bottom")
    c.set_ylim(-len(with_law) + 0.4, 1.0)
    c.set_xlabel("Phi^-1(P) D_s^(1/2) r^(1/4) / h_s from the measured accuracy P")
    c.text(0.02, 0.03, f"weighted mean {m:.3f} +/- {e:.3f}; pi^(1/4) = {LAW_CONSTANT:.3f}\n"
                       f"chi^2 = {s['law_constant_chi2']:.1f} for {s['law_constant_dof']} values; "
                       f"{s['trajectories_per_realization']} trajectories each",
           transform=c.transAxes, fontsize=6.2, color=INK2, va="bottom")
    _style(c)


def _law_panel_fold(c, reached, class_color, marker, s):
    from .codiscovery import DELAY_CONSTANT
    with_law = [r for r in reached if r.get("law_constant")]
    if not with_law:
        c.axis("off")
        _note(c, "the law was not simulated (--no-law)")
        return
    emax = 0.0
    for r in with_law:
        rows = r["law"]["rows"]
        e = np.array([row["rate"] for row in rows]) ** (1 / 3)
        v = np.array([row["constant"] for row in rows])
        emax = max(emax, float(e.max()))
        col = class_color[r["class"]]
        c.plot(e, v, marker[r["name"]], ms=3.2, color=col, lw=0, label=r["name"])
        f0, f1, f2 = r["law"]["fit"]
        ee = np.linspace(0.0, float(np.sort(e)[4]), 60)
        c.plot(ee, f0 + f1 * ee + f2 * ee ** 2, color=col, lw=0.7, alpha=0.8)
    c.plot([0.0], [DELAY_CONSTANT], "*", ms=9, color=INK, zorder=5)
    c.axhline(DELAY_CONSTANT, color=INK, lw=0.8, ls="--")
    c.text(0.76 * emax, DELAY_CONSTANT, "|a1'| = 1.01879", fontsize=6.5, color=INK, va="bottom")
    c.set_xlim(-0.003, emax * 1.03)
    c.set_xlabel("r^(1/3) (canonical sweep rate r)")
    c.set_ylabel("mu at the crossing of the static fold / r^(2/3)")
    c.text(0.02, 0.97, f"extrapolated to r = 0: {s['law_constant_mean']:.5f} +/- {s['law_constant_stderr']:.5f}\n"
                       f"largest deviation from |a1'|: {s['law_constant_max_deviation']:.1e}\n"
                       f"markers as in (b); lines: quadratic fits in r^(1/3) to the five slowest sweeps",
           transform=c.transAxes, fontsize=6.2, color=INK2, va="top")
    _style(c)


def codiscovery_figure(report: Dict, path: Path):
    """(a) the derivation of the target in every realization, as letters on the slots of a common chain;
    (b) the reduced drift at the write point in the canonical variable, against -x^3 (symmetric write) or x^2
    (threshold write, at the fold); (c) the constant of the law in each realization: the swept-write constant
    against pi^(1/4), or the delay of the switch against the sweep rate, extrapolated to |a1'|."""
    rows = report["rows"]
    s = report["summary"]
    key = s.get("target_key", "symmetric-write")
    if key == "phase-locking":
        return _phase_locking_figure(report, path)
    classes = list(s["derivation_classes"])
    class_color = {c: STATE_COLORS[i % len(STATE_COLORS)] for i, c in enumerate(classes)}
    reached = [r for c in classes for r in rows if r.get("class") == c and r["status"] != "obstructed"]
    marker = {}
    for c in classes:
        for i, r in enumerate([r for r in reached if r["class"] == c]):
            marker[r["name"]] = MARKERS[i % len(MARKERS)]
    ordered = reached + [r for r in rows if r["status"] == "obstructed"]
    fig = plt.figure(figsize=(12.0, 9.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0])
    a = fig.add_subplot(gs[0, :])
    _derivation_panel(a, ordered, class_color, key)
    title = "symmetric write" if key == "symmetric-write" else "threshold write"
    _title(a, "a", f"Derivation of the {title} in each realization")
    b = fig.add_subplot(gs[1, 0])
    if key == "symmetric-write":
        xm = max(max(abs(v) for v in r["canonical"]["x"]) for r in reached) if reached else 1.0
        xx = np.linspace(-xm, xm, 200)
        b.plot(xx, -xx ** 3, color=INK, lw=1.0, ls="--", label="-x^3", zorder=3)
        b.set_ylabel("g = |a3|^(1/2) ds/dt at the write point")
        b.set_xlabel("x = |a3|^(1/2) s (canonical coordinate)")
        note = (f"even part <= {s['even_part_max']:.0e} of the cubic term\ndeviations at large |x|: fifth-order terms, "
                f"specific to each realization")
    else:
        xm = 1.0
        xx = np.linspace(-xm, xm, 200)
        b.plot(xx, xx ** 2, color=INK, lw=1.0, ls="--", label="x^2", zorder=3)
        b.set_ylabel("g = a2 ds/dt at the fold")
        b.set_xlabel("x = a2 s (canonical coordinate)")
        b.set_xlim(-xm, xm)
        b.set_ylim(-0.25, 1.6)
        note = (f"linear part <= {s['linear_part_max']:.0e} of the slowest transverse rate\ndeviations at large |x|: "
                f"cubic terms, specific to each realization")
    for r in reached:
        cf = r["canonical"]
        x, g = np.asarray(cf["x"]), np.asarray(cf["g"])
        keep = np.abs(x) <= xm
        b.plot(x[keep], g[keep], marker[r["name"]], ms=2.4, color=class_color[r["class"]], alpha=0.85,
               label=f"{r['name']} ({r['field']})")
    b.legend(fontsize=5.4, loc="upper right" if key == "symmetric-write" else "upper center", frameon=False)
    b.text(0.02, 0.03, note, transform=b.transAxes, fontsize=6.2, color=INK2, va="bottom")
    _style(b)
    _title(b, "b", "Reduced drift at the write point in canonical units")
    c = fig.add_subplot(gs[1, 1])
    if key == "symmetric-write":
        _law_panel_pitchfork(c, reached, class_color, s)
        _title(c, "c", "Constant of the swept-write law in each realization")
    else:
        _law_panel_fold(c, reached, class_color, marker, s)
        _title(c, "c", "Delay of the switch against the sweep rate")
    fig.suptitle(f"Co-discovery by construction: {s['reached']} realizations from {len(s['fields_reached'])} fields "
                 f"reach the {title} by {len(classes)} classes of derivation", fontsize=10, color=INK,
                 x=0.01, ha="left", fontweight="bold")
    fig.savefig(path, dpi=130)
    plt.close(fig)


def _phase_locking_figure(report: Dict, path: Path):
    """(a) derivations of phase locking; (b) the phase response along the drive, different in every realization;
    (c) the phase gained per period in the driven realization, in canonical units, against -sin phi; (d) relaxation
    rate inside and slip frequency outside the locking range against the detuning, in units of K."""
    rows = report["rows"]
    s = report["summary"]
    classes = list(s["derivation_classes"])
    class_color = {c: STATE_COLORS[i % len(STATE_COLORS)] for i, c in enumerate(classes)}
    reached = [r for c in classes for r in rows if r.get("class") == c and r["status"] != "obstructed"]
    marker = {r["name"]: MARKERS[i % len(MARKERS)] for i, r in enumerate(reached)}
    ordered = reached + [r for r in rows if r["status"] == "obstructed"]
    fig = plt.figure(figsize=(13.0, 5.4 + 0.3 * len(ordered)), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[0.1 + 0.105 * len(ordered), 1.0])
    a = fig.add_subplot(gs[0, :])
    _derivation_panel(a, ordered, class_color, "phase-locking")
    _title(a, "a", "Derivation of phase locking in each realization")
    b = fig.add_subplot(gs[1, 0])
    for r in reached:
        c = r["canonical"]
        th = np.asarray(c["theta"]) / (2 * np.pi)
        z = np.asarray(c["Zp"]) - np.mean(c["Zp"])
        b.plot(th, z / max(np.max(np.abs(z)), 1e-300), "-", lw=1.0, color=class_color[r["class"]],
               marker=marker[r["name"]], ms=2.4, markevery=6, label=f"{r['name']} ({c['ratio']}:1)")
    b.set_xlabel("phase of the cycle theta / 2 pi")
    b.set_ylabel("Z_p - mean, divided by its largest value")
    b.legend(fontsize=5.2, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False)
    _style(b)
    _title(b, "b", "Phase response along the drive")
    c_ax = fig.add_subplot(gs[1, 1])
    ph = np.linspace(-np.pi, np.pi, 300)
    c_ax.plot(ph, -np.sin(ph), color=INK, lw=1.0, ls="--", zorder=3, label="-sin phi")
    for r in reached:
        c = r["canonical"]
        c_ax.plot(c["phi"], c["g"], marker[r["name"]], ms=3.2, color=class_color[r["class"]], lw=0)
    c_ax.set_xlabel("phi = n psi - phi_n - pi/2")
    c_ax.set_ylabel("(phase gained per period) / (T K)")
    dmax = s.get("deviation_max")
    c_ax.text(0.98, 0.97, f"largest deviation {dmax:.1e} of K\nat K = 0.01 min(omega_0, 2 kappa),\n"
                          f"first order in the drive;\nmarkers as in (b)" if dmax is not None else "",
              transform=c_ax.transAxes, fontsize=6.0, color=INK2, va="top", ha="right")
    _style(c_ax)
    _title(c_ax, "c", "Averaged drift of the phase difference")
    d = fig.add_subplot(gs[1, 2])
    with_law = [r for r in reached if r.get("law_constant")]
    if not with_law:
        d.axis("off")
        _note(d, "the law was not simulated (--no-law)")
    else:
        u = np.linspace(-1, 1, 200)
        d.plot(u, np.sqrt(1 - u ** 2), color=INK, lw=1.0, ls="--")
        for sgn in (-1, 1):
            v = sgn * np.linspace(1, 2.2, 100)
            d.plot(v, np.sqrt(v ** 2 - 1), color=INK, lw=1.0, ls="--")
        for r in with_law:
            fit = min(r["law"]["fits"], key=lambda f: f["eps"])
            col = class_color[r["class"]]
            for row in fit["rows"]:
                face = col if row["side"] == "inside" else "white"
                d.plot(row["nu"], row["canonical"], marker[r["name"]], ms=3.6, color=col, markerfacecolor=face, lw=0)
        d.text(0.0, 1.13, "relaxation rate\nlambda / (n K)", ha="center", va="bottom", fontsize=6.2, color=INK2)
        d.text(1.72, 0.55, "slip frequency\nOmega / K", ha="center", va="center", fontsize=6.2, color=INK2)
        d.set_xlim(-2.2, 2.2)
        d.text(0.5, 0.98, f"half-width of the locking range, extrapolated to zero drive:\n"
                          f"{s['law_constant_mean']:.4f} +/- {s['law_constant_stderr']:.4f} K "
                          f"(largest deviation from 1: {s['law_constant_max_deviation']:.1e})\n"
                          f"points at K = 0.01 min(omega_0, 2 kappa); filled: inside, open: outside",
               transform=d.transAxes, fontsize=5.8, color=INK2, va="top", ha="center")
        d.set_ylim(0, 2.45)
        d.set_xlabel("nu = Delta omega / K (detuning in units of K)")
        d.set_ylabel("rate or frequency in units of K")
        _style(d)
    _title(d, "d", "Relaxation and slips against the detuning")
    fig.suptitle(f"Co-discovery by construction: {s['reached']} realizations from {len(s['fields_reached'])} fields "
                 f"reach phase locking by {len(classes)} classes of derivation", fontsize=10, color=INK,
                 x=0.01, ha="left", fontweight="bold")
    fig.savefig(path, dpi=130)
    plt.close(fig)


def _grid_panel(ax, rows, unit):
    if not rows:
        ax.axis("off")
        return
    Ns = sorted({r["N"] for r in rows})
    rmax = max(Ns)
    code = np.full((len(Ns), rmax + 1), np.nan)
    for r in rows:
        k = r["holonomy"]
        nneg = sum(1 for x in r["signs"] if x < 0)
        i = Ns.index(r["N"])
        if unit == "gene":
            val = 0 if r["measured"] == "oscillation" else (1 if r["stable_states"] == 1 else 2)
            txt = {0: "osc", 1: "1", 2: str(r["stable_states"])}[val]
        else:
            val = 3 if r["measured"] == "frustrated" else 2
            txt = f"{r['unsatisfied_bonds']} / {r['ground_degeneracy']}"
        code[i, nneg] = val
        ax.text(nneg, i, txt, ha="center", va="center", fontsize=7.5, color="white" if val in (0, 3) else INK)
        if not k["satisfiable"]:
            ax.add_patch(plt.Rectangle((nneg - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK, lw=1.4, ls="-"))
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap([ORANGE, "#d1d5db", "#93c5fd", RED])
    ax.imshow(code, cmap=cmap, vmin=-0.5, vmax=3.5, aspect="auto")
    ax.set_yticks(range(len(Ns))); ax.set_yticklabels([f"N = {n}" for n in Ns], fontsize=7)
    ax.set_xticks(range(rmax + 1))
    cells = ("cells: stable states found, osc = oscillation" if unit == "gene"
             else "cells: unsatisfied bonds in the ground state / ground states found")
    ax.set_xlabel(("repressions" if unit == "gene" else "antiferromagnetic bonds") + " in the loop  (boxed: negative loop)\n"
                  + cells, fontsize=7.5)
    ax.tick_params(labelsize=7)


def _example_loops(ax, rot):
    picks = [r for r in rot if r.get("showcase")]
    colors = [GREEN, BLUE, RED]
    for k, r in enumerate(picks[:3]):
        pos = np.asarray(r["geometry"]["pos"])
        pos = pos - pos.mean(axis=0)
        pos = pos / max(np.abs(pos).max(), 1e-9) + np.array([2.6 * k, 0.0])
        for i, j in zip(r["geometry"]["src"], r["geometry"]["tgt"]):
            ax.plot(pos[[i, j], 0], pos[[i, j], 1], color=GRID, lw=1.6, zorder=0)
        ax.plot(pos[:, 0], pos[:, 1], "o", color=INK2, ms=3, zorder=1)
        _sticks(ax, pos, np.asarray(r["ground_state"]), 2, colors[k], 2.6, 0.32)
        hol = r["holonomy"]
        lab = r["showcase"] + (f"\n|Phi| = {abs(hol['flux']):.2f}" if hol["holonomy"] == "rotation" else "\nholonomy: a reflection")
        ax.text(pos[:, 0].mean(), -1.55, f"{lab}\nexcess energy per bond {r['measured']:.3f}", ha="center", va="top",
                fontsize=7, color=INK)
    ax.set_xlim(-1.4, 2.6 * max(len(picks), 1) - 1.2)
    ax.set_ylim(-2.6, 1.4)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[:].set_visible(False)


# ------------------------------------------------------------------------------------------------ HTML
def _img64(path: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(Path(path).read_bytes()).decode()


CSS = """
:root { --ground:#f5f7fa; --panel:#ffffff; --ink:#19202b; --muted:#586273; --rule:#dfe4eb; --accent:#1d5fd1;
  --ok:#0e7a55; --warn:#a55a0a; --bad:#b3261e; --plate:#ffffff; --chip-ink:#ffffff; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { color-scheme: dark;
  --ground:#0f141b; --panel:#161c25; --ink:#e6eaf0; --muted:#9aa4b2; --rule:#283141; --accent:#7fb0ff;
  --ok:#3fbf8f; --warn:#e0a458; --bad:#f07a70; --plate:#f3f5f8; --chip-ink:#0f141b; } }
:root[data-theme="dark"] { color-scheme: dark; --ground:#0f141b; --panel:#161c25; --ink:#e6eaf0; --muted:#9aa4b2;
  --rule:#283141; --accent:#7fb0ff; --ok:#3fbf8f; --warn:#e0a458; --bad:#f07a70; --plate:#f3f5f8; --chip-ink:#0f141b; }
body { background:var(--ground); color:var(--ink); font:15px/1.6 "Public Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  padding-inline:16px; padding-block:28px 72px; }
main { max-width:1120px; margin:0 auto; display:grid; gap:28px; }
h1, h2 { font-family:"Spectral", Georgia, "Times New Roman", serif; font-weight:600; text-wrap:balance; margin:0; }
h1 { font-size:clamp(28px, 4vw, 38px); line-height:1.15; }
h2 { font-size:22px; }
h3 { font-size:16px; margin:0; }
p { margin:0; max-width:72ch; }
.lead { color:var(--muted); }
.mono, code { font-family:"JetBrains Mono", ui-monospace, "SFMono-Regular", Menlo, monospace; font-size:0.86em; }
.stats { display:flex; flex-wrap:wrap; gap:8px 22px; color:var(--muted); font-size:14px; }
.stats b { color:var(--ink); font-variant-numeric:tabular-nums; }
.steps { display:grid; grid-template-columns:repeat(auto-fit, minmax(230px, 1fr)); gap:12px; }
.step { border-top:2px solid var(--accent); padding-top:8px; font-size:14px; color:var(--muted); }
.step b { color:var(--ink); display:block; font-size:14px; letter-spacing:0.02em; }
.scroll { overflow-x:auto; border:1px solid var(--rule); border-radius:8px; background:var(--panel); }
table { border-collapse:collapse; width:100%; font-size:13.5px; min-width:760px; }
th, td { text-align:left; vertical-align:top; padding:8px 10px; border-bottom:1px solid var(--rule); }
th { color:var(--muted); font-weight:600; font-size:12px; text-transform:uppercase; letter-spacing:0.06em; }
tr:last-child td { border-bottom:0; }
td a { color:var(--accent); text-decoration:none; }
td a:hover, td a:focus-visible { text-decoration:underline; }
.chip { display:inline-block; padding:1px 8px; border-radius:999px; font-size:12px; font-weight:600; color:var(--chip-ink); white-space:nowrap; }
.chip.ok { background:var(--ok); } .chip.warn { background:var(--warn); } .chip.bad { background:var(--bad); }
.chip.none { background:var(--muted); }
.chip.info { background:var(--accent); }
section.card { display:grid; gap:12px; padding-top:18px; border-top:1px solid var(--rule); }
.card-head { display:flex; flex-wrap:wrap; gap:6px 12px; align-items:baseline; }
.card-head .src { color:var(--muted); font-size:14px; }
.plate { background:var(--plate); border-radius:8px; padding:8px; }
.plate img { display:block; width:100%; height:auto; max-width:100%; }
.facts { display:grid; grid-template-columns:repeat(auto-fit, minmax(300px, 1fr)); gap:14px 28px; }
dl { display:grid; grid-template-columns:max-content 1fr; gap:4px 14px; margin:0; font-size:14px; }
dt { color:var(--muted); }
dd { margin:0; }
nav { display:flex; flex-wrap:wrap; gap:6px; }
nav a { color:var(--accent); border:1px solid var(--rule); border-radius:6px; padding:2px 8px; font-size:13px; text-decoration:none; background:var(--panel); }
nav a:focus-visible, td a:focus-visible { outline:2px solid var(--accent); outline-offset:2px; }
@media (max-width:640px) { dl { grid-template-columns:1fr; } dt { margin-top:6px; } }
"""

FONTS = ('<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>'
         '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&'
         'family=Public+Sans:wght@400;600&family=Spectral:wght@600&display=swap">')


VERDICT_NAMES = {"stores": "Stable states", "writes": "Write point", "constructs": "Writing mechanism",
                 "holds": "Retention law", "lock": "Writing protocols"}
SLOT_NAMES = {"Omega": "Operation (Omega)", "Xi": "Carrier (Xi)", "C": "Closure (C)", "R": "Observable (R)",
              "P": "Protocol (P)", "A": "Material parameters (A)"}


def _mechanism(card) -> str:
    """A descriptive classification of how a state is written; none of the classes is an error."""
    from .discovery import mechanism_class
    return mechanism_class(card)


def _chip(card) -> str:
    return f'<span class="chip info">{html.escape(_mechanism(card))}</span>'


def _short(text: str, n: int = 110) -> str:
    return text if len(text) <= n else text[: n - 1].rsplit(" ", 1)[0] + "…"


def html_report(cards: List[Dict], path: Path, loops: Optional[List[Dict]] = None,
                title: str = "Memory in Model Materials",
                heading: str = "Memory in model materials: a gallery of realizations"):
    """Writes the gallery twice: gallery.html (page content only, the form published as an Artifact) and
    index.html (the same content in a full document, for opening locally). `title` names the page; `heading` is
    its first heading."""
    path = Path(path)
    esc = html.escape
    head = f"<title>{esc(title)}</title>{FONTS}<style>{CSS}</style>"
    parts = ["<main>"]
    n_loops = len(loops) if loops else 0
    ok_loops = sum(r["consistent"] for r in loops) if loops else 0
    counts = {}
    for c in cards:
        counts[_mechanism(c)] = counts.get(_mechanism(c), 0) + 1
    parts.append(f"<header style='display:grid;gap:12px'><h1>{esc(heading)}</h1>"
                 "<p class='lead'>Each entry is a model material, called a realization: a set of state variables (the "
                 "carrier) with deterministic dynamics and thermal noise. For every realization the program determines "
                 "the stable states, the parameter value at which a new state can be written (the write point), the "
                 "law by which a stored state is lost, and how the local dynamics at the write point differs from the "
                 "symmetric pitchfork normal form &epsilon;s &minus; s&sup3; + h. The classifications are descriptive: "
                 "a fold, an oscillation or the absence of a bifurcation is a property of the material, not a failed "
                 "calculation.</p>"
                 f"<div class='stats'><span><b>{len(cards)}</b> realizations</span>"
                 + "".join(f"<span><b>{n}</b> {esc(k)}</span>" for k, n in sorted(counts.items(), key=lambda kv: -kv[1]))
                 + (f"<span><b>{ok_loops}</b> of <b>{n_loops}</b> loops as predicted from the holonomy</span>"
                    if loops else "") + "</div></header>")
    parts.append("<div class='steps'><div class='step'><b>Carrier and closure</b>The variables of the material and "
                 "what the dynamics conserves or exchanges with its environment.</div><div class='step'><b>Write "
                 "point</b>The parameter value at which a stable state loses stability, located by continuation along "
                 "the control parameter or along a uniform field.</div><div class='step'><b>Normal form</b>The drift "
                 "along the unstable direction at the write point. A quadratic term, a bias or a positive cubic "
                 "coefficient distinguishes it from the symmetric pitchfork; tuning a second parameter to a cusp "
                 "restores the symmetric case.</div></div>")
    rows = []
    for c in cards:
        v = c["verdict"]
        rows.append(f"<tr><td><a href='#{c['id']}'>{esc(_short(c['name'], 60))}</a></td><td>{_chip(c)}</td>"
                    f"<td>{esc(_short(v['stores'], 60))}</td><td>{esc(_short(v['writes'], 90))}</td>"
                    f"<td>{esc(_short(v['holds'], 60))}</td></tr>")
    parts.append("<section style='display:grid;gap:10px'><h2>At a glance</h2><div class='scroll'><table><thead><tr>"
                 "<th>Realization</th><th>Writing mechanism</th><th>Stable states</th><th>Write point</th>"
                 "<th>Retention law</th></tr>"
                 f"</thead><tbody>{''.join(rows)}</tbody></table></div></section>")
    if loops:
        img = path.parent / "loops.png"
        parts.append("<section class='card' id='loops'><div class='card-head'><h2>Loops across fields</h2></div>"
                     f"<p>For {ok_loops} of {n_loops} loops, the behaviour predicted from the composition of the bond "
                     "couplings around the loop (its holonomy) agrees with the calculated states. Gene loops follow "
                     "Thomas's rule (a negative feedback loop cannot sustain two stable states), loops of bistable spins "
                     "follow Toulouse's frustration criterion, and loops of rotors are frustrated by the mismatch "
                     "&Phi; of their bond directions, with minimum excess energy per bond 1 &minus; cos(&Phi;/N).</p>"
                     + (f"<div class='plate'><img alt='Loops of genes, spins and rotors: prediction and simulation' "
                        f"src='{_img64(img)}'></div>" if img.exists() else "") + "</section>")
    parts.append("<nav aria-label='Realizations'>" + "".join(
        f"<a href='#{c['id']}'>{esc(_short(c['name'].split(',')[0], 34))}</a>" for c in cards) + "</nav>")
    for c in cards:
        img = path.parent / f"{c['id']}.png"
        v = c["verdict"]
        ver = "".join(f"<dt>{esc(VERDICT_NAMES.get(k, k))}</dt><dd>{esc(str(val))}</dd>" for k, val in v.items())
        slots = "".join(f"<dt>{esc(SLOT_NAMES.get(k, k))}</dt><dd>{esc(str(val))}</dd>" for k, val in c["slots"].items())
        parts.append(f"<section class='card' id='{c['id']}'><div class='card-head'><h3>{esc(c['name'])}</h3>{_chip(c)}"
                     f"<span class='src'>{esc(c['provenance'])}</span></div>"
                     + (f"<div class='plate'><img alt='{esc(c['name'])}' src='{_img64(img)}'></div>" if img.exists() else "")
                     + f"<div class='facts'><dl>{ver}</dl><dl>{slots}</dl></div></section>")
    parts.append("</main>")
    body = "".join(parts)
    (path.parent / "gallery.html").write_text(head + body, encoding="utf-8")
    path.write_text("<!doctype html><html lang='en'><head><meta charset='utf-8'>"
                    "<meta name='viewport' content='width=device-width, initial-scale=1, viewport-fit=cover'>"
                    + head + "</head><body>" + body + "</body></html>", encoding="utf-8")
