"""Automatic evaluation of realizations: the states they store, where and how the states are written, how long
they are retained, how long a rewrite takes, and what a transfer to another carrier would require.

evaluate(real) returns a card:
  states       stable states from random starts, basin counts, kappa spectra, capacity (bits)
  loss_law     the law that loses stored information (from the spectra and the number of states)
  construct    write points along the control and along a write field, their normal forms, the obstruction
               to the pitchfork generator, its cancellation by a further material parameter, and the transferred
               swept-write law checked by simulation (construct.py)
  bifurcation  fixed points along the control, for the diagram
  barrier      minimum-energy-path barrier between the two most populated states (gradient systems)
  hold         Monte Carlo retention time at the realization's noise (censored at the simulated interval)
  writes       bounded writes below and above the restoring threshold at fixed control, and below it while the
               control is swept, from the end of its range where the barriers are low; accuracy, rewriting time,
               and the ratio of retention time to rewriting time (a lower bound when the retention is censored)
  kernel       memory kernel of an unobserved linear block, if the realization declares one
  transfers    intertwining defects of construction maps whose source is this realization
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from . import analysis as an
from .construct import construct
from .identity import Realization
from .transfer import Transfer


def _cplx(k: np.ndarray) -> List[List[float]]:
    return [[float(x.real), float(x.imag)] for x in k]


def _sweep_start(real: Realization, scan: List[Dict]) -> Optional[float]:
    """The end of the control range where the landscape retains least: fewer stable states, or weaker
    restoring rates (for a control that only scales the landscape, as a coupling or inverse temperature)."""
    if not scan:
        return None
    ends = [scan[0], scan[-1]]
    n_st = [sum(p["stable"] for p in e["points"]) for e in ends]
    if n_st[0] != n_st[1]:
        return ends[int(np.argmin(n_st))]["value"]
    k = [max((p["kappa_min"] for p in e["points"] if p["stable"]), default=np.inf) for e in ends]
    lo, hi = min(k), max(k)
    if np.isfinite(lo) and lo < 0.25 * hi:
        return ends[int(np.argmin(k))]["value"]
    return None


def evaluate(real: Realization, rng: np.random.Generator, transfers: Optional[List[Transfer]] = None,
             quick: bool = False) -> Dict:
    t0 = time.time()
    card: Dict = {"name": real.name, "slots": real.slots(), "provenance": real.provenance, "tags": list(real.tags),
                  "carrier": {"kind": real.carrier.kind, "dim": real.carrier.dim, "period": real.carrier.period,
                              "scale": real.carrier.scale},
                  "control": real.control, "control_range": list(real.control_range) if real.control_range else None,
                  "params": dict(real.params), "geometry": real.geometry}
    n_starts = 48
    states, counts, spectra, n_unconv, continuum = an.stored_states(real, rng, n_starts)
    n_null = 0 if real.null_modes is None else len(np.atleast_2d(real.null_modes))
    p = np.array(counts, float) / max(sum(counts), 1)
    card["states"] = {"count": len(states), "basin_counts": counts, "unconverged": n_unconv, "continuum": continuum,
                      "kappa": [_cplx(k) for k in spectra], "classes": [an.classify(k) for k in spectra],
                      "capacity_bits": float(np.log2(len(states))) if states else 0.0,
                      "basin_entropy_bits": float(-(p * np.log2(p)).sum()) if len(p) else 0.0,
                      "sampled_from": n_starts, "q": [s.tolist() for s in states]}
    card["loss_law"] = an.loss_law(len(states), spectra, n_unconv, n_starts, continuum, n_null)
    if real.carrier.dim == 1 or real.carrier.dim == 2:
        card["landscape"] = _landscape(real)
    stored = states if (len(states) >= 2 and not continuum) else None
    card["construct"] = construct(real, rng, states=stored, check_write=not quick)
    scan = []
    scale = "control_role" in card["construct"]
    if real.control and real.carrier.dim <= 12 and not scale:
        values = np.linspace(*real.control_range, 9 if real.carrier.dim > 4 else 13)
        scan = an.bifurcation_scan(real, rng, values, 12 if real.carrier.dim > 4 else 24, seeds=states)
        card["bifurcation"] = scan
    if not states and n_unconv > 0:
        q0 = real.carrier.sample(rng, 1)
        q, t, traj = an.integrate(real, q0, 200.0, 0.0, rng, record=lambda x: x[0].copy(), n_record=400)
        card["trajectory"] = {"t": t[200:].tolist(), "q": np.asarray(traj[200:]).tolist()}
    if stored is not None:
        a, b = states[0], states[1]
        if real.potential is not None:
            card["barrier"] = an.neb_barrier(real, a, b, iters=1500 if real.carrier.dim > 4 else 3000)
        wf = card["construct"].get("write_field", {})
        thr = float(wf["h_c"]) if "h_c" in wf else an.restoring_threshold(real, a, b)
        card["threshold"] = thr
        card["threshold_source"] = ("field at which the write-field continuation loses the stored state" if "h_c" in wf
                                    else "largest restoring force on the straight path")
        kmin = max(min(float(k.real[k.real > 1e-3].min()) for k in spectra), 1e-2)
        T_w, T_rel = float(np.clip(30 / kmin, 5, 150)), float(np.clip(10 / kmin, 2, 50))
        t_max = float(min(300.0, 30000 * real.dt))
        n_traj = (60 if quick else 150) if real.carrier.dim <= 4 else (40 if quick else 80)
        card["hold"] = an.hold_time(real, 0, states, real.noise, rng, n_traj=n_traj, t_max=t_max)
        writes = [an.write_test(real, 0, 1, states, f * thr, T_w, T_rel, real.noise, rng, n_traj=n_traj)
                  for f in (0.4, 1.5)]
        start = real.control_range[0] if scale else _sweep_start(real, scan)
        if start is not None:
            sweep = an.write_test(real, 0, 1, states, 0.4 * thr, T_w, T_rel, real.noise, rng, n_traj=n_traj,
                                  mode="sweep", sweep_from=start)
            sweep["sweep_from"], sweep["sweep_to"] = float(start), float(real.params[real.control])
            writes.append(sweep)
        for w in writes:
            w["h_over_threshold"] = w["h"] / thr
            w["lock_ratio"] = (card["hold"]["time"] / w["rewrite_time"]) if np.isfinite(w["rewrite_time"]) else None
            # a retention time censored by the simulated interval makes the ratio a lower bound
            w["lock_ratio_is_lower_bound"] = bool(card["hold"].get("censored")) and w["lock_ratio"] is not None
        card["writes"] = writes
    if getattr(real, "linear", None) is not None:
        tg = np.linspace(0, 2.0, 41)
        card["kernel"] = {"t": tg.tolist(), "K": an.memory_kernel(real, tg).tolist(),
                          "blocks": {k: np.asarray(v).tolist() for k, v in real.linear.items()}}
    if transfers:
        card["transfers"] = [{"name": tr.name, "target": tr.target.name, "note": tr.note, **tr.defect(rng)}
                             for tr in transfers]
    card["verdict"] = verdict(card)
    card["seconds"] = time.time() - t0
    return card


def _landscape(real: Realization) -> Dict:
    """Drift (and potential) on a grid, for the drawing of one- and two-dimensional carriers."""
    c = real.carrier
    lo = 0.0 if c.kind in ("orthant", "torus") else -c.scale
    hi = c.period if c.kind == "torus" else (1.1 * c.scale)
    if c.dim == 1:
        x = np.linspace(lo, hi, 400)[:, None]
        out = {"x": x[:, 0].tolist(), "F": real.F(x)[:, 0].tolist()}
        if real.potential is not None:
            out["V"] = real.V(x).tolist()
        return out
    g = np.linspace(lo, hi, 40)
    X, Y = np.meshgrid(g, g)
    F = real.F(np.stack([X.ravel(), Y.ravel()], axis=1))
    return {"grid": g.tolist(), "U": F[:, 0].reshape(X.shape).tolist(), "V": F[:, 1].reshape(X.shape).tolist()}


def sweep_label(card: Dict, w: Dict) -> str:
    """How the control moved during a writing protocol."""
    if not w["mode"].startswith("sw"):
        return "fixed control"
    if card["construct"].get("control_role") and "sweep_from" in w:
        return f"control raised from {w['sweep_from']:.3g} to {w['sweep_to']:.3g}"
    return "control swept through the write point"


def mechanism_class(card: Dict) -> str:
    """A descriptive class of how a state is written, from the calculated write points; no class is an error."""
    st = card["states"]
    if st["count"] == 0:
        return "limit cycle: phase memory" if "limit cycle" in card["loss_law"] else "no stable state"
    if st.get("continuum"):
        return "zero mode: diffusive retention"
    if st["count"] == 1:
        return "single state: relaxation"
    ctrl = [e for e in card["construct"]["events"] if e["source"] == "control"]
    obs = [e.get("obstruction", {}) for e in ctrl]
    if any(o.get("realizes_generator") for o in obs):
        return "symmetric write (supercritical pitchfork)"
    if any(e.get("cancellation") for e in ctrl):
        return "one-sided write (fold); symmetric at a cusp"
    if obs and all(o.get("subcritical") and not o.get("asymmetric") and not o.get("biased") for o in obs):
        return "symmetric write to a distant state (subcritical pitchfork)"
    if ctrl:
        return "asymmetric write"
    return "write by a uniform field"


def verdict(card: Dict) -> Dict[str, str]:
    """Short statements for the gallery: stable states, write point, writing mechanism, retention, writing protocols."""
    st = card["states"]
    if st["count"] == 0:
        stores = card["loss_law"] if card["loss_law"].startswith("no stable state") else (
            "no stable state: " + card["loss_law"])
    elif st["continuum"]:
        stores = "a continuum of neutral states (zero mode)"
    else:
        stores = f"{st['count']} stable state(s)" + (" (sampled)" if st["count"] >= 8 else "")
    ev = card["construct"]["events"]
    ctrl = [e for e in ev if e["source"] == "control"]
    field = [e for e in ev if e["source"] == "write field"]
    writes = "; ".join(sorted({f"{e['kind'].split(':')[0]} at {e['param']} = {e['value']:.3g}" for e in ctrl})) or (
        "no control parameter" if not card.get("control") else
        "no bifurcation along the control, which only rescales the energy landscape" if card["construct"].get("control_role")
        else "no bifurcation along the control")
    if field:
        writes += f"; a uniform field removes the stable state at h = {field[0]['value']:.3g}"
    gen = [e for e in ctrl if e.get("obstruction", {}).get("realizes_generator")]
    canc = [e["cancellation"] for e in ctrl if e.get("cancellation")]
    if gen:
        constructs = "supercritical pitchfork: symmetric write"
        law = gen[0].get("write_law")
        if law:
            constructs += (f"; accuracy of a swept write predicted {law['predicted']:.2f}, measured {law['measured']:.2f} "
                           f"(Gamma = {law['gamma']:.2g}); {law['at_material_noise']['measured']:.2f} at the material's "
                           f"own noise (Gamma = {law['at_material_noise']['gamma']:.2g})")
    elif canc:
        c = canc[0]
        constructs = (f"fold: one-sided write; a symmetric write is recovered at the cusp {c['second']} = "
                      f"{c['value']:.4g}, {card['control']} = {c['control']:.4g} "
                      f"({c['normal_form_along_sweep']['kind'].split(':')[0]} along the sweep)")
    elif ctrl and all(e["obstruction"].get("hopf") for e in ctrl):
        constructs = "Hopf bifurcation: the state oscillates and no stable state is selected; information is kept in the phase"
    elif ctrl and all(e["obstruction"].get("subcritical") and not e["obstruction"].get("asymmetric")
                      and not e["obstruction"].get("biased") for e in ctrl):
        constructs = ("subcritical pitchfork: symmetric write, but the state jumps to a distant branch; "
                      "a saturating term is needed for a local write")
    elif ctrl:
        constructs = ("differs from the pitchfork normal form by: "
                      + "; ".join(sorted({k for e in ctrl for k in e["obstruction"]["terms"]})))
    elif st["count"] >= 2:
        constructs = "no bifurcation along the control; writing requires a field"
    else:
        constructs = "a single stable state: no write point" + (
            "; the unobserved variables give the observed one a memory kernel" if card.get("kernel") else "")
    lock = "-"
    if card.get("writes"):
        parts = []
        for w in card["writes"]:
            r = w.get("lock_ratio")
            lab = f"field {w['h_over_threshold']:.1f} x threshold, {sweep_label(card, w)}"
            bound = "> " if w.get("lock_ratio_is_lower_bound") else ""
            parts.append(f"{lab}: accuracy {w['accuracy']:.2f}" + (f", retention/rewriting time ratio {bound}{r:.3g}"
                                                                    if r else ", not rewritten"))
        lock = "; ".join(parts)
    return {"stores": stores, "writes": writes, "constructs": constructs, "holds": card["loss_law"], "lock": lock}


def run(realizations, out_dir: Path, seed: int = 20260923, transfers: Optional[Dict[str, List[Transfer]]] = None,
        figures: bool = True, quick: bool = False, loops: Optional[List[Dict]] = None) -> List[Dict]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cards = []
    for i, real in enumerate(realizations):
        rng = np.random.default_rng(seed + 101 * i)
        card = evaluate(real, rng, (transfers or {}).get(real.name), quick=quick)
        card["id"] = f"card{i:02d}"
        cards.append(card)
        (out_dir / f"{card['id']}.json").write_text(json.dumps(card, indent=1, default=float), encoding="utf-8")
        print(f"[{card['id']}] {real.name}: {card['verdict']['stores']}; {card['verdict']['writes']} "
              f"({card['seconds']:.0f} s)", flush=True)
    if figures:
        from .visual import card_figure, html_report, loops_figure
        for card in cards:
            card_figure(card, out_dir / f"{card['id']}.png")
        if loops:
            loops_figure(loops, out_dir / "loops.png")
        html_report(cards, out_dir / "index.html", loops=loops)
    return cards
