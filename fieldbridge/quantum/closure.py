"""The closure rule as a construction: index, closed counterpart, excess and a check of the rule
(docs/constructor.md).

Rule. Under a Hamiltonian H the expectation of an observable O moves inside the closure of O: the smallest space of
operators that contains O and is invariant under i[H, .] (Chapter 11, tutorial Module 24). The generator of the
closure, M_ab = <B_a, i[H, B_b]>, is real and antisymmetric, and its eigenvalues +-i omega are the only frequencies
that <O>(t) can show from any prepared state. The rule is sufficient and not necessary: a frequency of the closure
appears with the amplitude the prepared state gives it, and that amplitude can vanish.

Index. The smallest number of terms of the Hamiltonian, as the specification writes them, whose removal leaves a
closure of two or three operators with a rotation frame (the letters O K of the quantum language): 0 for a
realization that already rotates, 1 when one term (an interaction, an anisotropy, a second field) carries every
further frequency, and so on. Every minimal set of terms is listed; the first is used for the excess.

Counterpart. The Hamiltonian without those terms: the closed counterpart, which rotates with one frequency Omega_c
and follows the Rabi law.

Excess. From the prepared state, a top eigenstate of the observable, the exact signal is a sum over the Bohr
frequencies of H. The excess is the spread of these lines, the weighted standard deviation of the Bohr frequencies,
which the counterpart's single line lacks; with it the shift of their weighted mean from Omega_c, the weight outside
the line that continues the counterpart's rotation, and the closure frequencies that the prepared state leaves silent. For a realization of index 0 the excess is the amplitude offset of the signal from sin^2 theta, which is
zero when the prepared state lies along the observable in the closure. The check of the rule: every line of the
signal lies in the closure spectrum.

Scaling, for comparison. With the removed terms multiplied by lambda, the spread falls as lambda^p near zero, p = 1
when the removed terms split the frequency at first order; the exponent between the two smallest lambda is reported. It is a comparison value, not a prediction for a device.
"""
from __future__ import annotations

import itertools
from math import degrees, log
from typing import Dict, List, Optional, Tuple

import numpy as np

from . import language as ql

TOL_FREQ = 1e-7       # relative, for merging Bohr frequencies and matching them to the closure spectrum
TOL_WEIGHT = 1e-9     # relative, for a line to count
RULE = ("the expectation of an observable moves with the frequencies of the closure of the observable under the "
        "Hamiltonian and with no other; sufficient, not necessary")


def prepared_state(Os: np.ndarray) -> np.ndarray:
    """A top eigenstate of the observable, as the Rabi law of the language uses it."""
    _, U = np.linalg.eigh(Os)
    return U[:, -1]


def signal_spectrum(Hs: np.ndarray, Os: np.ndarray, psi0: np.ndarray, tol: float = TOL_FREQ):
    """The exact signal <O>(t) = c + sum_omega 2 Re(a_omega e^{i omega t}) from psi0: the Bohr frequencies
    omega > 0 (merged within the tolerance) with their complex amplitudes, and the constant c."""
    E, V = np.linalg.eigh(Hs)
    c = V.conj().T @ psi0
    Ot = V.conj().T @ Os @ V
    scale = max(1.0, float(E[-1] - E[0]))
    const = 0.0
    pairs: List[Tuple[float, complex]] = []
    for m in range(len(E)):
        for k in range(len(E)):
            w = float(E[m] - E[k])
            amp = complex(np.conj(c[m]) * c[k] * Ot[m, k])
            if abs(w) <= tol * scale:
                const += amp.real
            elif w > 0:
                pairs.append((w, amp))
    pairs.sort(key=lambda p: p[0])
    lines: List[List] = []
    for w, amp in pairs:
        if lines and w - lines[-1][0] <= tol * scale:
            lines[-1][1] += amp
            lines[-1][0] = (lines[-1][0] * lines[-1][2] + w) / (lines[-1][2] + 1)
            lines[-1][2] += 1
        else:
            lines.append([w, amp, 1])
    return [(float(w), complex(a)) for w, a, _ in lines], float(const)


def _lines(Hs, Os, psi0):
    lines, const = signal_spectrum(Hs, Os, psi0)
    weights = [(w, 2.0 * abs(a)) for w, a in lines]
    total = sum(x for _, x in weights)
    return [(w, x) for w, x in weights if x > TOL_WEIGHT * max(1.0, total)], total, const


def closure_spectrum(Hs: np.ndarray, O0: np.ndarray):
    """The closure of the traceless observable, its dimension and its frequencies."""
    basis = ql.closure_basis(Hs, O0)
    n = len(basis)
    freqs = ql._frequencies(ql.closure_generator(Hs, basis)) if 1 < n <= ql.CLOSURE_MAX else []
    return basis, n, freqs


def rotation_frame(Hk: np.ndarray, O0: np.ndarray):
    """The closure of O0 under Hk and its rotation frame when the closure has two or three operators."""
    cb = ql.closure_basis(Hk, O0)
    if 2 <= len(cb) <= 3:
        return cb, ql.canonical_closure(cb, Hk, O0)
    return cb, None


def minimal_removals(real: ql.Realization, V: np.ndarray, O0: np.ndarray, max_removed: Optional[int] = None):
    """The index: the smallest number of Hamiltonian terms whose removal leaves a rotating closure, with every set of
    terms that attains it. (None, []) when no removal short of the whole Hamiltonian does."""
    T = len(real.terms)
    limit = T - 1 if max_removed is None else min(max_removed, T - 1)
    for r in range(0, limit + 1):
        found = []
        for drop in itertools.combinations(range(T), r):
            keep = [m for i, (_, m) in enumerate(real.terms) if i not in drop]
            Hk = V.conj().T @ sum(keep) @ V
            cb, canon = rotation_frame(Hk, O0)
            if canon is not None:
                found.append({"removed": [real.terms[i][0] for i in drop], "removed_indices": list(drop), "H": Hk,
                              "canon": canon, "closure_dimension": len(cb)})
        if found:
            return r, found
    return None, []


def _excess(lines, total, rate):
    """The excess of the signal over one rotation at the counterpart's frequency: the spread of its lines (the
    weighted standard deviation of the Bohr frequencies), the shift of their weighted mean from the counterpart's
    frequency, the line nearest to that frequency with its shift, and the weight outside that line."""
    if not lines:
        return {"lines": 0, "spread": 0.0, "relative_spread": 0.0, "mean_shift": -rate, "counterpart_line_frequency": 0.0,
                "frequency_shift": -rate, "weight_outside_counterpart_line": 0.0}
    w = np.array([x for x, _ in lines]); x = np.array([y for _, y in lines]) / total
    mean = float(np.sum(x * w))
    spread = float(np.sqrt(max(np.sum(x * (w - mean) ** 2), 0.0)))
    near = min(range(len(lines)), key=lambda i: abs(lines[i][0] - rate))
    return {"lines": len(lines), "spread": spread, "relative_spread": spread / rate if rate > 0 else float("nan"),
            "mean_shift": mean - rate, "counterpart_line_frequency": float(lines[near][0]),
            "frequency_shift": float(lines[near][0] - rate),
            "weight_outside_counterpart_line": float(sum(y for i, (_, y) in enumerate(lines) if i != near) / total)}


def construct(real: ql.Realization, max_removed: Optional[int] = None,
              scales: Tuple[float, ...] = (0.5, 0.25, 0.125, 0.0625)) -> Dict[str, object]:
    """Index, closed counterpart, excess, the check of the rule and the scaling, for one realization."""
    V, sec = ql.restrict(real)
    Hs, Os = V.conj().T @ real.H @ V, V.conj().T @ real.O @ V
    O0 = ql._traceless(Os)
    out: Dict[str, object] = {"name": real.name, "field": real.field, "carrier": real.carrier.description,
                              "sector_dimension": int(Hs.shape[0]), "terms": [t for t, _ in real.terms], "rule": RULE}
    _, n, freqs = closure_spectrum(Hs, O0)
    out["closure"] = {"dimension": n, "frequencies": freqs}
    psi0 = prepared_state(Os)
    lines, total, const = _lines(Hs, Os, psi0)
    scale = max(1.0, max(freqs) if freqs else 1.0, max((w for w, _ in lines), default=1.0))
    dist = [min((abs(w - f) for f in freqs), default=float("inf")) for w, _ in lines]
    silent = [f for f in freqs if min((abs(w - f) for w, _ in lines), default=float("inf")) > 1e-6 * scale]
    out["signal"] = {"frequencies": [w for w, _ in lines], "weights": [x for _, x in lines], "constant": const,
                     "lines_outside_closure": int(sum(d > 1e-6 * scale for d in dist)),
                     "largest_distance_to_closure": float(max(dist, default=0.0)), "silent_closure_frequencies": silent}
    out["rule_holds"] = out["signal"]["lines_outside_closure"] == 0
    if n <= 1:
        out.update(index=None, counterparts=[], excess=None,
                   obstruction="the observable is conserved: nothing moves" if n == 1 else
                   "the observable is a multiple of the identity: nothing moves")
        out["statement"] = statement(out)
        return out
    index, found = minimal_removals(real, V, O0, max_removed)
    out["index"] = index
    counterparts = []
    for f in found:
        canon = f["canon"]
        law = ql.rabi_law(f["H"], Os, canon)
        counterparts.append({"removed": f["removed"], "closure_dimension": f["closure_dimension"],
                             "rate": float(canon["rate"]), "theta_deg": float(degrees(canon["theta"])),
                             "amplitude_offset": float(canon["amplitude_offset"]), "law_residual": float(law["residual"])})
    out["counterparts"] = counterparts
    if index is None:
        out.update(excess=None, obstruction="no removal of terms short of the whole Hamiltonian leaves a rotating closure")
        out["statement"] = statement(out)
        return out
    cp0 = counterparts[0]
    excess = _excess(lines, total, cp0["rate"])
    excess["counterpart_frequency"] = cp0["rate"]
    excess["silent_closure_frequencies"] = len(silent)
    if index == 0:
        excess["amplitude_offset"] = cp0["amplitude_offset"]
    out["excess"] = excess
    if index > 0:
        drop = set(found[0]["removed_indices"])
        points = []
        for lam in (1.0,) + tuple(scales):
            Hl = V.conj().T @ sum((lam if i in drop else 1.0) * m for i, (_, m) in enumerate(real.terms)) @ V
            ls, tot, _ = _lines(Hl, Os, psi0)
            e = _excess(ls, tot, cp0["rate"])
            points.append({"lambda": lam, "spread": e["spread"], "mean_shift": e["mean_shift"],
                           "weight_outside": e["weight_outside_counterpart_line"]})
        a, b = points[-2]["spread"], points[-1]["spread"]
        p = log(a / b) / log(points[-2]["lambda"] / points[-1]["lambda"]) if a > 1e-12 and b > 1e-12 else None
        out["scaling"] = {"points": points, "exponent": p}
    out["statement"] = statement(out)
    return out


def _fmt(values: List[float], limit: int = 4) -> str:
    return ", ".join(f"{v:.4g}" for v in values[:limit]) + (", ..." if len(values) > limit else "")


def statement(r: Dict[str, object]) -> str:
    c = r["closure"]
    head = (f"The closure of the observable has {c['dimension']} operator{'s' if c['dimension'] != 1 else ''}"
            + (f" and {len(c['frequencies'])} frequenc{'y' if len(c['frequencies']) == 1 else 'ies'} ({_fmt(c['frequencies'])})"
               if c["frequencies"] else "") + ".")
    if r.get("obstruction"):
        return head + " " + r["obstruction"].capitalize() + "."
    cp = r["counterparts"][0]
    e = r["excess"]
    if r["index"] == 0:
        s = (head + f" The realization rotates at {cp['rate']:.4g} (index 0). From the top eigenstate of the observable "
             f"the signal has {e['lines']} line{'s' if e['lines'] != 1 else ''}")
        if abs(e.get("amplitude_offset", 0.0)) > ql.ALONG_TOL:
            s += (f", with an amplitude offset of {e['amplitude_offset']:.4g} from sin^2 theta: the prepared state does "
                  f"not lie along the observable in the closure")
        else:
            s += " with the Rabi amplitude"
    else:
        rem = " and ".join(cp["removed"])
        s = (head + f" Removing {r['index']} term{'s' if r['index'] != 1 else ''} ({rem}) leaves a rotation at "
             f"{cp['rate']:.4g}, the closed counterpart (index {r['index']}"
             + (f"; {len(r['counterparts'])} minimal sets of terms" if len(r["counterparts"]) > 1 else "") + "). "
             f"From the top eigenstate of the observable the signal has {e['lines']} lines, spread over "
             f"{e['spread']:.4g} ({e['relative_spread']:.3f} of the counterpart's frequency) around a mean shifted by "
             f"{e['mean_shift']:+.4g}; a share {e['weight_outside_counterpart_line']:.3f} of the oscillating amplitude "
             f"lies outside the line that continues the counterpart's rotation")
        if r.get("scaling") and r["scaling"]["exponent"] is not None:
            s += f"; with the removed terms scaled by lambda the spread falls as lambda^{r['scaling']['exponent']:.2f}"
    s += (f". Every line lies in the closure spectrum" if r["rule_holds"] else
          f". {r['signal']['lines_outside_closure']} lines lie outside the closure spectrum: the rule is violated")
    if e.get("silent_closure_frequencies"):
        s += f"; {e['silent_closure_frequencies']} closure frequenc{'y is' if e['silent_closure_frequencies'] == 1 else 'ies are'} silent for this prepared state"
    return s + "."


def markdown(rows: List[Dict[str, object]]) -> str:
    out = ["# The closure rule as a construction", "",
           "Index: the smallest number of Hamiltonian terms removed so that the closure of the observable rotates. "
           "Counterpart: the Hamiltonian without them, with one frequency. Excess: the spread of the signal's lines "
           "(weighted standard deviation of its Bohr frequencies, relative to the counterpart's frequency) from a top "
           "eigenstate of the observable, and the shift of their weighted mean from that frequency. The rule holds when every line of the signal lies in the "
           "closure spectrum; it is sufficient and not necessary, since closure frequencies can be silent.", "",
           "| Realization | Dimension | Closure | Frequencies | Index | Removed | Counterpart rate | Relative spread | "
           "Mean shift | Exponent | Rule |", "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        c = r["closure"]
        e = r.get("excess") or {}
        cps = r.get("counterparts") or []
        cp = cps[0] if cps else {}
        sc = r.get("scaling") or {}
        index = "" if r.get("index") is None else str(r["index"])
        removed = "; ".join(cp.get("removed", []))
        rate = f"{cp['rate']:.4g}" if cp else ""
        weight = f"{e['relative_spread']:.3f}" if e else ""
        shift = f"{e['mean_shift']:+.3g}" if e else ""
        expo = f"{sc['exponent']:.2f}" if sc.get("exponent") is not None else ""
        rule = "holds" if r["rule_holds"] else "violated"
        out.append(f"| {r['name']} | {r['sector_dimension']} | {c['dimension']} | {len(c['frequencies'])} | {index} | "
                   f"{removed} | {rate} | {weight} | {shift} | {expo} | {rule} |")
    out.append("")
    for r in rows:
        out.extend([f"## {r['name']}", "", r["statement"], ""])
    return "\n".join(out)
