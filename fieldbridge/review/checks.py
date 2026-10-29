"""Four checks of a claim in a manuscript (``fieldbridge review``).

The nested description assigns every quantity to a level: the operation, the carrier (the states), the closure (the
constraints), the observable, the protocol (preparation and observation) and the realization (the parameters). A
claim holds only at the level at which it was computed. The four checks test the places where a claim most often
moves from one level to another, or outside what supports it.

1. same_system: the quantity was computed for the system it is claimed for (same carrier, free and fixed parameters,
   closure and realization). Where the review gives equation and unknown counts, it also compares them with the
   parameters fixed, which shows whether the claimed state exists for a prescribed system.
2. protocol_independence: a property of the state does not change with the protocol (run length, rate, preparation),
   and a thermodynamic statement does not rest on a kinetic observable alone. Includes the observation-window test of
   a fitted divergence.
3. theorem_scope: each theorem used is applied inside its hypotheses and supports the inference drawn
   (theorems.py).
4. prior_art: novelty claims and the statements a review relies on carry a recorded search with a verdict; a
   published result that contradicts the claim fails it.

Each check returns a status, pass, warn, fail or not applicable, and an explanation in plain sentences. The checks
read what the review file states; they do not read the manuscript, and a pass is a necessary condition, not a proof.
"""
from __future__ import annotations

import math
from typing import Dict, List

from . import theorems

NA = "not applicable"


def _result(name, status, explanation, **extra):
    out = {"check": name, "status": status, "explanation": explanation}
    out.update(extra)
    return out


def _evaluate_count(expr, N):
    import sympy as sp
    n = sp.Symbol("N")                                      # N is also SymPy's numeric-evaluation function
    return float(sp.sympify(str(expr), locals={"N": n}).subs(n, N))


# 1 --------------------------------------------------------------------------------------------------------------
def same_system(claim: Dict, systems: Dict) -> Dict:
    a, b = claim.get("computed_on"), claim.get("claimed_for")
    if not a or not b:
        return _result("same_system", NA, "The claim does not say on which system it was computed and for which it is "
                                          "claimed.")
    if a not in systems or b not in systems:
        return _result("same_system", "fail", f"Unknown system: {a if a not in systems else b}.")
    sa, sb = systems[a], systems[b]
    diffs: List[str] = []
    for key, label in (("variables", "state variables"), ("free_parameters", "free parameters"),
                       ("fixed_parameters", "fixed parameters")):
        xa, xb = set(sa.get("carrier", {}).get(key, [])), set(sb.get("carrier", {}).get(key, []))
        if xa != xb:
            only_a, only_b = sorted(xa - xb), sorted(xb - xa)
            parts = []
            if only_a:
                parts.append(f"only in the computed system: {', '.join(only_a)}")
            if only_b:
                parts.append(f"only in the claimed system: {', '.join(only_b)}")
            diffs.append(f"{label} differ ({'; '.join(parts)})")
    ca, cb = set(sa.get("closure", [])), set(sb.get("closure", []))
    if ca != cb:
        diffs.append("constraints differ (" + "; ".join(
            ([f"only in the computed system: {', '.join(sorted(ca - cb))}"] if ca - cb else []) +
            ([f"only in the claimed system: {', '.join(sorted(cb - ca))}"] if cb - ca else [])) + ")")
    ra, rb = sa.get("realization", {}), sb.get("realization", {})
    rdiff = sorted(k for k in set(ra) | set(rb) if ra.get(k) != rb.get(k))
    if rdiff:
        diffs.append("realization differs in " + ", ".join(f"{k} ({ra.get(k)} against {rb.get(k)})" for k in rdiff))
    counts = []
    for name, s in (((a, sa), (b, sb)) if a != b else ((a, sa),)):
        c = s.get("counts")
        if c and "equations" in c and "unknowns" in c:
            rows = []
            for N in c.get("sizes", [64, 256, 1024]):
                e, u = _evaluate_count(c["equations"], N), _evaluate_count(c["unknowns"], N)
                rows.append((N, e, u))
            over = [r for r in rows if r[1] > r[2]]
            text = (f"In system {name} the counts are {c['equations']} equations for {c['unknowns']} unknowns"
                    + (f"; the equations exceed the unknowns by {c['equations']} - ({c['unknowns']}), so the state "
                       "exists only for parameters on a set of measure zero" if over and len(over) == len(rows) else ""))
            counts.append({"system": name, "rows": rows, "overdetermined": bool(over), "text": text + "."})
    if not diffs:
        status, text = "pass", f"The quantity is computed on the system it is claimed for ({a})."
    else:
        status = "fail"
        text = (f"The quantity is computed on system '{a}' and claimed for system '{b}', and they differ: "
                + "; ".join(diffs) + ". A property of one system is not a property of the other unless the "
                "difference cannot change the quantity, which the claim must show.")
    for c in counts:
        text += " " + c["text"]
        if c["system"] == b and c["overdetermined"]:
            if a != b:
                status = "fail"
            elif status == "pass":
                status = "warn"
                text += (" The state therefore exists only for tuned parameters, and the claim holds for the tuned "
                         "system, not for a prescribed one.")
    return _result("same_system", status, text, counts=counts)


# 2 --------------------------------------------------------------------------------------------------------------
def _trend(estimates):
    """Weighted least-squares slope of the estimate against the logarithm of the protocol setting, in units of its
    standard error."""
    xs = [math.log(float(e["setting"])) for e in estimates]
    ys = [float(e["value"]) for e in estimates]
    ws = [1.0 / max(float(e.get("error", 0.0)) ** 2, 1e-300) for e in estimates]
    W = sum(ws)
    xm = sum(w * x for w, x in zip(ws, xs)) / W
    ym = sum(w * y for w, y in zip(ws, ys)) / W
    sxx = sum(w * (x - xm) ** 2 for w, x in zip(ws, xs))
    if sxx == 0:
        return 0.0, float("inf")
    slope = sum(w * (x - xm) * (y - ym) for w, x, y in zip(ws, xs, ys)) / sxx
    return slope, slope / math.sqrt(1.0 / sxx)


def protocol_independence(claim: Dict) -> Dict:
    p = claim.get("protocol")
    kind = claim.get("kind", "")
    obs = claim.get("observable_type", "")
    msgs: List[str] = []
    status = "pass"
    if kind == "thermodynamic" and obs in ("relaxation time", "retention time", "escape time", "fitted divergence"):
        if not claim.get("equilibrium_evidence"):
            status = "fail"
            msgs.append(f"A thermodynamic statement rests on a kinetic observable ({obs}) without equilibrium "
                        "evidence (a free energy, an entropy or a coexistence). How long a prepared state survives "
                        "depends on its preparation and on the time of observation; a thermodynamic property does "
                        "not.")
    if not p:
        if not msgs:
            return _result("protocol_independence", NA, "No protocol-dependent estimate is given for this claim.")
        return _result("protocol_independence", status, " ".join(msgs))
    if p.get("estimates") and len(p["estimates"]) >= 3:
        slope, z = _trend(p["estimates"])
        quantity = p.get("quantity", "the quantity")
        setting = p.get("setting", "the protocol setting")
        if abs(z) > 2:
            status = "fail"
            msgs.append(f"The estimate of {quantity} changes with {setting}: slope {slope:.3g} per e-fold, "
                        f"{abs(z):.1f} standard errors. It is a property of the protocol as well as of the state.")
        else:
            msgs.append(f"The estimate of {quantity} does not change significantly with {setting} "
                        f"({abs(z):.1f} standard errors).")
    w = p.get("window")
    if w:
        est = float(w["estimate"])
        lo = w.get("warmest_without_relaxation")
        hi = w.get("coldest_with_relaxation")
        err = w.get("error")
        if hi is not None and est <= float(hi) and float(hi) - est < 0.05 * float(hi):
            status = "fail"
            msgs.append(f"The fitted divergence {est:g} lies within 5 % below the coldest temperature at which "
                        f"relaxation was still seen ({float(hi):g}): the fit places it at the edge of the observation "
                        "window, which a longer run would move.")
        elif lo is not None and hi is not None and float(lo) <= est <= float(hi):
            status = "warn" if status == "pass" else status
            msgs.append(f"The fitted divergence {est:g} lies in the interval {float(lo):g}-{float(hi):g} where no "
                        "relaxation time was measured; its position there is an extrapolation of the fitted form.")
        if err is not None and lo is not None and hi is not None and float(err) < 0.5 * (float(hi) - float(lo)):
            status = "fail"
            msgs.append(f"The quoted uncertainty {float(err):g} is smaller than half the unmeasured interval "
                        f"{float(lo):g}-{float(hi):g}; it is conditional on the fitted form, not set by the data.")
        if w.get("at_fit_bound"):
            status = "fail"
            msgs.append("The fitted divergence sits on a bound of the fit.")
        prof = w.get("profile")
        if prof:
            vals = [float(v) for _, v in prof]
            span = [float(t) for t, _ in prof]
            vmin = min(vals)
            ok = [t for t, v in zip(span, vals) if v <= 2 * vmin]
            if ok and (max(ok) - min(ok)) > 0.5 * (max(span) - min(span)):
                status = "fail" if status != "fail" else status
                msgs.append(f"The fit residual stays within a factor two of its minimum for divergence points from "
                            f"{min(ok):g} to {max(ok):g}; the data do not determine the divergence.")
        if w.get("alternative_without_divergence_fits"):
            status = "fail"
            msgs.append("A form without any divergence at positive temperature fits the same data comparably.")
    return _result("protocol_independence", status, " ".join(msgs) if msgs else "Nothing to test.")


# 3 --------------------------------------------------------------------------------------------------------------
def theorem_scope(claim: Dict) -> Dict:
    uses = claim.get("theorems", [])
    if not uses:
        return _result("theorem_scope", NA, "The claim cites no theorem.")
    results = [theorems.check_use(u["name"], u.get("case", {}), u.get("inference", "")) for u in uses]
    order = {"fail": 3, "unknown": 2, "warn": 1, "pass": 0}
    worst = max(results, key=lambda r: order[r["status"]])["status"]
    status = {"unknown": "warn"}.get(worst, worst)
    text = " ".join(f"{r['theorem']}: {r['explanation']}" for r in results)
    return _result("theorem_scope", status, text, uses=results)


# 4 --------------------------------------------------------------------------------------------------------------
VERDICTS = ("published", "partly published", "not found")


def prior_art(claim: Dict) -> Dict:
    items = claim.get("prior_art", [])
    novel = bool(claim.get("novelty_claimed"))
    if not items:
        if novel:
            return _result("prior_art", "fail", "Novelty is claimed, but no search of the literature is recorded.")
        return _result("prior_art", NA, "No novelty is claimed and no search is recorded.")
    msgs, status = [], "pass"
    for it in items:
        st = it.get("statement", "(no statement)")
        q = it.get("queries", [])
        v = it.get("verdict")
        if v not in VERDICTS:
            status = "fail"
            msgs.append(f"'{st}': the verdict must be one of {', '.join(VERDICTS)}.")
            continue
        if not q:
            status = "fail"
            msgs.append(f"'{st}': no search queries are recorded, so the verdict cannot be repeated.")
        if v == "not found" and len(q) < 2:
            status = "fail" if status == "fail" else "warn"
            msgs.append(f"'{st}': 'not found' rests on fewer than two recorded queries.")
        if v in ("published", "partly published") and not it.get("sources"):
            status = "fail"
            msgs.append(f"'{st}': the verdict '{v}' gives no source.")
        if it.get("contradicts") and v in ("published", "partly published"):
            status = "fail"
            msgs.append(f"'{st}': a published result contradicts the claim: " + "; ".join(it.get("sources", [])) + ".")
            continue
        if novel and v == "published" and it.get("bears_on_novelty", True):
            status = "fail"
            msgs.append(f"'{st}' is claimed as new but is published: " + "; ".join(it.get("sources", [])) + ".")
        elif v == "partly published" and novel:
            status = "warn" if status == "pass" else status
            msgs.append(f"'{st}' is partly published (" + "; ".join(it.get("sources", [])) + "); the claim must "
                        "state what is new against these sources.")
    if not msgs:
        msgs.append("Every statement has a recorded search and a verdict.")
    return _result("prior_art", status, " ".join(msgs))


CHECKS = (("same_system", lambda c, s: same_system(c, s)),
          ("protocol_independence", lambda c, s: protocol_independence(c)),
          ("theorem_scope", lambda c, s: theorem_scope(c)),
          ("prior_art", lambda c, s: prior_art(c)))


def review(doc: Dict) -> Dict:
    systems = doc.get("systems", {})
    out = []
    for claim in doc.get("claims", []):
        results = [fn(claim, systems) for _, fn in CHECKS]
        out.append({"id": claim.get("id"), "statement": claim.get("statement"), "results": results,
                    "failed": [r["check"] for r in results if r["status"] == "fail"]})
    tally = {k: sum(1 for c in out for r in c["results"] if r["status"] == k)
             for k in ("pass", "warn", "fail", NA)}
    return {"paper": doc.get("paper", {}), "claims": out, "tally": tally}
