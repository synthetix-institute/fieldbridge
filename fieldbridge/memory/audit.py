"""Audit of a derived writing law: the specific ways it could be wrong, each settled by a test of its own.

A law that the constructor derives rests on assumptions that the derivation itself does not test: the class of the
reduced equation read on one window of the fit, a law measured with the bias toward one state only, a noise ratio
chosen for the test rather than the element's own, a time step, a write that the description of the element cannot
represent, a write at a value of the control that the published model does not allow. For each assumption the audit
states the error that would follow if it failed, runs the one test that decides it, and records the result:

  confirmed       the test finds the error; the law as derived does not hold for this element
  qualified       the law holds under a condition that the test states (a range of sweep rates, for example)
  refuted         the test finds no effect beyond its statistical error
  not applicable  the law does not rest on this assumption
  not testable    the test cannot be run on this specification

Proposing a specific error and deciding it is easier than deciding a whole law at once. The verdict is "ready for
review" when no error is confirmed: the law then goes to a person or to the literature, and it is not certified by
the audit. Otherwise the verdict lists the confirmed errors, each with the assumption it breaks and its evidence.
"""
from __future__ import annotations

from dataclasses import replace
from math import pi, sqrt
from statistics import NormalDist
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from . import analysis as an
from . import codiscovery as cd
from . import noise_ratio as nr
from .construct import along_path, cancel_asymmetry, obstruction, swept_write_check, sweep_normal_form
from .identity import Realization

CONFIRMED, REFUTED, QUALIFIED = "confirmed", "refuted", "qualified"
NOT_APPLICABLE, NOT_TESTABLE = "not applicable", "not testable"
SIGMA = 3.0                 # a difference counts when it exceeds this many combined standard errors ...
FLOOR = 0.03                # ... and this absolute amount in units of the law constant
FEW_QUANTA = 10.0           # below this many quanta the semiclassical description of a write is not reliable
OWN_NOISE_RATIO = 1.0       # above this noise ratio the write constant departs from pi^(1/4) (noise_ratio.write_constant)


# ------------------------------------------------------------------------------------------------ the law under audit
def symmetric_write_target(real: Realization, rng) -> Tuple[Optional[Realization], Optional[Dict], str]:
    """The realization and normal form on which the symmetric-write derivation measures its law (the model itself, or
    the sweep through the cusp found by a second parameter), as codiscovery chooses them; with a note."""
    events = an.locate_writes(real, rng, n_starts=24 if real.carrier.dim <= 4 else 12)
    nfs = []
    for ev in events:
        nf = an.normal_form(real, ev["q"], real.control, ev["v"])
        nf["obstruction"] = obstruction(nf)
        nfs.append(nf)
    pick = next((nf for nf in nfs if nf["obstruction"]["realizes_generator"]), None)
    if pick is not None:
        return real, pick, "the model itself"
    fold = next((nf for nf in nfs if (nf["obstruction"]["asymmetric"] or nf["obstruction"]["biased"])
                 and not nf.get("hopf")), None)
    if fold is None:
        return None, None, "no write point that a law can be measured on"
    for second in real.material:
        with np.errstate(all="ignore"):
            try:
                c = cancel_asymmetry(real, fold, second)
            except np.linalg.LinAlgError:
                continue
        if c.get("success") and c.get("obstruction_along_sweep", {}).get("realizes_generator"):
            path = along_path(real, {real.control: c["control"], second: c["value"]}, c["direction"])
            nfp = sweep_normal_form(path, np.asarray(c["state"]))
            nfp["obstruction"] = obstruction(nfp)
            nfp["_second"], nfp["_second_value"], nfp["_control_at_write"] = second, float(c["value"]), float(c["control"])
            return path, nfp, f"the sweep through the cusp found along {second}"
    return None, None, "no cusp along the material parameters"


def _constant(law: Dict) -> Dict:
    return cd.law_constant(law) if law else {}


def _differs(a: Dict, b: Dict) -> Tuple[bool, float]:
    d = abs(a["constant"] - b["constant"])
    e = sqrt(a.get("stderr", 0.0) ** 2 + b.get("stderr", 0.0) ** 2)
    return d > max(SIGMA * e, FLOOR * abs(a["constant"])), d


def _item(error: str, assumption: str, test: str, status: str, **evidence) -> Dict[str, object]:
    return {"error": error, "assumption": assumption, "test": test, "status": status, "evidence": evidence}


# ------------------------------------------------------------------------------------------------ the checks
def check_window(real, target, nf, rng, base) -> Dict:
    """The class of the reduced equation is read on one window of the fit; read on half that window, closer to the
    write point, it must not change. A wider window takes up terms of higher order (at a cusp, the quartic term reads as
    a quadratic one), so its class is reported but does not decide."""
    error = "the class of the reduced equation is an artefact of the window on which the drift is fitted"
    assumption = "the reduced equation keeps its class on a narrower window"
    test = "refit the reduced equation on half the window and compare the classes (twice the window is reported)"
    half = float(nf.get("half", 0.0) or 0.0)
    if half <= 0 or nf.get("hopf"):
        return _item(error, assumption, test, NOT_APPLICABLE)
    param = target.control
    q = np.asarray(nf["state"], float)
    kinds = {}
    for label, h in (("half", 0.5 * half), ("used", half), ("double", 2.0 * half)):
        try:
            kinds[label] = an.normal_form(target, q, param, float(nf["value"]), half=h)["kind"].split(":")[0]
        except (ValueError, np.linalg.LinAlgError) as e:
            kinds[label] = f"not read ({type(e).__name__})"
    used = nf["kind"].split(":")[0]
    changed = kinds["half"] != used
    return _item(error, assumption, test, CONFIRMED if changed else REFUTED, classes=kinds, reported=used,
                 window=half)


def check_one_sign(real, target, nf, rng, base) -> Dict:
    """A law measured with the bias toward one state only contains any offset of the write; with the bias toward
    either state the response and the offset separate."""
    error = "the constant mixes the response to the bias with an offset of the write"
    assumption = "the write has no offset of its own"
    test = "measure the law with the bias toward either state and separate the response from the offset"
    const = base["constant"]
    if const.get("by_sign"):
        law, two = base["law"], const
    else:
        law = swept_write_check(target, nf, rng, n_traj=base["n_traj"], replicates=base["replicates"], signs=(1, -1))
        two = _constant(law)
    if not two or "offset" not in two:
        return _item(error, assumption, test, NOT_TESTABLE, note="the two-sign law was not measured")
    # stricter than the derivation's own report (two standard errors): an audit must rarely reject a valid law
    resolved = abs(two["offset"]) > max(SIGMA * two["offset_stderr"], FLOOR * abs(two["constant"]))
    status = CONFIRMED if resolved else REFUTED
    return _item(error, assumption, test, status, response=two["constant"], offset=two["offset"],
                 offset_stderr=two["offset_stderr"], one_sign_constant=const.get("constant"))


def check_own_noise(real, target, nf, rng, base) -> Dict:
    """The law is tested at a noise ratio Gamma = 0.05, where the constant is pi^(1/4). At the element's own noise the
    ratio can be large, and the constant is then c(Gamma) of the noise-ratio law."""
    error = "the constant pi^(1/4) does not apply at the element's own noise and the sweep rate of the test"
    assumption = "the element's own noise ratio is small (Gamma <= 1)"
    test = ("compute the noise ratio at the element's own noise, the constant c(Gamma) there, and the sweep rate above "
            "which pi^(1/4) applies")
    own = (base["law"] or {}).get("at_material_noise") or {}
    if "gamma" not in own:
        return _item(error, assumption, test, NOT_TESTABLE)
    g = float(own["gamma"])
    rate = float((base["law"] or {}).get("rate", 0.02))
    c = nr.write_constant(g)
    # Gamma = D_s |a3| / r, so pi^(1/4) applies at the element's own noise for sweeps faster than r = D_s |a3| / 1
    status = QUALIFIED if g > OWN_NOISE_RATIO else REFUTED
    return _item(error, assumption, test, status, own_noise_ratio=g, constant_at_own_noise=c,
                 regime=nr.regime(g), rate_of_test=rate, rate_above_which_pi_quarter_applies=g * rate / OWN_NOISE_RATIO,
                 measured_at_own_noise=own.get("measured"), stderr=own.get("stderr"))


def check_time_step(real, target, nf, rng, base) -> Dict:
    """The stochastic integration must not set the outcome: halving the time step must leave the constant."""
    error = "the constant depends on the time step of the integration"
    assumption = "the time step is small enough"
    test = "measure the law again with half the time step"
    const = base["constant"]
    if not const:
        return _item(error, assumption, test, NOT_TESTABLE)
    finer = replace(target)
    try:
        object.__setattr__(finer, "_safe_dt", an.safe_dt(target) / 2)
    except (AttributeError, TypeError):
        return _item(error, assumption, test, NOT_TESTABLE, note="the time step cannot be set")
    signs = (1, -1) if const.get("by_sign") else (1,)
    law = swept_write_check(finer, nf, rng, n_traj=base["n_traj"], replicates=base["replicates"], signs=signs)
    other = _constant(law)
    if not other:
        return _item(error, assumption, test, NOT_TESTABLE)
    differs, d = _differs(const, other)
    return _item(error, assumption, test, CONFIRMED if differs else REFUTED, constant=const["constant"],
                 constant_half_step=other["constant"], difference=d, dt=an.safe_dt(target), dt_half=an.safe_dt(finer))


def _quantum_description(real: Realization) -> bool:
    text = " ".join(str(real.spec.get(k, "")) for k in ("closure", "kind")).lower() if hasattr(real, "spec") else ""
    prov = str((getattr(real, "spec", {}) or {}).get("provenance", {})).lower()
    return "vacuum" in text or "truncated wigner" in text or "vacuum noise" in prov or "truncated wigner" in prov


def check_few_quanta(real, target, nf, rng, base) -> Dict:
    """For a semiclassical description of a quantum element (mean field with the vacuum noise of its loss), the
    written states must hold many quanta; at a few quanta the quantum selection differs from the semiclassical one."""
    error = "the write lies outside the range in which the semiclassical description holds"
    assumption = f"the written states hold more than {FEW_QUANTA:g} quanta"
    test = ("count the quanta that separate the written states from the unwritten one at the end of the sweep, in "
            "units of the vacuum noise")
    if not _quantum_description(real):
        return _item(error, assumption, test, NOT_APPLICABLE)
    q0 = np.asarray(nf["state"], float)
    rate = float((base["law"] or {}).get("rate", 0.02))
    g = abs(float(nf.get("kappa_slope", 0.0)))
    if g <= 0:
        return _item(error, assumption, test, NOT_TESTABLE, note="the growth rate of the critical mode is not resolved")
    # the end of the sweep of the law test (construct.swept_write_check): 4/sqrt(rate) past the write point
    end = {target.control: float(nf["value"]) + float(np.sign(-nf["kappa_slope"])) * (rate / g) * 4.0 / sqrt(rate)}
    from .completion import stable_states
    try:
        states = stable_states(target, rng, end, n_starts=24)
    except (ValueError, np.linalg.LinAlgError):
        states = []
    if not states:
        return _item(error, assumption, test, NOT_TESTABLE, note="no written state found past the write point")
    D = float(real.noise)
    quanta = [float(np.sum(np.asarray(target.carrier.diff(q0, s)) ** 2) / (8.0 * D)) for s in states]
    n = min(quanta)
    return _item(error, assumption, test, CONFIRMED if n < FEW_QUANTA else REFUTED, quanta=n,
                 note="the noise of the description is taken as the vacuum noise, D = kappa/8 per quadrature")


def periodic_partner(real, name: str, period: float = 2.0 * np.pi, points: int = 6) -> Optional[int]:
    """The variable whose shift by `period` compensates a shift of the parameter `name` by the same amount, so that
    the parameter is defined modulo the period (an external flux phase with the phase of the loop); None if there is
    none. Tested at random states, with a relative tolerance."""
    rng = np.random.default_rng(0)
    qs = real.carrier.sample(rng, points)
    shifted = {name: float(real.params[name]) + period}
    for k in range(real.carrier.dim):
        dq = np.zeros(real.carrier.dim)
        dq[k] = period
        a = real.F(qs)
        b = real.F(qs + dq, **shifted)
        if np.all(np.isfinite(a)) and np.allclose(a, b, rtol=1e-9, atol=1e-9 * max(1.0, float(np.abs(a).max()))):
            return k
    return None


def check_published_sign(real, target, nf, rng, base) -> Dict:
    """The write must occur at values of the control, and of a second parameter that unfolds a cusp, of the same sign
    as their published values: a rate, a mass or a loss made negative has no physical meaning. A parameter defined
    modulo 2 pi (a shift by 2 pi is compensated by a shift of a variable) is compared at its value nearest the
    published one."""
    error = "the write occurs at a value of a parameter that the published model does not allow"
    assumption = "the control, and the second parameter at a cusp, keep the signs of their published values"
    test = ("compare the signs of the parameters at the write point with the signs of their published values, a "
            "parameter defined modulo 2 pi at its value nearest the published one")
    published = dict((getattr(real, "spec", {}) or {}).get("parameters") or {})
    at_write = {real.control: float(nf.get("_control_at_write", nf["value"]))}
    if nf.get("_second"):
        at_write[nf["_second"]] = float(nf["_second_value"])
    pairs = {k: (float(published[k]), v) for k, v in at_write.items() if k in published and float(published[k]) != 0.0}
    if not pairs:
        return _item(error, assumption, test, NOT_APPLICABLE)
    evidence, opposite = {}, []
    for k, (pub, at) in pairs.items():
        row = [pub, at]
        if pub * at < 0 and k in real.params and periodic_partner(real, k) is not None:
            at = at + 2.0 * np.pi * round((pub - at) / (2.0 * np.pi))
            row.append(at)
            evidence[f"{k} (published, at the write, equivalent value modulo 2 pi)"] = row
        else:
            evidence[f"{k} (published, at the write)"] = row
        if pub * at < 0:
            opposite.append(k)
    return _item(error, assumption, test, CONFIRMED if opposite else REFUTED, **evidence)


CHECKS: Dict[str, Callable] = {"window": check_window, "one sign": check_one_sign, "own noise": check_own_noise,
                               "time step": check_time_step, "few quanta": check_few_quanta,
                               "published sign": check_published_sign}


# ------------------------------------------------------------------------------------------------ the audit
def audit_symmetric_write(real: Realization, rng, n_traj: int = 4000, replicates: int = 4,
                          checks: Optional[List[str]] = None) -> Dict[str, object]:
    """Derive the symmetric write of `real` with its law, then run every check on it. Returns the derivation's status
    and law, the ledger of checks and the verdict."""
    derivation = cd.TARGETS["symmetric-write"](real, rng, check_law=True, n_traj=n_traj, replicates=replicates)
    out: Dict[str, object] = {"target": "symmetric-write", "status": derivation.get("status"),
                              "word": derivation.get("word"), "law_constant": derivation.get("law_constant") or {}}
    if not str(derivation.get("status", "")).startswith("reached"):
        out.update(verdict="nothing to audit: the symmetric write is not reached",
                   obstruction=derivation.get("obstruction"))
        return out
    target, nf, where = symmetric_write_target(real, np.random.default_rng(int(rng.integers(2 ** 31))))
    if target is None:
        out.update(verdict="not audited: " + where)
        return out
    base = {"law": derivation.get("law") or {}, "constant": derivation.get("law_constant") or {},
            "n_traj": n_traj, "replicates": replicates}
    ledger = []
    for name, fn in CHECKS.items():
        if checks and name not in checks:
            continue
        try:
            item = fn(real, target, nf, rng, base)
        except (ValueError, np.linalg.LinAlgError) as e:
            item = _item(name, "", "", NOT_TESTABLE, note=f"{type(e).__name__}: {str(e)[:120]}")
        ledger.append({"check": name, **item})
    out.update(measured_on=where, ledger=ledger, verdict=verdict(ledger, out["law_constant"]),
               qualifications=[x["error"] for x in ledger if x["status"] == QUALIFIED])
    return out


def verdict(ledger: List[Dict[str, object]], law_constant: Dict[str, object]) -> str:
    """Errors found, if any check confirms one; otherwise ready for review, but only if the law was measured: a write
    whose law was not measured cannot be reviewed, whatever the other checks say."""
    confirmed = [x for x in ledger if x["status"] == CONFIRMED]
    if confirmed:
        return "errors found: " + "; ".join(x["error"] for x in confirmed)
    if (law_constant or {}).get("constant") is None:
        return "not testable: the law of the write was not measured"
    return "ready for review"


def markdown(name: str, res: Dict[str, object]) -> str:
    lines = [f"# Audit of the symmetric-write law: {name}", ""]
    c = res.get("law_constant") or {}
    lines.append(f"Derivation: {res.get('status')} ({res.get('word') or '-'}); law constant "
                 + (f"{c['constant']:.3f} ± {c['stderr']:.3f}" if c.get("constant") is not None else "not measured")
                 + (f", measured on {res['measured_on']}" if res.get("measured_on") else "") + ".")
    lines += ["", f"Verdict: {res['verdict']}.", ""]
    if res.get("ledger"):
        lines += ["| check | possible error | test | result |", "| --- | --- | --- | --- |"]
        for x in res["ledger"]:
            lines.append(f"| {x['check']} | {x['error']} | {x['test']} | {x['status']} |")
        lines.append("")
        for x in res["ledger"]:
            if x["status"] in (CONFIRMED, QUALIFIED):
                ev = ", ".join(f"{k} = {v:.4g}" if isinstance(v, float) else f"{k} = {v}" for k, v in x["evidence"].items())
                lines.append(f"- {x['check']} ({x['status']}): {x['error']}. Evidence: {ev}.")
    lines += ["", "A law that passes the audit is ready for review by a person or against the literature; the audit does "
              "not certify it.", ""]
    return "\n".join(lines)
