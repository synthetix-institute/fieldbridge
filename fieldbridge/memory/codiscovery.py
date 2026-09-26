"""Co-discovery by construction: one mechanism derived in models from different fields, each by its own chain of
transformations, with a check that the chains end on the same mechanism.

Different fields often reach the same mechanism by different derivations. Across the arXiv corpus such convergences
are not more frequent than chance (derivation chains of the Hyperion V2.1 language). The constructor makes the
convergence deliberate: given a target mechanism, it derives the target in every realization, records the chain of
transformations as a word, and certifies the end by invariants that do not depend on the field.

Three targets are implemented (TARGETS):

  symmetric-write   a supercritical pitchfork, ds/dt = eps s - s^3 + h: at the write point a weak bias chooses one of
                    two symmetric states
  threshold-write   a fold (saddle-node), ds/dt = mu + s^2: the occupied state disappears at a threshold of the
                    control or of a write field, and the material switches to another state
  phase-locking     the Adler equation, dpsi/dt = Delta omega - K sin phi: a periodic modulation of the control fixes
                    the phase of a limit cycle, whose phase is otherwise a flat direction (phase_locking.py)

Transformations (the letters of a derivation word):

  S  symmetry      a linear symmetry of the drift reverses the critical mode, which forces a pitchfork (Module 3); for
                   the threshold write it is the reason why the control gives no threshold
  C  continuation  a stable state is followed along the control to the write point, where kappa reaches zero
  W  write field   a bounded field toward another stored state is raised until the occupied state disappears
  R  reduction     the drift is reduced to the critical direction on the slow manifold, a0 + a1 s + a2 s^2 + a3 s^3;
                   the other directions are eliminated adiabatically
  U  unfolding     a one-sided write (fold) is made symmetric by tuning a second material parameter to the cusp; the
                   control is replaced by the sweep through the cusp that carries no bias (followed by R again)
  K  canonical     the critical coordinate is rescaled: x = |a3|^(1/2) s gives eps x - x^3 + h, and x = a2 s gives
                   mu + x^2 with mu = a2 b (p - p_f), b the push of the parameter p along the critical mode
  L  law           the law of the canonical form is tested by simulating the realization itself: for the symmetric
                   write the accuracy of a swept write, P = Phi(pi^(1/4) h_s / (D_s^(1/2) r^(1/4))); for the threshold
                   write the delay of the switch under a sweep mu = r t, which crosses the static fold at
                   mu = |a1'| r^(2/3), with a1' = -1.01879 the first zero of the derivative of the Airy function; for
                   phase locking the relaxation rate inside and the slip frequency outside the locking range, which
                   give its half-width in units of K

For phase locking the letters S, C, R and K act on the same slots with the meanings given in TARGET_INFO: S is a
symmetry that maps the cycle onto itself a fraction of a period later and sets the locking ratio, C moves the control
into the oscillating range, R reduces the state to the phase of the cycle, and K averages over the drive.

Every letter acts on named slots of the identity ((Omega, Xi); C, R, P; A) (SLOTS below), so a derivation is a chain of
transformations of the mechanism in the same terms as the derivation chains of the corpus.

A derivation ends as reached or obstructed, with the obstruction named. Two invariants certify that the ends are the
same mechanism: the canonical form of the reduced drift at the write point (its units fix the leading coefficient, so
the content of the check is in the other terms: no even part for the symmetric write, no linear part at the fold), and
the constant of the law, estimated independently in every realization from a simulation of the full realization.

Realizations are compared by their derivation class: the letters with the kind of symmetry and of each reduction,
without the names of the variables. A new target is a function with the signature of derive_symmetric_write, added
to TARGETS with its description in TARGET_INFO.
"""
from __future__ import annotations

from math import pi, sqrt
from statistics import NormalDist
from typing import Dict, List, Optional, Tuple

import numpy as np

from . import analysis as an
from .construct import along_path, cancel_asymmetry, obstruction, swept_write_check
from .identity import Realization

LAW_CONSTANT = pi ** 0.25
DELAY_CONSTANT = 1.0187929716474710  # -a1', the first zero of Ai' (checked against scipy.special.ai_zeros in the tests)
DELAY_RATES = tuple(float(r) for r in 10.0 ** -np.arange(3.0, 7.01, 0.5))
DELAY_SOLVER = "LSODA"  # switches to an implicit method where the transverse modes are stiff (rotor patches)
LETTERS = {"S": "symmetry", "C": "continuation", "W": "write field", "R": "reduction", "U": "unfolding",
           "K": "canonical form", "L": "law"}
# the slots of the identity ((Omega, Xi); C, R, P; A) on which each transformation acts, and how
SLOTS = {"S": ("Xi, Omega", "a linear map of the carrier that leaves the drift invariant and reverses the critical mode"),
         "C": ("A", "the control is moved to the value where the linear rate kappa of the critical mode vanishes"),
         "R": ("Xi, Omega", "the carrier is reduced to the critical line; the drift to a0 + a1 s + a2 s^2 + a3 s^3"),
         "U": ("A, P", "a second material parameter is set to the cusp; the protocol sweeps through the cusp"),
         "W": ("P", "a bounded field toward another stored state is raised until the occupied state disappears"),
         "K": ("Xi", "the coordinate of the critical line is rescaled to the canonical form"),
         "L": ("P, R", "the law of the canonical form is tested by simulating the realization")}
TARGET_INFO = {
    "symmetric-write": {
        "title": "the symmetric write", "letters": "SCRUKL", "canonical": "-x^3", "check": "even part",
        "target": "symmetric write (supercritical pitchfork, ds/dt = eps s - s^3 + h)",
        "law": "accuracy of a swept write, P = Phi(pi^(1/4) h_s / (D_s^(1/2) r^(1/4)))",
        "constant_name": "pi^(1/4)", "constant": LAW_CONSTANT,
        "K": "x = |a3|^(1/2) s, so that the reduced drift is eps x - x^3 + h",
        "L": "a swept write with a bias is simulated; the observable is the sign of the critical coordinate"},
    "threshold-write": {
        "title": "the threshold write", "letters": "SCWRKL", "canonical": "x^2", "check": "linear part",
        "target": "threshold write (fold, ds/dt = mu + s^2)",
        "law": "delay of the switch under a sweep mu = r t: the static fold is crossed at mu = |a1'| r^(2/3)",
        "constant_name": "|a1'|", "constant": DELAY_CONSTANT,
        "K": "x = a2 s and mu = a2 b (p - p_f), so that the reduced drift is mu + x^2",
        "L": "the parameter is swept through the fold without noise at rates r = 1e-3 to 1e-7 (canonical units); "
             "the delay of the crossing is extrapolated to r = 0"},
    "phase-locking": {
        "title": "phase locking", "letters": "SCRKL", "canonical": "nu - sin phi", "check": "averaged drift",
        "target": "phase locking (Adler equation, dpsi/dt = Delta omega - K sin phi)",
        "law": "relaxation rate inside and slip frequency outside the locking range: "
               "lambda / (n K) = (1 - nu^2)^(1/2), Omega / K = (nu^2 - 1)^(1/2)",
        "constant_name": "1 (half-width in units of K)", "constant": 1.0,
        "S": "a symmetry of the drift that holds for every value of the control maps the cycle onto itself 1/m of a "
             "period later: the drive acts through the harmonics m, 2m, ... and locks at m:1 with m locked phases",
        "C": "the control is moved into the range where the realization oscillates",
        "R": "the carrier is reduced to the phase of the cycle; the phase response Z is the periodic solution of the "
             "adjoint equation, normalized by Z . F = omega_0",
        "K": "averaging over a drive p1 + eps cos(omega_f t), omega_f ~ n omega_0: dpsi/dt = Delta omega - K sin phi, "
             "phi = n psi - phi_n - pi/2, K = eps |Z_n| / 2",
        "L": "the driven realization is simulated at seven detunings and two amplitudes; the half-width of the locking "
             "range in units of K is extrapolated to zero amplitude"},
}
for _key, _check, _point in (("symmetric-write", "even_part", "write point"), ("threshold-write", "linear_part",
                                                                                  "write point"),
                             ("phase-locking", "deviation", "operating point")):
    TARGET_INFO[_key].update(check_key=_check, point=_point)
TARGET_INFO["phase-locking"]["examples"] = ("examples/memory/oscillators", "examples/memory")
SHORT_OBSTRUCTIONS = (("no control parameter", "no control parameter"),
                      ("rate of time", "the control only rescales time"),
                      ("multiplies the whole drift", "the control only rescales the drift"),
                      ("no write point", "no write point along the control"),
                      ("no stable state", "no stable state: the realization oscillates"),
                      ("single stable state", "a single stable state: nothing to switch to"),
                      ("exchange of stability", "the field leaves the state an equilibrium: no fold"),
                      ("Hopf", "Hopf: an oscillation, not a write"),
                      ("subcritical", "subcritical: jump to a distant branch"),
                      ("fold", "fold; no cusp with a supercritical sweep"),
                      ("neutral cycles", "a family of neutral cycles: no isolated phase"),
                      ("no oscillation", "no oscillation in the control range"),
                      ("does not act on the phase", "the drive does not act on the phase"),
                      ("only through harmonics above", "the drive acts only through high harmonics"))


def _field(real: Realization) -> str:
    spec = getattr(real, "spec", None) or {}
    return str(spec.get("field", "unspecified"))


def reversing_symmetry(real: Realization, rng, mode: np.ndarray) -> Optional[Dict[str, object]]:
    """A linear symmetry g of the drift with g(mode) = -mode: it forbids even terms along the mode."""
    from .predict import symmetries
    mode = np.asarray(mode, float)
    for g in symmetries(real, rng):
        if "permutation" not in g:
            continue
        P = np.eye(mode.size)[list(g["permutation"])] * np.asarray(g["signs"], float)[:, None]
        # the drift is equivariant under q -> sgn * q[perm]; the image of the mode is sgn * mode[perm]
        image = np.asarray(g["signs"], float) * mode[list(g["permutation"])]
        if np.allclose(image, -mode, atol=1e-6):
            return {"permutation": g["permutation"], "signs": g["signs"], "order": g["order"],
                    "matrix": P.tolist()}
    return None


def describe_symmetry(real: Realization, sym: Dict[str, object]) -> str:
    """A symmetry in words, in the variables of the realization."""
    names = list(getattr(real, "spec", {}).get("carrier", {}).get("variables", [])) or \
        [f"q{i}" for i in range(real.carrier.dim)]
    perm, sgn = list(sym["permutation"]), list(sym["signs"])
    if all(p == i for i, p in enumerate(perm)):
        flipped = [names[i] for i, g in enumerate(sgn) if g < 0]
        return "reflection " + ", ".join(f"{v} -> -{v}" for v in flipped)
    if sym["order"] == 2 and all(g > 0 for g in sgn):
        pairs = sorted({tuple(sorted((i, p))) for i, p in enumerate(perm) if i != p})
        return "exchange " + ", ".join(f"{names[i]} <-> {names[j]}" for i, j in pairs)
    return f"cyclic permutation of order {sym['order']}"


def route_string(steps: List[Dict[str, object]]) -> str:
    parts = []
    for st in steps:
        L = st["letter"]
        if L == "S":
            parts.append(f"S[{st.get('meaning', 'symmetry')}]")
        elif L == "C":
            parts.append(f"C[{st['control']}" + (f" = {st['value']:.4g}" if "value" in st else "") + "]")
        elif L == "R":
            k = st.get("eliminated_directions")
            parts.append(f"R[{st.get('kind', '')}" + (f"; {k} direction{'s' if k > 1 else ''} eliminated" if k else "")
                         + "]")
        elif L == "U":
            parts.append(f"U[{st['second']} to the cusp]")
        elif L == "W":
            parts.append(f"W[field toward stored state {st['to_state']}]" if st.get("to_state") is not None else "W")
        elif L == "K" and st.get("ratio"):
            parts.append(f"K[{st['ratio']}:1]")
        else:
            parts.append(L)
    return " > ".join(parts)


def canonical_form(nf: Dict) -> Dict[str, object]:
    """The reduced drift at the write point in the canonical variable x = |a3|^(1/2) s, g = |a3|^(1/2) ds/dt.

    The coefficients are local: a polynomial of degree 5 is fitted over the sampled window, so that the fifth-order
    term does not bias the cubic one. In these units g = sum c_n x^n with c_n = a_n |a3|^((1-n)/2), and c_3 = -1 for
    a supercritical pitchfork. The even part (c_0, c_2, c_4) must vanish for a symmetric write; the linear part
    vanishes at the write point. Both are reported relative to the cubic term at the edge of the window."""
    s = np.asarray(nf.get("s", []), float)
    f = np.asarray(nf.get("f", []), float)
    if s.size < 7:
        return {}
    a = np.polyfit(s, f, 5)[::-1]
    if a[3] == 0.0:
        return {}
    k = sqrt(abs(a[3]))
    c = [float(a[n] * k ** (1 - n)) for n in range(6)]
    x, g = k * s, k * f
    xm = float(np.max(np.abs(x)))
    return {"scale": k, "a3_local": float(a[3]), "coefficients": c, "cubic_sign": int(np.sign(c[3])),
            "even_part": float(np.max(np.abs(g + g[::-1])) / 2 / xm ** 3),
            "linear_part": float(abs(c[1]) / xm ** 2), "quintic": c[5], "x": x.tolist(), "g": g.tolist(),
            "eps_per_unit_control": float(-nf.get("kappa_slope", 0.0))}


def derivation_class(steps: List[Dict[str, object]], coarse: bool = False) -> str:
    """The derivation without the names of the variables: letters, the kind of symmetry and of each reduction.
    With coarse=True the kind of symmetry and the side of a pitchfork are dropped: for the threshold write they
    only explain why the control gives no threshold."""
    parts = []
    for st in steps:
        L = st["letter"]
        if L == "S":
            m = str(st.get("meaning", ""))
            parts.append("S" if coarse else
                         "S(" + (f"cyclic, order {st['order']}" if m.startswith("cyclic") else m.split(" ")[0]) + ")")
        elif L == "R":
            kind = str(st.get("kind", ""))
            if "fold" in kind:
                parts.append("R(fold)")
            elif "pitchfork" in kind:
                parts.append("R(pitchfork)" if coarse or not kind.startswith("subcritical")
                             else "R(subcritical pitchfork)")
            else:
                parts.append(f"R({kind.split(':')[0]})")
        elif L == "K" and st.get("ratio"):
            parts.append(f"K({st['ratio']}:1)")
        else:
            parts.append(L)
    return " ".join(parts)


def _short(reason: str) -> str:
    return next((short for key, short in SHORT_OBSTRUCTIONS if key in reason), "other terms in the reduced drift")


def law_constant(law: Dict) -> Dict[str, float]:
    """The constant of the swept-write law inferred from a measured accuracy: Phi^-1(P) D_s^(1/2) r^(1/4) / h_s."""
    if not law or "measured" not in law:
        return {}
    p, se = float(law["measured"]), float(law.get("stderr", 0.0))
    p_c = min(max(p, 1e-4), 1 - 1e-4)
    nd = NormalDist()
    z = nd.inv_cdf(p_c)
    scale = sqrt(law["D_along_mode"]) * law["rate"] ** 0.25 / law["h_along_mode"]
    dz = se / max(nd.pdf(z), 1e-12)
    return {"constant": float(z * scale), "stderr": float(dz * scale), "expected": LAW_CONSTANT}


def derive_symmetric_write(real: Realization, rng, check_law: bool = True, n_traj: int = 400) -> Dict[str, object]:
    """Derive the symmetric write in one realization; returns the word, the steps and the certificate."""
    word: List[str] = []
    steps: List[Dict[str, object]] = []
    out: Dict[str, object] = {"name": real.name, "field": _field(real), "status": "obstructed"}
    if not real.control:
        out["obstruction"] = "no control parameter: nothing can be swept through a write point"
        out["obstruction_short"] = _short(out["obstruction"])
        return {**out, "word": "", "route": "", "class": "", "steps": steps}
    if an.is_scale_control(real, rng):
        out["obstruction"] = ("the control multiplies the whole drift: it only rescales the landscape, so no state "
                              "changes along it; a state is written by a field instead")
        out["obstruction_short"] = _short(out["obstruction"])
        return {**out, "word": "", "route": "", "class": "", "steps": steps}
    events = an.locate_writes(real, rng, n_starts=24 if real.carrier.dim <= 4 else 12)
    if not events:
        out["obstruction"] = "no write point along the control"
        out["obstruction_short"] = _short(out["obstruction"])
        return {**out, "word": "", "route": "", "class": "", "steps": steps}
    nfs = []
    for ev in events:
        nf = an.normal_form(real, ev["q"], real.control, ev["v"])
        nf["obstruction"] = obstruction(nf)
        nfs.append(nf)
    word += ["C", "R"]
    steps.append({"letter": "C", "control": real.control,
                  "write_points": [{"value": float(nf["value"]), "kind": nf["kind"].split(":")[0]} for nf in nfs]})
    pick = next((nf for nf in nfs if nf["obstruction"]["realizes_generator"]), None)
    target_real, target_nf = real, pick
    if pick is not None:
        steps.append({"letter": "R", "value": float(pick["value"]), "a0": pick["a0"], "a2": pick["a2"],
                      "a3": pick["a3"], "eliminated_directions": real.carrier.dim - 1,
                      "kind": pick["kind"].split(":")[0]})
        sym = reversing_symmetry(real, rng, np.asarray(pick["mode"]))
        if sym is not None:
            word.insert(0, "S")
            steps.insert(0, {"letter": "S", "meaning": describe_symmetry(real, sym),
                             **{k: sym[k] for k in ("permutation", "signs", "order")}})
    else:
        first = nfs[0]
        steps.append({"letter": "R", "value": float(first["value"]), "kind": first["kind"].split(":")[0],
                      "eliminated_directions": real.carrier.dim - 1,
                      **({"a2": first["a2"], "a3": first["a3"]} if not first.get("hopf") else
                         {"omega": first.get("omega")})})
        if not first.get("hopf") and "mode" in first:
            sym = reversing_symmetry(real, rng, np.asarray(first["mode"]))
            if sym is not None:
                word.insert(0, "S")
                steps.insert(0, {"letter": "S", "meaning": describe_symmetry(real, sym),
                                 **{k: sym[k] for k in ("permutation", "signs", "order")}})
        fold = next((nf for nf in nfs if (nf["obstruction"]["asymmetric"] or nf["obstruction"]["biased"])
                     and not nf.get("hopf")), None)
        unfolded = None
        if fold is not None:
            for second in real.material:
                with np.errstate(all="ignore"):
                    try:
                        c = cancel_asymmetry(real, fold, second)
                    except np.linalg.LinAlgError:
                        continue
                if c.get("success") and c.get("obstruction_along_sweep", {}).get("realizes_generator"):
                    unfolded = c
                    break
        if unfolded is None:
            ob = nfs[0]["obstruction"]
            if nfs[0].get("hopf"):
                reason = "a complex pair crosses: the state oscillates (Hopf), a phase memory rather than a write"
            elif ob.get("subcritical"):
                reason = "positive cubic term: the write is symmetric but jumps to a distant branch (subcritical)"
            elif fold is not None:
                reason = ("one-sided write (fold) and no material parameter reaches a cusp with a supercritical "
                          "sweep")
            else:
                reason = "the reduced drift differs from the pitchfork: " + ", ".join(ob.get("terms", {}))
            out.update(obstruction=reason, obstruction_short=_short(reason))
            return {**out, "word": "".join(word), "route": route_string(steps), "class": derivation_class(steps),
                    "steps": steps}
        path = along_path(real, {real.control: unfolded["control"], unfolded["second"]: unfolded["value"]},
                          unfolded["direction"])
        nf_path = an.normal_form(path, np.asarray(unfolded["state"]), "t", 0.0)
        nf_path["obstruction"] = obstruction(nf_path)
        word += ["U", "R"]
        steps.append({"letter": "U", "second": unfolded["second"], "cusp": {unfolded["second"]: unfolded["value"],
                                                                           real.control: unfolded["control"]},
                      "sweep_direction": unfolded["direction"]})
        steps.append({"letter": "R", "value": 0.0, "a0": nf_path["a0"], "a2": nf_path["a2"], "a3": nf_path["a3"],
                      "kind": nf_path["kind"].split(":")[0], "along": "the sweep through the cusp",
                      "eliminated_directions": real.carrier.dim - 1})
        target_real, target_nf = path, nf_path
        out["status"] = "reached after unfolding"
    canon = canonical_form(target_nf)
    word.append("K")
    steps.append({"letter": "K", "scale": canon.get("scale"), "a3_local": canon.get("a3_local"),
                  "cubic_sign": canon.get("cubic_sign"),
                  "even_part": canon.get("even_part"), "linear_part": canon.get("linear_part"),
                  "quintic": canon.get("quintic")})
    law, const = {}, {}
    if check_law:
        law = swept_write_check(target_real, target_nf, rng, n_traj=n_traj)
        const = law_constant(law)
        if law:
            word.append("L")
            steps.append({"letter": "L", "predicted": law["predicted"], "measured": law["measured"],
                          "stderr": law["stderr"], "gamma": law["gamma"], **{f"law_{k}": v for k, v in const.items()}})
    if out["status"] == "obstructed":
        out["status"] = "reached"
    return {**out, "word": "".join(word), "route": route_string(steps), "class": derivation_class(steps),
            "steps": steps, "canonical": canon, "law": law, "law_constant": const,
            "write_point": {"param": target_nf["param"], "value": float(target_nf["value"])}}


# ------------------------------------------------------------------------------------------------ threshold write
def _shift(real: Realization, q: np.ndarray, d: np.ndarray) -> np.ndarray:
    return real.carrier.wrap(q + d) if real.carrier.kind == "torus" else q + d


def _is_fold(nf: Dict) -> bool:
    ob = nf.get("obstruction") or obstruction(nf)
    return bool(not nf.get("hopf") and ob["asymmetric"] and ob["biased"])


def refine_fold(real: Realization, param: str, q, v: float) -> Optional[Tuple[np.ndarray, float]]:
    """The fold point from F = 0 and det J = 0, started from a state next to it on its branch. The determinant is
    divided by the product of the other eigenvalues, so that its residual is the critical eigenvalue."""
    from scipy.optimize import root
    q = np.asarray(q, float)
    n = q.size
    lam = np.sort(np.abs(np.linalg.eigvals(an.jacobian(real, q, {param: v}))))
    others = max(float(np.prod(lam[1:])), 1e-300) if n > 1 else 1.0

    def eqs(y):
        over = {param: y[n]}
        with np.errstate(all="ignore"):
            res = np.concatenate([real.F(y[:n], **over), [np.linalg.det(an.jacobian(real, y[:n], over)) / others]])
        return res if np.all(np.isfinite(res)) else np.full(n + 1, 1e3)

    sol = root(eqs, np.concatenate([q, [v]]), method="hybr", options={"xtol": 1e-13})
    if not np.all(np.isfinite(sol.x)) or np.max(np.abs(sol.fun)) > 1e-8 \
            or real.carrier.distance(sol.x[:n], q) > 0.25 * an.base_length(real):
        return None
    return sol.x[:n], float(sol.x[n])


def canonical_fold(real: Realization, param: str, q_f: np.ndarray, p_f: float, nf: Dict) -> Dict[str, object]:
    """The reduced drift at the fold in the canonical variable x = a2 s, g = a2 ds/dt: mu + x^2 + c3 x^3 + ...

    a2 is half the second derivative of the drift along the critical mode (the slow manifold enters only at third
    order), b the push of the parameter along the mode, mu = a2 b (p - p_f). The checks are the linear part at the
    fold relative to the slowest transverse rate, and c2 = 1 for the quadratic coefficient of the slow-manifold drift
    fitted near the fold. c3 is the leading correction and differs between realizations."""
    over = {param: p_f}
    vec, w = np.asarray(nf["mode"], float), np.asarray(nf["left_mode"], float)
    h = 1e-3 * an.base_length(real)
    d2 = (real.F(_shift(real, q_f, h * vec), **over) + real.F(_shift(real, q_f, -h * vec), **over)
          - 2 * real.F(q_f, **over)) / h ** 2
    a2 = 0.5 * float(w @ d2)
    J = an.jacobian(real, q_f, over)
    a1 = float(w @ J @ vec)
    rates = np.sort(np.abs(np.linalg.eigvals(J).real))
    half = float(nf.get("half", 0.1 * an.base_length(real)))
    ref = float(rates[1]) if rates.size > 1 and rates[1] > 0 else abs(a2) * half
    s_near = np.linspace(-1.0, 1.0, 21) * half / 8
    f_near, _ = an._slow_manifold_drift(real, q_f, vec, w, s_near, over)
    loc = np.polyfit(s_near, f_near, 4)[::-1]
    s, f = np.asarray(nf.get("s", []), float), np.asarray(nf.get("f", []), float)
    return {"a2": a2, "b": float(nf["bias_along_mode"]), "mu_per_unit_parameter": a2 * float(nf["bias_along_mode"]),
            "a1": a1, "linear_part": abs(a1) / max(ref, 1e-300), "c2": float(loc[2] / a2),
            "cubic": float(loc[3] / a2 ** 2), "x": (a2 * s).tolist(), "g": (a2 * f).tolist(),
            "mode": vec.tolist(), "left_mode": w.tolist()}


def fold_delay_law(real: Realization, param: str, q_f: np.ndarray, p_f: float, canon: Dict,
                   rates=DELAY_RATES, tau: float = 5.0) -> Dict[str, object]:
    """Sweep the parameter through the fold without noise, p = p_f + r t / (a2 b), and record the time t0 at which the
    state crosses the position of the static fold (s = 0). In canonical units mu = r t, and the crossing is at
    mu = |a1'| r^(2/3) plus corrections in powers of r^(1/3). The constant t0 r^(1/3) is extrapolated to r = 0 by a
    quadratic in r^(1/3) over the five slowest rates; its uncertainty is the difference from a linear fit to the
    three slowest. Each sweep starts on the occupied branch at mu = -tau r^(2/3)."""
    from scipy.integrate import solve_ivp
    vec, w = np.asarray(canon["mode"], float), np.asarray(canon["left_mode"], float)
    a2, k = canon["a2"], canon["mu_per_unit_parameter"]
    rows = []
    for r in rates:
        rho = r / k
        t_start = -tau * r ** (-1 / 3)
        p_start = p_f + rho * t_start
        q0, ok = an.newton(real, _shift(real, q_f, (-sqrt(tau) * r ** (1 / 3) / a2) * vec), {param: p_start})
        s0 = float(w @ real.carrier.diff(q_f, q0))
        if not ok or a2 * s0 >= 0 or an.n_unstable(an.kappa_spectrum(real, q0, {param: p_start})) > 0:
            continue  # at this rate the start lies beyond the branch of the occupied state

        def rhs(t, q, rho=rho):
            return real.F(q, **{param: p_f + rho * t})

        def cross(t, q):
            return float(w @ real.carrier.diff(q_f, q))

        cross.terminal, cross.direction = True, float(np.sign(a2))
        sol = solve_ivp(rhs, (t_start, 8.0 * r ** (-1 / 3)), q0, method=DELAY_SOLVER, rtol=1e-12, atol=1e-14,
                        events=cross)
        if sol.t_events[0].size:
            rows.append({"rate": float(r), "parameter_rate": float(rho),
                         "constant": float(sol.t_events[0][0] * r ** (1 / 3))})
    out: Dict[str, object] = {"rows": rows, "expected": DELAY_CONSTANT}
    if len(rows) >= 5:
        e = np.array([row["rate"] for row in rows]) ** (1 / 3)
        c = np.array([row["constant"] for row in rows])
        i5, i3 = np.argsort(e)[:5], np.argsort(e)[:3]
        quad, lin = np.polyfit(e[i5], c[i5], 2), np.polyfit(e[i3], c[i3], 1)
        out.update(constant=float(quad[2]), stderr=float(abs(quad[2] - lin[1])), correction_slope=float(quad[1]),
                   fit=[float(quad[2]), float(quad[1]), float(quad[0])])
    return out


def _field_fold(real: Realization, states: List[np.ndarray], max_pairs: int = 12):
    """A bounded write field toward another stored state, raised until the occupied state disappears. Pairs are
    tried in order of the force the field exerts on the occupied state: a field that exerts none leaves the state an
    equilibrium at every strength, and the state is then lost by an exchange of stability, not by a fold."""
    pairs = sorted(((float(np.linalg.norm(real.write_force(states[i], states[j], 1.0))), i, j)
                    for i in range(len(states)) for j in range(len(states)) if i != j), key=lambda t: -t[0])
    tried = []
    for force, i, j in pairs[:max_pairs]:
        if force < 1e-8:
            tried.append({"from": i, "to": j, "result": "the field exerts no force on the occupied state"})
            continue
        h_max, ev = 3.0 * an.restoring_threshold(real, states[i], states[j]), None
        for _ in range(14):
            rw = an.with_write_field(real, states[j], (0.0, h_max))
            _, ev = an.track(rw, "h", states[i], 0.0, h_max)
            if ev is not None:
                break
            h_max *= 4.0
        if ev is None:
            tried.append({"from": i, "to": j, "result": "the occupied state persists"})
            continue
        nf = an.normal_form(rw, np.asarray(ev["q"], float), "h", ev["v"])
        nf["obstruction"] = obstruction(nf)
        tried.append({"from": i, "to": j, "h": float(ev["v"]), "result": nf["kind"].split(":")[0]})
        if _is_fold(nf):
            return rw, ev, nf, i, j, tried
    return None, None, None, None, None, tried


def derive_threshold_write(real: Realization, rng, check_law: bool = True, n_traj: int = 0) -> Dict[str, object]:
    """Derive the threshold write in one realization: first along the control (C, R); if the control gives no fold,
    by a write field toward another stored state (W, R). The fold is then refined, rescaled (K) and its delay law
    tested (L). n_traj is not used: the law is deterministic."""
    word: List[str] = []
    steps: List[Dict[str, object]] = []
    notes: List[str] = []
    out: Dict[str, object] = {"name": real.name, "field": _field(real), "status": "obstructed"}
    fold, route = None, None
    scale = bool(real.control) and an.is_scale_control(real, rng)
    if real.control and not scale:
        events = an.locate_writes(real, rng, n_starts=24 if real.carrier.dim <= 4 else 12)
        found = []
        for ev in events:
            nf = an.normal_form(real, ev["q"], real.control, ev["v"])
            nf["obstruction"] = obstruction(nf)
            found.append((ev, nf))
        if found:
            word += ["C", "R"]
            steps.append({"letter": "C", "control": real.control,
                          "write_points": [{"value": float(nf["value"]), "kind": nf["kind"].split(":")[0]}
                                           for _, nf in found]})
            pick = next(((ev, nf) for ev, nf in found if _is_fold(nf)), None)
            ev, nf = pick if pick is not None else found[0]
            steps.append({"letter": "R", "value": float(nf["value"]), "kind": nf["kind"].split(":")[0],
                          "eliminated_directions": real.carrier.dim - 1})
            if pick is not None:
                fold, route = (real, real.control, ev, nf), "control"
            else:
                sym = None if nf.get("hopf") or "mode" not in nf else reversing_symmetry(real, rng, np.asarray(nf["mode"]))
                if sym is not None:
                    word.insert(0, "S")
                    steps.insert(0, {"letter": "S", "meaning": describe_symmetry(real, sym),
                                     **{k: sym[k] for k in ("permutation", "signs", "order")}})
                notes.append(f"along the control the write point is a {nf['kind'].split(':')[0]}"
                             + ("; the symmetry forbids a quadratic term" if sym is not None else ""))
        else:
            notes.append("no write point along the control")
    else:
        notes.append("the control multiplies the whole drift and changes no state" if scale else "no control parameter")
    if fold is None:
        states = an.stored_states(real, rng, n_starts=48 if real.carrier.dim <= 8 else 24)[0]
        if len(states) < 2:
            reason = ("no stable state (the realization oscillates): there is no state to leave" if not states else
                      "a single stable state: there is no second state to switch to")
            out.update(obstruction=reason, obstruction_short=_short(reason))
            return {**out, "word": "".join(word), "route": route_string(steps),
                    "class": derivation_class(steps, coarse=True), "steps": steps, "notes": notes}
        rw, ev, nf, i, j, tried = _field_fold(real, states)
        word.append("W")
        steps.append({"letter": "W", "from_state": i, "to_state": j, "states": len(states), "tried": tried,
                      "coercive_field": None if ev is None else float(ev["v"])})
        if nf is None:
            reason = ("along every write field the occupied state is lost by an exchange of stability or not at all, "
                      "not by a fold")
            out.update(obstruction=reason, obstruction_short=_short(reason))
            return {**out, "word": "".join(word), "route": route_string(steps),
                    "class": derivation_class(steps, coarse=True), "steps": steps, "notes": notes}
        word.append("R")
        steps.append({"letter": "R", "value": float(ev["v"]), "kind": nf["kind"].split(":")[0],
                      "eliminated_directions": real.carrier.dim - 1})
        fold, route = (rw, "h", ev, nf), "field"
    target_real, param, ev, nf = fold
    refined = refine_fold(target_real, param, ev["q"], ev["v"])
    q_f, p_f = refined if refined is not None else (np.asarray(ev["q"], float), float(ev["v"]))
    if refined is None:
        notes.append("the fold was located by continuation only; the law was not tested")
    canon = canonical_fold(target_real, param, q_f, p_f, an.normal_form(target_real, q_f, param, p_f))
    word.append("K")
    steps.append({"letter": "K", **{k: canon[k] for k in ("a2", "b", "mu_per_unit_parameter", "linear_part", "c2",
                                                          "cubic")}})
    law, const = {}, {}
    if check_law and refined is not None:
        law = fold_delay_law(target_real, param, q_f, p_f, canon)
        w = np.asarray(canon["left_mode"], float)
        law["noise_rate"] = float(canon["a2"] ** 2 * real.noise * float(w @ w))
        if "constant" in law:
            const = {"constant": law["constant"], "stderr": law["stderr"], "expected": DELAY_CONSTANT}
            word.append("L")
            steps.append({"letter": "L", "law_constant": law["constant"], "law_stderr": law["stderr"],
                          "correction_slope": law["correction_slope"], "rates": len(law["rows"])})
    out["status"] = "reached" if route == "control" else "reached by a write field"
    return {**out, "word": "".join(word), "route": route_string(steps), "class": derivation_class(steps, coarse=True),
            "steps": steps, "canonical": canon, "law": law, "law_constant": const,
            "write_point": {"param": param, "value": float(p_f)}, "notes": notes}


# ------------------------------------------------------------------------------------------------ phase locking
def derive_phase_locking(real: Realization, rng, check_law: bool = True, n_traj: int = 0) -> Dict[str, object]:
    """Derive phase locking in one realization (phase_locking.derive). n_traj is not used: the law is deterministic."""
    from .phase_locking import derive
    return derive(real, rng, check_law=check_law)


# ------------------------------------------------------------------------------------------------ comparison
TARGETS = {"symmetric-write": derive_symmetric_write, "threshold-write": derive_threshold_write,
           "phase-locking": derive_phase_locking}


def codiscover(reals: List[Realization], rng, target: str = "symmetric-write", check_law: bool = True,
               n_traj: int = 400) -> Dict[str, object]:
    """Derive the target in every realization and compare the ends: derivation classes, fields and invariants."""
    info = TARGET_INFO[target]
    rows = [TARGETS[target](r, rng, check_law=check_law, n_traj=n_traj) for r in reals]
    reached = [r for r in rows if r["status"].startswith("reached")]
    classes: Dict[str, List[str]] = {}
    for r in reached:
        classes.setdefault(r["class"], []).append(r["name"])
    summary = {"target_key": target, "target": info["target"], "reached": len(reached),
               "obstructed": len(rows) - len(reached), "fields_reached": sorted({r["field"] for r in reached}),
               "distinct_words": sorted({r["word"] for r in reached}), "derivation_classes": classes,
               "law_constant_expected": info["constant"], "constant_name": info["constant_name"]}
    canon = [r["canonical"] for r in reached if r.get("canonical")]
    if target == "symmetric-write":
        summary.update(even_part_max=max((c.get("even_part", np.nan) for c in canon), default=None),
                       cubic_signs=sorted({c.get("cubic_sign") for c in canon}),
                       trajectories_per_realization=n_traj if check_law else 0)
        floor = 1e-9
    elif target == "threshold-write":
        summary.update(linear_part_max=max((c["linear_part"] for c in canon), default=None),
                       c2_range=[min((c["c2"] for c in canon), default=None), max((c["c2"] for c in canon), default=None)],
                       rates=list(DELAY_RATES) if check_law else [])
        floor = 1e-5  # integration tolerance and extrapolation
    else:
        from .phase_locking import K_REL
        ratios: Dict[str, List[str]] = {}
        for r in reached:
            ratios.setdefault(f"{r['canonical']['ratio']}:1", []).append(r["name"])
        summary.update(ratios=ratios, deviation_max=max((c["deviation"] for c in canon), default=None),
                       amplitudes=list(K_REL) if check_law else [])
        floor = 1e-4  # integration tolerance and extrapolation
    consts = [r["law_constant"] for r in reached if r.get("law_constant")]
    if consts:
        c = np.array([k["constant"] for k in consts])
        e = np.array([max(k["stderr"], floor) for k in consts])
        wmean = float(np.sum(c / e ** 2) / np.sum(1 / e ** 2))
        summary.update(law_constant_mean=wmean, law_constant_stderr=float(1 / np.sqrt(np.sum(1 / e ** 2))),
                       law_constant_chi2=float(np.sum(((c - info["constant"]) / e) ** 2)), law_constant_dof=int(c.size),
                       law_constant_max_deviation=float(np.max(np.abs(c - info["constant"]))))
    return {"summary": summary, "rows": rows}


def _sci(v: Optional[float]) -> str:
    if v is None or not np.isfinite(v):
        return "-"
    return "0" if v == 0 else f"{v:.1e}"


def markdown(report: Dict[str, object]) -> str:
    s = report["summary"]
    key = s.get("target_key", "symmetric-write")
    info = TARGET_INFO[key]
    classes = s["derivation_classes"]
    digits = {"symmetric-write": 3, "threshold-write": 5}.get(key, 4)
    lines = [f"# Co-discovery by construction: {info['title']}", "",
             f"Target: {s['target']}. Reached in {s['reached']} realizations from {len(s['fields_reached'])} fields "
             f"({', '.join(s['fields_reached'])}) by {len(classes)} classes of derivation; obstructed in "
             f"{s['obstructed']}.", "",
             f"| realization | field | derivation | status | {info['point']} | {info['check']} | law constant |",
             "|---|---|---|---|---|---:|---:|"]
    for r in report["rows"]:
        wp = r.get("write_point")
        wp_s = f"{wp['param']} = {round(wp['value'], 6):.6g}" if wp else "-"
        if wp and key == "phase-locking":
            wp_s += f"; {r['canonical']['ratio']}:1"
        check = r.get("canonical", {}).get(info["check_key"])
        lc = r.get("law_constant", {})
        lc_s = f"{lc['constant']:.{digits}f} ± {lc['stderr']:.{digits}f}" if lc else "-"
        status = r["status"] if r["status"] != "obstructed" else "obstructed: " + r["obstruction_short"]
        lines.append(f"| {r['name']} | {r['field']} | {r['word'] or '-'} | {status} | {wp_s} | {_sci(check)} | "
                     f"{lc_s} |")
    lines += ["", "## Classes of derivation", ""]
    for cls, names in classes.items():
        lines.append(f"- {cls}: {', '.join(names)}")
    lines += ["", "## Derivations", ""]
    for r in report["rows"]:
        if r["status"] != "obstructed":
            note = f" ({'; '.join(r['notes'])})" if r.get("notes") else ""
            lines.append(f"- {r['name']} ({r['field']}): {r['route']}{note}")
    lines += ["", "## Obstructions", ""]
    for r in report["rows"]:
        if r["status"] == "obstructed":
            note = f" ({'; '.join(r['notes'])})" if r.get("notes") else ""
            lines.append(f"- {r['name']} ({r['field']}): {r.get('route') or 'no derivation'}; {r['obstruction']}{note}")
    lines += ["", "## Invariants", ""]
    if key == "symmetric-write":
        lines.append(f"- Canonical form at the write point: cubic sign {', '.join(str(v) for v in s['cubic_signs'])} "
                     f"in every reached realization; largest even part {_sci(s['even_part_max'])} of the cubic term "
                     f"at the edge of the window.")
        if "law_constant_mean" in s:
            lines.append(f"- Constant of the swept-write law ({s['trajectories_per_realization']} trajectories per "
                         f"realization): {s['law_constant_mean']:.3f} ± {s['law_constant_stderr']:.3f} over the "
                         f"reached realizations; pi^(1/4) = {s['law_constant_expected']:.3f}; chi^2 = "
                         f"{s['law_constant_chi2']:.1f} for {s['law_constant_dof']} values.")
    elif key == "phase-locking":
        lines.append("- Locking ratios: " + "; ".join(f"{k} in {', '.join(v)}" for k, v in s["ratios"].items())
                     + ". A ratio above 1:1 comes from a symmetry that maps the cycle onto itself a fraction of a "
                       "period later (S).")
        if s.get("deviation_max") is not None:
            lines.append(f"- Canonical form: the phase gained over one period in the driven realization deviates from "
                         f"T K cos(n psi - phi_n) by at most {s['deviation_max']:.1e} of K at K = 0.01 min(omega_0, "
                         f"2 kappa); the deviation is of first order in the drive amplitude.")
        if "law_constant_mean" in s:
            lines.append(f"- Half-width of the locking range in units of K, extrapolated to zero amplitude from "
                         f"K = {' and '.join(f'{a:g}' for a in s['amplitudes'])} min(omega_0, 2 kappa): "
                         f"{s['law_constant_mean']:.4f} ± {s['law_constant_stderr']:.4f}; expected 1; largest "
                         f"deviation {s['law_constant_max_deviation']:.1e}.")
    else:
        lines.append(f"- Canonical form at the fold, mu + x^2: largest linear part {_sci(s['linear_part_max'])} of the "
                     f"slowest transverse rate; quadratic coefficient of the slow-manifold drift c2 = "
                     f"{s['c2_range'][0]:.4f} to {s['c2_range'][1]:.4f} (1 by the choice of units).")
        if "law_constant_mean" in s:
            lines.append(f"- Delay of the switch: crossing of the static fold at mu = C r^(2/3), extrapolated to r = 0 "
                         f"from {len(s['rates'])} rates per realization: C = {s['law_constant_mean']:.5f} ± "
                         f"{s['law_constant_stderr']:.5f}; |a1'| = {s['law_constant_expected']:.5f}; largest deviation "
                         f"{s['law_constant_max_deviation']:.1e}.")
    lines += ["", "## Letters and the slots they act on", "", "| letter | transformation | slots | action |",
              "|---|---|---|---|"]
    for k in info["letters"]:
        action = info.get(k, SLOTS[k][1])
        lines.append(f"| {k} | {LETTERS[k]} | {SLOTS[k][0]} | {action} |")
    return "\n".join(lines) + "\n"
