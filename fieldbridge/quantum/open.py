"""An open-system carrier and the quantum write: a parametric oscillator swept through its threshold with a bias.

Specification schema ``fieldbridge-quantum-open/1``:

    carrier      {"kind": "mode", "max_quanta": N}                       one bosonic mode, Fock truncation N
    parameters   {name: value}                                           real numbers
    hamiltonian  [{"coefficient": expr, "operator": "ad ad a a", "hc": bool, "imaginary": bool}]
                 expr is restricted arithmetic in the parameters and the protocol variables; "imaginary" multiplies
                 the coefficient by i; "hc" adds the Hermitian conjugate of the term
    dissipators  [{"rate": expr, "operator": "a"}]                      Lindblad terms rate * D[operator]
    protocol     {"sweep": {"parameter": "eps2", "from": 0, "to": 1.2, "rate": 0.01}, "hold": 10,
                  "bias": {"parameter": "h", "value": 0.05}, "initial": {"thermal": expr}}
    observable   {"sign_of": "a + ad"}                                   the branch is the sign of this quadrature

The derivation of the target ``quantum-write`` has the letters G (Gaussian stage: the generator is quadratic plus
linear loss while the amplitude is small, the sweep crosses the threshold, the amplified quadrature is the measured
one and the bias pushes it), K (the canonical linear equation x' = (eps(t) - kappa/2) x + h + sqrt(2D) xi with the
noise fixed by the loss and the temperature) and L (the law P = Phi(h I1 / sqrt(sigma0^2 + 2 D I2)) over the protocol,
compared with the exact evolution of the density operator). Where the law fails, the biased steady state and the
switching gap of the Lindbladian decide whether the write is an equilibrium write.

Nothing in a specification is evaluated: coefficients go through the restricted parser of the memory module.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple, Union

import numpy as np

from ..memory.spec import SpecError as MemorySpecError, parse_expression
from .carriers import Carrier, CarrierError

SCHEMA = "fieldbridge-quantum-open/1"
REQUIRED = ("name", "question", "assumptions", "provenance", "carrier", "hamiltonian", "dissipators", "protocol",
            "observable")
LETTERS = {"G": "Gaussian stage", "K": "canonical form", "L": "law"}
SLOTS = {"G": ("Omega, P, R", "the generator is quadratic plus linear loss while the amplitude is small; the sweep "
                              "crosses the threshold; the amplified quadrature is the measured one; the bias pushes it"),
         "K": ("Omega", "x' = (eps(t) - kappa/2) x + h + sqrt(2D) xi, with 2D = kappa (2 nbar + 1)/4 fixed by the "
                        "loss and the temperature"),
         "L": ("P, R", "the density operator is evolved exactly over the protocol and the probability of the favoured "
                       "state is compared with the law")}
REACHED = "reached"
MAX_QUANTA = 160
STEADY_MAX_DIM = 48
TOL_LAW = 2e-3          # the law is a calibration: the exact evolution has no sampling error
TOL_DIRECTION = 1e-6


class SpecError(ValueError):
    pass


@dataclass
class Realization:
    name: str
    field: str
    spec: Dict
    carrier: Carrier
    params: Dict[str, float]
    protocol: Dict
    hamiltonian: Callable[[float, float], np.ndarray]     # H(eps2, h)
    terms: List[Dict]
    dissipators: List[Tuple[float, np.ndarray, str]]
    quadrature: np.ndarray
    initial: np.ndarray


# ------------------------------------------------------------------------------------------------ specifications
def _mode_carrier(spec: Dict) -> Carrier:
    if spec.get("kind") != "mode":
        raise SpecError("the open carrier must be one bosonic mode: {\"kind\": \"mode\", \"max_quanta\": N}")
    nmax = int(spec.get("max_quanta", 0))
    if not 4 <= nmax <= MAX_QUANTA:
        raise SpecError(f"max_quanta must be between 4 and {MAX_QUANTA}")
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1).astype(complex)
    tokens = {"I": np.eye(nmax + 1, dtype=complex), "a": a, "ad": a.conj().T, "na": a.conj().T @ a}
    return Carrier("mode", nmax + 1, tokens, f"one bosonic mode, at most {nmax} quanta", spec)


def _expression(text, names: Dict[str, object]):
    try:
        return parse_expression(str(text), names)
    except MemorySpecError as err:
        raise SpecError(str(err)) from err


def _real_number(text, params: Dict[str, float]) -> float:
    import sympy as sp
    names = {k: sp.Symbol(k) for k in params}
    val = complex(_expression(text, names).subs({names[k]: v for k, v in params.items()}).evalf())
    if abs(val.imag) > 1e-12:
        raise SpecError(f"'{text}' is not real")
    return float(val.real)


def _operator(carrier: Carrier, text: str) -> np.ndarray:
    try:
        return carrier.operator(text)
    except CarrierError as err:
        raise SpecError(str(err)) from err


def load(source: Union[str, Path, Dict]) -> Realization:
    import sympy as sp
    spec = source if isinstance(source, dict) else json.loads(Path(source).read_text(encoding="utf-8"))
    if spec.get("schema") != SCHEMA:
        raise SpecError(f"schema must be {SCHEMA}")
    missing = [k for k in REQUIRED if not spec.get(k)]
    if missing:
        raise SpecError("missing: " + ", ".join(missing))
    carrier = _mode_carrier(spec["carrier"])
    params = {k: float(v) for k, v in (spec.get("parameters") or {}).items()}
    proto = spec["protocol"]
    sweep, bias = proto.get("sweep") or {}, proto.get("bias") or {}
    for key in ("parameter", "from", "to", "rate"):
        if key not in sweep:
            raise SpecError(f"protocol.sweep needs '{key}'")
    if "parameter" not in bias or "value" not in bias:
        raise SpecError("protocol.bias needs 'parameter' and 'value'")
    if sweep["parameter"] in params or bias["parameter"] in params:
        raise SpecError("the protocol variables must not also be declared as parameters")
    if float(sweep["rate"]) <= 0 or float(sweep["to"]) <= float(sweep["from"]):
        raise SpecError("the sweep must increase its parameter at a positive rate")
    names = {k: sp.Symbol(k) for k in params}
    s_eps, s_h = sp.Symbol(sweep["parameter"]), sp.Symbol(bias["parameter"])
    names_all = {**names, sweep["parameter"]: s_eps, bias["parameter"]: s_h}
    subs = {names[k]: v for k, v in params.items()}
    # Hamiltonian terms: coefficient functions of (eps2, h) and operator matrices
    terms = []
    if not isinstance(spec["hamiltonian"], list) or not spec["hamiltonian"]:
        raise SpecError("hamiltonian must be a non-empty list of terms")
    for e in spec["hamiltonian"]:
        if not isinstance(e, dict) or "coefficient" not in e or "operator" not in e:
            raise SpecError("every Hamiltonian term needs 'coefficient' and 'operator'")
        expr = _expression(e["coefficient"], names_all).subs(subs)
        fn = sp.lambdify((s_eps, s_h), expr, "math")
        M = _operator(carrier, e["operator"])
        if e.get("imaginary"):
            M = 1j * M
        if e.get("hc"):
            M = M + M.conj().T
        terms.append({"text": f"{e['coefficient']}{' i' if e.get('imaginary') else ''} [{e['operator']}"
                              f"{' + h.c.' if e.get('hc') else ''}]", "fn": fn, "matrix": M, "operator": e["operator"],
                      "degree": _degree(e["operator"])})
    dim = carrier.dim

    def hamiltonian(eps2: float, h: float) -> np.ndarray:
        H = np.zeros((dim, dim), dtype=complex)
        for t in terms:
            H = H + float(t["fn"](eps2, h)) * t["matrix"]
        return H

    for probe in ((float(sweep["from"]), 0.0), (float(sweep["to"]), float(bias["value"]))):
        H = hamiltonian(*probe)
        if not np.allclose(H, H.conj().T, atol=1e-10 * max(1.0, float(np.abs(H).max()))):
            raise SpecError("the Hamiltonian is not Hermitian; mark non-Hermitian products with \"hc\": true")
    dissipators = []
    if not isinstance(spec["dissipators"], list) or not spec["dissipators"]:
        raise SpecError("dissipators must be a non-empty list")
    for d in spec["dissipators"]:
        if not isinstance(d, dict) or "rate" not in d or "operator" not in d:
            raise SpecError("every dissipator needs 'rate' and 'operator'")
        rate = _real_number(d["rate"], params)
        if rate < 0:
            raise SpecError(f"the rate '{d['rate']}' is negative")
        dissipators.append((rate, _operator(carrier, d["operator"]), str(d["operator"])))
    obs = spec["observable"]
    if not isinstance(obs, dict) or "sign_of" not in obs:
        raise SpecError("observable must be {\"sign_of\": \"a + ad\"}")
    quadrature = _operator(carrier, obs["sign_of"]) / 2.0
    nbar0 = _real_number((proto.get("initial") or {}).get("thermal", 0), params)
    if nbar0 < 0:
        raise SpecError("the initial thermal occupation is negative")
    protocol = {"sweep_parameter": sweep["parameter"], "from": float(sweep["from"]), "to": float(sweep["to"]),
                "rate": float(sweep["rate"]), "hold": float(proto.get("hold", 0.0)),
                "bias_parameter": bias["parameter"], "bias": float(bias["value"]), "nbar0": nbar0}
    return Realization(name=spec["name"], field=str(spec.get("field", "unspecified")), spec=spec, carrier=carrier,
                       params=params, protocol=protocol, hamiltonian=hamiltonian, terms=terms,
                       dissipators=dissipators, quadrature=quadrature, initial=thermal_state(dim, nbar0))


def _degree(operator_text: str) -> int:
    """Largest number of ladder operators in a product of the operator string."""
    return max(len([t for t in p.split() if t in ("a", "ad")]) + 2 * p.split().count("na") for p in operator_text.split("+"))


def thermal_state(dim: int, nbar: float) -> np.ndarray:
    if nbar <= 0:
        rho = np.zeros((dim, dim), dtype=complex)
        rho[0, 0] = 1.0
        return rho
    q = nbar / (nbar + 1.0)
    p = (1.0 - q) * q ** np.arange(dim)
    return np.diag(p / p.sum()).astype(complex)


# ------------------------------------------------------------------------------------------------ evolution
def evolve(real: Realization, rtol: float = 1e-8, atol: float = 1e-11) -> Dict[str, object]:
    """Exact Lindblad evolution over the protocol; returns the final state and the state at the end of the ramp."""
    from scipy.integrate import solve_ivp
    dim = real.carrier.dim
    pr = real.protocol
    T = (pr["to"] - pr["from"]) / pr["rate"]
    jumps = [(g, L, L.conj().T @ L) for g, L, _ in real.dissipators if g > 0]

    def rhs_factory(eps_of_t):
        def rhs(t, y):
            rho = y.reshape(dim, dim)
            H = real.hamiltonian(eps_of_t(t), pr["bias"])
            d = -1j * (H @ rho - rho @ H)
            for g, L, LdL in jumps:
                d += g * (L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL))
            return d.ravel()
        return rhs

    sol = solve_ivp(rhs_factory(lambda t: pr["from"] + pr["rate"] * t), (0.0, T), real.initial.ravel(),
                    method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    rho_ramp = sol.y[:, -1].reshape(dim, dim)
    rho_end, nfev = rho_ramp, int(sol.nfev)
    if pr["hold"] > 0:
        sol2 = solve_ivp(rhs_factory(lambda t: pr["to"]), (0.0, pr["hold"]), rho_ramp.ravel(), method="DOP853",
                         rtol=rtol, atol=atol)
        if not sol2.success:
            raise RuntimeError(sol2.message)
        rho_end, nfev = sol2.y[:, -1].reshape(dim, dim), nfev + int(sol2.nfev)
    return {"rho_end": rho_end, "rho_ramp": rho_ramp, "T_ramp": T, "nfev": nfev}


def hermite_functions(q: np.ndarray, N: int) -> np.ndarray:
    out = np.empty((N, q.size))
    out[0] = math.pi ** -0.25 * np.exp(-q ** 2 / 2.0)
    if N > 1:
        out[1] = math.sqrt(2.0) * q * out[0]
    for k in range(1, N - 1):
        out[k + 1] = math.sqrt(2.0 / (k + 1)) * q * out[k] - math.sqrt(k / (k + 1)) * out[k - 1]
    return out


def prob_positive(rho: np.ndarray) -> float:
    """Born probability that the quadrature (a + a†)/2 is positive, exact for the truncated state."""
    N = rho.shape[0]
    q_max = math.sqrt(2.0 * N + 1.0) + 6.0
    nodes, weights = np.polynomial.legendre.leggauss(6 * N + 200)
    q = 0.5 * q_max * (nodes + 1.0)
    w = 0.5 * q_max * weights
    psi = hermite_functions(q, N)
    return float(np.dot(w, np.einsum("mq,mn,nq->q", psi, rho, psi).real))


def moments(rho: np.ndarray, real: Realization) -> Dict[str, float]:
    x = real.quadrature
    n = real.carrier.tokens["na"]
    xm = float(np.trace(rho @ x).real)
    x2 = float(np.trace(rho @ x @ x).real)
    eig = np.linalg.eigvalsh((rho + rho.conj().T) / 2.0)
    return {"trace": float(np.trace(rho).real), "x_mean": xm, "x_std": math.sqrt(max(x2 - xm ** 2, 0.0)),
            "n_mean": float(np.trace(rho @ n).real), "min_eigenvalue": float(eig.min()),
            "top_population": float(np.diag(rho).real[-5:].sum())}


# ------------------------------------------------------------------------------------------------ the letters
def linear_stage(real: Realization, eps2: float, h: float) -> Dict[str, object]:
    """Mean-field drift of <a> from the quadratic part of H at fixed drive: d alpha/dt = -i w alpha - 2 i g alpha* + c."""
    dim = real.carrier.dim
    H = np.zeros((dim, dim), dtype=complex)
    for t in real.terms:
        if t["degree"] <= 2:
            H = H + float(t["fn"](eps2, h)) * t["matrix"]
    w = float((H[1, 1] - H[0, 0]).real)
    g = complex(H[2, 0]) / math.sqrt(2.0)
    c = -1j * complex(H[1, 0])            # push on alpha from the linear terms c a† + c* a
    kappa = sum(g_ for g_, _, op in real.dissipators if op.strip() == "a") - \
        sum(g_ for g_, _, op in real.dissipators if op.strip() == "ad")
    nbar = sum(g_ for g_, _, op in real.dissipators if op.strip() == "ad") / kappa if kappa > 0 else 0.0
    # linear map on (x, y): alpha' = -(kappa/2 + i w) alpha - 2 i g alpha*
    A = np.array([[-kappa / 2 + 2 * g.imag, w - 2 * g.real], [-w - 2 * g.real, -kappa / 2 - 2 * g.imag]])
    vals, vecs = np.linalg.eig(A)
    i = int(np.argmax(vals.real))
    v = vecs[:, i].real if abs(vecs[:, i].imag).max() < 1e-9 else None
    return {"w": w, "g": g, "push": c, "kappa": kappa, "nbar": nbar, "growth": float(vals[i].real),
            "direction": v, "gain": 2.0 * abs(g)}


def derive_quantum_write(real: Realization, law: bool = True) -> Dict[str, object]:
    pr = real.protocol
    out = {"name": real.name, "field": real.field, "target": "quantum-write", "word": "", "steps": [], "status": ""}
    steps, word = out["steps"], []

    # G: the generator in the linear stage, the threshold crossing, the amplified quadrature, the bias
    lo, hi = linear_stage(real, pr["from"], 0.0), linear_stage(real, pr["to"], 0.0)
    kappa, nbar = hi["kappa"], hi["nbar"]
    if kappa <= 0:
        return _obstructed(out, "G", "no loss of the mode: the open carrier needs a dissipator on a", "no loss")
    if lo["growth"] >= 0:
        return _obstructed(out, "G", f"the mode is already unstable at the start of the sweep (growth rate {lo['growth']:.3g})",
                           "no threshold")
    if hi["growth"] <= 0:
        return _obstructed(out, "G", f"the sweep does not reach the threshold (growth rate {hi['growth']:.3g} at the end)",
                           "no threshold")
    # the threshold drive by bisection on the growth rate
    a, b = pr["from"], pr["to"]
    for _ in range(80):
        m = 0.5 * (a + b)
        if linear_stage(real, m, 0.0)["growth"] > 0:
            b = m
        else:
            a = m
    threshold = 0.5 * (a + b)
    direction = hi["direction"]
    if direction is None or abs(direction[1]) > TOL_DIRECTION * max(1.0, abs(direction[0])):
        return _obstructed(out, "G", "the amplified quadrature is not the measured one (a detuning or a pump phase "
                                     "rotates it)", "wrong quadrature")
    push = linear_stage(real, pr["to"], 1.0)["push"]
    if abs(push.imag) > TOL_DIRECTION * max(1.0, abs(push.real)) or push.real == 0:
        return _obstructed(out, "G", "the bias does not push the amplified quadrature: the selection would come "
                                     "only through the Kerr rotation", "bias off the quadrature")
    nonlinear = [t["text"] for t in real.terms if t["degree"] > 2]
    steps.append({"letter": "G", "text": f"quadratic generator with loss kappa = {kappa:.4g} (nbar = {nbar:.3g}); "
                                         f"threshold at {pr['sweep_parameter']} = {threshold:.6g}; the amplified "
                                         f"quadrature is the measured one; the bias pushes it with {push.real:.4g} per "
                                         f"unit of {pr['bias_parameter']}; nonlinear terms: "
                                         + (", ".join(nonlinear) or "none")})
    word.append("G")

    # K: the canonical linear equation
    def gain(eps2):
        return linear_stage(real, eps2, 0.0)["gain"]
    r_eff = (gain(pr["to"]) - gain(pr["from"])) / (pr["to"] - pr["from"]) * pr["rate"]
    h_eff = push.real * pr["bias"]
    two_d = kappa * (2.0 * nbar + 1.0) / 4.0
    sigma0_sq = (2.0 * pr["nbar0"] + 1.0) / 4.0
    steps.append({"letter": "K", "text": f"x' = (eps(t) - kappa/2) x + h + sqrt(2D) xi with eps(t) = {gain(pr['from']):.4g} + "
                                         f"{r_eff:.4g} t, h = {h_eff:.5g}, 2D = {two_d:.5g}, sigma0^2 = {sigma0_sq:.4g}",
                  "r": r_eff, "h": h_eff, "two_d": two_d, "sigma0_sq": sigma0_sq, "kappa": kappa, "nbar": nbar})
    word.append("K")
    canon = {"r": r_eff, "h": h_eff, "two_d": two_d, "sigma0_sq": sigma0_sq, "kappa": kappa, "nbar": nbar,
             "threshold": threshold, "gain_from": gain(pr["from"]), "gain_to": gain(pr["to"])}
    out["canonical"] = canon
    p_law = protocol_probability(canon, pr)
    out["law"] = {"P_law": p_law}
    if not law:
        out["word"], out["status"], out["class"] = "".join(word), "derived without the law", "linear-stage write"
        return out

    # L: exact evolution against the law
    ev = evolve(real)
    m_end = moments(ev["rho_end"], real)
    p_exact = prob_positive(ev["rho_end"])
    trunc = max(m_end["top_population"] * 10.0, abs(m_end["trace"] - 1.0) * 10.0, 1e-5)
    out["law"].update({"P_exact": p_exact, "difference": p_exact - p_law, "truncation": trunc,
                       "n_end": m_end["n_mean"], "x_end": m_end["x_mean"], "top_population": m_end["top_population"],
                       "min_eigenvalue": m_end["min_eigenvalue"], "T_ramp": ev["T_ramp"], "nfev": ev["nfev"]})
    out["law_constant"] = {"constant": p_exact, "expected": p_law, "stderr": trunc}
    if abs(p_exact - p_law) <= max(TOL_LAW, 4.0 * trunc):
        steps.append({"letter": "L", "text": f"P_exact = {p_exact:.5f} against the law {p_law:.5f} "
                                             f"(truncation {trunc:.1e}): the linear-stage write", "law_constant": p_exact})
        word.append("L")
        out["word"], out["status"], out["class"] = "".join(word), REACHED, "linear-stage write"
        return out
    # the law fails: is it the equilibrium write?
    eq = equilibrium(real, pr["to"], pr["bias"])
    out["equilibrium"] = eq
    if eq is not None:
        near = equilibrium(real, min(threshold + 0.5 * kappa, pr["to"]), pr["bias"])
        eq["gap_near_threshold"] = near["gap"] if near else None
    if eq is not None and abs(p_exact - eq["P_eq"]) <= 0.02:
        return _obstructed(out, "L", f"P_exact = {p_exact:.4f} leaves the law ({p_law:.4f}) and equals the selection "
                                     f"of the biased steady state P_eq = {eq['P_eq']:.4f}; the switching gap just "
                                     f"above threshold, {eq['gap_near_threshold']:.3g}, is not small against the sweep "
                                     f"scale sqrt(r) = {math.sqrt(r_eff):.3g}, so the wells equilibrate during the "
                                     f"passage: the write is decided by the balance of the wells", "equilibrium write")
    return _obstructed(out, "L", f"P_exact = {p_exact:.4f} leaves the law ({p_law:.4f})"
                                 + (f" and the biased steady state ({eq['P_eq']:.4f})" if eq else "")
                                 + ": neither the linear stage nor the balance of the wells decides alone",
                       "intermediate regime")


def _obstructed(out: Dict, letter: str, reason: str, short: str) -> Dict:
    out["steps"].append({"letter": letter, "text": reason, "obstruction": True})
    out["word"] = "".join(s["letter"] for s in out["steps"] if not s.get("obstruction"))
    out["status"] = f"obstructed at {letter}"
    out["obstruction"], out["obstruction_short"] = reason, short
    out["class"] = short
    return out


def protocol_probability(canon: Dict, pr: Dict) -> float:
    """P = Phi(h I1 / sqrt(sigma0^2 + 2 D I2)) with the amplification exponent integrated over the ramp and the hold."""
    from scipy.integrate import quad
    from scipy.special import ndtr
    T = (pr["to"] - pr["from"]) / pr["rate"]
    k2 = canon["kappa"] / 2.0

    def eps(t):
        return canon["gain_from"] + canon["r"] * min(t, T)

    def phi(t):
        if t <= T:
            return (canon["gain_from"] - k2) * t + 0.5 * canon["r"] * t ** 2
        return (canon["gain_from"] - k2) * T + 0.5 * canon["r"] * T ** 2 + (eps(T) - k2) * (t - T)

    t_end = T + pr["hold"]
    t_star = min(max((k2 - canon["gain_from"]) / canon["r"], 0.0), t_end)
    phi_min = phi(t_star)

    def seg(k, a, b):
        if b <= a:
            return 0.0
        pts = [t_star] if a < t_star < b else None
        return quad(lambda s: math.exp(-k * (phi(s) - phi_min)), a, b, points=pts, limit=400)[0]

    log_i1 = math.log(seg(1, 0.0, min(T, t_end)) + seg(1, T, t_end)) - phi_min
    log_i2 = math.log(seg(2, 0.0, min(T, t_end)) + seg(2, T, t_end)) - 2.0 * phi_min
    log_var = float(np.logaddexp(math.log(canon["sigma0_sq"]), math.log(canon["two_d"]) + log_i2)) \
        if canon["two_d"] > 0 else math.log(canon["sigma0_sq"])
    if canon["h"] == 0:
        return 0.5
    z = math.copysign(math.exp(math.log(abs(canon["h"])) + log_i1 - 0.5 * log_var), canon["h"])
    return float(ndtr(z))


def equilibrium(real: Realization, eps2: float, h: float) -> Optional[Dict[str, float]]:
    """The biased steady state at a fixed drive and the smallest nonzero decay rate (dense; small carriers only)."""
    dim = real.carrier.dim
    if dim > STEADY_MAX_DIM:
        return None
    I = np.eye(dim, dtype=complex)
    H = real.hamiltonian(eps2, h)
    L = -1j * (np.kron(H, I) - np.kron(I, H.T))
    for g, J, _ in real.dissipators:
        if g > 0:
            JdJ = J.conj().T @ J
            L += g * (np.kron(J, J.conj()) - 0.5 * (np.kron(JdJ, I) + np.kron(I, JdJ.T)))
    vals, vecs = np.linalg.eig(L)
    order = np.argsort(-vals.real)
    rho = vecs[:, order[0]].reshape(dim, dim)
    rho = (rho + rho.conj().T) / 2.0
    rho /= np.trace(rho).real
    return {"P_eq": prob_positive(rho), "gap": float(-vals[order[1]].real), "n_eq": float(np.trace(rho @ real.carrier.tokens["na"]).real)}


# ------------------------------------------------------------------------------------------------ reports
def codiscover(reals: List[Realization], law: bool = True) -> Dict[str, object]:
    rows = [derive_quantum_write(r, law=law) for r in reals]
    reached = [r for r in rows if r["status"] == REACHED]
    summary = {"models": len(rows), "reached": len(reached), "obstructed": len(rows) - len(reached),
               "fields_reached": sorted({r["field"] for r in reached}),
               "classes": {c: sum(1 for r in rows if r.get("class") == c) for c in sorted({r.get("class", "") for r in rows})},
               "law_difference_max": max((abs(r["law"].get("difference", 0.0)) for r in reached), default=0.0)}
    return {"rows": rows, "summary": summary}


def markdown(report: Dict[str, object]) -> str:
    s = report["summary"]
    lines = ["# The quantum write on open carriers", "",
             f"{s['models']} realizations; {s['reached']} reach the linear-stage write; classes: "
             + ", ".join(f"{k}: {v}" for k, v in s["classes"].items()), "",
             "| realization | field | word | class | P exact | P law | status |", "|---|---|---|---|---:|---:|---|"]
    for r in report["rows"]:
        law = r.get("law", {})
        lines.append(f"| {r['name']} | {r['field']} | {r['word']} | {r.get('class', '')} | "
                     f"{law.get('P_exact', float('nan')):.4f} | {law.get('P_law', float('nan')):.4f} | {r['status']} |")
    lines.append("")
    for r in report["rows"]:
        lines.append(f"## {r['name']}")
        for st in r["steps"]:
            lines.append(f"- **{st['letter']}** {st['text']}")
        lines.append("")
    return "\n".join(lines)
