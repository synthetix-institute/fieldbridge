"""Open quantum networks: several modes and multilevel systems, time-dependent controls, losses, and pulses that
enter and leave through transmission lines. The model for which it was built is a quantum memory: a photon arrives
in a readout resonator, a transmon passes it to a storage resonator, and later the transmon passes it back and the
readout resonator emits it into the line.

Specification schema ``fieldbridge-quantum-network/1``:

    subsystems   [{"name": "r", "kind": "mode", "max_quanta": N},          tokens r, rd (creation), nr (number)
                  {"name": "q", "kind": "levels", "levels": d}]            tokens q (lowering, sqrt(k) |k-1><k|), qd,
                                                                           nq, and q_j_k = |j><k| for 0 <= j, k < d
    parameters   {name: value}                                             real numbers
    controls     {name: [{"from": t0, "to": t1, "value": expr}]}           piecewise in time, zero outside the
                                                                           segments; expr in t and the parameters
    hamiltonian  [{"coefficient": expr, "operator": "rd q", "hc": bool, "imaginary": bool}]
                                                                           expr in t, the parameters and controls
    dissipators  [{"rate": expr, "operator": "q", "name": "transmon decay"}]   rate * D[operator], constant rates
    ports        [{"name": "line", "operator": "r", "rate": expr,          a transmission line coupled to the
                   "input": {...}, "output": {...}}]                       operator at the given rate
                 input  {"kind": "fock", "photons": n, "shape": expr, "from": t0, "to": t1}
                        {"kind": "coherent", "amplitude": expr, "from": t0, "to": t1}  (photon-flux amplitude)
                 output {"shape": expr, "from": t0, "to": t1}             the temporal mode that is measured
    initial      {"r": 0, "q": 1, "s": {"thermal": expr}}                  number states, levels or thermal modes;
                                                                           subsystems not named start in 0
    duration     T                                                         the evolution runs from 0 to T
    observables  {"name": "operator"}                                      expectation values, recorded in time
    windows      [{"port": "line", "from": t0, "to": t1}]                  photons emitted into a line in a window

Input-output follows Gardiner and Collett: a port with operator c and rate kappa contributes the line operator
L = sqrt(kappa) c, the field leaving the system is b_out = b_in + sqrt(kappa) c, and a port without input or output
is the loss D[L]. Pulses are attached by the series product of the SLH formalism (Gough and James 2009). A coherent
input is the displaced vacuum b_in = beta(t). A Fock input of n photons in the temporal mode u(t) is emitted by a
virtual source mode whose coupling g_u(t) = u*(t) / sqrt(1 - int_0^t |u|^2) releases exactly that mode; the measured
output mode v(t) is absorbed by a virtual detector mode with g_v(t) = -v*(t) / sqrt(int_0^t |v|^2), whose final
occupation is the number of photons the line carries in v (Kiilerich and Moelmer, Phys. Rev. Lett. 123, 123604
(2019)). The shapes are normalized on their windows.

The measurement reports, for every port, the photons that arrived, the photons found in the output mode, the
photons emitted into the line in each window, the excitations each loss removes and the excitations left in the
system: a budget that closes when the Hamiltonian conserves the number of excitations.

Nothing in a specification is evaluated: expressions go through the restricted parser of the memory module.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple, Union

import numpy as np

from ..memory.spec import SpecError as MemorySpecError, parse_expression
from .carriers import Carrier, CarrierError, kron_all

SCHEMA = "fieldbridge-quantum-network/1"
REQUIRED = ("name", "question", "assumptions", "provenance", "subsystems", "duration")
MAX_DIM = 400                 # density operator of at most 400 x 400
MAX_SEGMENTS = 400
SHAPE_POINTS = 20001          # grid for the normalization and the running integral of a pulse shape
NAME = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")
RESERVED = {"t", "pi", "I", "E", "S", "N", "O", "Q"}


class SpecError(ValueError):
    pass


@dataclass
class Shape:
    """A normalized temporal mode on a window, with its running integral."""
    t: np.ndarray
    amp: np.ndarray
    cum: np.ndarray
    start: float
    stop: float

    def at(self, t: float, tm: Optional[float] = None) -> complex:
        """The amplitude at t; tm, the middle of the interval being integrated, decides whether the window is on."""
        c = t if tm is None else tm
        if c < self.start or c > self.stop:
            return 0.0
        return complex(np.interp(t, self.t, self.amp.real) + 1j * np.interp(t, self.t, self.amp.imag))

    def integral(self, t: float) -> float:
        if t <= self.start:
            return 0.0
        if t >= self.stop:
            return 1.0
        return float(np.interp(t, self.t, self.cum))


@dataclass
class Port:
    name: str
    operator: str
    rate: float
    L: np.ndarray                         # sqrt(rate) * operator on the full space
    source: Optional[Dict] = None         # {"kind", "photons" or "shape", "mode" (index of the virtual mode)}
    detector: Optional[Dict] = None


@dataclass
class Network:
    name: str
    field: str
    spec: Dict
    carrier: Carrier
    subsystems: List[Dict]
    params: Dict[str, float]
    controls: Dict[str, List[Tuple[float, float, Callable[[float], float]]]]
    terms: List[Dict]
    dissipators: List[Dict]
    ports: List[Port]
    initial: np.ndarray
    duration: float
    breakpoints: List[float]
    observables: Dict[str, np.ndarray]
    windows: List[Dict]
    excitations: np.ndarray
    notes: List[str] = field(default_factory=list)


# ------------------------------------------------------------------------------------------------ expressions
def _expression(text, names: Dict[str, object]):
    try:
        return parse_expression(str(text), names)
    except MemorySpecError as err:
        raise SpecError(str(err)) from err


def _function(text, params: Dict[str, float], variables: Tuple[str, ...]) -> Callable:
    """A numerical function of the variables (time, controls) with the parameters substituted."""
    import sympy as sp
    names = {k: sp.Symbol(k) for k in params}
    names.update({v: sp.Symbol(v) for v in variables})
    expr = _expression(text, names).subs({names[k]: v for k, v in params.items()})
    return sp.lambdify([names[v] for v in variables], expr, "math")


def _real_number(text, params: Dict[str, float]) -> float:
    val = complex(_function(text, params, ())())
    if abs(val.imag) > 1e-12:
        raise SpecError(f"'{text}' is not real")
    return float(val.real)


def _shape(spec: Dict, params: Dict[str, float], what: str, key: str) -> Shape:
    for k in (key, "from", "to"):
        if k not in spec:
            raise SpecError(f"{what} needs '{k}'")
    t0, t1 = float(spec["from"]), float(spec["to"])
    if not t1 > t0 >= 0:
        raise SpecError(f"{what}: the window must satisfy 0 <= from < to")
    fn = _function(spec[key], params, ("t",))
    t = np.linspace(t0, t1, SHAPE_POINTS)
    amp = np.array([complex(fn(x)) for x in t])
    if not np.all(np.isfinite(amp)):
        raise SpecError(f"{what}: the shape is not finite on its window")
    w = np.abs(amp) ** 2
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(t))])
    if cum[-1] <= 0:
        raise SpecError(f"{what}: the shape vanishes on its window")
    return Shape(t, amp / math.sqrt(cum[-1]), cum / cum[-1], t0, t1)


# ------------------------------------------------------------------------------------------------ carrier
def _subsystem_tokens(sub: Dict) -> Tuple[int, Dict[str, np.ndarray]]:
    name, kind = sub["name"], sub["kind"]
    if kind == "mode":
        n = int(sub.get("max_quanta", 0))
        if not 1 <= n <= 160:
            raise SpecError(f"mode '{name}': max_quanta must be between 1 and 160")
        a = np.diag(np.sqrt(np.arange(1, n + 1)), 1).astype(complex)
        return n + 1, {name: a, name + "d": a.conj().T, "n" + name: a.conj().T @ a}
    if kind == "levels":
        d = int(sub.get("levels", 0))
        if not 2 <= d <= 12:
            raise SpecError(f"levels '{name}': levels must be between 2 and 12")
        a = np.diag(np.sqrt(np.arange(1, d)), 1).astype(complex)
        tok = {name: a, name + "d": a.conj().T, "n" + name: a.conj().T @ a}
        for j in range(d):
            for k in range(d):
                P = np.zeros((d, d), dtype=complex)
                P[j, k] = 1.0
                tok[f"{name}_{j}_{k}"] = P
        return d, tok
    raise SpecError(f"subsystem '{name}': kind must be 'mode' or 'levels'")


def _carrier(subs: List[Dict]) -> Tuple[Carrier, List[int]]:
    dims, local = [], []
    for sub in subs:
        d, tok = _subsystem_tokens(sub)
        dims.append(d)
        local.append(tok)
    dim = int(np.prod(dims))
    if dim > MAX_DIM:
        raise SpecError(f"the product space has dimension {dim}; at most {MAX_DIM}")
    tokens = {"I": np.eye(dim, dtype=complex)}
    for k, tok in enumerate(local):
        for name, M in tok.items():
            if name in tokens:
                raise SpecError(f"operator name '{name}' is defined twice; rename a subsystem")
            tokens[name] = kron_all([M if i == k else np.eye(dims[i], dtype=complex) for i in range(len(dims))])
    desc = ", ".join(f"{s['name']} ({s['kind']}, {d} states)" for s, d in zip(subs, dims))
    return Carrier("network", dim, tokens, desc, {"subsystems": subs}), dims


def _operator(carrier: Carrier, text: str) -> np.ndarray:
    try:
        return carrier.operator(text)
    except CarrierError as err:
        raise SpecError(str(err)) from err


# ------------------------------------------------------------------------------------------------ loading
def load(source: Union[str, Path, Dict]) -> Network:
    spec = source if isinstance(source, dict) else json.loads(Path(source).read_text(encoding="utf-8"))
    if spec.get("schema") != SCHEMA:
        raise SpecError(f"schema must be {SCHEMA}")
    missing = [k for k in REQUIRED if not spec.get(k)] + (["hamiltonian"] if "hamiltonian" not in spec else [])
    if missing:
        raise SpecError("missing: " + ", ".join(missing))
    params = {k: float(v) for k, v in (spec.get("parameters") or {}).items()}
    controls_spec = spec.get("controls") or {}
    for nm in list(params) + list(controls_spec):
        if not NAME.match(nm) or nm in RESERVED:
            raise SpecError(f"'{nm}' cannot name a parameter or a control (reserved: {', '.join(sorted(RESERVED))})")
    clash = set(params) & set(controls_spec)
    if clash:
        raise SpecError("names used both as parameter and control: " + ", ".join(sorted(clash)))
    duration = float(spec["duration"])
    if duration <= 0:
        raise SpecError("duration must be positive")

    # subsystems, with a virtual mode for every pulse a port sends or receives
    subs = [dict(s) for s in spec["subsystems"]]
    names = [s.get("name") for s in subs]
    for nm in names:
        if not isinstance(nm, str) or not NAME.match(nm):
            raise SpecError(f"subsystem name '{nm}' must be letters and digits, starting with a letter")
    if len(set(names)) != len(names):
        raise SpecError("subsystem names must be distinct")
    ports_spec = spec.get("ports") or []
    virtual = []
    for p in ports_spec:
        inp, outp = p.get("input"), p.get("output")
        if inp and inp.get("kind") == "fock":
            n = int(inp.get("photons", 1))
            if not 1 <= n <= 4:
                raise SpecError(f"port '{p.get('name')}': a Fock input carries 1 to 4 photons")
            virtual.append({"name": f"{p['name']}In", "kind": "mode", "max_quanta": n, "virtual": True})
        if outp:
            n = int(outp.get("max_photons", (inp or {}).get("photons", 1) if (inp or {}).get("kind") == "fock" else 1))
            virtual.append({"name": f"{p['name']}Out", "kind": "mode", "max_quanta": n, "virtual": True})
    carrier, dims = _carrier(subs + virtual)
    dim = carrier.dim
    all_subs = subs + virtual

    # controls: piecewise functions of time
    controls, breaks = {}, {0.0, duration}
    for nm, segs in controls_spec.items():
        if not isinstance(segs, list) or not segs:
            raise SpecError(f"control '{nm}' must be a list of segments")
        rows = []
        for s in segs:
            for k in ("from", "to", "value"):
                if k not in s:
                    raise SpecError(f"control '{nm}': every segment needs '{k}'")
            a, b = float(s["from"]), float(s["to"])
            if not b > a:
                raise SpecError(f"control '{nm}': a segment must end after it starts")
            rows.append((a, b, _function(s["value"], params, ("t",))))
            breaks.update((a, b))
        rows.sort(key=lambda r: r[0])
        if any(r2[0] < r1[1] for r1, r2 in zip(rows, rows[1:])):
            raise SpecError(f"control '{nm}': segments overlap")
        controls[nm] = rows
    cnames = tuple(controls)

    def control_values(t: float, tm: Optional[float] = None) -> List[float]:
        """Control values at t, from the segment that contains tm (the middle of the interval being integrated),
        so that an integration step that ends on a switch does not see the next segment."""
        c = t if tm is None else tm
        out = []
        for nm in cnames:
            v = 0.0
            for a, b, fn in controls[nm]:
                if a <= c < b or (c == b == duration):
                    v = float(fn(t))
                    break
            out.append(v)
        return out

    # Hamiltonian terms
    if not isinstance(spec["hamiltonian"], list):
        raise SpecError("hamiltonian must be a list of terms")
    terms = []
    for e in spec["hamiltonian"]:
        if not isinstance(e, dict) or "coefficient" not in e or "operator" not in e:
            raise SpecError("every Hamiltonian term needs 'coefficient' and 'operator'")
        fn = _function(e["coefficient"], params, ("t",) + cnames)
        M = _operator(carrier, e["operator"])
        if e.get("imaginary"):
            M = 1j * M
        if e.get("hc"):
            M = M + M.conj().T
        terms.append({"text": f"{e['coefficient']}{' i' if e.get('imaginary') else ''} [{e['operator']}"
                              f"{' + h.c.' if e.get('hc') else ''}]", "fn": fn, "matrix": M})
    for tp in sorted(breaks):
        cv = control_values(min(tp, duration))
        H = sum((complex(t["fn"](tp, *cv)) * t["matrix"] for t in terms), np.zeros((dim, dim), dtype=complex))
        if not np.allclose(H, H.conj().T, atol=1e-10 * max(1.0, float(np.abs(H).max(initial=0.0)))):
            raise SpecError("the Hamiltonian is not Hermitian; mark non-Hermitian products with \"hc\": true")

    # dissipators
    dissipators = []
    for k, d in enumerate(spec.get("dissipators") or []):
        if not isinstance(d, dict) or "rate" not in d or "operator" not in d:
            raise SpecError("every dissipator needs 'rate' and 'operator'")
        rate = _real_number(d["rate"], params)
        if rate < 0:
            raise SpecError(f"the rate '{d['rate']}' is negative")
        L = math.sqrt(rate) * _operator(carrier, d["operator"])
        dissipators.append({"name": str(d.get("name", f"{d['operator']} at rate {d['rate']}")), "L": L})
    # the excitation number of the physical subsystems, and the rate at which each loss changes it: D^dagger[L](N)
    exc = sum((carrier.tokens["n" + s["name"]] for s in subs), np.zeros((dim, dim), dtype=complex))
    for d in dissipators:
        L, Ld = d["L"], d["L"].conj().T
        d["dN"] = Ld @ exc @ L - 0.5 * (Ld @ L @ exc + exc @ Ld @ L)

    # ports, with their pulses
    ports = []
    for p in ports_spec:
        for k in ("name", "operator", "rate"):
            if k not in p:
                raise SpecError(f"every port needs '{k}'")
        rate = _real_number(p["rate"], params)
        if rate < 0:
            raise SpecError(f"port '{p['name']}': the rate is negative")
        port = Port(p["name"], p["operator"], rate, math.sqrt(rate) * _operator(carrier, p["operator"]))
        inp = p.get("input")
        if inp:
            kind = inp.get("kind")
            if kind == "fock":
                shape = _shape(inp, params, f"port '{p['name']}' input", "shape")
                port.source = {"kind": "fock", "photons": int(inp.get("photons", 1)), "shape": shape,
                               "a": carrier.tokens[f"{p['name']}In"]}
            elif kind == "coherent":
                for k in ("amplitude", "from", "to"):
                    if k not in inp:
                        raise SpecError(f"port '{p['name']}' coherent input needs '{k}'")
                a0, a1 = float(inp["from"]), float(inp["to"])
                fn = _function(inp["amplitude"], params, ("t",))
                tt = np.linspace(a0, a1, SHAPE_POINTS)
                flux = np.abs(np.array([complex(fn(x)) for x in tt])) ** 2
                port.source = {"kind": "coherent", "fn": fn, "from": a0, "to": a1,
                               "photons": float(np.sum(0.5 * (flux[1:] + flux[:-1]) * np.diff(tt)))}
            else:
                raise SpecError(f"port '{p['name']}': input kind must be 'fock' or 'coherent'")
            breaks.update((float(inp["from"]), float(inp["to"])))
        outp = p.get("output")
        if outp:
            shape = _shape(outp, params, f"port '{p['name']}' output", "shape")
            port.detector = {"shape": shape, "a": carrier.tokens[f"{p['name']}Out"]}
            breaks.update((shape.start, shape.stop))
        ports.append(port)

    # initial product state
    init = spec.get("initial") or {}
    unknown = set(init) - set(names)
    if unknown:
        raise SpecError("initial names unknown subsystems: " + ", ".join(sorted(unknown)))
    factors = []
    for s, d in zip(all_subs, dims):
        v = init.get(s["name"], 0)
        rho = np.zeros((d, d), dtype=complex)
        if isinstance(v, dict) and "thermal" in v and s["kind"] == "mode":
            nbar = _real_number(v["thermal"], params)
            if nbar < 0:
                raise SpecError(f"initial '{s['name']}': negative thermal occupation")
            q = nbar / (nbar + 1.0)
            p_ = q ** np.arange(d)
            rho = np.diag(p_ / p_.sum()).astype(complex)
        elif isinstance(v, int) and 0 <= v < d:
            rho[v, v] = 1.0
        else:
            raise SpecError(f"initial '{s['name']}' must be a state index below {d}" +
                            (" or {\"thermal\": nbar}" if s["kind"] == "mode" else ""))
        if s.get("virtual") and s["name"].endswith("In"):
            port = next(pt for pt in ports if f"{pt.name}In" == s["name"])
            rho = np.zeros((d, d), dtype=complex)
            rho[port.source["photons"], port.source["photons"]] = 1.0
        factors.append(rho)
    initial = kron_all(factors)

    observables = {k: _operator(carrier, v) for k, v in (spec.get("observables") or {}).items()}
    windows = []
    pnames = [p.name for p in ports]
    for w in spec.get("windows") or []:
        if w.get("port") not in pnames or "from" not in w or "to" not in w:
            raise SpecError("every window needs a known 'port', 'from' and 'to'")
        if not 0.0 <= float(w["from"]) < float(w["to"]) <= duration:
            raise SpecError("a window must lie inside the evolution and end after it starts")
        windows.append({"port": w["port"], "from": float(w["from"]), "to": float(w["to"])})
        breaks.update((float(w["from"]), float(w["to"])))
    bp = sorted(b for b in breaks if 0.0 <= b <= duration)
    if len(bp) > MAX_SEGMENTS:
        raise SpecError(f"more than {MAX_SEGMENTS} time segments")
    net = Network(name=spec["name"], field=str(spec.get("field", "unspecified")), spec=spec, carrier=carrier,
                  subsystems=all_subs, params=params, controls=controls, terms=terms, dissipators=dissipators,
                  ports=ports, initial=initial, duration=duration, breakpoints=bp, observables=observables,
                  windows=windows, excitations=exc)
    net.control_values = control_values  # type: ignore[attr-defined]
    return net


# ------------------------------------------------------------------------------------------------ dynamics
def hamiltonian(net: Network, t: float, tm: Optional[float] = None) -> np.ndarray:
    cv = net.control_values(t, tm)  # type: ignore[attr-defined]
    dim = net.carrier.dim
    H = np.zeros((dim, dim), dtype=complex)
    for term in net.terms:
        H = H + complex(term["fn"](t, *cv)) * term["matrix"]
    return H


def port_components(port: Port, eye: np.ndarray) -> List[Tuple[Callable[[float, Optional[float]], complex], np.ndarray]]:
    """The components of a port in cascade order, source -> system -> detector, each as (a(t, tm), A): the line
    operator of the cascade is L = sum_i a_i(t) A_i. The couplings of the virtual modes diverge where the pulse is
    nearly all emitted (source) or has only begun (detector); they are set to zero where less than 1e-9 of the pulse
    remains or has passed."""
    comps = []
    src = port.source
    if src is not None:
        if src["kind"] == "coherent":
            def beta(t, tm, src=src):
                c = t if tm is None else tm
                return complex(src["fn"](t)) if src["from"] <= c <= src["to"] else 0.0
            comps.append((beta, eye))
        else:
            def g_in(t, tm, sh=src["shape"]):
                left = 1.0 - sh.integral(t)
                return sh.at(t, tm) / math.sqrt(left) if left > 1e-9 else 0.0     # conj(g_u) with g_u = u*/sqrt(left)
            comps.append((g_in, src["a"]))
    comps.append((lambda t, tm, k=math.sqrt(port.rate): k, port.L / math.sqrt(port.rate) if port.rate > 0 else port.L))
    if port.detector is not None:
        def g_out(t, tm, sh=port.detector["shape"]):
            seen = sh.integral(t)
            return -sh.at(t, tm) / math.sqrt(seen) if seen > 1e-9 else 0.0       # conj(g_v) with g_v = -v*/sqrt(seen)
        comps.append((g_out, port.detector["a"]))
    return comps


def evolve(net: Network, samples: int = 0, rtol: float = 1e-8, atol: float = 1e-10) -> Dict[str, object]:
    """Lindblad evolution from 0 to the duration, integrated between breakpoints so that no step straddles a switch.
    Alongside the density operator it integrates the photons each port emits into its line and the excitations each
    loss removes.

    The right-hand side is -i (H_eff rho - rho H_eff^dagger) + sum_k L_k rho L_k^dagger with
    H_eff = H - (i/2) sum_k L_k^dagger L_k; both products with rho are computed, so rho is never assumed Hermitian
    during the integration. A port with pulses is the SLH series product of its components (Gough and James 2009),
    L = sum_i a_i A_i in cascade order; the Hamiltonian of the series product cancels every backward term, and
    H_eff of the cascade is -(i/2) sum_i |a_i|^2 A_i^dagger A_i - i sum_{i > j} a_i* a_j A_i^dagger A_j.
    Operators are sparse in the product basis and built once."""
    from scipy import sparse
    from scipy.integrate import solve_ivp
    dim = net.carrier.dim
    csr = lambda M: sparse.csr_matrix(M)
    nd, npt = len(net.dissipators), len(net.ports)

    # every operator acts on rho from the left (M @ rho) or from the right (rho @ M = (M^T @ rho^T)^T)
    def left_right(M):
        return csr(M), csr(np.asarray(M).T)

    hterms = [(t["fn"], *left_right(t["matrix"])) for t in net.terms]
    damping = sum((d["L"].conj().T @ d["L"] for d in net.dissipators), np.zeros((dim, dim), dtype=complex))
    damp_l, damp_r = left_right(damping)
    jumps = [(csr(d["L"]), csr(d["L"].conj())) for d in net.dissipators]          # L @ rho and (conj(L) @ X^T)^T
    dN = [csr(d["dN"]).T.tocsr() for d in net.dissipators]                        # Tr(E rho) = sum(E^T * rho)
    eye = np.eye(dim, dtype=complex)
    cascades = []
    for port in net.ports:
        comps = port_components(port, eye)
        ops = [(a, csr(A), csr(A.conj())) for a, A in comps]
        pairs = {}
        for i in range(len(comps)):
            for j in range(i + 1):
                P = comps[i][1].conj().T @ comps[j][1]
                pairs[(i, j)] = (csr(P), csr(P.conj()))                       # P @ rho and rho @ P^dagger
        cascades.append((ops, pairs))

    def rhs(t, y, tm):
        rho = y[:dim * dim].reshape(dim, dim)
        rhoT = rho.T
        cv = net.control_values(t, tm)  # type: ignore[attr-defined]
        d = np.zeros((dim, dim), dtype=complex)
        for fn, Ml, Mr in hterms:                                         # -i [H, rho]
            c = complex(fn(t, *cv))
            if c != 0:
                d -= 1j * c * (Ml @ rho - (Mr @ rhoT).T)
        d -= 0.5 * (damp_l @ rho + (damp_r @ rhoT).T)                     # static losses, anticommutator
        rates = []
        for (L, Lc), E in zip(jumps, dN):
            X = L @ rho
            d += (Lc @ X.T).T
            rates.append(-np.sum(E.multiply(rho)))
        for ops, pairs in cascades:
            a = [f(t, tm) for f, _, _ in ops]
            for (i, j), (P, Pc) in pairs.items():                         # H_eff of the cascade
                c = -0.5j * abs(a[i]) ** 2 if i == j else -1j * np.conj(a[i]) * a[j]
                if c != 0:
                    d -= 1j * (c * (P @ rho) - np.conj(c) * (Pc @ rhoT).T)
            X = sum((ai * (A @ rho) for ai, (_, A, _) in zip(a, ops) if ai != 0), np.zeros((dim, dim), dtype=complex))
            J = sum((np.conj(ai) * (Ac @ X.T).T for ai, (_, _, Ac) in zip(a, ops) if ai != 0),
                    np.zeros((dim, dim), dtype=complex))
            d += J
            rates.append(np.trace(J))
        return np.concatenate([d.ravel(), np.array(rates, dtype=complex)])

    y = np.concatenate([net.initial.ravel(), np.zeros(nd + npt, dtype=complex)])
    grid = list(np.linspace(0.0, net.duration, samples + 1)) if samples else []
    rows, cum_at = [], {}
    bp = net.breakpoints
    nfev = 0

    def record(t, yy):
        rho = yy[:dim * dim].reshape(dim, dim)
        row = {"t": float(t), "excitations": float(np.trace(net.excitations @ rho).real)}
        for k, M in net.observables.items():
            row[k] = float(np.trace(M @ rho).real)
        for j, port in enumerate(net.ports):
            row[f"emitted_{port.name}"] = float(yy[dim * dim + nd + j].real)
        rows.append(row)

    for a, b in zip(bp, bp[1:]):
        inside = [t for t in grid if a <= t < b]
        t_eval = inside + [b]                     # the end of every interval is evaluated: the state carried on
        sol = solve_ivp(rhs, (a, b), y, method="DOP853", rtol=rtol, atol=atol, t_eval=t_eval, args=(0.5 * (a + b),))
        if not sol.success:
            raise RuntimeError(sol.message)
        nfev += int(sol.nfev)
        for k in range(sol.y.shape[1]):
            tk = float(sol.t[k])
            on_grid = any(abs(tk - g) <= 1e-12 * max(1.0, net.duration) for g in grid)
            if on_grid and (not rows or tk > rows[-1]["t"]):
                record(tk, sol.y[:, k])
        y = sol.y[:, -1]
        cum_at[b] = y[dim * dim + nd:].real.copy()
    cum_at[bp[0]] = np.zeros(npt)
    rho = y[:dim * dim].reshape(dim, dim)
    rho = 0.5 * (rho + rho.conj().T)
    return {"rho": rho, "absorbed": y[dim * dim:dim * dim + nd].real, "emitted": y[dim * dim + nd:].real,
            "emitted_at": cum_at, "trace": rows, "nfev": nfev}


# ------------------------------------------------------------------------------------------------ measurement
def measure(net: Network, samples: int = 200) -> Dict[str, object]:
    ev = evolve(net, samples=samples)
    rho = ev["rho"]
    tr = float(np.trace(rho).real)
    eig_min = float(np.linalg.eigvalsh(rho).min())
    ports = []
    supplied = 0.0
    detected_total = 0.0
    for j, port in enumerate(net.ports):
        row = {"port": port.name, "operator": port.operator, "rate": port.rate,
               "emitted_into_line": float(ev["emitted"][j])}
        if port.source is not None:
            row["input_kind"] = port.source["kind"]
            row["input_photons"] = float(port.source["photons"])
            supplied += row["input_photons"]
            if port.source["kind"] == "fock":
                row["left_in_source"] = float(np.trace(port.source["a"].conj().T @ port.source["a"] @ rho).real)
        if port.detector is not None:
            n_det = float(np.trace(port.detector["a"].conj().T @ port.detector["a"] @ rho).real)
            row["photons_in_output_mode"] = n_det
            detected_total += n_det
            if port.source is not None:
                row["efficiency"] = n_det / row["input_photons"]
        ports.append(row)
    windows = []
    for w in net.windows:
        j = [p.name for p in net.ports].index(w["port"])
        at = ev["emitted_at"]
        windows.append({**w, "photons": float(at[w["to"]][j] - at[w["from"]][j])})
    remaining = float(np.trace(net.excitations @ rho).real)
    absorbed = {d["name"]: float(a) for d, a in zip(net.dissipators, ev["absorbed"])}
    accounted = remaining + sum(absorbed.values()) + sum(p["emitted_into_line"] for p in ports) + detected_total \
        + sum(p.get("left_in_source", 0.0) for p in ports)
    initial_exc = float(np.trace(net.excitations @ net.initial).real)
    return {"name": net.name, "field": net.field, "carrier": net.carrier.description, "dimension": net.carrier.dim,
            "duration": net.duration, "ports": ports, "windows": windows, "removed_by_losses": absorbed,
            "remaining_in_system": remaining,
            "budget": {"supplied": supplied + initial_exc, "accounted": accounted,
                       "difference": accounted - supplied - initial_exc},
            "final": {k: float(np.trace(M @ rho).real) for k, M in net.observables.items()},
            "trace": ev["trace"], "checks": {"trace": tr, "smallest_eigenvalue": eig_min}, "nfev": ev["nfev"]}


# ------------------------------------------------------------------------------------------------ the swept write
def from_open(open_spec: Dict) -> Dict:
    """The one-mode swept write of ``fieldbridge-quantum-open/1`` written as a network: the sweep and the bias
    become controls, the mode keeps the name a (tokens a, ad, na) and every term and dissipator is unchanged."""
    pr = open_spec["protocol"]
    sw, bias = pr["sweep"], pr["bias"]
    T = (float(sw["to"]) - float(sw["from"])) / float(sw["rate"])
    hold = float(pr.get("hold", 0.0))
    total = T + hold
    off = bias.get("off_at")
    eps = [{"from": 0.0, "to": T, "value": f"{float(sw['from'])!r} + {float(sw['rate'])!r}*t"}]
    if hold > 0:
        eps.append({"from": T, "to": total, "value": repr(float(sw["to"]))})
    hval = repr(float(bias["value"]))
    h = [{"from": 0.0, "to": float(off) if off is not None and float(off) < total else total, "value": hval}]
    initial = {}
    th = (pr.get("initial") or {}).get("thermal")
    if th not in (None, 0, "0"):
        initial["a"] = {"thermal": th}
    return {"schema": SCHEMA, "name": open_spec["name"], "field": open_spec.get("field", "unspecified"),
            "question": open_spec["question"], "assumptions": open_spec["assumptions"],
            "provenance": open_spec["provenance"],
            "subsystems": [{"name": "a", "kind": "mode", "max_quanta": int(open_spec["carrier"]["max_quanta"])}],
            "parameters": open_spec.get("parameters") or {},
            "controls": {sw["parameter"]: eps, bias["parameter"]: h},
            "hamiltonian": open_spec["hamiltonian"], "dissipators": open_spec["dissipators"],
            "initial": initial, "duration": total}


# ------------------------------------------------------------------------------------------------ report
def markdown(report: Dict[str, object]) -> str:
    lines = [f"# {report['name']}", "", f"Field: {report['field']}. Carrier: {report['carrier']} "
             f"(dimension {report['dimension']}). Evolution from 0 to {report['duration']:g}.", ""]
    for p in report["ports"]:
        txt = f"- Port **{p['port']}** ({p['operator']}, rate {p['rate']:g}): "
        parts = []
        if "input_photons" in p:
            parts.append(f"{p['input_photons']:.4g} photons sent in ({p['input_kind']})")
        if "photons_in_output_mode" in p:
            parts.append(f"{p['photons_in_output_mode']:.4f} found in the output mode")
        if "efficiency" in p:
            parts.append(f"efficiency {p['efficiency']:.4f}")
        parts.append(f"{p['emitted_into_line']:.4f} emitted into the line"
                     + (" outside the output mode" if "photons_in_output_mode" in p else ""))
        lines.append(txt + "; ".join(parts) + ".")
    for w in report["windows"]:
        lines.append(f"- Photons emitted into {w['port']} between {w['from']:g} and {w['to']:g}: {w['photons']:.4f}.")
    if report["removed_by_losses"]:
        lines += ["", "Excitations removed by each loss (zero for pure dephasing, negative for a loss that adds "
                      "excitations, such as absorption from a thermal bath):"]
        lines += [f"- {k}: {v:.4f}" for k, v in report["removed_by_losses"].items()]
    b = report["budget"]
    lines += ["", f"Excitations left in the system: {report['remaining_in_system']:.4f}. Budget: "
              f"{b['supplied']:.4f} supplied, {b['accounted']:.4f} accounted for (difference {b['difference']:.2e}; "
              "it vanishes when the Hamiltonian conserves the number of excitations).",
              f"Checks: trace {report['checks']['trace']:.10f}, smallest eigenvalue "
              f"{report['checks']['smallest_eigenvalue']:.2e}."]
    if report["final"]:
        lines += ["", "Final expectation values: " + ", ".join(f"{k} = {v:.4f}" for k, v in report["final"].items()) + "."]
    return "\n".join(lines) + "\n"
