"""Data of the web page: the realizations of site_registry, computed with FieldBridge.

build() resolves every node of the registry to a specification, computes it with the functions of the quantum and
memory packages, checks that every edge changes only the component it names (check_edit), and fills the texts of the
edges and the sequences with the computed facts. The page (fieldbridge/web) recomputes the dynamics in the browser;
the mechanism, the derivation words and the law constants shown for a realization come from here.

Law constants need long simulations (twelve minutes for phase locking on a laptop). build(law=True) computes them
and write_law_record() stores them in docs/site/law_constants.json with the hash of every specification; a build
without law takes each constant from that record when the specification is unchanged, and shows none otherwise.
"""
from __future__ import annotations

import copy
import hashlib
import json
import platform
import re
import zlib
from math import atan2, degrees, pi
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional

import numpy as np

from . import site_references as refs
from . import site_registry as reg
from .quantum.language import closure_basis, observable_frequencies
from .web_models import arithmetic, model_record

ROOT = Path(__file__).resolve().parents[1]
LAW_RECORD = ROOT / "docs" / "site" / "law_constants.json"
SCHEMA = "fieldbridge-site/1"
TARGETS = ("symmetric-write", "threshold-write", "phase-locking")
PREFIX = {"symmetric-write": "sym", "threshold-write": "thr", "phase-locking": "lock"}


class SiteError(ValueError):
    """The registry and the calculations disagree: an edge does not change the component it names, a text refers to
    a fact that was not computed, or a specification does not load."""


# ------------------------------------------------------------------------------------------------ helpers
def _rng(nid: str, salt: str) -> np.random.Generator:
    """A generator seeded by the node and the calculation, so that a subset of nodes gives the same results."""
    return np.random.default_rng(zlib.crc32(f"{nid}:{salt}".encode()))


def _strip(html: str) -> str:
    return re.sub(r"<[^>]+>", "", html).replace("&gt;", ">").replace("&lt;", "<").replace("&amp;", "&")


def _mat(M: np.ndarray) -> Dict[str, List[float]]:
    M = np.asarray(M, complex)
    return {"n": int(M.shape[0]), "re": [float(f"{v:.12g}") for v in M.real.ravel()],
            "im": [float(f"{v:.12g}") for v in M.imag.ravel()]}


def _num(v, spec: str = "") -> str:
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        v = float(v)
        if spec:
            return format(v, spec)
        a = abs(v)
        if a < 5e-7:  # a write point located at zero to the tolerance of the search
            return "0"
        text = f"{v:.0f}" if a >= 100 else f"{v:.1f}" if a >= 10 else f"{v:#.3g}".rstrip(".")
        return text.replace("-", "−")
    if isinstance(v, (list, tuple)):
        return ", ".join(_num(x, spec) for x in v)
    return str(v)


PLACEHOLDER = re.compile(r"\{(from\.)?([a-z_0-9]+)(?::([^}]+))?\}")


def fill(text: str, facts: Dict, facts_from: Optional[Dict] = None, strict: bool = True) -> str:
    """Replace {key}, {from.key} and {key:spec} by computed facts."""
    def one(m):
        src = facts_from if m.group(1) else facts
        key = m.group(2)
        if src is None or key not in src or src[key] is None:
            if strict:
                raise SiteError(f"the text refers to '{m.group(0)}', which was not computed")
            return "—"
        return _num(src[key], m.group(3) or "")
    return PLACEHOLDER.sub(one, text)


def _memory_hash() -> str:
    here = Path(__file__).resolve().parent / "memory"
    return hashlib.sha256(b"".join(p.read_bytes() for p in sorted(here.glob("*.py")))).hexdigest()


def spec_hash(spec: Dict) -> str:
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()


# ------------------------------------------------------------------------------------------------ resolution
def _modify(spec: Dict, d: Dict, family: str) -> None:
    if family == "unitary":
        if "params" in d:
            spec.setdefault("parameters", {}).update(d["params"])
        if "observable" in d:
            spec["observable"] = copy.deepcopy(d["observable"])
        for k, term in (d.get("replace_terms") or {}).items():
            spec["hamiltonian"][int(k)] = copy.deepcopy(term)
        for k, value in (d.get("coefficients") or {}).items():
            spec["hamiltonian"][int(k)]["coefficient"] = value
        spec["hamiltonian"] = spec["hamiltonian"] + copy.deepcopy(d.get("add_terms", []))
        return
    if "params" in d:
        spec["parameters"].update(d["params"])
    if "drift" in d:
        spec["drift"].update(d["drift"])
    if "potential" in d:
        if d["potential"] is None:
            spec.pop("potential", None)
        else:
            spec["potential"] = d["potential"]
    if "control" in d:
        spec["control"] = copy.deepcopy(d["control"])


def _clean(spec: Dict) -> Dict:
    """Coefficients of order 1e-16 left by cos(90 degrees) in an attached specification are set to zero; the term is
    kept, so that the operators of the realization do not change."""
    for key in ("hamiltonian", "observable"):
        for term in spec.get(key, []):
            c = term.get("coefficient")
            if isinstance(c, (int, float)) and abs(c) < 1e-12:
                term["coefficient"] = 0.0
    return spec


def auto_nodes(root: Path = ROOT) -> List[Dict]:
    """Specifications in the folders of AUTO_DIRS that no node of the registry names: each becomes a node of its own,
    so that a contributed material appears on the page without an edit of the registry."""
    named = {d["spec"] for d in reg.NODES if "spec" in d}
    ids = {d["id"] for d in reg.NODES}
    out = []
    for folder, family in reg.AUTO_DIRS.items():
        for path in sorted((root / folder).glob("*.json")):
            rel = str(path.relative_to(root))
            if rel in named:
                continue
            spec = json.loads(path.read_text(encoding="utf-8"))
            fam = "field" if spec.get("kind") == "field" else family
            nid = path.stem if path.stem not in ids else f"{path.stem}_{Path(folder).name}"
            out.append({"id": nid, "family": fam, "spec": rel, "auto": True,
                        "tutorial": ("docs/tutorial/28_regulation_set_point.md#8-a-model-from-your-field"
                                     if fam == "regulation" else "docs/tutorial/21_memory_new_material.md"
                                     if fam != "unitary" else "docs/tutorial/24_spin_language.md#2-writing-a-realization")})
            ids.add(nid)
    return out


def all_defs(root: Path = ROOT) -> List[Dict]:
    return list(reg.NODES) + auto_nodes(root)


def resolve(root: Path = ROOT, only: Optional[Iterable[str]] = None) -> Dict[str, Dict]:
    """Every node as {"id", "def", "family", "spec", "path", "attach"}, bases before the nodes built from them."""
    defs = {d["id"]: d for d in all_defs(root)}
    out: Dict[str, Dict] = {}

    def visit(nid: str) -> Dict:
        if nid in out:
            return out[nid]
        d = defs[nid]
        attach = None
        if "spec" in d:
            spec = json.loads((root / d["spec"]).read_text(encoding="utf-8"))
            family, path = d["family"], d["spec"]
        elif "attach" in d:
            from .quantum import language as ql
            src_id, to, size = d["attach"]
            src = visit(src_id)
            attach = ql.attach(ql.load(src["spec"]), to, size)
            if not attach["attached"]:
                raise SiteError(f"{nid}: the rotation of {src_id} does not attach to {to}")
            spec, family, path = _clean(copy.deepcopy(attach["spec"])), "unitary", src["path"]
        else:
            base = visit(d["base"])
            spec, family, path = copy.deepcopy(base["spec"]), base["family"], base["path"]
            _modify(spec, d, family)
        if "name" in d:
            spec["name"] = _strip(d["name"])
        if "question" in d:
            spec["question"] = d["question"]
        out[nid] = {"id": nid, "def": d, "family": family, "spec": spec, "path": path, "attach": attach,
                    "derived": "spec" not in d, "name_html": d.get("name") or spec.get("name", nid)}
        return out[nid]

    for nid in (list(only) if only is not None else defs):
        visit(nid)
    return out


# ------------------------------------------------------------------------------------------------ unitary
# the closure of the observable (Chapter 11) and the frequencies with which it moves: the letter O of the derivation
def observable_closure(H: np.ndarray, O: np.ndarray, tol: float = 1e-9, max_dim: int = 256) -> int:
    """The dimension of the closure of O under i[H, .]."""
    return len(closure_basis(H, O, tol, max_dim))


def observable_frame(canon: Dict) -> List[np.ndarray]:
    """The rotating operators in the frame of the observable: F3 along the observable, the rotation axis in the plane
    of F1 and F3 at the angle theta from F3. canonical_su2 and canonical_closure give J with J3 along the axis and
    the observable in the plane of J1 and J3."""
    J, th = canon["J"], float(canon["theta"])
    c, s = np.cos(th), np.sin(th)
    return [-c * J[0] + s * J[2], -J[1], s * J[0] + c * J[2]]


def _operator_matrix(real, e: Dict) -> np.ndarray:
    M = real.carrier.operator(e["operator"])
    return M + M.conj().T if e.get("hc") else M


def _product_state(product: List[str]) -> np.ndarray:
    one = {"+z": np.array([1, 0], complex), "-z": np.array([0, 1], complex),
           "+x": np.array([1, 1], complex) / np.sqrt(2), "-x": np.array([1, -1], complex) / np.sqrt(2),
           "+y": np.array([1, 1j], complex) / np.sqrt(2), "-y": np.array([1, -1j], complex) / np.sqrt(2)}
    psi = np.array([1], complex)
    for p in product:
        psi = np.kron(psi, one[p])
    return psi


def _slider(name: str, value: float) -> List[float]:
    span = max(2.0 * abs(value), 1.0)
    return [-span, span, span / 100.0]


def unitary_record(node: Dict, parent: Optional[Dict] = None) -> Dict:
    from .quantum import language as ql
    spec = node["spec"]
    real = ql.load(spec)
    V, sec = ql.restrict(real)
    Hs, Os = V.conj().T @ real.H @ V, V.conj().T @ real.O @ V
    row = ql.derive_bloch_rotation(real)
    closure = observable_closure(Hs, Os)
    basis, closed = ql.lie_closure([Hs, Os])
    params = {k: float(v) for k, v in (spec.get("parameters") or {}).items()}
    in_terms = set()
    terms = []
    for k, e in enumerate(spec["hamiltonian"]):
        tree = arithmetic(str(e["coefficient"]))
        in_terms |= {m for m in re.findall(r"[A-Za-z_]\w*", str(e["coefficient"])) if m in params}
        terms.append({"id": f"t{k}", "operator": e["operator"], "hc": bool(e.get("hc")),
                      "coefficient": str(e["coefficient"]), "tree": tree,
                      "matrix": _mat(V.conj().T @ _operator_matrix(real, e) @ V)})
    frame, j_top, frame_from = None, None, None
    reached = ql.reached(row)  # through the algebra, or through the closure of the observable
    if reached:
        H0, O0 = ql._traceless(Hs), ql._traceless(Os)
        canon = (ql.canonical_su2(basis, H0, O0) if row["status"] == ql.REACHED
                 else ql.canonical_closure(closure_basis(Hs, O0), Hs, O0))
        frame, j_top, frame_from = observable_frame(canon), float(row["law"]["j_top"]), node["id"]
    elif parent is not None and parent["engine"].get("frame") and parent["engine"]["dim"] == Hs.shape[0]:
        frame, j_top, frame_from = parent["_frame"], parent["engine"]["j_top"], parent["id"]
    # the term without which H and the observable generate su(2): named also when the rotation is reached through
    # the closure, where it is the reason for the larger algebra and not an obstruction
    cause = ql._single_term_cause(real, V) if row["status"] != ql.REACHED and len(real.terms) > 1 else None
    sig = row.get("signature", {})
    freqs = observable_frequencies(Hs, Os)
    facts = {"status": row["status"], "word": row["word"], "dim": int(row["algebra_dimension"]),
             "closure": closure, "frequencies": len(freqs), "frequencies_text": NUMBER_WORDS.get(len(freqs), str(len(freqs))),
             "frequency_list": [float(f"{f:.6g}") for f in freqs[:6]],
             "carrier": real.carrier.description, "hilbert": int(real.carrier.dim), "sector": int(Hs.shape[0])}
    if reached:
        facts.update(rate=float(sig["rate"]), theta=float(sig["theta_deg"]))
        if row["representation"]:  # a closure that is not an su(2) carries no representation
            facts["rep"] = ql._rep(row["representation"]).replace(" x ", " × ").replace("j=", "j = ")
        facts["residual"] = float(row["law"]["residual"])
    if cause:
        facts["cause"] = _cause_html(cause)
    if node.get("attach"):
        native = node["attach"]["native"]
        for key in ("tunnelling", "energy_difference", "pairing_amplitude", "single_particle_energy",
                    "field_x", "field_z", "rabi_frequency", "detuning", "ising_coupling_g", "transverse_field_h"):
            if key in native:
                facts[key] = float(native[key])
        if "exchange_couplings" in native:
            facts["couplings"] = [float(x) for x in native["exchange_couplings"]]
    klass = "rotation" if reached else "conserved" if closure == 1 else "obstructed"
    preps = []
    for p in node["def"].get("preparations", []):
        if "product" in p:
            psi = V.conj().T @ _product_state(p["product"])
            preps.append({"id": p["id"], "label": p["label"], "state": {"re": psi.real.tolist(), "im": psi.imag.tolist()}})
        else:
            preps.append({"id": p["id"], "label": p["label"]})
    sector_params = set(params) - in_terms
    engine = {"kind": "unitary", "dim": int(Hs.shape[0]), "terms": terms,
              "observable": {"text": " + ".join(f"{t['coefficient']} {t['operator']}" for t in spec["observable"]),
                             "terms": [{"coefficient": str(t["coefficient"]), "operator": t["operator"]}
                                       for t in spec["observable"]], "matrix": _mat(Os)},
              "params": {k: v for k, v in params.items() if k not in sector_params},
              "fixed": {k: params[k] for k in sector_params},
              "ranges": {k: _slider(k, v) for k, v in params.items() if k not in sector_params},
              "frame": [_mat(F) for F in frame] if frame is not None else None, "frame_from": frame_from,
              "j_top": j_top, "preparations": preps,
              "carrier": {**real.carrier.spec, "description": real.carrier.description,
                          "sector": {"operator": sec["operator"], "value": sec["value"], "dimension": sec["dimension"],
                                     "full_dimension": sec["full_dimension"]} if sec else None}}
    law = None
    if reached:
        L = row["law"]
        law = {"times": L["times"][::2], "f": [float(f"{v:.6g}") for v in L["f_exact"][::2]]}
    record = {"id": node["id"], "family": "unitary", "name": node["name_html"], "field": real.field,
              "class": klass, "facts": facts, "engine": engine, "law_curve": law,
              "derivation": {"word": row["word"], "steps": [{k: v for k, v in s.items() if k in ("letter", "dimension", "closed", "operator", "value", "rate", "theta_deg", "residual")}
                                                            for s in row["steps"]],
                             "obstruction": row.get("obstruction")},
              "slots": {"Omega": "H = " + " + ".join(f"{t['coefficient']} [{t['operator']}{' + h.c.' if t.get('hc') else ''}]" for t in spec["hamiltonian"]),
                        "Xi": real.carrier.description + (f"; sector {sec['operator']} = {sec['value']:g} ({sec['dimension']} states)" if sec else ""),
                        "C": "closed unitary evolution, ħ = 1" + ("; the sector is conserved" if sec else ""),
                        "R": engine["observable"]["text"],
                        "P": "the top eigenstate of R is prepared; H acts from t = 0",
                        "A": ", ".join(f"{k} = {_num(v)}" for k, v in params.items()) or "numerical coefficients"}}
    record["_frame"] = frame
    return record


def _cause_html(label: str) -> str:
    """'h [X1]' -> 'h X<sub>1</sub>' for the texts."""
    m = re.match(r"(.*?) \[(.*)\]", label)
    coef, op = (m.group(1), m.group(2)) if m else ("", label)
    op = op.replace(" + h.c.", "")
    for k, v in GREEK.items():
        coef = re.sub(rf"\b{k}\b", v, coef)
    body = operator_html(op)
    return (coef + " " if coef else "") + (f"({body})" if " + " in op and coef else body)


NUMBER_WORDS = {0: "no", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six"}
GREEK = {"lam": "λ", "Delta": "Δ", "delta": "δ", "eps": "ε", "mu": "μ", "alpha": "α", "gamma": "γ", "sigma": "σ",
         "kappa": "κ", "theta": "θ", "psi": "ψ", "phi": "φ", "rabi": "Ω<sub>R</sub>", "detuning": "δ", "rf": "Ω<sub>rf</sub>"}


def operator_html(text: str) -> str:
    def token(t: str) -> str:
        m = re.fullmatch(r"([XYZ])(\d)", t)
        if m:
            return f"{m.group(1)}<sub>{m.group(2)}</sub>"
        m = re.fullmatch(r"J([xyzpm])", t)
        if m:
            return "J<sub>" + {"p": "+", "m": "−"}.get(m.group(1), m.group(1)) + "</sub>"
        m = re.fullmatch(r"cd(\d)", t)
        if m:
            return f"c<sup>†</sup><sub>{m.group(1)}</sub>"
        m = re.fullmatch(r"[cn](\d)", t)
        if m:
            return f"{t[0]}<sub>{m.group(1)}</sub>"
        m = re.fullmatch(r"n([a-z])", t)
        if m:
            return f"n<sub>{m.group(1)}</sub>"
        m = re.fullmatch(r"([a-z])d", t)
        if m:
            return f"{m.group(1)}<sup>†</sup>"
        return t

    def product(p: str) -> str:
        toks, out, k = p.split(), [], 0
        while k < len(toks):
            n = 1
            while k + n < len(toks) and toks[k + n] == toks[k]:
                n += 1
            out.append(token(toks[k]) + (f"<sup>{n}</sup>" if n > 1 else ""))
            k += n
        return "".join(out)

    return " + ".join(product(p.strip()) for p in text.split("+"))


# ------------------------------------------------------------------------------------------------ dissipative
# trajectories of the swept-write check when the law constants of the record are computed: the standard error of the
# constant is then 0.014 to 0.018, about 1% of pi^(1/4); with the 400 of a quick check it is 0.11
LAW_TRAJECTORIES = 25600
# independent runs into which those trajectories are split: the error of a constant is the larger of their scatter and
# the binomial error
LAW_REPLICATES = 8


def _derive(target: str, real, nid: str, law: bool) -> Dict:
    from .memory import codiscovery
    return codiscovery.TARGETS[target](real, _rng(nid, target), check_law=law, n_traj=LAW_TRAJECTORIES,
                                       replicates=LAW_REPLICATES)


def _row_facts(prefix: str, row: Dict) -> Dict:
    out = {f"{prefix}_status": row["status"], f"{prefix}_word": row.get("word") or "—",
           f"{prefix}_class": row.get("class", "")}
    if row.get("obstruction_short") or row.get("obstruction"):
        out[f"{prefix}_obstruction"] = row.get("obstruction_short") or row.get("obstruction")
    const = row.get("law_constant") or {}
    if const.get("constant") is not None:
        out[f"{prefix}_law"] = float(const["constant"])
        out[f"{prefix}_law_err"] = float(const["stderr"])
    if prefix == "lock" and row["status"].startswith("reached"):
        canon = row["canonical"]
        out.update(lock_ratio=int(canon["ratio"]), lock_K=float(canon["K"]), lock_period=float(canon["period"]))
        if const.get("constant") is not None:
            out["lock_width"], out["lock_width_err"] = out["lock_law"], out["lock_law_err"]
    return out


def dissipative_record(node: Dict, rows: Dict[str, Dict]) -> Dict:
    from .memory import analysis as an
    from .memory import spec as mspec
    spec = node["spec"]
    real = mspec.load(spec)
    nid = node["id"]
    facts: Dict[str, object] = {}
    for target, row in rows.items():
        facts.update(_row_facts(PREFIX[target], row))
    over: Dict[str, float] = {}
    lock = rows.get("phase-locking")
    if lock and lock["status"].startswith("reached") and lock.get("write_point"):
        wp = lock["write_point"]
        if abs(float(wp["value"]) - float(real.params[wp["param"]])) > 1e-9:
            over = {wp["param"]: float(wp["value"])}
            facts["operating_param"], facts["operating_value"] = wp["param"], float(wp["value"])
    states, counts, spectra, n_unconv, continuum = an.stored_states(real, _rng(nid, "states"), 48, over)
    facts["states"] = len(states)
    facts["states_text"] = NUMBER_WORDS.get(len(states), str(len(states)))
    facts["loss"] = an.loss_law(len(states), spectra, n_unconv, 48, continuum)
    scale = an.is_scale_control(real, _rng(nid, "scale")) if real.control else False
    events = []
    if real.control and not scale and not over and len(states) >= 1:
        lo, hi = real.control_range
        for ev in an.locate_writes(real, _rng(nid, "writes"), values=np.linspace(lo, hi, 7)):
            try:
                kind = an.normal_form(real, np.asarray(ev["q"], float), real.control, ev["v"]).get("kind", "")
            except (np.linalg.LinAlgError, ValueError):
                kind = ""
            events.append({"value": float(ev["v"]), "type": ev["type"], "kind": kind.split(":")[0], "q": ev["q"]})
    if events:
        facts.update(write_point=events[0]["value"], write_param=real.control, write_kind=events[0]["kind"])
    facts["scale_control"] = bool(scale)
    klass = _dissipative_class(states, n_unconv, events, scale, lock)
    if klass == "neutral-cycles":  # the wording of the memory card (memory.discovery) once the cycles are neutral
        facts["loss"] = ("no stable state: a family of neutral cycles, with no restoring force on the amplitude or the "
                         "phase (Law 2)")
    engine = model_record(spec, nid, node["path"])
    for k in ("tutorial", "saved", "source_url", "id", "name", "question", "field", "specification",
              "assumptions", "source", "closure", "observable"):
        engine.pop(k, None)
    engine["kind_dynamics"] = engine["kind"]  # equations, gene or rotor; model-physics.js reads "kind" as well
    engine["params"] = {**engine["params"], **over}
    for k, v in over.items():
        lo, hi, st = engine["ranges"][k]
        engine["ranges"][k] = [min(lo, v), max(hi, v), st]
    engine["control"] = {"name": real.control, "range": list(real.control_range)} if real.control else None
    engine["events"] = [{k: v for k, v in e.items() if k != "q"} for e in events]
    engine["states"] = [np.asarray(s).tolist() for s in states]
    engine["potential"] = arithmetic(spec["potential"]) if spec.get("potential") else None
    engine["noise"] = float(real.noise)
    if lock and lock["status"].startswith("reached"):
        canon = lock["canonical"]
        engine["drive"] = {"param": lock["write_point"]["param"], "p1": float(lock["write_point"]["value"]),
                           "ratio": int(canon["ratio"]), "omega": float(canon["omega"]), "K": float(canon["K"]),
                           "eps": float(canon["eps"]), "period": float(canon["period"])}
    carrier = spec["carrier"] if spec.get("kind") == "equations" else {"kind": real.carrier.kind, "name": real.carrier.name}
    record = {"id": nid, "family": "dissipative", "name": node["name_html"], "field": spec.get("field", ""),
              "class": klass, "facts": facts, "engine": engine,
              "question": spec.get("question"), "assumptions": spec.get("assumptions", []),
              "rows": {t: _row_summary(r) for t, r in rows.items()},
              "slots": {"Omega": spec.get("name") if spec.get("kind") == "network" else
                        "; ".join(f"d{v}/dt = {spec['drift'][v]}" for v in spec["carrier"]["variables"]),
                        "Xi": real.carrier.describe(),
                        "C": spec.get("closure") or real.closure or "thermal bath",
                        "R": spec.get("observable") or real.observable,
                        "P": ("bounded write toward a stored state; sweep of " + real.control) if real.control
                        else "bounded write toward a stored state",
                        "A": ", ".join(f"{k} = {_num(float(v))}" for k, v in {**real.params, **over}.items())
                        + f"; noise D = {_num(float(real.noise))}"},
              "carrier_spec": carrier}
    return record


def _short_name(name: str) -> str:
    words = _strip(name).split()
    out = ""
    for w in words:
        if len(out) + len(w) + 1 > 18:
            break
        out = (out + " " + w).strip()
    return out or _strip(name)[:18]


def _row_summary(row: Dict) -> Dict:
    out = {k: row.get(k) for k in ("status", "word", "class", "obstruction_short") if row.get(k) is not None}
    canon = row.get("canonical") or {}
    if "x" in canon and "g" in canon:
        x, g = np.asarray(canon["x"], float), np.asarray(canon["g"], float)
        stride = max(1, x.size // 48)
        out["canonical"] = {"x": [float(f"{v:.5g}") for v in x[::stride]], "g": [float(f"{v:.5g}") for v in g[::stride]]}
    if "phi" in canon and "g" in canon:
        # the phase gained per period, in units of K, against the canonical phase: -sin(phi) for the Adler form
        phi, g = np.asarray(canon["phi"], float), np.asarray(canon["g"], float)
        stride = max(1, phi.size // 64)
        out["canonical"] = {"phi": [float(f"{v:.5g}") for v in phi[::stride]], "g": [float(f"{v:.5g}") for v in g[::stride]],
                            "ratio": int(canon["ratio"])}
    const = row.get("law_constant") or {}
    if const.get("constant") is not None:
        out["law"] = {"constant": float(const["constant"]), "stderr": float(const["stderr"]),
                      "expected": float(const.get("expected", float("nan"))),
                      **({"signature": const["signature"]} if const.get("signature") else {})}
    if row.get("write_point"):
        out["write_point"] = {"param": row["write_point"]["param"], "value": float(row["write_point"]["value"])}
    return out


def _dissipative_class(states, n_unconv, events, scale, lock) -> str:
    if lock and "neutral" in str(lock.get("obstruction", "")):
        return "neutral-cycles"
    if not states:
        return "oscillation" if n_unconv else "single-state"
    if len(states) == 1:
        return "single-state"
    if scale:
        return "field-write"
    kind = events[0]["kind"] if events else ""
    if kind.startswith("supercritical pitchfork"):
        return "symmetric-write"
    if kind.startswith("saddle-node"):
        return "threshold-write"
    if kind.startswith("subcritical"):
        return "subcritical-write"
    return "field-write"


# ------------------------------------------------------------------------------------------------ fields, stochastic
def field_record(node: Dict) -> Dict:
    from .memory import fields as mf
    spec = node["spec"]
    model = mf.from_spec(spec)
    law = mf.predicted_law(model)
    times = np.geomspace(0.1, 3000.0, 48)
    curve = mf.spectral_snr(model, times)
    facts = {"law": law["law"]}
    if law["law"] == "power":
        from fractions import Fraction
        e = Fraction(float(law["exponent"])).limit_denominator(4)
        facts["exponent"] = ("−" if e < 0 else "") + (str(abs(e.numerator)) if e.denominator == 1
                                                     else f"{abs(e.numerator)}/{e.denominator}")
        facts["exponent_value"] = float(law["exponent"])
        facts["exponent_fit"] = float(mf.exponent_check(model)["exponent"])
    else:
        facts["rate"] = float(law["rate"])
    f = spec["field"]
    return {"id": node["id"], "family": "field", "name": node["name_html"], "field": "statistical physics",
            "class": "power-loss" if law["law"] == "power" else "exponential-loss", "facts": facts,
            "question": spec.get("question"), "assumptions": spec.get("assumptions", []),
            "engine": {"kind": "field", **{k: f.get(k) for k in ("d", "L", "conserved", "write", "amplitude", "M", "T",
                                                                    "g", "kappa0", "Dgrad")},
                       "kappa0": model.kappa0, "Dgrad": model.Dgrad,
                       "curve": {"t": [float(f"{v:.6g}") for v in curve["times"]],
                                 "snr": [float(f"{v:.6g}") for v in curve["snr_profile"]]}},
            "slots": {"Omega": "relaxation in F = Σ φ²/2 + g φ⁴/4" + (" with a gradient term" if not model.conserved else ""),
                      "Xi": f"lattice field in d = {model.d}, {model.L}^{model.d} sites",
                      "C": "conserved density (Model B)" if model.conserved else "no conservation law (Model A)",
                      "R": "the profile of the write, measured against a field written without it",
                      "P": "a write that adds material at one site" if model.write == "charge"
                      else "a write that moves material between two sites",
                      "A": ", ".join(f"{k} = {_num(float(getattr(model, k)))}" for k in
                                     (("M", "T", "g", "amplitude") if model.conserved
                                      else ("M", "T", "g", "amplitude", "kappa0", "Dgrad")))}}


def hysteron_record(node: Dict) -> Dict:
    """Interacting hysterons under a slow drive: the prediction from the signs of the couplings and of the drive, the
    subloops measured as by ``memory hysterons`` (same seed and count), the major loop and one excursion."""
    from .memory import hysterons as hy
    spec = node["spec"]
    model = hy.from_spec(spec)
    c = hy.card(model, np.random.default_rng(HYSTERON_SEED), count=HYSTERON_SUBLOOPS)
    result, cls, loop, shown = c["check"], c["class"], c["loop"], c["excursion"]
    pred = result["prediction"]
    lo, hi = c["window"]
    moving = [H for branch in (loop["up"], loop["down"]) for (H, r), (_, r0) in zip(branch[1:], branch) if r != r0]
    span = max(moving) - min(moving) if moving else 1.0
    xr = [min(moving) - 0.08 * span, max(moving) + 0.08 * span] if moving else [lo - 1, hi + 1]
    keep = _thin(len(shown["trace"]), 500)
    facts = {"elements": model.n, "guaranteed": "yes" if pred["guaranteed"] else "no",
             "failed": result["failed"], "total": result["total"], "bonds": pred["bonds"],
             "bonds_frustrated": pred["bonds_frustrated_with_drive"],
             "couplings_frustrated": "no" if pred["couplings_balanced"] else "yes",
             "plaquettes": "—" if pred["frustrated_plaquettes"] is None else _num(pred["frustrated_plaquettes"]),
             "max_differ": max((r["max_differ"] for r in result["subloops"].values()), default=0),
             "loss": "no thermal loss on the time of the drive (Law 3 neglected)",
             "statement": pred["statement"], "verdict": result["verdict"]}
    return {"id": node["id"], "family": "hysterons", "name": node["name_html"], "field": spec.get("field", ""),
            "class": cls, "facts": facts, "question": spec.get("question"), "assumptions": spec.get("assumptions", []),
            "engine": {"kind": "hysterons", "elements": model.n, "xr": [float(f"{v:.6g}") for v in xr], "window": [lo, hi],
                       "major": {k: _trace_points(loop[k], 500) for k in ("up", "down")},
                       "excursion": {"H1": shown["H1"], "H2": shown["H2"], "differ": shown["differ"],
                                     "trace": [_pair(shown["trace"][i]) for i in keep],
                                     "progress": [float(f"{shown['progress'][i]:.5g}") for i in keep],
                                     "unlike": [shown["unlike"][i] for i in keep]}},
            "slots": hysteron_slots(spec, model)}


def _thin(n: int, most: int) -> List[int]:
    """At most `most` indices of n, the first and the last kept."""
    if n <= most:
        return list(range(n))
    return sorted(set(np.linspace(0, n - 1, most).round().astype(int).tolist()))


def _pair(p) -> List[float]:
    return [float(f"{p[0]:.6g}"), float(f"{p[1]:.6g}")]


def _trace_points(trace: List, most: int) -> List[List[float]]:
    return [_pair(trace[i]) for i in _thin(len(trace), most)]


def _amount(item, what: str) -> str:
    if item is None or item == 0:
        return f"no {what}"
    if isinstance(item, (int, float)):
        return f"{what} {_num(float(item))}"
    if isinstance(item, list):
        return f"{what} given element by element"
    if item.get("distribution") == "gaussian":
        return f"Gaussian {what} of width {_num(float(item['width']))}"
    return f"{what} uniform between {_num(float(item['low']))} and {_num(float(item['high']))}"


def hysteron_slots(spec: Dict, model) -> Dict[str, str]:
    hs = spec["hysterons"]
    lat, cp = hs.get("lattice"), hs.get("couplings")
    if lat:
        xi = (f"{model.n} elements on a periodic square lattice of {lat['L']} × {lat['L']}" if lat["shape"] == "square"
              else f"{model.n} elements on a periodic chain")
        p = float((cp or {}).get("negative_fraction", 0.0))
        bonds = ("couplings of strength " + _num(float((cp or {}).get("value", 1.0))) + " between neighbours, "
                 + ("all positive" if p == 0 else "all negative" if p == 1 else f"a fraction {_num(p)} of them negative"))
    else:
        xi = f"{model.n} elements"
        bonds = "no couplings" if not cp else f"{len(cp)} couplings" + ("" if model.reciprocal else ", directed")
    drive = hs.get("drive", "uniform")
    drive_text = {"uniform": "a uniform drive", "staggered": "a drive of opposite signs on the two sublattices",
                  "random": "a drive of random sign on each element"}.get(drive if isinstance(drive, str) else "",
                                                                            "a drive with the given sign on each element")
    return {"Omega": "switching at thresholds of fᵢ = Σⱼ Jᵢⱼsⱼ + hᵢ + ηᵢH, with " + bonds,
            "Xi": xi + ", each in one of two states",
            "C": "a quasi-static drive without thermal activation; each avalanche relaxed at a fixed drive",
            "R": "the response Σ ηᵢsᵢ / Σ |ηᵢ|, and the configuration at each turning point",
            "P": drive_text + ", from a large negative value, with excursions inside its turning points",
            "A": "; ".join([_amount(hs.get("fields"), "random fields"), _amount(hs.get("half_widths"), "half-widths"),
                            f"seed {hs.get('seed', 0)}"])}


# the parameter points of a regulated node: fewer than ``regulation card`` (32), so that the page builds quickly
REGULATION_SAMPLES = 8
REGULATION_CLASSES = ("perfect-adaptation", "fine-tuned-adaptation", "partial-adaptation", "no-adaptation")


def regulation_record(node: Dict) -> Dict:
    """A regulated realization: the card of ``regulation card`` (with fewer parameter points) and the response of the
    output to the first step of the input."""
    from .regulation import card as rcard
    from .regulation import spec as rspec
    from .regulation.cli import _short_phi
    spec = node["spec"]
    real = rspec.load(spec)
    res = rcard.card(real, samples=REGULATION_SAMPLES, seed=0, first_step_only=True, trace=True)
    if res["class"] not in REGULATION_CLASSES:
        raise SiteError(f"{node['id']}: {res['class']}: a regulated realization on the page needs a stable steady "
                        f"state at its reference input")
    step = res["step_responses"][0]
    rob = res["robustness"]
    att = [a for a in res.get("attenuation") or [] if a.get("remaining_fraction") is not None]
    clamp = min(att, key=lambda a: abs(a["remaining_fraction"])) if att else None
    integ = res.get("integrator") or {}
    ratio = res.get("gain_ratio") or 0.0
    # without a variable to clamp the reference gain is the gain itself, and the ratio says nothing
    shown_ratio = "0" if ratio < 1e-9 else _num(ratio) if att else "— (no variable to clamp)"
    facts = {"integrator": _short_phi(integ), "gain_ratio": shown_ratio,
             "set_point": _num(res["y0_steady"]), "u0": _num(real.u0), "u1": _num(step["u1"]),
             "input": real.input, "output": real.output_text if len(real.output_text) <= 24 else "the observable",
             "peak": _num(step["peak"]), "final": "0" if abs(step["final"]) < 1e-9 else _num(step["final"]),
             "final_over_peak": "0" if step.get("final_over_peak") is None or abs(step["final_over_peak"]) < 1e-9
             else _num(step["final_over_peak"]),
             "return_time": _num(step["return_time"]),
             "robust": f"{rob['gain_zero']} of {rob['stable']}",
             "clamp": f"{clamp['variable']} ({clamp['role']}): {_num(clamp['remaining_fraction'])}" if clamp else "—",
             "calibration": "—" if step.get("calibration_ratio") is None
             else f"{abs(step['calibration_ratio'] - 1):.1e}"}
    tr = step["trace"]
    return {"id": node["id"], "family": "regulation", "name": node["name_html"], "field": spec.get("field", ""),
            "class": res["class"], "facts": facts, "question": spec.get("question"),
            "assumptions": spec.get("assumptions", []),
            "engine": {"kind": "regulation", "t": tr["t"], "y": tr["y"], "phi": tr["phi"],
                       "y_ref": step["y_ref"], "y_end": step["y_end"], "set_point": res["y0_steady"],
                       "u0": real.u0, "u1": step["u1"], "input": real.input},
            "slots": regulation_slots(spec, real, step)}


OPEN_CLASSES = ("linear-stage-write", "equilibrium-write")
OPEN_TRACE = 160


def open_record(node: Dict) -> Dict:
    """An open realization: the quantum write of ``quantum write`` with a trace of the amplified quadrature."""
    from .quantum import open as qo
    spec = node["spec"]
    real = qo.load(spec)
    row = qo.derive_quantum_write(real, law=True, trace=OPEN_TRACE)
    cls = (row.get("class") or "").replace(" ", "-")
    if cls not in OPEN_CLASSES:
        raise SiteError(f"{node['id']}: {row.get('class')}: an open realization on the page reaches the linear-stage "
                        f"or the equilibrium write")
    law, canon, eq = row["law"], row["canonical"], row.get("equilibrium") or {}
    facts = {"word": row["word"], "P_exact": _num(law["P_exact"]), "P_law": _num(law["P_law"]),
             "difference": f"{law['difference']:+.1e}", "photons": _num(law["n_end"]),
             "kappa": _num(canon["kappa"]), "nbar": _num(canon["nbar"]), "two_d": _num(canon["two_d"]),
             "r": _num(canon["r"]), "h": _num(canon["h"]), "threshold": _num(canon["threshold"]),
             "truncation": f"{law['truncation']:.1e}",
             "P_eq": _num(eq["P_eq"]) if eq.get("P_eq") is not None else "—",
             "gap": _num(eq["gap_near_threshold"]) if eq.get("gap_near_threshold") is not None else "—"}
    tr = row["trace"]
    return {"id": node["id"], "family": "open", "name": node["name_html"], "field": spec.get("field", ""),
            "class": cls, "facts": facts, "question": spec.get("question"),
            "assumptions": spec.get("assumptions", []),
            "engine": {"kind": "open", "t": tr["t"], "x_mean": tr["x_mean"], "x_std": tr["x_std"], "n": tr["n"],
                       "T_ramp": law["T_ramp"], "t_threshold": canon.get("t_threshold"),
                       "P_exact": law["P_exact"], "P_law": law["P_law"]},
            "slots": open_slots(spec, real)}


def open_slots(spec: Dict, real) -> Dict[str, str]:
    pr = real.protocol
    params = real.params
    return {"Omega": "the Lindbladian: the Kerr term, the two-photon drive, the one-photon bias and the loss",
            "Xi": f"one bosonic mode, at most {spec['carrier']['max_quanta']} quanta",
            "C": "the frame rotating at half the pump frequency; the truncation of the Fock space",
            "R": f"the sign of the quadrature ({spec['observable']['sign_of']})/2 at the end",
            "P": f"{pr['sweep_parameter']} from {_num(pr['from'])} to {_num(pr['to'])} at the rate {_num(pr['rate'])}, "
                 f"held for {_num(pr['hold'])}, with the bias {pr['bias_parameter']} = {_num(pr['bias'])}",
            "A": ", ".join(f"{k} = {_num(v)}" for k, v in list(params.items())[:6])}


def regulation_slots(spec: Dict, real, step: Dict) -> Dict[str, str]:
    reg = spec["regulation"]
    params = [k for k in real.params if k != real.input]
    return {"Omega": f"the drift of {', '.join(real.variables)} ({spec.get('kind', 'equations')})",
            "Xi": f"{len(real.variables)} variables on the {spec['carrier'].get('kind', 'euclid')} carrier",
            "C": spec.get("closure") or "; ".join(spec.get("assumptions", [])[:1]) or "as in the source",
            "R": f"the output {reg['output']}" if len(reg["output"]) <= 40 else (spec.get("observable") or "the output"),
            "P": f"a step of {real.input} from {_num(real.u0)} to {_num(step['u1'])}",
            "A": ", ".join(f"{k} = {_num(real.params[k])}" for k in params[:6]) + (" …" if len(params) > 6 else "")}


def stochastic_record(node: Dict) -> Dict:
    from .verification import verify_construction
    spec = node["spec"]
    report = verify_construction(spec)
    growth = report["target"]["drift"]
    return {"id": node["id"], "family": "stochastic", "name": node["name_html"] if node["derived"] else
            ("log of a geometric Brownian motion, " + spec["convention"].capitalize() + " convention"),
            "field": "stochastic processes", "class": "convention",
            "facts": {"convention": spec["convention"].capitalize(), "growth": sympy_html(growth),
                      "correction": sympy_html(report["source_convention_correction"])},
            "question": spec.get("question"), "assumptions": spec.get("assumptions", []),
            "engine": {"kind": "stochastic", "map": "log", "convention": spec["convention"],
                       "params": {"mu": 0.3, "sigma": 0.6}, "ranges": {"mu": [-1, 1, 0.01], "sigma": [0, 1.5, 0.01]},
                       "report": {"target_drift": growth, "variance": report["target"]["variance"]}},
            "slots": {"Omega": f"dX = {spec['drift']} dt + {spec['noise']} dW", "Xi": "X > 0; Y = log X",
                      "C": spec["convention"].capitalize() + " convention for the noise term",
                      "R": "the mean of log X", "P": "X(0) > 0", "A": "μ, σ"}}


def sympy_html(text: str) -> str:
    t = str(text)
    for k, v in GREEK.items():
        t = re.sub(rf"\b{k}\b", v, t)
    return t.replace("**2", "²").replace("*", "·").replace("-", "−")


# ------------------------------------------------------------------------------------------------ edges
def _numeric_drift(spec: Dict) -> Callable[[np.ndarray, Dict[str, float]], np.ndarray]:
    from .memory import spec as mspec
    if spec.get("schema") == "fieldbridge-regulation/1":      # the body of a regulated realization
        spec = {**{k: v for k, v in spec.items() if k != "regulation"}, "schema": mspec.SCHEMA}
    real = mspec.load(spec)
    return lambda q, p: real.drift(np.atleast_2d(q), {**real.params, **p})


def _same_carrier(a: Dict, b: Dict) -> bool:
    if a["family"] != b["family"]:
        return False
    if a["family"] == "unitary":
        ca, cb = a["spec"]["carrier"], b["spec"]["carrier"]
        from .quantum.carriers import make
        A, B = make(ca), make(cb)
        return A.kind == B.kind and A.dim == B.dim and json.dumps(a["spec"].get("sector"), sort_keys=True) == \
            json.dumps(b["spec"].get("sector"), sort_keys=True)
    if a["family"] == "hysterons":
        ha, hb = a["spec"]["hysterons"], b["spec"]["hysterons"]
        return ha.get("lattice") == hb.get("lattice") and ha.get("units") == hb.get("units")
    if a["family"] in ("dissipative", "regulation"):
        ca, cb = a["spec"].get("carrier"), b["spec"].get("carrier")
        if a["spec"].get("kind") == "network" or b["spec"].get("kind") == "network":
            return a["spec"].get("network") == b["spec"].get("network")
        return ca["kind"] == cb["kind"] and ca["variables"] == cb["variables"] and \
            str(ca.get("period", "2*pi")) == str(cb.get("period", "2*pi"))
    return True


def _operators(spec: Dict) -> List[str]:
    return sorted(f"{t['operator']}{'+hc' if t.get('hc') else ''}" for t in spec["hamiltonian"])


def _observable_matrix(spec: Dict) -> np.ndarray:
    from .quantum import language as ql
    real = ql.load(spec)
    return real.O


def _same_up_to_scale(A: np.ndarray, B: np.ndarray) -> bool:
    k = np.vdot(A.ravel(), B.ravel()).real / max(np.vdot(A.ravel(), A.ravel()).real, 1e-300)
    return k > 0 and np.allclose(k * A, B, atol=1e-9 * max(1.0, float(np.abs(B).max())))


def _points(spec: Dict, n: int = 24) -> np.ndarray:
    c = spec["carrier"]
    rng = np.random.default_rng(7)
    dim = len(c["variables"])
    if c["kind"] == "orthant":
        return rng.uniform(0.05, 2.0 * float(c.get("scale", 1.0)), (n, dim))
    if c["kind"] == "torus":
        return rng.uniform(0.0, 2 * np.pi, (n, dim))
    return rng.uniform(-1.5, 1.5, (n, dim))


def classical_limit(q_spec: Dict, d_spec: Dict, j: Optional[float] = None) -> float:
    """The largest deviation, over the magnetization angle phi, of the coherent-state energy of the quantum spin from
    the classical energy of the dissipative realization, after the known finite-spin correction:

        <H>/(2|D| j^2) = e(phi) - 1/2 - cos(phi)^2 / (4j),   m = -<J>/j = (cos phi, 0, sin phi),

    with h = |b| / (2 |D| j) and the field at the angle psi to the easy axis z. For D = 0 the energy is compared
    with e(phi) = -h sin(phi + psi) up to its scale |b| j / h."""
    from fractions import Fraction
    from .quantum.carriers import spin_matrices
    from .memory import spec as mspec
    jj = float(Fraction(str(q_spec["carrier"]["j"]))) if j is None else j
    S = spin_matrices(jj)
    params = {k: float(v) for k, v in (q_spec.get("parameters") or {}).items()}
    from .quantum.language import _coefficient
    coef = {"Jx": 0.0, "Jz": 0.0, "Jz Jz": 0.0}
    for t in q_spec["hamiltonian"]:
        if t["operator"] not in coef:
            raise SiteError(f"the classical limit is defined here for Jx, Jz and Jz Jz, not {t['operator']}")
        coef[t["operator"]] += _coefficient(t["coefficient"], params)
    bx, bz, D = coef["Jx"], coef["Jz"], coef["Jz Jz"]
    real = mspec.load(d_spec)
    h, psi = float(real.params["h"]), float(real.params["psi"])
    field = np.hypot(bx, bz)
    if abs(atan2(bx, bz) - psi) > 1e-9:
        raise SiteError(f"the field of the spin is at {degrees(atan2(bx, bz)):.6g} degrees to the easy axis, "
                        f"the realization has psi = {degrees(psi):.6g}")
    if D != 0 and (D > 0 or abs(field / (2 * abs(D) * jj) - h) > 1e-9):
        raise SiteError("the reduced field h = |b| / (2|D| j) of the spin differs from the realization, or D > 0")
    H = D * S["Jz"] @ S["Jz"] + bx * S["Jx"] + bz * S["Jz"]
    worst = 0.0
    for phi in np.linspace(0.0, 2 * np.pi, 73):
        m = np.array([np.cos(phi), 0.0, np.sin(phi)])
        n = -m
        w, U = np.linalg.eigh(n[0] * S["Jx"] + n[2] * S["Jz"])
        psi_state = U[:, -1]
        E = float(np.real(psi_state.conj() @ H @ psi_state))
        e = float(real.V(np.array([phi]))) if real.potential is not None else float("nan")
        if D != 0:
            dev = E / (2 * abs(D) * jj ** 2) - (e - 0.5 - np.cos(phi) ** 2 / (4 * jj))
        else:
            dev = E / (field * jj) - e / h
        worst = max(worst, abs(dev))
    return worst


def check_edit(edge: Dict, nodes: Dict[str, Dict], records: Dict[str, Dict]) -> str:
    """Raise SiteError unless the edge changes only the component it names; return how it was checked."""
    a, b = nodes[edge["from"]], nodes[edge["to"]]
    kind, slot = edge["kind"], edge["slot"]
    expected = {"param": "A", "term": "Omega", "observable": "R", "attach": "Xi", "codiscovery": "Xi",
                "closure": "C", "protocol": "P"}[kind]
    if slot != expected:
        raise SiteError(f"{edge['id']}: an edge of kind {kind} changes {expected}, not {slot}")
    fail = lambda why: SiteError(f"{edge['id']}: {why}")  # noqa: E731
    if kind == "codiscovery":
        ra, rb = records[edge["from"]], records[edge["to"]]  # only a co-discovery needs the computed mechanisms
        if _same_carrier(a, b):
            raise fail("a co-discovery joins two different carriers")
        target = edge.get("target")
        if target == "rotation":
            if not (str(ra["facts"].get("status", "")).startswith("reached")
                    and str(rb["facts"].get("status", "")).startswith("reached")):
                raise fail("the rotation is not reached in both realizations")
            return "both realizations reach the Bloch rotation (derive_bloch_rotation)"
        if target:
            p = PREFIX[target]
            if not (str(ra["facts"].get(f"{p}_status", "")).startswith("reached")
                    and str(rb["facts"].get(f"{p}_status", "")).startswith("reached")):
                raise fail(f"the target {target} is not reached in both realizations")
            return f"both realizations reach the {target} (codiscovery.{target})"
        if ra["class"] != rb["class"]:
            raise fail(f"the classes differ: {ra['class']} and {rb['class']}")
        return f"both realizations are of the class {ra['class']}"
    if kind == "attach":
        if a["family"] == "hysterons" and b["family"] == "hysterons":
            ha, hb = a["spec"]["hysterons"], b["spec"]["hysterons"]
            diff = {k for k in set(ha) | set(hb) if ha.get(k) != hb.get(k)}
            if diff != {"lattice"}:
                raise fail(f"a change of the lattice changes the lattice only, not {sorted(diff)}")
            return "the hysterons differ only in the lattice"
        if a["family"] == "field":
            fa, fb = a["spec"]["field"], b["spec"]["field"]
            diff = {k for k in set(fa) | set(fb) if fa.get(k) != fb.get(k)}
            if not diff or not diff <= {"d", "L"}:
                raise fail(f"a change of the lattice changes d and L only, not {sorted(diff)}")
            return "the fields differ only in the dimension and the size of the lattice"
        rep = b.get("attach")
        if not rep or rep["source"] != a["spec"]["name"]:
            raise fail("the target was not attached from the source")
        for k, v in rep["preserved"].items():
            if v["target"] is None or abs(v["source"] - v["target"]) > 1e-8 * max(1.0, abs(v["source"])):
                raise fail(f"the attachment does not keep {k}")
        return "the detached rotation is kept: rate, weight and angle (quantum.language.attach)"
    if kind == "closure":
        if a["family"] == "unitary" and b["family"] == "dissipative":
            dev = classical_limit(a["spec"], b["spec"])
            if dev > 1e-9:
                raise fail(f"the classical energy differs from the coherent-state energy by {dev:.3g}")
            return ("the coherent-state energy of the spin equals the classical energy of the realization with the "
                    "finite-spin correction cos²φ/(4j) (site_data.classical_limit)")
        if a["family"] == "field":
            fa, fb = a["spec"]["field"], b["spec"]["field"]
            diff = {k for k in set(fa) | set(fb) if fa.get(k) != fb.get(k)}
            if "conserved" not in diff or not diff <= {"conserved", "kappa0", "Dgrad"}:
                raise fail(f"a change of the conservation law changes conserved (and the non-conserved mass), not "
                           f"{sorted(diff)}")
            return "the fields differ only in the conservation law"
        if a["family"] == "stochastic":
            keys = {k for k in set(a["spec"]) | set(b["spec"]) if a["spec"].get(k) != b["spec"].get(k)}
            if not keys <= {"convention", "question", "assumptions", "provenance"} or "convention" not in keys:
                raise fail(f"a change of convention changes only the convention, not {sorted(keys)}")
            return "the specifications differ only in the convention"
        raise fail("no check for this change of closure")
    if a["family"] != b["family"] or not _same_carrier(a, b):
        raise fail("the carrier changes; name the edge Xi")
    if a["family"] == "open":
        if kind != "param":
            raise fail(f"no check for {kind} on an open realization")
        strip = lambda sp: {k: v for k, v in sp.items() if k not in ("parameters", "carrier", "name", "question",  # noqa: E731
                                                                     "assumptions", "provenance")}
        if strip(a["spec"]) != strip(b["spec"]):
            raise fail("a change of parameters keeps the generator, the protocol and the observable")
        return "the specifications differ only in their parameters and truncation"
    if a["family"] == "regulation":
        ra_, rb_ = a["spec"]["regulation"], b["spec"]["regulation"]
        for key in ("input", "output", "steps"):
            if ra_.get(key) != rb_.get(key):
                raise fail(f"a change of {slot} keeps the {key} of the regulation block")
        if kind not in ("param", "term"):
            raise fail(f"no check for {kind} on a regulated realization")
    if a["family"] == "hysterons":
        ha, hb = a["spec"]["hysterons"], b["spec"]["hysterons"]
        diff = {k for k in set(ha) | set(hb) if ha.get(k) != hb.get(k)}
        allowed = {"term": {"couplings", "directed"}, "protocol": {"drive"},
                   "param": {"fields", "half_widths", "seed"}}.get(kind, set())
        if diff and diff <= allowed:
            return f"the hysterons differ only in {', '.join(sorted(diff))}"
        raise fail(f"the hysterons differ in {sorted(diff)}")
    if a["family"] == "field":
        fa, fb = a["spec"]["field"], b["spec"]["field"]
        diff = {k for k in set(fa) | set(fb) if fa.get(k) != fb.get(k)}
        if kind == "protocol" and diff == {"write"}:
            return "the fields differ only in the shape of the write"
        raise fail(f"the fields differ in {sorted(diff)}")
    if a["family"] == "unitary":
        sa, sb = a["spec"], b["spec"]
        same_ops = _operators(sa) == _operators(sb)
        same_obs = _same_up_to_scale(_observable_matrix(sa), _observable_matrix(sb))
        if kind == "param":
            if not (same_ops and same_obs):
                raise fail("a parameter change keeps the operators and the observable")
            return "same operators and observable; coefficients differ"
        if kind == "term":
            if same_ops or not same_obs:
                raise fail("a change of Omega changes the operators and keeps the observable")
            return "same carrier and observable; the operators of H differ"
        if kind == "observable":
            ha = [(str(t["coefficient"]), t["operator"], bool(t.get("hc"))) for t in sa["hamiltonian"]]
            hb = [(str(t["coefficient"]), t["operator"], bool(t.get("hc"))) for t in sb["hamiltonian"]]
            if ha != hb or sa.get("parameters") != sb.get("parameters") or same_obs:
                raise fail("a change of R keeps the Hamiltonian and changes the observable")
            return "same Hamiltonian; the observable differs"
        raise fail(f"no check for {kind} on a unitary realization")
    # dissipative
    sa, sb = a["spec"], b["spec"]
    Fa, Fb = _numeric_drift(sa), _numeric_drift(sb)
    pa, pb = dict(sa.get("parameters", {})), dict(sb.get("parameters", {}))
    q = _points(sa)
    if kind == "param":
        if pa == pb:
            raise fail("no parameter changes")
        new = set(pb) - set(pa)
        if new:
            # a parameter that the source does not have: at the value named by the edge the target is the source
            reduces = edge.get("reduces") or {}
            if set(reduces) != new:
                raise fail(f"name the values of {sorted(new)} at which the target reduces to the source")
            if not np.allclose(Fb(q, {**pb, **{k: v for k, v in pa.items() if k in pb}, **reduces}), Fa(q, pa),
                               rtol=1e-9, atol=1e-12):
                raise fail(f"the target with {reduces} is not the source")
            return f"the target with {', '.join(f'{k} = {v:g}' for k, v in reduces.items())} is the source"
        mapped = {**pa, **{k: v for k, v in pb.items() if k in pa}}
        if not np.allclose(Fa(q, mapped), Fb(q, pb), rtol=1e-9, atol=1e-12):
            raise fail("the drift of the source with the parameters of the target is not the drift of the target")
        return "the drift of the source with the parameters of the target equals the drift of the target"
    if kind == "term":
        shared = set(pa) & set(pb)
        if any(pa[k] != pb[k] for k in shared):
            raise fail("a change of Omega keeps the shared parameters")
        if np.allclose(Fa(q, pa), Fb(q, pb), rtol=1e-9, atol=1e-12):
            raise fail("the drift does not change")
        return "same carrier and shared parameters; the drift differs"
    if kind == "protocol":
        if not np.allclose(Fa(q, pa), Fb(q, pb), rtol=1e-9, atol=1e-12):
            raise fail("a change of protocol keeps the dynamics")
        if (sa.get("control") or {}).get("name") == (sb.get("control") or {}).get("name"):
            raise fail("the protocol acts on the same parameter")
        return "the same drift; the protocol modulates another parameter"
    raise fail(f"no check for {kind} on a dissipative realization")


# ------------------------------------------------------------------------------------------------ law record
def read_law_record(path: Path = LAW_RECORD) -> Dict:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def write_law_record(data: Dict, path: Path = LAW_RECORD) -> Path:
    nodes = {}
    for nid, rec in data["nodes"].items():
        if rec["family"] != "dissipative":
            continue
        laws = {t: {**r["law"], "word": r.get("word"), "class": r.get("class")}
                for t, r in rec.get("rows", {}).items() if r.get("law")}
        if laws:
            nodes[nid] = {"spec_sha256": rec["spec_sha256"], **laws}
    out = {"schema": "fieldbridge-site-laws/1", "implementation_sha256": _memory_hash(),
           "versions": _versions(), "nodes": nodes}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _versions() -> Dict[str, Optional[str]]:
    v = {"python": platform.python_version()}
    for mod in ("numpy", "scipy", "sympy"):
        try:
            v[mod] = __import__(mod).__version__
        except ImportError:  # pragma: no cover
            v[mod] = None
    return v


def _merge_laws(rec: Dict, record: Dict) -> Optional[str]:
    """Law constants from the record, if the specification is unchanged. Returns the status of the merge."""
    entry = (record.get("nodes") or {}).get(rec["id"])
    if not entry:
        return None
    if entry.get("spec_sha256") != rec["spec_sha256"]:
        return "specification changed since the record"
    for target, law in entry.items():
        if target not in PREFIX or target not in rec["rows"]:
            continue
        row = rec["rows"][target]
        row["law"] = {k: law[k] for k in ("constant", "stderr", "expected", "signature") if k in law}
        if law.get("word"):
            row["word"], row["class"] = law["word"], law.get("class", row.get("class"))
            rec["facts"][f"{PREFIX[target]}_word"] = law["word"]
        p = PREFIX[target]
        rec["facts"][f"{p}_law"], rec["facts"][f"{p}_law_err"] = law["constant"], law["stderr"]
        if p == "lock":
            rec["facts"]["lock_width"], rec["facts"]["lock_width_err"] = law["constant"], law["stderr"]
    return "from the record"


# ------------------------------------------------------------------------------------------------ build
# the mechanisms of the memory materials, in the order of the section on memory: the writes, then no memory, then phase
MEMORY_ORDER = ["symmetric-write", "threshold-write", "subcritical-write", "field-write", "return-point", "no-return",
                "single-state", "oscillation", "neutral-cycles"]
# the subloops of a hysteron node: those of ``memory hysterons`` with its default seed and count
HYSTERON_SEED, HYSTERON_SUBLOOPS = 20260923, 24


def _targets_for(node: Dict) -> List[str]:
    path = node["path"]
    # the oscillators of Module 11 and the ring of three repressors: the realizations in which the tutorial derives
    # phase locking
    oscill = "oscillators" in path or node["id"] == "repressilator"
    if node["family"] != "dissipative":
        return []
    if oscill:
        return ["phase-locking"]
    return ["symmetric-write", "threshold-write"]


def build(root: Path = ROOT, law: bool = False, only: Optional[Iterable[str]] = None, strict: bool = True,
          derive: bool = True, log: Callable[..., None] = print) -> Dict:
    """The data of the page. derive=False skips the co-discovery derivations of the memory realizations (for fast
    tests of the page); the texts that name their results are then left with a dash, which strict=True refuses."""
    nodes = resolve(root, only)
    records: Dict[str, Dict] = {}
    record_file = read_law_record() if not law else {}
    order = sorted(nodes, key=lambda n: (nodes[n]["def"].get("parent") is not None, n))
    for nid in order:
        node = nodes[nid]
        fam = node["family"]
        log(f"  {nid} ({fam})")
        if fam == "unitary":
            pid = node["def"].get("parent")
            parent = None
            if pid:
                if pid not in records:
                    records[pid] = unitary_record(resolve(root, [pid])[pid])
                parent = records[pid]
            records[nid] = unitary_record(node, parent)
        elif fam == "dissipative":
            from .memory import spec as mspec
            real = mspec.load(node["spec"])
            rows = {t: _derive(t, real, nid, law) for t in _targets_for(node)} if derive else {}
            records[nid] = dissipative_record(node, rows)
        elif fam == "field":
            records[nid] = field_record(node)
        elif fam == "hysterons":
            records[nid] = hysteron_record(node)
        elif fam == "regulation":
            records[nid] = regulation_record(node)
        elif fam == "open":
            records[nid] = open_record(node)
        else:
            records[nid] = stochastic_record(node)
        rec = records[nid]
        rec["spec_sha256"] = spec_hash(node["spec"])
        rec["spec"] = node["path"] if not node["derived"] else None
        rec["base_spec"] = node["path"]
        rec["tutorial"] = node["def"].get("tutorial")
        rec["source"] = (node["spec"].get("provenance") or {}).get("source") if isinstance(node["spec"].get("provenance"), dict) else None
        rec["derived"] = node["derived"]
        rec["short"] = reg.SHORT.get(nid) or _short_name(rec["name"])
        rec["auto"] = bool(node["def"].get("auto"))
        rec["universal"] = bool(node["def"].get("universal"))
        if rec["family"] == "dissipative" and not law:
            rec["law_status"] = _merge_laws(rec, record_file)
    records = {k: v for k, v in records.items() if k in nodes}
    edges, checks = [], {}

    def check(e: Dict) -> str:
        # a co-discovery of a memory target is checked on the derivations, which a fast build (derive=False) skips
        if not derive and e["kind"] == "codiscovery" and e.get("target") in PREFIX:
            return "not checked: this build does not derive the targets"
        return check_edit(e, nodes, records)

    for e in reg.EDGES:
        if e["from"] not in nodes or e["to"] not in nodes:
            continue
        checks[e["id"]] = check(e)
        edges.append({"id": e["id"], "from": e["from"], "to": e["to"], "slot": e["slot"], "kind": e["kind"],
                      "target": e.get("target"), "change": e["change"], "reduces": e.get("reduces"),
                      "text": fill(e["text"], records[e["to"]]["facts"], records[e["from"]]["facts"], strict),
                      "check": checks[e["id"]]})
    for e in auto_edges(nodes, records):
        try:
            checks[e["id"]] = check_edit(e, nodes, records)
        except SiteError:
            continue
        edges.append({**e, "text": fill(e["text"], records[e["to"]]["facts"], records[e["from"]]["facts"], strict),
                      "reduces": None, "check": checks[e["id"]]})
    by_id = {e["id"]: e for e in edges}
    sequences = []
    for s in reg.SEQUENCES:
        steps, ok = [], True
        for st in s["steps"]:
            if "node" in st:
                if st["node"] not in records:
                    ok = False
                    break
                steps.append({"kind": "start", "node": st["node"],
                              "text": fill(st["text"], records[st["node"]]["facts"], strict=strict)})
            elif "edge" in st:
                e = by_id.get(st["edge"])
                if e is None:
                    ok = False
                    break
                at = e["from"] if st.get("reverse") else e["to"]
                text = fill(st["text"], records[at]["facts"], strict=strict) if "text" in st else e["text"]
                steps.append({"kind": "edge", "edge": e["id"], "reverse": bool(st.get("reverse")), "text": text})
            else:
                steps.append({"kind": "prepare", "preparation": st["prepare"], "text": st["text"]})
        if ok:
            sequences.append({"id": s["id"], "title": s["title"], "steps": steps})
    for rec in records.values():
        rec.pop("_frame", None)
    present = [c for c in reg.CLASSES if any(r["class"] == c for r in records.values())]
    # the mechanism of each class present, opened on its canonical realization or, without one, on its first node
    mechanisms = {}
    for c in present:
        m = dict(reg.MECHANISMS.get(c, {}), **refs.reading(c))  # empty for a class not yet in site_references
        if m.get("node") not in records or records[m["node"]]["class"] != c:
            m["node"] = sorted(r["id"] for r in records.values() if r["class"] == c)[0]
        mechanisms[c] = m
    start = reg.START if reg.START in records else sorted(records)[0]
    # memory in model materials: every specification file of examples/memory (not the fields), each a material asked
    # whether and how it stores a state; ordered by the mechanism that writes it, the writes first
    order = MEMORY_ORDER + [c for c in present if c not in MEMORY_ORDER]
    materials = sorted((r["id"] for r in records.values()
                        if r["family"] in ("dissipative", "hysterons") and (r.get("spec") or "").startswith("examples/memory/")),
                       key=lambda i: (order.index(records[i]["class"]), _strip(records[i]["name"]).lower()))
    atlas = {"columns": [{"class": c, "label": reg.CLASSES[c],
                          "nodes": sorted((r["id"] for r in records.values() if r["class"] == c),
                                          key=lambda i: (records[i]["family"], records[i]["field"], i))}
                         for c in present]}
    from .memory.cli import BOUNDARY
    law_info = {"computed": bool(law)}
    if not law and record_file:
        law_info.update(record_implementation_sha256=record_file.get("implementation_sha256"),
                        current_implementation_sha256=_memory_hash(), record_versions=record_file.get("versions"))
    return {"schema": SCHEMA, "slots": reg.SLOTS, "classes": reg.CLASSES, "classes_short": reg.CLASS_SHORT,
            "start": start, "mechanisms": mechanisms, "materials": materials, "absent": reg.CLASS_ABSENT,
            "retention_laws": {"link": reg.M4_RETENTION, "laws": reg.RETENTION_LAWS,
                               "reading": {k[-1]: refs.reading(k) for k in refs.RETENTION}},
            "planned": [dict(p, **refs.reading(p["id"])) for p in reg.PLANNED],
            "law_reading": {k[len("law-"):]: refs.reading(k) for k in refs.CERTIFIED}, "mechanisms_doc": refs.DOC,
            "nodes": records, "edges": edges, "sequences": sequences, "atlas": atlas,
            "codiscovery": codiscovery_summary(records), "boundary": BOUNDARY, "law": law_info,
            "provenance": {"memory_implementation_sha256": _memory_hash(), "versions": _versions()}}


AUTO_TEXT = {"symmetric-write": "In the {name} the constructor reaches the symmetric write as well (derivation {sym_word}).",
             "threshold-write": "In the {name} the constructor reaches the one-sided write as well (derivation {thr_word}).",
             "phase-locking": "The {name} locks to a drive as well, at {lock_ratio}:1 (derivation {lock_word}).",
             "rotation": "In the {name} the constructor reaches the Bloch rotation as well (derivation {word})."}


def auto_edges(nodes: Dict[str, Dict], records: Dict[str, Dict]) -> List[Dict]:
    """For every automatic node, one edge from the first registered node that reaches the same target (or, without
    a target, has the same class) on another carrier."""
    out = []
    registered = [n for n in nodes if not nodes[n]["def"].get("auto")]
    for nid, node in nodes.items():
        if not node["def"].get("auto"):
            continue
        rec = records[nid]
        targets = [t for t in TARGETS if str(rec["facts"].get(f"{PREFIX[t]}_status", "")).startswith("reached")]
        if rec["family"] == "unitary" and rec["class"] == "rotation":
            targets = ["rotation"]
        for target in targets + [None]:
            for rid in registered:
                r = records[rid]
                same = (str(r["facts"].get("status", "")).startswith("reached") if target == "rotation" else
                        str(r["facts"].get(f"{PREFIX[target]}_status", "")).startswith("reached") if target else
                        r["class"] == rec["class"])
                if same and not _same_carrier(nodes[rid], node):
                    text = AUTO_TEXT.get(target, "The calculation places the {name} in the class "
                                         + reg.CLASSES.get(rec["class"], rec["class"])
                                         + "; the specification has no entry in the registry of the page.")
                    out.append({"id": f"auto_{nid}", "from": rid, "to": nid, "slot": "Xi", "kind": "codiscovery",
                                "target": target, "change": f"{r['short']} → {rec['short']}",
                                "text": text.replace("{name}", _strip(rec["name"]))})
                    break
            if any(e["to"] == nid for e in out):
                break
    return out


def codiscovery_summary(records: Dict[str, Dict]) -> Dict:
    """The data of the collapse plots: every realization that reaches a target, in the canonical units of the
    target."""
    rot = [{"id": r["id"], "name": r["name"], "field": r["field"], "rate": r["facts"]["rate"],
            "theta": r["facts"]["theta"], "residual": r["facts"]["residual"], **r["law_curve"]}
           for r in records.values() if r["family"] == "unitary" and r["class"] == "rotation" and r.get("law_curve")]
    out = {"rotation": rot}
    for target in TARGETS:
        rows = []
        for r in records.values():
            row = (r.get("rows") or {}).get(target)
            if row and str(row.get("status", "")).startswith("reached"):
                rows.append({"id": r["id"], "name": r["name"], "field": r["field"], "word": row.get("word"),
                             **({"canonical": row["canonical"]} if row.get("canonical") else {}),
                             **({"law": {k: v for k, v in row["law"].items() if k != "signature"}}
                                if row.get("law") else {})})
        # a realization that simulates the equations of another (the unfolded unequal toggle and the toggle) is a
        # replicate of that model: the page draws it open and names the model
        from .memory.construct import same_equations
        first = same_equations([((records[x["id"]].get("rows") or {}).get(target) or {}).get("law", {}).get("signature")
                                for x in rows])
        for x, f in zip(rows, first):
            if f != rows.index(x):
                x["same_as"] = rows[f]["id"]
        out[target] = rows
    return out


def write(data: Dict, out: Path) -> Path:
    """The data as a script, so that the page also works when opened from the file system."""
    out.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(data, ensure_ascii=False, separators=(",", ":"), default=float).replace("</", "<\\/")
    out.write_text("window.FIELDBRIDGE_SITE = " + text + ";\n", encoding="utf-8")
    return out
