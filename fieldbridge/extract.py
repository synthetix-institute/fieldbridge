from __future__ import annotations

import re
from typing import Dict, List, Tuple

from .models import Comparison, Fingerprint, MechanismSheet
from .routes import FIBERS, ROUTES, fingerprint_text
from .search import cosine


DISPLAY_MATH_PATTERNS = [
    re.compile(r"\\begin\{(?:equation\*?|align\*?|gather\*?|multline\*?)\}([\s\S]*?)\\end\{(?:equation\*?|align\*?|gather\*?|multline\*?)\}"),
    re.compile(r"\\\[([\s\S]*?)\\\]"),
    re.compile(r"\$\$([\s\S]*?)\$\$"),
]


def compact(text: str, limit: int = 500) -> str:
    clean = re.sub(r"\s+", " ", text).strip()
    return clean if len(clean) <= limit else f"{clean[: limit - 3].strip()}..."


def extract_equations(text: str, limit: int = 8) -> List[str]:
    equations: List[str] = []
    seen = set()

    def add(value: str) -> None:
        clean = compact(value.replace("\\\\", "\n"), 420)
        if len(clean) < 3:
            return
        key = clean.lower()
        if key in seen:
            return
        seen.add(key)
        equations.append(clean)

    for pattern in DISPLAY_MATH_PATTERNS:
        for match in pattern.finditer(text):
            add(match.group(1))
            if len(equations) >= limit:
                return equations

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not (4 <= len(line) <= 420):
            continue
        word_count = len(line.split())
        if word_count > 18:
            continue
        if word_count > 8 and re.search(r"[.!?]", line):
            continue
        has_operator = bool(re.search(r"=|\\partial|partial_t|∂|\\nabla|nabla|\\int|\\sum|_\{?n\+1\}?", line))
        has_symbol = bool(re.search(r"[A-Za-z][A-Za-z0-9_{}()\\^]*", line))
        if has_operator and has_symbol:
            add(line)
            if len(equations) >= limit:
                return equations
    return equations


def active_routes(fp: Fingerprint, threshold: float = 0.22) -> List[str]:
    return [route for route, score in sorted(fp.routes.items(), key=lambda item: item[1], reverse=True) if score >= threshold]


def guess_roles(text: str, fp: Fingerprint) -> Tuple[str, str, str, str]:
    lower = text.lower()
    if re.search(r"cell|tissue|wound|morphology|regeneration|immune|neural|bioelectric|organoid", lower):
        state = "q(x,t): retained cellular, tissue, bioelectric, immune, or developmental state"
        input_signal = "u(t): morphogen, metabolite, electrical, mechanical, immune, or neural cue"
        boundary = "B: membrane, wound edge, niche, organoid, tissue compartment, or assay boundary"
        output = "y(t): phenotype, migration, secretion, repair, morphology, or fate response"
    elif re.search(r"hydrogel|material|droplet|interface|polymer|membrane|porous|crystal|colloid", lower):
        state = "m(x,t): retained material state"
        input_signal = "u(t): write protocol, field, load, fuel, light, voltage, or chemical exposure"
        boundary = "B: sample boundary, interface, membrane, pore geometry, wall, or device body"
        output = "y(t): transport, conductance, release, shape, force, motion, or actuation"
    elif re.search(r"agent|collective|swarm|network|institution|consensus|allocation|policy|trace", lower):
        state = "m_t: retained trace, record, belief, policy, or environmental mark"
        input_signal = "u_t: signal, task cue, message, local interaction, reward, or event"
        boundary = "B: network, platform, arena, institution, task, or communication topology"
        output = "a_t: route choice, consensus, allocation, cascade, decision, or group response"
    else:
        state = "q: retained internal state"
        input_signal = "u: input, perturbation, write condition, or operation"
        boundary = "B: boundary, interface, context, domain, or realization condition"
        output = "y: measured output or later response"

    if fp.routes.get("discrete_protocol_route", 0) >= 0.5:
        state = state.replace("q(x,t)", "q_n").replace("q:", "q_n:").replace("m(x,t)", "m_n").replace("m_t", "m_t")
    return state, input_signal, boundary, output


def default_equations(fp: Fingerprint) -> List[str]:
    routes = set(active_routes(fp))
    equations = []
    if "discrete_protocol_route" in routes:
        equations.append("q_{n+1}=P(q_n,u_n,B)")
    else:
        equations.append("partial_t q = W(u,t)-q/tau_q")
    equations.append("C(q,u,B)=0")
    if "transport_flow_route" in routes:
        equations.append("partial_t rho + nabla cdot J(rho,q,B)=S")
    if "boundary_weak_form_route" in routes:
        equations.append("B[q]=b or int_Omega phi Lq = int_Omega phi f")
    if "spectral_operator_route" in routes:
        equations.append("L(q,B) psi = lambda psi")
    if "commutator_incompatibility_route" in routes:
        equations.append("[A,B]q != 0")
    equations.append("y=A(q,B)")
    return equations


def route_explanation(fp: Fingerprint) -> List[str]:
    rows = []
    for route_id in active_routes(fp):
        rows.append(f"{route_id}: {ROUTES[route_id][0]} ({fp.routes[route_id]:.2f})")
    return rows or ["No route crossed the active threshold; extraction is low confidence."]


def extract_mechanism(text: str, title: str = "input") -> MechanismSheet:
    fp = fingerprint_text(text)
    equations = extract_equations(text)
    if not equations:
        equations = default_equations(fp)
    state, input_signal, boundary, output = guess_roles(text, fp)
    active = set(active_routes(fp))
    invariant_parts = ["A mechanism is extracted as state, input, boundary, output, and controls."]
    if {"transport_flow_route", "constraint_closure_route"} & active:
        invariant_parts.append("The central invariant is a state transition constrained by an admissibility or residual law.")
    if "boundary_weak_form_route" in active:
        invariant_parts.append("The boundary is treated as part of the mechanism, not as passive scenery.")
    if "discrete_protocol_route" in active:
        invariant_parts.append("The order of operations is part of the mechanism.")
    if "commutator_incompatibility_route" in active:
        invariant_parts.append("Non-commuting operations or history effects must be tested directly.")

    return MechanismSheet(
        source_title=title,
        fingerprint=fp,
        invariant=" ".join(invariant_parts),
        state=state,
        input_signal=input_signal,
        boundary=boundary,
        output=output,
        equations=equations,
        variables=[state, input_signal, boundary, output],
        measurements=[
            "measure the retained state before and after input removal",
            "measure the later output under the same read condition",
            "fit the residual using state, input, boundary, and output roles",
        ],
        controls=[
            "no-write or no-perturbation baseline",
            "erased or reset state",
            "shuffled write/read order",
            "same final condition with different history",
            "boundary or context swap",
        ],
        route_explanation=route_explanation(fp),
        evidence_boundary=(
            "This sheet is an extracted mechanism hypothesis. It is not a validation; "
            "the named state, equations, and controls must be checked in the source field."
        ),
    )


def _overlap(left: Dict[str, float], right: Dict[str, float]) -> Dict[str, float]:
    return {
        key: min(left_value, right.get(key, 0.0))
        for key, left_value in left.items()
        if min(left_value, right.get(key, 0.0)) > 0.05
    }


def _changed(left: Dict[str, float], right: Dict[str, float]) -> Dict[str, float]:
    rows = {
        key: abs(left.get(key, 0.0) - right.get(key, 0.0))
        for key in sorted(set(left) | set(right))
    }
    return {key: value for key, value in rows.items() if value > 0.15}


def compare_mechanisms(source_text: str, target_text: str, source_title: str = "source", target_title: str = "target") -> Comparison:
    source = extract_mechanism(source_text, source_title)
    target = extract_mechanism(target_text, target_title)
    preserved_routes = _overlap(source.fingerprint.routes, target.fingerprint.routes)
    preserved_fibers = _overlap(source.fingerprint.fibers, target.fingerprint.fibers)
    changed_routes = _changed(source.fingerprint.routes, target.fingerprint.routes)
    changed_fibers = _changed(source.fingerprint.fibers, target.fingerprint.fibers)
    score = cosine(source.fingerprint.vector(), target.fingerprint.vector())
    analogous = []
    for left in source.equations[:4]:
        for right in target.equations[:4]:
            pair_fp_left = fingerprint_text(left)
            pair_fp_right = fingerprint_text(right)
            if cosine(pair_fp_left.vector(), pair_fp_right.vector()) >= 0.62:
                analogous.append(f"{left}  <->  {right}")
                break
    if score >= 0.72:
        interpretation = "The two papers share a strong mechanism-level analogy."
    elif score >= 0.45:
        interpretation = "The two papers share a partial mechanism analogy; inspect changed routes before translating."
    else:
        interpretation = "The analogy is weak at route/fiber level; any translation should be treated as speculative."
    return Comparison(
        source=source,
        target=target,
        score=round(score, 4),
        preserved_routes=preserved_routes,
        changed_routes=changed_routes,
        preserved_fibers=preserved_fibers,
        changed_fibers=changed_fibers,
        analogous_equations=analogous[:6],
        interpretation=interpretation,
    )
