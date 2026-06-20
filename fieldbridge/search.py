from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from .database import field_by_id, load_all
from .models import AnalogyMatch, Fingerprint, MechanismRecord, Translation
from .routes import fingerprint_text


def cosine(left: List[float], right: List[float]) -> float:
    dot = sum(a * b for a, b in zip(left, right))
    left_norm = math.sqrt(sum(a * a for a in left))
    right_norm = math.sqrt(sum(b * b for b in right))
    if left_norm == 0 or right_norm == 0:
        return 0.0
    return dot / (left_norm * right_norm)


def tokens(text: str) -> List[str]:
    return re.findall(r"[a-zA-Z][a-zA-Z0-9_-]{2,}", text.lower())


def overlap(left: Dict[str, float], right: Dict[str, float]) -> Dict[str, float]:
    return {
        key: min(value, right.get(key, 0.0))
        for key, value in left.items()
        if min(value, right.get(key, 0.0)) > 0.05
    }


def keyword_hits(query: str, record: MechanismRecord) -> List[str]:
    query_tokens = set(tokens(query))
    return sorted({word for word in record.keywords if word.lower() in query_tokens})


def score_record(query: str, query_fp: Fingerprint, record: MechanismRecord) -> AnalogyMatch:
    record_fp = record.fingerprint()
    route_overlap = overlap(query_fp.routes, record_fp.routes)
    fiber_overlap = overlap(query_fp.fibers, record_fp.fibers)
    hits = keyword_hits(query, record)
    vector_score = cosine(query_fp.vector(), record_fp.vector())
    route_score = sum(route_overlap.values()) / max(1, len(route_overlap))
    fiber_score = 0.35 * sum(fiber_overlap.values()) / max(1, len(fiber_overlap))
    keyword_score = min(0.18, 0.035 * len(hits))
    score = max(0.0, min(1.0, 0.62 * vector_score + 0.2 * route_score + fiber_score + keyword_score))
    return AnalogyMatch(
        record=record,
        score=round(score, 4),
        route_overlap=route_overlap,
        fiber_overlap=fiber_overlap,
        keyword_hits=hits,
    )


def is_hyperion_record(record: MechanismRecord) -> bool:
    return record.field_id == "hyperion_equation" or record.source == "hyperion_equation_witnesses"


def record_matches_field(record: MechanismRecord, target_field: Optional[str]) -> bool:
    if target_field is None:
        return True
    return record.field_id == target_field or target_field in (record.target_fields or [])


def find_analogs(
    text: str,
    target_field: Optional[str] = None,
    data_dir: Path | None = None,
    top_k: int = 5,
    include_hyperion: bool = False,
) -> List[AnalogyMatch]:
    _, records = load_all(data_dir, include_static_index=include_hyperion)
    query_fp = fingerprint_text(text)
    scored = [
        score_record(text, query_fp, record)
        for record in records
        if record_matches_field(record, target_field)
    ]
    return sorted(scored, key=lambda item: item.score, reverse=True)[:top_k]


def _compose_invariant(matches: Iterable[AnalogyMatch]) -> str:
    for match in matches:
        if match.record.invariant:
            return match.record.invariant
    return (
        "A source operation writes or constrains an internal state; a later operation "
        "reads that state through transport, closure, operator response, and boundary conditions."
    )


def _active_route_names(fp: Fingerprint, threshold: float = 0.22) -> List[str]:
    return [key.replace("_route", "") for key, value in sorted(fp.routes.items(), key=lambda item: item[1], reverse=True) if value >= threshold]


def _first(values: List[str], fallback: str) -> str:
    return next((value for value in values if value), fallback)


def _target_role_words(target) -> tuple[str, str, str, str]:
    if target.field_id == "material_intelligence":
        return (
            "a retained material state",
            "a write protocol or environmental exposure",
            "a material boundary, interface, sample geometry, or device condition",
            "a later transport, conductance, release, shape, force, motion, or actuation response",
        )
    if target.field_id == "biological_intelligence":
        return (
            "a retained cellular, tissue, immune, bioelectric, or developmental state",
            "a morphogen, metabolite, mechanical, electrical, immune, or neural cue",
            "a membrane, niche, wound, organoid, tissue compartment, or assay boundary",
            "a later phenotype, migration, secretion, repair, morphology, or fate response",
        )
    if target.field_id == "collective_intelligence":
        return (
            "a retained trace, memory, belief, record, policy, or environmental mark",
            "a signal, task cue, message, reward, event, or local interaction",
            "a network, platform, institution, task, communication, or incentive boundary",
            "a later route choice, consensus, allocation, cascade, decision, or collective response",
        )
    return (
        _first(target.state_words, "retained state q"),
        _first(target.input_words, "input u"),
        _first(target.boundary_words, "boundary B"),
        _first(target.output_words, "output y"),
    )


def _sentence_after_intro(text: str) -> str:
    text = (text or "").strip()
    if not text:
        return text
    return text[:1].lower() + text[1:]


def _target_formulation(text: str, source_fp: Fingerprint, target, best: MechanismRecord | None) -> str:
    state, input_signal, boundary, output = _target_role_words(target)
    routes = _active_route_names(source_fp)
    route_clause = f" Active source routes: {', '.join(routes)}." if routes else ""
    receptor_clause = ""
    if best:
        receptor_clause = f" Closest concrete receptor in the pack: {best.title}. In that receptor, {_sentence_after_intro(best.summary)}"
    return (
        f"Render the source mechanism in {target.label} as follows: {input_signal} writes or selects "
        f"{state}; after the writing condition is removed, changed, or read under a new condition, "
        f"{state} changes {output} under {boundary}. The transfer is supported only if the named "
        f"state reduces the output residual against erased-state, shuffled-history, and boundary-swap controls."
        f"{route_clause}{receptor_clause}"
    )


def translate_mechanism(
    text: str,
    target_field: str,
    data_dir: Path | None = None,
    top_k: int = 4,
    include_hyperion: bool = True,
) -> Translation:
    field_packs, _ = load_all(data_dir)
    target = field_by_id(field_packs, target_field)
    native_matches = find_analogs(
        text,
        target_field=target_field,
        data_dir=data_dir,
        top_k=max(1, top_k),
        include_hyperion=False,
    )
    evidence_matches = (
        [
            match
            for match in find_analogs(
                text,
                target_field=target_field,
                data_dir=data_dir,
                top_k=max(top_k * 2, top_k),
                include_hyperion=True,
            )
            if is_hyperion_record(match.record)
        ]
        if include_hyperion
        else []
    )
    matches = (native_matches[: max(1, top_k - 2)] + evidence_matches[:top_k])[:top_k]
    source_fp = fingerprint_text(text)
    best = native_matches[0].record if native_matches else None
    if native_matches:
        variables = best.variables
        equations = best.equations
        measurements = best.measurements
        controls = best.controls
    else:
        variables = [
            f"q: {', '.join(target.state_words[:4])}",
            f"u: {', '.join(target.input_words[:4])}",
            f"B: {', '.join(target.boundary_words[:4])}",
            f"y: {', '.join(target.output_words[:4])}",
        ]
        equations = ["partial_t q = W(u,t)-q/tau", "C(q,u,B)=0", "y=A(q,B)"]
        measurements = target.output_words[:4]
        controls = target.validation_words[:4]
    formulation = _target_formulation(text, source_fp, target, best)

    return Translation(
        source_fingerprint=source_fp,
        target_field=target,
        matches=matches,
        invariant=_compose_invariant(matches),
        target_formulation=formulation,
        variables=variables,
        equations=equations,
        measurements=measurements,
        controls=controls,
        evidence_boundary=(
            "FieldBridge proposes a structural analogue. It is evidence-backed only by "
            "the shipped examples, deterministic route/fiber fingerprint, and optional "
            "Hyperion static witness index. Raw arXiv witnesses are audit evidence, not "
            "the translated mechanism."
        ),
    )
