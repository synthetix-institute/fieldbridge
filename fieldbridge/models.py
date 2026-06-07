from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass(frozen=True)
class Fingerprint:
    routes: Dict[str, float]
    fibers: Dict[str, float]
    hits: Dict[str, List[str]] = field(default_factory=dict)

    def vector(self) -> List[float]:
        route_order = [
            "transport_flow_route",
            "constraint_closure_route",
            "spectral_operator_route",
            "boundary_weak_form_route",
            "commutator_incompatibility_route",
            "discrete_protocol_route",
        ]
        fiber_order = ["structure", "spectral", "geometry", "syntax", "entropy"]
        return [self.routes.get(key, 0.0) for key in route_order] + [
            self.fibers.get(key, 0.0) for key in fiber_order
        ]


@dataclass(frozen=True)
class FieldPack:
    field_id: str
    label: str
    description: str
    state_words: List[str]
    input_words: List[str]
    boundary_words: List[str]
    output_words: List[str]
    validation_words: List[str]


@dataclass(frozen=True)
class MechanismRecord:
    record_id: str
    title: str
    field_id: str
    summary: str
    invariant: str
    equations: List[str]
    variables: List[str]
    measurements: List[str]
    controls: List[str]
    references: List[str]
    keywords: List[str]
    routes: Dict[str, float]
    fibers: Dict[str, float]
    source: Optional[str] = None
    target_fields: List[str] = field(default_factory=list)
    hyperion: Optional[Dict[str, Any]] = None

    def fingerprint(self) -> Fingerprint:
        return Fingerprint(routes=self.routes, fibers=self.fibers, hits={})


@dataclass(frozen=True)
class AnalogyMatch:
    record: MechanismRecord
    score: float
    route_overlap: Dict[str, float]
    fiber_overlap: Dict[str, float]
    keyword_hits: List[str]


@dataclass(frozen=True)
class MechanismSheet:
    source_title: str
    fingerprint: Fingerprint
    invariant: str
    state: str
    input_signal: str
    boundary: str
    output: str
    equations: List[str]
    variables: List[str]
    measurements: List[str]
    controls: List[str]
    route_explanation: List[str]
    evidence_boundary: str


@dataclass(frozen=True)
class Comparison:
    source: MechanismSheet
    target: MechanismSheet
    score: float
    preserved_routes: Dict[str, float]
    changed_routes: Dict[str, float]
    preserved_fibers: Dict[str, float]
    changed_fibers: Dict[str, float]
    analogous_equations: List[str]
    interpretation: str


@dataclass(frozen=True)
class Translation:
    source_fingerprint: Fingerprint
    target_field: FieldPack
    matches: List[AnalogyMatch]
    invariant: str
    target_formulation: str
    variables: List[str]
    equations: List[str]
    measurements: List[str]
    controls: List[str]
    evidence_boundary: str
