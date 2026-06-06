from __future__ import annotations

import math
import re
from typing import Dict, Iterable, List, Tuple

from .models import Fingerprint


ROUTES: Dict[str, Tuple[str, List[str]]] = {
    "transport_flow_route": (
        "flow, balance, propagation, diffusion, force, current, or evolution",
        [
            r"\\partial_t",
            r"\bpartial_t\b",
            r"\bnabla\b",
            r"\\nabla",
            r"\bdiv\b",
            r"\\cdot",
            r"\bflow\b",
            r"\btransport\b",
            r"\bdiffusion\b",
            r"\bcurrent\b",
            r"\bflux\b",
            r"\bforce\b",
            r"\bmotion\b",
            r"\bevolution\b",
            r"\bpropagation\b",
        ],
    ),
    "constraint_closure_route": (
        "constraint, admissibility, constitutive law, normalization, or residual",
        [
            r"\bconstraint\b",
            r"\bclosure\b",
            r"\bresidual\b",
            r"\badmissib",
            r"\bconstitutive\b",
            r"\bbalance\b",
            r"\bnormalization\b",
            r"\bgauge\b",
            r"\bconserved\b",
            r"\bno-flux\b",
            r"\bforce balance\b",
            r"C\(",
        ],
    ),
    "spectral_operator_route": (
        "operator, mode, kernel, spectrum, eigenvalue, response, or transfer function",
        [
            r"\boperator\b",
            r"\bspectral\b",
            r"\bspectrum\b",
            r"\beigen",
            r"\blambda\b",
            r"\\lambda",
            r"\bkernel\b",
            r"\bmode\b",
            r"\bfrequency\b",
            r"\bresponse curve\b",
            r"\btransfer function\b",
            r"\bHamiltonian\b",
            r"\bLaplacian\b",
        ],
    ),
    "boundary_weak_form_route": (
        "boundary, interface, domain, weak form, surface, membrane, or context",
        [
            r"\bboundary\b",
            r"\binterface\b",
            r"\bdomain\b",
            r"\bsurface\b",
            r"\bmembrane\b",
            r"\bwall\b",
            r"\bweak form\b",
            r"\bcontext\b",
            r"\bcompartment\b",
            r"\bniche\b",
            r"\bgeometry\b",
            r"\\partial_n",
        ],
    ),
    "commutator_incompatibility_route": (
        "non-commuting order, hysteresis, curvature, bracket, or compatibility failure",
        [
            r"\[[^\]]+,[^\]]+\]",
            r"\bcommutator\b",
            r"\bnoncommut",
            r"\bhysteresis\b",
            r"\border-dependent\b",
            r"\bcurvature\b",
            r"\bincompatib",
            r"\bobstruction\b",
            r"\bpath dependent\b",
            r"\bmemory\b",
        ],
    ),
    "discrete_protocol_route": (
        "algorithm, recurrence, update, protocol, feedback, learning, or self-modification",
        [
            r"\bn\+1\b",
            r"_\{?n\+1\}?",
            r"\bt\+1\b",
            r"\bupdate\b",
            r"\biteration\b",
            r"\bprotocol\b",
            r"\bfeedback\b",
            r"\blearning\b",
            r"\btraining\b",
            r"\bcontroller\b",
            r"\bpolicy\b",
            r"\bself-modif",
            r"\balgorithm\b",
        ],
    ),
}

FIBERS: Dict[str, Tuple[str, List[str]]] = {
    "structure": (
        "symbolic architecture: brackets, equations, nesting, dependencies, and graph form",
        [
            r"=",
            r":=",
            r"\\frac",
            r"\\int",
            r"\\sum",
            r"\bdepends\b",
            r"\bgraph\b",
            r"\bstate variable\b",
            r"\balgebra\b",
        ],
    ),
    "spectral": (
        "global operator behavior: modes, spectra, transfer curves, and stability",
        [
            r"\bspectral\b",
            r"\beigen",
            r"\bmode\b",
            r"\bfrequency\b",
            r"\bstability\b",
            r"\bdispersion\b",
            r"\bresponse\b",
            r"\btransfer\b",
        ],
    ),
    "geometry": (
        "realization layer: shape, coordinates, dimensions, domains, surfaces, and boundaries",
        [
            r"\bgeometry\b",
            r"\bcoordinate\b",
            r"\bdimension\b",
            r"\bdomain\b",
            r"\bsurface\b",
            r"\bboundary\b",
            r"\bmanifold\b",
            r"\binterface\b",
        ],
    ),
    "syntax": (
        "local names and topic words that may change across fields",
        [
            r"\bfluid\b",
            r"\bcell\b",
            r"\btissue\b",
            r"\bagent\b",
            r"\bdroplet\b",
            r"\bhydrogel\b",
            r"\bnetwork\b",
            r"\bparticle\b",
        ],
    ),
    "entropy": (
        "uncertainty, distribution, density, probability, information, or free energy",
        [
            r"\bentropy\b",
            r"\bprobability\b",
            r"\bdensity\b",
            r"\bdistribution\b",
            r"\bfree energy\b",
            r"\binformation\b",
            r"\bnoise\b",
            r"\bensemble\b",
        ],
    ),
}


def _compact(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def _hits(text: str, patterns: Iterable[str]) -> List[str]:
    found: List[str] = []
    for pattern in patterns:
        if re.search(pattern, text, flags=re.IGNORECASE):
            found.append(pattern)
    return found


def _score(hit_count: int, bonus: float = 0.0) -> float:
    return max(0.0, min(1.0, 1.0 - math.exp(-(0.34 * hit_count + bonus))))


def _equation_bonus(text: str, route_id: str) -> float:
    bonus = 0.0
    if (
        route_id == "constraint_closure_route"
        and "=" in text
        and re.search(r"\bC\s*\(|constraint|closure|residual|conserved|balance|admissib", text, flags=re.IGNORECASE)
    ):
        bonus += 0.2
    if re.search(r"\\partial_t|partial_t|∂_t", text) and route_id == "transport_flow_route":
        bonus += 0.35
    if re.search(r"\\lambda|lambda|eigen|spectrum|operator", text, flags=re.IGNORECASE) and route_id == "spectral_operator_route":
        bonus += 0.25
    if re.search(r"\\partial_n|boundary|interface|domain", text, flags=re.IGNORECASE) and route_id == "boundary_weak_form_route":
        bonus += 0.25
    if re.search(r"_\{?n\+1\}?|n\+1|update|protocol|feedback", text, flags=re.IGNORECASE) and route_id == "discrete_protocol_route":
        bonus += 0.25
    return bonus


def fingerprint_text(text: str) -> Fingerprint:
    """Compute a transparent route/fiber fingerprint for a text or equation block."""
    clean = _compact(text)
    route_scores: Dict[str, float] = {}
    fiber_scores: Dict[str, float] = {}
    all_hits: Dict[str, List[str]] = {}

    for route_id, (_, patterns) in ROUTES.items():
        hits = _hits(clean, patterns)
        route_scores[route_id] = _score(len(hits), _equation_bonus(clean, route_id))
        all_hits[route_id] = hits

    for fiber_id, (_, patterns) in FIBERS.items():
        hits = _hits(clean, patterns)
        fiber_scores[fiber_id] = _score(len(hits))
        all_hits[fiber_id] = hits

    return Fingerprint(routes=route_scores, fibers=fiber_scores, hits=all_hits)
