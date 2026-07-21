from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

from .extract import active_routes, extract_mechanism
from .models import ConstructorTransfer
from .routes import ROUTES
from .search import translate_mechanism


def _operator_proxy(routes: List[str]) -> str:
    operator_routes = [
        route
        for route in routes
        if route in {
            "transport_flow_route",
            "spectral_operator_route",
            "commutator_incompatibility_route",
            "discrete_protocol_route",
        }
    ]
    if not operator_routes:
        return "operator apparatus unresolved from the public fingerprint"
    labels = [ROUTES[route][0] for route in operator_routes[:3]]
    return "; ".join(labels)


def _first(values: List[str], fallback: str) -> str:
    return next((value for value in values if value), fallback)


def construct_transfer(
    text: str,
    target_field: str,
    data_dir: Path | None = None,
    top_k: int = 4,
    include_hyperion: bool = True,
) -> ConstructorTransfer:
    """Build a reviewable mechanism-preserving transfer contract.

    The public fingerprint supplies transparent proxies for the learned
    operator/substrate language. It does not assign private Omega or Xi tokens.
    """
    source = extract_mechanism(text, title="source")
    translation = translate_mechanism(
        text,
        target_field=target_field,
        data_dir=data_dir,
        top_k=top_k,
        include_hyperion=include_hyperion,
    )
    routes = active_routes(source.fingerprint)
    operator = _operator_proxy(routes)
    target_state = _first(translation.variables, "target carrier unresolved")
    target_boundary = next(
        (value for value in translation.variables if value.startswith("B:")),
        "target closure or boundary unresolved",
    )
    target_input = next(
        (value for value in translation.variables if value.startswith("u:")),
        next(
            (value for value in translation.variables if value.startswith(("T:", "c:"))),
            "target protocol or input unresolved",
        ),
    )

    source_identity: Dict[str, Any] = {
        "M": {
            "Omega_proxy": operator,
            "Xi_proxy": source.state,
        },
        "F": {
            "C": source.boundary,
            "R": source.output,
            "P": source.input_signal,
        },
        "A": source.source_title,
    }
    target_identity: Dict[str, Any] = {
        "M": {
            "Omega_proxy": operator,
            "Xi_proxy": target_state,
        },
        "F": {
            "C": target_boundary,
            "R": translation.measurements,
            "P": [target_input, *translation.controls],
        },
        "A": translation.target_field.label,
    }

    preserved = [
        f"operator contract: {operator}",
        f"invariant: {translation.invariant}",
    ]
    changed = [
        f"substrate/carrier: {source.state} -> {target_state}",
        f"closure: {source.boundary} -> {target_boundary}",
        f"realization: {source.source_title} -> {translation.target_field.label}",
    ]
    attachments = {
        "C_closure": [target_boundary],
        "R_readout": translation.measurements,
        "P_protocol": [target_input],
        "A_realization": translation.variables,
        "falsifiers": translation.controls,
    }
    moves = [
        {"road": "preserve", "acts_on": "Omega", "operation": "retain the qualified operator contract"},
        {"road": "reattach", "acts_on": "Xi", "operation": "replace the source carrier with a target-field carrier"},
        {"road": "close", "acts_on": "C", "operation": "supply target admissibility, boundary, or normalization"},
        {"road": "observe", "acts_on": "R", "operation": "define a measurable consequence"},
        {"road": "execute", "acts_on": "P", "operation": "define the intervention or computational protocol"},
        {"road": "falsify", "acts_on": "I_op", "operation": "break one retained clause with a negative control"},
    ]
    predictions = [
        "The target equations must preserve the stated invariant within the declared closure.",
        "The proposed readout must respond to the target protocol through the retained operator contract.",
        "At least one clause-breaking control must abolish or materially reduce that response.",
    ]
    structural_gates = {
        "source_equation_present": bool(source.equations),
        "target_receptor_present": bool(translation.matches),
        "candidate_equation_present": bool(translation.equations),
        "closure_present": bool(attachments["C_closure"]),
        "readout_present": bool(attachments["R_readout"]),
        "protocol_present": bool(attachments["P_protocol"]),
        "falsifier_present": bool(attachments["falsifiers"]),
    }
    gates = {
        **structural_gates,
        "field_blind_retrieval": False,
        "independent_target_validation": False,
        "prospective_prediction": False,
    }
    readiness = (
        "structurally_complete_constructor_proposal"
        if all(structural_gates.values())
        else "incomplete_constructor_proposal"
    )
    boundary = (
        "This is a deterministic public constructor contract. Route/fiber evidence is a proxy for the "
        "learned language, target records are receptors rather than proof, and the candidate becomes a "
        "scientific result only after derivation, dimensional and closure checks, residual tests, and "
        "independent empirical or computational validation."
    )
    return ConstructorTransfer(
        source=source,
        translation=translation,
        source_identity=source_identity,
        target_identity=target_identity,
        preserved_contract=preserved,
        changed_clauses=changed,
        constructor_moves=moves,
        required_attachments=attachments,
        predictions=predictions,
        validation_gates=gates,
        readiness=readiness,
        evidence_boundary=boundary,
    )
