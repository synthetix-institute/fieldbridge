"""Bind a retrieved source annotation to an explicit mathematical construction.

The annotation is authored or supplied by a future extraction stage. Its
canonical equation must match the retrieved record. This proves linkage to
that local record, not alignment to an original paper or physical validity.
The target map/observable is supplied, never inferred from retrieval scores.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, replace
import hashlib
import json

from .models import ConstructorTransfer


class CalculationRefused(ValueError):
    pass


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _strings(values, label):
    if not isinstance(values, list) or not values or not all(isinstance(v, str) and v.strip() for v in values):
        raise CalculationRefused(f"{label} must contain explicit nonempty assumptions")
    return values


def _keys(value, allowed, required, label):
    if not isinstance(value, dict):
        raise CalculationRefused(f"{label} must be an object")
    missing, extra = required - value.keys(), value.keys() - allowed
    if missing or extra:
        raise CalculationRefused(f"{label}: missing {sorted(missing)}; unsupported {sorted(extra)}")


def canonical_source_equation(source):
    """A narrow transcription format, not a general LaTeX-to-model parser."""
    if source.get("family") in {"scalar_ito", "scalar_sde"}:
        return f"dX = ({source['drift']}) dt + ({source['noise']}) dW"
    if source.get("family") == "finite_hamiltonian":
        return "H = " + json.dumps(source["hamiltonian"])
    raise CalculationRefused("Retrieved source has no supported calculation family")


def emit_spec(transfer: ConstructorTransfer, request):
    common = {"schema", "source_record_id", "target_field", "kind", "question", "assumptions"}
    optional = {"state_map", "inverse_map", "target_domain", "candidate", "observable"}
    _keys(request, common | optional, common, "Correspondence")
    if request["schema"] != "fieldbridge-correspondence/1":
        raise CalculationRefused("Unsupported correspondence schema")
    if request["target_field"] != transfer.translation.target_field.field_id:
        raise CalculationRefused("Correspondence target_field differs from --to")
    if not isinstance(request["question"], str) or not request["question"].strip():
        raise CalculationRefused("A physical question is required")
    _strings(request["assumptions"], "Correspondence assumptions")
    if not isinstance(request["source_record_id"], str) or not request["source_record_id"]:
        raise CalculationRefused("source_record_id must identify a retrieved record")
    matches = [m for m in transfer.translation.matches if m.record.record_id == request["source_record_id"]]
    if len(matches) != 1:
        raise CalculationRefused("The chosen source must occur exactly once in the actual retrieved matches")
    match = matches[0]
    source = match.record.calculation_source
    if not isinstance(source, dict):
        raise CalculationRefused("Retrieved record lacks a typed source model; no equation is guessed from its keywords")
    sde_fields = {"domain", "drift", "noise", "convention"}
    families = {"scalar_ito": ("ito_transfer", sde_fields),
                "scalar_sde": ("ito_transfer", sde_fields),
                "finite_hamiltonian": ("quantum_closure", {"hamiltonian"})}
    if source.get("family") not in families:
        raise CalculationRefused("Retrieved source has no supported calculation family")
    kind, fields = families[source["family"]]
    header = {"schema", "family", "parameters", "equation_index", "assumptions"}
    _keys(source, header | fields, header | fields, "Source annotation")
    if source["schema"] != "fieldbridge-source-model/1" or request["kind"] != kind:
        raise CalculationRefused("Source model and correspondence calculation kinds do not agree")
    if source["family"] == "scalar_ito" and source["convention"] != "ito":
        raise CalculationRefused("scalar_ito requires convention=ito; use scalar_sde for other declared conventions")
    required = {"state_map", "inverse_map", "target_domain"} if kind == "ito_transfer" else {"observable"}
    allowed = common | required | ({"candidate"} if kind == "ito_transfer" else set())
    _keys(request, allowed, common | required, "Correspondence")
    _strings(source["assumptions"], "Source assumptions")
    index = source["equation_index"]
    if type(index) is not int or index < 0 or index >= len(match.record.equations):
        raise CalculationRefused("Source equation_index is outside the retrieved record")
    stored = match.record.equations[index]
    canonical = canonical_source_equation(source)
    if "".join(stored.split()) != "".join(canonical.split()):
        raise CalculationRefused("Source annotation does not match its stored canonical equation")
    spec = {"schema": "fieldbridge-construction/1", "kind": kind,
            "question": request["question"], "parameters": deepcopy(source["parameters"]),
            "assumptions": source["assumptions"] + request["assumptions"],
            "provenance": {"origin": "retrieved_record_with_explicit_correspondence",
                "record_id": match.record.record_id, "record_source": match.record.source,
                "references": match.record.references, "equation_ids": [f"{match.record.record_id}:{index}"],
                "source_alignment": "not_established_by_this_adapter",
                "source_record_sha256": _hash(asdict(match.record)),
                "correspondence_sha256": _hash(request),
                "novelty": "not_established"}}
    spec.update({key: deepcopy(source[key]) for key in fields})
    spec.update({key: deepcopy(request[key]) for key in required})
    if "candidate" in request:
        spec["candidate"] = deepcopy(request["candidate"])
    binding = {"record_id": match.record.record_id, "source_equation": stored,
               "equation_index": index, "retrieval_score": match.score,
               "source_record_sha256": spec["provenance"]["source_record_sha256"],
               "source_model_sha256": _hash(source), "correspondence_sha256": _hash(request),
               "binding_scope": "canonical mathematical annotation agrees with the local retrieved equation; not original-paper alignment"}
    return spec, binding


def attach_calculation(transfer: ConstructorTransfer, request):
    gates = {**transfer.validation_gates, "source_record_linked": False,
             "local_generator_identity": False, "finite_quantum_closure": False,
             "omission_control_detected": False}
    try:
        spec, binding = emit_spec(transfer, request)
        from .verification import verify_construction
        report = verify_construction(spec)
        if report["calculation"] == "scalar_ito_transfer":
            if report["residual_coefficients"] != ["0", "0"]:
                raise CalculationRefused("The local generator coefficients did not verify")
            gates["local_generator_identity"] = True
            gates["omission_control_detected"] = report["omission_control"]["detects_omission"]
            predictions = [f"Derived target drift: {report['target']['drift']}.",
                           f"Derived target quadratic variation rate: {report['target']['variance']}."]
            if report["mean_prediction"]:
                predictions.append("Mean response: " + report["mean_prediction"]["expression"]
                                   + "; " + report["mean_prediction"]["condition"] + ".")
        else:
            if not report["closed_identities"] or not all(report["closed_identities"]):
                raise CalculationRefused("Observable evolution did not close")
            gates["finite_quantum_closure"] = True
            gates["omission_control_detected"] = not report["omission_control"]["single_observable_closed"]
            predictions = [report["prediction"], "Derived evolution matrix G: " + json.dumps(report["evolution_matrix"])]
        gates["source_record_linked"] = True
        return replace(transfer, calculation={"status": "calculated", "binding": binding,
                        "specification": spec, "report": report,
                        "verified_scope": "selected_retrieved_source_model_only"},
                       predictions=predictions, validation_gates=gates,
                       readiness="calculated_retrieved_source_model",
                       evidence_boundary="The selected retrieved equations and explicit correspondence produce the stated calculation. The query text and other retrieved target equations are not thereby verified. Source-paper alignment, independent target validation and novelty remain unestablished.")
    except ModuleNotFoundError as error:
        if error.name != "sympy":
            raise
        reason = "Install the calculation dependency: pip install -e '.[construction]'"
    except (CalculationRefused, ValueError, KeyError, TypeError) as error:
        reason = str(error)
    return replace(transfer, calculation={"status": "refused", "reason": reason},
                   predictions=[], validation_gates=gates, readiness="calculation_refused",
                   evidence_boundary="No calculated prediction was issued. The retrieved material remains an unverified proposal.")
