"""Investigate signed prediction departures independently of task acceptance.

An accepted approximation may still have a reproducible correction. Conversely,
a failed approximation is not evidence of a new phenomenon. This module keeps
the measured probabilities, sampling precision, numerical checks and repeated
calculations separate, and proposes bounded follow-up calculations.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import chi2, fisher_exact


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def adjusted_probabilities(pvalues):
    """Benjamini-Hochberg adjusted p-values, preserving the original order."""
    p = np.asarray(pvalues, float)
    if len(p) == 0:
        return []
    if not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError("Invalid comparison probabilities")
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(1., q)
    return result.tolist()


def probability_variance(p, n):
    if not (n > 0 and 0 <= p <= 1):
        raise ValueError("Invalid probability or sample count")
    # A half-count keeps the precision finite for unobserved states.
    smoothed = (p*n + .5)/(n+1)
    return smoothed*(1-smoothed)/n


def read_case(comparison_path):
    path = Path(comparison_path)
    prediction_path = path.with_name("prediction.json")
    comparison = json.loads(path.read_text(encoding="utf-8"))
    prediction = json.loads(prediction_path.read_text(encoding="utf-8"))
    if file_hash(prediction_path) != comparison["prediction_sha256"]:
        raise ValueError("The comparison belongs to another prediction")
    model = prediction.get("probability_prediction", prediction.get("repaired_prediction"))
    if model is None:
        raise ValueError("The record has no probability prediction")
    reference = comparison["reference"]
    probabilities = model["probabilities"]
    if len(probabilities) != len(reference["probabilities"]):
        raise ValueError("Different state catalogues")
    identity = {
        "task_sha256": comparison["task_sha256"],
        "states": prediction["states"]["states"],
        "implementation_sha256": prediction["implementation_sha256"],
        "policy": prediction["policy"],
    }
    group = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    rows = []
    for i, (p, q) in enumerate(zip(probabilities, reference["probabilities"])):
        variance = probability_variance(p, model["n"])+probability_variance(q, reference["n"])
        numerical = abs(p-model.get("coarse_probabilities", probabilities)[i])
        numerical += abs(q-reference.get("coarse_probabilities", reference["probabilities"])[i])
        rows.append({"state_index": i, "state": prediction["states"]["states"][i],
                     "predicted": p, "observed": q, "signed_departure": q-p,
                     "prediction_count":int(round(p*model["n"])),"prediction_n":model["n"],
                     "reference_count":int(round(q*reference["n"])),"reference_n":reference["n"],
                     "sampling_variance": variance, "standard_error": math.sqrt(variance),
                     "step_refinement_difference": numerical})
    return {"group": group, "task_sha256": comparison["task_sha256"], "path": str(path.resolve()),
            "comparison_sha256": file_hash(path), "prediction_sha256": file_hash(prediction_path),
            "prediction_seed": model["seed"], "reference_seed": reference["seed"],
            "original_acceptance": bool(comparison.get("prediction_verified", comparison.get("repair_verified", False))),
            "numerically_valid": bool(model["valid"] and reference["valid"]),
            "rate": prediction.get("task",{}).get("protocol",{}).get("rate"), "states": rows}


def investigate(paths, *, alpha=.05, maximum_followup_samples=65536):
    if not 0 < alpha < 1 or maximum_followup_samples < 256:
        raise ValueError("Invalid follow-up policy")
    cases = []
    seen = set()
    for path in paths:
        case = read_case(path)
        # A copied result or a repeated seed pair is not an independent replication.
        key = (case["group"], case["prediction_seed"], case["reference_seed"])
        if key not in seen:
            cases.append(case)
            seen.add(key)
    grouped = {}
    for case in cases:
        grouped.setdefault(case["group"], []).append(case)
    questions = []
    for group, repeats in grouped.items():
        for i in range(len(repeats[0]["states"])):
            measurements = [r["states"][i] for r in repeats]
            variance = np.array([m["sampling_variance"] for m in measurements])
            residual = np.array([m["signed_departure"] for m in measurements])
            weights = 1/variance
            effect = float(np.dot(weights, residual)/weights.sum())
            se = float(1/np.sqrt(weights.sum()))
            z = effect/se
            heterogeneity = float(np.sum(weights*(residual-effect)**2))
            heterogeneous = len(repeats)>1 and chi2.sf(heterogeneity,len(repeats)-1)<.01
            same_sign = bool(np.all(residual>0) or np.all(residual<0))
            numerical = max(m["step_refinement_difference"] for m in measurements)
            valid = all(r["numerically_valid"] for r in repeats)
            npred=sum(m["prediction_n"] for m in measurements)
            nref=sum(m["reference_n"] for m in measurements)
            cpred=sum(m["prediction_count"] for m in measurements)
            cref=sum(m["reference_count"] for m in measurements)
            exact_p=float(fisher_exact([[cref,nref-cref],[cpred,npred-cpred]]).pvalue)
            questions.append({"group": group, "task_sha256": repeats[0]["task_sha256"],
                              "state_index": i, "state": measurements[0]["state"],
                              "replications": len(repeats), "signed_departures": residual.tolist(),
                              "estimated_departure": effect, "standard_error": se,
                              "z": z, "p_value": exact_p,"probability_test":"pooled exact two-sample test",
                              "same_sign": same_sign, "heterogeneous": bool(heterogeneous),
                              "step_refinement_difference": numerical, "numerically_valid": valid,
                              "accepted_original_cases": sum(r["original_acceptance"] for r in repeats)})
    adjusted = adjusted_probabilities([q["p_value"] for q in questions])
    for row, qvalue in zip(questions, adjusted):
        row["adjusted_p_value"] = qvalue
        resolved = (row["step_refinement_difference"]==0 or
                    abs(row["estimated_departure"]) > 2*row["step_refinement_difference"])
        if not row["numerically_valid"] or not resolved:
            status = "check_numerical_resolution"
        elif row["replications"]>=2 and row["same_sign"] and not row["heterogeneous"] and qvalue<=alpha:
            status = "repeatable_departure"
        elif qvalue<=alpha:
            status = "independent_replication_needed"
        else:
            status = "sampling_precision_limited"
        # No small effect is deleted merely because it satisfies task acceptance.
        effect = abs(row["estimated_departure"])
        requested = math.ceil(4.5/max(effect,1e-9)**2)
        row["status"] = status
        row["followup"] = {
            "calculation": "step refinement" if status=="check_numerical_resolution" else "fresh independent trajectories",
            "maximum_samples": min(maximum_followup_samples, max(256, requested)),
            "precision_target_exceeds_budget": requested>maximum_followup_samples,
            "reason": "test the signed departure; task acceptance does not decide physical interest",
        }
    questions.sort(key=lambda q:q["adjusted_p_value"])
    return {"schema": "fieldbridge-physical-inquiry/1", "cases": cases, "questions": questions,
            "all_outcomes_retained": True, "alpha": alpha,
            "repeatable_departures": sum(q["status"]=="repeatable_departure" for q in questions),
            "scope": "sampling and numerical evidence for investigating a discrepancy; not a novelty test"}


def write_report(paths, output):
    report = investigate(paths)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    (out/"investigation.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    lines = ["# Signed probability departures", "",
             f"Independent calculations: {len(report['cases'])}.",
             f"Repeatable departures: {report['repeatable_departures']}.", "",
             "Original acceptance and physical follow-up are separate. All states remain in the JSON record.", "",
             "| State | Replications | Signed difference | Standard error | Adjusted p | Follow-up |",
             "| --- | ---: | ---: | ---: | ---: | --- |"]
    for q in report["questions"]:
        lines.append(f"| {q['state']} | {q['replications']} | {q['estimated_departure']:.5g} | "
                     f"{q['standard_error']:.4g} | {q['adjusted_p_value']:.4g} | {q['status']} |")
    (out/"investigation.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    return report
