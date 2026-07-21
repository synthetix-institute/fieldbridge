from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple


FIELDS = (
    "next_move",
    "destination_omega",
    "destination_xi",
    "destination_completion",
)
CONTEXT_FIELDS = (
    "current_omega",
    "current_xi",
    "current_completion",
    "first_move",
)
REQUIRED_FIELDS = ("paper_id",) + CONTEXT_FIELDS + FIELDS


def _load_rows(path: Path) -> List[Dict[str, str]]:
    raw = path.read_text(encoding="utf-8")
    value: Any
    if path.suffix.lower() == ".jsonl":
        value = [json.loads(line) for line in raw.splitlines() if line.strip()]
    else:
        value = json.loads(raw)
        if isinstance(value, Mapping):
            value = value.get("transitions") or value.get("records") or []
    if not isinstance(value, list):
        raise ValueError("Continuation manifest must contain a list of transition records.")
    rows: List[Dict[str, str]] = []
    for index, row in enumerate(value):
        if not isinstance(row, Mapping):
            raise ValueError(f"Transition row {index} is not an object.")
        missing = [key for key in REQUIRED_FIELDS if str(row.get(key) or "").strip() == ""]
        if missing:
            raise ValueError(f"Transition row {index} is missing: {', '.join(missing)}")
        rows.append({key: str(row[key]) for key in REQUIRED_FIELDS})
    return rows


def _paper_fold(paper_id: str, seed: int) -> int:
    digest = hashlib.sha256(f"{seed}:{paper_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % 10


def _context(row: Mapping[str, str]) -> Tuple[str, ...]:
    return tuple(row[key] for key in CONTEXT_FIELDS)


def _fit_counts(
    rows: Sequence[Mapping[str, str]], target: str
) -> tuple[dict[Tuple[str, ...], Counter[str]], dict[str, Counter[str]], Counter[str]]:
    complete: dict[Tuple[str, ...], Counter[str]] = defaultdict(Counter)
    first_move: dict[str, Counter[str]] = defaultdict(Counter)
    global_counts: Counter[str] = Counter()
    for row in rows:
        value = row[target]
        complete[_context(row)][value] += 1
        first_move[row["first_move"]][value] += 1
        global_counts[value] += 1
    return complete, first_move, global_counts


def _distribution(
    row: Mapping[str, str],
    counts: tuple[dict[Tuple[str, ...], Counter[str]], dict[str, Counter[str]], Counter[str]],
    vocabulary: Sequence[str],
    alpha: float,
    baseline: bool,
) -> Dict[str, float]:
    complete, first_move, global_counts = counts
    selected = first_move.get(row["first_move"], global_counts) if baseline else complete.get(_context(row))
    if not selected:
        selected = first_move.get(row["first_move"], global_counts)
    total = sum(selected.values()) + alpha * len(vocabulary)
    return {value: (selected.get(value, 0) + alpha) / total for value in vocabulary}


def _score_target(
    training: Sequence[Mapping[str, str]],
    evaluation: Sequence[Mapping[str, str]],
    target: str,
    alpha: float,
) -> Dict[str, Any]:
    vocabulary = sorted({row[target] for row in training})
    counts = _fit_counts(training, target)
    correct = 0
    top3 = 0
    model_bits = 0.0
    baseline_bits = 0.0
    for row in evaluation:
        truth = row[target]
        if truth not in vocabulary:
            continue
        model = _distribution(row, counts, vocabulary, alpha, baseline=False)
        baseline = _distribution(row, counts, vocabulary, alpha, baseline=True)
        ranking = sorted(vocabulary, key=lambda value: (-model[value], value))
        correct += int(ranking[0] == truth)
        top3 += int(truth in ranking[:3])
        model_bits -= math.log2(model[truth])
        baseline_bits -= math.log2(baseline[truth])
    scored = sum(1 for row in evaluation if row[target] in vocabulary)
    return {
        "evaluation_transitions": scored,
        "accuracy": correct / scored if scored else 0.0,
        "top3_accuracy": top3 / scored if scored else 0.0,
        "bits_per_transition": model_bits / scored if scored else 0.0,
        "first_move_baseline_bits_per_transition": baseline_bits / scored if scored else 0.0,
        "gain_over_first_move_baseline_bits": (baseline_bits - model_bits) / scored if scored else 0.0,
        "vocabulary_size": len(vocabulary),
    }


def validate_future_state(
    manifest: Path,
    *,
    seed: int = 20260710,
    evaluation_fold: int = 9,
    alpha: float = 0.5,
    min_evaluation_transitions: int = 100,
) -> Dict[str, Any]:
    rows = _load_rows(manifest)
    evaluation_papers = {
        row["paper_id"] for row in rows if _paper_fold(row["paper_id"], seed) == evaluation_fold
    }
    training = [row for row in rows if row["paper_id"] not in evaluation_papers]
    evaluation = [row for row in rows if row["paper_id"] in evaluation_papers]
    scores = {target: _score_target(training, evaluation, target, alpha) for target in FIELDS}
    enough = len(evaluation) >= min_evaluation_transitions
    predictive = enough and all(scores[target]["gain_over_first_move_baseline_bits"] > 0 for target in FIELDS)
    return {
        "report_type": "fieldbridge_future_mechanism_state_validation",
        "readiness": "usable" if enough else "insufficient_evaluation_transitions",
        "question": (
            "Given the current operational state and first observed move, what move and "
            "operational state will the next equation in a withheld paper occupy?"
        ),
        "protocol": {
            "unit_of_holdout": "complete_paper",
            "query": "(current_omega, current_xi, current_completion, first_move)",
            "targets": list(FIELDS),
            "destination_definition": "frozen-language assignment of the unseen next equation",
            "current_state_classification": False,
            "evaluation_fold": evaluation_fold,
            "paper_hash_seed": seed,
            "additive_smoothing": alpha,
        },
        "transitions": len(rows),
        "training_transitions": len(training),
        "evaluation_transitions": len(evaluation),
        "training_papers": len({row["paper_id"] for row in training}),
        "evaluation_papers": len(evaluation_papers),
        "targets": scores,
        "predictive_gate": "pass" if predictive else "not_passed",
        "claim_scope": (
            "Complete-paper held-out prediction of a future symbolic move and destination state. "
            "It is not classification of the current equation, exact equation reconstruction, or "
            "validation of a physical mechanism."
        ),
    }


def render_markdown(report: Mapping[str, Any]) -> str:
    labels = {
        "next_move": "Next move",
        "destination_omega": "Future operator state",
        "destination_xi": "Future substrate state",
        "destination_completion": "Future completion state",
    }
    lines = [
        "# FieldBridge Future-State Validation",
        "",
        f"- Readiness: `{report['readiness']}`",
        f"- Held-out unit: `{report['protocol']['unit_of_holdout']}`",
        f"- Evaluation papers: `{report['evaluation_papers']}`",
        f"- Evaluation transitions: `{report['evaluation_transitions']}`",
        "",
        "## What Was Predicted",
        "",
        "| Target withheld from the model | Top-1 | Top-3 | Gain over first-move baseline |",
        "|---|---:|---:|---:|",
    ]
    for key in FIELDS:
        score = report["targets"][key]
        lines.append(
            f"| {labels[key]} | {score['accuracy']:.1%} | {score['top3_accuracy']:.1%} | "
            f"{score['gain_over_first_move_baseline_bits']:.3f} bits |"
        )
    lines.extend([
        "",
        "`Destination` is the frozen-language state assigned to an unseen next equation. "
        "It is not the classification of the current equation.",
        "",
        "## Scope",
        "",
        str(report["claim_scope"]),
    ])
    return "\n".join(lines) + "\n"
