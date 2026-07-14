from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import random
import re
from typing import Any, Dict, Iterable, List, Mapping, Sequence

from .pdf_sparse_builder import chunk_text, read_document
from .routes import fingerprint_text


TOKEN_RE = re.compile(r"[A-Za-z][A-Za-z0-9_-]{2,}")


def _load_manifest(path: Path) -> List[Dict[str, str]]:
    raw = path.read_text(encoding="utf-8")
    value: Any
    if path.suffix.lower() == ".jsonl":
        value = [json.loads(line) for line in raw.splitlines() if line.strip()]
    else:
        value = json.loads(raw)
        if isinstance(value, Mapping):
            value = value.get("papers") or value.get("records") or []
    if not isinstance(value, list):
        raise ValueError("Manifest must be a JSON list, a JSON object with papers/records, or JSONL.")
    required = ("paper_id", "path", "mechanism_id", "field_id")
    rows: List[Dict[str, str]] = []
    for index, row in enumerate(value):
        if not isinstance(row, Mapping):
            raise ValueError(f"Manifest row {index} is not an object.")
        missing = [key for key in required if not str(row.get(key) or "").strip()]
        if missing:
            raise ValueError(f"Manifest row {index} is missing: {', '.join(missing)}")
        resolved = Path(str(row["path"]))
        if not resolved.is_absolute():
            resolved = (path.parent / resolved).resolve()
        if not resolved.exists():
            raise FileNotFoundError(f"Paper file does not exist: {resolved}")
        rows.append({
            "paper_id": str(row["paper_id"]),
            "path": str(resolved),
            "mechanism_id": str(row["mechanism_id"]),
            "field_id": str(row["field_id"]),
        })
    ids = [row["paper_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("paper_id values must be unique; duplicate papers would invalidate the holdout.")
    return rows


def _cosine(left: Sequence[float], right: Sequence[float]) -> float:
    dot = sum(a * b for a, b in zip(left, right))
    left_norm = math.sqrt(sum(a * a for a in left))
    right_norm = math.sqrt(sum(b * b for b in right))
    return dot / (left_norm * right_norm) if left_norm and right_norm else 0.0


def _paper_fingerprint(text: str, max_chars: int) -> List[float]:
    chunks = chunk_text(text, max_chars=max_chars, overlap_sentences=0) or [text]
    vectors = [fingerprint_text(chunk).vector() for chunk in chunks if chunk.strip()]
    if not vectors:
        return []
    dimensions = len(vectors[0])
    # A full paper often contains long descriptive sections. Combining the
    # maximum with the mean retains a mechanism that is explicit in one section
    # without allowing one isolated cue to determine the entire paper vector.
    return [
        0.65 * max(vector[index] for vector in vectors)
        + 0.35 * sum(vector[index] for vector in vectors) / len(vectors)
        for index in range(dimensions)
    ]


def _tokens(text: str) -> Counter[str]:
    return Counter(token.lower() for token in TOKEN_RE.findall(text))


def _tfidf_vector(
    counter: Counter[str],
    document_frequency: Mapping[str, int],
    document_count: int,
) -> Dict[str, float]:
    total = max(sum(counter.values()), 1)
    return {
        token: (frequency / total)
        * (math.log((1 + document_count) / (1 + document_frequency.get(token, 0))) + 1.0)
        for token, frequency in counter.items()
    }


def _sparse_cosine(left: Mapping[str, float], right: Mapping[str, float]) -> float:
    common = set(left) & set(right)
    dot = sum(left[key] * right[key] for key in common)
    left_norm = math.sqrt(sum(value * value for value in left.values()))
    right_norm = math.sqrt(sum(value * value for value in right.values()))
    return dot / (left_norm * right_norm) if left_norm and right_norm else 0.0


def _metrics(query_rows: Sequence[Mapping[str, Any]], key: str, k: int) -> Dict[str, float]:
    if not query_rows:
        return {"queries": 0, "top1_accuracy": 0.0, "mrr": 0.0, f"precision_at_{k}": 0.0, f"recall_at_{k}": 0.0}
    reciprocal_ranks: List[float] = []
    precisions: List[float] = []
    recalls: List[float] = []
    top1: List[float] = []
    for row in query_rows:
        ranking = row[key]
        relevant_total = int(row["relevant_gallery_papers"])
        relevant_positions = [index for index, item in enumerate(ranking, start=1) if item["relevant"]]
        reciprocal_ranks.append(1.0 / relevant_positions[0] if relevant_positions else 0.0)
        top1.append(float(bool(ranking and ranking[0]["relevant"])))
        hits = sum(bool(item["relevant"]) for item in ranking[:k])
        precisions.append(hits / max(min(k, len(ranking)), 1))
        recalls.append(hits / max(relevant_total, 1))
    return {
        "queries": int(len(query_rows)),
        "top1_accuracy": sum(top1) / len(top1),
        "mrr": sum(reciprocal_ranks) / len(reciprocal_ranks),
        f"precision_at_{k}": sum(precisions) / len(precisions),
        f"recall_at_{k}": sum(recalls) / len(recalls),
    }


def _bootstrap_delta(
    rows: Sequence[Mapping[str, Any]],
    k: int,
    samples: int,
    seed: int,
) -> Dict[str, float]:
    if not rows:
        return {"mean": 0.0, "ci95_low": 0.0, "ci95_high": 0.0}
    rng = random.Random(seed)

    def precision(row: Mapping[str, Any], key: str) -> float:
        ranking = row[key]
        return sum(bool(item["relevant"]) for item in ranking[:k]) / max(min(k, len(ranking)), 1)

    differences = [precision(row, "operational_ranking") - precision(row, "lexical_ranking") for row in rows]
    draws = []
    for _ in range(max(samples, 1)):
        draws.append(sum(differences[rng.randrange(len(differences))] for _ in differences) / len(differences))
    draws.sort()
    low = draws[int(0.025 * (len(draws) - 1))]
    high = draws[int(0.975 * (len(draws) - 1))]
    return {"mean": sum(differences) / len(differences), "ci95_low": low, "ci95_high": high}


def validate_full_paper_zero_shot(
    manifest: Path,
    *,
    top_k: int = 10,
    max_chars: int = 2600,
    bootstrap_samples: int = 2000,
    seed: int = 17,
    min_eligible_queries: int = 100,
) -> Dict[str, Any]:
    rows = _load_manifest(manifest)
    papers: List[Dict[str, Any]] = []
    for row in rows:
        text = read_document(Path(row["path"]))
        papers.append({
            **row,
            "text_sha256": hashlib.sha256(text.encode("utf-8", errors="replace")).hexdigest(),
            "characters": len(text),
            "operational_vector": _paper_fingerprint(text, max_chars=max_chars),
            "token_counts": _tokens(text),
        })
    global_document_frequency: Counter[str] = Counter()
    for paper in papers:
        global_document_frequency.update(paper["token_counts"].keys())

    queries: List[Dict[str, Any]] = []
    skipped: List[Dict[str, str]] = []
    for index, query in enumerate(papers):
        gallery = [paper for paper in papers if paper["paper_id"] != query["paper_id"]]
        relevant = [
            paper for paper in gallery
            if paper["mechanism_id"] == query["mechanism_id"] and paper["field_id"] != query["field_id"]
        ]
        if not relevant:
            skipped.append({"paper_id": query["paper_id"], "reason": "no_same_mechanism_different_field_gallery_paper"})
            continue

        # Fit the lexical baseline on the gallery only. The query paper is
        # excluded from both retrieval and inverse-document-frequency fitting.
        gallery_document_frequency = global_document_frequency.copy()
        gallery_document_frequency.subtract(query["token_counts"].keys())
        gallery_document_frequency += Counter()
        lexical_query = _tfidf_vector(
            query["token_counts"], gallery_document_frequency, len(gallery)
        )
        lexical_gallery = {
            paper["paper_id"]: _tfidf_vector(
                paper["token_counts"], gallery_document_frequency, len(gallery)
            )
            for paper in gallery
        }

        def ranked(score_key: str) -> List[Dict[str, Any]]:
            scored = []
            for paper in gallery:
                score = (
                    _cosine(query["operational_vector"], paper["operational_vector"])
                    if score_key == "operational"
                    else _sparse_cosine(lexical_query, lexical_gallery[paper["paper_id"]])
                )
                scored.append({
                    "paper_id": paper["paper_id"],
                    "field_id": paper["field_id"],
                    "mechanism_id": paper["mechanism_id"],
                    "score": round(float(score), 8),
                    "relevant": bool(
                        paper["mechanism_id"] == query["mechanism_id"]
                        and paper["field_id"] != query["field_id"]
                    ),
                })
            return sorted(scored, key=lambda item: (-item["score"], item["paper_id"]))

        queries.append({
            "paper_id": query["paper_id"],
            "field_id": query["field_id"],
            "mechanism_id": query["mechanism_id"],
            "relevant_gallery_papers": len(relevant),
            "operational_ranking": ranked("operational"),
            "lexical_ranking": ranked("lexical"),
        })

    operational = _metrics(queries, "operational_ranking", top_k)
    lexical_metrics = _metrics(queries, "lexical_ranking", top_k)
    gain = _bootstrap_delta(queries, top_k, bootstrap_samples, seed)
    performance_claim_gate = {
        "passed": bool(len(queries) >= min_eligible_queries and gain["ci95_low"] > 0.0),
        "min_eligible_queries": int(min_eligible_queries),
        "eligible_queries_observed": int(len(queries)),
        "precision_gain_ci95_must_exclude_zero": True,
        "independent_labels_required": True,
    }
    return {
        "report_type": "fieldbridge_full_paper_zero_shot_validation",
        "readiness": "usable" if queries else "insufficient_cross_field_labels",
        "protocol": {
            "unit_of_holdout": "complete_paper",
            "query_excluded_from_gallery": True,
            "chunk_leakage_permitted": False,
            "relevance": "same mechanism_id, different field_id",
            "operational_representation": "0.65*chunkwise_max + 0.35*chunkwise_mean of public route/fiber fingerprint",
            "lexical_baseline": "full-corpus TF-IDF cosine",
            "top_k": int(top_k),
            "chunk_max_chars": int(max_chars),
            "bootstrap_samples": int(bootstrap_samples),
            "seed": int(seed),
            "lexical_fit": "gallery_only; query excluded from document-frequency estimation",
        },
        "papers": len(papers),
        "eligible_queries": len(queries),
        "skipped_queries": skipped,
        "operational": operational,
        "lexical_baseline": lexical_metrics,
        f"operational_minus_lexical_precision_at_{top_k}": gain,
        "performance_claim_gate": performance_claim_gate,
        "paper_provenance": [
            {key: paper[key] for key in ("paper_id", "path", "field_id", "mechanism_id", "characters", "text_sha256")}
            for paper in papers
        ],
        "queries": queries,
        "claim_scope": (
            "Leave-one-paper-out cross-field retrieval with complete-document holdout. "
            "The protocol establishes zero-shot evaluation only; a performance claim requires an independently "
            "labelled manifest, adequate eligible-query count, and confidence intervals excluding zero."
        ),
    }


def render_markdown(report: Mapping[str, Any]) -> str:
    k = int(report["protocol"]["top_k"])
    gain = report[f"operational_minus_lexical_precision_at_{k}"]
    return "\n".join([
        "# FieldBridge Full-Paper Zero-Shot Validation",
        "",
        f"- Readiness: `{report['readiness']}`",
        f"- Papers: `{report['papers']}`",
        f"- Eligible queries: `{report['eligible_queries']}`",
        f"- Holdout unit: complete paper",
        f"- Operational precision@{k}: `{report['operational'][f'precision_at_{k}']:.4f}`",
        f"- Lexical precision@{k}: `{report['lexical_baseline'][f'precision_at_{k}']:.4f}`",
        f"- Paired gain: `{gain['mean']:.4f}` (95% bootstrap CI `{gain['ci95_low']:.4f}`, `{gain['ci95_high']:.4f}`)",
        "",
        "A query is relevant only when another complete paper carries the same independently supplied mechanism label in a different field. No chunk from the query paper enters the gallery.",
        "",
        "## Scope",
        "",
        str(report["claim_scope"]),
    ]) + "\n"
