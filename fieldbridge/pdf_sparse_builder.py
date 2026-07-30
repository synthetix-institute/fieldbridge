from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .extract import extract_equations, extract_mechanism
from .routes import FIBERS, ROUTES, fingerprint_text


ROLE_PATTERNS: Dict[str, Tuple[str, List[str]]] = {
    "state": (
        "state, memory, trace, carrier, order parameter, or retained configuration",
        [
            r"\bstate\b",
            r"\bmemory\b",
            r"\btrace\b",
            r"\bcarrier\b",
            r"\border parameter\b",
            r"\bconfiguration\b",
            r"\bmicrostructure\b",
            r"\bcoverage\b",
            r"\bpotential\b",
            r"\bpolicy\b",
            r"\bbelief\b",
        ],
    ),
    "input": (
        "input, stimulus, perturbation, write operation, signal, treatment, or event",
        [
            r"\binput\b",
            r"\bstimulus\b",
            r"\bperturb",
            r"\bwrite\b",
            r"\bsignal\b",
            r"\btreatment\b",
            r"\bevent\b",
            r"\bexposure\b",
            r"\bload\b",
            r"\bfield\b",
            r"\bvoltage\b",
            r"\bgradient\b",
        ],
    ),
    "boundary": (
        "boundary, context, geometry, interface, compartment, network, or assay condition",
        [
            r"\bboundary\b",
            r"\bcontext\b",
            r"\bgeometry\b",
            r"\binterface\b",
            r"\bcompartment\b",
            r"\bnetwork\b",
            r"\bassay\b",
            r"\bsurface\b",
            r"\bmembrane\b",
            r"\bwall\b",
            r"\bniche\b",
            r"\bplatform\b",
        ],
    ),
    "output": (
        "output, response, phenotype, action, decision, transport, motion, or readout",
        [
            r"\boutput\b",
            r"\bresponse\b",
            r"\bphenotype\b",
            r"\baction\b",
            r"\bdecision\b",
            r"\btransport\b",
            r"\bmotion\b",
            r"\breadout\b",
            r"\brepair\b",
            r"\bconductance\b",
            r"\ballocation\b",
            r"\bconsensus\b",
        ],
    ),
    "validation": (
        "control, ablation, reset, washout, shuffled history, residual, or falsifier",
        [
            r"\bcontrol\b",
            r"\bablation\b",
            r"\berase",
            r"\breset\b",
            r"\bwashout\b",
            r"\bshuffle",
            r"\bresidual\b",
            r"\bfalsif",
            r"\bbaseline\b",
            r"\bnull\b",
            r"\bknock",
        ],
    ),
}

CONSTRUCTOR_ROLE_PATTERNS: Dict[str, Tuple[str, List[str]]] = {
    "state_or_carrier": (
        "the quantity, state space, memory, distribution, field, or configuration on which the mechanism acts",
        [
            r"\bstate\b",
            r"\bstates\b",
            r"\bstate space\b",
            r"\bvariable\b",
            r"\bfield\b",
            r"\bcarrier\b",
            r"\bmemory\b",
            r"\btrace\b",
            r"\border parameter\b",
            r"\bconfiguration\b",
            r"\bdensity\b",
            r"\bdistribution\b",
            r"\bconcentration\b",
            r"\bpotential\b",
            r"\bactivity\b",
            r"\bbelief\b",
            r"\bparameter\b",
        ],
    ),
    "operator_apparatus": (
        "the operator, generator, map, kernel, objective, force, or update rule that transforms the carrier",
        [
            r"\boperator\b",
            r"\bgenerator\b",
            r"\bHamiltonian\b",
            r"\bLaplacian\b",
            r"\bkernel\b",
            r"\bmap\b",
            r"\btransition matrix\b",
            r"\btransfer matrix\b",
            r"\bJacobian\b",
            r"\bpropagator\b",
            r"\bgradient\b",
            r"\bforce\b",
            r"\bfree energy\b",
            r"\bobjective\b",
            r"\boptimization\b",
            r"\bmessage passing\b",
            r"\btensor\b",
        ],
    ),
    "update_or_transport": (
        "the directed change: evolution, flow, inference, propagation, learning, transport, or recurrence",
        [
            r"\bevolution\b",
            r"\bflow\b",
            r"\btransport\b",
            r"\bdiffusion\b",
            r"\bpropagation\b",
            r"\bcurrent\b",
            r"\bflux\b",
            r"\bupdate\b",
            r"\brecurrence\b",
            r"\biteration\b",
            r"\blearning\b",
            r"\binference\b",
            r"\btrajectory\b",
            r"\bdynamics\b",
        ],
    ),
    "admissibility_logic": (
        "the constraints, closure laws, boundaries, normalizations, gates, and compatibility tests",
        [
            r"\bconstraint\b",
            r"\bclosure\b",
            r"\badmissib",
            r"\bnormalization\b",
            r"\bconservation\b",
            r"\bboundary condition\b",
            r"\bthreshold\b",
            r"\bgate\b",
            r"\bcriterion\b",
            r"\bcompatib",
            r"\bfeasib",
            r"\bresidual\b",
            r"\bregularization\b",
            r"\bselection\b",
            r"\bloss\b",
        ],
    ),
    "readout_rule": (
        "the observable, spectrum, response, decision, phenotype, output, or measurement used to test the mechanism",
        [
            r"\breadout\b",
            r"\bobservable\b",
            r"\bmeasurement\b",
            r"\bmeasure\b",
            r"\bresponse\b",
            r"\boutput\b",
            r"\bprediction\b",
            r"\bspectrum\b",
            r"\beigenvalue\b",
            r"\bmode\b",
            r"\bphenotype\b",
            r"\bbehavior\b",
            r"\bdecision\b",
            r"\bconductance\b",
            r"\berror\b",
        ],
    ),
    "protocol_execution": (
        "the preparation, intervention, ordered procedure, or computational schedule used to execute the mechanism",
        [
            r"\bprotocol\b",
            r"\bprocedure\b",
            r"\bpreparation\b",
            r"\bintervention\b",
            r"\bsequence\b",
            r"\bworkflow\b",
            r"\btime step\b",
            r"\bupdate order\b",
            r"\btraining schedule\b",
            r"\bdosage\b",
            r"\bpulse\b",
            r"\biteration\b",
            r"\balgorithm\b",
        ],
    ),
    "falsifier": (
        "the control, perturbation, ablation, reset, swap, or validation test that could break the proposed mechanism",
        [
            r"\bcontrol\b",
            r"\bvalidation\b",
            r"\bfalsif",
            r"\bablation\b",
            r"\bknockout\b",
            r"\bknockdown\b",
            r"\breset\b",
            r"\berase",
            r"\bwashout\b",
            r"\bshuffle",
            r"\bbaseline\b",
            r"\bnull\b",
            r"\bperturbation\b",
            r"\bswap\b",
            r"\bcross-validation\b",
        ],
    ),
}

SUBSTRATE_PATTERNS: Dict[str, Tuple[str, List[str]]] = {
    "coordinate_domain": (
        "ordinary spatial or spatiotemporal coordinates, domains, and physical fields",
        [
            r"\bcoordinate\b",
            r"\bspatial\b",
            r"\bposition\b",
            r"\bdomain\b",
            r"\bregion\b",
            r"\bvolume\b",
            r"\bsurface\b",
            r"\binterface\b",
            r"\bx\s*,\s*y\b",
            r"\bx\s*,\s*t\b",
            r"\\Omega",
            r"\bOmega\b",
            r"\bR\^?d\b",
            r"\\mathbb\{R\}",
        ],
    ),
    "metric_manifold": (
        "metric, manifold, curvature, geodesic, and differential-geometric substrate",
        [
            r"\bmanifold\b",
            r"\bmetric\b",
            r"\bgeodesic\b",
            r"\bcurvature\b",
            r"\bRiemann",
            r"\bRicci\b",
            r"\bconnection\b",
            r"\bLaplace-Beltrami\b",
            r"\bdifferential geometry\b",
        ],
    ),
    "inner_product_space": (
        "linear vector, Hilbert, basis, norm, inner-product, and projection substrate",
        [
            r"\bHilbert\b",
            r"\bvector space\b",
            r"\binner product\b",
            r"\bnorm\b",
            r"\bbasis\b",
            r"\bprojection\b",
            r"\borthogonal\b",
            r"\beigenvector\b",
            r"\bket\b",
            r"\bbra\b",
            r"\\langle",
            r"\\psi",
        ],
    ),
    "phase_space": (
        "Hamiltonian, symplectic, canonical, phase-space, and state-trajectory substrate",
        [
            r"\bphase space\b",
            r"\bHamiltonian\b",
            r"\bsymplectic\b",
            r"\bcanonical\b",
            r"\bPoisson\b",
            r"\bmomentum\b",
            r"\btrajectory\b",
            r"\bq\s*,\s*p\b",
        ],
    ),
    "probability_space": (
        "probability, distribution, stochastic, measure, Bayesian, and information substrate",
        [
            r"\bprobability\b",
            r"\bdistribution\b",
            r"\bdensity\b",
            r"\bstochastic\b",
            r"\bMarkov\b",
            r"\bBayes",
            r"\bmeasure\b",
            r"\bentropy\b",
            r"\bexpectation\b",
            r"\bposterior\b",
            r"\blikelihood\b",
        ],
    ),
    "graph_topology": (
        "graph, network, adjacency, edge, node, connectivity, and topological substrate",
        [
            r"\bgraph\b",
            r"\bnetwork\b",
            r"\badjacency\b",
            r"\bnode\b",
            r"\bedge\b",
            r"\bconnectivity\b",
            r"\btopolog",
            r"\bLaplacian matrix\b",
            r"\bdegree\b",
        ],
    ),
    "lattice_site_space": (
        "lattice, grid, site, chain, spin, and discrete spatial substrate",
        [
            r"\blattice\b",
            r"\bsite\b",
            r"\bgrid\b",
            r"\bspin chain\b",
            r"\bcellular automaton\b",
            r"\bnearest-neighbor\b",
            r"\bIsing\b",
        ],
    ),
    "bundle_gauge_space": (
        "fiber bundle, gauge, connection, covariant derivative, and local-frame substrate",
        [
            r"\bbundle\b",
            r"\bgauge\b",
            r"\bcovariant\b",
            r"\bconnection\b",
            r"\bcurvature form\b",
            r"\bfiber\b",
            r"\blocal frame\b",
            r"\bholonomy\b",
        ],
    ),
    "configuration_quotient": (
        "configuration, quotient, orbit, symmetry-reduced, moduli, and equivalence-class substrate",
        [
            r"\bconfiguration space\b",
            r"\bquotient\b",
            r"\borbit\b",
            r"\bsymmetry\b",
            r"\bequivalence class\b",
            r"\bmoduli\b",
            r"\bgauge fixing\b",
            r"\bstate manifold\b",
        ],
    ),
    "stoichiometric_space": (
        "chemical species, stoichiometric matrix, reaction-network, and concentration-vector substrate",
        [
            r"\bstoichiometric\b",
            r"\breaction network\b",
            r"\bchemical species\b",
            r"\bmass action\b",
            r"\brate law\b",
            r"\bconcentration vector\b",
            r"\bmetabolic\b",
        ],
    ),
}

ROLE_NODE_TYPES = {
    "input": "Input",
    "state": "State",
    "boundary": "Boundary",
    "output": "Output",
    "validation": "Validation",
}

STRUCTURAL_ROLE_EDGES = (
    ("input", "state", "writes_state", "The input, perturbation, or writing condition creates or selects the retained state."),
    ("boundary", "state", "constrains_state", "The boundary, interface, context, or geometry limits which state is admissible."),
    ("state", "output", "drives_readout", "The retained state changes the later measurable response."),
    ("boundary", "output", "gates_readout", "The boundary or context changes how the response is expressed."),
    ("state", "validation", "tested_by", "Controls test whether the retained state is really causal."),
    ("output", "validation", "falsified_by", "Controls or residuals can kill the proposed output mechanism."),
)

TOKEN_RE = re.compile(r"[A-Za-z][A-Za-z0-9_+-]{2,}")
SENTENCE_RE = re.compile(r"(?<=[.!?])\s+")

STOPWORDS = {
    "about",
    "after",
    "against",
    "also",
    "and",
    "are",
    "because",
    "been",
    "before",
    "between",
    "both",
    "but",
    "can",
    "could",
    "does",
    "each",
    "from",
    "have",
    "into",
    "may",
    "more",
    "not",
    "only",
    "over",
    "paper",
    "such",
    "that",
    "the",
    "then",
    "there",
    "these",
    "this",
    "those",
    "through",
    "under",
    "using",
    "when",
    "where",
    "which",
    "with",
}


@dataclass(frozen=True)
class Chunk:
    doc_id: str
    path: str
    chunk_id: str
    text: str
    title: str
    routes: Dict[str, float]
    fibers: Dict[str, float]
    role_scores: Dict[str, float]
    constructor_scores: Dict[str, float]
    substrate_scores: Dict[str, float]
    equations: List[str]


@dataclass
class KGNode:
    id: str
    node_type: str
    label: str
    properties: Dict[str, Any] = field(default_factory=dict)
    vector_id: str = ""
    source_file: str = ""
    content: str = ""
    original_ids: List[str] = field(default_factory=list)


@dataclass
class KGEdge:
    source: str
    target: str
    relationship_type: str
    weight: float = 1.0
    vector_id: str = ""
    source_file: str = ""
    evidence: Dict[str, Any] = field(default_factory=dict)


def clean_text(value: Any, limit: int = 240) -> str:
    text = re.sub(r"\s+", " ", str(value or "")).strip()
    return text if len(text) <= limit else text[: limit - 3].rstrip() + "..."


def slugify(value: str, fallback: str = "field") -> str:
    text = re.sub(r"[^a-zA-Z0-9]+", "_", value.lower()).strip("_")
    return text or fallback


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def tokens(text: str) -> List[str]:
    return [
        token.lower()
        for token in TOKEN_RE.findall(text)
        if token.lower() not in STOPWORDS and not token.lower().startswith("http")
    ]


def top_terms(texts: Iterable[str], limit: int = 12) -> List[str]:
    counts: Counter = Counter()
    for text in texts:
        counts.update(tokens(text))
    return [term for term, _ in counts.most_common(limit)]


def pattern_score(text: str, patterns: Sequence[str]) -> Tuple[float, List[str]]:
    hits = [pattern for pattern in patterns if re.search(pattern, text, re.IGNORECASE)]
    return max(0.0, min(1.0, 1.0 - math.exp(-0.42 * len(hits)))), hits


def extract_pdf_text(path: Path) -> str:
    """Extract text from a PDF with optional local tools.

    The concept-economics prototype used ``pdfplumber`` because it preserves
    page text better for many scientific PDFs. FieldBridge keeps every parser
    optional: it tries ``pdfplumber`` first, then ``pypdf``/``PyPDF2``, then the
    command-line ``pdftotext`` fallback.
    """

    try:
        import pdfplumber  # type: ignore

        pages = []
        with pdfplumber.open(str(path)) as pdf:
            for page in pdf.pages:
                pages.append(page.extract_text() or "")
        text = "\n\n".join(pages).strip()
        if text:
            return text
    except Exception:
        pass

    for module_name in ("pypdf", "PyPDF2"):
        try:
            module = __import__(module_name)
            reader_cls = getattr(module, "PdfReader")
            reader = reader_cls(str(path))
            pages = []
            for page in reader.pages:
                pages.append(page.extract_text() or "")
            text = "\n\n".join(pages).strip()
            if text:
                return text
        except Exception:
            pass

    pdftotext = shutil.which("pdftotext")
    if pdftotext:
        result = subprocess.run(
            [pdftotext, "-layout", str(path), "-"],
            check=False,
            capture_output=True,
            text=True,
            errors="replace",
        )
        if result.stdout.strip():
            return result.stdout

    raise RuntimeError(
        f"Could not extract text from {path}. Install pypdf/PyPDF2 or poppler pdftotext."
    )


def read_document(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".pdf":
        return extract_pdf_text(path)
    return path.read_text(encoding="utf-8", errors="replace")


def iter_documents(root: Path, extensions: Sequence[str]) -> List[Path]:
    suffixes = {ext if ext.startswith(".") else f".{ext}" for ext in extensions}
    return sorted(path for path in root.rglob("*") if path.is_file() and path.suffix.lower() in suffixes)


def chunk_text(text: str, max_chars: int = 2600, overlap_sentences: int = 1) -> List[str]:
    paragraphs = [part.strip() for part in re.split(r"\n\s*\n", text) if part.strip()]
    chunks: List[str] = []
    carry: List[str] = []
    current = ""
    for paragraph in paragraphs:
        candidate = f"{current}\n\n{paragraph}".strip() if current else paragraph
        if len(candidate) <= max_chars:
            current = candidate
            continue
        if current:
            chunks.append(current)
            sentences = SENTENCE_RE.split(current)
            carry = sentences[-overlap_sentences:] if overlap_sentences else []
        current = " ".join(carry + [paragraph]).strip()
        if len(current) > max_chars:
            for start in range(0, len(current), max_chars):
                piece = current[start : start + max_chars].strip()
                if piece:
                    chunks.append(piece)
            current = ""
            carry = []
    if current:
        chunks.append(current)
    return chunks


def make_chunks(
    paths: Sequence[Path],
    max_chars: int,
    max_chunks_per_doc: int,
    failures: Optional[List[Dict[str, str]]] = None,
) -> List[Chunk]:
    chunks: List[Chunk] = []
    for doc_index, path in enumerate(paths):
        try:
            text = read_document(path)
        except Exception as exc:
            if failures is not None:
                failures.append({"path": str(path), "error": clean_text(exc, 500)})
            continue
        if not text.strip():
            if failures is not None:
                failures.append({"path": str(path), "error": "No extractable text."})
            continue
        doc_id = f"D{doc_index:04d}"
        for chunk_index, piece in enumerate(chunk_text(text, max_chars=max_chars)):
            if chunk_index >= max_chunks_per_doc:
                break
            fp = fingerprint_text(piece)
            role_scores = {
                role: pattern_score(piece, patterns)[0]
                for role, (_, patterns) in ROLE_PATTERNS.items()
            }
            constructor_scores = {
                role: pattern_score(piece, patterns)[0]
                for role, (_, patterns) in CONSTRUCTOR_ROLE_PATTERNS.items()
            }
            substrate_scores = {
                substrate: pattern_score(piece, patterns)[0]
                for substrate, (_, patterns) in SUBSTRATE_PATTERNS.items()
            }
            chunks.append(
                Chunk(
                    doc_id=doc_id,
                    path=str(path),
                    chunk_id=f"{doc_id}:C{chunk_index:03d}",
                    title=path.stem,
                    text=piece,
                    routes=fp.routes,
                    fibers=fp.fibers,
                    role_scores=role_scores,
                    constructor_scores=constructor_scores,
                    substrate_scores=substrate_scores,
                    equations=extract_equations(piece, limit=4),
                )
            )
    return chunks


def score_chunk(chunk: Chunk) -> float:
    route_values = sorted(chunk.routes.values(), reverse=True)
    fiber_values = sorted(chunk.fibers.values(), reverse=True)
    role_values = sorted(chunk.role_scores.values(), reverse=True)
    constructor_values = sorted(chunk.constructor_scores.values(), reverse=True)
    substrate_values = sorted(chunk.substrate_scores.values(), reverse=True)
    equation_bonus = 0.16 if chunk.equations else 0.0
    return (
        0.28 * sum(route_values[:3])
        + 0.18 * sum(fiber_values[:2])
        + 0.14 * sum(role_values[:2])
        + 0.19 * sum(constructor_values[:3])
        + 0.12 * sum(substrate_values[:2])
        + equation_bonus
    )


def top_chunks(chunks: Sequence[Chunk], scorer, limit: int) -> List[Chunk]:
    return sorted(chunks, key=scorer, reverse=True)[:limit]


def sparse_attention(chunks: Sequence[Chunk], top_k: int = 10) -> Dict[str, Any]:
    heads: Dict[str, Any] = {}
    for route_id, (meaning, _) in ROUTES.items():
        selected = top_chunks(chunks, lambda chunk, rid=route_id: chunk.routes.get(rid, 0.0), top_k)
        heads[f"route:{route_id}"] = {
            "meaning": meaning,
            "top_chunks": chunk_rows(selected, route_id=route_id),
            "terms": top_terms(chunk.text for chunk in selected),
        }
    for fiber_id, (meaning, _) in FIBERS.items():
        selected = top_chunks(chunks, lambda chunk, fid=fiber_id: chunk.fibers.get(fid, 0.0), top_k)
        heads[f"fiber:{fiber_id}"] = {
            "meaning": meaning,
            "top_chunks": chunk_rows(selected, fiber_id=fiber_id),
            "terms": top_terms(chunk.text for chunk in selected),
        }
    for role, (meaning, _) in ROLE_PATTERNS.items():
        selected = top_chunks(chunks, lambda chunk, role_id=role: chunk.role_scores.get(role_id, 0.0), top_k)
        heads[f"role:{role}"] = {
            "meaning": meaning,
            "top_chunks": chunk_rows(selected, role=role),
            "terms": top_terms(chunk.text for chunk in selected),
        }
    for role, (meaning, _) in CONSTRUCTOR_ROLE_PATTERNS.items():
        selected = top_chunks(chunks, lambda chunk, role_id=role: chunk.constructor_scores.get(role_id, 0.0), top_k)
        heads[f"constructor:{role}"] = {
            "meaning": meaning,
            "top_chunks": chunk_rows(selected, constructor_role=role),
            "terms": top_terms(chunk.text for chunk in selected),
        }
    for substrate, (meaning, _) in SUBSTRATE_PATTERNS.items():
        selected = top_chunks(chunks, lambda chunk, substrate_id=substrate: chunk.substrate_scores.get(substrate_id, 0.0), top_k)
        heads[f"substrate:{substrate}"] = {
            "meaning": meaning,
            "top_chunks": chunk_rows(selected, substrate=substrate),
            "terms": top_terms(chunk.text for chunk in selected),
        }
    return heads


def chunk_rows(
    chunks: Sequence[Chunk],
    route_id: str = "",
    fiber_id: str = "",
    role: str = "",
    constructor_role: str = "",
    substrate: str = "",
) -> List[Dict[str, Any]]:
    rows = []
    for chunk in chunks:
        score = (
            chunk.routes.get(route_id, 0.0)
            if route_id
            else chunk.fibers.get(fiber_id, 0.0)
            if fiber_id
            else chunk.role_scores.get(role, 0.0)
            if role
            else chunk.constructor_scores.get(constructor_role, 0.0)
            if constructor_role
            else chunk.substrate_scores.get(substrate, 0.0)
            if substrate
            else score_chunk(chunk)
        )
        rows.append(
            {
                "chunk_id": chunk.chunk_id,
                "document": chunk.title,
                "path": chunk.path,
                "score": round(float(score), 4),
                "excerpt": clean_text(chunk.text, 340),
            }
        )
    return rows


def average_scores(chunks: Sequence[Chunk], attr: str, keys: Sequence[str]) -> Dict[str, float]:
    if not chunks:
        return {key: 0.0 for key in keys}
    rows = [getattr(chunk, attr) for chunk in chunks]
    return {
        key: round(sum(float(row.get(key, 0.0)) for row in rows) / len(rows), 4)
        for key in keys
    }


def active_keys(scores: Mapping[str, float], limit: int = 3, threshold: float = 0.18) -> List[str]:
    return [
        key
        for key, value in sorted(scores.items(), key=lambda item: item[1], reverse=True)
        if value >= threshold
    ][:limit]


def _phrase_candidates_from_patterns(
    chunks: Sequence[Chunk],
    key: str,
    pattern_map: Mapping[str, Tuple[str, List[str]]],
    score_attr: str,
    fallback: Sequence[str],
    limit: int = 10,
) -> List[str]:
    selected = top_chunks(chunks, lambda chunk: getattr(chunk, score_attr).get(key, 0.0), max(8, limit))
    phrases: List[str] = []
    for chunk in selected:
        sentences = SENTENCE_RE.split(chunk.text)
        for sentence in sentences:
            if pattern_score(sentence, pattern_map[key][1])[0] <= 0:
                continue
            clean = clean_text(sentence, 110)
            if 16 <= len(clean) <= 110:
                phrases.append(clean)
    phrases.extend(fallback)
    out = []
    seen = set()
    for phrase in phrases:
        key = phrase.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(phrase)
        if len(out) >= limit:
            break
    return out


def phrase_candidates(chunks: Sequence[Chunk], role: str, fallback: Sequence[str], limit: int = 10) -> List[str]:
    return _phrase_candidates_from_patterns(chunks, role, ROLE_PATTERNS, "role_scores", fallback, limit)


def constructor_phrase_candidates(chunks: Sequence[Chunk], role: str, fallback: Sequence[str], limit: int = 10) -> List[str]:
    return _phrase_candidates_from_patterns(
        chunks, role, CONSTRUCTOR_ROLE_PATTERNS, "constructor_scores", fallback, limit
    )


def field_pack_from_attention(field_id: str, label: str, description: str, chunks: Sequence[Chunk]) -> Dict[str, Any]:
    return {
        "field_id": field_id,
        "label": label,
        "description": description or f"Sparse-attention field pack built from {len(chunks)} document chunks.",
        "state_words": constructor_phrase_candidates(chunks, "state_or_carrier", ["state or carrier q"], 10),
        "input_words": phrase_candidates(chunks, "input", ["input u", "write operation", "perturbation"], 10),
        "boundary_words": constructor_phrase_candidates(chunks, "admissibility_logic", ["admissibility condition", "boundary condition B"], 10),
        "output_words": constructor_phrase_candidates(chunks, "readout_rule", ["output y", "later response", "readout"], 10),
        "validation_words": constructor_phrase_candidates(chunks, "falsifier", ["erased-state control", "shuffled-history control", "boundary swap"], 10),
        "protocol_words": constructor_phrase_candidates(
            chunks,
            "protocol_execution",
            ["preparation or execution protocol"],
            10,
        ),
    }


def evidence_rows(
    chunks: Sequence[Chunk],
    score_attr: str,
    key: str,
    limit: int = 5,
) -> List[Dict[str, Any]]:
    selected = top_chunks(chunks, lambda chunk: getattr(chunk, score_attr).get(key, 0.0), limit)
    return [
        {
            "chunk_id": chunk.chunk_id,
            "document": chunk.title,
            "path": chunk.path,
            "score": round(float(getattr(chunk, score_attr).get(key, 0.0)), 4),
            "excerpt": clean_text(chunk.text, 320),
            "equations": chunk.equations[:3],
        }
        for chunk in selected
        if getattr(chunk, score_attr).get(key, 0.0) > 0
    ]


def profile_rows(
    scores: Mapping[str, float],
    definitions: Mapping[str, Tuple[str, List[str]]],
    chunks: Sequence[Chunk],
    score_attr: str,
    threshold: float = 0.05,
) -> List[Dict[str, Any]]:
    rows = []
    for key, score in sorted(scores.items(), key=lambda item: item[1], reverse=True):
        if score < threshold:
            continue
        meaning = definitions[key][0]
        rows.append(
            {
                "id": key,
                "score": round(float(score), 4),
                "definition": meaning,
                "terms": top_terms(
                    (chunk.text for chunk in top_chunks(chunks, lambda c, k=key: getattr(c, score_attr).get(k, 0.0), 8)),
                    10,
                ),
                "evidence": evidence_rows(chunks, score_attr, key, limit=4),
            }
        )
    return rows


ROUTE_TO_CONSTRUCTOR_ROLES: Dict[str, List[str]] = {
    "transport_flow_route": ["state_or_carrier", "update_or_transport", "readout_rule"],
    "constraint_closure_route": ["state_or_carrier", "admissibility_logic", "falsifier"],
    "spectral_operator_route": ["operator_apparatus", "readout_rule", "state_or_carrier"],
    "boundary_weak_form_route": ["admissibility_logic", "state_or_carrier", "readout_rule"],
    "commutator_incompatibility_route": ["operator_apparatus", "admissibility_logic", "falsifier"],
    "discrete_protocol_route": [
        "update_or_transport",
        "protocol_execution",
        "admissibility_logic",
        "readout_rule",
    ],
}


def route_adapter_rows(
    chunks: Sequence[Chunk],
    route_priors: Mapping[str, float],
    substrate_priors: Mapping[str, float],
    limit: int = 6,
) -> List[Dict[str, Any]]:
    rows = []
    active_substrates = [
        key for key, value in sorted(substrate_priors.items(), key=lambda item: item[1], reverse=True) if value >= 0.05
    ][:4]
    for route_id, score in sorted(route_priors.items(), key=lambda item: item[1], reverse=True):
        if score < 0.05:
            continue
        selected = top_chunks(chunks, lambda chunk, rid=route_id: chunk.routes.get(rid, 0.0), limit)
        roles = ROUTE_TO_CONSTRUCTOR_ROLES.get(route_id, [])
        rows.append(
            {
                "hyperion_route": route_id,
                "score": round(float(score), 4),
                "route_definition": ROUTES[route_id][0],
                "constructor_roles": roles,
                "field_terms": top_terms((chunk.text for chunk in selected), 14),
                "substrate_links": active_substrates,
                "evidence": chunk_rows(selected, route_id=route_id),
            }
        )
    return rows


def build_field_adapter(
    field_id: str,
    label: str,
    description: str,
    chunks: Sequence[Chunk],
    selected: Sequence[Chunk],
    paths: Sequence[Path],
    generated_at: str,
) -> Dict[str, Any]:
    route_priors = average_scores(selected, "routes", list(ROUTES.keys()))
    fiber_priors = average_scores(selected, "fibers", list(FIBERS.keys()))
    constructor_priors = average_scores(selected, "constructor_scores", list(CONSTRUCTOR_ROLE_PATTERNS.keys()))
    substrate_priors = average_scores(selected, "substrate_scores", list(SUBSTRATE_PATTERNS.keys()))
    constructor_profiles = profile_rows(
        constructor_priors,
        CONSTRUCTOR_ROLE_PATTERNS,
        selected,
        "constructor_scores",
        threshold=0.04,
    )
    substrate_profiles = profile_rows(
        substrate_priors,
        SUBSTRATE_PATTERNS,
        selected,
        "substrate_scores",
        threshold=0.025,
    )
    route_map = route_adapter_rows(selected, route_priors, substrate_priors)
    active_routes = [row["hyperion_route"] for row in route_map[:4]]
    active_substrates = [row["id"] for row in substrate_profiles[:4]]
    return {
        "schema_version": 1,
        "artifact_type": "fieldbridge_field_adapter",
        "field_id": field_id,
        "label": label,
        "description": description or f"Adapter built from {len(paths)} field documents.",
        "generated_at": generated_at,
        "corpus": {
            "documents": len(paths),
            "chunks": len(chunks),
            "anchor_chunks": len(selected),
            "source_paths": [str(path) for path in paths],
        },
        "adapter_contract": {
            "purpose": (
                "Map a portable Hyperion route/fiber fingerprint into field-native "
                "constructor roles while keeping the substrate evidence explicit."
            ),
            "input": "mechanism sheet or route/fiber fingerprint",
            "output": (
                "field-native state/carrier, operator, update, admissibility, "
                "readout, protocol, and falsifier candidates"
            ),
            "review_rule": (
                "Promote a translation only when evidence exists for carrier, operator/update, "
                "admissibility, readout, execution protocol, and a falsifying control "
                "in the field corpus."
            ),
        },
        "route_profile": route_priors,
        "fiber_profile": fiber_priors,
        "constructor_role_profile": constructor_priors,
        "substrate_profile": substrate_priors,
        "active_routes": active_routes,
        "active_substrates": active_substrates,
        "constructor_roles": constructor_profiles,
        "universal_substrates": substrate_profiles,
        "route_to_field_adapter": route_map,
        "field_native_receivers": {
            "state_or_carrier": constructor_phrase_candidates(selected, "state_or_carrier", [], 12),
            "operator_apparatus": constructor_phrase_candidates(selected, "operator_apparatus", [], 12),
            "update_or_transport": constructor_phrase_candidates(selected, "update_or_transport", [], 12),
            "admissibility_logic": constructor_phrase_candidates(selected, "admissibility_logic", [], 12),
            "readout_rule": constructor_phrase_candidates(selected, "readout_rule", [], 12),
            "protocol_execution": constructor_phrase_candidates(selected, "protocol_execution", [], 12),
            "falsifier": constructor_phrase_candidates(selected, "falsifier", [], 12),
        },
        "gap_report": {
            "missing_constructor_roles": [
                key for key, value in constructor_priors.items() if value < 0.04
            ],
            "weak_substrate_evidence": [
                key for key, value in substrate_priors.items() if value < 0.025
            ],
            "interpretation": (
                "Missing or weak entries are not failures of the field. They mark where "
                "additional PDFs, equations, experiments, or manual curation are needed "
                "before a transfer should be treated as evidence-backed."
            ),
        },
        "claim_scope": (
            "Deterministic adapter built from sparse attention over a field PDF folder. "
            "It is a receptor map for mechanism translation, not a proof that a proposed "
            "translation is physically true."
        ),
    }


def render_adapter_markdown(adapter: Mapping[str, Any]) -> str:
    def score_list(rows: Mapping[str, float], limit: int = 8) -> str:
        items = sorted(rows.items(), key=lambda item: item[1], reverse=True)[:limit]
        return "\n".join(f"- `{key}`: `{value:.3f}`" for key, value in items)

    lines = [
        f"# FieldBridge Adapter: {adapter['label']}",
        "",
        "## Corpus",
        f"- Documents: `{adapter['corpus']['documents']}`",
        f"- Chunks: `{adapter['corpus']['chunks']}`",
        f"- Anchor chunks: `{adapter['corpus']['anchor_chunks']}`",
        "",
        "## Active Hyperion Routes",
        score_list(adapter["route_profile"]),
        "",
        "## Universal Substrate Evidence",
        score_list(adapter["substrate_profile"]),
        "",
        "## Constructor Role Evidence",
        score_list(adapter["constructor_role_profile"]),
        "",
        "## Adapter Contract",
        str(adapter["adapter_contract"]["purpose"]),
        "",
        "A translation should name a field-native carrier, operator or update rule, "
        "admissibility condition, readout, execution protocol, and falsifying control "
        "before it is promoted.",
        "",
        "## Route To Field Adapter",
    ]
    for row in adapter["route_to_field_adapter"]:
        lines.extend(
            [
                f"### `{row['hyperion_route']}`",
                f"- Score: `{row['score']:.3f}`",
                f"- Route definition: {row['route_definition']}",
                f"- Constructor roles: {', '.join(row['constructor_roles'])}",
                f"- Substrate links: {', '.join(row['substrate_links']) or 'none detected'}",
                f"- Field terms: {', '.join(row['field_terms'][:10])}",
            ]
        )
        if row["evidence"]:
            first = row["evidence"][0]
            lines.append(f"- Evidence: `{first['document']}` - {first['excerpt']}")
        lines.append("")
    lines.extend(
        [
            "## Gaps",
            f"- Missing constructor roles: {', '.join(adapter['gap_report']['missing_constructor_roles']) or 'none above threshold'}",
            f"- Weak substrate evidence: {', '.join(adapter['gap_report']['weak_substrate_evidence']) or 'none above threshold'}",
            "",
            "## Scope",
            str(adapter["claim_scope"]),
            "",
        ]
    )
    return "\n".join(lines)


def mechanism_record_from_chunk(field_id: str, chunk: Chunk, rank: int) -> Dict[str, Any]:
    sheet = extract_mechanism(chunk.text, title=chunk.title)
    routes = active_keys(chunk.routes)
    fibers = active_keys(chunk.fibers)
    invariant = (
        f"{field_id} sparse-attention anchor is indexed by "
        f"{', '.join(routes) or 'no strong route'} over "
        f"{', '.join(fibers) or 'no strong fiber'}."
    )
    equations = chunk.equations or sheet.equations

    def grounded_role(role: str, limit: int = 3) -> List[str]:
        return constructor_phrase_candidates([chunk], role, [], limit)

    variables: List[str] = []
    for role in (
        "state_or_carrier",
        "operator_apparatus",
        "update_or_transport",
        "admissibility_logic",
        "protocol_execution",
    ):
        variables.extend(f"{role}: {value}" for value in grounded_role(role, 2))
    measurements = grounded_role("readout_rule", 4)
    controls = grounded_role("falsifier", 4)
    return {
        "record_id": f"{field_id}:pdf_anchor:{rank:04d}",
        "title": f"{chunk.title} mechanism anchor {rank:03d}",
        "field_id": field_id,
        "summary": clean_text(chunk.text, 420),
        "invariant": invariant,
        "equations": equations[:6],
        "variables": variables,
        "measurements": measurements,
        "controls": controls,
        "references": [
            f"{Path(chunk.path).suffix.lower().lstrip('.') or 'document'}:"
            f"{Path(chunk.path).name}"
        ],
        "keywords": top_terms([chunk.text], 24),
        "routes": chunk.routes,
        "fibers": chunk.fibers,
        "source": "document_sparse_attention",
    }


def node_row(node: KGNode) -> Dict[str, Any]:
    return {
        "id": node.id,
        "label": node.label,
        "node_type": node.node_type,
        "properties": node.properties,
        "vector_id": node.vector_id,
        "source_file": node.source_file,
        "content": node.content,
        "original_ids": node.original_ids,
    }


def edge_row(edge: KGEdge) -> Dict[str, Any]:
    return {
        "source": edge.source,
        "target": edge.target,
        "relationship": edge.relationship_type,
        "relationship_type": edge.relationship_type,
        "weight": round(float(edge.weight), 4),
        "vector_id": edge.vector_id,
        "source_file": edge.source_file,
        "evidence": edge.evidence,
    }


def compute_graph_metrics(nodes: List[Dict[str, Any]], edges: List[Dict[str, Any]]) -> Dict[str, Any]:
    ids = [node["id"] for node in nodes]
    node_set = set(ids)
    in_degree: Counter = Counter()
    out_degree: Counter = Counter()
    weighted_out: Dict[str, float] = defaultdict(float)
    incoming: Dict[str, List[Tuple[str, float]]] = defaultdict(list)
    for edge in edges:
        source = edge["source"]
        target = edge["target"]
        if source not in node_set or target not in node_set:
            continue
        weight = float(edge.get("weight", 1.0))
        out_degree[source] += 1
        in_degree[target] += 1
        weighted_out[source] += max(weight, 0.001)
        incoming[target].append((source, max(weight, 0.001)))

    if ids:
        rank = {node_id: 1.0 / len(ids) for node_id in ids}
        damping = 0.85
        for _ in range(24):
            base = (1.0 - damping) / len(ids)
            next_rank = {node_id: base for node_id in ids}
            dangling = sum(rank[node_id] for node_id in ids if weighted_out.get(node_id, 0.0) <= 0)
            for node_id in ids:
                next_rank[node_id] += damping * dangling / len(ids)
                for source, weight in incoming.get(node_id, []):
                    next_rank[node_id] += damping * rank[source] * weight / weighted_out[source]
            rank = next_rank
    else:
        rank = {}

    for node in nodes:
        node_id = node["id"]
        node["in_degree"] = int(in_degree[node_id])
        node["out_degree"] = int(out_degree[node_id])
        node["degree"] = int(in_degree[node_id] + out_degree[node_id])
        node["pagerank"] = round(float(rank.get(node_id, 0.0)), 6)

    relationship_counts = Counter(edge["relationship"] for edge in edges)
    node_type_counts = Counter(node["node_type"] for node in nodes)
    return {
        "nodes": len(nodes),
        "edges": len(edges),
        "node_types": dict(sorted(node_type_counts.items())),
        "relationship_types": dict(sorted(relationship_counts.items())),
    }


def detect_receptor_gaps(nodes: List[Dict[str, Any]], edges: List[Dict[str, Any]], role_terms: Mapping[str, Sequence[str]]) -> Dict[str, Any]:
    node_ids = {node["id"] for node in nodes}
    edge_types = {edge["relationship"] for edge in edges}
    gaps: List[Dict[str, Any]] = []

    for role in ROLE_PATTERNS:
        if not role_terms.get(role):
            gaps.append(
                {
                    "gap_type": "missing_role_language",
                    "severity": 0.7,
                    "description": f"No strong field phrase was found for the {role} role.",
                    "involved_nodes": [f"role:{role}"] if f"role:{role}" in node_ids else [],
                    "suggested_fixes": [f"Add papers or notes that name the field-specific {role} variable."],
                }
            )

    required_links = {
        "writes_state": "The corpus names inputs and states, but does not clearly connect the input as the writing operation.",
        "constrains_state": "The corpus names states and boundaries, but the admissibility link is weak.",
        "drives_readout": "The corpus names retained states and outputs, but the causal readout link is weak.",
        "tested_by": "The corpus names mechanisms, but explicit controls or falsifiers are weak.",
    }
    for relationship, description in required_links.items():
        if relationship not in edge_types:
            gaps.append(
                {
                    "gap_type": "missing_mechanism_edge",
                    "severity": 0.82,
                    "description": description,
                    "involved_nodes": [],
                    "suggested_fixes": ["Add or inspect documents with explicit experiment-to-state-to-readout language."],
                }
            )

    isolated = [
        node
        for node in nodes
        if int(node.get("degree", 0)) == 0 and not node["id"].startswith(("role:", "route:", "fiber:"))
    ][:40]
    similar = []
    for index, left in enumerate(nodes):
        for right in nodes[index + 1 :]:
            if left["node_type"] != right["node_type"]:
                continue
            score = SequenceMatcher(None, left["label"].lower(), right["label"].lower()).ratio()
            if 0.72 <= score < 0.99:
                similar.append(
                    {
                        "node1": left["id"],
                        "node2": right["id"],
                        "label1": left["label"],
                        "label2": right["label"],
                        "similarity": round(score, 3),
                    }
                )

    return {
        "gaps": gaps,
        "isolated_nodes": [
            {"node_id": node["id"], "label": node["label"], "type": node["node_type"]}
            for node in isolated
        ],
        "similar_nodes": sorted(similar, key=lambda row: row["similarity"], reverse=True)[:40],
    }


def receptor_graph(chunks: Sequence[Chunk], max_nodes: int = 120, max_edges: int = 260) -> Dict[str, Any]:
    role_terms = {role: phrase_candidates(chunks, role, [], 18) for role in ROLE_PATTERNS}
    nodes: Dict[str, KGNode] = {}
    edge_weights: Dict[Tuple[str, str, str], KGEdge] = {}

    def add_node(node_id: str, node_type: str, label: str, **properties: Any) -> None:
        if node_id in nodes:
            nodes[node_id].properties.update({key: value for key, value in properties.items() if value is not None})
            return
        nodes[node_id] = KGNode(
            id=node_id,
            node_type=node_type,
            label=label,
            properties={key: value for key, value in properties.items() if value is not None},
            vector_id=str(properties.get("vector_id", "")),
            source_file=str(properties.get("source_file", "")),
            content=str(properties.get("content", "")),
            original_ids=[node_id],
        )

    def add_edge(
        source: str,
        target: str,
        relationship: str,
        weight: float = 1.0,
        vector_id: str = "",
        source_file: str = "",
        evidence: Optional[Dict[str, Any]] = None,
    ) -> None:
        if source == target or source not in nodes or target not in nodes:
            return
        key = (source, target, relationship)
        if key not in edge_weights:
            edge_weights[key] = KGEdge(
                source=source,
                target=target,
                relationship_type=relationship,
                weight=0.0,
                vector_id=vector_id,
                source_file=source_file,
                evidence={},
            )
        edge_weights[key].weight += weight
        if evidence:
            supports = edge_weights[key].evidence.setdefault("supporting_chunks", [])
            evidence_row = dict(evidence)
            if len(supports) < 8 and evidence_row not in supports:
                supports.append(evidence_row)

    for role, (meaning, _) in ROLE_PATTERNS.items():
        add_node(f"role:{role}", "MechanismRole", role.replace("_", " "), role=role, meaning=meaning)
        for rank, phrase in enumerate(role_terms[role][:10]):
            term_id = f"term:{role}:{slugify(phrase)[:72]}"
            add_node(term_id, ROLE_NODE_TYPES[role], clean_text(phrase, 120), role=role, rank=rank)
            add_edge(f"role:{role}", term_id, "has_field_phrase", weight=max(0.2, 1.0 - 0.06 * rank))

    route_priors = average_scores(chunks, "routes", list(ROUTES.keys()))
    fiber_priors = average_scores(chunks, "fibers", list(FIBERS.keys()))
    for route_id in active_keys(route_priors, limit=len(ROUTES), threshold=0.08):
        meaning = ROUTES[route_id][0]
        add_node(f"route:{route_id}", "Route", route_id.replace("_route", "").replace("_", " "), score=route_priors[route_id], meaning=meaning)
    for fiber_id in active_keys(fiber_priors, limit=len(FIBERS), threshold=0.08):
        meaning = FIBERS[fiber_id][0]
        add_node(f"fiber:{fiber_id}", "Fiber", fiber_id.replace("_", " "), score=fiber_priors[fiber_id], meaning=meaning)

    concept_counts: Counter = Counter()
    concept_sources: Dict[str, str] = {}
    for chunk in chunks:
        chunk_terms = [term for term in tokens(chunk.text) if 4 <= len(term) <= 32]
        for term, count in Counter(chunk_terms).most_common(16):
            concept_counts[term] += count
        for term in chunk_terms[:60]:
            concept_sources.setdefault(term, chunk.path)
    for term, count in concept_counts.most_common(max(20, max_nodes // 2)):
        concept_id = f"concept:{slugify(term)}"
        add_node(concept_id, "Concept", term, count=int(count), source_file=concept_sources.get(term, ""))

    role_coactivation: Counter = Counter()
    concept_coactivation: Counter = Counter()
    for chunk in chunks:
        active_roles = [
            role
            for role, score in chunk.role_scores.items()
            if score >= 0.05 or any(pattern_score(chunk.text, ROLE_PATTERNS[role][1])[1])
        ]
        route_ids = active_keys(chunk.routes, limit=len(ROUTES), threshold=0.2)
        fiber_ids = active_keys(chunk.fibers, limit=len(FIBERS), threshold=0.2)
        evidence = {
            "chunk_id": chunk.chunk_id,
            "document": chunk.title,
            "path": chunk.path,
            "excerpt": clean_text(chunk.text, 180),
        }
        for left, right, relationship, explanation in STRUCTURAL_ROLE_EDGES:
            if left in active_roles and right in active_roles:
                role_coactivation[(left, right, relationship)] += 1
                add_edge(
                    f"role:{left}",
                    f"role:{right}",
                    relationship,
                    weight=1.0 + chunk.role_scores.get(left, 0.0) + chunk.role_scores.get(right, 0.0),
                    vector_id="mechanism_role_order",
                    source_file=chunk.path,
                    evidence={**evidence, "explanation": explanation},
                )
        for route_id in route_ids:
            route_node = f"route:{route_id}"
            for role in active_roles:
                add_edge(route_node, f"role:{role}", "attends_to_role", weight=chunk.routes[route_id], vector_id=route_id, source_file=chunk.path, evidence=evidence)
        for fiber_id in fiber_ids:
            fiber_node = f"fiber:{fiber_id}"
            for role in active_roles:
                add_edge(fiber_node, f"role:{role}", "realizes_role", weight=chunk.fibers[fiber_id], vector_id=fiber_id, source_file=chunk.path, evidence=evidence)

        important = [term for term, _ in Counter(tokens(chunk.text)).most_common(8) if f"concept:{slugify(term)}" in nodes]
        for role in active_roles:
            for term in important[:6]:
                add_edge(f"role:{role}", f"concept:{slugify(term)}", "mentions_concept", weight=0.25, source_file=chunk.path, evidence=evidence)
        for index, left in enumerate(important[:6]):
            for right in important[index + 1 : 6]:
                concept_coactivation[tuple(sorted((left, right)))] += 1

    for (left, right), count in concept_coactivation.most_common(80):
        add_edge(f"concept:{slugify(left)}", f"concept:{slugify(right)}", "co_occurs", weight=float(count))

    node_rows = [node_row(node) for node in nodes.values()]
    edge_rows = [edge_row(edge) for edge in edge_weights.values()]
    protected = {node["id"] for node in node_rows if node["id"].startswith(("role:", "route:", "fiber:"))}
    edge_rows = sorted(
        edge_rows,
        key=lambda edge: (edge["relationship"] in {"writes_state", "constrains_state", "drives_readout", "gates_readout", "tested_by", "falsified_by"}, edge["weight"]),
        reverse=True,
    )[:max_edges]
    referenced = {edge["source"] for edge in edge_rows} | {edge["target"] for edge in edge_rows} | protected
    node_rows = [node for node in node_rows if node["id"] in referenced][:max_nodes]
    referenced = {node["id"] for node in node_rows}
    edge_rows = [edge for edge in edge_rows if edge["source"] in referenced and edge["target"] in referenced]
    statistics = compute_graph_metrics(node_rows, edge_rows)
    gaps = detect_receptor_gaps(node_rows, edge_rows, role_terms)

    return {
        "schema_version": 1,
        "model": "fieldbridge_sparse_attention_kg",
        "role_terms": role_terms,
        "nodes": node_rows,
        "edges": edge_rows,
        "statistics": statistics,
        "gap_report": gaps,
    }


def build_pdf_field_pack(
    pdf_dir: Path,
    field_id: str,
    label: str,
    out_dir: Path,
    description: str = "",
    max_docs: int = 200,
    max_chunks_per_doc: int = 40,
    max_chars: int = 2600,
    max_anchors: int = 80,
    extensions: Sequence[str] = (".pdf", ".txt", ".tex", ".md"),
) -> Dict[str, Any]:
    paths = iter_documents(pdf_dir, extensions)[:max_docs]
    if not paths:
        raise FileNotFoundError(f"No input documents found in {pdf_dir}")
    extraction_failures: List[Dict[str, str]] = []
    chunks = make_chunks(
        paths,
        max_chars=max_chars,
        max_chunks_per_doc=max_chunks_per_doc,
        failures=extraction_failures,
    )
    if not chunks:
        detail = extraction_failures[0]["error"] if extraction_failures else "No text found."
        raise RuntimeError(f"No readable chunks extracted from {pdf_dir}: {detail}")
    successful_path_names = {chunk.path for chunk in chunks}
    successful_paths = [path for path in paths if str(path) in successful_path_names]

    selected = top_chunks(chunks, score_chunk, max_anchors)
    pack = field_pack_from_attention(field_id, label, description, selected)
    records = [mechanism_record_from_chunk(field_id, chunk, index) for index, chunk in enumerate(selected)]
    route_priors = average_scores(selected, "routes", list(ROUTES.keys()))
    fiber_priors = average_scores(selected, "fibers", list(FIBERS.keys()))
    constructor_priors = average_scores(selected, "constructor_scores", list(CONSTRUCTOR_ROLE_PATTERNS.keys()))
    substrate_priors = average_scores(selected, "substrate_scores", list(SUBSTRATE_PATTERNS.keys()))
    attention = sparse_attention(chunks, top_k=10)
    graph = receptor_graph(selected)
    generated_at = datetime.now(timezone.utc).isoformat()
    adapter = build_field_adapter(
        field_id=field_id,
        label=label,
        description=description,
        chunks=chunks,
        selected=selected,
        paths=successful_paths,
        generated_at=generated_at,
    )
    adapter_report = render_adapter_markdown(adapter)

    evidence = {
        "schema_version": 1,
        "field_id": field_id,
        "label": label,
        "generated_at": generated_at,
        "source_artifacts": {
            "pdf_dir": str(pdf_dir),
            "documents": [str(path) for path in successful_paths],
            "extraction_failures": extraction_failures,
        },
        "claim_boundary": (
            "This sidecar is built by deterministic sparse attention over a field "
            "document folder. "
            "It identifies receptor vocabulary and mechanism anchors; it does not validate "
            "that the proposed mechanisms are experimentally true."
        ),
        "route_priors": route_priors,
        "fiber_priors": fiber_priors,
        "constructor_role_priors": constructor_priors,
        "substrate_priors": substrate_priors,
        "mechanism_anchors": records[:24],
        "sparse_attention": attention,
        "field_adapter_path": f"field_adapters/{field_id}.json",
        "field_adapter_report_path": f"reports/{field_id}_adapter.md",
        "receptor_graph": graph,
        "receptor_graph_path": f"kg/{field_id}_knowledge_graph.json",
        "gap_report": graph.get("gap_report", {}),
        "telemetry": {
            "document_count": len(successful_paths),
            "documents_discovered": len(paths),
            "documents_processed": len(successful_paths),
            "documents_failed": len(extraction_failures),
            "chunk_count": len(chunks),
            "anchor_count": len(records),
            "max_docs": max_docs,
            "max_chunks_per_doc": max_chunks_per_doc,
            "max_chars": max_chars,
        },
    }

    write_json(out_dir / "field_packs" / f"{field_id}.json", pack)
    write_json(out_dir / "field_pack_evidence" / f"{field_id}.json", evidence)
    write_json(out_dir / "field_adapters" / f"{field_id}.json", adapter)
    write_text(out_dir / "reports" / f"{field_id}_adapter.md", adapter_report)
    write_json(out_dir / "kg" / f"{field_id}_knowledge_graph.json", graph)
    write_json(out_dir / "index" / "core_examples.json", records)
    write_json(out_dir / "manifest.json", {
        "schema_version": 1,
        "generated_at": generated_at,
        "source": "FieldBridge PDF sparse-attention export",
        "fields": [
            {
                "field_id": field_id,
                "field_pack": f"field_packs/{field_id}.json",
                "evidence": f"field_pack_evidence/{field_id}.json",
                "field_adapter": f"field_adapters/{field_id}.json",
                "field_adapter_report": f"reports/{field_id}_adapter.md",
                "knowledge_graph": f"kg/{field_id}_knowledge_graph.json",
                "native_mechanism_anchors": len(records),
                "document_count": len(successful_paths),
                "documents_discovered": len(paths),
                "documents_failed": len(extraction_failures),
                "chunk_count": len(chunks),
            }
        ],
        "index": {
            "core_examples": "index/core_examples.json",
            "hyperion_static_index": None,
        },
    })
    return {
        "out_dir": str(out_dir),
        "field_id": field_id,
        "documents": len(successful_paths),
        "documents_discovered": len(paths),
        "documents_failed": len(extraction_failures),
        "chunks": len(chunks),
        "anchors": len(records),
        "route_priors": route_priors,
        "fiber_priors": fiber_priors,
        "constructor_role_priors": constructor_priors,
        "substrate_priors": substrate_priors,
        "field_adapter": str(out_dir / "field_adapters" / f"{field_id}.json"),
        "field_adapter_report": str(out_dir / "reports" / f"{field_id}_adapter.md"),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build a FieldBridge field pack from a folder of PDFs/texts.")
    parser.add_argument("pdf_dir", help="Folder containing PDFs, .txt, .tex, or .md files.")
    parser.add_argument("--field-id", required=True, help="Stable field id, e.g. material_intelligence.")
    parser.add_argument("--label", default="", help="Human label. Defaults to title-cased field id.")
    parser.add_argument("--description", default="", help="Human description for the generated field pack.")
    parser.add_argument("--out-dir", default="build/fieldbridge_pdf_export", help="Output data tree.")
    parser.add_argument("--max-docs", type=int, default=200)
    parser.add_argument("--max-chunks-per-doc", type=int, default=40)
    parser.add_argument("--max-chars", type=int, default=2600)
    parser.add_argument("--max-anchors", type=int, default=80)
    parser.add_argument("--extensions", default=".pdf,.txt,.tex,.md")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    field_id = slugify(args.field_id)
    label = args.label or field_id.replace("_", " ").title()
    extensions = [item.strip().lower() for item in args.extensions.split(",") if item.strip()]
    summary = build_pdf_field_pack(
        pdf_dir=Path(args.pdf_dir),
        field_id=field_id,
        label=label,
        description=args.description,
        out_dir=Path(args.out_dir),
        max_docs=args.max_docs,
        max_chunks_per_doc=args.max_chunks_per_doc,
        max_chars=args.max_chars,
        max_anchors=args.max_anchors,
        extensions=extensions,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
