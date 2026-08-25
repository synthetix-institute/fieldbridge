from __future__ import annotations

import json
from dataclasses import fields
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

from .models import FieldPack, MechanismRecord


RECORD_FIELDS = {field.name for field in fields(MechanismRecord)}


def default_data_dir() -> Path:
    repo_data = Path.cwd() / "data"
    if repo_data.exists():
        return repo_data
    return Path(__file__).resolve().parent.parent / "data"


def load_field_packs(data_dir: Path | None = None) -> List[FieldPack]:
    root = data_dir or default_data_dir()
    packs: List[FieldPack] = []
    for path in sorted((root / "field_packs").glob("*.json")):
        data = json.loads(path.read_text(encoding="utf-8"))
        packs.append(FieldPack(**data))
    return packs


def load_records(data_dir: Path | None = None) -> List[MechanismRecord]:
    root = data_dir or default_data_dir()
    path = root / "index" / "core_examples.json"
    rows = json.loads(path.read_text(encoding="utf-8"))
    return [mechanism_record_from_row(row) for row in rows]


BAD_EQUATION_MARKERS = (
    "rotatebox",
    "begintabular",
    "__HYPERION",
    "bfseries",
    "multicolumn",
    "includegraphics",
)


def _is_hyperion_row(row: Dict[str, Any]) -> bool:
    return row.get("field_id") == "hyperion_equation" or row.get("source") == "hyperion_equation_witnesses"


def _clean_equation(value: Any) -> str:
    text = str(value or "").replace("\\n", " ").replace("\n", " ")
    text = " ".join(text.split())
    return text[:420].strip()


def _safe_public_equation(value: Any) -> bool:
    text = _clean_equation(value)
    lower = text.lower()
    if len(text) < 8 or len(text) > 420:
        return False
    if any(marker.lower() in lower for marker in BAD_EQUATION_MARKERS):
        return False
    if lower.count("begin") > 1 or lower.count("end") > 1:
        return False
    if text.count("{") + text.count("}") > 80:
        return False
    if not any(marker in text for marker in ("=", "partial", "\\partial", "nabla", "\\nabla", "lambda", "\\lambda", "int", "\\int")):
        return False
    return True


def _sanitize_equations(row: Dict[str, Any]) -> List[str]:
    equations = row.get("equations") or []
    if not isinstance(equations, list):
        return []
    if _is_hyperion_row(row):
        return []
    cleaned = [_clean_equation(item) for item in equations if _safe_public_equation(item)]
    return cleaned or [_clean_equation(item) for item in equations if _clean_equation(item)]


def _sanitize_summary(row: Dict[str, Any]) -> str:
    summary = str(row.get("summary") or "")
    if _is_hyperion_row(row) and "Clean source witness:" in summary:
        return summary.split("Clean source witness:", 1)[0].strip()
    return summary


def mechanism_record_from_row(row: Dict[str, Any]) -> MechanismRecord:
    payload = {key: value for key, value in row.items() if key in RECORD_FIELDS}
    payload["equations"] = _sanitize_equations(row)
    payload["summary"] = _sanitize_summary(row)
    return MechanismRecord(**payload)


def load_static_index_records(data_dir: Path | None = None, max_records: int | None = None) -> List[MechanismRecord]:
    root = data_dir or default_data_dir()
    path = root / "index" / "hyperion_static_index.json"
    if not path.exists():
        return []
    rows = (json.loads(path.read_text(encoding="utf-8")).get("records") or [])
    if max_records is not None:
        rows = rows[:max_records]
    return [mechanism_record_from_row(row) for row in rows]


def load_field_evidence(data_dir: Path | None = None) -> Dict[str, Dict[str, Any]]:
    root = data_dir or default_data_dir()
    evidence_dir = root / "field_pack_evidence"
    if not evidence_dir.exists():
        return {}
    evidence: Dict[str, Dict[str, Any]] = {}
    for path in sorted(evidence_dir.glob("*.json")):
        row = json.loads(path.read_text(encoding="utf-8"))
        field_id = row.get("field_id") or path.stem
        evidence[field_id] = row
    return evidence


def load_all(
    data_dir: Path | None = None,
    *,
    include_static_index: bool = False,
    max_static_records: int | None = None,
) -> Tuple[List[FieldPack], List[MechanismRecord]]:
    records = load_records(data_dir)
    if include_static_index:
        records = records + load_static_index_records(data_dir, max_records=max_static_records)
    return load_field_packs(data_dir), records


def field_by_id(field_packs: Iterable[FieldPack], field_id: str) -> FieldPack:
    for pack in field_packs:
        if pack.field_id == field_id:
            return pack
    known = ", ".join(pack.field_id for pack in field_packs)
    raise KeyError(f"Unknown field '{field_id}'. Known fields: {known}")
