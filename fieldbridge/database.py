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
        data = json.loads(path.read_text())
        packs.append(FieldPack(**data))
    return packs


def load_records(data_dir: Path | None = None) -> List[MechanismRecord]:
    root = data_dir or default_data_dir()
    path = root / "index" / "core_examples.json"
    rows = json.loads(path.read_text())
    return [mechanism_record_from_row(row) for row in rows]


def mechanism_record_from_row(row: Dict[str, Any]) -> MechanismRecord:
    payload = {key: value for key, value in row.items() if key in RECORD_FIELDS}
    return MechanismRecord(**payload)


def load_static_index_records(data_dir: Path | None = None, max_records: int | None = None) -> List[MechanismRecord]:
    root = data_dir or default_data_dir()
    path = root / "index" / "hyperion_static_index.json"
    if not path.exists():
        return []
    rows = (json.loads(path.read_text()).get("records") or [])
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
        row = json.loads(path.read_text())
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
