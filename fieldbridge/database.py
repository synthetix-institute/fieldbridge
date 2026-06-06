from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, List, Tuple

from .models import FieldPack, MechanismRecord


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
    return [MechanismRecord(**row) for row in rows]


def load_all(data_dir: Path | None = None) -> Tuple[List[FieldPack], List[MechanismRecord]]:
    return load_field_packs(data_dir), load_records(data_dir)


def field_by_id(field_packs: Iterable[FieldPack], field_id: str) -> FieldPack:
    for pack in field_packs:
        if pack.field_id == field_id:
            return pack
    known = ", ".join(pack.field_id for pack in field_packs)
    raise KeyError(f"Unknown field '{field_id}'. Known fields: {known}")
