"""Provenance and writing of reports for the families of mechanisms added after memory and quantum.

A report records the hash of the package that computed it, the versions of the numerical libraries, the hash of the
specification it read, and the boundary of what the calculation shows. The memory and quantum packages keep their own
copies of this code: their implementation hashes are recorded in the law record and in published reports.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Dict, Optional


def implementation_hash(package_dir: Path) -> str:
    return hashlib.sha256(b"".join(p.read_bytes() for p in sorted(Path(package_dir).glob("*.py")))).hexdigest()


def provenance(package_dir: Path, boundary: str, spec: Optional[Dict] = None,
               libraries=("numpy", "scipy", "sympy", "matplotlib")) -> Dict[str, object]:
    versions = {"python": platform.python_version()}
    for mod in libraries:
        try:
            versions[mod] = __import__(mod).__version__
        except ImportError:
            versions[mod] = None
    out = {"implementation_sha256": implementation_hash(package_dir), "versions": versions,
           "novelty_established": False, "evidence_boundary": boundary}
    if spec is not None:
        out.update(question=spec.get("question"),
                   input_sha256=hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest(),
                   assumptions=spec.get("assumptions"), provenance=spec.get("provenance"))
    return out


def write_report(out_dir: Path, name: str, report: Dict, markdown: str) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{name}.json").write_text(json.dumps(report, indent=2, default=float) + "\n", encoding="utf-8")
    (out_dir / f"{name}.md").write_text(markdown, encoding="utf-8")
