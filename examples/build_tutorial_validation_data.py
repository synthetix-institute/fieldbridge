"""Build small, deterministic manifests for the FieldBridge tutorial.

The records exercise the public validation contracts. They are synthetic
software examples, not scientific performance evidence.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "build" / "tutorial_data"


def build_retrieval_manifest() -> Path:
    paper_dir = OUT / "papers"
    paper_dir.mkdir(parents=True, exist_ok=True)
    papers = [
        ("continuum-diffusion", "continuum", "diffusion", "partial_t u = kappa nabla^2 u. No-flux closure conserves mass."),
        ("graph-diffusion", "graphs", "diffusion", "dot u = -kappa L_G u. Graph Laplacian closure conserves total state."),
        ("mechanical-oscillation", "mechanics", "oscillation", "partial_t^2 q + omega^2 q = 0. Periodic motion conserves energy."),
        ("network-oscillation", "networks", "oscillation", "dot z = A z. Complex eigenmodes produce network oscillation."),
    ]
    manifest = []
    for paper_id, field_id, mechanism_id, sentence in papers:
        path = paper_dir / f"{paper_id}.txt"
        path.write_text((sentence + "\n\n") * 30, encoding="utf-8")
        manifest.append({
            "paper_id": paper_id,
            "path": str(path.relative_to(OUT)),
            "field_id": field_id,
            "mechanism_id": mechanism_id,
        })
    path = OUT / "papers.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return path


def build_continuation_manifest() -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    # demo-00 and demo-01 fall in the default evaluation fold; the remaining
    # complete paper groups fit the conditional tables.
    for paper_index in range(10):
        for step in range(6):
            rows.append({
                "paper_id": f"demo-{paper_index:02d}",
                "current_omega": f"Omega{step % 3:02d}",
                "current_xi": f"Xi{step % 2:02d}",
                "current_completion": "C+R" if step % 2 else "C",
                # Keep T1 deliberately uninformative so the tutorial shows the
                # extra information supplied by the complete current state.
                "first_move": "T00",
                "next_move": f"T{(step + 1) % 3:02d}",
                "destination_omega": f"Omega{step % 3:02d}",
                "destination_xi": f"Xi{(step + 1) % 2:02d}",
                "destination_completion": "C+R+P" if step % 2 else "C+R",
            })
    path = OUT / "transitions.jsonl"
    path.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")
    return path


if __name__ == "__main__":
    retrieval = build_retrieval_manifest()
    continuation = build_continuation_manifest()
    print(f"retrieval manifest: {retrieval}")
    print(f"continuation manifest: {continuation}")
