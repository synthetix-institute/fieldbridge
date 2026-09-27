"""Build the offline constructor studio without running discovery sweeps."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

ASSETS = Path(__file__).with_name("web")
ROOT = ASSETS.parents[1]
REPO = "https://github.com/synthetix-institute/fieldbridge"

SOURCES = {
    "ito": {
        "title": "Stochastic Integral", "authors": "K. Ito", "year": 1944,
        "journal": "Proceedings of the Imperial Academy 20, 519–524",
        "url": "https://doi.org/10.3792/pia/1195572786",
        "supports": "Stochastic integration underlying the change-of-variable calculation. The numerical parameters and transformations below are worked examples.",
    },
    "toggle": {
        "title": "Construction of a genetic toggle switch in Escherichia coli",
        "authors": "T. S. Gardner, C. R. Cantor and J. J. Collins", "year": 2000,
        "journal": "Nature 403, 339–342", "url": "https://doi.org/10.1038/35002131",
        "supports": "An experimentally constructed bistable gene-regulatory network. The scalar normal form in this demo is a separate model; a reduction from the two-concentration equations requires its own derivation.",
    },
    "kramers": {
        "title": "Brownian motion in a field of force and the diffusion model of chemical reactions",
        "authors": "H. A. Kramers", "year": 1940, "journal": "Physica 7, 284–304",
        "url": "https://doi.org/10.1016/S0031-8914(40)90098-2",
        "supports": "Thermally activated escape over an energy barrier. The displayed escape time uses the overdamped, weak-noise approximation.",
    },
    "snap": {
        "title": "Amplifying the response of soft actuators by harnessing snap-through instabilities",
        "authors": "J. T. B. Overvelde, T. Kloek, J. J. A. D'haen and K. Bertoldi", "year": 2015,
        "journal": "PNAS 112, 10863–10868", "url": "https://doi.org/10.1073/pnas.1504947112",
        "supports": "Experimental use of snap-through instability in interconnected fluidic actuators. Retention and sensing performance of a proposed event latch require additional measurements.",
    },
}

def build_gallery(out: Path) -> list[dict]:
    """Publish the existing model calculations, without recomputing their simulations."""
    destination = out / "gallery"
    destination.mkdir(exist_ok=True)
    filenames = ("colloid_patch", "compartments", "dipole_patch", "laser", "pitchfork",
                 "repressilator", "repressor_ring4", "schlogl", "toggle", "toggle_unequal",
                 "tubes", "tubes_unequal")
    citations = {0: "https://doi.org/10.1103/PhysRevE.62.5263",
                 3: "https://doi.org/10.1103/RevModPhys.47.67",
                 5: "https://doi.org/10.1038/35002125",
                 7: "https://doi.org/10.1007/BF01379769",
                 8: SOURCES["toggle"]["url"], 9: SOURCES["toggle"]["url"],
                 10: "https://doi.org/10.1016/j.jtbi.2006.07.015",
                 11: "https://doi.org/10.1016/j.jtbi.2006.07.015"}
    cards = []
    for i, filename in enumerate(filenames):
        card_id = f"card{i:02d}"
        path = ROOT / "examples/gallery" / (card_id + ".json")
        card = json.loads(path.read_text(encoding="utf-8"))
        spec = json.loads((ROOT / "examples/memory" / (filename + ".json")).read_text(encoding="utf-8"))
        for suffix in (".json", ".png"):
            shutil.copy2(path.with_suffix(suffix), destination / (card_id + suffix))
        cards.append({"id": card_id, "name": card["name"], "slots": card["slots"],
                      "verdict": card["verdict"], "tags": card["tags"],
                      "params": card["params"], "provenance": spec["provenance"]["source"],
                      "source_url": citations.get(i),
                      "image": f"gallery/{card_id}.png", "record": f"gallery/{card_id}.json",
                      "specification": f"examples/memory/{filename}.json"})
    (out / "gallery.html").write_text(
        '<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0;url=index.html#gallery">'
        '<a href="index.html#gallery">Memory in model materials: a gallery of realizations</a>\n', encoding="utf-8")
    return cards

def build_quantum(out: Path) -> dict:
    from .quantum.language import load, attach, codiscover
    source = load(ROOT / "examples/quantum/two_spins.json")
    targets = ("qubit", "spin", "bosons", "chain", "fermion-pair", "correlated-pair", "collective")
    attachments = {target: attach(source, target, size=4 if target in ("bosons", "chain") else 3)
                   for target in targets}
    for target, report in attachments.items():
        if not report["attached"] or report["law"]["target_residual"] > 1e-8:
            raise ValueError(f"Failed quantum attachment: {target}")
    discovery = codiscover([load(p) for p in sorted((ROOT / "examples/quantum").glob("*.json"))])
    records = {"source": source.name, "attachments": attachments, "examples": discovery,
               "tutorial": REPO + "/blob/main/docs/tutorial/24_spin_language.md"}
    (out / "quantum_examples.json").write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    return records

def build_studio(out_dir: str | Path) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    destination = out / "studio"
    destination.mkdir(exist_ok=True)
    for name in ("studio.css", "collections.css", "workspace.css", "physics.js", "studio.js", "collections.js", "lucide.min.js", "LUCIDE_LICENSE"):
        shutil.copy2(ASSETS / name, destination / name)
    config = {"repo": REPO, "sources": SOURCES, "schema": "fieldbridge-browser-construction/1"}
    from .verification import verify_construction
    examples = Path(__file__).resolve().parents[1] / "examples" / "construction"
    verified = {}
    for name in ("ito_square", "log_signal_ito", "log_signal_stratonovich"):
        spec = json.loads((examples / (name + ".json")).read_text(encoding="utf-8"))
        report = verify_construction(spec)
        if report["residual_coefficients"] != ["0", "0"]:
            raise ValueError(f"Failed generator identity: {name}")
        verified[name] = {"specification": spec, "calculation": report,
                          "tutorial": REPO + "/blob/main/docs/tutorial/10_stochastic_construction.md",
                          "test_file": REPO + "/blob/main/tests/test_stochastic_conventions.py"}
    (out / "verified_examples.json").write_text(json.dumps(verified, indent=2) + "\n", encoding="utf-8")
    config["verified_examples"] = {name: item["calculation"] for name, item in verified.items()}
    config["gallery"] = build_gallery(out)
    config["quantum"] = build_quantum(out)
    template = (ASSETS / "index.html").read_text(encoding="utf-8")
    for asset in ASSETS.iterdir():
        if asset.suffix in (".js", ".css"):
            version = hashlib.sha256(asset.read_bytes()).hexdigest()[:12]
            template = template.replace(f'studio/{asset.name}"', f'studio/{asset.name}?v={version}"')
    template = template.replace("__CONFIG__", json.dumps(config, ensure_ascii=False).replace("<", "\\u003c"))
    page = out / "index.html"
    page.write_text(template, encoding="utf-8")
    return page
