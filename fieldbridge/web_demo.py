"""Build the offline constructor studio without running discovery sweeps."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

ASSETS = Path(__file__).with_name("web")
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

def build_studio(out_dir: str | Path) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    destination = out / "studio"
    destination.mkdir(exist_ok=True)
    for name in ("studio.css", "physics.js", "studio.js", "lucide.min.js", "LUCIDE_LICENSE"):
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
    template = (ASSETS / "index.html").read_text(encoding="utf-8")
    template = template.replace("__CONFIG__", json.dumps(config, ensure_ascii=False).replace("<", "\\u003c"))
    page = out / "index.html"
    page.write_text(template, encoding="utf-8")
    return page
