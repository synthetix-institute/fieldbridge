"""Build the web page: the realizations of the language, the single-component edits that connect them, and their
mechanisms.

site_data.py computes every realization with FieldBridge and checks every edit; the page (fieldbridge/web) shows the
realization expression ((Omega, Xi); C, R, P; A), lets the visitor change one component at a time, and recomputes the
dynamics in the browser with engines that the tests compare with the Python calculations.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import time
from pathlib import Path
from typing import Callable

from . import site_data

ASSETS = Path(__file__).with_name("web")
ROOT = ASSETS.parents[1]
REPO = "https://github.com/synthetix-institute/fieldbridge"
SCRIPTS = ("model-physics.js", "unitary.js", "dissipative.js", "fields.js", "mathml.js", "views.js",
           "atlas.js", "expression.js", "site.js")
STYLES = ("site.css",)


def _gallery_cards(out: Path) -> dict:
    """The full memory card of each specification in examples/gallery, keyed by its name."""
    gallery = ROOT / "examples" / "gallery"
    index = json.loads((gallery / "gallery.json").read_text(encoding="utf-8"))
    dest = out / "gallery"
    dest.mkdir(parents=True, exist_ok=True)
    cards = {}
    for card in index["cards"]:
        png = gallery / f"{card['id']}.png"
        if png.exists():
            shutil.copy2(png, dest / png.name)
            cards[card["name"]] = {"image": f"gallery/{png.name}", "verdict": card["verdict"]}
    return cards


def publish_assets(out_dir: str | Path) -> Path:
    """Copy the page and its scripts next to an existing data.js; each file is versioned by the hash of its content,
    so that a browser does not keep an older script."""
    out = Path(out_dir)
    (out / "site").mkdir(parents=True, exist_ok=True)
    page = (ASSETS / "index.html").read_text(encoding="utf-8")
    for name in SCRIPTS + STYLES + ("data.js",):
        src = out / "site" / name if name == "data.js" else ASSETS / name
        if name != "data.js":
            shutil.copy2(src, out / "site" / name)
        version = hashlib.sha256(src.read_bytes()).hexdigest()[:12]
        page = page.replace(f'"site/{name}"', f'"site/{name}?v={version}"')
    (out / "index.html").write_text(page, encoding="utf-8")
    return out / "index.html"


def build_site(out_dir: str | Path, law: bool = False, save_law_record: bool = False,
               log: Callable[..., None] = print) -> Path:
    out = Path(out_dir)
    (out / "site").mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    data = site_data.build(law=law, log=log)
    if save_law_record:
        if not law:
            raise ValueError("the law record is written from a build with --law")
        log(f"law constants written to {site_data.write_law_record(data)}")
    data["repo"] = REPO
    cards = _gallery_cards(out)
    for rec in data["nodes"].values():
        if rec.get("spec") and not rec.get("derived"):
            name = json.loads((ROOT / rec["spec"]).read_text(encoding="utf-8")).get("name")
            if name in cards:
                rec["card"] = cards[name]
    site_data.write(data, out / "site" / "data.js")
    publish_assets(out)
    summary = {k: v for k, v in data.items() if k != "nodes"}
    summary["nodes"] = {k: {kk: vv for kk, vv in v.items() if kk != "engine"} for k, v in data["nodes"].items()}
    text = json.dumps(summary, indent=1, default=float, ensure_ascii=False) + "\n"
    (out / "site.json").write_text(text, encoding="utf-8")
    log(f"page with {len(data['nodes'])} realizations and {len(data['edges'])} edits in {time.time() - t0:.0f} s: "
        f"{out / 'index.html'}")
    return out / "index.html"
