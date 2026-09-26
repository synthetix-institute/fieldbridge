"""Figures for the chapter on the language of mechanisms (chapter 24).

Run the chapter's commands first; they write their results under build/tut. Then

    python3 -B docs/tutorial/figures/quantum/make_figures.py [--build build/tut]

copies the co-discovery and attachment figures here. No calculation is repeated.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
COPIES = {"q_codiscover/codiscover.png": "q_codiscovery.png", "q_attach_chain/attach.png": "q_attach_chain.png"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--build", default="build/tut", help="directory with the chapter's results")
    build = Path(ap.parse_args().build)
    for src, dst in COPIES.items():
        shutil.copyfile(build / src, HERE / dst)
    print("figures written to", HERE)


if __name__ == "__main__":
    main()
