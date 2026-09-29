"""The web page (python -m fieldbridge demo): realizations of the language, the single-component edits that connect
them, and their mechanisms. See web_demo.py and site_data.py."""
from __future__ import annotations

import argparse
from pathlib import Path

from .site_data import LAW_RECORD


def run(out_dir, law: bool = False, save_law_record: bool = False, log=print) -> Path:
    from .web_demo import build_site
    return build_site(out_dir, law=law, save_law_record=save_law_record, log=log)


def add_parser(sub) -> None:
    p = sub.add_parser("demo", help="Build the web page: realizations, the edits that connect them, their mechanisms.")
    p.add_argument("--out-dir", default="build/site", help="Output directory (default build/site)")
    p.add_argument("--law", action="store_true",
                   help="Compute the law constants (about fifteen minutes); otherwise they are read from "
                        f"{Path(LAW_RECORD).relative_to(Path(LAW_RECORD).parents[2])}")
    p.add_argument("--save-law-record", action="store_true", help="With --law, write the law constants to the record")
    p.add_argument("--studio-only", action="store_true", help=argparse.SUPPRESS)  # accepted from older workflows
    p.add_argument("--quiet", action="store_true", help="Print only the path of the page")

    def execute(args):
        page = run(args.out_dir, law=args.law, save_law_record=args.save_law_record,
                   log=(lambda *a: None) if args.quiet else print)
        print(page)
        return 0
    p.set_defaults(func=execute)
