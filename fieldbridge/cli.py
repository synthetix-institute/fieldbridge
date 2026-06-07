from __future__ import annotations

import argparse
from pathlib import Path

from .database import load_field_packs
from .extract import compare_mechanisms, extract_mechanism
from .render import render_comparison, render_fingerprint, render_mechanism_sheet, render_search, render_translation
from .routes import fingerprint_text
from .search import find_analogs, translate_mechanism


def read_input(path_or_text: str) -> str:
    path = Path(path_or_text)
    if path.exists():
        return path.read_text(errors="replace")
    return path_or_text


def cmd_fields(args: argparse.Namespace) -> int:
    packs = load_field_packs(Path(args.data_dir) if args.data_dir else None)
    for pack in packs:
        print(f"{pack.field_id}\t{pack.label}\t{pack.description}")
    return 0


def cmd_fingerprint(args: argparse.Namespace) -> int:
    text = read_input(args.input)
    print(render_fingerprint(fingerprint_text(text)))
    return 0


def cmd_search(args: argparse.Namespace) -> int:
    text = read_input(args.input)
    matches = find_analogs(
        text,
        target_field=args.target_field,
        data_dir=Path(args.data_dir) if args.data_dir else None,
        top_k=args.top_k,
        include_hyperion=args.include_hyperion,
    )
    print(render_search(matches))
    return 0


def cmd_extract(args: argparse.Namespace) -> int:
    text = read_input(args.input)
    title = args.title or (Path(args.input).name if Path(args.input).exists() else "input")
    print(render_mechanism_sheet(extract_mechanism(text, title=title)))
    return 0


def cmd_compare(args: argparse.Namespace) -> int:
    source = read_input(args.source)
    target = read_input(args.target)
    print(render_comparison(compare_mechanisms(source, target, source_title=args.source, target_title=args.target)))
    return 0


def cmd_translate(args: argparse.Namespace) -> int:
    text = read_input(args.input)
    translation = translate_mechanism(
        text,
        target_field=args.to,
        data_dir=Path(args.data_dir) if args.data_dir else None,
        top_k=args.top_k,
        include_hyperion=not args.no_hyperion,
    )
    print(render_translation(translation))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fieldbridge",
        description=(
            "No-LLM mechanism translation across fields: extract a portable "
            "mechanism, find analogues, or render it in a new field formulation."
        ),
    )
    parser.add_argument("--data-dir", default="", help="Optional data directory with field_packs/ and index/.")
    sub = parser.add_subparsers(required=True)

    fields = sub.add_parser("fields", help="List available field packs.")
    fields.set_defaults(func=cmd_fields)

    fingerprint = sub.add_parser("fingerprint", help="Compute a route/fiber fingerprint.")
    fingerprint.add_argument("input", help="Path or literal text.")
    fingerprint.set_defaults(func=cmd_fingerprint)

    search = sub.add_parser("search", help="Find existing analogous mechanisms in field packs.")
    search.add_argument("input", help="Path or literal text.")
    search.add_argument("--target-field", default=None, help="Restrict to a field id.")
    search.add_argument("--top-k", type=int, default=5)
    search.add_argument("--include-hyperion", action="store_true", help="Include Hyperion static witness records as audit evidence.")
    search.set_defaults(func=cmd_search)

    extract = sub.add_parser("extract", help="Extract a mechanism sheet from one paper or fragment.")
    extract.add_argument("input", help="Path or literal text.")
    extract.add_argument("--title", default="", help="Optional title shown in the rendered sheet.")
    extract.set_defaults(func=cmd_extract)

    compare = sub.add_parser("compare", help="Compare mechanisms extracted from two papers or fragments.")
    compare.add_argument("source", help="Source path or literal text.")
    compare.add_argument("target", help="Target path or literal text.")
    compare.set_defaults(func=cmd_compare)

    translate = sub.add_parser(
        "translate",
        help="Translate an extracted mechanism into a target field formulation.",
    )
    translate.add_argument("input", help="Path or literal text.")
    translate.add_argument("--to", required=True, help="Target field id.")
    translate.add_argument("--top-k", type=int, default=4)
    translate.add_argument("--no-hyperion", action="store_true", help="Use only field-pack anchors, without Hyperion witness evidence.")
    translate.set_defaults(func=cmd_translate)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
