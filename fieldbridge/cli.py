from __future__ import annotations

import argparse
from pathlib import Path

from .database import load_field_packs
from .constructor import construct_transfer
from .continuation import render_markdown as render_continuation_markdown
from .continuation import validate_future_state
from .extract import compare_mechanisms, extract_mechanism
from .pdf_sparse_builder import build_pdf_field_pack, read_document, slugify
from .render import render_comparison, render_constructor, render_fingerprint, render_mechanism_sheet, render_search, render_translation
from .routes import fingerprint_text
from .search import find_analogs, translate_mechanism
from .zero_shot import render_markdown as render_zero_shot_markdown
from .zero_shot import validate_full_paper_zero_shot


def read_input(path_or_text: str) -> str:
    path = Path(path_or_text)
    if path.is_file():
        return read_document(path)
    if path.exists():
        raise IsADirectoryError(f"Expected a document, received directory: {path}")
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


def cmd_construct(args: argparse.Namespace) -> int:
    import json
    from dataclasses import asdict

    calculate = args.calculate
    if calculate and (not args.correspondence or not args.out_dir):
        raise SystemExit("--calculate requires --correspondence and --out-dir")
    if not calculate and (args.correspondence or args.out_dir):
        raise SystemExit("--correspondence and --out-dir require --calculate")
    request = None
    destination = Path(args.out_dir) if calculate else None
    if calculate:
        destination.mkdir(parents=True, exist_ok=True)
        (destination / "run_status.json").write_text(json.dumps({"status": "running"}) + "\n", encoding="utf-8")
        try:
            request = json.loads(Path(args.correspondence).read_text(encoding="utf-8"))
            if not isinstance(request, dict):
                raise ValueError("Correspondence must be a JSON object")
        except (OSError, ValueError) as error:
            (destination / "run_status.json").write_text(json.dumps({"status": "failed", "reason": str(error)}) + "\n", encoding="utf-8")
            raise SystemExit(f"Invalid correspondence: {error}") from error
    try:
        text = read_input(args.input)
        transfer = construct_transfer(
            text,
            target_field=args.to,
            data_dir=Path(args.data_dir) if args.data_dir else None,
            top_k=args.top_k,
            include_hyperion=not args.no_hyperion,
            calculation_request=request,
        )
    except Exception as error:
        if calculate:
            (destination / "run_status.json").write_text(json.dumps({"status": "failed", "reason": str(error)}) + "\n", encoding="utf-8")
        raise
    rendered = render_constructor(transfer)
    if calculate:
        (destination / "transfer.json").write_text(json.dumps(asdict(transfer), indent=2) + "\n", encoding="utf-8")
        (destination / "transfer.md").write_text(rendered + "\n", encoding="utf-8")
        outcome = transfer.calculation
        # Replace stale success artifacts with null/refusal on unsuccessful reruns.
        (destination / "construction_spec.json").write_text(json.dumps(outcome.get("specification"), indent=2) + "\n", encoding="utf-8")
        (destination / "calculation.json").write_text(json.dumps(outcome.get("report", outcome), indent=2) + "\n", encoding="utf-8")
        (destination / "run_status.json").write_text(json.dumps({"status": outcome["status"]}) + "\n", encoding="utf-8")
    print(rendered)
    return 2 if calculate and transfer.calculation["status"] != "calculated" else 0


def cmd_mechanisms(args: argparse.Namespace) -> int:
    from .site_references import write_doc
    print(write_doc(args.out))
    return 0


def cmd_verify_construction(args: argparse.Namespace) -> int:
    import json
    try:
        from .verification import ConstructionError, run_file
    except ModuleNotFoundError as error:
        if error.name != "sympy":
            raise
        raise SystemExit("Install the calculation extra: pip install -e '.[construction]'") from error
    try:
        report = run_file(args.input, args.out_dir)
    except (ConstructionError, KeyError, TypeError, json.JSONDecodeError) as error:
        raise SystemExit(f"Construction input rejected: {error}") from error
    print(json.dumps(report, indent=2))
    return 0


def cmd_design_spin_cancellation(args: argparse.Namespace) -> int:
    import json
    try:
        from .inverse_spin import run_file
    except ModuleNotFoundError as error:
        if error.name not in {"sympy", "numpy"}:
            raise
        raise SystemExit("Install the calculation extra: pip install -e '.[construction]'") from error
    try:
        report = run_file(args.input, args.out_dir)
    except (ValueError, KeyError, TypeError, OSError) as error:
        raise SystemExit(f"Spin design rejected: {error}") from error
    print(json.dumps({"status": report["status"], "design": report["construction"],
                      "nearest_neighbor": report["nearest_neighbor_comparison"],
                      "prediction": report["prediction"], "checks": report["checks"]}, indent=2))
    return 0


def cmd_build_field_pack(args: argparse.Namespace) -> int:
    field_id = slugify(args.field_id)
    label = args.label or field_id.replace("_", " ").title()
    extensions = [item.strip().lower() for item in args.extensions.split(",") if item.strip()]
    summary = build_pdf_field_pack(
        pdf_dir=Path(args.pdf_dir),
        field_id=field_id,
        label=label,
        description=args.description,
        out_dir=Path(args.out_dir),
        max_docs=args.max_docs,
        max_chunks_per_doc=args.max_chunks_per_doc,
        max_chars=args.max_chars,
        max_anchors=args.max_anchors,
        extensions=extensions,
    )
    import json

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


def cmd_validate_zero_shot(args: argparse.Namespace) -> int:
    import json

    report = validate_full_paper_zero_shot(
        Path(args.manifest),
        top_k=args.top_k,
        max_chars=args.max_chars,
        bootstrap_samples=args.bootstrap_samples,
        seed=args.seed,
        min_eligible_queries=args.min_eligible_queries,
    )
    out_json = Path(args.out_json)
    out_md = Path(args.out_md)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(render_zero_shot_markdown(report), encoding="utf-8")
    print(json.dumps({
        "json": str(out_json),
        "markdown": str(out_md),
        "readiness": report["readiness"],
        "eligible_queries": report["eligible_queries"],
    }, indent=2))
    return 0


def cmd_validate_continuation(args: argparse.Namespace) -> int:
    import json

    report = validate_future_state(
        Path(args.manifest),
        seed=args.seed,
        evaluation_fold=args.evaluation_fold,
        alpha=args.alpha,
        min_evaluation_transitions=args.min_evaluation_transitions,
    )
    out_json = Path(args.out_json)
    out_md = Path(args.out_md)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(render_continuation_markdown(report), encoding="utf-8")
    print(json.dumps({
        "json": str(out_json),
        "markdown": str(out_md),
        "readiness": report["readiness"],
        "predictive_gate": report["predictive_gate"],
        "next_move_accuracy": report["targets"]["next_move"]["accuracy"],
    }, indent=2))
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

    construct = sub.add_parser(
        "construct",
        help="Propose a mechanism transfer or calculate a supported annotated source.",
    )
    construct.add_argument("input", help="Source paper, equation fragment, or literal text.")
    construct.add_argument("--to", required=True, help="Target field id.")
    construct.add_argument("--top-k", type=int, default=4)
    construct.add_argument("--no-hyperion", action="store_true", help="Use only public field-pack receptors.")
    construct.add_argument("--calculate", action="store_true", help="Calculate a supported retrieved model using an explicit correspondence.")
    construct.add_argument("--correspondence", help="fieldbridge-correspondence/1 JSON selecting a retrieved source and a map or observable.")
    construct.add_argument("--out-dir", help="Calculation output directory; required with --calculate.")
    construct.set_defaults(func=cmd_construct)

    verify = sub.add_parser("verify-construction", help="Derive and verify a supplied physical model; does not retrieve or infer source equations.")
    verify.add_argument("input", help="Construction JSON with equations, map and assumptions.")
    verify.add_argument("--out-dir", required=True, help="Directory for the calculation and its provenance.")
    verify.set_defaults(func=cmd_verify_construction)

    mechanisms = sub.add_parser("mechanisms", help="Write the definitions, derivations and original publications of the "
                                                   "mechanisms on the web page (docs/mechanisms.md).")
    mechanisms.add_argument("--out", default="docs/mechanisms.md")
    mechanisms.set_defaults(func=cmd_mechanisms)

    try:  # the memory constructor needs the memory extra; the rest of the program does not
        from .memory.cli import add_parser as add_memory_parser
        add_memory_parser(sub)
    except ImportError:
        pass

    try:  # the quantum language needs numpy and sympy (the construction extra)
        from .quantum.cli import add_parser as add_quantum_parser
        add_quantum_parser(sub)
    except ImportError:
        pass

    try:  # the regulation module reads its bodies with the memory parser (the memory extra)
        from .regulation.cli import add_parser as add_regulation_parser
        add_regulation_parser(sub)
    except ImportError:
        pass

    try:  # the demo calculates with the memory and quantum packages (the memory extra)
        from .demo import add_parser as add_demo_parser
        add_demo_parser(sub)
    except ImportError:
        pass

    design = sub.add_parser("design-spin-cancellation", help="Solve pairwise spin-interaction constraints and test the polarization.")
    design.add_argument("input", help="fieldbridge-spin-design/1 JSON with exchange bonds and collective coupling.")
    design.add_argument("--out-dir", required=True, help="Directory for design, input and direct-dynamics comparison.")
    design.set_defaults(func=cmd_design_spin_cancellation)

    def add_build_pack_parser(name: str, help_text: str) -> None:
        build_pack = sub.add_parser(name, help=help_text)
        build_pack.add_argument("pdf_dir", help="Folder containing PDFs, .txt, .tex, or .md files.")
        build_pack.add_argument("--field-id", required=True, help="Stable field id, e.g. material_intelligence.")
        build_pack.add_argument("--label", default="", help="Human label. Defaults to title-cased field id.")
        build_pack.add_argument("--description", default="", help="Human description for the generated field pack.")
        build_pack.add_argument("--out-dir", default="build/fieldbridge_pdf_export", help="Output data tree.")
        build_pack.add_argument("--max-docs", type=int, default=200)
        build_pack.add_argument("--max-chunks-per-doc", type=int, default=40)
        build_pack.add_argument("--max-chars", type=int, default=2600)
        build_pack.add_argument("--max-anchors", type=int, default=80)
        build_pack.add_argument("--extensions", default=".pdf,.txt,.tex,.md")
        build_pack.set_defaults(func=cmd_build_field_pack)

    add_build_pack_parser(
        "build-field-pack",
        "Build a FieldBridge field pack, evidence sidecar, adapter, and mechanism KG from a folder of PDFs/texts.",
    )
    add_build_pack_parser(
        "build-field-adapter",
        "Build a field-native mechanism adapter from a folder of PDFs/texts.",
    )
    validate = sub.add_parser(
        "validate-zero-shot",
        help="Evaluate cross-field mechanism retrieval while holding out each complete paper.",
    )
    validate.add_argument("manifest", help="JSON/JSONL rows with paper_id, path, mechanism_id, and field_id.")
    validate.add_argument("--top-k", type=int, default=10)
    validate.add_argument("--max-chars", type=int, default=2600)
    validate.add_argument("--bootstrap-samples", type=int, default=2000)
    validate.add_argument("--seed", type=int, default=17)
    validate.add_argument("--min-eligible-queries", type=int, default=100)
    validate.add_argument("--out-json", default="build/full_paper_zero_shot.json")
    validate.add_argument("--out-md", default="build/full_paper_zero_shot.md")
    validate.set_defaults(func=cmd_validate_zero_shot)

    continuation = sub.add_parser(
        "validate-continuation",
        help="Predict the next move and future mechanism state in complete held-out papers.",
    )
    continuation.add_argument(
        "manifest",
        help=(
            "JSON/JSONL transition records with paper_id, current state, first_move, "
            "next_move, and destination state."
        ),
    )
    continuation.add_argument("--seed", type=int, default=20260710)
    continuation.add_argument("--evaluation-fold", type=int, default=9, choices=range(10))
    continuation.add_argument("--alpha", type=float, default=0.5)
    continuation.add_argument("--min-evaluation-transitions", type=int, default=100)
    continuation.add_argument("--out-json", default="build/future_state_validation.json")
    continuation.add_argument("--out-md", default="build/future_state_validation.md")
    continuation.set_defaults(func=cmd_validate_continuation)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
