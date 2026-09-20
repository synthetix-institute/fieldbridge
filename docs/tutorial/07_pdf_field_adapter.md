# Build a searchable field from a paper collection

A field pack supplies the records that retrieval can return. Building one from
papers makes the source passages inspectable and allows a new field to enter
the workbench without changing the search algorithm.

For a small offline rehearsal, use the bundled text and TeX examples:

```bash
python3 -B -m fieldbridge build-field-adapter examples \
  --field-id tutorial_field --label "Tutorial field" \
  --extensions .txt,.tex --max-docs 4 --max-anchors 12 \
  --out-dir build/tutorial_field

python3 -B -m fieldbridge --data-dir build/tutorial_field \
  fields
python3 -B -m fieldbridge --data-dir build/tutorial_field \
  search examples/brownian_probability_flow.tex \
  --target-field tutorial_field --top-k 3
```

The global `--data-dir` option belongs **before** the subcommand. Keeping it
there avoids accidentally searching the default dataset.

## Inspect what was extracted

Start with `build/tutorial_field/reports/tutorial_field_adapter.md`.
The same output tree contains:

| Path under the output directory | Purpose |
| --- | --- |
| `field_packs/tutorial_field.json` | Target field description |
| `field_adapters/tutorial_field.json` | Recognized carrier and operation cues |
| `field_pack_evidence/tutorial_field.json` | Source artifacts and extraction failures |
| `index/core_examples.json` | Records used by retrieval |
| `kg/tutorial_field_knowledge_graph.json` | Relations between recognized roles |

[build_pdf_field_pack](../../fieldbridge/pdf_sparse_builder.py) discovers files,
extracts text, makes bounded passages, scores evidence and writes these
coordinated outputs. Its `sparse_attention` function is rule-based passage
scoring, not a trained attention model that proves relations between equations.

## Move to PDFs

```bash
python3 -m pip install -e '.[pdf]'
python3 -B -m fieldbridge build-field-adapter /path/to/papers \
  --field-id active_matter --label "Active Matter" \
  --max-docs 300 --max-chunks-per-doc 40 --max-anchors 120 \
  --out-dir build/active_matter
```

Scanned pages need OCR. An extracted PDF equation may lose indices, fractions
or signs. Compare a candidate passage with its original display and nearby
definitions before using it as a mathematical input. The evidence file lists
unreadable files under `source_artifacts.extraction_failures`.

## Connect a source to a calculation

The field builder generates retrieval records; it does not automatically
supply the typed `calculation_source` annotation used by the exact verifier.
Add that annotation only after checking the equation, symbols and domain.
Then supply a proposed map or observable using the
[retrieved-source walkthrough](13_retrieval_to_calculation.md).

**Exercise.** Select one output record and trace its equation back to the input
document. Identify one physical assumption that would be lost if only the
displayed formula, rather than its surrounding text, were retained.

[Next: first calculated construction](08_end_to_end_walkthrough.md) · [Tutorial](index.md)
