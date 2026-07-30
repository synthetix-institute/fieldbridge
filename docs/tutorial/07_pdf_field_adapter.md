# Chapter 7: Build a Field Adapter from PDFs

A field adapter turns a local paper collection into the target side of a
mechanism transfer. It identifies which carriers, operators, update laws,
closures, readouts, protocols, and falsifiers are supported by that collection.

## Install PDF Support

```bash
python3 -m pip install -e '.[pdf]'
```

## Build the Adapter

```bash
fieldbridge build-field-adapter /path/to/papers \
  --field-id active_matter \
  --label "Active Matter" \
  --out-dir build/active_matter \
  --max-docs 300 \
  --max-chunks-per-doc 40 \
  --max-anchors 120
```

The folder is searched recursively for `.pdf`, `.txt`, `.tex`, and `.md`
documents. Text-layer PDFs are parsed directly. Scanned PDFs require OCR first.
Unreadable files are listed in
`field_pack_evidence/active_matter.json` under
`source_artifacts.extraction_failures`.

## Read the Outputs

Start with:

```text
build/active_matter/reports/active_matter_adapter.md
```

Then inspect:

```text
field_adapters/active_matter.json
field_pack_evidence/active_matter.json
kg/active_matter_knowledge_graph.json
index/core_examples.json
```

The adapter implements the public identity:

```text
M = (Omega, Xi)
I_op = (M; C, R, P)
```

`operator_apparatus` and `update_or_transport` provide evidence for `Omega`;
`state_or_carrier` and the substrate profile provide evidence for `Xi`;
`admissibility_logic`, `readout_rule`, and `protocol_execution` provide
`C`, `R`, and `P`. Falsifiers remain a separate validation layer.

## Use the New Field

```bash
fieldbridge translate examples/bioelectric_regeneration.txt \
  --to active_matter \
  --data-dir build/active_matter

fieldbridge construct examples/bioelectric_regeneration.txt \
  --to active_matter \
  --data-dir build/active_matter
```

The generated adapter is an evidence index, not an accepted field model.
Inspect source passages and equations before promoting a transfer.
