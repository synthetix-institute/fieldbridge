# FieldBridge

FieldBridge extracts mechanisms from scientific papers and finds structurally
analogous equations in other fields.

It is built for one public task:

```text
paper A -> extracted mechanism -> analogous mechanism in field B -> testable translation
```

The comparison is not based on topic words alone. FieldBridge uses a small,
transparent operational fingerprint: transport, closure, spectral/operator
structure, boundary realization, incompatibility, discrete protocol, and a few
realization fibers.

The first public version is intentionally transparent. It does not require a trained
model or a private database. It ships with three starter field packs:

- material intelligence
- biological intelligence
- collective intelligence

FieldBridge can be used as:

- a mechanism extractor for one paper or LaTeX fragment;
- a paper-to-paper analogy checker;
- an equation analogy finder across fields;
- a mechanism translator from one field language into another;
- a field-pack format for community contributions;
- a small public workbench that can later connect to a larger private or
  institutional fingerprint database.

## Install

```bash
python -m pip install -e .
```

## Quick Start

```bash
fieldbridge fields
fieldbridge extract examples/bioelectric_regeneration.txt
fieldbridge compare examples/bioelectric_regeneration.txt examples/material_memory.tex
fieldbridge fingerprint examples/bioelectric_regeneration.txt
fieldbridge search examples/bioelectric_regeneration.txt --target-field material_intelligence
fieldbridge translate examples/bioelectric_regeneration.txt --to material_intelligence
```

## What It Does

FieldBridge reads a text or LaTeX fragment and returns a mechanism sheet:

- state;
- input or writing condition;
- boundary, interface, or context;
- output;
- candidate equations;
- measurements;
- falsifying controls;
- active route/fiber structure.

It can then compare that mechanism with another paper or search a field pack for
analogous mechanisms.

Example:

```text
Source field: biological intelligence
Target field: material intelligence

Invariant:
An input writes an internal state; after the input is removed, the state changes
a later response through a closure law and a boundary condition.

Material analogue:
A thermal, electrical, chemical, optical, or mechanical writing step imprints a
retained material state. Later transport, conductance, release, shape, or motion
must be predicted by that state and fail after erasure or shuffled-history controls.
```

## What It Does Not Claim

FieldBridge is not a proof engine and does not validate a physical claim by itself.
The public starter version produces testable hypotheses. A mechanism becomes real
only when its state variable, equations, measurements, and falsifying controls are
confirmed in the target system.

FieldBridge also does not expose a full private parser or a full large-corpus
embedding system. The starter repository uses a deterministic, inspectable
fingerprint so the community can understand, test, and extend the method without
needing private data.

## Optional Full Database

The public repo is deliberately small. Larger installations can attach:

- high-dimensional fingerprints;
- vector indexes over equations or patents;
- paper-level witnesses and references;
- validated analogy records;
- private field packs.

The command-line interface should remain the same. A stronger database should make
the results better, not change the public contract.

## Community Field Packs

A field pack adds the vocabulary, variables, observables, controls, and example
mechanism records for a field. See `docs/DATA_MODEL.md`.

Useful starter contributions:

- neuroscience and bioelectricity;
- synthetic biology;
- soft robotics;
- chemical reaction networks;
- active matter;
- swarm robotics;
- economics and institutions;
- patents and invention claims.

## License

MIT.
