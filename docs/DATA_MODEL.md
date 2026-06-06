# FieldBridge Data Model

FieldBridge has three public data objects: mechanism sheets, field packs, and
mechanism records.

## Mechanism Sheet

A mechanism sheet is extracted from one paper or fragment. It is the intermediate
object used by `extract`, `compare`, and `translate`.

```text
state q
input u
boundary B
output y
candidate equations
measurements
falsifying controls
route/fiber fingerprint
```

This is intentionally smaller than a full paper parser. The public tool only needs
enough structure to compare mechanisms and produce testable translations.

## Field Pack

A field pack names how a field talks about the same operational roles.

```json
{
  "field_id": "material_intelligence",
  "label": "Material intelligence",
  "description": "...",
  "state_words": ["retained material state m(x,t)", "..."],
  "input_words": ["thermal writing protocol", "..."],
  "boundary_words": ["sample boundary B", "..."],
  "output_words": ["transport", "..."],
  "validation_words": ["no-write baseline", "..."]
}
```

Use field packs to add a new scientific community without changing the engine.
For example, a neuroscience pack might call the state a membrane potential,
synaptic trace, latent policy, or cortical state; a patent pack might call it
an internal technical state.

## Mechanism Record

A mechanism record is a concrete analogy target. It has equations, variables,
measurements, controls, references, keywords, and an explicit route/fiber
profile.

```json
{
  "record_id": "mat_droplet_interface_001",
  "title": "Active droplet interfacial memory",
  "field_id": "material_intelligence",
  "summary": "...",
  "invariant": "...",
  "equations": ["partial_t Gamma = ..."],
  "variables": ["Gamma(s,t): interfacial coverage state"],
  "measurements": ["surface-tension map"],
  "controls": ["fresh droplet in written bath"],
  "references": ["doi/arxiv/community reference"],
  "keywords": ["droplet", "marangoni", "interface"],
  "routes": {
    "transport_flow": 0.96,
    "constraint_closure": 0.82,
    "spectral_operator": 0.16,
    "boundary_weak_form": 0.94,
    "commutator_incompatibility": 0.34,
    "discrete_protocol": 0.20
  },
  "fibers": {
    "structure": 0.84,
    "spectral": 0.14,
    "geometry": 0.82,
    "syntax": 0.60,
    "entropy": 0.18
  }
}
```

## Optional Full Fingerprint Database

The starter repo uses route/fiber records. A larger installation can add:

- high-dimensional fingerprint vectors;
- nearest-neighbor vector indexes;
- equation witnesses and arXiv identifiers;
- validated analogy records;
- learned embeddings or other coordinates.

The public rule is simple: if the database is not present, FieldBridge must say
that it is returning a deterministic starter fingerprint, not a full embedding.
