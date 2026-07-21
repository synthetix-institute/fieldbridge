# FieldBridge Data Model

FieldBridge has five public data objects: mechanism sheets, constructor
transfers, field packs, field adapters, and mechanism records.

## Constructor Transfer

A constructor transfer records the proposed identity change rather than only a
similarity score:

```text
I_op=(Omega, Xi; C, R, P) -> I'_op=(Omega, Xi'; C', R', P')
preserved contract
changed carrier and realization
required C/R/P attachments
candidate equations and predictions
falsifying controls and validation gates
```

The public tool reports interpretable route/fiber proxies for `Omega` and `Xi`.
It does not claim access to private learned token assignments. A structurally
complete constructor record is a testable proposal, not a validated transfer.

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

## Field Adapter

A field adapter is generated from a folder of PDFs or text files. It is the
reviewable layer between the public Hyperion route grammar and a field's own
experimental language.

See [`FIELD_ADAPTERS.md`](FIELD_ADAPTERS.md) for the full corpus workflow,
review criteria, and interpretation guide.

```json
{
  "artifact_type": "fieldbridge_field_adapter",
  "field_id": "neuroscience",
  "route_profile": {"spectral_operator_route": 0.31},
  "constructor_role_profile": {
    "state_or_carrier": 0.44,
    "operator_apparatus": 0.28,
    "admissibility_logic": 0.24
  },
  "substrate_profile": {
    "coordinate_domain": 0.22,
    "probability_space": 0.19,
    "graph_topology": 0.15
  },
  "field_native_receivers": {
    "state_or_carrier": ["...evidence phrase..."],
    "operator_apparatus": ["...evidence phrase..."],
    "readout_rule": ["...evidence phrase..."]
  }
}
```

The adapter keeps two layers separate:

- constructor roles: carrier, operator, update, admissibility, readout, falsifier;
- substrate evidence: coordinate domain, metric manifold, inner-product space,
  phase space, probability space, graph topology, lattice, bundle/gauge space,
  configuration quotient, or stoichiometric space.

Field-specific nouns are retained as evidence phrases, not as final substrate
classes. This prevents a neuroscience, biology, or chemistry corpus from
fragmenting into many local labels when the relevant transferable substrate may
be graph topology, probability space, coordinate domain, or stoichiometric
space.

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
