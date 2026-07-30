# Field Adapters From PDF Corpora

FieldBridge can build a field adapter from a folder of papers. The adapter is a
reviewable bridge between the public Hyperion route grammar and a field's own
experimental language.

The target use case is a folder with tens to hundreds of PDFs representing a
field or subfield: neuroscience, active matter, chemical reaction networks,
bioelectricity, materials ageing, soft robotics, collective behaviour, or a
focused disease/mechanism literature.

The adapter answers one practical question:

```text
If a portable mechanism arrives from another field, which state carriers,
operators, update laws, admissibility tests, readouts, execution protocols,
substrates, and falsifiers are available in this field corpus?
```

It is not a static glossary. It is a receptor map for mechanism transfer.

## Conceptual Contract

FieldBridge separates three layers that are usually mixed in prose:

1. **Public Hyperion route grammar.**  
   The transferable operation type: transport, closure, spectral/operator
   structure, boundary realization, incompatibility, or discrete protocol.

2. **Constructor roles.**  
   The field-native parts required to make a mechanism testable: carrier,
   operator, update, admissibility, readout, protocol, and falsifier.

3. **Universal substrate evidence.**  
   The type of space on which the construction acts: coordinate domain,
   metric manifold, inner-product space, phase space, probability space, graph
   topology, lattice/site space, bundle/gauge space, configuration quotient, or
   stoichiometric space.

Field-specific nouns remain evidence phrases. They are not treated as final
ontology. For example, in a neuroscience corpus, *synapse*, *connectome*,
*cortical state*, and *belief state* may appear as evidence. The adapter asks
whether those phrases instantiate a smaller substrate class such as graph
topology, probability space, coordinate domain, or inner-product space.

## Installation

```bash
cd /path/to/fieldbridge
python -m pip install -e ".[pdf]"
```

The optional `pdf` extra installs `pdfplumber`. If it is absent, the builder
tries `pypdf`, `PyPDF2`, and the command-line `pdftotext` fallback when
available.

## Input Contract

The input path is traversed recursively. The builder reads text-layer PDFs and
plain `.txt`, `.tex`, and `.md` documents. A corrupt or unreadable document is
recorded under `source_artifacts.extraction_failures` and the rest of the folder
continues. The command fails only when no readable document produces a text
chunk.

The builder does not perform optical character recognition. Run OCR before
ingestion when a PDF contains page images without an embedded text layer. A
quick preflight is:

```bash
pdftotext first_paper.pdf - | head
```

The same parser is used by `extract`, `compare`, `translate`, `construct`,
complete-paper validation, and the folder builder, so a PDF has one extraction
contract across the CLI.

`fieldbridge extract paper.pdf` produces one coarse mechanism sheet for the
whole document. Use `build-field-adapter` for collections and for papers that
contain several mechanisms: the adapter splits each document into source-linked
chunks before scoring roles, routes, and substrates.

## Basic Command

```bash
fieldbridge build-field-adapter /path/to/papers \
  --field-id neuroscience \
  --label "Neuroscience" \
  --description "Mechanism adapter built from a neuroscience PDF corpus." \
  --out-dir build/neuroscience_adapter \
  --max-docs 500 \
  --max-chunks-per-doc 40 \
  --max-anchors 160
```

The command accepts `.pdf`, `.txt`, `.tex`, and `.md` files by default. Use
`--extensions` to restrict the corpus:

```bash
fieldbridge build-field-adapter /path/to/papers \
  --field-id chemical_reaction_networks \
  --label "Chemical Reaction Networks" \
  --extensions .pdf,.txt \
  --out-dir build/chemical_reaction_networks
```

`build-field-pack` remains as a compatibility alias. The stronger output is the
field adapter.

## Output Tree

```text
build/neuroscience_adapter/
  field_packs/neuroscience.json
  field_adapters/neuroscience.json
  reports/neuroscience_adapter.md
  field_pack_evidence/neuroscience.json
  kg/neuroscience_knowledge_graph.json
  index/core_examples.json
  manifest.json
```

### `field_adapters/<field_id>.json`

The main adapter contract. It contains route profiles, constructor role profiles,
substrate profiles, evidence snippets, and route-to-field mappings.

### `reports/<field_id>_adapter.md`

A human-readable audit report. Use it first when deciding whether the field
corpus is mature enough for transfer.

### `field_packs/<field_id>.json`

A compact field pack that can be used by the existing `translate` command.

### `field_pack_evidence/<field_id>.json`

The detailed evidence sidecar: sparse-attention heads, priors, anchors, adapter
paths, and KG paths.

### `kg/<field_id>_knowledge_graph.json`

A typed mechanism graph. It links role nodes, route nodes, field phrases, and
concept nodes through relationships such as `writes_state`, `constrains_state`,
`drives_readout`, and `tested_by`.

### `index/core_examples.json`

Native mechanism anchors extracted from the corpus. These are searchable
examples for analog search and target-field translation.

## Adapter Schema

The adapter JSON has this top-level structure:

```json
{
  "schema_version": 1,
  "artifact_type": "fieldbridge_field_adapter",
  "field_id": "neuroscience",
  "label": "Neuroscience",
  "corpus": {
    "documents": 500,
    "chunks": 14320,
    "anchor_chunks": 160
  },
  "adapter_contract": {
    "input": "mechanism sheet or route/fiber fingerprint",
    "output": "field-native state/carrier, operator, update, admissibility, readout, protocol, and falsifier candidates"
  },
  "route_profile": {},
  "fiber_profile": {},
  "constructor_role_profile": {},
  "substrate_profile": {},
  "constructor_roles": [],
  "universal_substrates": [],
  "route_to_field_adapter": [],
  "field_native_receivers": {},
  "gap_report": {}
}
```

Scores are deterministic pattern evidence over selected mechanism-bearing
chunks. They are comparative signals, not calibrated probabilities.

## Constructor Roles

| Role | Question Answered | Typical Evidence |
| --- | --- | --- |
| `state_or_carrier` | What carries the mechanism? | state, field, density, distribution, parameter, potential, memory, configuration |
| `operator_apparatus` | What transforms the carrier? | operator, generator, kernel, Hamiltonian, Laplacian, force, objective, message passing |
| `update_or_transport` | How does it change? | flow, diffusion, propagation, inference, recurrence, learning, trajectory |
| `admissibility_logic` | What makes a state or step legal? | constraint, closure, normalization, threshold, gate, compatibility, residual |
| `readout_rule` | What is measured or predicted? | observable, measurement, response, spectrum, eigenvalue, phenotype, behaviour |
| `protocol_execution` | How is the mechanism prepared or executed? | preparation, intervention, sequence, workflow, update order, training schedule |
| `falsifier` | What can break the claim? | control, validation, ablation, reset, washout, baseline, null, perturbation |

A strong transfer candidate should have evidence for all seven roles. A field can
still be useful with weaker evidence, but the missing roles should become the
next curation targets.

## Universal Substrate Classes

| Substrate | Meaning |
| --- | --- |
| `coordinate_domain` | Physical or spatiotemporal domain: position, region, surface, interface, `x,t`, `Omega` |
| `metric_manifold` | Metric, manifold, curvature, geodesic, differential-geometric substrate |
| `inner_product_space` | Vector, Hilbert, basis, norm, projection, eigenvector, inner product |
| `phase_space` | Hamiltonian, symplectic, canonical, trajectory, position-momentum substrate |
| `probability_space` | Probability, distribution, stochastic process, Bayesian update, entropy, expectation |
| `graph_topology` | Graph, network, adjacency, node, edge, topology, connectivity |
| `lattice_site_space` | Lattice, grid, site, chain, spin, nearest-neighbour substrate |
| `bundle_gauge_space` | Bundle, gauge, connection, covariant derivative, local-frame substrate |
| `configuration_quotient` | Configuration, quotient, orbit, symmetry-reduced, moduli, equivalence classes |
| `stoichiometric_space` | Chemical species, stoichiometric matrix, reaction network, concentration vector |

These classes are deliberately few. The aim is to collapse local field nouns
into transferable substrates rather than multiply synonyms.

## How Sparse Attention Is Used

The builder performs four deterministic passes:

1. **Document extraction.**  
   PDFs and text files are converted to plain text. The text is split into
   bounded chunks.

2. **Route and fiber scoring.**  
   Each chunk is scored against the public FieldBridge route/fiber patterns.

3. **Constructor and substrate scoring.**  
   Each chunk is scored against constructor-role and universal-substrate
   patterns.

4. **Anchor selection and graph construction.**  
   High-scoring chunks become native mechanism anchors. Role co-activation and
   route-role co-activation build the typed mechanism graph.

The result is transparent enough to inspect: every active route, role, and
substrate has evidence snippets attached.

## Reading The Adapter Report

Start with `reports/<field_id>_adapter.md`.

Review these sections in order:

1. **Active Hyperion Routes.**  
   This tells which operation types are well represented in the corpus.

2. **Universal Substrate Evidence.**  
   This tells what kind of mathematical or physical space the field tends to use.

3. **Constructor Role Evidence.**  
   This tells whether the corpus has carriers, operators, updates,
   admissibility tests, readouts, protocols, and falsifiers.

4. **Route To Field Adapter.**  
   This is the operational mapping from Hyperion route to field-native evidence.

5. **Gaps.**  
   This gives the curation tasks: missing constructor roles or weak substrate
   evidence.

## Review Criteria For A Transfer

A transferred mechanism is ready for serious discussion when the adapter supplies:

- a named field-native carrier or admissible state space;
- an operator, generator, update rule, or objective acting on that carrier;
- an admissibility condition, boundary, normalization, threshold, or closure law;
- a readout that can be measured or computed;
- an intervention, preparation, or computational protocol that executes it;
- a falsifier or control that can remove, alter, or reduce the claimed effect;
- at least one substrate class that explains where the operator acts.

For experimental fields, the falsifier is essential. A translation without a
control is only a mechanism sketch.

## Example: Neuroscience Folder

Run:

```bash
fieldbridge build-field-adapter ~/corpora/neuroscience_pdfs \
  --field-id neuroscience \
  --label "Neuroscience" \
  --out-dir build/neuroscience_adapter \
  --max-docs 800 \
  --max-anchors 200
```

Expected useful signals:

- `probability_space` when the corpus emphasizes belief states, Bayesian
  inference, uncertainty, stochastic dynamics, or variational objectives;
- `graph_topology` when the corpus emphasizes connectivity, networks,
  adjacency, connectomes, or graph Laplacians;
- `coordinate_domain` when the corpus emphasizes fields over tissue, cortical
  sheets, spatial gradients, or anatomical domains;
- `inner_product_space` when the corpus emphasizes basis expansions, spectra,
  modes, projections, or linear state spaces.

The local terms may be synapse, neuron, cortical state, population code,
connectome, or policy. Those terms should be inspected as evidence phrases. The
adapter-level claim concerns the substrate and constructor roles.

## Example: Chemistry Folder

For a chemical reaction corpus, useful signals often include:

- `stoichiometric_space` for species vectors, stoichiometric matrices, reaction
  networks, mass-action laws, and concentration vectors;
- `probability_space` for stochastic reaction networks and master equations;
- `coordinate_domain` for reaction-diffusion or transport in physical space;
- `constraint_closure_route` for conservation, detailed balance, steady-state
  constraints, or thermodynamic feasibility.

This makes chemical corpora useful for translating mechanisms where the portable
object is not a molecule name but a constrained update over a species-state
substrate.

## Example: Materials Folder

For a material-intelligence corpus, useful signals often include:

- `coordinate_domain` for fields, boundaries, interfaces, and physical samples;
- `metric_manifold` when curvature, shape, deformation, or geometry is active;
- `phase_space` for dynamical systems and trajectories;
- `lattice_site_space` for spin, lattice, defect, or grid models;
- `transport_flow_route` and `boundary_weak_form_route` for flux, diffusion,
  interfaces, and weak-form descriptions.

This is the natural adapter target for turning a mechanism from a biological or
collective system into a material experiment.

## Using An Adapter In A Transfer Workflow

1. Build the adapter from a target-field corpus.
2. Inspect `reports/<field_id>_adapter.md`.
3. Run `fieldbridge extract` on the source paper or mechanism fragment.
4. Run `fieldbridge translate --to <field_id>` using the generated data directory.
5. Compare the translation against the adapter report.
6. Promote the translation only when carrier, operator/update, admissibility,
   readout, protocol, falsifier, and substrate evidence are all present.

Example:

```bash
fieldbridge translate examples/bioelectric_regeneration.txt \
  --to neuroscience \
  --data-dir build/neuroscience_adapter
```

## Corpus Design

The adapter quality depends on corpus design.

Recommended corpus types:

- a coherent field or subfield;
- papers with equations, experimental protocols, or explicit models;
- review papers plus primary mechanism papers;
- negative controls and validation papers when available;
- enough breadth to cover carriers, operators, readouts, protocols, and falsifiers.

Avoid mixing unrelated domains in the same adapter. If the folder contains
neuroscience, immunology, materials, and economics together, the substrate
profile will be less interpretable. Build separate adapters and compare them.

## Parameter Guidance

| Parameter | Use |
| --- | --- |
| `--max-docs` | Number of documents to read. Increase for large fields. |
| `--max-chunks-per-doc` | Upper bound on chunks per paper. Increase for long review papers. |
| `--max-chars` | Chunk size. Larger chunks preserve context; smaller chunks improve locality. |
| `--max-anchors` | Number of mechanism anchors exported. Increase for broad fields. |
| `--extensions` | File types included in the corpus. |

For hundreds of PDFs, start with:

```bash
--max-docs 500 --max-chunks-per-doc 40 --max-anchors 160
```

For a smaller focused corpus:

```bash
--max-docs 80 --max-chunks-per-doc 50 --max-anchors 80
```

## Quality Checks

After building an adapter, check:

- Do the top routes match the field's actual mechanism style?
- Are the top substrate classes few and interpretable?
- Are constructor roles supported by different documents rather than one repeated
  source?
- Do evidence snippets contain equations, protocols, or measurements?
- Are falsifiers present, or does the corpus mostly describe positive effects?
- Are the mechanism anchors readable enough to support translation?

Weak adapters are still useful. They tell which papers or measurements are
missing from the field corpus.

## Relation To Hyperion

The public FieldBridge adapter uses deterministic patterns and small route/fiber
fingerprints. A larger Hyperion installation can replace the scoring layer with
high-dimensional fingerprints, equation witnesses, Lagrangian roads, Noether
checks, or Gromov-Wasserstein neighborhoods.

The adapter contract should remain stable:

```text
portable route/fingerprint
  -> field-native carrier/operator/update/admissibility/readout/protocol/falsifier
  -> substrate evidence
  -> evidence snippets
  -> transfer review
```

This keeps the public repo usable while allowing private or institutional
databases to improve the evidence layer.

## Practical Interpretation

The strongest output is not a high score by itself. The strongest output is a
field-native transfer recipe:

```text
carrier:       what stores or carries the mechanism
operator:      what acts on it
update:        how it changes
admissibility: what makes the state legal
readout:       what is measured
protocol:      how the mechanism is prepared or executed
falsifier:     what breaks the proposed mechanism
substrate:     where the operator acts
```

That recipe is what allows FieldBridge to move from analogy to experimental
design.
