# FieldBridge

FieldBridge is a public companion to the Hyperion project at
[Synthetix Institute](https://synthetix.institute). It is motivated by the
Hyperion atlas result: mathematical-scientific texts are extremely diverse in
their objects and terminology, but their mechanisms collapse to a small
operational grammar of transport, closure, spectral/operator structure, boundary
realization, incompatibility, and discrete protocol. In that view, field-specific
nouns are often the dressing; the portable unit is the mechanism.

**The invariant is the transformation; discovery is a successful translation.**

FieldBridge does not look for papers that talk about similar topics. It converts
a paper into an operational mechanism language first, and only then asks whether
another field contains the same route of action.

## Guided Tutorial

The [FieldBridge tutorial](docs/tutorial/index.md) follows the codebase-
knowledge format developed by PocketFlow: it identifies the core abstractions,
maps their relationships, and introduces them in dependency order through
runnable examples. The six chapters cover fingerprints, mechanism sheets,
cross-field retrieval, constructor transfers, complete-paper validation and
future-state prediction.

That changes the comparison. A bioelectric memory, a material hysteresis, and a
collective trace may use unrelated words. They become comparable only after they
are written as a mechanism: a state is written, constrained by a boundary or
context, transported or transformed, and tested by a later output and a residual
control.

The important question is not *which papers sound similar?* It is:

> Do these systems implement the same operational route, and how would that route
> be rendered in a new field?

![FieldBridge workflow](docs/assets/fieldbridge-workflows.svg)

The public version exposes one main route, with two supporting routes:

```text
main route:       paper A -> extracted identity -> constructor transfer -> target tests
support route 1:  paper A -> extracted mechanism -> translated mechanism in field B
support route 2:  paper A -> extracted mechanism -> analog search in field B
```

The first public version is intentionally transparent. It does not require a trained
model or a private database. It ships with four starter field packs:

- material intelligence
- biological intelligence
- collective intelligence
- stochastic optimization

FieldBridge can be used as:

- a mechanism extractor for one paper or LaTeX fragment;
- a paper-to-paper analogy checker;
- an equation analogy finder that retrieves existing cross-field examples;
- a mechanism translator that renders a source mechanism in a new field's
  variables, equations, measurements, and falsifying controls;
- a constructor that states what is preserved, which substrate changes, which
  closure/readout/protocol clauses must be attached, and how the transfer fails;
- a field-pack format for community contributions;
- a small public workbench that can later connect to a larger private or
  institutional fingerprint database.

## The Hyperion Result Behind It

The website atlas makes the central point visible: scientific knowledge does not
only organize around nouns such as *cell*, *droplet*, *agent*, *field*, or
*particle*. Those nouns are local embodiments. The reusable structure is the
route by which a system transports, closes, constrains, resolves spectra, meets a
boundary, fails to commute, or updates by protocol.

That is the reason translation is the main product. If the route is universal,
then a mechanism can be moved into a new substrate. Analogy search is useful, but
it is secondary: it tells us where the route has already appeared. Translation
asks where the route could appear next.

## The New Language

FieldBridge uses a small public version of the Hyperion language. The private
atlas has richer symbols and fingerprints; the public workbench keeps the core
idea inspectable:

```text
mechanism := state q + input u + boundary/context B + output y
route     := transport + closure + spectral/operator + boundary
test      := measurement + residual + falsifying control
```

In ordinary language, the nouns dominate: cell, droplet, particle, robot,
network, patent claim. In the operational language, those nouns become
realization choices. The central object is the route by which the mechanism acts.

This is why cross-field translation is possible at all. FieldBridge is not
claiming that two topics are similar. It is asking whether two different
scientific descriptions can be expressed by the same operational route, and then
rendering that route in a target field's variables, equations, measurements, and
controls.

## Why This Is Different

This is not another semantic-search app, RAG demo, or LLM prompt wrapper.

Most tools ask whether two papers use similar words. FieldBridge asks whether two
systems implement the same mechanism under different names, geometries, and
experimental cultures.

```text
ordinary search:       "find papers about regeneration"
LLM rewriting:         "explain regeneration in material language"
FieldBridge:           "extract the route, then build its material form"
```

The output is not a poetic analogy. A useful translation must name:

- the state variable that carries memory or structure;
- the input that writes or perturbs it;
- the boundary or context that makes the mechanism admissible;
- the equation or equation skeleton;
- the measurement that would reveal it;
- the control experiment that would make the claim fail.

## Hyperion Philosophy

Hyperion treats scientific knowledge as a system of transformations rather than a
taxonomy of objects. The same mechanism can appear as a bioelectric state, a
hydrogel memory, a collective trace, an interface law, or a patentable technical
state. The surface nouns change; the operational form may persist.

FieldBridge exposes a small public version of that idea. It does not try to
reproduce the full private Hyperion parser, atlas, or high-dimensional
fingerprint database. Instead, it provides a minimal open workflow that others can
inspect, extend, and test.

## How It Works

1. Extract a mechanism sheet from a paper or LaTeX fragment.
2. Score a transparent route/fiber fingerprint.
3. Search field packs for already-known analogues.
4. Translate the mechanism into a selected field language.
5. Return variables, equations, measurements, and falsifying controls.

The public engine is deterministic and does not require an LLM. A larger
installation can attach high-dimensional fingerprints, arXiv witnesses, vector
indexes, patent corpora, or private field packs, while keeping the same public
commands.

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
fieldbridge construct examples/brownian_probability_flow.tex --to stochastic_optimization
```

## Constructor: From A Match To A Test

The Hyperion language is hierarchical:

```text
mechanism core       M = (Omega, Xi)
completion fiber     F = (C, R, P)
operational identity I_op = (M; F)
realized model       I = (I_op; A)
```

`Omega` is the transformation apparatus, `Xi` its carrier, `C` closes the
construction, `R` makes it observable, `P` makes it executable, and `A` binds
it to field-specific objects, parameters, units and provenance.

`fieldbridge construct` turns this hierarchy into a practical transfer record.
It reports six operations separately:

1. preserve a qualified operator contract;
2. replace or attach a substrate;
3. supply closure or admissibility conditions;
4. define a measurable readout;
5. define an intervention or computational protocol;
6. break a retained clause with a falsifying control.

A transfer is mechanism-preserving only when the retained contract survives
the new closure and its consequences pass the controls. Occupancy of a coarse
operator--substrate pair is not enough.

The Brownian example demonstrates the public workflow:

```bash
fieldbridge construct examples/brownian_probability_flow.tex \
  --to stochastic_optimization --no-hyperion
```

It preserves gradient drift plus diffusion, replaces particle position by a
parameter-space carrier, and attaches normalization, stationary-current,
readout and negative-control obligations. This is a controlled reconstruction
of a known Langevin/Fokker--Planck transfer, not a claim that FieldBridge
discovered stochastic optimization.

## Full-Paper Zero-Shot Validation

FieldBridge includes a leave-one-paper-out benchmark for cross-field mechanism
retrieval. The query is an entire PDF, TeX, Markdown or text paper. Every chunk
from that paper is excluded from the retrieval gallery. A match counts as
relevant only when an independently labelled paper implements the same
mechanism in a different field.

The benchmark compares the public operational route/fiber representation with
a full-paper TF--IDF baseline fitted on the gallery after removing the query
paper. It reports top-1 accuracy, mean reciprocal rank, precision and recall at
`k`, and a paired bootstrap interval for the precision gain. A manifest row
has four required fields:

```json
{"paper_id":"paper-001","path":"papers/a.pdf","mechanism_id":"conservative_diffusion","field_id":"continuum_physics"}
```

Run the benchmark with:

```bash
fieldbridge validate-zero-shot benchmark/papers.jsonl \
  --top-k 10 \
  --min-eligible-queries 100 \
  --out-json build/full_paper_zero_shot.json \
  --out-md build/full_paper_zero_shot.md
```

This command establishes a full-paper zero-shot evaluation protocol. It does
not establish successful transfer until an independently labelled benchmark
contains enough eligible cross-field queries and the paired confidence interval
supports a gain over the lexical baseline.

## Predict The Next Mechanism State

Cross-field retrieval asks whether another paper contains the same mechanism.
Continuation prediction asks a stronger and different question: given the
current operational state and the first observed move, what move and mechanism
state will the next equation occupy?

Here `destination` means the state assigned to an unseen future equation. It is
not classification of the current equation. In the corpus-scale Hyperion
evaluation, complete papers were excluded from fitting. The model identified
the next move with 57.4% accuracy and predicted the future operator, substrate
and completion states with 74.4%, 71.7% and 73.4% accuracy, respectively.
FieldBridge records these results as companion evidence and provides a public
reproduction command for labelled transition manifests:

```bash
fieldbridge validate-continuation benchmark/transitions.jsonl \
  --min-evaluation-transitions 100 \
  --out-json build/future_state_validation.json \
  --out-md build/future_state_validation.md
```

The command reports next-move and destination-state performance separately,
together with gains over first-move-only baselines. See
[`docs/FUTURE_STATE_VALIDATION.md`](docs/FUTURE_STATE_VALIDATION.md) for the
schema, exact interpretation and corpus-scale reference record.

## Build A Field Pack From PDFs

FieldBridge can also build a small public field pack from a folder of papers. The
builder follows the older KG-generation pattern used in the concept-economics
prototype: parse PDFs, split them into mechanism-bearing chunks, run sparse
attention heads over route/fiber/role language, and export both a field pack and
a typed mechanism knowledge graph.

Full documentation: [`docs/FIELD_ADAPTERS.md`](docs/FIELD_ADAPTERS.md).

```bash
python -m pip install -e ".[pdf]"  # optional; otherwise pypdf/PyPDF2/pdftotext can be used

fieldbridge build-field-adapter /path/to/papers \
  --field-id material_intelligence \
  --label "Material Intelligence" \
  --out-dir build/material_intelligence_export \
  --max-docs 200 \
  --max-anchors 80
```

The export contains:

- `field_packs/<field_id>.json`: target-field vocabulary for translation;
- `field_adapters/<field_id>.json`: route-to-field adapter with constructor
  roles, universal substrate evidence, and review gates;
- `reports/<field_id>_adapter.md`: a human-readable adapter report;
- `field_pack_evidence/<field_id>.json`: sparse-attention heads, anchors, priors,
  and evidence boundaries;
- `kg/<field_id>_knowledge_graph.json`: typed nodes and edges such as
  `input -> state` (`writes_state`), `boundary -> state`
  (`constrains_state`), and `state -> output` (`drives_readout`);
- `index/core_examples.json`: native mechanism anchors that can be searched or
  used as target-field examples.

The KG is not a validation claim. It is a receptor map: it says which variables,
boundaries, outputs, controls, and missing links the field literature appears to
make available for mechanism transfer.

The adapter is the stronger object for a large field corpus. It converts a folder
of papers into field-native receivers for the six public Hyperion routes:
state/carrier, operator apparatus, update or transport, admissibility logic,
readout rule, and falsifier. It also keeps substrate evidence separate from
field nouns by scoring universal substrate classes such as coordinate domains,
inner-product spaces, probability spaces, graph topology, metric manifolds,
lattices, bundles, quotient spaces, and stoichiometric spaces.

For example, a neuroscience folder should not become a list of local names such
as cortex, connectome, synapse, or Markov blanket. The adapter should ask which
substrate those names instantiate: coordinate domain, graph topology,
probability space, inner-product space, phase space, or another universal
substrate class. Field-specific terms remain evidence; they are not the ontology.

## Main Route

### Mechanism translation

`fieldbridge translate` asks: if the same mechanism were constructed in a new
field, what would the target-field version look like?

It does not merely rename keywords. It produces a target formulation containing:

- the conserved invariant;
- target-field state, input, boundary, and output variables;
- candidate equations or equation skeletons;
- measurements that would make the mechanism observable;
- falsifying controls that would kill the translation.

This mode is the more constructive workflow. It is meant for designing new
experiments, inventions, and theory candidates, while keeping the evidence
boundary explicit.

The intended outcome is concrete: not "biology is like materials," but a proposed
state variable, a candidate equation, a measurement, and the control experiment
that would make the translation fail.

## Supporting Route

### Analog search

`fieldbridge search` asks: where has this mechanism already appeared? It returns
records from the field packs whose route/fiber profile, equations, variables,
measurements, and controls match the source mechanism.

This mode is useful for literature navigation, paper-to-paper comparison, and
finding equations in another field that play the same operational role.

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

It can then compare that mechanism with another paper, search a field pack for
analogous mechanisms, or translate the mechanism into a new target field.

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
