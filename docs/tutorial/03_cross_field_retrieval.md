# Find a useful source record

A query about gradient drift and diffusion should bring back records whose
equations and assumptions are worth inspecting. FieldBridge ranks records
for that purpose; the ranking is the start of the correspondence, not its
derivation.

```bash
python3 -B -m fieldbridge search examples/brownian_probability_flow.tex \
  --target-field stochastic_optimization --top-k 5
```

Without optional atlas data, this searches the bundled field records. The
response shows each record's identity, score, overlapping route and fiber
evidence, and keyword matches. The field filter is supplied by the user; this
command is not a field-blind discovery test.

## Understand the score

[score_record](../../fieldbridge/search.py) combines cosine similarity of the
fingerprint with route overlap, fiber overlap and a capped keyword bonus.
Thus the result depends on both the rules and the stored examples. A score
of 0.8 is a ranking value, not an 80% probability of physical equivalence.

[database.load_all](../../fieldbridge/database.py) reads field packs and
mechanism records. The optional atlas witnesses are loaded with
`search --include-hyperion`. They can point to source material beyond the
small authored seed set, but their presence does not turn a target formulation
into a derived equation.

## Ask for a target formulation

```bash
python3 -B -m fieldbridge translate examples/brownian_probability_flow.tex \
  --to stochastic_optimization --no-hyperion
```

[translate_mechanism](../../fieldbridge/search.py) selects existing target
variables, equations, measurements and controls, with fallback formulations
when needed. For example, a gradient-diffusion query can retrieve a noisy
parameter-update example. The program is presenting that stored model in the
context of the query; it has not proved that the particle system and the
optimization process have corresponding trajectories.

To establish a relation, a researcher must specify the map between states,
the domains, the parameter correspondence and the observable to preserve.
The [calculated-source tutorial](13_retrieval_to_calculation.md) executes this
next step for supported annotations.

**Exercise.** Open the highest-ranked record in
[data/index/core_examples.json](../../data/index/core_examples.json).
Find its equation, assumptions and references. Which part of the proposed
correspondence is supported by that record, and which part is still a choice?

[Next: read a constructor proposal](04_constructor_transfers.md) · [Tutorial](index.md)
