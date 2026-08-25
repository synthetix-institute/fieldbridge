# Adding a Field

FieldBridge ships four fields plus a worked fifth, `epidemic_dynamics`, added
purely as a template. This page walks through that fifth field end to end. It
is not an intelligence domain, which is the point: the atlas is arXiv-wide, so
nothing about the method is specific to the domains that happen to ship.

## What a field is

Two files, and nothing else.

| File | Role |
|---|---|
| `data/field_packs/<field_id>.json` | Vocabulary. Five word lists that let the retriever recognise the field's state, drive, boundary, readout and controls. |
| `data/index/core_examples.json` | One or more seed records. Supplies the *target formulation* a transfer is written into. |

The field pack is about 1 KB. The seed record is the part that takes thought.

## Step 1 — the field pack

`data/field_packs/epidemic_dynamics.json`:

```json
{
  "field_id": "epidemic_dynamics",
  "label": "Epidemic dynamics",
  "description": "Compartmental population flow on contact networks: incidence, recovery, conserved population and threshold behaviour.",
  "state_words": ["compartment fractions S, I, R", "prevalence I(t)", "..."],
  "input_words": ["transmission rate beta", "vaccination or isolation schedule", "..."],
  "boundary_words": ["conserved total population S+I+R=N", "non-negativity of compartments", "..."],
  "output_words": ["epidemic peak height and timing", "final attack rate", "..."],
  "validation_words": ["zero-transmission control", "population-conservation residual", "..."]
}
```

Write the words a paper in this field would actually use. The five lists map
onto the operational identity `I_op = (M; C, R, P)`: `state_words` and
`input_words` describe the carrier and drive, `boundary_words` the closure `C`,
`output_words` the readout `R`, and `validation_words` the controls that would
falsify a transfer.

`validation_words` is the list people skip. Don't: a transfer with no
falsifier is a restatement, and the constructor's `falsifier_present` gate
reads this list.

## Step 2 — the seed record

Append to `data/index/core_examples.json`. The fields that carry weight:

- `invariant` — one sentence naming what survives the transfer. This is what
  the constructor prints as the preserved contract, so vague text here makes
  every transfer into this field vague.
- `equations` — plain ASCII, not LaTeX.
- `variables` — prefix-tagged. `q:` state, `u:` drive, `B:` boundary or
  closure, `y:` readout. The constructor parses these prefixes to fill the
  target clauses; untagged variables are not picked up.
- `measurements`, `controls` — what a reader would measure, and what would
  show the transfer failed.
- `routes`, `fibers` — the six route and five fiber weights. Estimate them from
  the equations; they bias retrieval but do not gate it.
- `target_fields` — include the field's own id.

## Step 3 — run it

```bash
fieldbridge construct examples/brownian_probability_flow.tex --to epidemic_dynamics
```

The identity block shows a preserved operator and a re-attached carrier:

```text
source core   M  = (Omega08 + Omega04 + Omega14, Xi01)   carrier: p(x,t) probability density
target core   M' = (Omega08 + Omega04 + Omega14, Xi_unassigned)
                                                  carrier: compartment fractions S, I, R
```

The operator is preserved because Brownian probability flow and mass-action
incidence are both conserved-density transport. The target carrier reports
`Xi_unassigned` because re-attachment *proposes* a carrier; the atlas never
assigned one for this field, and printing the source token there would claim
an assignment that was never made.

Retrieved witnesses mix the new seed with atlas evidence:

```text
0.730  epidemic_dynamics: SIR compartmental flow with conserved population
0.782  EW000003427: A00 Omega08 + Omega04 + Omega14; arXiv 1905.02221
```

## Step 4 — check it discriminates

A field pack that returns the same operator for every input is not working.
Feed it two mechanisms that should differ:

```bash
fieldbridge construct <a wave or Schroedinger equation> --to epidemic_dynamics
fieldbridge construct <a diffusion equation>            --to epidemic_dynamics
```

Spectral and transport sources should resolve to different operator tokens. In
the shipped atlas a wave equation and a Schrödinger equation resolve to the
*same* tokens, which is correct: both are normal-mode operators on different
carriers.

## What you get, and what you don't

The retrieval side is arXiv-wide: 2,633 witnesses, none of them written for any
of these fields. The target side is not. There are eight seed records in
total, one or two per field, so the target formulation a transfer is written
into is the seed you wrote.

This is the honest asymmetry. Adding a field gives you real atlas retrieval
against a target vocabulary you authored. The tool proposes a structured
correspondence and names what would falsify it; it does not supply target-field
knowledge you did not put in.

## Getting the atlas

The atlas index is **not tracked in git**. It is a ~31 MB generated artifact
(index plus shards) attached to a release, because tracking it would add that
much immutable history on every refresh.

```bash
python3 scripts/fetch_atlas.py
```

Without it everything still runs: the constructor reports a route-derived proxy
instead of an assignment, says so in its output, and atlas-dependent tests skip
rather than fail. A red suite therefore always means broken code, never a
missing download.

## Refreshing the atlas

The current index was generated on 2026-06-07 from
`discoveries/equation_witnesses.jsonl` (3,824 witnesses, dated 2026-05-15) and
holds 2,633 records with 897 analog equations. It is a snapshot, not the
full-corpus run.

To regenerate from a newer witness export, in the KnowledgeParser tree:

```bash
python3 scripts/export_fieldbridge_static_index.py \
  --witnesses discoveries/equation_witnesses.jsonl \
  --out discoveries/fieldbridge_static_index/hyperion_static_index.json \
  --shard-dir discoveries/fieldbridge_static_index/shards \
  --max-records 3000 --max-per-paper 2
```

Then copy `hyperion_static_index.json` into `data/index/` and the shards into
`data/shards/`, and attach the index to a new release so other clones can fetch
it. Pass `--max-records` explicitly: the script default is 1200, below the
2,633 currently shipped, so omitting it silently shrinks the atlas.
`fetch_atlas.py` warns when an index falls below that count.

Regenerating from the current witness file reproduces the same index. A real
refresh needs the upstream witness extraction re-run first.
