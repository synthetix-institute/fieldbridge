# Contributing to FieldBridge

FieldBridge is extended by adding material specifications, mechanisms and
analyses. A material specification gives the equations of one material with
their source and assumptions. From it the programs calculate the stable
states, the write points, the writing protocols and the law of loss, and
compare them with materials from other fields. The
[catalog of materials](docs/materials.md) lists the specifications in the
repository, and the [list of wanted materials](docs/wanted_materials.md) names classic models from fields
that are not covered yet.

## What a contribution returns

- **The memory card of the material.** For every specification a pull request adds, the checks of the pull
  request compute its card: the stable states, the write points and their normal forms, the writing protocols and
  the law of loss. The card is shown in the summary of the "material card" check; nothing has to be installed
  to read it.
- **Its connections to other fields.** The same summary lists the materials of the catalog that are written by
  the same mechanism, for example a threshold write at a fold or a symmetric write at a pitchfork, with their
  fields. A model from ecology can turn out to be written like a laser or a genetic switch; the derivations of
  [Module 9](docs/tutorial/23_memory_codiscovery.md) then compare the two routes step by step.
- **A place on the web page.** After the next build the material appears on the
  [web page](https://synthetix-institute.github.io/fieldbridge/): in the column of its mechanism on the map, joined to
  a realization from another field that FieldBridge finds written by the same mechanism, with its equations and its
  dynamics calculated in the browser.
- **Credit.** The catalog names the source of every material and its contributor (`provenance.contributor`).

## How contributions add up

Every material added to the catalog is compared with every material already in it. The more fields the catalog
covers, the more a new material connects to, and the more each mechanism is tested: a threshold write found in
ten fields, each with the same law constant, is a stronger result than one found in two. Mechanisms, carriers and
analyses added to the code apply at once to every material in the catalog.

## Adding a material

### 1. Installation

```bash
git clone https://github.com/synthetix-institute/fieldbridge.git
cd fieldbridge
python3 -m venv .venv && source .venv/bin/activate
python3 -m pip install -e '.[memory]' pytest
```

### 2. Template

```bash
python3 -B -m fieldbridge memory new liquid_crystal_cell --carrier torus
```

This writes `examples/memory/liquid_crystal_cell.json`. The file loads and
runs, but its dynamics and descriptive fields are placeholders. Choose the
carrier by what the variables are:

| `--carrier` | Variables | Examples |
| --- | --- | --- |
| `euclid` | real numbers | an order parameter, a position, a current |
| `orthant` | non-negative amounts | concentrations, populations, conductances |
| `torus` | angles | orientations of rods or spins, phases of oscillators |

### 3. Specification fields

| Field | What to write |
| --- | --- |
| `question` | The question the model answers, for example "Does a nematic cell anchored at both plates store two orientations, and how is one written?" |
| `provenance.source` | Where the equations come from: authors, journal, volume, page and year, a DOI or an arXiv identifier |
| `provenance.contributor` | Your name or GitHub handle |
| `field` | The field of the material, for example "soft matter" |
| `carrier.variables`, `drift` | One expression per variable, in `+ - * / **`, `sqrt exp log sin cos tanh abs` and `pi` |
| `parameters` | Every constant of the drift, with the values stated in the source |
| `control` | The parameter an experiment varies, with its range; write points are searched along it |
| `closure` | What is specified or discarded to obtain closed equations: the bath (`noise`), boundaries, imposed conservation |
| `observable` | What an experiment measures |
| `assumptions` | Every modelling choice a reader could dispute |

[Module 2](docs/tutorial/16_memory_specification.md) explains each field on
a worked example, and the [glossary](docs/tutorial/memory_glossary.md)
defines the terms.

### 4. Checks

```bash
python3 -B -m fieldbridge memory check examples/memory/liquid_crystal_cell.json
python3 -B -m fieldbridge memory card examples/memory/liquid_crystal_cell.json \
  --out-dir build/liquid_crystal_cell
```

`memory check` lists required checks (the file loads; the question,
source, field and assumptions are stated) and recommended ones (a citable
reference, a control parameter, the closure and the observable, and
agreement between the structural predictions and a quick calculation), and
the materials of the catalog written by the same mechanism. It exits with
status 1 until the required checks pass. `memory card` writes
`card.md` and `card.png`: the stable states, the write point and its normal
form, the writing protocols and the law of loss.
[Module 1](docs/tutorial/15_memory_first_card.md) shows how to read a card.

### 5. Pull request

```bash
python3 -B -m fieldbridge memory catalog --out docs/materials.md   # also writes docs/materials.json
python3 -B -m pytest -q -p no:cacheprovider
git checkout -b material/liquid-crystal-cell
git add examples/memory/liquid_crystal_cell.json docs/materials.md docs/materials.json
git commit -m "Add a nematic cell with two anchored orientations"
```

The branch is pushed to a fork and proposed as a pull request. The tests run
`memory check` on every specification in `examples/memory`, so a material
with placeholders or without a source does not pass.

## Requirements for a material

- The equations are those of a named source, and the parameter values are
  the ones it states. A value chosen for the example is listed under
  `assumptions`.
- The question is one that a measurement could answer, and the observable
  says what would be measured.
- What the programs find is a property of the model. The card does not show
  that the physical material behaves this way; that needs comparison with
  experiment. The pull request can say how to test it.
- Physical terms are used: observable and measurement, not readout; state
  and carrier, not latent variables.

A material can also be proposed without a specification, through the
[material proposal form](https://github.com/synthetix-institute/fieldbridge/issues/new?template=material.yml)
with its source; another contributor can then write the specification.

## Adding a mechanism

A mechanism enters the map of the web page when FieldBridge derives it in models from several fields, each by its
own chain of verified transformations, and certifies a law whose constant does not depend on the field (or, as for
return-point memory, a structural condition checked in each model). A new mechanism is a derivation target in
`fieldbridge/memory/codiscovery.py`; the module docstring describes the three targets that exist. A target needs:

1. **A canonical form and its check.** The reduced drift at the write point, rescaled so that its leading coefficient
   is fixed (letter K), and an invariant that the other coefficients must satisfy in every model, for example no even
   part at a pitchfork or no linear part at a fold. Give its entry in `TARGET_INFO`: `title`, `letters`,
   `canonical`, `check`, `target`, `law`, `constant_name`, `constant`, the meaning of the letters that differ from
   those of `LETTERS`, and `check_key` and `point`.
2. **A law with a field-independent constant**, estimated in every model by simulating the full realization (letter
   L), with an uncertainty. A constant from stochastic trajectories is split into independent runs on spawned
   random streams, as in `construct.swept_write_check`; a deterministic constant is extrapolated with a stated error
   floor, as in `fold_delay_law`.
3. **A derivation function** `derive_<target>(real, rng, check_law=True, n_traj=..., replicates=1)` that returns
   the fields of `derive_symmetric_write` (`status`, `word`, `steps`, `canonical`, `law`, `law_constant`,
   `write_point`) and, where a step fails, `status = "obstructed"` with the obstruction named (`obstruction`,
   `obstruction_short`). Add it to `TARGETS`.
4. **Tests** in `tests/test_memory_constructor.py`: models from at least two fields that reach the target by
   different words, a model whose derivation stops at a named step, and the law constant within its error of the
   expected value. Compare numbers differentially; do not pin solver residuals.

On the web page the target appears under *One mechanism in different fields* once `site_data.TARGETS` and `PREFIX`
include it and `site.js` (`collapse`) draws its figure. Its class on the map needs entries in `site_registry.py`
(`CLASSES`, `MECHANISMS`, `CLASS_SHORT`), a glyph in `web/mechanisms.js` (`GLYPHS`), and its definition, derivation
and original publications in `site_references.py` (`READING`), from which `docs/mechanisms.md` is regenerated
(`python3 -B -m fieldbridge mechanisms --out docs/mechanisms.md`). The law constants of the page come from
`docs/site/law_constants.json`, regenerated with `python3 -B -m fieldbridge demo --law --save-law-record`. Until its
law is certified, a mechanism is listed under *Mechanisms in preparation* (`PLANNED` in `site_registry.py` and
`docs/ROADMAP.md`).

## Other contributions

| Contribution | Where to start |
| --- | --- |
| A mechanism derived in models from different fields | [Adding a mechanism](#adding-a-mechanism), [Modules 9 to 11](docs/tutorial/23_memory_codiscovery.md) and `fieldbridge/memory/codiscovery.py` |
| A new carrier or analysis | `fieldbridge/memory/`, with tests in `tests/test_memory_constructor.py` |
| A quantum carrier for the language of mechanisms | [Chapter 24](docs/tutorial/24_spin_language.md), `fieldbridge/quantum/` and `examples/quantum/` |
| Retrieval, calculation and field packs | [Extensions](docs/tutorial/12_reproduction_and_discovery.md) and [adding a field](docs/NEW_FIELD.md) |
| A tutorial correction or figure | `docs/tutorial/`; the memory figures are produced by `docs/tutorial/figures/memory/make_figures.py` |

## Code conventions

- Read and write text files with `encoding="utf-8"`.
- Specifications are parsed by the restricted arithmetic parser in
  `fieldbridge/memory/spec.py`; never evaluate input with `eval` or
  `sympify`.
- Every calculated report carries the input hash, the hash of the
  implementation and the evidence boundary; keep them when adding a command.
- A new behaviour comes with a test that fails without it.

```bash
python3 -B -m pytest -q -p no:cacheprovider
```

By contributing you agree that your contribution is released under the
[MIT license](LICENSE).
