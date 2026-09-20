# Reproduce a construction and extend it to a new question

The supplied examples separate three questions that otherwise become mixed:
whether a transfer preserves the equations, whether a missing term can be
derived, and whether the resulting consequence is new. The first two are
mathematical calculations. The third requires a source search and comparison
with existing theory or experiment.

Run the three calculations independently:

```bash
python3 -m pip install -e '.[construction]'
python3 -B -m fieldbridge verify-construction examples/construction/affine_transfer.json --out-dir build/affine_transfer
python3 -B -m fieldbridge verify-construction examples/construction/ito_square.json --out-dir build/ito_square
python3 -B -m fieldbridge verify-construction examples/construction/quantum_correlations.json --out-dir build/quantum_correlations
python3 -B -m pytest -q -p no:cacheprovider
```

After installing dependencies, the calculations need no network, learned model,
source archive or cluster. Install
pytest separately to run the development tests. Each calculation writes JSON
and Markdown, preserving the question, assumptions, provenance and input hash.
Authored examples carry `source_alignment_verified: false` and
`novelty_established: false`, even when their equations verify exactly.

MorphWiki can run these same inputs through its
`scripts/build_construction_companion.py`. Its output combines the reports,
copies the input records, and records the calculation source hashes and
Python/SymPy versions. The repositories remain independently usable.

## From the example to a discovery candidate

Begin with a physical response to explain or produce. Recover the full source
equations and the assumptions that allow their use. A field-pack match or
atlas neighbour can locate a promising source neighbourhood, but it does not
establish a transformation between the two displayed endpoint equations.

Specify the proposed state correspondence or retained observables. For a
stochastic transfer, the calculation needs the actual noise law, not merely
the word "diffusion". For a quantum reduction, it needs the Hamiltonian and
the measured operator, not only their field labels. This is where a proposed
analogy becomes a mathematical question.

Calculate the target equation and its measurable consequence. Then remove
the assumption or term responsible for that consequence. In the supplied
examples, omitting quadratic variation changes the mean, and omitting the
spin correlation removes the information that fixes the initial slope. Such
controls identify the physical reason the construction works.

A prospective example should predict an unprovided response or a new regime
of a specified system. Its parameter choices, observable and comparison
method should be fixed before examining that response. Numerical agreement
with a supplied full model verifies the reduction; agreement with a new
experiment tests whether that full model describes the system. Neither test
alone establishes priority over existing literature.

## Extending the implementation

[The inverse-design example](14_inverse_construction.md) now makes one such
extension runnable. It solves for Ising couplings compatible with independently
variable exchange, calculates the polarization and reports an infeasible
nearest-neighbour restriction. It illustrates a research calculation; novelty
of a future construction remains a separate question.

The current handlers cover scalar Ito maps and finite closed Hamiltonian
dynamics. A new handler must supply a derivation, an independently checked
consequence, and a deliberately broken comparison. Boundary-value problems
need explicit domains; reduced open-system dynamics need the allowed initial
states and environmental couplings. Those problems are not handled by adding
another text field to the existing report.

The [retrieval adapter](13_retrieval_to_calculation.md) now joins a selected
retrieved source to a supplied correspondence and executes the calculation.
It still requires a typed source annotation and an explicit map or observable.
An automatic proposal stage could select these, while the calculation module
tests the candidate without seeing a reference answer. A useful comparison
would give alternative proposal methods the same sources and mathematical
tools and score their verified consequences. The current examples establish
the calculation path; they do not measure that proposal advantage.

## A practical research record

Use the following structure when extending a calculation to a new physical
question. Each entry should refer to an equation, an assumption or a result
that another researcher can inspect.

| Entry | Example from the coupled-spin problem |
| --- | --- |
| Question | Which expectations determine a transverse signal under a changed interaction? |
| Source model | The full Hamiltonian, admissible preparations and measured operator |
| Proposed change | Add a transverse field while retaining the measured operator |
| Derived consequence | An additional independent observable enters the closed span |
| Independent comparison | Evolve the full density matrix and compare the signal |
| Control | Remove the additional field and recover the two-observable result |
| Remaining research question | Does this change create a useful response in a specified physical implementation? |

The [transverse-field exercise](11_quantum_closure.md#exercise-let-a-field-compete-with-the-interaction)
executes the mathematical step. A discovery claim would require the additional
target-system result and literature comparison, not simply a larger span.

Record failures as carefully as successes. An incorrect state map should be
discarded; a numerical discrepancy should be checked under refinement; a
missing term should be derived and tested through its observable consequence.
These outcomes give different directions for the next calculation.

[Tutorial](index.md)
