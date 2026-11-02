# Learning a construction from a failed reduction

A cubic amplitude equation predicts saturation proportional to the square root
of the growth rate. If its cubic coefficient vanishes and a fifth-order term
provides saturation, that prediction changes. The constructor should calculate
the term responsible for the change rather than search for parameters that make
the original cubic equation appear successful.

This chapter runs a bounded learning calculation. It starts with complete source
equations, a proposed set of observables and equations for their evolution. The
program calculates what the proposed equations omit, derives a correction, and
tests it on preparations and parameter settings excluded from candidate selection.
It then retains the successful form of the correction for later constructions.

## 1. Derive a missing term

The first example supplies the source equation and a proposed cubic description:

$$
\dot x=ax-bx^5,\qquad u=x,\qquad \dot u=au-bu^3.
$$

Differentiating the observable along the source dynamics gives a residual
$b x^3-b x^5$. The constructor searches polynomial terms in $u$ and solves for
their coefficients by equality of the equations, without fitting a trajectory.
The correction $b u^3-b u^5$ removes the cubic term and supplies fifth-order
saturation. For $a,b>0$, the nonzero steady amplitudes satisfy
$|u|=(a/b)^{1/4}$; the cubic description would give $(a/b)^{1/2}$.

```bash
python3 -B -m fieldbridge memory learn examples/learning/quintic_saturation.json \
  --library build/constructor_learning/moves.sqlite \
  --out-dir build/constructor_learning/quintic
```

The JSON report contains the original residual, attempted polynomial supports,
derived equations and trajectory comparisons. `verified_extension` means that
the equation identity and the specified numerical tests succeeded. This example
reproduces established amplitude-equation mathematics from supplied equations.

## 2. Retain information that the observable needs

For coupled deviations from equilibrium,

$$
\dot x=-a x+b y,\qquad \dot y=-c y,
$$

a description retaining only $u=x$ omits $b y$. Different hidden preparations
$y(0)$ give different visible responses even though $x$ never drives $y$.
No polynomial equation for $u$ alone can reproduce all of them. The constructor
therefore attaches $y$ and derives its evolution as part of the completed model.

```bash
python3 -B -m fieldbridge memory learn examples/learning/hidden_exchange.json \
  --library build/constructor_learning/moves.sqlite \
  --out-dir build/constructor_learning/exchange
```

An attached observable need not be an original state coordinate. Replacing
$b y$ by $b y^2$ makes the residual point to $v=y^2$. Its equation is
$\dot v=-2cv$, so $(x,v)$ closes without resolving the sign of $y$.
The experiment below includes this case. It demonstrates an exact predictive
description, not an increase in the physical dimension of the source system.

## 3. Measure what improves with experience

```bash
python3 -B scripts/check_constructor_learning.py \
  --out-dir build/constructor_learning_comparison
```

Use a fresh output directory. The script first learns the fifth-order correction,
then applies its support to a different equation with both cubic and fifth-order
terms. A separate run starts with an empty library. Compare
`transfer_warm.candidate_evaluations` with
`transfer_without_library.candidate_evaluations` in `experiment.json`.
The warm run must use fewer candidate supports and pass the same tests.
Coefficients are solved again from the new source equations; the stored move is
not a copied answer.

The final example changes coordinates to a moving branch, $u=x-m$ with
$\dot m=v$. The residual gives the missing drift $-v$. The constructor keeps this
asymmetry in the predicted dynamics rather than forcing a symmetric-write law.

## 4. What qualifies a new move

For a source vector field $F$, observable map $h$ and proposed target vector
field $G$, the exact check is

$$
D h(x)F(x)-G(h(x))=0.
$$

The verifier reparses the candidate and recalculates every component of this
identity. Numerical evaluation then compares the full source dynamics with the
constructed equations, using different integration methods and tighter step
sizes. Removing the correction must change the prediction by a resolved amount.

The task file separates `development` from `withheld` interventions. Only the
development results select the candidate. The program freezes its hash before
evaluating withheld results. A withheld failure ends that campaign: those same
responses are not used to select another candidate. A subsequent study needs
fresh withheld interventions. At least one withheld parameter setting must be
absent from development whenever parameters are present.

The persistent library records which interventions have already been evaluated
for the source vector field. Reusing them in a later campaign returns
`withheld_already_used`, even if their labels or integration durations change.
Tests are marked as used before their responses are calculated, including tests
that fail. This guard applies within the chosen library; separate libraries do
not share a global experimental registry.

The SQLite library stores equation recipes, calculation records and unsuccessful
attempts. It contains no executable generated code. A previously successful
recipe is a proposal on a new model; the equation identity and numerical checks
are repeated before accepting the new application.

## 5. Keep the existing constructor's unresolved results

```bash
python3 -B -m fieldbridge memory departures path/to/records.jsonl \
  --library build/constructor_learning/moves.sqlite
```

The importer retains target mismatches, numerical failures, measured law
departures, and successful target classifications as different kinds of record.
Existing completion results retain their equations, source references and
reported discrepancies. Importing a record never makes it a learned move.
To collect outcomes during a new calculation, add
`--learning-library build/constructor_learning/moves.sqlite` to `memory record`.

An unexplained difference in a stochastic write probability still needs a
convergence study and a suitable stochastic derivation. The polynomial learner
does not reinterpret it as a deterministic term. Its current scope is autonomous,
noiseless polynomial equations on real Euclidean coordinates, with at most four
source coordinates. Boundary-constrained systems, stochastic calculus, periodic
coordinates and non-polynomial equations remain explicit unsupported cases.
Learning from those cases requires another calculation backend with its own
independent checks. Source fidelity, dimensional consistency and experimental
applicability also require evidence outside this algebraic calculation.

## 6. Repair a recorded failure

The exchange example can now enter through the same record library used by
`memory record`. Run the complete example in a fresh directory:

```bash
python3 -B scripts/check_repair_campaign.py \
  --out-dir build/repair_campaign
```

It records the incomplete equation for $x$, derives the additional coordinate
$y$, and retains the checked correction. A second system has twice the coupling
to $y$. Reusing the correction solves its coefficients again and tests new
preparations. `campaign_comparison.json` contains the equations, the original
error, the repaired prediction and the comparison with enumeration and direct
symbolic closure.

For an existing record library, use:

```bash
python3 -B -m fieldbridge memory repair \
  --library build/constructor_learning/moves.sqlite \
  --families build/constructor_learning/families.json \
  --out-dir build/constructor_learning/repairs
```

An eligible record contains `spec` and a target's `learning_task`. The task has
the same fields as the JSON examples above: the source model, retained
`coordinates`, `candidate_drift`, and separate `development` and `withheld`
interventions. Its source model must equal the recorded specification. A record
that says only "no write point found" has no proposed equation to correct; it
returns `needs_equation_task`. Stochastic equations return
`unsupported_or_invalid` for this backend. Neither result creates a recipe.

The library distinguishes an incorrect term in a closed set of coordinates
from a missing coordinate. For example, $\dot x=-ax+bxy$ with $\dot y=-y$
requires $y$ even though the defect contains both the retained and discarded
variables. An explicit closure check is repeated after every enlargement.
Search limits apply equally to recalled and newly derived proposals.

A counterexample can suppress a previously attempted hint for the same
equations and search policy. The exclusion is evaluated before proposing that
hint again. Increasing the coordinate budget permits a new attempt; exhausting
a budget never labels the underlying dynamics impossible. Family snapshots are
frozen during transfer tests, while a separate SQLite registry reserves each
withheld intervention before evaluating it. Repair campaigns have one learning
writer; parallel evaluation uses frozen snapshots and separate experiment
registries.

## 7. Separate memory from a stronger solver

```bash
python3 -B scripts/check_repair_transfer.py \
  --out-dir build/repair_family_comparison
```

This comparison trains on amplitude and exchange equations, freezes the
libraries, then evaluates cascade kinetics, a product interaction, moving
coordinates, crowding and mixed hidden-variable feedback. A four-stage cascade
exceeds the declared two-coordinate attachment budget and should remain
unresolved. The development revision is completed before the libraries are
frozen; withheld outcomes never revise the library being tested.

Read `summary.totals` in `transfer_experiment.json` for four strategies:

| Strategy | Source of its proposals |
|---|---|
| `blind` | Polynomial enumeration |
| `symbolic` | Direct algebraic closure, with no memory |
| `recipes` | Stored monomial supports, followed by enumeration |
| `families` | Verified supports under their applicability conditions, followed by the family solver and enumeration |

Compare verified outcomes, candidate counts, Groebner reductions and elapsed
time. The family algorithms are supplied mathematical procedures. A memory
advantage requires improvement over the same direct solver without memory,
not only over enumeration. Exact reproduction of supplied dynamics tests model
correction; a new physical mechanism additionally changes an interaction or
driving condition and predicts a consequence of that change.

## Implementation and tests

- [Learning engine](../../fieldbridge/memory/learning.py)
- [Executable comparison](../../scripts/check_constructor_learning.py)
- [Regression and negative-control tests](../../tests/test_constructor_learning.py)
- [Existing completion of memory models](../../fieldbridge/memory/completion.py)
- [Conditional repair library](../../fieldbridge/memory/repair_families.py)
- [Recorded-failure replay](../../scripts/check_repair_campaign.py)
- [Four-way transfer comparison](../../scripts/check_repair_transfer.py)

```bash
python3 -B -m pytest tests/test_constructor_learning.py tests/test_repair_families.py tests/test_periodic_construction.py
```

Next: [What counts as a discovery](12_reproduction_and_discovery.md).
