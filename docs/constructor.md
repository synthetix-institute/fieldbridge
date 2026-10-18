# The constructor: from a structural rule to a prediction

FieldBridge states structural rules: a property of a material follows from the signs or the algebra of its
equations. A rule on its own classifies. The constructor turns a rule into a measurement and a prediction in four
steps, the same for every rule.

1. **Index.** Count the smallest edit of the material that makes the rule's hypothesis hold: edges whose signs
   must change, writes that must be removed from a node, couplings that must be cut. The index is zero when the
   hypothesis holds and grows with the distance from it.
2. **Reference.** Build the material with that edit, keeping everything else: the graph, the magnitudes, the
   fields, the thresholds, the writes. Where the minimal edit is not unique, build the unique reference that
   removes the whole effect instead (every coupling cooperative; every coupling removed) and keep the minimal edit
   as a second comparison.
3. **Excess.** Measure the same quantity on the material and on its reference. The difference is the consequence
   of the violated rule in a unit an experiment can compare: bits, a failure rate, a number of distinguishable
   states.
4. **Law and prediction.** Measure the excess over an ensemble to get the typical value as a function of the
   index, then name the published systems in which the index is not zero and the excess has not been measured.

Three rules are implemented this way, and the write point as a construction.

## Return to a turning point (rule R8)

Hysterons under a scalar drive return to every turning point when the signed network of couplings with the drive
as a node is balanced. Index: the frustration index F of that network. Reference: the ferromagnetic reference, every
coupling positive and every element driven upward, unique; the balanced counterpart, the F unsatisfiable edges
reversed, as the minimal edit. Excess: the number of distinguishable states, the deep memory and the failures of the
exact return test, measured by exhaustive enumeration of the transition graph. Law, from the sequence-memory study
on random networks of up to 64 elements: the gain over the reference is about 0.11 to 0.17 bit per unsatisfiable
edge above a threshold that grows with the degree, almost all of the added states need three or more reversals of
the drive, and the first switching fields of a descending sweep identify them.

```bash
python3 -B -m fieldbridge memory counterpart examples/counterpart/frustrated_network.json --out-dir out/counterpart
python3 -B -m fieldbridge memory predictions examples/memory/hysterons --out-dir out/memory_predictions
```

`examples/counterpart/platforms_R8.md` lists eight published platforms, what each paper measured and what the rule
asks to measure. The rule is sufficient and not necessary: antiferromagnetic chains of 12 elements return exactly in
three of four cases with F = 6, gaining states over the reference without deep memory, because a chain has no
avalanches to lock; the return test, not the gain, decides whether a memory is beyond return-point.

## The order of writes (the holonomy rule)

Two writes leave a record of their order only where they fail to commute on the coordinates that retain the
state; non-commutation is necessary and not sufficient, since a device that overwrites records only the last
write. Material: bistable coordinates coupled through the Laplacian of a graph, written by force pulses on supports
and released. Index: the smallest number of edits, writes removed from a node or couplings cut, after which every
write acts on an independent component and all writes commute, H. Reference: the uncoupled network with the same
writes, in which each element records the last write that moved it, unique. Excess: the number of distinct final
states over all orders of the writes and the split, the number of last writes after which the final state still
depends on the earlier order. Law, from the composition census over every connected graph with three to seven
nodes: with attractive couplings the share of graphs that record the order falls from 91 % for trees to 16 % for
dense graphs, with antiferromagnetic couplings it is 26 to 27 % for trees, sparse and medium graphs and 63 % for
dense ones.

```bash
python3 -B -m fieldbridge memory order examples/counterpart/order/ring8.json --out-dir out/order
python3 -B -m fieldbridge memory predictions examples/counterpart/order --out-dir out/order_predictions
```

The four examples under `examples/counterpart/order/` show the measurement. The ring of eight (writes on five nodes,
coupling 0.13) reproduces the census control state for state: four distinct final states over the six orders and the
order retained after each of the three last writes, with the uncoupled reference at four states and two retained
orders, so the ring's record is that of its elements. The star of seven with attractive couplings records the order
(four states, three retained, one bit over its reference). The complete graph of six with attractive couplings ends
every order in the same state: consensus erases the order and the last write, one bit below the reference. The same
graph with antiferromagnetic couplings keeps all six orders apart, 1.6 bits over the reference. The index H is 10 to
12 edits in each case; it says that a record is possible, and the couplings' sign decides whether it is kept.

## The write point (construction and verification)

A write needs a point along the control where the stable state loses stability. The counterpart holds the control
short of that point: `write_point_counterpart` cuts the control's range by a margin on the side of the write point,
and `write_point_check` verifies that the derivation of the symmetric write finds write points on the realization and
none on the counterpart. The excess is the selection probability of the written state, the law constant that
`codiscover` measures on the realization, against one half on the counterpart; that measurement uses the existing
`memory codiscover` command on both specifications.

## The write rules over a directory of drift specifications

`memory predictions DIRECTORY` applies the rules to every specification of the directory (its subfolders with
`--recursive`). For hysteron networks it is rule
R8, for bistable networks of kind `order` the holonomy of the writes, and for drift equations and networks (kind
`equations` or `network`) the three write targets of `memory codiscover`: the symmetric write, the threshold write
and phase locking, each derived under a time budget (`--budget`, 120 s by default), with the law constant where the
law is checked. Where the symmetric write is obstructed, the completion rung of `memory complete` runs on the same
specification: it keeps what the obstruction names (a second parameter that cancels a bias, the full critical space,
the next order, the distant branch) and derives the law of the completed model with its certificate. Where the
certificate deviates, a hold at the end value of the control tells an equilibrium write (the noise exchanges the
choice between shallow written states) from a failure in the nonlinear stage. A sweep through a cusp found by the
second parameter is read on a window narrowed until the quartic term of the drift no longer appears as a quadratic one,
and its law is refused when the sweep is not a pitchfork; with `--curvature`, `memory complete` also finds the curvature
of the sweep at which its intrinsic offset vanishes (tutorial Module 15). The report lists, per specification,
the status of each target and the completion with its law, class, response constant and offset; `--workers` runs
the specifications in parallel processes.

The first set run this way is the thirty-one elements of quantum-computing hardware written as semiclassical drift
equations in [examples/quantum/semiclassical](../examples/quantum/semiclassical/README.md); their report is kept in
[docs/predictions/quantum_semiclassical.md](predictions/quantum_semiclassical.md).

## Audit of a derived law

`memory audit SPEC` derives the symmetric write of a specification with its law and then asks, for each assumption
that the derivation does not test, what error would follow if the assumption failed, and runs the one test that
decides it. Six checks are implemented: the class of the reduced equation read on half the window of the fit; the law
measured with the bias toward either state, so that an offset of the write separates from the response; the noise
ratio at the element's own noise, where the constant is c(Γ) of the noise-ratio law rather than π^(1/4); the law
measured again with half the time step; for a semiclassical description of a quantum element, the number of quanta
that separate the written states from the unwritten one (below ten the description is not reliable); and the signs of
the control and of a second parameter at the write point against their published values. Each check ends as
confirmed, qualified (the law holds under a stated condition, such as a range of sweep rates), refuted, not applicable
or not testable. The verdict is "ready for review" when no error is confirmed: the law then goes to a person or to the
literature, and the audit does not certify it.

```bash
python3 -B -m fieldbridge memory audit examples/memory/pitchfork.json --out-dir out/audit_pitchfork
```

On the canonical pitchfork the audit finds no error and one condition (π^(1/4) applies at its own noise for sweeps
faster than 0.08); on the cusp of a driven Kerr resonator it confirms the offset of a straight sweep and that the
written states are a few quanta apart; on a model written only at a negative rate it confirms the sign error.

## Changes of a specification and the closure of their pairs

`memory changes SPEC` changes one part of a drift specification at a time and lets the write rules judge the result
([tutorial Module 14](tutorial/34_memory_changes.md)). Four kinds of change are applied: one additive term of one
equation removed; another parameter made the control and swept from a quarter to four times its value; the range of
the control extended by ten times its width below it, above it or on both sides; and, with `--donor`, one term of a
donor specification added to the equation in the same position, in the specification's own variables, with the
donor's parameters set to one and a new coefficient. A changed specification is kept only when its states stay
bounded: short integrations from large amplitude, with the control at both ends of its range, must not run away or
still be growing at the end. The test rejects a loss rate made negative and a saturating term removed, and keeps a
widened detuning; it integrates because a static radial test fails both ways (a Kerr oscillator above threshold has
an outward radial drift and is bounded by its rotation; a frustrated triangle without one cubic term has an inward
drift on average over directions and runs away along one axis).

Each kept specification is derived for the symmetric and the threshold write, every derivation in a child process of
its own under a wall-clock limit (`--budget`), so that an exception inside a solver's callback costs one derivation
and not the run. A change is a step when it reaches a write the specification does not reach, or changes the class of
the reason why a derivation stops (a one-sided write, several critical directions, a vanishing cubic term, a distant
branch, an oscillation). With `--pairs`, the first changes are combined with every change of the changed
specification, and a pair is tested for closure: it closes when its outcome equals the outcome of the specification or
of a single change, or combines classes that single changes produce, and it is open when it reaches a write that no
single change reaches or a class that no single change produces. On the quantum-computing elements the open pairs
were of three kinds: a second control that the specification did not declare (the detuning of a Kerr oscillator,
swept past the parametric threshold), a lifted degeneracy (unequal losses on a frustrated triangle), and changes
without a physical realization, which the boundedness test now rejects. The first two are known physics; the changes
locate where a mechanism appears and do not decide whether it is new.

```bash
python3 -B -m fieldbridge memory changes examples/quantum/semiclassical/kpo_pair.json --pairs --out-dir out/changes_kpo
python3 -B -m fieldbridge memory changes examples/memory/pitchfork.json --donor examples/memory/laser.json --out-dir out/changes_pitchfork
```

## The frequencies of an observable (the closure rule)

Under a Hamiltonian the expectation of an observable moves inside the closure of the observable, the smallest space
of operators that contains it and is invariant under the commutator with the Hamiltonian (tutorial Module 24, §6).
The generator of that closure has pairs of eigenvalues ±iω, and these ω are the only frequencies the signal can show
from any prepared state. The rule is sufficient and not necessary: a frequency of the closure appears with the
amplitude the prepared state gives it, which can be zero. Index: the smallest number of terms of the Hamiltonian, as
the specification writes them, whose removal leaves a closure of two or three operators with a rotation frame.
Reference: the Hamiltonian without those terms, the closed counterpart, which rotates with one frequency and follows
the Rabi law. Excess, from a top eigenstate of the observable: the spread of the signal's lines, the weighted
standard deviation of its Bohr frequencies, which the counterpart's single line lacks; the shift of their weighted
mean from the counterpart's frequency; the weight outside the line that continues the counterpart's rotation; and
the closure frequencies the prepared state leaves silent. For a realization of index 0 the excess is the amplitude
offset of the signal from sin²θ, zero when the prepared state lies along the observable in the closure. Comparison:
with the removed terms scaled by λ the spread falls as λ^p; p = 1 when the terms split the frequency at first order.

```bash
python3 -B -m fieldbridge quantum closure --out-dir out/closure
python3 -B -m fieldbridge memory predictions examples/quantum --out-dir out/quantum_predictions
```

Over the thirteen quantum examples and controls every line of every signal lies in the closure spectrum. Eight
realizations rotate (index 0). The nitrogen-vacancy centre has a closure of seven operators with three frequencies;
removing the anisotropy D J_z² leaves a rotation at the Larmor rate, and the signal's three lines spread over 6.5
times that rate, falling as λ^0.9 when the anisotropy is scaled down. The Bose–Josephson junction with interaction
has a closure of 21 operators and ten frequencies; removing the interaction leaves the Josephson rotation, and the
ten lines spread over 0.31 of its frequency. The exchange chain of Module 14 and the control with a field on both
spins have index 1 with two lines each; two nuclear spins of different species need both Zeeman terms removed
(index 2). The control with a weighted observable has index 0 and one line, with an amplitude offset of 0.10: the
rule's necessary direction fails there, as the language reports.

`fieldbridge.memory.counterpart.RULES` lists the four rules with their indices, constructions, excesses and commands.

## What the constructor does not do

It does not establish that a model describes a physical system, and it does not establish novelty; both need a
comparison with experiment and with the literature. The ensemble laws are typical values in the ensembles in which
they were measured, and the reports quote them as comparisons, never as predictions for a lattice or a designed
network. A positive index makes the consequence possible; the measurement decides.
