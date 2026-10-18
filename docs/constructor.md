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

Two rules are implemented this way.

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

`memory predictions DIRECTORY` applies the rules to every specification it finds. For hysteron networks it is rule
R8, for bistable networks of kind `order` the holonomy of the writes, and for drift equations and networks (kind
`equations` or `network`) the three write targets of `memory codiscover`: the symmetric write, the threshold write
and phase locking, each derived under a time budget (`--budget`, 120 s by default), with the law constant where the
law is checked. Where the symmetric write is obstructed, the completion rung of `memory complete` runs on the same
specification: it keeps what the obstruction names (a second parameter that cancels a bias, the full critical space,
the next order, the distant branch) and derives the law of the completed model with its certificate. The report lists,
per specification, the status of each target and the completion with its law; `--workers` runs the specifications in
parallel processes.

The first set run this way is the thirty-one elements of quantum-computing hardware written as semiclassical drift
equations in [examples/quantum/semiclassical](../examples/quantum/semiclassical/README.md); their report is kept in
[docs/predictions/quantum_semiclassical.md](predictions/quantum_semiclassical.md).

## A rule with its construction named and not yet built

- **Closure.** An observable returns to itself under a sequence of operations only when the algebra it generates
  closes. Reference: the generators restricted to the closed subalgebra. Excess: the part of the observable's motion
  outside the closed subalgebra.

`fieldbridge.memory.counterpart.RULES` lists all four with their constructions and their status.

## What the constructor does not do

It does not establish that a model describes a physical system, and it does not establish novelty; both need a
comparison with experiment and with the literature. The ensemble laws are typical values in the ensembles in which
they were measured, and the reports quote them as comparisons, never as predictions for a lattice or a designed
network. A positive index makes the consequence possible; the measurement decides.
