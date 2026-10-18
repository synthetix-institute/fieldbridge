# Module 12. Return to a turning point

**Learning objectives.** After this module you can

1. describe a material as interacting hysterons under a slow drive, in a specification of kind `hysterons`;
2. state when such a material returns exactly to its state at a turning point of the drive, from the signs of its
   couplings and of the drive alone;
3. run `memory hysterons` and read the prediction against subloops of the dynamics;
4. explain why the frustration of the couplings does not decide the return, and why the structural condition is
   sufficient and not necessary;
5. build the balanced counterpart of a frustrated network with `memory counterpart` and measure the network's states
   against it;
6. follow a material given by equations under a slow drive, test its return on its own transition graph and reduce it
   to hysterons.

**Prerequisites.** Modules [3](17_memory_predictions.md) and [4](18_memory_writing_and_retention.md). **Time.** About
30 minutes.

## 1. The question

Many materials record the turning points of a scalar drive: magnets the extremes of the applied field, porous solids
those of the vapour pressure, sheared amorphous solids those of the strain. Their elements switch at thresholds of
the drive, behind barriers that are large compared with the thermal energy, so that the state depends on the
sequence of turning points and not on the rate of driving. The exact form of this record is return-point memory.
When the drive turns at $H_1$, makes an excursion that does not pass $H_1$ and comes back to it, the material returns
exactly to the state that it had at $H_1$. The excursion is erased, and the turning points that enclose it are
retained (Sethna et al., 1993; Keim et al., 2019).

Independent elements have this property by construction: each element with its own loop is the Preisach model. The
question of this module is what interactions do to it.

## 2. The model

Elements $s_i = \pm 1$ switch at thresholds of their local field

$$
f_i = \sum_j J_{ij} s_j + h_i + \eta_i H , \qquad \text{element } i \text{ keeps its state while } s_i f_i + b_i \ge 0 . \qquad (1)
$$

$h_i$ is a quenched random field, $\eta_i$ the sign (and strength) with which the drive acts on element $i$, and
$b_i \ge 0$ half the width of the element's own loop. With $b = 0$ this is the zero-temperature random-field Ising
model; with $J = 0$ and $b > 0$ it is the Preisach model. The drive changes quasi-statically and without thermal
noise. Between events nothing moves; at an event the element that has become unstable switches, and the avalanche
that follows is relaxed at fixed $H$, switching the least stable element first or, as a control, a randomly chosen
unstable one. The measured response is the observable conjugate to the drive,
$R = \sum_i \eta_i s_i / \sum_i \lvert\eta_i\rvert$: the magnetization in a uniform field, the staggered
magnetization in a staggered one.

In the language of the realization: $\Xi$ is the set of configurations, $\Omega$ the switching rule with the
couplings $J_{ij}$, $C$ the absence of thermal activation (the quasi-static, athermal closure), $R$ the response
above, $P$ the history of turning points of the drive, and $A$ the random fields, the half-widths and the strength of
the couplings.

## 3. The prediction from structure

Module 3 decided storage and frustration from the loops of the couplings: all couplings of a loop can be satisfied
if and only if the composition of their transports has a fixed point, and for couplings that act by sign this means a
positive product of signs around the loop. The return to a turning point is guaranteed by the same rule, applied to
a larger network. Count the drive as one more element, joined to every driven element $i$ with the sign of $\eta_i$.
Each coupling $ij$ then closes a loop through the drive with the product of signs $\eta_i J_{ij} \eta_j$.

If no loop of this extended network is frustrated, there are signs $g_i = \pm 1$ (Harary, 1953) such that the
relabeled states $\sigma_i = g_i s_i$ have only non-negative couplings $g_i J_{ij} g_j$ and a drive that pushes every
element the same way. The relabeled dynamics has the no-passing property (Middleton, 1992; Sethna et al., 1993): a
configuration that lies above another at every element stays above it under the same history of the drive. While
the drive rises only upward switches occur, so every avalanche ends with each element switched at most once, and its
final state does not depend on the order of the switches. The return follows from four comparisons. Let $A$ be the
state at $H_1$, $C$ the state passed at $H_2 < H_1$ on the way up, $B$ the state after the descent to $H_2$, and
$A'$ the state after the rise back to $H_1$. No passing gives $B \le A$ (compare with $A$ kept at $H_1$),
$B \ge C$ (compare with $C$ kept at $H_2$), $A' \ge A$ (the rise from $B \ge C$) and $A' \le A$ (the rise stays
below $H_1$, where $A$ is stable). Hence $A' = A$.

The same sign condition makes a system with an input monotone (Angeli and Sontag, 2003). It is a statement about the
realization extended by its protocol: the couplings alone do not decide it.

## 4. Running the command

```bash
python3 -B -m fieldbridge memory hysterons examples/memory/hysterons/rfim_antiferromagnet.json --out-dir out/afm
```

The command reads the specification, predicts from the signs, then measures 24 subloops for each of four procedures:
least stable element first and random order, each with simple excursions and with excursions that contain a nested
subloop ($H_1 \to H_2 \to H_3 \to H_2 \to H_1$ with $H_2 < H_3 < H_1$). The excursions of one procedure start from
one rising branch, followed from saturation and stopped at each $H_1$ in turn. $H_1$ is drawn inside the switching
window of the rising branch of the major loop, between 10 and 90 per cent of its switches, and $H_2$ inside that of
the falling branch, below $H_1$: the descent reaches the part of the loop in which elements switch back. These
excursions are a sample. All of them start from the rising branch of the major loop, and an excursion that starts from
another turning point of a longer history can fail where none of the sampled ones does: a count of zero shows that the
sampled excursions return, not that every excursion does. Where the network extended by the drive is balanced, the
structure guarantees the return of every excursion. The report
`hysterons.md` states the prediction, the measured returns and the largest number of elements whose state at the
return differs from that at $H_1$. The figure `hysterons.png` draws the major loop with one excursion, the elements
unlike the state at $H_1$ along it, the signs of the couplings with the drive and the subloops that do not return;
the web page shows its first panel for each material.

The six examples in [`examples/memory/hysterons`](../../examples/memory/hysterons) come from three fields
(default seed, 96 subloops each):

| Example | Field | Loops of couplings | Couplings whose loop through the drive is frustrated | Return guaranteed | Subloops that did not return |
| --- | --- | --- | --- | --- | --- |
| [random-field Ising ferromagnet](../../examples/memory/hysterons/rfim_ferromagnet.json) | magnetism | none frustrated | 0 of 512 | yes | 0 |
| [antiferromagnet, uniform field](../../examples/memory/hysterons/rfim_antiferromagnet.json) | magnetism | none frustrated | 512 of 512 | no | 18 |
| [antiferromagnet, staggered field](../../examples/memory/hysterons/rfim_antiferromagnet_staggered.json) | magnetism | none frustrated | 0 of 512 | yes | 0 |
| [antiferromagnetic chain](../../examples/memory/hysterons/antiferromagnetic_chain.json) | magnetism | none frustrated | 256 of 256 | no | 0 |
| [independent pores](../../examples/memory/hysterons/adsorption_pores.json) | porous media | no couplings | none | yes | 0 |
| [soft spots under shear](../../examples/memory/hysterons/soft_spots.json) | amorphous solids | half of the plaquettes frustrated | 138 of 288 | no | 12 |

Every realization whose extended network is balanced returned in every subloop, for both orders of relaxation and
for nested excursions. The failures are small: at most three elements differ at the return. They show that the state
at the turning point is no longer a function of the turning points alone.

## 5. Frustration of the couplings is not the criterion

The antiferromagnet on the square lattice frustrates no plaquette: the lattice is bipartite, every loop of couplings
is even, and the product of signs around it is positive. In a uniform field, however, every coupling closes a
frustrated loop through the drive, and 18 of 96 subloops failed. The same couplings in a field of opposite signs on
the two sublattices returned in every subloop, because relabeling one sublattice turns them into the ferromagnet in a
uniform field. Conversely, since every element is driven, the product of $\eta_i J_{ij} \eta_j$ around a loop of
couplings equals the product of the $J_{ij}$: a frustrated loop of couplings always contains a coupling whose loop
through the drive is frustrated. The soft spots have both.

## 6. Sufficient, not necessary

A frustrated loop through the drive allows the return to fail. It does not force a failure. Random antiferromagnetic
chains started from a large field return exactly (Deutsch, Dhar and Narayan, 2004), although every loop through the
drive is frustrated and the standard argument does not apply, and the command reproduces this: 0 of 96 subloops
failed. Chains whose couplings or drive have mixed signs can fail (Exercise 4). The command therefore reports what the dynamics
does next to what the structure guarantees, and it never turns a frustrated loop into a predicted failure.

## 7. The balanced counterpart

A frustrated loop through the drive says only that the return is no longer guaranteed. A counterpart turns this into
a measurement. Reverse the smallest set of signs, couplings or actions of the drive, that makes the network with the
drive balanced; their number is the frustration index $F$ of the network that counts the drive as an element. Keep
everything else: the graph, the magnitudes of the couplings, the fields, the half-widths and the strength of the
drive. By Section 3 this counterpart returns to every turning point. Measured against it, the states of the material
separate into two parts,

$$
C - C_0 = \underbrace{(C_{\mathrm{bal}} - C_0)}_{\text{loss}} + \underbrace{(C - C_{\mathrm{bal}})}_{\text{gain}} , \qquad (2)
$$

where $C$ is $\log_2$ of the number of states that some history of the drive reaches from negative saturation,
$C_{\mathrm{bal}}$ the same for the counterpart and $C_0$ for the uncoupled elements. The loss is the locking of the
balanced network below its elements; the gain counts the states that exist only because some couplings cannot be
satisfied. States that need three or more reversals of the drive are reported separately as deep memory.

```bash
python3 -B -m fieldbridge memory counterpart examples/counterpart/frustrated_network.json --out-dir out/counterpart
```

The command computes $F$ over all $2^n$ relabelings up to 16 elements (above that, the best of local searches, an
upper bound), builds the counterpart, and enumerates the transition graph of each material: its stable states and the
avalanche that leads from each state when the drive leaves the interval on which it is stable. The return is then
tested from every reachable state, not from a sample: every state on a rising or falling branch is a turning point,
and every excursion that does not pass the previous turning point is followed back. The number of states grows
quickly with the number of elements; above `--max-states` the enumeration stops and the report says so. For the
example, ten hysterons with competing couplings:

| Material | States | Bits | Return test: failed excursions |
| --- | --- | --- | --- |
| network ($F = 8$) | 20 | 4.32 | 65 of 425 |
| balanced counterpart | 5 | 2.32 | 0 of 7 |
| uncoupled elements | 31 | 4.95 | not tested |

The counterpart reverses six couplings and the drive on two elements. It returns in every excursion, as the structure
guarantees, while the material fails in 65 of 425 and reaches four times as many states: a gain of 2.00 bits against a
loss of 2.63 bits. The gain is a number an experiment can compare, with the same elements wired as the counterpart as
the control.

The balanced counterpart depends on which minimizing relabeling is taken: different minimizers reverse different
edges, and on random networks of up to 16 elements the counts of two minimizers' counterparts differ in about one
case in five, by up to 1.5 bit. The report therefore also measures the material against its ferromagnetic reference,
the same graph, magnitudes, fields and half-widths with every coupling positive and every element driven upward.
The reference is unique and balanced, and the identities $C = C_{\mathrm{ref}} + \text{gain}$ and
$C_{\mathrm{ref}} = C_0 + \text{loss}$ define the gain and the loss that the report quotes first. Beside the measured
gain the report prints the typical gain of a random network of the same size, mean degree and frustration index in
the ensemble of the sequence-memory study (random graphs, thresholds spread over 0.7, couplings of 0.3 to 0.9):
about 0.11 to 0.17 bit per unsatisfiable edge above a threshold that grows with the degree. That number is a
comparison for networks that resemble the ensemble, not a prediction for a lattice or a designed network.

```bash
python3 -B -m fieldbridge memory predictions examples/memory/hysterons --out-dir out/memory_predictions
```

lists, for every specification of kind `hysterons` in a directory, the frustration index with the drive, whether the
network is balanced, the typical gain of the ensemble, and, up to `--enumerate-up-to` elements, the measured gain,
deep memory and return test. The rule is sufficient and not necessary: in the one-dimensional antiferromagnetic
chain of Deutsch, Dhar and Narayan every negative coupling closes an unsatisfiable triangle with the drive, yet three
of four random chains of 12 elements return in every excursion (one fails 13 % of them), and all hold more states than
their reference without any deep memory, because a chain has no avalanches to lock. The gain over the reference
therefore has two sources, the absence of locking and the states beyond return-point memory, and the return test
tells them apart. `examples/counterpart/platforms_R8.md` lists the published platforms with the measurements that
would test the rule's consequence on each.

The same construction applies to any rule that decides a property from structure: the counterpart is the smallest
change of one layer that makes the rule hold, everything else is kept, and the consequence of the violation is
measured against it. The return to a turning point is the first rule for which FieldBridge builds it.

## 8. Materials given by equations

A material is often known by its equations of motion rather than by its thresholds: coupled flux loops,
magnetic grains, bistable chemical reactors. `memory hysterons` also accepts such a material: a specification of kind
`equations` (Module [2](16_memory_specification.md)) with a block `slow_drive`.

| Key | Meaning |
| --- | --- |
| `parameter` | the declared parameter through which the slow drive acts |
| `range` | `[low, high]`: a range of the drive with a single stable state at each end |
| `elements` | one expression per element whose sign is its state; for coupled bistable units, a coordinate of each unit |

The command follows the equations under a quasi-static drive, without noise. Each stable equilibrium is continued
in the drive until it disappears; past that point the state relaxes along the drift to another equilibrium. The
states that some history reaches from the low end of the range, with the intervals of the drive on which they exist
and the transitions at their ends, form the transition graph. The exhaustive return test of §7 runs on this graph.
Its result is a statement about the equations, with no model of elements in between.

The elements are needed only for the reduction to hysterons. At the end of a state's interval, the element whose
sign changes first during the relaxation is the one that lost stability. The drive at that point gives one linear
condition on the fields, half-widths and couplings of a hysteron model. The command finds every stable state on a
grid of the drive, reached or not, and solves the conditions of all of them by least squares. The signs of the
reduced couplings and of the drive then decide whether the return is guaranteed (§3). The reduced model is also
measured, and its states and its return are compared with those of the equations.

```bash
python3 -B -m fieldbridge memory hysterons examples/slow_drive/antiferromagnetic_pair.json --out-dir out/slow_drive
```

The example has two bistable coordinates, $\dot x_i = x_i - x_i^3 + J x_j + H \pm a$, with $J = -0.1$ and
$a = 0.05$. They switch one at a time, through the states $--$, $+-$ and $++$. The state $-+$ exists but no
history of the drive reaches it. Both excursions return. The reduced model has the coupling $-0.100$, the
fields $\pm 0.053$ and the half-widths $0.406$, in units of the drive, and the same three states. Its network has a
frustrated loop through the drive, which allows a failure; here none occurs. With $J = +0.1$, one avalanche switches
both coordinates, and the reduced coupling $+0.092$ is balanced with the drive.

The construction stops with its reason wherever the quasi-static dynamics is not that of hysterons:
- the range does not saturate the material;
- after a loss of stability the state oscillates or escapes;
- the drive does not decide the transition.

The last case occurs for identical elements. With $a = 0$, the symmetric state $--$ of the antiferromagnetic pair loses
stability through a pitchfork in the antisymmetric direction, so noise chooses which coordinate switches. A spread
of the elements' parameters lifts the degeneracy, as fabrication does in devices.

**A frustrated triangle of flux qubits** (`examples/slow_drive/rfsquid_triangle_unequal.json`). Three rf-SQUID loops of an annealing processor (Harris et al., 2010) are
coupled antiferromagnetically through their mutual inductance, and the flux bias $\varepsilon$ is swept slowly. Take
barrier parameters $\beta = (1.85, 1.46, 1.39)$, couplings $J_{12}, J_{13}, J_{23} = 0.05, 0.15, 0.26$ and flux
offsets $(0.14, -0.14, -0.13)$. The equations then reach six states. Two of nine excursions do not return to their
turning point. Rising from full negative flux, the triangle reaches $+-+$. The drive turns at $\varepsilon = 1.03$,
falls past $-0.37$ and rises again; the third loop switches back down, and on the way up the second loop switches
first. The triangle returns as $++-$: the excursion has changed which loop of the frustrated triangle is left
down. A direct integration of the equations under a slow ramp gives the same sequence of states.

**When the reduced model is wrong.** In the second example, `rfsquid_triangle_close_thresholds.json`, the equations
reach six states and every excursion returns. The hysteron model fitted to the same equations has the same six states,
but one transition is wrong.
- The first loop switches first only at folds where the other two loops have the same sign. Those folds fix only the
  sum of its two couplings; the report marks both as "not determined", and the fit takes them equal.
- In state $--+$ the two couplings act with opposite signs, so their unknown difference moves the first loop's
  threshold. The reduced model puts it 0.006 below that of the second loop.
- The equations switch the second loop there, the reduced model the first. Two of the reduced model's excursions then
  fail to return.

A coupling that the folds do not determine leaves the structural rule undecided, and the command says so. The return
is then decided on the transition graph of the equations, as here.

## 9. A model from your field

A specification of kind `hysterons` has the usual fields (`schema`, `kind`, `id`, `name`, `field`, `question`,
`assumptions`, `provenance`) and a `hysterons` object:

| Key | Meaning |
| --- | --- |
| `lattice` | `{"shape": "square" or "chain", "L": n}`: periodic nearest-neighbour bonds; or, instead, |
| `units`, `couplings` | the number of elements and a list of couplings `[i, j, J]` (reciprocal) or, with `"directed": true`, `{"from": j, "to": i, "J": J}` |
| `couplings` (lattice) | `{"value": J, "negative_fraction": p}`: bonds of strength J, a fraction p of them reversed at random |
| `fields` | the random fields $h_i$: a number, a list, `{"distribution": "gaussian", "width": w}` or `{"distribution": "uniform", "low": a, "high": b}` |
| `half_widths` | the half-widths $b_i \ge 0$ of the elements' own loops, in the same forms |
| `drive` | `"uniform"`, `"staggered"` (lattices), `"random"`, or a list of $\eta_i$ |
| `seed` | the seed of the random draws |

Write first the question that the model answers, then the assumptions that make the elements hysterons: which
barriers are large compared with the thermal energy, and on which time scale the drive is slow.

## 10. Exercises

1. Reverse the response of one element of the ferromagnet to the drive, with its couplings unchanged. How many loops
   through the drive become frustrated on the square lattice? Then reverse its four couplings as well. What does
   `memory hysterons` predict now, and why?
2. Couplings that are not reciprocal can make an avalanche cycle at fixed drive. Construct two elements for which this
   happens, and show that it cannot happen when the extended network is balanced.
3. In the soft-spot example, give every soft spot the same orientation to the strain (`"drive": "uniform"`). Is the
   return now guaranteed? Explain the answer from the plaquettes.
4. Copy the chain example, give it 512 elements, positive couplings (`"negative_fraction": 0`) and a drive of random
   sign (`"drive": "random"`), and run `memory hysterons` with its default seed. How many subloops fail, and how does
   this differ from the antiferromagnetic chain?
5. Run `memory counterpart` on the example with `--exact-up-to 0`, so that the frustration index comes from local
   search. Does the bound agree with the exact index? Why can a local search never report fewer unsatisfiable edges
   than the index?

<details>
<summary>Answers</summary>

1. Four, one for each coupling of the element. Reversing the couplings as well is the relabeling $s_0 \to -s_0$: the
   extended network is balanced again, and the return is guaranteed.
2. Let element 0 raise element 1 ($J_{10} > 0$) and element 1 lower element 0 ($J_{01} < 0$). At a drive where 0
   switches up, 1 follows, 0 switches back, 1 follows it down, and the cycle repeats (`test_an_avalanche_that_does_not_end_is_reported`).
   After a relabeling to non-negative couplings, a rising drive makes only upward switches, and each element switches
   at most once.
3. No. Half of the plaquettes are frustrated, and a frustrated loop of couplings always contains a coupling whose
   loop through the drive is frustrated, whatever the drive.
4. 4 of 96. About half of the couplings (240 of 512) close a frustrated loop through the drive, where in the
   antiferromagnetic chain every one does, yet here some subloops fail: the exact return of the antiferromagnetic
   chain is a property of that chain, not of chains in general.
5. Yes, both give $F = 8$ for this network. Every relabeling that a search visits leaves some number of unsatisfiable
   edges, and the index is the minimum of that number over all relabelings, so no search can report less; it can
   report more when it stops in a local minimum, which is why the report marks the value as an upper bound.

</details>

## Summary

- Interacting hysterons under a slow drive are specified by their couplings, random fields, half-widths and the
  signs with which the drive acts on them.
- If no loop of the network that counts the drive as an element is frustrated, a relabeling makes every coupling
  cooperative and the drive uniform, and the return to every turning point is exact.
- The frustration of the couplings alone does not decide the return: the antiferromagnet on a bipartite lattice fails
  in a uniform field and returns in a staggered one.
- The condition is sufficient and not necessary: antiferromagnetic chains return although every loop through the
  drive is frustrated.
- A material given by equations is followed under the slow drive itself. The return is decided on its transition graph, and its reduction to hysterons gives the signs for the structural rule. Where the folds do not determine a coupling, the reduction can get a transition wrong, and the return is decided on the equations.
- The balanced counterpart reverses the fewest signs that make the extended network balanced and keeps everything
  else; it returns to every turning point, and the states of the material are measured against it as a loss and a
  gain.

## Reference

| Result | Function | Test |
| --- | --- | --- |
| Specification and model | [`hysterons.from_spec`](../../fieldbridge/memory/hysterons.py), `Hysterons` | `test_specifications_outside_the_contract_are_refused` |
| Prediction from structure | `hysterons.predict`, `relabel`, `canonical_form` | `test_a_network_cooperative_after_relabeling_returns_to_every_turning_point`, `test_the_relabeling_maps_the_dynamics_onto_the_cooperative_form` |
| Dynamics and subloops | `hysterons.Run`, `subloops`, `check` | `test_a_frustrated_loop_through_the_drive_can_break_the_return_and_need_not` |
| Card and figure | `hysterons.card`, [`visual.hysteron_figure`](../../fieldbridge/memory/visual.py) | `test_the_card_of_a_material_draws_its_loop_its_excursion_and_the_signs_with_the_drive` |
| Report | [`cli.cmd_hysterons`](../../fieldbridge/memory/cli.py) | `test_command_line_writes_a_report_with_provenance` |
| Frustration index and counterpart | [`counterpart.frustration`](../../fieldbridge/memory/counterpart.py), `balanced_counterpart` | `test_the_counterpart_keeps_magnitudes_fields_and_thresholds_and_is_balanced`, `test_local_search_never_undercuts_the_exact_index` |
| States, deep memory and exhaustive return | `counterpart.reachable`, `reversals`, `exact_return`, `measure` | `test_the_enumeration_reproduces_the_preisach_counts`, `test_the_dynamics_of_run_stays_inside_the_enumerated_states` |
| Materials given by equations: transition graph, return and reduction | [`slow_drive.EquationGraph`](../../fieldbridge/memory/slow_drive.py), `transition_graph`, `survey`, `reduce`, `check` | `test_cooperative_pair_reduces_to_a_balanced_network`, `test_antiferromagnetic_pair_has_a_frustrated_loop_through_the_drive`, `test_identical_antiferromagnetic_elements_are_a_choice_the_drive_does_not_make` |
| Loss and gain against the counterpart | `counterpart.construct`, [`cli.cmd_counterpart`](../../fieldbridge/memory/cli.py) | `test_the_counterpart_returns_to_every_turning_point_and_the_parts_add_up`, `test_the_example_gains_states_over_a_counterpart_that_restores_the_return` |

Sources: F. Preisach, Z. Phys. 94, 277 (1935); F. Harary, Michigan Math. J. 2, 143 (1953); I. D. Mayergoyz,
*Mathematical Models of Hysteresis* (Springer, 1991); A. A. Middleton, Phys. Rev. Lett. 68, 670 (1992);
J. P. Sethna, K. Dahmen, S. Kartha, J. A. Krumhansl, B. W. Roberts and J. D. Shore, Phys. Rev. Lett. 70, 3347
(1993); D. Angeli and E. D. Sontag, IEEE Trans. Automat. Control 48, 1684 (2003); J. M. Deutsch, A. Dhar and
O. Narayan, Phys. Rev. Lett. 92, 227203 (2004); N. C. Keim, J. D. Paulsen, Z. Zeravcic, S. Sastry and S. R. Nagel,
Rev. Mod. Phys. 91, 035002 (2019); M. van Hecke, Phys. Rev. E 104, 054608 (2021).

[Previous: Module 11](26_memory_phase_locking.md) · [Tutorial index](index.md) · [Glossary](memory_glossary.md)
