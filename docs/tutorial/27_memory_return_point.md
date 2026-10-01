# Module 12. Return to a turning point

**Learning objectives.** After this module you can

1. describe a material as interacting hysterons under a slow drive, in a specification of kind `hysterons`;
2. state when such a material returns exactly to its state at a turning point of the drive, from the signs of its
   couplings and of the drive alone;
3. run `memory hysterons` and read the prediction against subloops of the dynamics;
4. explain why the frustration of the couplings does not decide the return, and why the structural condition is
   sufficient and not necessary.

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

The command reads the specification, predicts from the signs, then measures 24 subloops for each of four
procedures: least stable element first and random order, each with simple excursions and with excursions that contain
a nested subloop ($H_1 \to H_2 \to H_3 \to H_2 \to H_1$ with $H_2 < H_3 < H_1$). The excursions of one procedure
start from one rising branch, followed from saturation and stopped at each $H_1$ in turn. $H_1$ is drawn inside the
switching window of the rising branch of the major loop, between 10 and 90 per cent of its switches, and $H_2$ inside
that of the falling branch, below $H_1$: the descent reaches the part of the loop in which elements switch back. The report `hysterons.md` states the prediction, the measured returns and the largest number of
elements whose state at the return differs from that at $H_1$.

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

## 7. A model from your field

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

## 8. Exercises

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

## Reference

| Result | Function | Test |
| --- | --- | --- |
| Specification and model | [`hysterons.from_spec`](../../fieldbridge/memory/hysterons.py), `Hysterons` | `test_specifications_outside_the_contract_are_refused` |
| Prediction from structure | `hysterons.predict`, `relabel`, `canonical_form` | `test_a_network_cooperative_after_relabeling_returns_to_every_turning_point`, `test_the_relabeling_maps_the_dynamics_onto_the_cooperative_form` |
| Dynamics and subloops | `hysterons.Run`, `subloops`, `check` | `test_a_frustrated_loop_through_the_drive_can_break_the_return_and_need_not` |
| Report | [`cli.cmd_hysterons`](../../fieldbridge/memory/cli.py) | `test_command_line_writes_a_report_with_provenance` |

Sources: F. Preisach, Z. Phys. 94, 277 (1935); F. Harary, Michigan Math. J. 2, 143 (1953); I. D. Mayergoyz,
*Mathematical Models of Hysteresis* (Springer, 1991); A. A. Middleton, Phys. Rev. Lett. 68, 670 (1992);
J. P. Sethna, K. Dahmen, S. Kartha, J. A. Krumhansl, B. W. Roberts and J. D. Shore, Phys. Rev. Lett. 70, 3347
(1993); D. Angeli and E. D. Sontag, IEEE Trans. Automat. Control 48, 1684 (2003); J. M. Deutsch, A. Dhar and
O. Narayan, Phys. Rev. Lett. 92, 227203 (2004); N. C. Keim, J. D. Paulsen, Z. Zeravcic, S. Sastry and S. R. Nagel,
Rev. Mod. Phys. 91, 035002 (2019); M. van Hecke, Phys. Rev. E 104, 054608 (2021).

[Previous: Module 11](26_memory_phase_locking.md) · [Tutorial index](index.md) · [Glossary](memory_glossary.md)
