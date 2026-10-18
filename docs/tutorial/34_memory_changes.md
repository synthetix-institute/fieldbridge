# Module 14. Changes of a model and their composition

**Learning objectives.** After this module you can

1. apply the three changes of a model that FieldBridge implements and read what each one shows;
2. tell a pair of changes that closes from one that is open, and say why an open pair matters;
3. state what a write reached after a change says about the supplied model, and what it does not.

**Prerequisites.** [Module 9](23_memory_codiscovery.md) and [Module 10](25_memory_threshold_write.md).
**Time.** About 40 minutes, most of it waiting for the pairs.

## 1. The question

Module 9 carries one mechanism to models from different fields and stops where a derivation fails. Papers suggest
other routes. Consecutive equations in a paper differ by a few kinds of change: a term is added, removed or modified,
a parameter becomes the control, the range or the order of a protocol changes, a closure or an observable is replaced.
Applied to a published model, each change gives a new model whose mechanism the derivations of Modules 9 and 10
decide, so a change is a derivation step that needs no target chosen in advance. This module asks which mechanisms
such changes reveal, and when two changes together do something that neither does alone.

## 2. Three changes and how they are judged

`memory changes` applies three kinds of change to an equation specification:

| Change | What it does |
| --- | --- |
| remove | drops one additive term of one equation |
| control | makes another parameter the control, swept from a quarter to four times its value |
| range | extends the range of the control by ten times its width below it, above it, or on both sides |

A changed model whose trajectories run away from large starts is discarded: the change removed the saturation, and
there is no stable state to write. Every other changed model is judged by the symmetric-write and threshold-write
derivations. Its outcome is the set of writes it reaches together with, for each write it does not reach, the class
of the reason why the derivation stops: no write point along the control, a single stable state, an oscillation, a
one-sided write that needs a second parameter, a jump to a distant branch, several critical directions, or a
vanishing cubic term. A change is a **step** when it reaches a write that the model did not reach or changes one of
these classes.

A **pair** of changes is a change of the changed model. It **closes** when its outcome is the outcome of the model or
of one of its single changes, or combines classes that single changes produce. It is **open** when it reaches a write
that no single change reaches, or a reason for stopping that no single change produces. A derivation that fails
numerically or exceeds its time limit decides nothing and does not open a pair.

## 3. Removing terms of a Hopf normal form

The normal form of a Hopf bifurcation, written as in a model of olfactory oscillators
([examples/changes/hopf_normal_form.json](../../examples/changes/hopf_normal_form.json)),

$$
\dot x=\omega(\mu x-y-x\rho^2),\qquad \dot y=\omega(\mu y+x-y\rho^2),\qquad \rho^2=x^2+y^2, \qquad (1)
$$

begins to oscillate at $\mu = 0$. Neither write is reached: the symmetric write stops at the oscillation, and the
threshold write finds a single stable state.

```bash
python3 -m fieldbridge memory changes examples/changes/hopf_normal_form.json --out-dir results/changes/hopf
```

Of twelve single changes, two are discarded because they remove the cubic saturation and the states run away. Of the
other ten, five are steps:

| Change | Symmetric write | Threshold write |
| --- | --- | --- |
| remove $-\omega y$ from $\dot x$ (or $\omega x$ from $\dot y$) | several critical directions | single stable state |
| remove $\omega\mu x$ from $\dot x$ | oscillation (no step) | single stable state |
| $\omega$ as the control | the control scales the drift | single stable state |
| range of $\mu$ above or below the published one | no write point | single stable state |

Removing one rotation term turns the rotation into a shear. The linear part becomes
$\omega\bigl(\begin{smallmatrix}\mu & 0\\ 1 & \mu\end{smallmatrix}\bigr)$, whose two eigenvalues are equal, so both
directions lose stability together at $\mu = 0$. Removing the growth term $\omega\mu x$ alone leaves a rotation and
the oscillation.

```bash
python3 -m fieldbridge memory changes examples/changes/hopf_normal_form.json --pairs --save \
    --out-dir results/changes/hopf_pairs
```

Of 71 pairs, 16 are open. Removing $-\omega y$ and $\omega\mu x$ together, in either order, reaches the symmetric
write and the threshold write, which no single change reaches. The changed model is

$$
\dot x=-\omega x\rho^2,\qquad \dot y=\omega(\mu y+x-y\rho^2). \qquad (2)
$$

The direction $x$ is now neutral at linear order and decays through the cubic term, and $y$ is the only direction
that loses stability. With $x = 0$, $\dot y = \omega(\mu y - y^3)$: a pitchfork with the two stable states
$y = \pm\sqrt{\mu}$, related by the symmetry $(x, y) \to (-x, -y)$ that Eq. (1) already had. The oscillator has become
a symmetric switch, and neither removal alone produces it. The other open pairs reach the threshold write, for
example removing both growth terms after the range of $\mu$ is extended.

## 4. Changing what is controlled: the Lorenz model and the disk dynamo

A published study varies one parameter over the range it needs. The transition that a model contains can lie
outside that range or along another parameter.

The Lorenz model in the chaotic range ([examples/changes/lorenz.json](../../examples/changes/lorenz.json)) is studied
for $24.74 \le r \le 60$, where no state loses stability. Extending the range of $r$ below the published one reaches
the symmetric write at $r = 1$, where the conduction state gives way to the two convection rolls,
$x = y = \pm\sqrt{b(r-1)}$, $z = r - 1$, related by the reflection $(x, y, z) \to (-x, -y, z)$ that reverses the
sense of rotation of the roll. It is the textbook pitchfork of the model, found by a change of the protocol.

```bash
python3 -m fieldbridge memory changes examples/changes/lorenz.json --out-dir results/changes/lorenz
```

The modified Bullard disk dynamo ([examples/changes/bullard_dynamo.json](../../examples/changes/bullard_dynamo.json))
couples two components of the magnetic flux, $x$ and $y$, to the angular velocity $Z$ of the disk. The control of its
source, the velocity $Z_0$ of the state without field, appears in the stability analysis of the paper but in none of
the three equations, where that velocity is $g/f$. No state of the specified model loses stability along $Z_0$: the symmetric write finds no write point. Two of the twelve
single changes are steps, and both make another parameter the control: the torque $g$ or the friction $f$. Either
choice reaches the symmetric write. The equations are unchanged by $(x, y) \to (-x, -y)$, the reversal of the field,
so at the threshold of the dynamo the state without field gives way to two states of opposite polarity, and a weak
bias selects one of them. None of the 76 pairs is open: the two steps already contain everything the pairs reach.

```bash
python3 -m fieldbridge memory changes examples/changes/bullard_dynamo.json --pairs --out-dir results/changes/dynamo
```

## 5. Removing a term of an epidemic model and of a frustrated triangle

Two symmetric populations with pairwise and three-body infection
([examples/changes/two_population_sis.json](../../examples/changes/two_population_sis.json)) have four coexisting
stable endemic states in their source. Along the published control, the pairwise transmission rate $\tau$ over
$[0, 0.1]$, no state loses stability symmetrically, and the symmetric write finds no write point. Eight of the 17
single changes are steps. Making another rate the control, a coupling between the populations ($d_2$, $d_3$, $d_4$)
or the three-body strength of one population ($h_1$, $h_2$), reaches the symmetric write: the choice between two of
the endemic states is then made by a weak bias at a point the published protocol never crosses. Removing the
three-body term of one population instead leaves a one-sided write. None of the 98 pairs is open.

```bash
python3 -m fieldbridge memory changes examples/changes/two_population_sis.json --pairs --first 6 \
    --out-dir results/changes/sis
```

Three degenerate parametric oscillators of a coherent Ising machine on an antiferromagnetic triangle
([examples/hardware/cim_triangle.json](../../examples/hardware/cim_triangle.json), Module 15) write into a
two-dimensional critical plane. The uniform mode stays stable, and the two modes of the frustrated pattern lose
stability together, so the symmetric write stops at several critical directions. Removing any one of the three
coupling terms reaches the symmetric write. The triangle then loses the symmetry that made the two directions
equivalent; one of them loses stability first, and it does so alone, through a pitchfork. The change shows what the
coupling does in the published device: it is the coupling of all three edges that holds the two directions together
and makes the write choose among six configurations instead of two. Of 72 pairs, 7 are open, and none of them reaches
a further write.

```bash
python3 -m fieldbridge memory changes examples/hardware/cim_triangle.json --pairs --out-dir results/changes/triangle
```

## 6. What a change shows and what it does not

A write reached after a change is a calculated property of the changed model, not of the supplied one: removing the
rotation term of Eq. (1) does not make a Hopf oscillator bistable. What the change does show is the role of the
term in the supplied model. The rotation term keeps the two directions of Eq. (1) from losing stability together;
the coupling of a frustrated triangle keeps its two critical directions equivalent. A step names the term or the
parameter that was changed, and `--save` writes the changed specification, so that every other command of the
memory modules can be run on it.

The changes are limited to three kinds, judged by two writes. A change that adds a term, replaces a closure or moves
the model to another carrier is not implemented, and the derivations run under a time limit (`--budget`, 60 s by
default), so a slow derivation decides nothing. An open pair is a statement about this set of changes and these two
derivations; whether its mechanism is new to physics is decided by the literature, not by the program.

## 7. Exercises

1. Run the single changes of the epidemic model with `--kinds control` and explain why another rate, and not the
   published one, carries a symmetric write.
2. The pair that removes $\omega x$ from $\dot y$ and $\omega\mu x$ from $\dot x$ also reaches both writes. Find its
   two stable states by hand and the symmetry that relates them.
3. Take a model of your own field and find a pair of changes that is open. Which term does it show to be necessary
   for the mechanism of the published model?

Tests: `tests/test_memory_changes.py`.
