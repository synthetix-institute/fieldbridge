# Module 15. Mechanisms of quantum-computing hardware

**Learning objectives.** After this module you can

1. describe an element of quantum-computing hardware by the classical image of its Lindbladian, with the noise of
   its vacuum;
2. read which known laws of the memory language the hardware carries and where a derivation stops;
3. derive the law at an obstruction with `memory complete`, and read the intrinsic offset of a sweep through a cusp
   and the write into a critical plane;
4. find with `memory changes` which term of a device makes its mechanism.

**Prerequisites.** [Module 9](23_memory_codiscovery.md), [Module 10](25_memory_threshold_write.md),
[Module 13](33_memory_noise_ratio.md) and [Module 14](34_memory_changes.md); the
[quantum write](30_quantum_write.md) for the exact dynamics of one oscillator.
**Time.** About 45 minutes of reading; the full survey of the 31 models takes about two hours on eight cores.

## 1. The question

Quantum computers are built from driven, damped and nonlinear elements: resonators that measure qubits, parametric
oscillators, flux and phase qubits, trapped ions, Rydberg ensembles, optomechanical and magnonic elements. Many of them are
operated by sweeping a drive through an instability, so that a weak signal, often the state of a qubit, decides
which of two states the element ends in. Which mechanisms of the memory language does this hardware carry, and
which laws appear where a derivation stops?

## 2. The elements and their models

The elements were not chosen by hand. A census of reviews of superconducting, semiconductor, trapped-ion,
neutral-atom, photonic and optomechanical platforms lists 60 hardware elements (`examples/quantum/semiclassical/survey/census.json`). Of
these, 21 have a classical write, 15 have one in part, and 24 have none: gates, error-correcting codes and linear
amplifiers store or process information without a write. A search of the abstracts of quantum-computing papers added
two classes that the reviews omit, cavity magnonics and spin ensembles in a cavity. Every element ends as a model, as
a gap with its reason (equations not in the sources, a rate model, a chaotic image, a running state), or as an
element without a classical write.

Each of the 31 models in `examples/quantum/semiclassical/` is the classical image of the element's Lindbladian: the equations of
the mean field, with the noise of the vacuum of each loss channel, $D = \kappa(2\bar n + 1)/8$ per quadrature, and
the parameters of the source named in its provenance. The image holds when the stored states contain many photons;
it leaves out the quantum dynamics of a few photons, which [the quantum write](30_quantum_write.md) computes for one
oscillator.

## 3. Running the survey

```bash
python3 -B examples/quantum/semiclassical/survey/survey.py examples/quantum/semiclassical/jba_resonator.json --out results/hardware
python3 -B examples/quantum/semiclassical/survey/report.py --results results/hardware
```

For each model, `survey.py` runs the three derivations of Modules 9–11 with their law constants
(`memory codiscover`), the law at every obstruction (`memory complete`, with 8 000 samples of the linear stage and a
certificate of 4 000 trajectories per bias), and the single changes of Module 14 (`memory changes`). Each of the
three derivations is stopped after 30 minutes and then recorded as not derived. The results of the full survey are
saved in [examples/quantum/semiclassical/survey/results](../../examples/quantum/semiclassical/survey/results), each with the stamp of the code that
produced it, and `report.py` collects them into [REPORT.md](../../examples/quantum/semiclassical/survey/REPORT.md), whose table this module reads. The
[records](../records/quantum_semiclassical.md) of `memory record` measure the same quantities at another seed, and the laws at the obstructions with
fewer trajectories; this module says so where they are the source. Each step is also a command of its own:

```bash
python3 -m fieldbridge memory codiscover examples/quantum/semiclassical/jba_resonator.json --target threshold-write --out-dir results/jba
python3 -m fieldbridge memory complete examples/quantum/semiclassical/jba_resonator.json --out-dir results/jba_complete
python3 -m fieldbridge memory changes examples/quantum/semiclassical/cim_triangle.json --out-dir results/cim_changes
```

## 4. Known laws carried to the hardware

In most models a derivation of Modules 9–11 reaches its mechanism, and the law carries its constant. The symmetric
write is reached in 18 of the 31 models. In 15 of them its constant lies within three standard errors or 3 % of
$\pi^{1/4}$, among them the persistent-current flux qubit, the rf-SQUID annealer qubit, the quantum flux parametron
latch and the SNAIL parametric oscillator. These constants carry standard errors of 0.11 to 0.22 at the default 400
trajectories and 8 replicates, so the count cannot see a deviation of a few per cent. The threshold write is reached
in 23 models, and in 19 of them the derivation measures the Airy constant $1.019$. Three of the other four, the
Kerr-cat qubit and the pair and the triangle of Kerr parametric oscillators, carry that constant in the records: the
derivation starts from a random state, and from some starts the branch of its slowest sweep does not end at the fold.
In the fourth, the Jaynes-Cummings resonator at high power, no stable state exists at the start of that sweep.

Phase locking is reached in the three self-sustained oscillators: the trapped-ion phonon laser, the quantum van der
Pol oscillator and the optomechanical self-oscillator. The first two have the half-width $1.000$. In the third the
fit is not a measurement, because the slip frequency outside the locking range cannot be read and the fitted centre
lies outside the range, and no constant is recorded. The phase-locking derivations of the driven Jaynes-Cummings
resonator at high power and of the spin ensemble did not finish within their 30 minutes and are recorded as not
derived.

These results are expected. A universal law returns its constant on every carrier whose dynamics reduce to its
normal form; what the survey adds is the device's own scales, which set how weak a signal the element can write.

The law record (`docs/site/law_constants.json`), with 25,600 trajectories, resolves the symmetric write to about 0.015.
Thirteen of the 15 models stay within 2.2 standard errors of $\pi^{1/4}$ there, and two do not. The dissipative cat
qubit gives $1.217 \pm 0.017$, for a reason that is not identified. The rf-SQUID phase qubit gives $1.297 \pm 0.014$,
and six more runs of the same size bring its mean to $1.307 \pm 0.006$, which is 4.2 standard errors low. This
deviation comes from the way the constant is measured and not from the write. The check starts the element from rest,
with the bias and the noise on, a time $4/\sqrt{r} = 28$ before the crossing, and a phase particle with the damping
$\gamma = 0.2$ of this element rings for a time $2/\gamma = 10$, which is not short against that time. The linear stage
integrated exactly from this start gives 1.311, and from a start twice as early 1.324; the noise on the phase itself
accounts for the remaining 0.5 %. At the lower damping $\gamma = 0.1$ the check gives $1.154 \pm 0.008$ against 1.169
from the linear stage, the value does not change when the noise ratio falls from 0.05 to 0.0008, and in the full
dynamics a start twice as early restores $1.33$. At $\gamma \geq 0.5$ the start has no measurable effect
([phase_qubit_start.py](../../examples/quantum/semiclassical/survey/phase_qubit_start.py)). The phase qubit therefore
carries the law, and its recorded constant is low by a known amount.

## 5. Laws where a derivation stops

Where a derivation stops, `memory complete` keeps what the obstruction names and derives the law of the completed
model from the sweep itself, with a certificate from the full stochastic dynamics. Four kinds of law appear.

**A sweep through a cusp writes with an offset of its own.** The driven Kerr resonator of a Josephson bifurcation
amplifier, the strongly driven transmon and the driven trapped ion lose their state at a fold: the write is
one-sided. A second
parameter (the dispersive shift, the detuning or the anharmonicity) unfolds the fold at a cusp, and a sweep through
the cusp writes symmetrically. Its response to a bias keeps the constant of the pitchfork, 1.32 to 1.37 in the full
dynamics, but the sweep also writes with an intrinsic offset, in units of the smallest bias that the noise lets the
write resolve. Two write points of an element unfold at the same cusp and give the two numbers of a row:

| Element | Offset, full dynamics | Offset, reduced law |
| --- | ---: | ---: |
| Kerr resonator of a bifurcation amplifier | 0.77, 0.75 | 1.02, 1.03 |
| strongly driven transmon | 0.81, 0.81 | 0.98, 0.96 |
| driven trapped ion | 0.78, 0.79 | 0.95, 0.92 |

The three elements share the Kerr normal form and therefore the size of the offset and the state it favours: averaged
over the two signs of a bias of one noise floor, the state of smaller amplitude is written with probability 0.68 to
0.69 in all three. The report prints these offsets with a sign, which is only the orientation that the eigensolver
gave the critical direction and has no meaning of its own. When the resonator measures a qubit, the dispersive shift
of the qubit is the bias of this write, so the offset is a systematic error of the measurement that the shape of the
sweep must cancel. The reduced law overestimates the offset by about a quarter and is not certified at the strict
tolerance; the full dynamics are the reference. Other cusps have offsets of their own: 2.2 to 2.3 for the dynamically
polarised nuclear-spin bath, where the law is certified, and 0.11 for the Rydberg ensemble. The magnonic element has
an offset of 1.03 toward its state of smaller amplitude, and its reduced law fails: it gives a response of −0.08
where the full dynamics give 1.30. Both of its write points along the drive are folds and unfold at the same cusp. In
the optomechanical element the law departs from the full dynamics by up to 0.30, and 13 % of the written states
change during a hold, so only the full dynamics give the write.

The flux-biased rf-SQUID phase qubit shows the case without an offset. Along the flux it loses a well at folds, where
$\cos\varphi = -\ell$, and it carries the threshold write there. The inductance ratio $\ell$ unfolds each fold at a
cusp at $\ell = 1$ and a flux of an odd multiple of $\pi$. The reflection of the potential about that flux makes the
sweep through the cusp a pitchfork, so the offset in the full dynamics stays at 0.02 or less, with a response of 1.31
to 1.34. The reduced law of one coordinate gives only 0.86 to 0.89 here and is not certified.

**The curvature of the sweep that cancels the offset.** `memory complete --curvature` curves the sweep through the
cusp, p = cusp + t d + c2 t² n with n normal to the sweep direction d in the plane of the two parameters, reads the
offset from the full dynamics at five curvatures, places its zero between the curvatures that bracket it, and refines
it by secant steps until the offset is below 0.1 or three standard errors:

```bash
python3 -m fieldbridge memory complete examples/quantum/semiclassical/jba_resonator.json --curvature --out-dir results/jba_curvature
```

For the resonator of the bifurcation amplifier the offset vanishes at c2 = −0.124, −0.122, −0.104, −0.070 and −0.022
for noise ratios Γ = 0.0125, 0.025, 0.05, 0.1 and 0.2, with the response staying at 1.27 to 1.40. The curvature has a
part that does not depend on the noise, the curvature of the line on which the two written states are mirror images,
and a part that grows with the noise ratio, which compensates the drift that the quadratic term of the reduced
equation draws from the fluctuations away from the cusp. The transmon and the trapped ion are cancelled at 0.081 to
0.085 and at 0.067 to 0.069 in the coordinates of their own second parameters, and the magnon element at 0.205, where
an offset of $0.06 \pm 0.03$ remains. For the nuclear-spin bath no curvature within reach cancels the offset: it is
still 0.62 at c2 = 0.38, and its zero lies beyond the curvatures that the sweep can reach. The equation of this
element near the cusp differs from the Kerr form. These runs are collected in
[curvature.json](../../examples/quantum/semiclassical/survey/curvature.json).

**Frustrated triangles write into a critical plane.** Three coherent-Ising-machine nodes, three Kerr parametric
oscillators or three rf-SQUID annealer qubits on an antiferromagnetic triangle have two directions that lose
stability together. The law of the plane, the linear-stage write applied in both directions and read as the sector
in which the state lands, gives the probabilities of the six frustrated configurations to within 0.024 of the full
dynamics on the coherent Ising machine, where it is certified, and to within 0.039 and 0.044 at the two upper write
points of the parametric-oscillator triangle, where these small differences still exceed the statistical error.

At the write point where the parametric-oscillator triangle leaves the state of zero amplitude, the law gives the six
configurations probabilities between 0.08 and 0.29 at a bias of half a noise floor, but the full dynamics end in none
of them: as the amplitude grows, the Kerr shift brings the uniform mode into resonance, and the network writes a state
that is not frustrated. The plane law holds where the shift moves every other mode away from resonance. On the
rf-SQUID triangle the law misses two uniform configurations, which are stable there and which the law of the plane
does not reach: the full dynamics end in one of them with a probability of up to 0.35, which is also the largest
deviation from the law. During a hold of 56.6 time units 20 % of the written states change in the survey and 8 % in
the record, on either side of the 10 % that separates an equilibrium write from the nonlinear stage, so the class of
this write is not settled.

**A jump to a distant branch keeps the law of the linear stage.** A parametric oscillator that measures a qubit
through its amplitude loses its state through a subcritical pitchfork and jumps to a distant branch. The choice is
made while the deviation still grows linearly, before the jump, so the response keeps the constant of the
pitchfork: 1.35 from the law and 1.32 from the full dynamics, certified. Module 13 explains why a first-order switch obeys the same law.

**Where the cubic term vanishes, the vacuum decides.** The Kerr-cat qubit at zero detuning, and the in-phase mode of
two coupled Kerr parametric oscillators, have no cubic term in the reduced equation, and the leading term is of fifth
order. The response to a bias then almost vanishes: the full dynamics give 0.08 and 0.07, and the law of the completed
model, which keeps the leading term, gives 0.04 and 0.00. The written states do not keep the choice either, since in
the [records](../records/quantum_semiclassical.md) about half of them change during a hold of 56.6 time units. Here the
classical image is not enough, and the quantum dynamics decide ([the quantum write](30_quantum_write.md)). The
phase-locked parametric oscillator used to measure qubits has a vanishing cubic term too, with the leading term
$-5.76 \times 10^{-4}\, s^5 = -K^2 s^5/\kappa$, but at many photons, where the written states are kept; there the
[record](../records/quantum_semiclassical.md) gives 0.78 from the law and 1.27 from the full dynamics.

In all, 17 laws at obstructions are certified and 25 are derived but not certified; the remaining
write points have no cusp, a continuous family of written states (the non-degenerate parametric converter) or an
oscillation.

## 6. Single changes of the hardware

The single changes of Module 14 show which term of a device makes its mechanism. In 10 of the 31 models a single
change reaches a write the device does not have.

- Removing one coupling of a frustrated triangle (coherent Ising machine, Kerr parametric oscillators, rf-SQUID
  qubits) reaches the symmetric write: without the symmetry of the triangle the two critical directions separate, and
  the first one to lose stability does so alone.
- Removing the loss of the parametric amplifier below threshold reaches the threshold write: the amplifier is then
  above its threshold.
- Removing a loss term or the coupling of two Kerr parametric oscillators reaches the symmetric write that the
  vanishing cubic of the in-phase mode prevents.
- Removing one coupling term of the non-degenerate parametric converter, whose written states form a continuous
  family of phases, reaches the symmetric write: the change breaks the phase symmetry that made every phase
  equivalent, and two states remain.
- Removing the Kerr term from one equation of the parametric oscillator that measures through its amplitude
  reaches the symmetric write.

## 7. What the survey does not compute

The models are classical images. Where a device holds only a few photons, the noise of the vacuum still acts after
the linear stage and the quantum dynamics can differ from the image; the Kerr-cat qubit at zero detuning is such a
case, and [the quantum write](30_quantum_write.md) treats it exactly. The survey computes laws of writing, not the
fidelity of a measurement or the weights with which an annealer samples its ground states. A law derived at an
obstruction is a calculated statement about the model, and whether it is new is decided by the literature of the
device.

## 8. Exercises

1. Choose an element of the survey that measures a qubit and convert its law into a measurement error for a given
   dispersive shift, sweep rate and loss.
2. Remove the coupling of one edge of the Kerr-parametric-oscillator triangle with `memory changes --save`, and run
   `memory complete` on the changed model. Which of the six configurations survive as written states?
3. Add an element of your own platform as a classical image, with the vacuum noise of its loss channels, and run the
   survey on it.
