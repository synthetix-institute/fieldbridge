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
**Time.** About 45 minutes of reading; the full survey of the 31 models takes about two hours on 16 cores.

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
saved in [examples/quantum/semiclassical/survey/results](../../examples/quantum/semiclassical/survey/results), and `report.py` collects them into
[REPORT.md](../../examples/quantum/semiclassical/survey/REPORT.md), whose table this module reads. Each step is also a command of its own:

```bash
python3 -m fieldbridge memory codiscover examples/quantum/semiclassical/jba_resonator.json --target threshold-write --out-dir results/jba
python3 -m fieldbridge memory complete examples/quantum/semiclassical/jba_resonator.json --out-dir results/jba_complete
python3 -m fieldbridge memory changes examples/quantum/semiclassical/cim_triangle.json --out-dir results/cim_changes
```

## 4. Known laws carried to the hardware

In most models a derivation of Modules 9–11 reaches its mechanism, and the law carries its constant. The symmetric
write is reached with a constant consistent with $\pi^{1/4}$ in 14 of the 31 models, among them the dissipative cat
qubit ($1.33 \pm 0.13$), the persistent-current flux qubit, the rf-SQUID annealer qubit, the quantum flux parametron
latch and the SNAIL parametric oscillator. The threshold write is reached with the Airy constant $1.019$ in
21 models. Phase locking is reached with the half-width $1.000$ in the three self-sustained oscillators: the
trapped-ion phonon laser, the quantum van der Pol oscillator and the optomechanical self-oscillator. The
phase-locking derivation of the driven Jaynes-Cummings resonator at high power did not finish within its 30 minutes
and is recorded as not derived. The constants of the symmetric write carry standard errors of 0.11 to 0.22 at the
default 400 trajectories and 8 replicates.

These results are expected. A universal law returns its constant on every carrier whose dynamics reduce to its
normal form; what the survey adds is the device's own scales, which set how weak a signal the element can write.

## 5. Laws where a derivation stops

Where a derivation stops, `memory complete` keeps what the obstruction names and derives the law of the completed
model from the sweep itself, with a certificate from the full stochastic dynamics. Four kinds of law appear.

**A sweep through a cusp writes with an offset of its own.** The driven Kerr resonator of a Josephson bifurcation
amplifier, the strongly driven transmon and the driven trapped ion lose their state at a fold: the write is
one-sided. A second
parameter (the dispersive shift, the detuning or the anharmonicity) unfolds the fold at a cusp, and a sweep through
the cusp writes symmetrically. Its response to a bias keeps the constant of the pitchfork, 1.33 to 1.35 in the full
dynamics, but the sweep also writes with an intrinsic offset, in units of the smallest bias that the noise lets the
write resolve:

| Element | Offset, full dynamics | Offset, reduced law |
| --- | ---: | ---: |
| Kerr resonator of a bifurcation amplifier | −0.80 | −1.03 |
| strongly driven transmon | +0.82, +0.84 | +1.02, +1.01 |
| driven trapped ion | −0.85, −0.83 | −0.99, −0.98 |

The three elements share the Kerr normal form and therefore the size of the offset; the side of the unfolding sets
its sign. When the resonator measures a qubit, the dispersive shift of the qubit is the bias of this write, so the
offset is a systematic error of the measurement that the shape of the sweep must cancel. The reduced law
overestimates the offset by about a quarter and is not certified at the strict tolerance; the full dynamics are the
reference. Other cusps have offsets
of their own: −2.3 for the dynamically polarised nuclear-spin bath, where the law is certified, and +0.08 for the
Rydberg ensemble. For the optomechanical and magnonic elements the reduced law fails, and only the full dynamics give
the write.

**Frustrated triangles write into a critical plane.** Three coherent-Ising-machine nodes, three Kerr parametric
oscillators or three rf-SQUID annealer qubits on an antiferromagnetic triangle have two directions that lose
stability together. The law of the plane, the linear-stage write applied in both directions and read as the sector
in which the state lands, gives the probabilities of the six frustrated configurations to within 0.038 to 0.056 of
the full dynamics. The exception is the write point at which the parametric-oscillator triangle leaves the state of
zero amplitude. There the law gives the six configurations probabilities between 0.09 and 0.28, but the full dynamics
end in none of them: as the amplitude grows, the Kerr shift brings the uniform mode into resonance, and the network
writes a state that is not frustrated. The plane law holds where the shift moves every other mode away from
resonance. The rf-SQUID phase qubit has a plane write with four states, certified.

**A jump to a distant branch keeps the law of the linear stage.** A parametric oscillator that measures a qubit
through its amplitude loses its state through a subcritical pitchfork and jumps to a distant branch. The choice is
made while the deviation still grows linearly, before the jump, so the response keeps the constant of the
pitchfork: 1.35 from the law and 1.33 from the full dynamics, certified. Module 13 explains why a first-order switch obeys the same law.

**Where the cubic term vanishes, the vacuum decides.** The Kerr-cat qubit at zero detuning, and the in-phase mode of
two coupled Kerr parametric oscillators, have no cubic term in the reduced equation, and the leading term is of fifth
order. The law of the completed model gives 1.35, but the full dynamics give 0.04 to 0.06: with the noise of the
vacuum the choice is not complete at the end of the linear stage. Here the classical image is not enough, and the
quantum dynamics decide ([the quantum write](30_quantum_write.md)). The phase-locked parametric oscillator used to
measure qubits has a degenerate leading order too, but at many photons, and its law is certified.

In all, 20 laws at obstructions are certified and 22 are derived but not certified; the remaining
write points have no cusp, no choice at the end of the sweep, a continuous family of written states (the
non-degenerate parametric converter) or an oscillation.

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
