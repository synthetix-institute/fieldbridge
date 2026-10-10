# Decision, Module 1. When a population acts as one

**Learning objectives.** After this module you can

1. describe a population in a specification of schema `fieldbridge-decision/1`: limit-cycle units with a spread
   parameter (`oscillators`), the mean-field equations of a population with a noise that falls as 1/N (`collective`),
   or units with a diverse parameter coupled through their mean (`units`);
2. compute, from the equations, the coupling at which limit-cycle units synchronize and the rate at which synchrony
   grows above it;
3. compute the probability that a population swept through its collective threshold follows a weak bias, and say
   when the composition of a random sample sets the outcome instead;
4. compute the window of passage times after a step that makes the collective state unstable, from the rate of the
   unstable mode along the mean trajectory;
5. run `decision predict`, `decision card` and `decision survey`, and read the class and the laws.

**Prerequisites.** Modules [1](15_memory_first_card.md), [2](16_memory_specification.md) and
[4](18_memory_writing_and_retention.md) for the realization, its specification and the swept write;
[Module 11](26_memory_phase_locking.md) for the phase response of a limit cycle. **Time.** About 60 minutes.

## 1. The question

A population acts as one when its units, each of which could do something else, do the same thing. Oscillators with
different natural frequencies lock to a common rhythm; a group of animals, spins, molecules or neurons chooses one of
two collective states; after a cue the whole group leaves the state it was in. Three questions follow. At which
coupling does a population synchronize? How accurately does it follow a weak bias when it is driven through its
collective threshold, for N units and a given speed? How long, and with what spread, does it take to leave a state
that a step has made unstable? This module answers them from the equations of the units, in published models from
electronics, neuroscience, chemistry, chronobiology, ecology, statistical physics, animal behaviour, photonics,
nonlinear dynamics and magnetism.

## 2. Synchronization

Limit-cycle units dx<sub>i</sub>/dt = F(x<sub>i</sub>; p<sub>i</sub>) with a spread of the parameter p are coupled
through an observable c(x) acting on one variable, K(c̄ − c(x<sub>i</sub>)). At weak coupling each unit is reduced to
its phase (Kuramoto 1984): the frequency map ω(p) turns the density of p into the density g(ω) of natural frequencies
(a change of variables, not a fit), and the phase response Z of the unit gives the coupling function

$$H(\phi) = \frac{1}{2\pi}\int_0^{2\pi} Z_k(\psi)\,\bigl[c(x(\psi+\phi)) - c(x(\psi))\bigr]\,d\psi
= a_1\cos\phi + b_1\sin\phi + \dots$$

The incoherent state of the phases loses stability where a perturbation e<sup>λt</sup> first grows,

$$1 = \frac{K}{2}\,(b_1 - i a_1)\int \frac{g(\omega)}{\lambda + i\omega}\,d\omega,$$

which at the onset (λ → −iΩ + 0) gives two real conditions for K<sub>c</sub> and the frequency Ω of the emerging
rhythm. For a<sub>1</sub> = 0 and a symmetric unimodal g, K<sub>c</sub> = 2/(π g(0) b<sub>1</sub>) (Kuramoto 1975);
a<sub>1</sub> ≠ 0 shifts K<sub>c</sub> and Ω (Sakaguchi and Kuramoto 1986). Above K<sub>c</sub> the coherent
perturbation grows at the rate μ, the real part of the unstable root (Strogatz and Mirollo 1991). K<sub>c</sub> itself
is textbook; what the card computes is the reduction of the full unit to g(ω) and H, and with `--simulate` it checks
the growth of synchrony of the full units against μ.

```bash
python -m fieldbridge decision predict examples/decision/synchronization/van_der_pol.json
```

For van der Pol units at μ = 1 with a 1% spread of frequencies, coupled through y on y, the reduction gives the
frequency 0.943, the Floquet rate κ = 1.06, a<sub>1</sub> = −0.093 and b<sub>1</sub> = 0.499; the incoherent state
loses stability at K<sub>c</sub> = 0.0301 with Ω = 0.941, and at 1.5 K<sub>c</sub> synchrony grows at
μ = 5.44·10⁻³.

The full units start on the cycles of the unit forced by the self-consistent mean field of the incoherent state,
each group of equal parameters spread evenly in time, so that the mean field is constant until a seed of 10⁻⁷ grows.
The measured rate is that of the amplitude of the mean field at the reference frequency, period by period. It exceeds
μ by an amount of order K/κ, the correction from the amplitude of the units, and tends to μ as the spread, and with it
K/κ, goes to zero.

## 3. The swept collective write

A population whose units choose between two states, with the symmetry that exchanges them, has a collective
pitchfork: at the control c* the symmetric state loses stability. With the critical vectors v and w of the mean-field
equations, the rate a = dλ/dc, the bias projected on the critical mode h<sub>s</sub>, and the noise of the critical
mode D<sub>s</sub> = w<sup>T</sup>Dw, a sweep dc/dt = r selects the state favoured by the bias h with

$$P = \Phi\!\left(\frac{\pi^{1/4}\,h_s\,h}{D_s^{1/2}\,(a r)^{1/4}}\right),$$

the law of the symmetric write of the memory module. In a population the noise of the collective mode comes from its
units and falls as 1/N, so the probit grows as h N<sup>1/2</sup> r<sup>−1/4</sup>. The noise is computed from the
population itself: for units that move between states it is the chemical Langevin noise of the transitions,
D = (1/2N) Σ ν ν<sup>T</sup> a. The law needs Λ = a r/(|b| D<sub>s</sub>) ≫ 1, b the cubic coefficient.

The closed form is the limit of a slow sweep. The card also gives the law along the actual passage: the critical
eigenvalue, the bias and the noise integrated along the mean trajectory of the sweep, from the biased equilibrium at
the start. The two differ when the crossing window (r/a)<sup>1/2</sup> is not narrow against the range on which these
change, or when the mean state lags the control because the symmetric modes relax slowly.

```bash
python -m fieldbridge decision card examples/decision/write/decision_network_wong_wang2006.json --out-dir out/ww
```

In the decision network of two neural pools of Wong and Wang (2006) the spontaneous state loses stability at the
stimulus μ<sub>0</sub>* = 10.68 Hz in a subcritical pitchfork (b = −101). Ramped at 5 Hz/s, the symmetric mode relaxes
at only about 4/s, the mean state lags the ramp by about 1 Hz, and the antisymmetric mode crosses later, where the
bias (proportional to μ<sub>0</sub>) is larger: at the bias z = 0.4 in units of the closed form's spread, the law along
the passage gives P = 0.670 against the closed form's 0.655.

## 4. Passage from a seed

A step of the control, or a cue beyond the coercive field, makes a mode of the collective state unstable. The mode has
d components: d = 1 when one real eigenvalue leads (its symmetry is a reflection: the difference of two pools, an
excess of one enantiomer, the in-phase amplitude of parametric oscillators), d = 2 when a complex pair leads (a
rotation: the transverse magnetization of a macrospin). Along the mean trajectory the mode grows as
s(t) = e<sup>Λ(t)</sup>(ξ + …), Λ = ∫λ dt, from a Gaussian seed ξ of d components whose spread falls as N<sup>−1/2</sup>
(the stationary fluctuations before the step and the noise of the growth). A threshold far from the start is crossed
when Λ(τ) = ln(L/|ξ|), so the 10–90% window of passage times is

$$\Lambda(\tau_{90}) - \Lambda(\tau_{10}) = \ln q_d,\qquad q_d = \frac{\chi_d^{-1}(0.9)}{\chi_d^{-1}(0.1)},\quad
q_1 = 13.09,\ q_2 = 4.675,$$

whatever the noise and the threshold, and Λ(τ<sub>50</sub>) rises by 1/2 per factor e in N. For a constant rate the
window is λΔτ = ln q<sub>d</sub>. These passage statistics are known (Haake, Haus and Glauber 1981); the card takes the
rate from the linearization of the model along its own mean trajectory, with nothing fitted. The statistics need the
seed to lie far below the threshold; in practice Λ(τ<sub>50</sub>) ≳ 2.

```bash
python -m fieldbridge decision predict examples/decision/passage/coherent_ising_machine_step.json
```

For parametric oscillators pumped from p* − 0.5 to p* + 0.4 the in-phase mode grows at λ = 0.4 (d = 1): the window is
ln 13.09/0.4 = 6.4 in time. For a macrospin whose field is reversed to twice the anisotropy field, the transverse
magnetization rotates while it grows at λ = α(h − 1) = 0.1 (d = 2): the window is ln 4.675/0.1 = 15.4.

A bias gives the seed a mean z (in units of its spread): the population then leaves its state towards the favoured
side with Φ(z), and the window narrows to the folded-normal ratio of the seed magnitudes.

## 5. Running the command

A specification has the header of every example, a `kind`, the equations payload of the memory module where the kind
needs one (`carrier`, `parameters`, `drift`), a block `population` and a block `protocol`:

| Kind | `population` | Target |
| --- | --- | --- |
| `oscillators` | `heterogeneity` {parameter, density: gaussian, centre, spread: the frequency spread σ<sub>ω</sub>/ω}; `coupling` {through: c(x), acting_on, factors: [K/K<sub>c</sub>, …]} | synchronization |
| `collective` | `size` (a declared parameter: N), `control` {name, range}, `bias` {name}, `noise` {transitions: [{change, rate}]} or {variance: {variable: expr}}, `symmetric_state` | collective-write, seeded-passage |
| `units` | `unit` {variable, mean, drift}, `heterogeneity` {parameter, density: gaussian, centre, width}, `control`, `bias`, `noise` (a number D) | collective-write |

The protocol of a write is `{"sweep": {"from", "to", "rates"}, "sizes", "z", "sample"}` (z are biases in units of the
spread of the closed form; `sample` is `quantiles` or `random` for units); that of a passage is
`{"step": {"from", "to"}, "threshold", "duration", "sizes", "dt"}`. A collective may omit its drift when it lists its
transitions: the drift is then F = Σ ν a.

```bash
python -m fieldbridge decision card examples/decision/write/ising_glauber1963.json --simulate --replicas 2000 --out-dir out/ising
python -m fieldbridge decision survey examples/decision
```

With `--simulate` the card runs the stochastic check of its target: full units for synchronization; replicas of the
mean-field equations with their noise (or of the diverse units themselves) through the sweep, within
|P − P<sub>law</sub>| < 4·stderr + 0.02; passage times after the step at every size, the window in the growth exponent
against ln q<sub>d</sub> within 4·se + 3%, and the slope of the median against ln N. The collective checks run the
mean-field equations with a noise that falls as 1/N, not the individual-based models of the sources.

## 6. Controls

The first two controls change one component and change the class; the third is a contrast between two bodies after a
step:

| Control | Change | Class |
| --- | --- | --- |
| `controls/josephson_no_onset.json` | C: coupling through sin φ, whose first harmonic has no sine part (b<sub>1</sub> = 0) | synchronizes → no onset |
| `controls/diverse_units_random_sample.json` | P: a random sample of the diversity for every population | follows the bias → set by the sample |
| `passage/macrospin_stoner_wohlfarth.json` against `passage/coherent_ising_machine_step.json` | a complex pair instead of one real eigenvalue leads after the step | one-component seed → two-component seed |

An overdamped Josephson junction coupled through its supercurrent has a coupling function with a cosine part only
(a<sub>1</sub> = 0.382, b<sub>1</sub> = 0): the coupling shifts the frequencies and no coupling makes the incoherent
state unstable, as for the neutral in-phase state of a resistively loaded array (Tsang, Mirollo, Strogatz and Wiesenfeld
1991).

Bistable units with a Gaussian diversity (Tessone et al. 2006), swept in their coupling, follow the swept law when the
sample of the diversity is symmetric. A random sample is not: at the reference profile x*(a) of the infinite
population every unit feels the drift ∂f/∂X (X<sub>sample</sub> − X<sub>∞</sub>), so the sample carries a frozen bias
of standard deviation s<sub>q</sub> = |Σw ∂f/∂X / Σw ∂f/∂h| (Var<sub>a</sub> x*/N)<sup>1/2</sup>. At N = 200 it is
0.074 against the thermal spread 0.0021: the outcome is set by the sample, P = Φ(h/(σ<sub>th</sub><sup>2</sup> +
s<sub>q</sub><sup>2</sup>)<sup>1/2</sup>) hardly depends on the sweep rate, and a slower sweep does not average the
frozen bias out.

## 7. The published populations

```bash
python -m fieldbridge decision survey examples/decision
```

| Population | Field | Target | Class | From the equations |
| --- | --- | --- | --- | --- |
| Population of van der Pol oscillators | electronics | synchronization | synchronizes above K_c | a₁ = −0.0933, b₁ = 0.499; K_c = 0.0301, μ(1.5 K_c) = 0.00544 |
| Population of FitzHugh–Nagumo neurons | neuroscience | synchronization | synchronizes above K_c | a₁ = 0.0763, b₁ = 0.422; K_c = 0.00589, μ(1.5 K_c) = 0.000897 |
| Population of Brusselator reactors | chemical kinetics | synchronization | synchronizes above K_c | a₁ = −0.595, b₁ = 0.705; K_c = 0.0177, μ(1.5 K_c) = 0.00473 |
| Population of Goodwin clocks | chronobiology | synchronization | synchronizes above K_c | a₁ = −0.916, b₁ = 0.364; K_c = 0.0107, μ(1.5 K_c) = 0.00171 |
| Population of Rosenzweig–MacArthur predator–prey patches | ecology | synchronization | synchronizes above K_c | a₁ = −0.821, b₁ = 0.518; K_c = 0.00568, μ(1.5 K_c) = 0.00119 |
| Population of overdamped Josephson junctions coupled through the supercurrent | superconductivity | synchronization | no onset | a₁ = 0.382, b₁ = 0 (10⁻¹⁰ numerically); no K_c |
| Mean-field Ising model with Glauber dynamics | statistical physics | collective-write | follows the bias | threshold 1; P = 0.655 at z = 0.4, N = 200 (closed form 0.655); Λ = 12 |
| Honeybee nest-site choice with stop signals | animal behaviour | collective-write | follows the bias | threshold 0.8681; P = 0.648 at z = 0.4, N = 200 (closed form 0.655); Λ = 6.23 |
| Chiral autocatalysis with recycling | chemistry | collective-write | follows the bias | threshold 39.53; P = 0.657 at z = 0.4, N = 16 (closed form 0.655); Λ = 33.7 |
| Decision network of two neural pools | neuroscience | collective-write | follows the bias | threshold 10.68; P = 0.670 at z = 0.4, N = 4000 (closed form 0.655); Λ = 35.9 |
| Coherent Ising machine of degenerate parametric oscillators | photonics | collective-write | follows the bias | threshold 0.4; P = 0.655 at z = 0.4, N = 16 (closed form 0.655); Λ = 4.8e+03 |
| Globally coupled bistable units with Gaussian diversity | nonlinear dynamics | collective-write | follows the bias | threshold 1.358; P = 0.655 at z = 0.4, N = 200; frozen bias — against the thermal spread 0.0021 |
| Bistable units with Gaussian diversity, a random sample per population | nonlinear dynamics | collective-write | set by the sample | threshold 1.358; P = 0.509 at z = 0.8, N = 200; frozen bias 0.0742 against the thermal spread 0.0021 |
| Coherent Ising machine after a step of the pump | photonics | seeded-passage | one-component seed | d = 1; rate at the end 0.4; window ln 13.09 |
| Chiral autocatalysis after a step of the substrate | chemistry | seeded-passage | one-component seed | d = 1; rate at the end 0.4644; window ln 13.09 |
| Decision network after the stimulus is switched on | neuroscience | seeded-passage | one-component seed | d = 1; rate at the end 4.361; window ln 13.09 |
| Uniaxial macrospin after the field is reversed | magnetism | seeded-passage | two-component seed | d = 2; rate at the end 0.1; window ln 4.675 |

The stochastic checks of `decision card --simulate` (2000 replicas for the collectives, 400 for the diverse units):

Synchronization (full units at 1.5 K<sub>c</sub>, 16 parameter groups of 8 phases):

| Units | Field | K/κ at 1.5 K_c | μ (reduction) | μ (128 full units) | ratio |
| --- | --- | --- | --- | --- | --- |
| Population of van der Pol oscillators | electronics | 0.0426 | 0.005437 | 0.005844 | 1.075 |
| Population of FitzHugh–Nagumo neurons | neuroscience | 0.0176 | 0.0008972 | 0.000984 | 1.097 |
| Population of Brusselator reactors | chemical kinetics | 0.023 | 0.004727 | 0.005007 | 1.059 |
| Population of Goodwin clocks | chronobiology | 0.269 | 0.001706 | 0.002944 | 1.725 |
| Population of Rosenzweig–MacArthur predator–prey patches | ecology | 0.0905 | 0.001194 | 0.0014 | 1.172 |

The excess over μ is the amplitude correction of order K/κ: 6–10% where K/κ is 0.02–0.04, 17% for the predator–prey
patches (K/κ = 0.09) and 73% for the Goodwin clocks, whose amplitude relaxes slowly (κ = 0.06, K/κ = 0.27 even at a
spread of 0.25%). The reduction is the limit of a small spread: as the spread and with it K/κ go to zero, the ratio
tends to 1 (Exercise 2).

The swept write (2000 replicas of the collective, 400 populations of diverse units):

| Population | Field | Conditions within the gate | Largest |P − P_law| / stderr |
| --- | --- | --- | --- |
| Mean-field Ising model with Glauber dynamics | statistical physics | 18 of 18 (2000 replicas) | 2.6 |
| Honeybee nest-site choice with stop signals | animal behaviour | 18 of 18 (2000 replicas) | 2.6 |
| Chiral autocatalysis with recycling | chemistry | 18 of 18 (2000 replicas) | 1.6 |
| Decision network of two neural pools | neuroscience | 18 of 18 (2000 replicas) | 3.6 |
| Coherent Ising machine of degenerate parametric oscillators | photonics | 18 of 18 (2000 replicas) | 2.7 |
| Globally coupled bistable units with Gaussian diversity | nonlinear dynamics | 8 of 8 (400 replicas) | 3.8 |
| Bistable units with Gaussian diversity, a random sample per population | nonlinear dynamics | 6 of 6 (400 replicas) | 1.3 |

For the diverse units the thermal noise D = 0.01 also softens the units with a near zero, whose restoring force is weak;
this raises P slightly above the deterministic reduction (by a fraction of order D), the largest deviation in the table.
In the random samples the biases are 0.8, 17.5 and 35 thermal spreads: the frozen bias, 35 thermal spreads at any N
(both fall as N<sup>−1/2</sup>), sets P = Φ(h/(σ<sub>th</sub><sup>2</sup> + s<sub>q</sub><sup>2</sup>)<sup>1/2</sup>),
close to ½ at the smallest bias, where a symmetric sample would give 0.79.

The passage from a seed (8000 replicas at each size):

| Population | Field | d | Λ window / ln q_d at the three sizes | Λ(τ₅₀) | Λ(τ₅₀) per ln N |
| --- | --- | --- | --- | --- | --- |
| Coherent Ising machine after a step of the pump | photonics | 1 | 1.017, 0.988, 1.013 | 5.28, 6.30, 7.34 | 0.495 ± 0.005 |
| Chiral autocatalysis after a step of the substrate | chemistry | 1 | 0.974, 0.966, 1.004 | 2.15, 2.85, 3.55 | 0.506 ± 0.006 |
| Decision network after the stimulus is switched on | neuroscience | 1 | 0.974, 0.990, 1.000 | 1.63, 2.34, 3.06 | 0.513 ± 0.007 |
| Uniaxial macrospin after the field is reversed | magnetism | 2 | 0.992, 0.972, 0.970 | 1.82, 2.51, 3.16 | 0.483 ± 0.004 |

Every window lies within 3.5% of ln q<sub>d</sub>, with the rate taken from the linearization. The windows fall below
the law where Λ(τ<sub>50</sub>) is near 2: the seed is then not far below the threshold, and the noise of the growth has
not yet accumulated when the earliest replicas cross. In the decision network the rate must be taken along the mean
trajectory: with the final rate at 30 Hz, 4.36/s, the measured windows would be 16–20% short of ln 13.09 (0.84, 0.81 and
0.80 of it).

## 8. A model from your field

Write the equations of your units or of your population's mean field, name its size N as a parameter (for a
collective), its control and its bias, and give the interval in which the threshold lies. Run `decision predict`
first: it tells you where the symmetric state loses stability, whether the pitchfork is supercritical, how large Λ is
and, after a step, how many components the unstable mode has. Then run `decision card --simulate`. A specification in
`examples/decision` appears on the web page in the column of its class.

## 9. Exercises

1. Show that for a Lorentzian density of half-width γ and H = sin φ the growth rate above onset is μ = K/2 − γ, and
   check it with `onset.growth_rate`.
2. Halve the spread of the Goodwin clocks. How do K<sub>c</sub> and K/κ change, and what do you expect for the ratio
   of the measured growth rate to μ?
3. In the Ising example, raise N by a factor 16 at a fixed field. By how much does the probit of the favoured state
   grow, and by how much would a 16 times slower sweep raise it?
4. Ramp the stimulus of the decision network at 1.25 Hz/s instead of 5. Does the law along the passage approach the
   closed form? Why?
5. Reverse the field of the macrospin only to 1.5 times the anisotropy field. Which condition of the window law is then
   hardest to meet at Δ = 40?

## Summary

- Limit-cycle units reduced to the density of their frequencies and the first harmonic of their coupling function
  synchronize above the K<sub>c</sub> of the dispersion relation; a coupling with no sine part synchronizes no
  population.
- A population swept through its collective pitchfork follows a weak bias with
  P = Φ(π<sup>1/4</sup>h<sub>s</sub>h/(D<sub>s</sub><sup>1/2</sup>(ar)<sup>1/4</sup>)), D<sub>s</sub> ∝ 1/N; along the
  actual passage where the crossing window is wide or the mean state lags. A random sample of diverse units carries a
  frozen bias that can set the outcome whatever the sweep rate.
- After a step the unstable mode leaves its state within a window of ln 13.09 (one component) or ln 4.675 (two) in
  the growth exponent along the mean trajectory, and the median moves by 1/2 per factor e in N.

## Reference

`fieldbridge.decision`: `spec` (the schema and its kinds), `oscillators` (the reduction to g(ω) and H, the full
units), `onset` (the dispersion relation and the growth rate), `collective` (the reduction at the threshold, the laws
of the swept write, the stochastic collective, the passage after a step), `units` (diverse units and the frozen bias),
`passage` (the window laws, the fold normal form), `predict`, `card`, `cli`. Examples in `examples/decision`:
synchronization (five units), write (six populations), passage (four), controls (two). Sources: Kuramoto (1975, 1984),
Sakaguchi and Kuramoto (1986), Strogatz and Mirollo (1991), Winfree (1967), Glauber (1963), Pais et al. (2013), Saito,
Sugimori and Hyuga (2007), Wong and Wang (2006), Martí et al. (2008), Wang et al. (2013), Tessone et al. (2006), Haake,
Haus and Glauber (1981), Stoner and Wohlfarth (1948), Brown (1963), Kondepudi and Nelson (1983); see docs/mechanisms.md.
