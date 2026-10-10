# Module 11. Phase locking in eight oscillators

**Learning objectives.** After this module you can

1. reduce an oscillator under a weak periodic drive to the Adler equation for its phase;
2. run `memory codiscover --target phase-locking` and read the derivations in eight oscillators from eight fields;
3. explain how a symmetry of the cycle sets the ratio at which an oscillator locks to a drive;
4. test the locking law with relaxation rates inside the locking range and slip frequencies outside it.

**Prerequisites.** Modules [6](20_memory_phase.md) and [9](23_memory_codiscovery.md). **Time.** About 45 minutes.

## 1. The target

Module 6 showed that a limit cycle stores a phase. Time-translation symmetry makes the phase a flat direction: a
pulse shifts it, the shift persists, and noise spreads it by diffusion. The target of this module gives this flat direction
a restoring force. A weak periodic modulation of the control, $p(t) = p_1 + \varepsilon\cos\omega_f t$, acts on the
phase $\theta$ through the phase response along the drive, $Z_p(\theta) = Z(\theta)\cdot\partial F/\partial p$. Here
$Z$ is the gradient of the phase on the cycle: the periodic solution of the adjoint equation
$\dot Z = -J^{\mathsf T} Z$, normalized by $Z\cdot F = \omega_0$. Averaged over the drive, only the harmonic of $Z_p$
that matches the ratio $\omega_f \approx n\,\omega_0$ survives, and the phase difference $\psi = \theta - \omega_f t/n$
obeys the Adler equation (Adler, 1946)

$$
\dot\psi = \Delta\omega - K \sin\phi, \qquad \phi = n\psi - \phi_n - \tfrac{\pi}{2}, \qquad
K = \tfrac12\,\varepsilon\,\lvert Z_n\rvert, \qquad \Delta\omega = \omega_0 - \omega_f/n  \qquad (1)
$$

where $Z_n = a_n - i b_n$ is the $n$-th Fourier coefficient of $Z_p$ and $\phi_n = \mathrm{atan2}(b_n, a_n)$. In
the units $\tau = nKt$ and $\nu = \Delta\omega/K$ it reads $d\phi/d\tau = \nu - \sin\phi$. For $\lvert\nu\rvert < 1$
the phase difference relaxes to one of $n$ values, $2\pi/n$ apart, at the rate $\lambda = nK(1-\nu^2)^{1/2}$: the
drive writes the phase and restores it after a perturbation. For $\lvert\nu\rvert > 1$ the phase slips at the
frequency $\Omega = K(\nu^2-1)^{1/2}$.

## 2. The transformations

The letters act on the same components as in Module 9, with these meanings:

| Letter | Acts on | What is computed | The derivation stops if |
| --- | --- | --- | --- |
| S | $\Xi$, $\Omega$ | a symmetry of the drift that holds for every value of the control maps the cycle onto itself $1/m$ of a period later; then $Z_p(\theta + 2\pi/m) = Z_p(\theta)$, and the drive acts only through the harmonics $m, 2m, \dots$ | (not required) |
| C | $A$ | the control is moved into the range where the model oscillates: the middle of the widest interval of oscillating values on a grid of 11 | no value oscillates, or the control multiplies the whole drift |
| R | $\Xi$, $\Omega$ | the state is reduced to the phase; $Z$ follows from the adjoint equation, and the Floquet multipliers of the cycle are checked | a transverse multiplier equals 1 (a family of neutral cycles) |
| K | $\Xi$ | averaging over the drive gives Eq. (1), with $n$ the smallest harmonic present; the phase gained over one period in the driven model is compared with $TK\cos(n\psi - \phi_n)$ | no harmonic up to 8 |
| L | $P$, $R$ | relaxation rates inside the locking range and slip frequencies outside it, simulated at $\nu = 0, \pm 0.7, \pm 1.3, \pm 2$ and two drive amplitudes; the half-width of the locking range in units of $K$ is extrapolated to zero amplitude | — |

The drive is weak against both the frequency and the attraction of the cycle: $K = 0.02$ and $0.01$ times
$\min(\omega_0, 2\kappa)$, where $\kappa$ is the relaxation rate of the amplitude, the slowest transverse Floquet
exponent.

## 3. Running the constructor

```bash
python3 -B -m fieldbridge memory codiscover --target phase-locking --out-dir build/tut/m9_phase
```

Without arguments the command reads `examples/memory/oscillators` and `examples/memory`. The first directory holds
eight oscillators from seven fields; the second adds the ring of three repressors and the models of the write
targets. The command takes about 15 minutes on a laptop, or 2 minutes with `--no-law`.

| Realization | Field | Derivation | Operating point | Ratio | $K/\varepsilon$ | Half-width in units of $K$ |
| --- | --- | --- | --- | --- | --- | --- |
| van der Pol oscillator | electronics | RKL | bias $= 0$ | 1:1 | 0.2674 | 0.9998 ± 0.0015 |
| Brusselator | chemical kinetics | RKL | $b = 3$ | 1:1 | 0.5120 | 1.0000 ± 0.0025 |
| Goodwin clock | chronobiology | RKL | $\mathrm{tr} = 10$ | 1:1 | 0.1579 | 0.9998 ± 0.0042 |
| FitzHugh–Nagumo neuron | neuroscience | CRKL | current $= 0.7$ | 1:1 | 0.2128 | 0.9999 ± 0.0032 |
| Rosenzweig–MacArthur predator and prey | ecology | CRKL | capacity $= 3.25$ | 1:1 | 0.0591 | 0.9992 ± 0.0073 |
| overdamped Josephson junction | superconductivity | CRKL | current $= 1.6$ | 1:1 | 0.4003 | 0.9997 ± 0.0002 |
| pumped oscillator (parametron principle) | computing hardware | SRKL | $k = 1$ | 2:1 | 0.2542 | 0.9998 ± 0.0015 |
| ring of 3 repressors | synthetic biology | SRKL | $\alpha = 10$ | 3:1 | 0.0165 | 0.997 ± 0.019 |
| Lotka–Volterra predator and prey | ecology | R | stops at R: a family of neutral cycles | — | — | — |
| capillary rotors, dipoles | soft matter, magnetism | — | stop: the control only rescales time | — | — | — |
| three compartments | compartment models | — | stops: no control | — | — | — |
| models of the write targets | six fields | — | stop at C: no oscillation in the control range | — | — | — |

Eight models from eight fields reach the target by four classes of derivation. Three oscillate as specified (R K L).
Three must first be moved into oscillation (C R K L). The FitzHugh–Nagumo neuron is excitable without current and
fires periodically for currents between about 0.33 and 1.42; the constructor takes 0.7, the middle of the
oscillating grid values. The Rosenzweig–MacArthur equilibrium loses stability when the carrying capacity exceeds
$1 + 2x^*$, where $x^* = 2/3$ is the prey density at equilibrium (enrichment destabilizes the equilibrium; Rosenzweig,
1971). The junction carries no voltage below its critical current and rotates above it. Two models reach the target
through a symmetry and lock at a higher ratio (S R K L).

![Derivations of phase locking, the averaged drift of the phase and the locking law](figures/memory/m9_phase.png)

*Figure 1. Phase locking. (a) The derivation in each model; C moves the control into the oscillation, and S sets the
ratio. (b) The phase response along the drive, different in every model. (c) The phase gained per period in the
driven model, in units of $TK$, against $-\sin\phi$. (d) The relaxation rate inside the locking range (filled) and
the slip frequency outside it (open), in units of $K$, against the circle $(1-\nu^2)^{1/2}$ and the hyperbola
$(\nu^2-1)^{1/2}$.*

## 4. A symmetry sets the ratio

The van der Pol oscillator and the pumped oscillator have the same reflection,
$(x, y) \to (-x, -y)$, which maps the cycle onto itself half a period later, so that $Z(\theta+\pi) = -Z(\theta)$. A
drive that reverses under the reflection, such as the bias of the van der Pol oscillator, has
$Z_p(\theta+\pi) = -Z_p(\theta)$. Its phase response contains only odd harmonics: the constructor finds
$\lvert Z_2\rvert/\lvert Z_1\rvert < 10^{-6}$ and $\lvert Z_3\rvert/\lvert Z_1\rvert = 0.12$, and the oscillator
locks at 1:1. A drive that the reflection leaves unchanged, such as the stiffness $k$ of the pumped oscillator, has
$Z_p(\theta+\pi) = Z_p(\theta)$ and only even harmonics. A modulation at the free frequency then does not act at
first order, and a modulation at twice the frequency locks the oscillator with two phases half a cycle apart. The
drive holds either phase, and the history decides which one is occupied. This is the storage principle of the
parametron, in which a binary digit is the phase of an oscillator pumped at twice its frequency (von Neumann, 1957;
Goto, 1959). The ring of three repressors has the cyclic symmetry of Module 3, which maps its cycle onto itself a
third of a period later. The common promoter strength $\alpha$ respects this symmetry, so $Z_\alpha$ contains only
the harmonics 3, 6, …, and the ring locks at 3:1 with three phases. The symmetry that forbids a pitchfork in the ring
of three (Module 9, Section 5) thus fixes its locking ratio.

## 5. The junction

For the overdamped junction, $\dot\varphi = I - \sin\varphi$, the phase response along the current
has a single harmonic,
$Z(\theta) = (I^2 + \cos\theta + \omega_0\sin\theta)/(\omega_0 I)$ with $\omega_0 = (I^2-1)^{1/2}$, up to the origin
of $\theta$. The constructor finds $K/\varepsilon = 0.40032 = 1/(2\omega_0)$ at $I = 1.6$, and no higher harmonic
above $10^{-9}$. A microwave current of amplitude $\varepsilon$ therefore locks the junction over the frequencies
$\omega_0 \pm \varepsilon/(2\omega_0)$. Since $d\omega_0/dI = I/\omega_0$, the locking range is a step of half-width
$\varepsilon/(2I)$ in the bias current, the first Shapiro step at small amplitude (Shapiro, 1963). The other models
have phase responses of different shapes (Figure 1b), but averaging keeps only their $n$-th harmonic.

## 6. Where the derivation stops

In the Lotka–Volterra model the quantity $x - \ln x + y - g\ln y$ is conserved, so
every orbit around the centre is a cycle. A transverse Floquet multiplier equals 1, and the amplitude is as flat as
the phase. A periodic drive moves the orbit across the family without a restoring force; there is no isolated phase
to fix. In the Rosenzweig–MacArthur model of the same field, logistic prey growth and saturating predation make the
cycle attracting ($\kappa = 0.12$), and it locks. The capillary rotors and the in-plane dipoles stop because their
control multiplies the whole drift: modulating it changes only the rate of time, so $Z_p$ is constant and has no
harmonic. The models of the write targets stop at C, since none oscillates in its control range.

## 7. Two invariants

The phase gained over one period in the driven model agrees with $TK\cos(n\psi-\phi_n)$ to 3.5%
of $K$ at $K = 0.01\min(\omega_0, 2\kappa)$ (Figure 1c); the deviation is of first order in the drive. In units of
$K$ the measured relaxation rates and slip frequencies lie on the circle and the hyperbola of Eq. (1) (Figure 1d). A
joint fit of both gives the half-width and the centre of the locking range. Extrapolated linearly in the amplitude
from the two drives, the half-width is $1.0003 \pm 0.0006$ in units of $K$ ($\chi^2 = 0.86$ for 8 values), and no
model deviates from 1 by more than $3.2 \times 10^{-3}$. The largest uncertainty belongs to the ring of three, whose
third harmonic is small ($\lvert Z_3\rvert = 0.033$), so that its drive must be strong and the second-order
corrections are large. The centre moves in proportion to the amplitude, from $0.012K$ to $0.006K$ for the van der Pol
oscillator when $K$ is halved: a frequency pull of second order in the drive.

## 8. Memory in the locked phase

In Module 6 a written phase was retained only because nothing restored it, and its
variance grew linearly in time. Within the locking range the drive restores the phase at the rate $\lambda$, so for
weak noise the variance saturates, and the phase is lost only by slips over the barriers of the tilted periodic
potential of Eq. (1) (Pikovsky, Rosenblum and Kurths, 2001). For a ratio $n:1$ the drive holds $n$ phases. The three
memory targets of Modules 9 to 11 thus sort the models by their states: the write models reach only the writes, the oscillators reach
only phase locking, and the ring of three, which stops at R for the symmetric write and has no stable state for the
threshold write, locks at 3:1.

## 9. A model from your field

Add the specification of your oscillator to the command of Section 3,

```bash
python3 -B -m fieldbridge memory codiscover --target phase-locking my_oscillator.json \
  examples/memory/oscillators/van_der_pol.json --out-dir build/my_phase
```

to find the ratio at which a periodic modulation of its control fixes its phase, and the range of frequencies over
which it does. If the derivation stops at C, the control range of the specification contains no oscillation; if it
stops at R, the cycles form a neutral family.

## 10. Exercises

1. The ring of three repressors does not lock at 1:1 when $\alpha$ is modulated. Which modulation locks it at 1:1?
2. Show that a microwave current of amplitude $\varepsilon$ gives the overdamped junction a first Shapiro step of
   half-width $\varepsilon/(2I)$ in the bias current.
3. The Lotka–Volterra model stops at R. What does a weak periodic drive do to its cycles instead, and which
   change of the model makes the cycle lock?

<details><summary>Answers</summary>

1. A modulation that breaks the cyclic symmetry, for example of the promoter of one gene only. Give each gene its
   own strength $\alpha_i$ and modulate $\alpha_0$: the drive no longer commutes with the shift by one gene, $Z_p$
   acquires a first harmonic, and the ring locks at 1:1. With the common $\alpha$ it locks only at 3:1.
2. With $K = \varepsilon\lvert Z_1\rvert/2 = \varepsilon/(2\omega_0)$, the junction is locked where
   $\lvert\omega_0(I) - \omega_f\rvert < \varepsilon/(2\omega_0)$. Near the step $\omega_0(I)$ changes at the
   rate $d\omega_0/dI = I/\omega_0$, so the step extends over $\lvert I - I_f\rvert < \varepsilon/(2I)$, where
   $\omega_0(I_f) = \omega_f$. The time-averaged voltage, proportional to the frequency, is constant on the step.
3. The drive moves the state from one cycle of the family to another. Because the period depends on the cycle, this
   changes the frequency without a restoring force, and no isolated locked phase appears. A term that makes the cycle
   attracting removes the family: logistic prey growth and saturating predation, as in the Rosenzweig–MacArthur
   model.

</details>

## Summary

- Phase locking, the Adler equation for the phase difference to a periodic drive, is reached in eight oscillators
  from eight fields by four classes of derivation.
- A symmetry that maps the cycle onto itself $1/m$ of a period later sets the ratio $m:1$ (the parametron, the ring
  of three).
- The half-width of the locking range is $1.0003 \pm 0.0006$ in units of $K = \varepsilon\lvert Z_n\rvert/2$.
- A family of neutral cycles, a control that only rescales time, and a control range without oscillation stop the
  derivation, each at a named letter.
- Within the locking range the drive restores the phase, and the phase is lost only by slips.

## Reference

| Result | Function | Test |
| --- | --- | --- |
| Phase locking: cycle, phase response, ratio, Adler form | [`phase_locking.derive`](../../fieldbridge/memory/phase_locking.py), `find_cycle`, `floquet`, `phase_response`, `drive_symmetry`, `averaged_drift` | `test_the_symmetry_of_the_cycle_sets_the_locking_ratio`, `test_phase_locking_obstructions_name_their_reason` |
| Locking law | `phase_locking.adler_law`, `locked_rate`, `slip_frequency` | `test_the_junction_has_a_one_harmonic_phase_response_and_the_adler_width` |
| Report and figure | [`cli.cmd_codiscover`](../../fieldbridge/memory/cli.py), [`visual.codiscovery_figure`](../../fieldbridge/memory/visual.py) | `test_codiscover_command_writes_report_and_figure` |

Sources: R. Adler, Proc. IRE 34, 351 (1946); B. van der Pol, Phil. Mag. 2, 978 (1926); I. Prigogine and R. Lefever,
J. Chem. Phys. 48, 1695 (1968); R. FitzHugh, Biophys. J. 1, 445 (1961); J. Nagumo, S. Arimoto and S. Yoshizawa,
Proc. IRE 50, 2061 (1962); M. L. Rosenzweig and R. H. MacArthur, Am. Nat. 97, 209 (1963); M. L. Rosenzweig, Science
171, 385 (1971); B. C. Goodwin, Adv. Enzyme Regul. 3, 425 (1965); W. C. Stewart, Appl. Phys. Lett. 12, 277 (1968);
D. E. McCumber, J. Appl. Phys. 39, 3113 (1968); S. Shapiro, Phys. Rev. Lett. 11, 80 (1963); J. von Neumann, US
patent 2,815,488 (1957); E. Goto, Proc. IRE 47, 1304 (1959); M. B. Elowitz and S. Leibler, Nature 403, 335 (2000);
A. J. Lotka, *Elements of Physical Biology* (1925); V. Volterra, Nature 118, 558 (1926); A. Pikovsky, M. Rosenblum
and J. Kurths, *Synchronization* (Cambridge University Press, 2001).

[Previous: Module 10](25_memory_threshold_write.md) · [Next: Module 12](27_memory_return_point.md) · [Tutorial index](index.md) · [Glossary](memory_glossary.md)
