# Computation, Module 1. Which functions of an input history a body represents

**Learning objectives.** After this module you can

1. describe a body driven by an input, and the observables measured on it, in a specification of schema
   `fieldbridge-computation/1` (or `fieldbridge-springs/1` for a network of springs and masses);
2. state from the equations alone whether the body has fading memory, how many independent functions of the input
   history a linear combination of its observables can reproduce, and which degrees are excluded by symmetry;
3. run `computation predict`, `computation card` and `computation survey`, and read the class, the profile and the
   capacities by degree;
4. state a capacity together with the measurement noise and the amplitude at which it holds.

**Prerequisites.** Modules [1](15_memory_first_card.md) and [2](16_memory_specification.md) of the memory part, for
the realization and its specification. **Time.** About 40 minutes.

## 1. The question

A body that is driven by an input carries traces of the input's past in its state. A neuron's membrane potential
depends on the currents of the last milliseconds; the methylation of a bacterium's receptors on the attractant of the
last seconds; the resistance of a magnetic tunnel junction on the voltages of the last pulses. An experimenter who
measures a few observables and adds them with constant weights can reproduce some functions of the input history and
not others. The question of this module is which functions: how many, of which degree in the inputs, and from how far
back.

## 2. The measure

The input $u_t$ is drawn independently in every interval of length $\Delta t$ and held within it. For a target $z$,
a function of the past inputs, the **capacity** of the measured signals $y$ is the fraction of the variance of $z$
that the best linear combination reproduces (Dambre et al., 2012):

$$C[z] = 1 - \min_w \frac{\langle (z - w \cdot y)^2 \rangle}{\langle z^2 \rangle}.$$

The targets are products of normalized Legendre polynomials of delayed inputs, $\prod_j P_{d_j}(u_{t-k_j}/A)$, which
are orthonormal for inputs uniform on $[-A, A]$ (Hermite polynomials for Gaussian inputs; products of distinct
delayed inputs for binary ones). The **degree** of a target is $\sum_j d_j$. With fading memory and without noise the
capacities of all targets together equal the number of independent measured signals (Dambre et al., 2012); the
module shows how that total is shared between degrees and delays.

A specification has the `equations` body of the memory module (carrier, parameters and drift, read by the same
restricted parser) and a block `computation`:

| Key | Meaning |
| --- | --- |
| `input` | a declared parameter, redrawn for every interval and held within it |
| `law`, `amplitude`, `offset` | `uniform` on $[-A, A]$, `gaussian` with deviation $A$, or `binary` $\pm A$; a constant added |
| `hold` | the length of an interval |
| `observables` | expressions in the variables and parameters |
| `virtual_nodes` | measurement times per interval (default 1) |
| `washout`, `length`, `substeps` | intervals discarded and used; Runge–Kutta steps per interval |
| `solver`, `form` | `rk4` or `stiff`; `flow` (a rate) or `map` (the next state) |
| `initial`, `exact_box`, `expect` | a starting state, the range of the state for the exact capacities, the expected class |

In the language of the realization: $\Xi$ is the set of states, $\Omega$ the drift, $C$ the closure, $R$ the
observables, $P$ the input law and its hold, and $A$ the parameters.

## 3. The prediction from structure

`computation predict` linearizes the body at its steady state and reads four things from the equations, without a
simulation:

1. **Fading memory.** Over one interval the deviation obeys $s_t = \Phi s_{t-1} + \Gamma u_t$. The class is decided
   by the largest $|\lambda(\Phi)|$ over the modes that the input reaches and the observables see: below 1 the
   memory fades; equal to 1 a mode keeps the running sum of the input (class *no fading memory*). Conserved amounts
   that the input cannot change are not among these modes.
2. **Rank.** The impulse responses of the measured signals, $g_k$, span $n_\mathrm{lin}$ directions: the linear
   capacity at small amplitude, at most the number of measured signals and the number of modes reached and seen
   (Gonon, Grigoryeva and Ortega, 2020).
3. **Profile.** The capacity for the input $k$ intervals back is $C(k) = g_k^\top (\sum_j g_j g_j^\top)^{+} g_k$. For
   one mode, $C(k) = (1 - a^2)a^{2k}$ with $a = e^{-\kappa\Delta t}$.
4. **Symmetry.** If the drift and the observables are odd about the steady state and the input law is symmetric,
   every even degree has zero capacity. If they are affine, only degree 1 has capacity.

These decide the class of the card: *linear memory*, *odd degrees only*, *nonlinear capacity* or *no fading memory*.

```bash
python -m fieldbridge computation predict examples/computation/benchmarks/one_mode_map.json
```

For one linear mode with $a = 0.8$ the profile is $0.36 \cdot 0.64^k$ and sums to 1. For three linear stages in a
chain, measured at every stage, $n_\mathrm{lin} = 3$; measured at the last stage only, $n_\mathrm{lin} = 1$; measured
at the last stage three times per interval, $n_\mathrm{lin} = 3$ again. The oscillator with a cubic restoring force,
$\ddot x = -g\dot x - \omega^2 x - \beta x^3 + u$, is odd, and only odd degrees carry capacity. A linear stage that
drives the square of its state into a second stage ends at degree 2: the first stage holds the linear memory and the
second only products of two inputs, so the capacities are 1 at degree 1 and 1 at degree 2, summed over all delays,
at every amplitude. Computed on a grid of the state they approach these values as the grid is refined (degree 2:
0.976 on 41 × 41 points, 0.994 on 101 × 101, delays up to 8).

## 4. Running the command

`computation card` adds capacities to the predictions. For a body with one or two state variables they are computed
exactly from the equations: driven by independent inputs, the state after each interval is a Markov chain, and the
covariance of a signal with a target follows by a backward recursion of conditional expectations over the delays, on
a grid of the state, with the stationary distribution of the discretized transfer operator. There is no expansion in
the amplitude and no sampling noise. With `--simulate` the card also estimates the capacities from a simulation of
the full body (the estimator of Dambre et al., 2012, with their threshold: twice the $\chi^2$ quantile).

```bash
python -m fieldbridge computation card examples/computation/chemotaxis_methylation_tu2008.json --out-dir out/chemotaxis
```

The chemotaxis model of Tu, Shimizu and Berg (2008) is driven by $u = \ln([L]/[L]_0)$ around 158 μM with amplitude
1.151, held for $0.05/F(0)$, and measured through the kinase activity and the methylation level. Its memory fades at
the rate 24 $F(0)$ of the steep band of the methylation rate; it is not odd. From the equations: degree 1 holds 1.731,
degree 2 0.162 and degree 3 0.031, against 1.743, 0.161 and 0.029 estimated from a simulation of 95 000 intervals
(`--simulate --length 6000`).

For larger bodies the capacities come from the simulation alone:

```bash
python -m fieldbridge computation card examples/computation/hodgkin_huxley1952.json --simulate --amplitudes 1,4 --out-dir out/hh
```

The membrane of Hodgkin and Huxley (1952), driven by a current held for 1 ms and measured through its potential,
fades at 0.12 per ms over its four modes. Its predicted profile returns at delays of 5–8 ms (0.024, 0.034, 0.032):
the damped subthreshold oscillation of the membrane. Below threshold (1 μA/cm²) all of the capacity is linear: 1.009
at degree 1 and nothing above. At 4 μA/cm² the membrane fires, and only 0.316 of the capacity lies within degree 3
and delay 8 (0.211, 0.093, 0.012); the spikes spread the rest over higher degrees and longer delays.

## 5. Controls

Each control changes one component and changes the answer:

| Control | Change | Class |
| --- | --- | --- |
| `controls/odd_oscillator_bias.json` | $P$: offset of the input 0 → 0.5 | odd degrees only → nonlinear capacity |
| `controls/linear_chain_no_decay.json` | $A$: decay of the first stage 0.5 → 0 | linear memory → no fading memory |
| `benchmarks/running_sum.json` | a mode without decay | no fading memory |
| the linear chain measured at its last stage | $R$: three observables → one | $n_\mathrm{lin}$ 3 → 1 |

With the bias, degree 2 takes about 0.1 of the capacity of the cubic oscillator, against exactly 0 without it.

## 6. Measurement noise

Without noise, the capacity counts every independent direction of the measured signals, however weak. Measurement
at many times within an interval, impulse responses that fall quickly, and nonlinear distortion of weak directions
produce directions whose strength is $10^{-4}$ to $10^{-8}$ of the largest. A measured signal resolves them only
above its noise. The card therefore states the capacity at a relative measurement noise $\varepsilon$: the covariance
of the signals gains $\varepsilon^2$ times its diagonal.

```bash
python -m fieldbridge computation card examples/computation/spin_torque_furuta2018.json --noise 1e-3 --delays 1:30,2:6,3:4
```

The macrospin junction of Furuta et al. (2018) is driven by binary pulses of ±44 mV held for 20 ns; its resistance is
measured 50 times per pulse. With the reference layer perpendicular to the plane, the resistance depends on $m_z$
alone, the azimuth is hidden from it, and the body reduces to one variable. Its linear response has two directions
(singular values 1 and 0.4: the state at the start of a pulse and the pulse itself). The 50 measured signals are
nonlinear functions of the same two quantities and span many more directions, whose strengths fall geometrically:
the eigenvalues of their covariance are 1, 0.23, $2.9 \times 10^{-3}$, $2.4 \times 10^{-4}$, $1.2 \times 10^{-5}$, $8.8 \times 10^{-7}$, … of the largest. From the equations:

| Relative noise | Degree 1 | Degree 2 | Degree 3 | Degree 1, delays 1–30 |
| --- | --- | --- | --- | --- |
| $10^{-2}$ | 2.341 | 1.638 | 0.802 | 1.344 |
| $10^{-3}$ | 2.456 | 2.254 | 1.181 | 1.458 |
| $10^{-4}$ | 2.855 | 2.687 | 1.508 | 1.857 |

The source reports a short-term memory (delays 1–30) of about 2.3 from simulations without noise; in the table, that
value is approached only for a resistance measured to better than $10^{-5}$.

## 7. The published bodies

```bash
python -m fieldbridge computation survey examples/computation
```

The survey gives the predictions of every specification in a folder, and the exact capacities of those with one or
two state variables (without measurement noise, delays up to 40, 8 and 5 for degrees 1, 2 and 3, on a grid of 2001
points for one variable and 161 × 161 for two):

| Body | Field | Class | Memory, slowest rate | Modes reached and seen | n_lin / signals | Odd | Exact capacity by degree |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Three linear stages in a chain | — | linear memory only | fading, 0.5 | 3 of 3 | 3 / 3 | yes | — |
| Damped oscillator with a cubic restoring force (Duffing) | — | odd degrees only | fading, 0.25 | 2 of 2 | 2 / 2 | yes | 1: 1.9, 2: 0, 3: 0.0435 |
| One linear mode, x_t = a x_{t-1} + u_t | — | linear memory only | fading, 0.223 | 1 of 1 | 1 / 1 | yes | 1: 1, 2: 0, 3: 0 |
| A running sum of the input | — | no fading memory | none (integrating) | 1 of 1 | 1 / 1 | yes | — |
| A linear stage that drives the square of its state into a second stage | — | nonlinear capacity | fading, 0.5 | 1 of 2 | 1 / 2 | no | 1: 1, 2: 0.997, 3: 0 |
| Three linear stages, the first without decay | — | no fading memory | none (integrating) | 3 of 3 | 3 / 3 | yes | — |
| Duffing oscillator with a constant bias of the input | — | nonlinear capacity | fading, 0.25 | 2 of 2 | 2 / 2 | no | 1: 1.83, 2: 0.108, 3: 0.0169 |
| E. coli chemotaxis: receptor methylation driven by a held log-ligand input | bacterial chemotaxis | nonlinear capacity | fading, 24 | 1 of 1 | 2 / 2 | no | 1: 1.73, 2: 0.162, 3: 0.0307 |
| Echo state network of 50 variables (Dambre et al. 2012 parameterization) | recurrent networks | odd degrees only | fading, 0.622 | 50 of 50 | 18 / 50 | yes | — |
| Hodgkin-Huxley space-clamped membrane driven by a held current | neurophysiology | nonlinear capacity | fading, 0.121 | 4 of 4 | 1 / 1 | no | — |
| MAPK cascade (Huang and Ferrell 1996) driven by a held log-input of its activator E1 | cell signalling | nonlinear capacity | fading, 0.0317 | 4 of 18 | 1 / 1 | no | — |
| Random network of 30 masses and 78 nonlinear springs (Hauser et al. 2011 construction) | soft robotics and morphological computation | nonlinear capacity | fading, 0.0492 | 112 of 112 | 30 / 78 | no | — |
| Macrospin magnetic tunnel junction with a perpendicular reference layer, reduced to m_z | spintronics | nonlinear capacity | fading, 0.158 | 1 of 1 | 2 / 50 | no | 1: 3.07, 2: 3.05, 3: 2.16 |

The junction values hold without measurement noise and therefore include directions that no measurement resolves
(Section 6).

The mass–spring network is built by the procedure of Hauser et al. (2011): 30 masses placed at random, 78 springs
along the edges of their Delaunay triangulation, cubic stiffness and damping, a horizontal force on six nodes. It is
given as data in a specification of schema `fieldbridge-springs/1`, since its equations, written out, would exceed what
the expression parser accepts. A spring length is not an odd function of the node positions, so even degrees appear.

## 8. A model from your field

Write the equations of your body as an `equations` specification, add a block `computation` with the input you can
apply and the observables you can measure, and run `computation predict` first: it takes seconds and tells you the
class, the rank and the profile. Then run `computation card`, with `--noise` set to the precision of your
measurement, and with `--simulate` for a body of more than two variables. A specification in
`examples/computation` appears on the web page in the column of its class.

## 9. Exercises

1. In `benchmarks/one_mode_map.json`, set $a = 0.5$. Predict the profile and its sum before running the command.
2. Measure the linear chain at its last stage two times per interval. What is $n_\mathrm{lin}$?
3. Give the cubic oscillator a bias of 0.2 instead of 0.5. Is degree 2 smaller or larger than with 0.5?
4. Run the junction card at $\varepsilon = 10^{-1}$. Which degree loses the most?
5. Add a linear observable $x_1 + x_2$ to the square cascade. Does the class change? Does the rank?

## Summary

- The capacity of measured signals for a function of the input history is the fraction of its variance that a
  linear combination reproduces (Dambre et al., 2012).
- From the equations: fading memory over the modes reached and seen, the rank of the linear response, the degree-1
  profile, and the degrees excluded by symmetry.
- For one or two state variables, the capacities of every degree and delay follow exactly from the equations through
  the driven Markov chain; for larger bodies, from a simulation.
- A capacity is stated at an amplitude and a measurement noise: without noise it counts directions no measurement
  resolves.

## Reference

`fieldbridge.computation`: `spec` (the two schemas), `predict` (the linearization and the symmetry checks), `exact`
(the driven Markov chain), `ipc` (the estimator), `springs` (networks of springs and masses), `card`, `cli`. Examples
in `examples/computation`: six published bodies, five benchmarks, two controls. Sources: Dambre et al. (2012), Jaeger
(2001, 2002), Boyd and Chua (1985), Gonon, Grigoryeva and Ortega (2020), Kubota, Takahashi and Nakajima (2021),
Hodgkin and Huxley (1952), Tu, Shimizu and Berg (2008), Huang and Ferrell (1996), Hauser et al. (2011), Furuta et al.
(2018); see docs/mechanisms.md.
