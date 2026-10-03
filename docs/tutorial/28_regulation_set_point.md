# Regulation, Module 1. Return to a set point

**Learning objectives.** After this module you can

1. describe a cell, a machine or a material whose output is held at a set point, in a specification of schema
   `fieldbridge-regulation/1`;
2. state from the equations alone when the output returns exactly after a step of the input: an integrator and a
   stable steady state;
3. run `regulation card` and read the class, the integrator, the clamps and the steps;
4. tell a leaky integrator from an exact one, and check a certificate without repeating the search.

**Prerequisites.** Modules [1](15_memory_first_card.md) and [2](16_memory_specification.md) of the memory part, for
the realization and its specification. **Time.** About 30 minutes.

## 1. The question

Many systems hold an output while an input changes. Bacteria keep the activity of their chemoreceptors after the
concentration of an attractant has changed (Barkai and Leibler, 1997). A car keeps its speed on a hill. A cell keeps
the level of a protein when its production is disturbed. After a step of the input the output moves, and then it
returns exactly, part of the way, or not at all. The question of this module is which of the three occurs, through
which variable the output returns, and how fast.

## 2. The model

A specification has the `equations` body of the memory module (carrier, parameters and drift, read by the same
restricted parser) and a block `regulation`:

| Key | Meaning |
| --- | --- |
| `input` | a declared parameter that the protocol steps |
| `output` | an expression in the variables and parameters: the measured observable $y = h(q, u)$ |
| `steps` | input values; the first is the reference value $u_0$, every other one a value after a step |
| `initial` | a state from which the steady state at $u_0$ is reached |
| `fixed` | parameters that are not varied when parameter points are sampled |
| `set_point`, `expect` | the set point and the class that the source states, compared with the result and never used |

In the language of the realization: $\Xi$ is the set of states, $\Omega$ the drift, $C$ the closure (for example a
well-mixed cell), $R$ the output, $P$ the step of the input and $A$ the parameters.

## 3. The prediction from structure

At a steady state $F(q^*, u) = 0$ the output changes with a constant input as

$$
G = \frac{dy}{du} = h_u - \nabla h \cdot J^{-1} F_u = \frac{\det \begin{pmatrix} J & F_u \\ \nabla h & h_u \end{pmatrix}}{\det J} . \qquad (1)
$$

A variable can be clamped, that is, held at its steady value as a reservoir would hold it. The gain with that
variable clamped, $G_{\rm open}$, is the response without the part of the loop that passes through the variable, and
$G/G_{\rm open}$ is the fraction of the step that remains. The output adapts perfectly when $G = 0$ at every parameter
point, and it is fine-tuned when $G = 0$ only at the stated parameters.

A robust return needs an integrator (the internal model principle: Francis and Wonham, 1976; Sontag, 2003; for
chemotaxis, Yi et al., 2000). This is a function $\phi$ of the state with

$$
\frac{d\phi}{dt} = g(q)\,\bigl(y - y_0\bigr), \qquad g \text{ of one sign}. \qquad (2)
$$

At a steady state $d\phi/dt = 0$, so $y = y_0$ whatever the input and whatever the other parameters: an integrator and
a stable steady state make the return exact. The command searches for $\phi$ by linear algebra on the drift at sampled
states and inputs, in three stages:

1. a constant gain, with $\phi = w \cdot q + v \cdot \ln q$ (logarithms for concentrations);
2. a gain linear in the state;
3. a rate $w \cdot F$ that depends on the state only through the output.

A combination with $w \cdot F = 0$ is a conservation law, not an integrator. An integrator in other coordinates is not
found; a return without a detected integrator is reported as such.

A leaky integrator, $d\phi/dt = k(y - y_0) - \delta \phi$, leaves the fraction $1/(1 + G_{\rm loop})$ of the step,
with loop gain $G_{\rm loop} = k g/\delta$. This is standard loop analysis (Åström and Murray, 2021); the command uses
it as a calibration, not as a law.

## 4. Running the command

```bash
python3 -B -m fieldbridge regulation card examples/regulation/chemotaxis_tu2008.json --out-dir build/chemotaxis
python3 -B -m fieldbridge regulation check build/chemotaxis/regulation.json
```

The card finds the steady state, the gains and the clamps, and the integrator. It repeats the gain at 32 parameter
points, each free parameter multiplied by a factor between 1/2 and 2. It then follows the output after each step of the
input. The integral of the rate of $\phi$ over the response equals the change of $\phi$; the report gives their ratio
as a calibration of the simulation. The report `regulation.md` states the class, the integrator, the clamps and the
steps, and `regulation.json` holds a certificate (Section 7).

The published models in [`examples/regulation`](../../examples/regulation) come from five fields (32 parameter points,
default seed):

| Model | Field | Class | Integrator | G/G_open | Final / peak after the first step |
| --- | --- | --- | --- | --- | --- |
| antithetic integral feedback on gene expression | synthetic biology | perfect-adaptation | -z1 +z2 (constant gain 1) | 0 | 0 |
| E. coli chemotaxis: receptor methylation (MWC clusters) | bacterial chemotaxis | perfect-adaptation | -m (rate a function of the output) | 0 | 0 |
| cruise control of a car on a hill (PI controller) | engineering | perfect-adaptation | -z (constant gain 1) | 0 | 0 |
| EnvZ-OmpR two-component signalling | chemical reaction networks | perfect-adaptation | XD -XT -6 Xp -6 XpY -6 Yp -4 XDYp (gain 0.5 XD) | 0 | 0 |
| incoherent feedforward loop with a proportioner node (IFFLP) | biochemical networks | partial-adaptation | none found | 0.00295 | -0.00314 |
| negative feedback loop with a buffer node (NFBLB) | biochemical networks | partial-adaptation | none found | 0.308 | 0.812 |
| negative feedback loop with a buffer node, enzymes on B closer to saturation (NFBLB) | biochemical networks | partial-adaptation | none found | 0.0497 | 0.161 |
| quasi-integral control in a growing cell (type II) | synthetic biology | partial-adaptation | none found | 0.0216 | 0.167 |

Every model with an integrator returned exactly at every sampled parameter point with a stable steady state. In the
antithetic controller one of the 32 points lies where the deterministic model oscillates, as Briat et al. (2016)
describe. The methylation level of chemotaxis is found by the third stage: its rate $F(a)$ depends on the receptor
activity alone and decreases through $a_0 = 1/3$. In EnvZ–OmpR the integrator has a gain proportional to the
concentration of EnvZ-ADP. This is the linear constrained integrator of networks with absolute concentration
robustness (Cappelletti, Gupta and Khammash, 2020).

## 5. Controls

Each published model has a control that must change the class:

| Control | Change | Class |
| --- | --- | --- |
| Hill-type static controller on gene expression (control, the source's comparison) | the antithetic pair replaced by the static controller f(x2) = alpha K^n/(K^n + x2^n) of Eq. (12), n = 1, alpha = 8.22, K = 3 (Fig. 3) | no-adaptation |
| chemotaxis with a methylation rate that also depends on the methylation level (control) | a turnover -kd (m - m0) added to F(a): the rate no longer depends on the activity alone | partial-adaptation |
| cruise control with a proportional controller only (control) | the integral term removed: u = kp (vr - v) | no-adaptation |
| EnvZ-OmpR with a second phosphatase of OmpR-P (control) | an added reaction Yp -> Y (k12 = 0.05) that does not pass through XD | partial-adaptation |
| leaky integral control in a growing cell (type II, eps = 1; control) | controller reactions as slow as in the leaky controller of Eq. (2.4): eps = 1 | partial-adaptation |
| ideal antithetic integral control without dilution (type II, gamma = 0; control) | dilution of the controller species removed (Eq. 2.2) | perfect-adaptation |

A clamp names the variable through which the output returns, and the remaining fraction measures how much of the step
that variable removes. In the static Hill controller the feedback removes a third of the step, but the output moves
to its new value without returning: feedback that is proportional reduces the deviation and leaves it.

## 6. Leaks

Dilution by growth makes the memory of a cell leak (Qian and Del Vecchio, 2018). Their antithetic controller with
reactions $1/\varepsilon$ times faster than dilution leaves the error $\varepsilon \gamma (z_1 - z_2)/k$ at steady
state (their Eq. 2.8). With $\varepsilon = 0.02$ the error is a few per cent of the step, with $\varepsilon = 1$ more
than half, and without dilution the return is exact. In the three-node networks of Ma et al. (2009) the Michaelis
constants of the enzymes acting on the buffer node play the role of the leak: the closer these enzymes are to
saturation, the smaller the remaining fraction.

A leak does not always leave a residual. When the controller is autocatalytic, $dz/dt = k z (y - y_0) - \delta z$, the
logarithm $\ln z$ integrates $y - y_0 - \delta/k$: dilution moves the set point to $y_0 + \delta/k$ and the return is
still exact (benchmark `autocatalytic_diluted`; compare Ruoff et al., 2019).

## 7. Certificates

The search for an integrator explores sampled states; checking a recorded one does not. `regulation check` takes a
report and its specification and verifies, in seconds and without a search:

1. that the specification is the one the report read;
2. that the recorded steady state solves $F = 0$ and is stable;
3. that the recorded integrator satisfies its identity at fresh states and inputs, and symbolically where the drift
   and the output are rational functions;
4. that the gain ratio agrees with the class.

For the antithetic pair, EnvZ–OmpR, the network with absolute concentration robustness and the multiplicative
feedforward loop, the identity simplifies to zero with rational coefficients: the integrator is exact at every state,
not only at the sampled ones. A check that passes shows that the report follows from the specification. It does not
show that the specification describes a material.

## 8. A model from your field

Write first the question that the model answers, then the input that the experiment steps and the observable that it
measures. Add the source with the location of every equation and value in it, and list under `assumptions` every value
that the source does not state. A control is a copy with one change that must change the class: a leak, a removed
variable, a parameter off a cancellation. `regulation survey` compares each specification with the class it expects.

## 9. Exercises

1. Copy the PI loop benchmark and add a leak $-\delta z$ with $\delta = 0.05$ to the integral. Predict the remaining
   fraction from $k_I = 0.5$, $\lambda = 1$ and $\delta$, and compare it with the clamp of $z$ in the card.
2. Why does the search report the antithetic integrator as $z_2 - z_1$ and not as either species?
3. Run the search on the network with absolute concentration robustness with the logarithmic columns switched off
   (`find_integrator(..., logs=False, stage2=False)`). What is found, and why?
4. In the subtractive feedforward loop, set $k_1 = 1.2$. What is the class now, and what has happened to the
   integrator?

<details>
<summary>Answers</summary>

1. The closed-loop gain is $1/(\lambda + k_I/\delta) = 1/11$ and the gain with $z$ clamped is $1/\lambda$, so the
   remaining fraction is $\delta\lambda/(\delta\lambda + k_I) = 1/11 \approx 0.091$, the value the clamp of $z$ reports
   (`test_clamping_the_leaky_integrator_gives_the_remaining_fraction`).
2. Each species has the sequestration term $\eta z_1 z_2$ in its rate, which is not a function of the error. The term
   cancels in the difference, whose rate is $\theta x_2 - \mu = \theta (x_2 - \mu/\theta)$
   (`test_antithetic_pair_is_found_as_a_difference`).
3. Only the conservation law $A + B + C$. No integrator linear in the concentrations exists for this network
   (Cappelletti, Gupta and Khammash, 2020); $\ln B$ integrates $A - k_2/k_1$ with gain $k_1$
   (`test_concentration_robustness_carries_a_logarithmic_integrator_and_a_conservation_law`).
4. Partial or no adaptation. The combination $(k_1/k_3)\,x - y$ integrates $y - y_0$ only when $k_1 k_4 = k_2 k_3$; off
   this surface no integrator exists (`test_subtractive_feedforward_has_an_integrator_only_where_it_cancels`).

</details>

## Summary

- A regulated realization is a body with an input that is stepped and an output that is measured.
- An integrator, a function of the state whose rate is the error times a gain of one sign, together with a stable
  steady state, makes the return exact at every parameter point where that steady state exists.
- The command finds integrators in linear and logarithmic coordinates, with a gain that depends on the state, and
  with a rate that is a function of the output; it separates them from conservation laws.
- A leak leaves a fraction of the step, unless it only moves the set point, as for an autocatalytic controller.
- A certificate lets a reader check a report against its specification without repeating the search.

## Reference

| Result | Function | Test |
| --- | --- | --- |
| Specification | [`spec.load`](../../fieldbridge/regulation/spec.py) | `test_refusals` |
| Static gain and clamps | [`gains.static_gain`, `clamped_gain`, `reference_gain`](../../fieldbridge/regulation/gains.py) | `test_static_gain_equals_bordered_determinant_and_finite_difference` |
| Integrator | [`integrator.find_integrator`](../../fieldbridge/regulation/integrator.py) | `test_logarithmic_integrator_needs_the_logarithmic_columns`, `test_methylation_level_is_an_integrator_of_a_function_of_the_activity` |
| Card | [`card.card`](../../fieldbridge/regulation/card.py) | `test_benchmark_classes`, `test_published_models_and_controls_get_their_class` |
| Steps | [`step.step_response`](../../fieldbridge/regulation/step.py) | `test_response_integral_equals_the_change_of_the_integrator` |
| Certificate | [`certificate.check`](../../fieldbridge/regulation/certificate.py) | `test_a_certificate_fails_against_a_changed_specification` |
| Report | [`cli.cmd_card`](../../fieldbridge/regulation/cli.py) | `test_command_line_writes_a_card_whose_certificate_passes` |

The tables of Sections 4 and 5 are those of
`python3 -B -m fieldbridge regulation survey examples/regulation examples/regulation/controls --out-dir build/regulation`.

Sources: B. A. Francis and W. M. Wonham, Automatica 12, 457 (1976); N. Barkai and S. Leibler, Nature 387, 913
(1997); T.-M. Yi, Y. Huang, M. I. Simon and J. Doyle, PNAS 97, 4649 (2000); E. D. Sontag, Syst. Control Lett. 50, 119
(2003); Y. Tu, T. S. Shimizu and H. C. Berg, PNAS 105, 14855 (2008); W. Ma, A. Trusina, H. El-Samad, W. A. Lim and
C. Tang, Cell 138, 760 (2009); G. Shinar and M. Feinberg, Science 327, 1389 (2010); D. F. Anderson, G. A. Enciso and
M. D. Johnston, J. R. Soc. Interface 11, 20130943 (2014); C. Briat, A. Gupta and M. Khammash, Cell Syst. 2, 15
(2016); R. P. Araujo and L. A. Liotta, Nat. Commun. 9, 1757 (2018); Y. Qian and D. Del Vecchio, J. R. Soc. Interface
15, 20170902 (2018); P. Ruoff, O. Agafonov, D. M. Tveit, K. Thorsen and T. Drengstig, PLoS ONE 14, e0207831 (2019);
D. Cappelletti, A. Gupta and M. Khammash, J. R. Soc. Interface 17, 20200437 (2020); K. J. Åström and R. M. Murray,
*Feedback Systems*, 2nd ed. (Princeton, 2021).
