# From a discrepancy to a physical question

The writing calculation in [the previous chapter](38_stochastic_self_correction.md)
finds a rate that produces the desired retained pattern. Reaching that target
leaves a physical question: which part of the intervening motion determines
the pattern? The final equilibria alone cannot answer it, since every writing
protocol ends with the same equations.

## Keep the signed difference

For each final state, the probability difference is

$$
e_i=p_i^{\mathrm{reference}}-p_i^{\mathrm{prediction}}.
$$

The `memory investigate` command retains the sign, sampling precision and
integration-step difference for every state, including states from accepted
predictions. Independent repetitions are grouped by the same equations,
preparation, protocol, numerical policy and state catalogue. Copied records and
reused seed pairs do not count as new repetitions.

```bash
python3 -B -m fieldbridge memory investigate \
  build/stochastic_passage/*/comparison.json \
  --out-dir build/probability_inquiry
```

The JSON contains every state and its proposed follow-up. A numerical
resolution check takes precedence when the integration-step difference is
comparable to the measured departure. Repeated departures are assessed using
their sampling uncertainty and a multiple-comparison correction. A small
effect can therefore become a question even when the original approximation
was accepted. Opposite signs in repeated calculations do not establish a
systematic correction.

This examination distinguishes a discrepancy that needs explanation from the
question of originality. It can identify an ordinary finite-rate correction,
a missed interaction or a numerical problem; the mathematical and physical
explanation comes from the next calculation.

## Change the interaction during writing

Consider the three-amplitude source model,

$$
\mathrm d q_i=\left[(\mu-g_i)q_i-q_i^3+
\kappa\sum_j(q_j-q_i)+h_i\right]\mathrm dt+
\sqrt{2D}\,\mathrm dW_i.
$$

Increase the coupling by the same amount for the same duration, but place that
pulse at different stages of the sweep. Every pulse ends before the final
control value. The bias is then removed and every trajectory undergoes the
same relaxation and noisy hold. The final drift and its stable equilibria are
unchanged.

The pulse has an exact property derived from the equations. Define

$$
m=\frac{q_1+q_2+q_3}{3},\qquad
V=\frac12\sum_i(q_i-m)^2.
$$

Here $V$ measures differences among the three amplitudes within one
trajectory; it is not the variance across an ensemble of experiments. The
coupling contribution sums to zero, so the pulse makes no direct contribution
to $\dot m$. Its contribution to the evolution of $V$ is

$$
\left.\dot V\right|_{\mathrm{pulse}}=-6\,\delta\kappa V.
$$

It therefore suppresses differences among the amplitudes while preserving
their instantaneous mean. The full mean dynamics still depend on those
differences. With $\delta_i=q_i-m$, they satisfy

$$
\mathrm d m=\left[\mu m-\frac13\sum_i g_iq_i-m^3-2mV
-\frac13\sum_i\delta_i^3+\frac13\sum_i h_i\right]\mathrm dt
+\sqrt{\frac{2D}{3}}\,\mathrm dW_m.
$$

Changing the differences can thus change subsequent mean evolution and final
selection, despite the zero direct mean force. The source calculation also
evaluates the Lie bracket between mixing and the original drift. Its action on
the mean is nonzero: the two evolutions cannot generally be exchanged in time.
This is a concrete explanation to test, rather than a conclusion inferred from
the visual arrangement of an atlas.

## Select a timing contrast, then test it afresh

```bash
python3 -B scripts/investigate_writing_path.py \
  --task examples/learning/writing_path_source.json \
  --out-dir build/writing_path_inquiry \
  --pilot-n 1024 --n 4096 --workers 2
```

The exploratory calculation compares five pulse timings and a no-pulse
reference. The nonzero pulses have the same integrated excess coupling at a
given rate. The timings with the largest and smallest predicted all-positive
occupation are selected. Full probability predictions for these timings are
then saved before fresh reference trajectories are evaluated.

The reference calculations use a different integration scheme, paired
Brownian step refinement and common noise across pulse schedules. Common noise
allows the timing contrast to be measured as a paired difference between
trajectory outcomes. A second test changes both the rate and the initial
preparation, without changing the selected pulse windows or final equations.

The report includes the final state catalogue, probability contrasts, the
evolution of the mean and internal dispersion, and equilibria found during the
pulse. Temporary loss of a stable branch and redistribution among basins are
different explanations; the changing intermediate equations must be inspected
before deciding between them. The final landscape is fixed, not the landscape
at every point during the pulse.

The question and pulse family are supplied. Selection of the timing contrast,
probability prediction and fresh trajectory checks are automated. The inquiry
continues after a successful construction, without requiring a new fundamental
law as its outcome.

## Use a small threshold split to construct selective erasure

At the selected stage, $\mu=0.16$, the source equations have two mixed-sign
states. They disappear at different coupling strengths,

$$
\kappa_-=0.02795253,\qquad \kappa_+=0.02929186.
$$

The constructor finds these folds by solving the equilibrium equations
together with a zero Jacobian determinant. Each solution has one zero
eigenvalue, two negative eigenvalues, and nonzero parameter and quadratic
coefficients of its local normal form. These checks distinguish a fold from
an unresolved root or an accidentally small eigenvalue.

The interval between the folds has a direct use. Raising the coupling to its
midpoint, $\kappa=0.02862219$, removes one mixed-sign state while the other
survives. The pulse duration is calculated from the local passage time near
the lower fold. A pause at this stage allows the targeted state to leave its
former basin. The sweep then resumes, and the original coupling is restored.
All protocols end with the same eight checked stable states.

```bash
python3 -B scripts/construct_selective_erasure.py \
  --task examples/learning/writing_path_source.json \
  --parent build/writing_path_inquiry \
  --out-dir build/selective_erasure --n 256
```

The constructed operating regime is deliberately low-noise: the local barrier
estimate sets $D=1.07\times10^{-7}$, rather than the earlier writing diffusion
$1.25\times10^{-4}$. The chosen long pulse lasts about 483 time units. These
settings are computed and saved before fresh reference trajectories are used.
At this noise level, all 256 reference trajectories prepared in the targeted
$(+,-)$ pattern end in $(+,+)$, while all 256 prepared in the protected
$(-,+)$ pattern retain it. These signs refer to the first two amplitudes;
the complete three-coordinate final states are checked against the shared
equilibrium catalogue.

A weaker pulse preserves both patterns. A pulse above both fold thresholds
erases both mixed-sign patterns. A short pulse at the selective amplitude
preserves the target because it does not allow enough time to pass through
the bottleneck. Two deterministic solvers also reproduce the state mapping.
The stored recipe records the fold interval, timing rule, tested noise level
and controls; applying it to another model requires calculating that model's
folds and testing the proposed pulse.

## Follow an effect outside the original target

The late pulse leaves the original all-positive occupation unchanged in the
tested sample, but five of 4,096 paired trajectories exchange two other final
patterns. This hint is retained even though the probability predictions pass.

```bash
python3 -B scripts/confirm_writing_window.py \
  --parent build/writing_path_inquiry --out-dir build/late_pulse --n 4096
python3 -B scripts/follow_small_occupation_change.py \
  --parent build/late_pulse --out-dir build/rare_pattern_followup \
  --maximum-n 16384
```

The follow-up chooses its larger sample budget from the measured precision,
then freezes it before drawing new trajectories. Exact paired tests assess
rare changes, and zero observed changes retain a finite uncertainty interval.
The recorded fresh calculation finds five changes in one direction and one
in the other among 16,384 trajectories. That does not resolve a repeatable
redistribution after correction across the tested states. The hint remains
an open question, alongside the resolved timing and selective-erasure results.

Implementation: [signed departures](../../fieldbridge/memory/investigation.py),
[source-derived mixing and trajectories](../../fieldbridge/memory/pulse_inquiry.py),
[timing study](../../scripts/investigate_writing_path.py),
[selective erasure](../../scripts/construct_selective_erasure.py),
[rare-effect follow-up](../../scripts/follow_small_occupation_change.py),
[tests](../../tests/test_physical_investigation.py).
