# Correcting a stochastic writing prediction

A weak force acts on the first of three coupled bistable amplitudes while a
control parameter increases. The force is then removed. The final system has
several stable patterns, so predicting the sign of the first amplitude is not
enough: the calculation must predict the probability of each retained pattern.

The first model is

$$
\mathrm d q_i=\left[(\mu-g_i)q_i-q_i^3+
\kappa\sum_{j=1}^{3}(q_j-q_i)+h_i\right]\mathrm dt
+\sqrt{2D}\,\mathrm dW_i.
$$

The independent Brownian motions have the Itô convention. Here
$g=(0,0.02,0.25)$ and $\kappa=0.01$. After writing, every calculation ends at
$\mu=0.5$, removes $h$, relaxes without noise for 80 time units and then holds
for another 80 units at $D_{\mathrm{hold}}=10^{-5}$. The final drift is therefore
the same at every writing rate. Its stable equilibria are calculated from the
equations and checked through their Jacobian spectra. A trajectory far from
every checked equilibrium remains unassigned rather than entering a directional
bin.

## Why the first prediction can fail

Linearization about the zero state gives a Gaussian distribution during the
early part of the sweep. Its mean and covariance satisfy

$$
\dot m=J(t)m+h,\qquad
\dot C=J(t)C+CJ(t)^{\mathsf T}+2D\mathbb I.
$$

Assigning the direction of a Gaussian sample to a final equilibrium assumes
that its direction already determines the stored pattern. Continued coupling,
saturation and noise can change that assignment. In the coupled amplitudes,
growth of the first two components supplies a later force on the third through
$\kappa(q_1+q_2)$. The final pattern depends on that subsequent evolution.

The constructor evaluates the generator on coordinate and quadratic
observables. For a constant isotropic diffusion,

$$
\mathcal L\phi=F\cdot\nabla\phi+D\nabla^2\phi.
$$

Coordinate observables expose the nonlinear drift missing from the linear
description. Quadratic observables also expose missing diffusion:
$\mathcal L(q_i^2)$ contains $2D$. The corrected generator retains both terms.
The Gaussian prefix ends at the first sampling-grid crossing where the rms
nonlinear remainder reaches 0.05 times $\sqrt r$ times the rms amplitude. Its
samples then evolve under the complete stochastic generator. This threshold
defines the approximation being tested; the separate trajectory comparison
measures its accuracy.

## Freeze the correction, then test another system

With the `stochastic` dependencies installed, run:

```bash
python3 -B scripts/run_stochastic_passage_loop.py \
  --out-dir build/stochastic_passage \
  --n 4096 --workers 2
```

The two development calculations use coupled amplitudes. A correction is
retained only after a resolved failure of the angular prediction and successful
checks of the repaired predictions. The resulting `frozen_rule.json` fixes the
approximation and acceptance policy before evaluating new preparations and
rates. The second drift family is a cubic-anisotropy vector model,

$$
\dot q_i=(\mu-g_i)q_i-|q|^2q_i-\gamma q_i^3,
\qquad \gamma=0.5.
$$

Its coefficients are read from its own specification. They are not copied from
the first system. The source equations determine a new remainder and a new set
of final equilibria. Both numerical families are supplied constructions from
the earlier multi-mode passage investigation, not newly extracted experimental
models.

Each `prediction.json` is written before the corresponding reference outcomes
are evaluated. Repaired predictions use Euler–Maruyama evolution; references
use stochastic Heun evolution. Both have paired Brownian step refinement.
The final noisy hold tests whether the pattern survives after writing, and
unassigned or non-finite trajectories prevent acceptance.

The finite-rate rescaling is derived symbolically. For homogeneous cubic
drifts, $q=r^{1/4}v$ and $\tau=\mu/\sqrt r$, with the physical time origin
translated to $\mu=0$, give scaled signal $h/r^{3/4}$, diffusion $D/r$, and linear
splittings proportional to $G/\sqrt r$. The comparison varies rate at matched
scaled signal and noise. The final landscape and holding conditions remain
fixed.

## Construct a writing protocol

The first search asks for all three amplitudes to be positive with probability
at least 0.60. It chooses among new rates using only committed predictions.
When those predictions miss the target, run the bounded refinement:

```bash
python3 -B scripts/refine_stochastic_protocol.py \
  --parent build/stochastic_passage \
  --out-dir build/stochastic_protocol_design --workers 2
```

The target remains 0.60. At most two additional proposals halve the best
predicted rate, down to a fixed lower limit. Signal and writing noise follow the
derived scaling. A candidate is selected only when its predicted 95% lower
bound reaches the target. Its reference simulation then uses a fresh noise
sample; a failed reference stops the campaign.

Two additional controls change the rate alone or the writing noise and signal
alone. These distinguish a rate effect from the effect of rescaling the drive.
All controls preserve the final drift. The reports retain a failed design as
well as a successful one; a changed outcome is never promoted merely because
it is rare.

## Retain the change that mattered

The completed calculation first selected a rate of 0.0025. Its prediction and
reference probability were 0.568 and 0.575, both below the unchanged target of
0.60. Refining the rate to 0.00125 gave a prediction of 0.650 and a reference
probability of 0.667. Every reference uses 4,096 independent trajectories.

The controls distinguish two possible explanations for this gain. Changing the
rate alone gave 0.660; changing the writing noise and force alone gave 0.570.
The retained proposal therefore halves the rate while leaving the force, noise,
final drift and holding conditions unchanged. It is selected from the controls,
not assumed to work for every initial preparation.

```bash
python3 -B scripts/learn_verified_protocol.py \
  --parent build/stochastic_passage \
  --design build/stochastic_protocol_design \
  --out-dir build/stochastic_protocol_reuse
```

This command writes the conditional rate-change recipe before testing it on a
new base rate and a new initial Gaussian distribution. In the recorded run, the
recipe changes the new rate from 0.002 to 0.001. The committed prediction is
0.682; the fresh reference gives 0.702, with a 95% Wilson lower bound of 0.688.
The same goal is reached without another search over writing settings. No
trajectory changes its assigned pattern during the finite 80-unit noisy hold.
This finite observation establishes retention over the tested interval, not an
infinite lifetime.

An independent cluster replay uses the unchanged task and correction with
fresh trajectories. Its prediction is 0.693 and its reference probability is
0.690, with a 95% Wilson lower bound of 0.675. Both the probability comparison
and the integration-step refinement pass on the second machine.

The reusable knowledge has two parts: a stochastic prediction method that was
checked on both supplied model families, and a rate-only proposal supported by
controls and a new application within the coupled-amplitude family. The
matching method and available rate changes are supplied algorithms. The
experiment tests their conditional selection and reuse; it does not establish
that the constructor invented a new stochastic law.

## Apply the rule to a new declared task

The public command separates prediction from verification:

```bash
python3 -B -m fieldbridge memory passage --stage predict \
  --task path/to/passage_task.json \
  --rule build/stochastic_passage/frozen_rule.json \
  --out-dir build/new_passage

python3 -B -m fieldbridge memory passage --stage verify \
  --out-dir build/new_passage
```

The task uses schema `fieldbridge-stochastic-passage/1`, an explicit Itô
convention, the source model, sweep, diffusion, bias, initial mean and covariance,
and holding conditions. The numerical example generates complete task records
in `plan.json`. This backend currently accepts up to four real polynomial
coordinates with a fixed zero reference branch, symmetric linearization and
constant isotropic diffusion. Other equations need another derivation backend.

The verified consequence is a basin-probability prediction for specified
equations and a proposed protocol. Experimental realization and originality
require their own evidence. The reusable result is the calculated correction
and the conditions under which it predicts successfully.

Implementation: [stochastic generator and prediction](../../fieldbridge/memory/stochastic_passage.py),
[frozen study](../../scripts/run_stochastic_passage_loop.py),
[protocol refinement](../../scripts/refine_stochastic_protocol.py),
[conditional reuse](../../scripts/learn_verified_protocol.py),
[tests](../../tests/test_stochastic_passage.py).
