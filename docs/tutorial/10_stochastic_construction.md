# Why a change of stochastic coordinate adds a drift

[Visual constructor](https://synthetix-institute.github.io/fieldbridge/#stochastic):
detach the derived drift, build the target equation, and compare its mean with
the source mean. The equations and Python checks below supply that demonstration.

Consider a positive fluctuating coordinate X with restoring drift. In
dimensionless variables its Ito equation is

$$
dX_t=\left(\frac{\theta+1/2}{X_t}-\alpha X_t\right)dt+dW_t,
\qquad \theta>0,\quad\alpha>0.
$$

Here W is standard Brownian motion. Suppose a detector responds to the squared
coordinate Y=X squared. To predict that signal, the governing equation must be
transformed along with the coordinate. Brownian increments have quadratic
variation dW squared equal to dt, so the second derivative of the state map
contributes to the drift. This is the ordinary Ito rule, not an additional
empirical force. A derivation of the rule is provided in
[MIT's lecture on the Ito process and formula](https://live.ocw.mit.edu/courses/15-070j-advanced-stochastic-processes-fall-2013/d9d7372cbf65d56aa8aa9d59ba0ab2e8_MIT15_070JF13_Lec17.pdf).

For a source equation dX=a(X)dt+b(X)dW and a twice differentiable map h, the
generator acting on a smooth function f is L_X f=a f'+b squared f''/2. For a
target test function phi, direct differentiation gives

$$
L_X(\phi\circ h)
=\left(ah'+\frac{b^2}{2}h''\right)\phi'(h)
 +\frac{b^2(h')^2}{2}\phi''(h).
$$

The coefficients of the first and second derivatives determine the target
generator. The implementation obtains both from the source expressions,
substitutes the declared inverse map, and substitutes back to verify the two
coefficients. Checking only phi(y)=y would miss an incorrect diffusion term.

For h(x)=x squared, h'=2x and h''=2. Consequently,

$$
dY_t=(2\theta+2-2\alpha Y_t)dt+2\sqrt{Y_t}\,dW_t.
$$

Omitting quadratic variation would give 2 theta + 1 rather than 2 theta + 2.
On phi(y)=y, the difference between the source generator and the pulled-back
naive target generator is exactly 1. On a general phi the missing operator
is the first derivative with respect to y. The required correction is
therefore determined by the generator identity rather than fitted to data.

## Execute the calculation

From the repository root:

```bash
python3 -m pip install -e '.[construction]'
python3 -B -m fieldbridge verify-construction \
  examples/construction/ito_square.json --out-dir build/ito_square
```

The [input](../../examples/construction/ito_square.json) supplies the source
equation, map and deliberately incomplete candidate. The JSON output contains
the derived coefficients, zero residuals for the corrected equation, residual
1 for the incomplete candidate, and the generator applied to y and y squared.
It also records a hash of the exact input and the original assumptions.

For finite moments and a boundary realization with no additional contribution
to the mean equation, writing m(t)=E[Y_t] and m(0)=m0 gives

$$
m(t)=\frac{\theta+1}{\alpha}
 +\left(m_0-\frac{\theta+1}{\alpha}\right)e^{-2\alpha t}.
$$

With the same initial mean, the corrected and naive predictions differ by
(1-exp(-2 alpha t))/(2 alpha). This is an observable bias in the supplied
model, not evidence of a newly discovered law. No experimental comparison is
contained in the example.

The proof concerns interior differential expressions. Global stochastic
evolution also depends on boundary behavior and existence and uniqueness.
The program records this distinction rather than labelling the calculation
a global equivalence of stochastic processes.

For comparison, [the affine example](../../examples/construction/affine_transfer.json)
has h''=0. It verifies an exact rescaling without an extra drift. This control
shows why a correction is required by a particular transformation rather than
added to every transfer.

## The convention changes a measurable consequence

For the displayed coefficients $dX=\mu X\,dt+\sigma X\,dW$ and $Y=\log X$,
an Ito interpretation gives drift $\mu-\sigma^2/2$, whereas a Stratonovich
interpretation gives drift $\mu$. Both give variance rate $\sigma^2$.
For a fixed initial state the variance at time $t$ is $\sigma^2t$ in both
cases, so a variance measurement alone cannot detect this convention error.

```bash
python3 -B -m fieldbridge verify-construction \
  examples/construction/log_signal_ito.json --out-dir build/log_ito
python3 -B -m fieldbridge verify-construction \
  examples/construction/log_signal_stratonovich.json --out-dir build/log_stratonovich
```

The inputs contain no proposed target coefficients. Each declares the source
convention, positive source domain and real target domain. The verifier first
converts a Stratonovich drift $a$ to the Ito drift $a+bb'/2$ and then derives
the target generator, following the standard [conversion formula in
Pavliotis, Section 3.2](https://www.ma.imperial.ac.uk/~pavl/PavliotisBook.pdf).
Missing conventions are refused rather than guessed.
Changing only the convention changes the physical process; changing the
drift by the conversion formula preserves it. Consistently renaming a
coefficient must also preserve the prediction.

The [convention tests](../../tests/test_stochastic_conventions.py) integrate
the source coordinate using Euler-Maruyama for Ito and stochastic Heun for
Stratonovich, then compare logarithmic moments with the symbolic predictions.
These are authored numerical controls, not an evaluation of extraction from
held-out papers. The output key `target.variance` denotes the local quadratic
variation rate, not the variance of the state at a finite time. Only in this
constant-drift, constant-noise logarithmic example does multiplying it by
elapsed time give that variance directly for a fixed initial state.

Source: [`ito_transfer`](../../fieldbridge/verification.py).
Tests: [`test_verification.py`](../../tests/test_verification.py).
Next: [quantum observable closure](11_quantum_closure.md).

## Exercise: change the noise, predict the correction

Keep the source drift and map fixed, but double the Brownian amplitude.
Before running anything, use the factor $b^2 h''/2$ to predict the additional
drift. It should become 4, not 2: quadratic variation scales with the square
of the noise amplitude.

This creates a separate input and leaves the original example unchanged:

```bash
python3 - <<'PY'
import json
from pathlib import Path

spec = json.loads(Path('examples/construction/ito_square.json').read_text(encoding='utf-8'))
spec['noise'] = '2'
spec['question'] = 'How does doubling the Brownian amplitude change the squared signal?'
spec.pop('candidate')
spec['assumptions'] = [text.replace('standard scalar Ito Brownian motion',
    'standard scalar Ito Brownian motion with noise coefficient 2')
    for text in spec['assumptions']]
spec['provenance']['origin'] = 'tutorial perturbation of an authored benchmark'
spec['provenance']['equation_ids'] = ['tutorial_double_noise:source', 'tutorial_double_noise:map']
out = Path('build/tutorial_double_noise.json')
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(spec, indent=2), encoding='utf-8')
PY
python3 -B -m fieldbridge verify-construction \
  build/tutorial_double_noise.json --out-dir build/tutorial_double_noise
```

The derived drift is $2\theta+5-2\alpha y$, the variance is $16y$, and
the omission residual is 4. This is an interior generator calculation; any
global process interpretation must reconsider the boundary with the changed
noise. The experiment shows that the correction follows the supplied noise
law rather than a fixed example answer.

[Tutorial](index.md)
