# Turn a match into a physical question

The Brownian-motion query and a noisy optimization rule share a possible
gradient-diffusion description. To investigate it, specify the correspondence:
particle position becomes a parameter coordinate, the potential becomes a
loss, and the noise law and time convention must be matched. Boundary behavior
and the observable remain part of the question.

```bash
python3 -B -m fieldbridge construct examples/brownian_probability_flow.tex \
  --to stochastic_optimization --no-hyperion
```

The ordinary constructor assembles the extracted description with the
retrieved target example. Read “Preserved Contract” in this output as a
**proposed** preserved relation. The code's content-presence checks do not
establish conservation or equality of predictions.

## Decide what has to be calculated

For a stochastic coordinate map $y=h(x)$, a precise question is whether the
two generators agree on every smooth target observable:

```math
L_X(\phi\circ h)=(L_Y\phi)\circ h .
```

Here $L_X$ and $L_Y$ generate expectations in the two descriptions.
The map and its domain are specified before calculating the difference.
A mismatch then refers to a definite relation, rather than to a general
difference between responses.

In the [squared-coordinate example](10_stochastic_construction.md), direct
differentiation determines the target coefficients. The omitted Itô drift
leaves residual 1 on $\phi(y)=y$, so the error has an immediate consequence
for the measured mean. In contrast, ordinary translation only supplies a
possible target equation.

## Two constructor outputs

| Mode | Input beyond the query | Result |
| --- | --- | --- |
| `construct` | Target field and stored field records | A proposal, attachments and suggested tests |
| `construct --calculate` | Retrieved mathematical annotation and explicit correspondence | Derived equation or quantum closure, with calculated residuals |

The calculation adapter reads the selected source from the actual retrieval
results. It checks the annotation against the stored canonical equation and
then invokes [verification.py](../../fieldbridge/verification.py).
See [the complete command](13_retrieval_to_calculation.md).

For a new problem, write the intended relation before deciding whether the
calculation passed. A changing response may be the correct prediction of a
changed parameter. An obstruction is a failure of the relation that was
claimed to survive that change.

[Run a calculated construction](08_end_to_end_walkthrough.md) · [Tutorial](index.md)
