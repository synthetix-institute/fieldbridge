# Changes of a model: coherent Ising machine, three DOPOs on an antiferromagnetic triangle

The supplied model: symmetric-write several critical directions; threshold-write reached.

10 single changes keep the states bounded (3 discarded because states run away); 9 of them are steps.

| change | symmetric write | threshold write | step |
| --- | --- | --- | --- |
| control = e over [0.05, 0.8] | no write point | reached | yes |
| range of a wider: [-21, 21] | several critical directions | reached |  |
| range of a above: [1, 21] | no write point | reached | yes |
| range of a below: [-21, -1] | no write point | reached | yes |
| remove 'a*x1' from dx1/dt | reached | reached | reaches symmetric-write |
| remove '- e*(x2 + x3)' from dx1/dt | reached | reached | reaches symmetric-write |
| remove 'a*x2' from dx2/dt | one-sided write | reached | yes |
| remove '- e*(x1 + x3)' from dx2/dt | reached | reached | reaches symmetric-write |
| remove 'a*x3' from dx3/dt | reached | reached | reaches symmetric-write |
| remove '- e*(x1 + x2)' from dx3/dt | reached | reached | reaches symmetric-write |

## Pairs of changes

72 pairs; 7 open.

- range of a below: [-21, -1], then remove '- x1**3' from dx1/dt: a reason for stopping that no single change produces (threshold-write: oscillation).
- range of a below: [-21, -1], then remove '- x2**3' from dx2/dt: a reason for stopping that no single change produces (threshold-write: oscillation).
- range of a below: [-21, -1], then remove '- x3**3' from dx3/dt: a reason for stopping that no single change produces (threshold-write: oscillation).
- remove '- e*(x2 + x3)' from dx1/dt, then remove '- e*(x1 + x3)' from dx2/dt: a reason for stopping that no single change produces (threshold-write: oscillation).
- remove '- e*(x2 + x3)' from dx1/dt, then remove '- e*(x1 + x2)' from dx3/dt: a reason for stopping that no single change produces (threshold-write: oscillation).
- remove '- e*(x1 + x3)' from dx2/dt, then remove '- e*(x2 + x3)' from dx1/dt: a reason for stopping that no single change produces (threshold-write: oscillation).
- remove '- e*(x1 + x3)' from dx2/dt, then remove '- e*(x1 + x2)' from dx3/dt: a reason for stopping that no single change produces (threshold-write: oscillation).

A write reached after a change is a property of the changed model, not of the supplied one.

Calculated and structural statements about a supplied model. They do not establish that the model describes a physical system, and they do not establish novelty; both need comparison with experiment and with the literature.
