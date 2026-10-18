# Changes of a model: three Kerr parametric oscillators on an antiferromagnetic triangle

The supplied model: symmetric-write several critical directions; threshold-write reached.

19 single changes keep the states bounded (6 discarded because states run away); 17 of them are steps.

| change | symmetric write | threshold write | step |
| --- | --- | --- | --- |
| control = K over [1.25, 20] | no write point | reached | yes |
| control = Delta over [0.125, 2] | no write point | reached | yes |
| control = kappa over [0.25, 4] | no write point | reached | yes |
| control = J over [-2, -0.125] | no write point | reached | yes |
| range of p wider: [-120, 132] | several critical directions | reached |  |
| range of p above: [12, 132] | no write point | reached | yes |
| range of p below: [-120, 0] | several critical directions | reached |  |
| remove '- 0.5*kappa*x1' from dx1/dt | reached | reached | reaches symmetric-write |
| remove '- J*(y2 + y3)' from dx1/dt | reached | reached | reaches symmetric-write |
| remove '- 0.5*kappa*x2' from dx2/dt | reached | reached | reaches symmetric-write |
| remove '- J*(y3 + y1)' from dx2/dt | reached | reached | reaches symmetric-write |
| remove '- 0.5*kappa*x3' from dx3/dt | reached | reached | reaches symmetric-write |
| remove '- J*(y1 + y2)' from dx3/dt | reached | reached | reaches symmetric-write |
| remove '- 0.5*kappa*y1' from dy1/dt | reached | reached | reaches symmetric-write |
| remove '+ J*(x2 + x3)' from dy1/dt | vanishing cubic | reached | yes |
| remove '- 0.5*kappa*y2' from dy2/dt | reached | reached | reaches symmetric-write |
| remove '+ J*(x3 + x1)' from dy2/dt | vanishing cubic | reached | yes |
| remove '- 0.5*kappa*y3' from dy3/dt | reached | reached | reaches symmetric-write |
| remove '+ J*(x1 + x2)' from dy3/dt | vanishing cubic | reached | yes |

## Pairs of changes

144 pairs; 42 open.

- control = K over [1.25, 20], then range of K wider: [-186.2, 207.5]: a reason for stopping that no single change produces (symmetric-write: oscillation).
- control = K over [1.25, 20], then range of K below: [-186.2, 1.25]: a reason for stopping that no single change produces (symmetric-write: oscillation).
- control = K over [1.25, 20], then remove '- 0.5*kappa*x1' from dx1/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = K over [1.25, 20], then remove '- 0.5*kappa*x2' from dx2/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = K over [1.25, 20], then remove '- 0.5*kappa*x3' from dx3/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = K over [1.25, 20], then remove '- 0.5*kappa*y1' from dy1/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = K over [1.25, 20], then remove '- 0.5*kappa*y2' from dy2/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = K over [1.25, 20], then remove '- 0.5*kappa*y3' from dy3/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = Delta over [0.125, 2], then range of Delta below: [-18.62, 0.125]: a reason for stopping that no single change produces (symmetric-write: distant branch).
- control = Delta over [0.125, 2], then remove '- 0.5*kappa*x1' from dx1/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = Delta over [0.125, 2], then remove '- 0.5*kappa*x2' from dx2/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = Delta over [0.125, 2], then remove '- 0.5*kappa*x3' from dx3/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = Delta over [0.125, 2], then remove '- 0.5*kappa*y1' from dy1/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = Delta over [0.125, 2], then remove '- 0.5*kappa*y2' from dy2/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = Delta over [0.125, 2], then remove '- 0.5*kappa*y3' from dy3/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = kappa over [0.25, 4], then remove '- 0.5*kappa*x1' from dx1/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = kappa over [0.25, 4], then remove '- 0.5*kappa*x2' from dx2/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = kappa over [0.25, 4], then remove '- 0.5*kappa*x3' from dx3/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = kappa over [0.25, 4], then remove '- 0.5*kappa*y1' from dy1/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = kappa over [0.25, 4], then remove '- 0.5*kappa*y2' from dy2/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = kappa over [0.25, 4], then remove '- 0.5*kappa*y3' from dy3/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = J over [-2, -0.125], then remove '- 0.5*kappa*x1' from dx1/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = J over [-2, -0.125], then remove '- 0.5*kappa*x2' from dx2/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = J over [-2, -0.125], then remove '- 0.5*kappa*x3' from dx3/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = J over [-2, -0.125], then remove '- 0.5*kappa*y1' from dy1/dt: a reason for stopping that no single change produces (symmetric-write: distant branch; threshold-write: one-sided write).
- control = J over [-2, -0.125], then remove '+ J*(x2 + x3)' from dy1/dt: a reason for stopping that no single change produces (symmetric-write: distant branch).
- control = J over [-2, -0.125], then remove '- 0.5*kappa*y2' from dy2/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = J over [-2, -0.125], then remove '- 0.5*kappa*y3' from dy3/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove '- 0.5*kappa*x1' from dx1/dt, then control = K over [1.25, 20]: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove '- 0.5*kappa*x1' from dx1/dt, then control = Delta over [0.125, 2]: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove '- 0.5*kappa*x1' from dx1/dt, then control = kappa over [0.25, 4]: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove '- 0.5*kappa*x1' from dx1/dt, then control = J over [-2, -0.125]: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove '- 0.5*kappa*x1' from dx1/dt, then remove '+ J*(x2 + x3)' from dy1/dt: a reason for stopping that no single change produces (threshold-write: single stable state).
- remove '- 0.5*kappa*x1' from dx1/dt, then remove '- 0.5*kappa*y2' from dy2/dt: a reason for stopping that no single change produces (symmetric-write: one-sided write).
- remove '- 0.5*kappa*x1' from dx1/dt, then remove '- 0.5*kappa*y3' from dy3/dt: a reason for stopping that no single change produces (symmetric-write: one-sided write).
- remove '- 0.5*kappa*x2' from dx2/dt, then control = K over [1.25, 20]: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove '- 0.5*kappa*x2' from dx2/dt, then control = Delta over [0.125, 2]: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove '- 0.5*kappa*x2' from dx2/dt, then control = kappa over [0.25, 4]: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove '- 0.5*kappa*x2' from dx2/dt, then control = J over [-2, -0.125]: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove '- 0.5*kappa*x2' from dx2/dt, then remove '- 0.5*kappa*y1' from dy1/dt: a reason for stopping that no single change produces (symmetric-write: one-sided write).
- remove '- 0.5*kappa*x2' from dx2/dt, then remove '+ J*(x3 + x1)' from dy2/dt: a reason for stopping that no single change produces (threshold-write: oscillation).
- remove '- 0.5*kappa*x2' from dx2/dt, then remove '- 0.5*kappa*y3' from dy3/dt: a reason for stopping that no single change produces (symmetric-write: one-sided write).

A write reached after a change is a property of the changed model, not of the supplied one.

Calculated and structural statements about a supplied model. They do not establish that the model describes a physical system, and they do not establish novelty; both need comparison with experiment and with the literature.
