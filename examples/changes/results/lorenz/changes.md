# Changes of a model: Lorenz model of convection in the chaotic range

The supplied model: symmetric-write no write point; threshold-write oscillation.

6 single changes keep the states bounded (4 discarded because states run away); 2 of them are steps.

| change | symmetric write | threshold write | step |
| --- | --- | --- | --- |
| control = sigma over [2.5, 40] | no write point | oscillation |  |
| control = b over [0.6667, 10.67] | no write point | oscillation |  |
| range of r below: [-327.9, 24.74] | reached | oscillation | reaches symmetric-write |
| remove 'r*x' from dy/dt | no write point | single stable state | yes |
| remove '-y' from dy/dt | no write point | oscillation |  |
| remove '-b*z' from dz/dt | no write point | oscillation |  |

## Pairs of changes

33 pairs; 7 open.

- range of r below: [-327.9, 24.74], then remove '-b*z' from dz/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove 'r*x' from dy/dt, then remove '-b*z' from dz/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- control = sigma over [2.5, 40], then range of sigma above: [40, 415]: a reason for stopping that no single change produces (symmetric-write: oscillation).
- control = sigma over [2.5, 40], then remove '-b*z' from dz/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove '-b*z' from dz/dt, then control = sigma over [2.5, 40]: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove '-b*z' from dz/dt, then range of r below: [-327.9, 24.74]: a reason for stopping that no single change produces (threshold-write: one-sided write).
- remove '-b*z' from dz/dt, then remove 'r*x' from dy/dt: a reason for stopping that no single change produces (threshold-write: one-sided write).

A write reached after a change is a property of the changed model, not of the supplied one.

Calculated and structural statements about a supplied model. They do not establish that the model describes a physical system, and they do not establish novelty; both need comparison with experiment and with the literature.
