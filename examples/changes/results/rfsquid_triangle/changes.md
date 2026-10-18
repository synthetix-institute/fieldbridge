# Changes of a model: three rf-SQUID annealer qubits on an antiferromagnetic triangle

The supplied model: symmetric-write several critical directions; threshold-write reached.

18 single changes keep the states bounded (0 discarded because states run away); 14 of them are steps.

| change | symmetric write | threshold write | step |
| --- | --- | --- | --- |
| control = gam over [0.125, 2] | no write point | reached | yes |
| control = eps over [-1, 1] | distant branch | reached | yes |
| control = J over [0.0125, 0.2] | no write point | reached | yes |
| range of beta wider: [-19.5, 22.5] | several critical directions | reached |  |
| range of beta above: [2.5, 22.5] | no write point | reached | yes |
| range of beta below: [-19.5, 0.5] | no write point | reached | yes |
| remove '-gam*v1' from dv1/dt | several critical directions | oscillation |  |
| remove '- (th1 - eps)' from dv1/dt | no write point | reached | yes |
| remove '+ beta*sin(th1)' from dv1/dt | reached | reached | reaches symmetric-write |
| remove '- J*(th2 + th3)' from dv1/dt | reached | reached | reaches symmetric-write |
| remove '-gam*v2' from dv2/dt | several critical directions | oscillation |  |
| remove '- (th2 - eps)' from dv2/dt | no write point | reached | yes |
| remove '+ beta*sin(th2)' from dv2/dt | reached | reached | reaches symmetric-write |
| remove '- J*(th1 + th3)' from dv2/dt | reached | reached | reaches symmetric-write |
| remove '-gam*v3' from dv3/dt | several critical directions | oscillation |  |
| remove '- (th3 - eps)' from dv3/dt | no write point | reached | yes |
| remove '+ beta*sin(th3)' from dv3/dt | reached | reached | reaches symmetric-write |
| remove '- J*(th1 + th2)' from dv3/dt | reached | reached | reaches symmetric-write |

## Pairs of changes

118 pairs; 3 open.

- control = eps over [-1, 1], then remove '+ beta*sin(th1)' from dv1/dt: a reason for stopping that no single change produces (symmetric-write: one-sided write).
- remove '+ beta*sin(th1)' from dv1/dt, then control = eps over [-1, 1]: a reason for stopping that no single change produces (symmetric-write: one-sided write).
- remove '- J*(th2 + th3)' from dv1/dt, then remove '-gam*v1' from dv1/dt: a reason for stopping that no single change produces (symmetric-write: oscillation).

A write reached after a change is a property of the changed model, not of the supplied one.

Calculated and structural statements about a supplied model. They do not establish that the model describes a physical system, and they do not establish novelty; both need comparison with experiment and with the literature.
