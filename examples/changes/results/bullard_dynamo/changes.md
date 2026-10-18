# Changes of a model: modified Bullard disk dynamo

The supplied model: symmetric-write no write point; threshold-write reached.

12 single changes keep the states bounded (0 discarded because states run away); 2 of them are steps.

| change | symmetric write | threshold write | step |
| --- | --- | --- | --- |
| control = q over [0.25, 4] | no write point | reached |  |
| control = m over [0.25, 4] | no write point | reached |  |
| control = g over [0.25, 4] | reached | reached | reaches symmetric-write |
| control = f over [0.25, 4] | reached | reached | reaches symmetric-write |
| range of Z0 wider: [-19, 23] | no write point | reached |  |
| range of Z0 above: [3, 23] | no write point | reached |  |
| range of Z0 below: [-19, 1] | no write point | reached |  |
| remove 'x * Z' from dy/dt | no write point | single stable state |  |
| remove '+ m * x' from dy/dt | no write point | single stable state |  |
| remove '- (m + 1) * y' from dy/dt | no write point | reached |  |
| remove 'g * (1 - (m + 1) * x * y + m * x ** 2)' from dZ/dt | no write point | single stable state |  |
| remove '- f * Z' from dZ/dt | no write point | one-sided write |  |

## Pairs of changes

76 pairs; 0 open.


A write reached after a change is a property of the changed model, not of the supplied one.

Calculated and structural statements about a supplied model. They do not establish that the model describes a physical system, and they do not establish novelty; both need comparison with experiment and with the literature.
