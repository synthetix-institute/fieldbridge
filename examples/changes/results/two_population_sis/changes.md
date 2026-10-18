# Changes of a model: two symmetric populations with pairwise and three-body infection

The supplied model: symmetric-write no write point; threshold-write reached.

17 single changes keep the states bounded (0 discarded because states run away); 8 of them are steps.

| change | symmetric write | threshold write | step |
| --- | --- | --- | --- |
| control = beta over [0.9625, 15.4] | no write point | reached |  |
| control = gamma over [0.25, 4] | no write point | reached |  |
| control = d_1 over [1, 16] | one-sided write | reached | yes |
| control = d_2 over [0.75, 12] | reached | reached | reaches symmetric-write |
| control = d_3 over [0.75, 12] | reached | reached | reaches symmetric-write |
| control = d_4 over [1, 16] | reached | reached | reaches symmetric-write |
| control = h_1 over [0.25, 4] | reached | reached | reaches symmetric-write |
| control = h_2 over [0.25, 4] | reached | reached | reaches symmetric-write |
| range of tau wider: [-1, 1.1] | no write point | reached |  |
| range of tau above: [0.1, 1.1] | no write point | reached |  |
| range of tau below: [-1, 0] | no write point | reached |  |
| remove 'tau*(1-y_1)*(d_1*y_1+d_2*y_2)' from dy_1/dt | no write point | reached |  |
| remove '+beta*h_1*(1-y_1)*y_1**2' from dy_1/dt | one-sided write | reached | yes |
| remove '-gamma*y_1' from dy_1/dt | no write point | single stable state |  |
| remove 'tau*(1-y_2)*(d_3*y_1+d_4*y_2)' from dy_2/dt | no write point | reached |  |
| remove '+beta*h_2*(1-y_2)*y_2**2' from dy_2/dt | one-sided write | reached | yes |
| remove '-gamma*y_2' from dy_2/dt | no write point | single stable state |  |

## Pairs of changes

98 pairs; 0 open.


A write reached after a change is a property of the changed model, not of the supplied one.

Calculated and structural statements about a supplied model. They do not establish that the model describes a physical system, and they do not establish novelty; both need comparison with experiment and with the literature.
