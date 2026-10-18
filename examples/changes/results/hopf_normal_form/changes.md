# Changes of a model: Hopf normal form of an oscillator

The supplied model: symmetric-write oscillation; threshold-write single stable state.

10 single changes keep the states bounded (2 discarded because states run away); 5 of them are steps.

| change | symmetric write | threshold write | step |
| --- | --- | --- | --- |
| control = omega over [0.25, 4] | the control scales the drift | single stable state | yes |
| range of mu wider: [-10.5, 10.5] | oscillation | single stable state |  |
| range of mu above: [0.5, 10.5] | no write point | single stable state | yes |
| range of mu below: [-10.5, -0.5] | no write point | single stable state | yes |
| remove '- x * y**2 * omega' from dx/dt | oscillation | single stable state |  |
| remove '+ mu * omega * x' from dx/dt | oscillation | single stable state |  |
| remove '- omega * y' from dx/dt | several critical directions | single stable state | yes |
| remove '- x**2 * y * omega' from dy/dt | oscillation | single stable state |  |
| remove '+ mu * omega * y' from dy/dt | oscillation | single stable state |  |
| remove '+ omega * x' from dy/dt | several critical directions | single stable state | yes |

## Pairs of changes

71 pairs; 16 open.

- range of mu above: [0.5, 10.5], then remove '+ mu * omega * x' from dx/dt: reaches the threshold-write, which no single change reaches.
- range of mu above: [0.5, 10.5], then remove '- omega * y' from dx/dt: a reason for stopping that no single change produces (symmetric-write: one-sided write).
- range of mu above: [0.5, 10.5], then remove '+ mu * omega * y' from dy/dt: reaches the threshold-write, which no single change reaches.
- remove '- omega * y' from dx/dt, then range of mu above: [0.5, 10.5]: a reason for stopping that no single change produces (symmetric-write: one-sided write).
- remove '- omega * y' from dx/dt, then remove '+ mu * omega * x' from dx/dt: reaches the symmetric-write and the threshold-write, which no single change reaches.
- remove '- omega * y' from dx/dt, then remove '+ mu * omega * y' from dy/dt: reaches the threshold-write, which no single change reaches.
- remove '+ omega * x' from dy/dt, then remove '+ mu * omega * x' from dx/dt: reaches the symmetric-write and the threshold-write, which no single change reaches.
- remove '+ omega * x' from dy/dt, then remove '+ mu * omega * y' from dy/dt: reaches the threshold-write, which no single change reaches.
- range of mu wider: [-10.5, 10.5], then remove '+ mu * omega * x' from dx/dt: reaches the threshold-write, which no single change reaches.
- range of mu wider: [-10.5, 10.5], then remove '- x**2 * y * omega' from dy/dt: reaches the threshold-write, which no single change reaches.
- range of mu wider: [-10.5, 10.5], then remove '+ mu * omega * y' from dy/dt: reaches the threshold-write, which no single change reaches.
- remove '+ mu * omega * x' from dx/dt, then range of mu wider: [-10.5, 10.5]: reaches the threshold-write, which no single change reaches.
- remove '+ mu * omega * x' from dx/dt, then range of mu above: [0.5, 10.5]: reaches the threshold-write, which no single change reaches.
- remove '+ mu * omega * x' from dx/dt, then remove '- omega * y' from dx/dt: reaches the symmetric-write and the threshold-write, which no single change reaches.
- remove '+ mu * omega * x' from dx/dt, then remove '+ mu * omega * y' from dy/dt: a reason for stopping that no single change produces (threshold-write: oscillation).
- remove '+ mu * omega * x' from dx/dt, then remove '+ omega * x' from dy/dt: reaches the symmetric-write and the threshold-write, which no single change reaches.

A write reached after a change is a property of the changed model, not of the supplied one.

Calculated and structural statements about a supplied model. They do not establish that the model describes a physical system, and they do not establish novelty; both need comparison with experiment and with the literature.
