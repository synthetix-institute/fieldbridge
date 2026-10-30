# Co-discovery by construction: the symmetric write

Target: symmetric write (supercritical pitchfork, ds/dt = eps s - s^3 + h). Reached in 1 realizations from 1 fields (chemical reaction dynamics) by 1 classes of derivation; obstructed in 0.

| realization | field | derivation | status | write point | even part | law constant |
|---|---|---|---|---|---:|---:|
| post-transition-state bifurcation, transverse coordinate | chemical reaction dynamics | SCRKL | reached | s = 0.25 | 2.8e-16 | 1.289 ± 0.112 |

## Classes of derivation

- S(reflection) C R(pitchfork) K L: post-transition-state bifurcation, transverse coordinate

## Derivations

- post-transition-state bifurcation, transverse coordinate (chemical reaction dynamics): S[reflection y -> -y] > C[s] > R[supercritical pitchfork] > K > L

## Obstructions


## Invariants

- Canonical form at the write point: cubic sign -1 in every reached realization; largest even part 2.8e-16 of the cubic term at the edge of the window.
- Constant of the swept-write law (400 trajectories per realization): 1.289 ± 0.112 over 1 models from 1 realizations; pi^(1/4) = 1.331; chi^2 = 0.1 for 1 values.

## Letters and the slots they act on

| letter | transformation | slots | action |
|---|---|---|---|
| S | symmetry | Xi, Omega | a linear map of the carrier that leaves the drift invariant and reverses the critical mode |
| C | continuation | A | the control is moved to the value where the linear rate kappa of the critical mode vanishes |
| R | reduction | Xi, Omega | the carrier is reduced to the critical line; the drift to a0 + a1 s + a2 s^2 + a3 s^3 |
| U | unfolding | A, P | a second material parameter is set to the cusp; the protocol sweeps through the cusp |
| K | canonical form | Xi | x = |a3|^(1/2) s, so that the reduced drift is eps x - x^3 + h |
| L | law | P, R | a swept write with a bias is simulated; the observable is the sign of the critical coordinate |

Calculated and structural statements about a supplied model. They do not establish that the model describes a physical system, and they do not establish novelty; both need comparison with experiment and with the literature.
