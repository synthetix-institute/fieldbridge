# Mechanisms of quantum-computing hardware

Written by `python3 -B examples/quantum/semiclassical/survey/report.py` from the results of `survey.py`.

The census lists 60 hardware elements described in reviews of superconducting, semiconductor, trapped-ion, neutral-atom, photonic and optomechanical platforms: 21 have a classical write, 15 have one in part and 24 have none (gates, codes and linear amplifiers store or process information without a write). 31 models represent the elements with a write; the status of every element is in `census.json`.

## Derivations and laws at the obstructions

| model | symmetric write | threshold write | phase locking | laws at the obstructions | a single change that reaches a write |
| --- | --- | --- | --- | --- | --- |
| bright_state_resonator | fold; no cusp with a supercritical sweep | reached, no law | time limit of 1800 s | second parameter (dq): no cusp found; second parameter (dc): no cusp found; second parameter (dq): no cusp found; second parameter (dc): no cusp found |  |
| cim_dopo | known law, 1.195 ± 0.150 | a single stable state: nothing to switch to | no oscillation in the control range | none: response 1.35, offset -0.02; full dynamics 1.3, -0.01; certified | remove '-0.5*s' from ds/dt (threshold write) |
| cim_pair | known law, 1.389 ± 0.144 | known law, 1.019 ± 0.000 | no oscillation in the control range | none: response 1.35, offset -0.02; full dynamics 1.31, -0.02; certified; distant branch: response 1.34, offset -0.03; full dynamics 1.39, 0.0; certified; distant branch: response 1.33, offset -0.02; full dynamics 1.37, 0.01; certified |  |
| cim_triangle | degenerate critical mode: 2 directions | known law, 1.019 ± 0.000 | no oscillation in the control range | critical space: 6 states, largest deviation 0.038, not certified | remove '- e*(x2 + x3)' from dx1/dt (symmetric write); remove 'a*x2' from dx2/dt (symmetric write) |
| dissipative_cat | known law, 1.331 ± 0.126 | known law, 1.019 ± 0.000 | no oscillation in the control range | none: response 1.35, offset -0.02; full dynamics 1.24, -0.0; certified |  |
| flux_qubit_tilt | known law, 1.346 ± 0.139 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (beta): response 1.21, offset -0.01; full dynamics 1.32, -0.02; not certified; second parameter (beta): response 1.19, offset -0.0; full dynamics 1.3, -0.02; not certified |  |
| fluxonium | known law, 1.171 ± 0.112 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (ell): response 1.27, offset -0.0; full dynamics 1.32, -0.02; certified; second parameter (ell): response 1.27, offset -0.0; full dynamics 1.33, 0.03; certified |  |
| ion_duffing | known law, 1.474 ± 0.163 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (k): response 1.33, offset -0.99; full dynamics 1.34, -0.85; not certified; second parameter (k): response 1.33, offset -0.98; full dynamics 1.34, -0.83; not certified |  |
| ion_phonon_laser | law differs, 0.216 ± 0.299 | known law, 1.019 ± 0.000 | known law, 1.000 ± 0.000 | second parameter (Delta): response 0.64, offset 0.45; full dynamics 0.78, 0.83; not certified |  |
| jba_resonator | known law, 1.267 ± 0.133 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (chi): response 1.32, offset -1.03; full dynamics 1.33, -0.8; not certified; second parameter (chi): response 1.31, offset -1.01; full dynamics 1.33, -0.8; not certified |  |
| jpc_nondegenerate | degenerate critical mode: 2 directions | a single stable state: nothing to switch to | no oscillation in the control range | critical space: continuous family of written states | remove '+ chi*(cx*bx + cy*by)' from dax/dt (symmetric write); remove '+ chi*(cy*bx - cx*by)' from day/dt (symmetric write) |
| kerr_cat_kpo | vanishing cubic: leading term of order 5 | reached, no law | no oscillation in the control range | leading order: response 1.35, offset -0.02; full dynamics 0.04, 0.01; not certified |  |
| kpo_pair | vanishing cubic: leading term of order 5 | known law, 1.019 ± 0.000 | no oscillation in the control range | leading order: response 1.34, offset -0.01; full dynamics 0.06, -0.01; not certified; distant branch: response 1.36, offset -0.02; full dynamics 1.34, 0.01; certified; distant branch: response 1.33, offset -0.01; full dynamics 1.36, 0.0; certified | remove '- 0.5*kappa*x1' from dx1/dt (symmetric write); remove '- J*y2' from dx1/dt (symmetric write) |
| kpo_triangle | degenerate critical mode: 2 directions | known law, 1.019 ± 0.000 | no oscillation in the control range | critical space: 6 states, largest deviation 0.69, not certified; critical space: 6 states, largest deviation 0.045, not certified; critical space: 6 states, largest deviation 0.044, not certified | remove '- 0.5*kappa*x1' from dx1/dt (symmetric write); remove '- J*(y2 + y3)' from dx1/dt (symmetric write) |
| magnon_kerr | known law, 1.780 ± 0.220 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (Dm): response 0.0, offset -3.72; full dynamics 1.3, -1.06; not certified; critical space: 2 states, largest deviation 0.0, certified |  |
| nuclear_bath | reached, no law | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (p): response 1.34, offset -2.28; full dynamics 1.32, -2.3; certified; second parameter (p): response 1.27, offset -2.22; full dynamics 1.28, -2.21; certified |  |
| optomech_bistable | reached, no law | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (Dl): response 0.0, offset 3.72; full dynamics 0.35, 2.46; not certified |  |
| optomech_selfosc | Hopf: an oscillation, not a write | no stable state: the realization oscillates | known law, 1.000 ± 0.005 | phase: phase |  |
| parametric_amplifier | no write point along the control | a single stable state: nothing to switch to | no oscillation in the control range | no write point along the control | remove '- 0.5*kappa*x' from dx/dt (threshold write); remove '- 0.5*kappa*y' from dy/dt (threshold write) |
| phase_qubit_rfsquid | degenerate critical mode: 2 directions | known law, 1.019 ± 0.000 | no oscillation in the control range | critical space: 4 states, largest deviation 0.0, certified; critical space: 4 states, largest deviation 0.012, certified; critical space: 4 states, largest deviation 0.004, certified |  |
| po_amplitude_oscillator | subcritical: jump to a distant branch | a single stable state: nothing to switch to | no oscillation in the control range | distant branch: response 1.35, offset -0.02; full dynamics 1.33, 0.01; certified | remove '+ alpha*(x**2 + y**2)*x' from dy/dt (symmetric write) |
| pplo_oscillator | a line of states: the drift vanishes along the critical direction | a single stable state: nothing to switch to | no oscillation in the control range | leading order: response 1.35, offset -0.02; full dynamics 1.32, -0.02; certified | remove '-0.5*kappa*qx' from dqx/dt (threshold write); remove '-0.5*kappa*qy' from dqy/dt (threshold write) |
| qfp_latch | known law, 1.623 ± 0.121 | a single stable state: nothing to switch to | no oscillation in the control range | none: response 1.35, offset -0.02; full dynamics 1.34, -0.04; certified |  |
| quantum_vdp | no write point along the control | no stable state: the realization oscillates | known law, 1.000 ± 0.000 | no write point along the control | range of F wider: [-20, 22] (symmetric write, threshold write); range of F below: [-20, 0] (symmetric write, threshold write) |
| rfsquid_annealer | known law, 1.248 ± 0.111 | known law, 1.019 ± 0.000 | no oscillation in the control range | none: response 1.34, offset -0.02; full dynamics 1.31, -0.01; certified |  |
| rfsquid_pair | known law, 1.289 ± 0.112 | known law, 1.019 ± 0.000 | no oscillation in the control range | none: response 1.23, offset -0.01; full dynamics 1.32, 0.01; not certified; distant branch: response 1.38, offset -0.01; full dynamics 1.54, -0.0; not certified; distant branch: response 1.4, offset 0.01; full dynamics 1.53, 0.01; not certified |  |
| rfsquid_triangle | degenerate critical mode: 2 directions | known law, 1.019 ± 0.000 | no oscillation in the control range | critical space: 6 states, largest deviation 0.056, not certified | remove '+ beta*sin(th1)' from dv1/dt (symmetric write); remove '- J*(th2 + th3)' from dv1/dt (symmetric write) |
| rydberg_bistability | known law, 1.465 ± 0.128 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (Om): response 1.31, offset 0.2; full dynamics 1.36, 0.08; not certified |  |
| snail_po | known law, 1.169 ± 0.129 | known law, 1.016 ± 0.004 | no oscillation in the control range | none: response 1.34, offset -0.02; full dynamics 1.32, -0.02; certified |  |
| spin_ensemble_cavity | reached, no law | known law, 1.012 ± 0.033 | no oscillation in the control range | second parameter (dc): no choice at the end of the sweep; second parameter (ds): no cusp found |  |
| transmon_driven | known law, 1.542 ± 0.150 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (Delta): response 1.35, offset 1.02; full dynamics 1.35, 0.82; not certified; second parameter (Delta): response 1.37, offset 1.01; full dynamics 1.35, 0.84; not certified |  |

## Counts

- law at an obstruction: certified: 20
- law at an obstruction: continuous family of written states: 1
- law at an obstruction: no choice at the end of the sweep: 1
- law at an obstruction: no cusp found: 5
- law at an obstruction: not certified: 22
- law at an obstruction: phase: 1
- phase locking: known law: 3
- phase locking: not derived within the time limit: 1
- phase locking: stops at an obstruction: 27
- symmetric write: known law: 14
- symmetric write: reached: 4
- symmetric write: stops at an obstruction: 13
- threshold write: known law: 21
- threshold write: reached: 2
- threshold write: stops at an obstruction: 8

## Sweeps through a cusp: response to the bias and intrinsic offset

Response constant and offset at a bias of one noise floor, from the law of the completed model and from the full stochastic dynamics.

| model | second parameter | law | full dynamics |
| --- | --- | --- | --- |
| flux_qubit_tilt | beta | 1.21, -0.01 | 1.32, -0.02 |
| flux_qubit_tilt | beta | 1.19, -0.0 | 1.3, -0.02 |
| fluxonium | ell | 1.27, -0.0 | 1.32, -0.02 |
| fluxonium | ell | 1.27, -0.0 | 1.33, 0.03 |
| ion_duffing | k | 1.33, -0.99 | 1.34, -0.85 |
| ion_duffing | k | 1.33, -0.98 | 1.34, -0.83 |
| ion_phonon_laser | Delta | 0.64, 0.45 | 0.78, 0.83 |
| jba_resonator | chi | 1.32, -1.03 | 1.33, -0.8 |
| jba_resonator | chi | 1.31, -1.01 | 1.33, -0.8 |
| magnon_kerr | Dm | 0.0, -3.72 | 1.3, -1.06 |
| nuclear_bath | p | 1.34, -2.28 | 1.32, -2.3 |
| nuclear_bath | p | 1.27, -2.22 | 1.28, -2.21 |
| optomech_bistable | Dl | 0.0, 3.72 | 0.35, 2.46 |
| rydberg_bistability | Om | 1.31, 0.2 | 1.36, 0.08 |
| transmon_driven | Delta | 1.35, 1.02 | 1.35, 0.82 |
| transmon_driven | Delta | 1.37, 1.01 | 1.35, 0.84 |

## Writes into a critical plane

| model | written states | largest deviation | certified |
| --- | ---: | ---: | --- |
| cim_triangle | 6 | 0.038 | no |
| kpo_triangle | 6 | 0.69 | no |
| kpo_triangle | 6 | 0.045 | no |
| kpo_triangle | 6 | 0.044 | no |
| magnon_kerr | 2 | 0.0 | yes |
| phase_qubit_rfsquid | 4 | 0.0 | yes |
| phase_qubit_rfsquid | 4 | 0.012 | yes |
| phase_qubit_rfsquid | 4 | 0.004 | yes |
| rfsquid_triangle | 6 | 0.056 | no |

Calculated statements about the classical images of the hardware elements, not about the devices; the quantum dynamics at a few photons can differ (Module 15).
