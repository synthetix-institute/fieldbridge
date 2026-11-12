# Mechanisms of quantum-computing hardware

Written by `python3 -B examples/quantum/semiclassical/survey/report.py` from the results of `survey.py` (code stamp cd3149c7, versions {"python": "3.11.3", "numpy": "1.26.0", "scipy": "1.14.1"}).

The census lists 60 hardware elements described in reviews of superconducting, semiconductor, trapped-ion, neutral-atom, photonic and optomechanical platforms: 21 have a classical write, 15 have one in part and 24 have none (gates, codes and linear amplifiers store or process information without a write). 31 models represent the elements with a write; the status of every element is in `census.json`.

## Derivations and laws at the obstructions

| model | symmetric write | threshold write | phase locking | laws at the obstructions | a single change that reaches a write |
| --- | --- | --- | --- | --- | --- |
| bright_state_resonator | fold; no cusp with a supercritical sweep | reached, no law | no oscillation in the control range | second parameter (dq): no cusp found; second parameter (dc): no cusp found; second parameter (dq): no cusp found; second parameter (dc): no cusp found |  |
| cim_dopo | known law, 1.195 ± 0.150 | a single stable state: nothing to switch to | no oscillation in the control range | none: response 1.35, offset -0.02; full dynamics 1.36, -0.01; certified | remove '-0.5*s' from ds/dt (threshold write) |
| cim_pair | known law, 1.389 ± 0.144 | known law, 1.019 ± 0.000 | no oscillation in the control range | none: response 1.35, offset -0.02; full dynamics 1.34, -0.03; certified; distant branch: response 1.38, offset -0.01; full dynamics 1.37, 0.01; certified; distant branch: response 1.4, offset 0.02; full dynamics 1.39, 0.02; certified |  |
| cim_triangle | degenerate critical mode: 2 directions | known law, 1.019 ± 0.000 | no oscillation in the control range | critical space: 6 states, largest deviation 0.026, certified | remove '- e*(x2 + x3)' from dx1/dt (symmetric write); remove 'a*x2' from dx2/dt (symmetric write) |
| dissipative_cat | known law, 1.331 ± 0.126 | known law, 1.019 ± 0.000 | no oscillation in the control range | none: 2 states, largest deviation 0.031, not certified |  |
| flux_qubit_tilt | known law, 1.346 ± 0.139 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (beta): response 1.33, offset -0.02; full dynamics 1.33, 0.01; certified; second parameter (beta): response 1.31, offset -0.0; full dynamics 1.3, 0.0; certified |  |
| fluxonium | known law, 1.277 ± 0.133 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (ell): response 1.32, offset -0.02; full dynamics 1.34, 0.02; certified; second parameter (ell): response 1.31, offset 0.01; full dynamics 1.32, -0.0; certified |  |
| ion_duffing | law differs, 0.904 ± 0.137 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (k): response 1.36, offset -0.8; full dynamics 1.32, -0.78; certified; second parameter (k): response 1.34, offset -0.8; full dynamics 1.31, -0.78; certified |  |
| ion_phonon_laser | law differs, 0.124 ± 0.275 | known law, 1.019 ± 0.000 | known law, 1.000 ± 0.000 | second parameter (Delta): response 0.74, offset 0.31; full dynamics 0.72, 0.8; not certified |  |
| jba_resonator | known law, 1.379 ± 0.153 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (chi): response 1.36, offset -0.8; full dynamics 1.33, -0.79; certified; second parameter (chi): response 1.37, offset -0.81; full dynamics 1.33, -0.79; certified |  |
| jpc_nondegenerate | degenerate critical mode: 2 directions | a single stable state: nothing to switch to | no oscillation in the control range | critical space: continuous family of written states | remove '+ chi*(cx*bx + cy*by)' from dax/dt (symmetric write); remove '+ chi*(cy*bx - cx*by)' from day/dt (symmetric write) |
| kerr_cat_kpo | vanishing cubic: leading term of order 5 | reached, no law | no oscillation in the control range | leading order: response 1.35, offset -0.02; full dynamics 0.08, 0.01; not certified |  |
| kpo_pair | vanishing cubic: leading term of order 5 | reached, no law | no oscillation in the control range | leading order: response 1.35, offset -0.0; full dynamics 0.07, -0.01; not certified; distant branch: response 1.36, offset -0.02; full dynamics 1.37, -0.02; certified; distant branch: response 1.36, offset -0.01; full dynamics 1.3, 0.01; certified | remove '- 0.5*kappa*x1' from dx1/dt (symmetric write); remove '- J*y2' from dx1/dt (symmetric write) |
| kpo_triangle | degenerate critical mode: 2 directions | reached, no law | no oscillation in the control range | critical space: 6 states, largest deviation 0.694, not certified; critical space: 6 states, largest deviation 0.039, not certified; critical space: 6 states, largest deviation 0.044, not certified | remove '- 0.5*kappa*x1' from dx1/dt (symmetric write); remove '- J*(y2 + y3)' from dx1/dt (symmetric write) |
| magnon_kerr | known law, 1.202 ± 0.185 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (Dm): response 1.37, offset -1.15; full dynamics 1.3, -1.03; not certified; second parameter (Dm): response 1.37, offset -1.13; full dynamics 1.31, -1.03; not certified |  |
| nuclear_bath | reached, no law | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (p): response 1.35, offset -2.28; full dynamics 1.32, -2.3; certified; second parameter (p): response 1.29, offset -2.2; full dynamics 1.27, -2.14; certified |  |
| optomech_bistable | reached, no law | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (Dl): 2 states, largest deviation 0.297, not certified |  |
| optomech_selfosc | Hopf: an oscillation, not a write | no stable state: the realization oscillates | reached, no law | phase: phase |  |
| parametric_amplifier | no write point along the control | a single stable state: nothing to switch to | no oscillation in the control range | no write point along the control | remove '- 0.5*kappa*x' from dx/dt (threshold write); remove '- 0.5*kappa*y' from dy/dt (threshold write) |
| phase_qubit_rfsquid | known law, 1.363 ± 0.167 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (ell): response 1.31, offset -0.02; full dynamics 1.3, 0.01; certified; second parameter (ell): response 1.34, offset -0.01; full dynamics 1.33, 0.0; certified; second parameter (ell): response 1.31, offset 0.0; full dynamics 1.33, -0.03; certified |  |
| po_amplitude_oscillator | subcritical: jump to a distant branch | a single stable state: nothing to switch to | no oscillation in the control range | distant branch: response 1.35, offset -0.02; full dynamics 1.32, -0.03; certified | remove '+ alpha*(x**2 + y**2)*x' from dy/dt (symmetric write) |
| pplo_oscillator | vanishing cubic: leading term of order 5 | a single stable state: nothing to switch to | no oscillation in the control range | leading order: response 1.35, offset -0.02; full dynamics 1.33, -0.02; certified | remove '-0.5*kappa*qx' from dqx/dt (threshold write); remove '-0.5*kappa*qy' from dqy/dt (threshold write) |
| qfp_latch | known law, 1.623 ± 0.121 | a single stable state: nothing to switch to | no oscillation in the control range | none: response 1.35, offset -0.02; full dynamics 1.33, -0.0; certified |  |
| quantum_vdp | no write point along the control | no stable state: the realization oscillates | known law, 1.000 ± 0.000 | no write point along the control | range of F wider: [-20, 22] (symmetric write, threshold write); range of F below: [-20, 0] (symmetric write, threshold write) |
| rfsquid_annealer | known law, 1.248 ± 0.111 | known law, 1.019 ± 0.000 | no oscillation in the control range | none: response 1.33, offset -0.02; full dynamics 1.33, -0.02; certified |  |
| rfsquid_pair | known law, 1.289 ± 0.112 | known law, 1.019 ± 0.000 | no oscillation in the control range | none: response 1.34, offset -0.0; full dynamics 1.32, 0.01; not certified; distant branch: response 1.48, offset -0.0; full dynamics 1.51, -0.0; certified; distant branch: response 1.5, offset 0.0; full dynamics 1.56, -0.05; certified |  |
| rfsquid_triangle | degenerate critical mode: 2 directions | known law, 1.019 ± 0.000 | no oscillation in the control range | critical space: 8 states, largest deviation 0.024, not certified | remove '+ beta*sin(th1)' from dv1/dt (symmetric write); remove '- J*(th2 + th3)' from dv1/dt (symmetric write) |
| rydberg_bistability | known law, 1.329 ± 0.127 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (Om): response 1.36, offset 0.14; full dynamics 1.33, 0.06; not certified |  |
| snail_po | known law, 1.389 ± 0.114 | known law, 1.016 ± 0.004 | no oscillation in the control range | none: response 1.34, offset -0.02; full dynamics 1.28, -0.03; certified |  |
| spin_ensemble_cavity | fold; no cusp with a supercritical sweep | known law, 1.012 ± 0.033 | no oscillation in the control range | second parameter (dc): no cusp found; second parameter (ds): no cusp found |  |
| transmon_driven | known law, 1.542 ± 0.150 | known law, 1.019 ± 0.000 | no oscillation in the control range | second parameter (Delta): response 1.37, offset 0.78; full dynamics 1.3, 0.78; certified; second parameter (Delta): response 1.34, offset 0.78; full dynamics 1.36, 0.8; certified |  |

## Counts

- law at an obstruction: certified: 29
- law at an obstruction: continuous family of written states: 1
- law at an obstruction: no cusp found: 6
- law at an obstruction: not certified: 13
- law at an obstruction: phase: 1
- phase locking: known law: 2
- phase locking: reached: 1
- phase locking: stops at an obstruction: 28
- symmetric write: known law: 14
- symmetric write: reached: 4
- symmetric write: stops at an obstruction: 13
- threshold write: known law: 19
- threshold write: reached: 4
- threshold write: stops at an obstruction: 8

## Sweeps through a cusp: response to the bias and intrinsic offset

Response constant and offset at a bias of one noise floor, from the law of the completed model and from the full stochastic dynamics.

| model | second parameter | law | full dynamics |
| --- | --- | --- | --- |
| flux_qubit_tilt | beta | 1.33, -0.02 | 1.33, 0.01 |
| flux_qubit_tilt | beta | 1.31, -0.0 | 1.3, 0.0 |
| fluxonium | ell | 1.32, -0.02 | 1.34, 0.02 |
| fluxonium | ell | 1.31, 0.01 | 1.32, -0.0 |
| ion_duffing | k | 1.36, -0.8 | 1.32, -0.78 |
| ion_duffing | k | 1.34, -0.8 | 1.31, -0.78 |
| ion_phonon_laser | Delta | 0.74, 0.31 | 0.72, 0.8 |
| jba_resonator | chi | 1.36, -0.8 | 1.33, -0.79 |
| jba_resonator | chi | 1.37, -0.81 | 1.33, -0.79 |
| magnon_kerr | Dm | 1.37, -1.15 | 1.3, -1.03 |
| magnon_kerr | Dm | 1.37, -1.13 | 1.31, -1.03 |
| nuclear_bath | p | 1.35, -2.28 | 1.32, -2.3 |
| nuclear_bath | p | 1.29, -2.2 | 1.27, -2.14 |
| phase_qubit_rfsquid | ell | 1.31, -0.02 | 1.3, 0.01 |
| phase_qubit_rfsquid | ell | 1.34, -0.01 | 1.33, 0.0 |
| phase_qubit_rfsquid | ell | 1.31, 0.0 | 1.33, -0.03 |
| rydberg_bistability | Om | 1.36, 0.14 | 1.33, 0.06 |
| transmon_driven | Delta | 1.37, 0.78 | 1.3, 0.78 |
| transmon_driven | Delta | 1.34, 0.78 | 1.36, 0.8 |

## Writes into a critical plane

| model | written states | largest deviation | certified |
| --- | ---: | ---: | --- |
| cim_triangle | 6 | 0.026 | yes |
| dissipative_cat | 2 | 0.031 | no |
| kpo_triangle | 6 | 0.694 | no |
| kpo_triangle | 6 | 0.039 | no |
| kpo_triangle | 6 | 0.044 | no |
| optomech_bistable | 2 | 0.297 | no |
| rfsquid_triangle | 8 | 0.024 | no |

Calculated statements about the classical images of the hardware elements, not about the devices; the quantum dynamics at a few photons can differ (Module 15).
