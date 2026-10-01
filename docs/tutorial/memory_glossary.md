# Glossary for the memory modules

The memory modules (15–23) use the terms below with one meaning each. Symbols
follow the companion manuscript *Principles of material memory*. For each
mechanism, the section that derives it and the publications in which it was first
found are listed in [mechanisms](../mechanisms.md).

## Description of a material

| Term | Definition |
| --- | --- |
| Realization | A complete model of a material, written $I_{\mathrm{real}} = ((\Omega,\Xi);\,C,\,R,\,P;\,A)$. Memory is a property of a realization, not of a substance. |
| Carrier $\Xi$ | The space of states that can hold information: concentrations, angles, positions, conductances. Its kind is Euclidean space, a torus (angles) or an orthant (non-negative concentrations). |
| Operation $\Omega$ | The deterministic dynamics on the carrier, given as a drift $\dot q = F(q)$ or as a network of pairwise couplings. |
| Closure $C$ | What is specified or discarded to obtain closed equations: constitutive relations, admissible states, boundaries, imposed conservation laws and eliminated degrees of freedom, such as the thermal bath (noise intensity $D$). |
| Observable $R$ | The quantity that is measured. |
| Protocol $P$ | How the material is driven: the control parameter, a writing field, a temperature history. |
| Material parameters $A$ | The material or apparatus that implements the model, entered as parameter values. |
| Control parameter | The parameter an experiment varies, for example a promoter strength or a coupling strength. |
| Role of the control | *Scale*: the control multiplies the whole drift and acts as an inverse temperature. *Bias*: it adds a constant tilt to the landscape. *Shape*: it changes the form of the landscape and can create states. |

## Writing

| Term | Definition |
| --- | --- |
| Stable state | A state to which the noise-free dynamics returns after small perturbations. |
| Local relaxation rate $\kappa$ | $\kappa = -\mathrm{Re}\,\lambda$ for an eigenvalue $\lambda$ of the linearized drift. $\kappa > 0$: perturbations decay; $\kappa = 0$: a direction without restoring force; $\kappa < 0$: perturbations grow. |
| Write point | A parameter value at which a stable state loses stability ($\kappa$ crosses zero along one direction), so that a small bias can select a new state. |
| Normal form | The drift along the unstable direction at the write point, $\dot s = a_0 + a_1 s + a_2 s^2 + a_3 s^3$. |
| Supercritical pitchfork | $a_0 = a_2 = 0$, $a_3 < 0$: one state splits into two symmetric states; a weak bias selects one. This is the reference write. |
| Subcritical pitchfork | $a_2 = 0$, $a_3 > 0$: the selected state is far from the symmetric one. |
| Fold (saddle-node) | $a_2 \neq 0$: a state appears or disappears at a threshold; writing is one-sided. |
| Hopf bifurcation | A complex pair of eigenvalues crosses zero; the state starts to oscillate and no stable state is selected. |
| Cusp | A point in a plane of two parameters where $a_2$ vanishes; a sweep through it gives a symmetric write in a material whose single-parameter write is a fold. |
| Obstruction | The terms by which a normal form differs from the supercritical pitchfork: $a_2$, a bias $a_0$, or a positive $a_3$. |
| Swept write | Writing while the control parameter crosses the write point at rate $r$. For a pitchfork its accuracy is $\Phi(\pi^{1/4} h_s / (D_s^{1/2} r^{1/4}))$. |

## Retention

| Term | Definition |
| --- | --- |
| Retention time | The time after the writing field is removed at which the stored state is still identified with probability 0.9. |
| Writing time | The time under a field of bounded strength after which the new state is identified with probability 0.9. |
| Loss laws | Law 1, relaxation in a curved minimum: information decays as $e^{-2\kappa t}$. Law 2, diffusion along a direction with $\kappa = 0$: as $1/t$. Law 3, activation over a barrier: exponentially, at the Kramers rate. See [Module 4, Section 1.3](18_memory_writing_and_retention.md#13-retention). |
| Relation between retention and writing times | For a state protected by a barrier and a field of bounded strength, the ratio of retention time to writing time depends only on the work $E_w$ of the field: $\propto e^{E_w/k_BT}$. The barrier and the mobility cancel. The relation does not apply when the landscape changes between writing and retention, or when the state is written at an instability. |
| Zero mode | A direction along which the energy does not change (a continuous symmetry). A state written along it is not protected by a barrier and is lost by diffusion (Law 2). |
| Phase memory | Memory in the phase of a limit cycle, which is a zero mode created by the time-translation symmetry of the autonomous dynamics. |
| Phase response curve | The phase shift produced by a short kick, as a function of the phase at which the kick arrives. |
| Phase diffusion constant $D_\phi$ | The rate at which noise spreads the phase of a limit cycle, $\langle\delta\phi^2\rangle = 2D_\phi t$. |
| Information about a write | The mutual information between the write and the state at a later time, in nats. For the linear stage, $I(t) = \tfrac12\ln[1 + (\kappa s_0^2/D)/(e^{2\kappa t}-1)]$. |
| $W$ | The ratio $\lvert\kappa\rvert s_0^2/D$ of expansion to noise at an instability; it sets the information retained, $\tfrac12\ln(1+W)$, and the probability $\Phi(W^{1/2})$ that the sign of a binary write is kept. |

## Structure and transfer

| Term | Definition |
| --- | --- |
| Structural prediction | A statement derived from symmetries, loop signs or couplings before any simulation; each names its rule and what would falsify it. |
| Holonomy | The composition of the bond transports around a closed loop. A loop can satisfy all its bonds only if the holonomy has a fixed point. |
| Frustration | A loop whose bonds cannot all be satisfied. For $N$ equal rotor bonds with mismatch $\Phi$, the minimum excess energy per bond is $1 - \cos(\Phi/N)$. |
| Memory signature | The properties of a memory that follow from the structure of the operation: reciprocity, kind of write, retention law, oscillation, multistability, role of the control, zero modes and satisfiability of loops. |
| Transfer (attach) | Carrying a memory mechanism to another carrier, recording which entries of the signature are kept and which change. |
| Conserved field | A field whose total is fixed by the closure (Model B). Its slow modes have $\kappa(k) \propto k^2$, and a written profile is lost as $t^{-(d+2n)/z}$, with $n = 0$ for a write that adds material, $n = 1$ for one that moves it, and $z = 2$ for diffusive relaxation. |
| Non-conserved field | A field whose total is not fixed (Model A). Every mode relaxes at a rate of at least $M\kappa_0$, and a written profile is lost exponentially. |
| Signal-to-noise ratio (SNR) | For a written profile, the squared signal of the write divided by the equilibrium variance, summed over the modes of the field. |

## Derivations across fields

| Term | Definition |
| --- | --- |
| Derivation | A chain of transformations that takes a realization to a target mechanism. Each transformation acts on named components of the realization and is verified by a calculation. |
| Derivation word | The letters of a derivation: S symmetry, C continuation, W write field, R reduction, U unfolding, K canonical form, L law. For phase locking, R is the reduction to the phase and K the averaging over the drive. |
| Class of a derivation | The word together with the kind of symmetry (reflection, exchange, cyclic permutation) and the kind of the first reduction (pitchfork or fold), without the names of the variables. |
| Co-discovery by construction | Deriving one target mechanism in realizations from different fields, recording each derivation, naming the step at which a derivation stops, and testing invariants of the end points that do not depend on the field. |
| Canonical form | The reduced drift at the write point in the coordinate $x = \lvert a_3\rvert^{1/2} s$; for the symmetric write it is $\varepsilon x - x^3 + h$, with no even part. |
| Constant of the write law | $\Phi^{-1}(P)\,D_s^{1/2} r^{1/4}/h_s$ estimated from the accuracy $P$ of a simulated swept write; for the symmetric write it equals $\pi^{1/4} = 1.331$ in every realization. |
| Threshold write | A write at a fold (saddle-node), $\dot s = \mu + s^2$: the occupied state disappears when a parameter crosses a threshold, and the material switches to another state. It needs no symmetry. |
| Write field (W) | A bounded field toward another stored state, raised until the occupied state disappears; the value at which it disappears is the coercive field. |
| Delay of a switch | Under a sweep $\mu = r t$ the state crosses the position of the static fold at $\mu = \lvert a_1'\rvert r^{2/3}$, with $a_1' = -1.01879$ the first zero of the derivative of the Airy function; the switch lags behind the static threshold by an amount proportional to $r^{2/3}$. |
| Phase locking | A periodic drive of the control fixes the phase of a limit cycle: the phase difference $\psi = \theta - \omega_f t/n$ obeys the Adler equation $\dot\psi = \Delta\omega - K\sin\phi$ and relaxes to one of $n$ locked values when $\lvert\Delta\omega\rvert < K$. |
| Phase response along the drive | $Z_p(\theta) = Z(\theta)\cdot\partial F/\partial p$, where $Z$ is the gradient of the phase on the cycle (the periodic solution of the adjoint equation, normalized by $Z\cdot F = \omega_0$). Only its harmonic $n$ acts on the phase when $\omega_f \approx n\omega_0$. |
| Locking range | The detunings $\lvert\Delta\omega\rvert < K = \varepsilon\lvert Z_n\rvert/2$ over which the drive fixes the phase. Inside it the phase relaxes at $\lambda = n(K^2 - \Delta\omega^2)^{1/2}$; outside it slips at $\Omega = (\Delta\omega^2 - K^2)^{1/2}$. |
| Locking ratio | The ratio $n:1$ of drive to oscillator frequency. A symmetry that maps the cycle onto itself $1/m$ of a period later and that the drive respects leaves only the harmonics $m, 2m, \dots$, so the ratio is $m:1$ and the drive holds $m$ phases (two for a parametron). |
| Neutral cycles | A family of periodic orbits, one through every point near the cycle, as in the Lotka–Volterra model with its conserved quantity. A transverse Floquet multiplier equals 1, and a drive cannot fix an isolated phase. |

[Tutorial](index.md)
