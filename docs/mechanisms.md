# Mechanisms: definitions, derivations and original publications

For each mechanism of the [web page](https://synthetix-institute.github.io/fieldbridge/), each law of retention, each law whose constant FieldBridge certifies and each mechanism in preparation: where it is defined, where the tutorial derives it (with the functions that calculate it), and the publications in which it was first found, derived or measured. The terms of the memory modules are defined in the [glossary](tutorial/memory_glossary.md), those of the quantum chapters in [Chapter 24, §1](tutorial/24_spin_language.md#1-the-words-of-the-language). The realizations of each mechanism cite their own sources in their specification files and in the [catalog of materials](materials.md).

This page is written from [`fieldbridge/site_references.py`](../fieldbridge/site_references.py) by `python3 -B -m fieldbridge mechanisms --out docs/mechanisms.md`; a test checks that it is current.

| Mechanism | Canonical form | Defined | Derived | Original publications |
| --- | --- | --- | --- | --- |
| [Bloch rotation](#bloch-rotation) | <b>ṁ</b> = <b>Ω</b> × <b>m</b> | [Chapter 24, §1](tutorial/24_spin_language.md#1-the-words-of-the-language) | [Chapter 24, §3](tutorial/24_spin_language.md#3-detaching-the-rotation-from-the-spins-of-chapter-11), [Chapter 24, §6](tutorial/24_spin_language.md#6-a-larger-algebra-the-closure-of-the-observable) | Rabi, 1937; Bloch, 1946; Feynman, Vernon and Hellwarth, 1957 |
| [Conserved observable](#conserved-observable) | [H, R] = 0 | [Chapter 24, §6](tutorial/24_spin_language.md#6-a-larger-algebra-the-closure-of-the-observable) | [Chapter 11](tutorial/11_quantum_closure.md#how-the-program-finds-the-additional-observable) | Born and Jordan, 1925; Born, Heisenberg and Jordan, 1926 |
| [Several frequencies](#several-frequencies) | closure of R larger than three operators | [Chapter 24, §6](tutorial/24_spin_language.md#6-a-larger-algebra-the-closure-of-the-observable) | [Chapter 24, §7](tutorial/24_spin_language.md#7-obstructions) | Mori, 1965; Kitagawa and Ueda, 1993 |
| [Single stable state](#single-stable-state) | ẋ = −κx, κ &gt; 0 | [Glossary, Writing](tutorial/memory_glossary.md#writing) | [Module 4, §1.1](tutorial/18_memory_writing_and_retention.md#11-write-points-and-normal-forms), [Module 8, §2](tutorial/22_memory_time.md#2-one-expression-for-the-three-regimes-of-kappa) | Lyapunov, 1892 |
| [Symmetric write (supercritical pitchfork)](#symmetric-write-supercritical-pitchfork) | ẋ = εx − x<sup>3</sup> + h | [Glossary, Writing](tutorial/memory_glossary.md#writing) | [Module 4, §1.1–1.2](tutorial/18_memory_writing_and_retention.md#12-swept-writes), [Module 9](tutorial/23_memory_codiscovery.md#4-four-classes-of-derivation) | Poincaré, 1885; Landau, 1937; Kondepudi and Nelson, 1983; Kondepudi and Nelson, 1985; Kondepudi, Moss and McClintock, 1986 |
| [One-sided write (fold)](#one-sided-write-fold) | ẋ = μ + x<sup>2</sup> | [Glossary, Derivations across fields](tutorial/memory_glossary.md#derivations-across-fields) | [Module 10, §1](tutorial/25_memory_threshold_write.md#1-the-target), [Module 10, §3](tutorial/25_memory_threshold_write.md#3-the-delay-law) | Poincaré, 1885; Stoner and Wohlfarth, 1948; Haberman, 1979; Jung, Gray, Roy and Mandel, 1990 |
| [Write to a distant state (subcritical pitchfork)](#write-to-a-distant-state-subcritical-pitchfork) | ẋ = εx + ax<sup>3</sup> − x<sup>5</sup>, a &gt; 0 | [Glossary, Writing](tutorial/memory_glossary.md#writing) | [Module 4, §1.1](tutorial/18_memory_writing_and_retention.md#11-write-points-and-normal-forms) | Griffiths, 1970 |
| [Write by a uniform field](#write-by-a-uniform-field) | g U(q) − h·m(q) | [Glossary, Description of a material](tutorial/memory_glossary.md#description-of-a-material) | [Module 4, §2](tutorial/18_memory_writing_and_retention.md#2-worked-example-anisotropic-colloids-at-a-fluid-interface), [Module 4, §1.4](tutorial/18_memory_writing_and_retention.md#14-relation-between-retention-and-writing-times) | Stoner and Wohlfarth, 1948; Néel, 1949; Brown, 1963 |
| [Return-point memory](#return-point-memory) | f<sub>i</sub> = Σ<sub>j</sub> K<sub>ij</sub>σ<sub>j</sub> + h<sub>i</sub> + H, K<sub>ij</sub> ≥ 0 | [Module 12, §1](tutorial/27_memory_return_point.md#1-the-question) | [Module 12, §3](tutorial/27_memory_return_point.md#3-the-prediction-from-structure), [Module 12, §4](tutorial/27_memory_return_point.md#4-running-the-command) | Preisach, 1935; Harary, 1953; Barker et al., 1983; Middleton, 1992; Sethna et al., 1993; Angeli and Sontag, 2003 |
| [Return not exact](#return-not-exact) | η<sub>i</sub>J<sub>ij</sub>η<sub>j</sub> &lt; 0 on some coupling | [Module 12, §5](tutorial/27_memory_return_point.md#5-frustration-of-the-couplings-is-not-the-criterion) | [Module 12, §6](tutorial/27_memory_return_point.md#6-sufficient-not-necessary) | Deutsch, Dhar and Narayan, 2004; van Hecke, 2021 |
| [Limit cycle](#limit-cycle) | driven: φ̇ = ν − K sin φ | [Glossary, Retention](tutorial/memory_glossary.md#retention) | [Module 6](tutorial/20_memory_phase.md#1-concepts), [Module 11](tutorial/26_memory_phase_locking.md#1-the-target) | Poincaré, 1881; van der Pol, 1926; Andronov, 1929; Adler, 1946; Winfree, 1967; Lax, 1967; Guckenheimer, 1975 |
| [Neutral cycles](#neutral-cycles) | dI/dt = 0 | [Glossary, Derivations across fields](tutorial/memory_glossary.md#derivations-across-fields) | [Module 11, §6](tutorial/26_memory_phase_locking.md#6-where-the-derivation-stops) | Lotka, 1920; Volterra, 1926 |
| [Exponential loss](#exponential-loss) | SNR ∝ e<sup>−2Mκ<sub>0</sub>t</sup> | [Glossary, Structure and transfer](tutorial/memory_glossary.md#structure-and-transfer) | [Module 8, §3](tutorial/22_memory_time.md#3-fields-conservation-dimension-and-the-shape-of-the-write) | Hohenberg and Halperin, 1977; Allen and Cahn, 1979 |
| [Power-law loss](#power-law-loss) | SNR ∝ t<sup>−d/2−n</sup> | [Glossary, Structure and transfer](tutorial/memory_glossary.md#structure-and-transfer) | [Module 8, §3](tutorial/22_memory_time.md#3-fields-conservation-dimension-and-the-shape-of-the-write) | Fick, 1855; Cahn, 1961; Hohenberg and Halperin, 1977 |
| [Stochastic calculus convention](#stochastic-calculus-convention) | Itô μ − σ<sup>2</sup>/2, Stratonovich μ | [Chapter 10](tutorial/10_stochastic_construction.md) | [Chapter 10](tutorial/10_stochastic_construction.md#the-convention-changes-a-measurable-consequence) | Itô, 1944; Wong and Zakai, 1965; Stratonovich, 1966; van Kampen, 1981 |

## Mechanisms of the map

### Bloch rotation

Canonical form: <b>ṁ</b> = <b>Ω</b> × <b>m</b>. A three-vector of expectations rotates about a fixed axis, and the measured signal follows the Rabi law.

- **Defined:** [Chapter 24, §1](tutorial/24_spin_language.md#1-the-words-of-the-language), Eqs. (1) and (2): su(2) and the Rabi law.
- **Derived:** [Chapter 24, §3](tutorial/24_spin_language.md#3-detaching-the-rotation-from-the-spins-of-chapter-11), the rotation detached from two spins and derived on six carriers; [Chapter 24, §6](tutorial/24_spin_language.md#6-a-larger-algebra-the-closure-of-the-observable), the rotation reached through the closure of the observable. Code: `quantum.language.derive_bloch_rotation`, `rabi_law`, `closure_basis`.
- **Original publications:**
  - I. I. Rabi, Phys. Rev. 51, 652 (1937), [doi:10.1103/PhysRev.51.652](https://doi.org/10.1103/PhysRev.51.652): the probability of a transition of a moment in a rotating field (the Rabi law).
  - F. Bloch, Phys. Rev. 70, 460 (1946), [doi:10.1103/PhysRev.70.460](https://doi.org/10.1103/PhysRev.70.460): the precession of a nuclear magnetization in a field.
  - R. P. Feynman, F. L. Vernon and R. W. Hellwarth, J. Appl. Phys. 28, 49 (1957), [doi:10.1063/1.1722572](https://doi.org/10.1063/1.1722572): every two-level system as a vector that rotates in three dimensions.

### Conserved observable

Canonical form: [H, R] = 0. The observable commutes with the generator: its measured value does not change.

- **Defined:** [Chapter 24, §6](tutorial/24_spin_language.md#6-a-larger-algebra-the-closure-of-the-observable), a closure of the observable that contains one operator.
- **Derived:** [Chapter 11](tutorial/11_quantum_closure.md#how-the-program-finds-the-additional-observable), the closure of an observable under the Hamiltonian. Code: `quantum.language.closure_basis`, `observable_frequencies`.
- **Original publications:**
  - M. Born and P. Jordan, Z. Phys. 34, 858 (1925), [doi:10.1007/BF01328531](https://doi.org/10.1007/BF01328531): the equation of motion of an observable in matrix mechanics.
  - M. Born, W. Heisenberg and P. Jordan, Z. Phys. 35, 557 (1926), [doi:10.1007/BF01379806](https://doi.org/10.1007/BF01379806): quantities that commute with the Hamiltonian are constants of the motion.

### Several frequencies

Canonical form: closure of R larger than three operators. The commutators of the generator with the observable do not close on three operators: the observable moves with several frequencies, and the measured signal is not that of one rotation.

- **Defined:** [Chapter 24, §6](tutorial/24_spin_language.md#6-a-larger-algebra-the-closure-of-the-observable), a closure of more than three operators moves with several frequencies.
- **Derived:** [Chapter 24, §7](tutorial/24_spin_language.md#7-obstructions), the obstructions, and the term whose removal restores one rotation. Code: `quantum.language.observable_frequencies`, `_single_term_cause`.
- **Original publications:**
  - H. Mori, Prog. Theor. Phys. 34, 399 (1965), [doi:10.1143/PTP.34.399](https://doi.org/10.1143/PTP.34.399): the time dependence of an observable from its successive commutators with the generator.
  - M. Kitagawa and M. Ueda, Phys. Rev. A 47, 5138 (1993), [doi:10.1103/PhysRevA.47.5138](https://doi.org/10.1103/PhysRevA.47.5138): one-axis twisting, the term J<sub>z</sub><sup>2</sup> named for the junction with interaction.
- **Reviews and textbooks:** V. S. Viswanath and G. Müller, The Recursion Method (Springer, Berlin, 1994), [doi:10.1007/978-3-540-48651-0](https://doi.org/10.1007/978-3-540-48651-0).

### Single stable state

Canonical form: ẋ = −κx, κ &gt; 0. Every preparation relaxes to one state, and nothing of the preparation is kept.

- **Defined:** [Glossary, Writing](tutorial/memory_glossary.md#writing), stable state and local relaxation rate κ.
- **Derived:** [Module 4, §1.1](tutorial/18_memory_writing_and_retention.md#11-write-points-and-normal-forms), the state below the write point; [Module 8, §2](tutorial/22_memory_time.md#2-one-expression-for-the-three-regimes-of-kappa), the loss of a write at κ &gt; 0. Code: `memory.analysis.locate_writes`.
- **Original publications:**
  - A. M. Lyapunov, The general problem of the stability of motion (Kharkov, 1892); English translation in Int. J. Control 55, 531 (1992), [doi:10.1080/00207179208934253](https://doi.org/10.1080/00207179208934253): the stability of an equilibrium decided by the linearized motion.
- **Reviews and textbooks:** S. H. Strogatz, Nonlinear Dynamics and Chaos (Addison-Wesley, Reading, 1994).

### Symmetric write (supercritical pitchfork)

Canonical form: ẋ = εx − x<sup>3</sup> + h. Two stable states appear together at a supercritical pitchfork; a weak bias during the crossing selects the state that is written.

- **Defined:** [Glossary, Writing](tutorial/memory_glossary.md#writing), supercritical pitchfork, swept write.
- **Derived:** [Module 4, §1.1–1.2](tutorial/18_memory_writing_and_retention.md#12-swept-writes), the normal form and the law of the swept write, with its derivation; [Module 9](tutorial/23_memory_codiscovery.md#4-four-classes-of-derivation), the same canonical form and law constant in models from different fields. Code: `memory.analysis.normal_form`, `construct.swept_write_check`, `codiscovery.derive_symmetric_write`.
- **Original publications:**
  - H. Poincaré, Acta Math. 7, 259 (1885), [doi:10.1007/BF02402204](https://doi.org/10.1007/BF02402204): the bifurcation of a family of equilibria as a parameter changes.
  - L. D. Landau, Zh. Eksp. Teor. Fiz. 7, 19 (1937); English translation in Collected Papers of L. D. Landau (Pergamon, Oxford, 1965), p. 193, [doi:10.1016/B978-0-08-010586-4.50034-1](https://doi.org/10.1016/B978-0-08-010586-4.50034-1): the order parameter and the symmetric expansion of the free energy at a continuous transition.
  - D. K. Kondepudi and G. W. Nelson, Phys. Rev. Lett. 50, 1023 (1983), [doi:10.1103/PhysRevLett.50.1023](https://doi.org/10.1103/PhysRevLett.50.1023): symmetry breaking at a pitchfork far from equilibrium and its sensitivity to a small bias.
  - D. K. Kondepudi and G. W. Nelson, Nature 314, 438 (1985), [doi:10.1038/314438a0](https://doi.org/10.1038/314438a0): the probability that a weak bias selects the state when the bifurcation is crossed at a finite rate.
  - D. K. Kondepudi, F. Moss and P. V. E. McClintock, Physica D 21, 296 (1986), [doi:10.1016/0167-2789(86)90006-0](https://doi.org/10.1016/0167-2789(86)90006-0): this selection measured in a noisy electronic circuit.
- **Reviews and textbooks:** M. Golubitsky, I. Stewart and D. G. Schaeffer, Singularities and Groups in Bifurcation Theory, Vol. II (Springer, New York, 1988), [doi:10.1007/978-1-4612-4574-2](https://doi.org/10.1007/978-1-4612-4574-2); S. H. Strogatz, Nonlinear Dynamics and Chaos (Addison-Wesley, Reading, 1994).

### One-sided write (fold)

Canonical form: ẋ = μ + x<sup>2</sup>. A stored state disappears at a fold; a field past the threshold switches the state after a delay set by the Airy law.

- **Defined:** [Glossary, Derivations across fields](tutorial/memory_glossary.md#derivations-across-fields), threshold write, write field, delay of a switch.
- **Derived:** [Module 10, §1](tutorial/25_memory_threshold_write.md#1-the-target), the target and the derivation in seven fields; [Module 10, §3](tutorial/25_memory_threshold_write.md#3-the-delay-law), the delay law. Code: `memory.codiscovery.derive_threshold_write`, `fold_delay_law`, `delay_constant`.
- **Original publications:**
  - H. Poincaré, Acta Math. 7, 259 (1885), [doi:10.1007/BF02402204](https://doi.org/10.1007/BF02402204): the bifurcation of a family of equilibria as a parameter changes.
  - E. C. Stoner and E. P. Wohlfarth, Phil. Trans. R. Soc. A 240, 599 (1948), [doi:10.1098/rsta.1948.0007](https://doi.org/10.1098/rsta.1948.0007): switching of a single-domain particle when the field removes the occupied minimum (the coercive field).
  - R. Haberman, SIAM J. Appl. Math. 37, 69 (1979), [doi:10.1137/0137006](https://doi.org/10.1137/0137006): the delay of a slowly swept fold, given by the Airy function.
  - P. Jung, G. Gray, R. Roy and P. Mandel, Phys. Rev. Lett. 65, 1873 (1990), [doi:10.1103/PhysRevLett.65.1873](https://doi.org/10.1103/PhysRevLett.65.1873): the scaling of hysteresis with the rate of the sweep.
- **Reviews and textbooks:** S. H. Strogatz, Nonlinear Dynamics and Chaos (Addison-Wesley, Reading, 1994); Y. A. Kuznetsov, Elements of Applied Bifurcation Theory, 3rd edn (Springer, New York, 2004), [doi:10.1007/978-1-4757-3978-7](https://doi.org/10.1007/978-1-4757-3978-7).

### Write to a distant state (subcritical pitchfork)

Canonical form: ẋ = εx + ax<sup>3</sup> − x<sup>5</sup>, a &gt; 0. The state loses stability at a subcritical pitchfork and jumps to a distant state.

- **Defined:** [Glossary, Writing](tutorial/memory_glossary.md#writing), subcritical pitchfork.
- **Derived:** [Module 4, §1.1](tutorial/18_memory_writing_and_retention.md#11-write-points-and-normal-forms), the normal forms and the kind of each write. Code: `memory.analysis.normal_form`.
- **Original publications:**
  - R. B. Griffiths, Phys. Rev. Lett. 24, 715 (1970), [doi:10.1103/PhysRevLett.24.715](https://doi.org/10.1103/PhysRevLett.24.715): the tricritical point, at which the quartic coefficient changes sign and a continuous transition becomes discontinuous.
- **Reviews and textbooks:** S. H. Strogatz, Nonlinear Dynamics and Chaos (Addison-Wesley, Reading, 1994).

### Write by a uniform field

Canonical form: g U(q) − h·m(q). The control only rescales the energy; a uniform field writes a state, and the barriers between states keep it.

- **Defined:** [Glossary, Description of a material](tutorial/memory_glossary.md#description-of-a-material), role of the control: scale, bias or shape.
- **Derived:** [Module 4, §2](tutorial/18_memory_writing_and_retention.md#2-worked-example-anisotropic-colloids-at-a-fluid-interface), anisotropic colloids at a fluid interface, written by a field; [Module 4, §1.4](tutorial/18_memory_writing_and_retention.md#14-relation-between-retention-and-writing-times), the relation between retention and writing times. Code: `memory.analysis.hold_time`, `write_test`.
- **Original publications:**
  - E. C. Stoner and E. P. Wohlfarth, Phil. Trans. R. Soc. A 240, 599 (1948), [doi:10.1098/rsta.1948.0007](https://doi.org/10.1098/rsta.1948.0007): hysteresis of single-domain particles written by a uniform field.
  - L. Néel, Ann. Géophys. 5, 99 (1949): thermally activated loss of the magnetization of fine grains.
  - W. F. Brown, Phys. Rev. 130, 1677 (1963), [doi:10.1103/PhysRev.130.1677](https://doi.org/10.1103/PhysRev.130.1677): the thermal fluctuations of a single-domain particle.

### Return-point memory

Canonical form: f<sub>i</sub> = Σ<sub>j</sub> K<sub>ij</sub>σ<sub>j</sub> + h<sub>i</sub> + H, K<sub>ij</sub> ≥ 0. Elements switch at thresholds of a slow drive. With cooperative couplings and a drive that pushes every element the same way, the state at a turning point is recovered exactly after any excursion inside it.

- **Defined:** [Module 12, §1](tutorial/27_memory_return_point.md#1-the-question), return-point memory.
- **Derived:** [Module 12, §3](tutorial/27_memory_return_point.md#3-the-prediction-from-structure), the drive counted as an element: the relabeling and no passing; [Module 12, §4](tutorial/27_memory_return_point.md#4-running-the-command), six realizations from three fields. Code: `memory.hysterons.predict`, `relabel`, `check`.
- **Original publications:**
  - F. Preisach, Z. Phys. 94, 277 (1935), [doi:10.1007/BF01349418](https://doi.org/10.1007/BF01349418): hysteresis as a population of elementary loops (hysterons).
  - F. Harary, Michigan Math. J. 2, 143 (1953), [doi:10.1307/mmj/1028989917](https://doi.org/10.1307/mmj/1028989917): the balance of a signed graph, here with the drive as one of its vertices.
  - J. A. Barker, D. E. Schreiber, B. G. Huth and D. H. Everett, Proc. R. Soc. Lond. A 386, 251 (1983), [doi:10.1098/rspa.1983.0035](https://doi.org/10.1098/rspa.1983.0035): return-point memory and minor loops in magnets.
  - A. A. Middleton, Phys. Rev. Lett. 68, 670 (1992), [doi:10.1103/PhysRevLett.68.670](https://doi.org/10.1103/PhysRevLett.68.670): the order of states under a monotonic drive (no passing).
  - J. P. Sethna, K. Dahmen, S. Kartha, J. A. Krumhansl, B. W. Roberts and J. D. Shore, Phys. Rev. Lett. 70, 3347 (1993), [doi:10.1103/PhysRevLett.70.3347](https://doi.org/10.1103/PhysRevLett.70.3347): return-point memory in a disordered model of first-order transitions.
  - D. Angeli and E. D. Sontag, IEEE Trans. Automat. Control 48, 1684 (2003), [doi:10.1109/TAC.2003.817920](https://doi.org/10.1109/TAC.2003.817920): the same sign condition for a monotone system with an input.
- **Reviews and textbooks:** N. C. Keim, J. D. Paulsen, Z. Zeravcic, S. Sastry and S. R. Nagel, Rev. Mod. Phys. 91, 035002 (2019), [doi:10.1103/RevModPhys.91.035002](https://doi.org/10.1103/RevModPhys.91.035002).

### Return not exact

Canonical form: η<sub>i</sub>J<sub>ij</sub>η<sub>j</sub> &lt; 0 on some coupling. A coupling closes a frustrated loop through the drive: after an excursion the state at a turning point can differ.

- **Defined:** [Module 12, §5](tutorial/27_memory_return_point.md#5-frustration-of-the-couplings-is-not-the-criterion), the frustration of the couplings is not the criterion.
- **Derived:** [Module 12, §6](tutorial/27_memory_return_point.md#6-sufficient-not-necessary), a frustrated loop through the drive allows a failure and does not force one. Code: `memory.hysterons.predict`, `check`.
- **Original publications:**
  - J. M. Deutsch, A. Dhar and O. Narayan, Phys. Rev. Lett. 92, 227203 (2004), [doi:10.1103/PhysRevLett.92.227203](https://doi.org/10.1103/PhysRevLett.92.227203): random antiferromagnetic chains return exactly, although no passing fails.
  - M. van Hecke, Phys. Rev. E 104, 054608 (2021), [doi:10.1103/PhysRevE.104.054608](https://doi.org/10.1103/PhysRevE.104.054608): most transition graphs of three interacting hysterons violate return-point memory.

### Limit cycle

Canonical form: driven: φ̇ = ν − K sin φ. The preparations settle on a limit cycle whose phase is neutral; a periodic drive locks the phase inside the Adler range.

- **Defined:** [Glossary, Retention](tutorial/memory_glossary.md#retention), phase memory, phase response curve, phase diffusion.
- **Derived:** [Module 6](tutorial/20_memory_phase.md#1-concepts), memory in the phase of an oscillator; [Module 11](tutorial/26_memory_phase_locking.md#1-the-target), phase locking in eight oscillators and the locking ratio. Code: `memory.phase.response_curve`, `hold`, `lock`, `memory.phase_locking.derive`, `adler_law`.
- **Original publications:**
  - H. Poincaré, J. Math. Pures Appl. (3) 7, 375 (1881) and 8, 251 (1882): the limit cycle.
  - B. van der Pol, Phil. Mag. 2, 978 (1926), [doi:10.1080/14786442608564127](https://doi.org/10.1080/14786442608564127): relaxation oscillations of a triode circuit.
  - A. A. Andronov, C. R. Acad. Sci. Paris 189, 559 (1929): self-sustained oscillations as Poincaré's limit cycles.
  - R. Adler, Proc. IRE 34, 351 (1946), [doi:10.1109/JRPROC.1946.229930](https://doi.org/10.1109/JRPROC.1946.229930): the locking of an oscillator to an injected signal (the Adler equation).
  - A. T. Winfree, J. Theor. Biol. 16, 15 (1967), [doi:10.1016/0022-5193(67)90051-3](https://doi.org/10.1016/0022-5193(67)90051-3): the phase response and the synchronization of populations of oscillators.
  - M. Lax, Phys. Rev. 160, 290 (1967), [doi:10.1103/PhysRev.160.290](https://doi.org/10.1103/PhysRev.160.290): the diffusion of the phase of a self-sustained oscillator under noise.
  - J. Guckenheimer, J. Math. Biol. 1, 259 (1975), [doi:10.1007/BF01273747](https://doi.org/10.1007/BF01273747): the isochrons, which give the phase of a state near a limit cycle.
- **Reviews and textbooks:** Y. Kuramoto, Chemical Oscillations, Waves, and Turbulence (Springer, Berlin, 1984), [doi:10.1007/978-3-642-69689-3](https://doi.org/10.1007/978-3-642-69689-3); A. Pikovsky, M. Rosenblum and J. Kurths, Synchronization (Cambridge University Press, 2001), [doi:10.1017/CBO9780511755743](https://doi.org/10.1017/CBO9780511755743).

### Neutral cycles

Canonical form: dI/dt = 0. A conserved quantity fills the plane with closed orbits, none of which attracts its neighbours.

- **Defined:** [Glossary, Derivations across fields](tutorial/memory_glossary.md#derivations-across-fields), neutral cycles.
- **Derived:** [Module 11, §6](tutorial/26_memory_phase_locking.md#6-where-the-derivation-stops), why a drive cannot fix an isolated phase. Code: `memory.phase_locking.floquet`.
- **Original publications:**
  - A. J. Lotka, Proc. Natl. Acad. Sci. USA 6, 410 (1920), [doi:10.1073/pnas.6.7.410](https://doi.org/10.1073/pnas.6.7.410): undamped oscillations in a model of interacting species.
  - V. Volterra, Nature 118, 558 (1926), [doi:10.1038/118558a0](https://doi.org/10.1038/118558a0): the predator-prey cycles and their conserved quantity.

### Exponential loss

Canonical form: SNR ∝ e<sup>−2Mκ<sub>0</sub>t</sup>. Every mode relaxes at a finite rate: the trace of a write is lost exponentially.

- **Defined:** [Glossary, Structure and transfer](tutorial/memory_glossary.md#structure-and-transfer), non-conserved field.
- **Derived:** [Module 8, §3](tutorial/22_memory_time.md#3-fields-conservation-dimension-and-the-shape-of-the-write), conservation, dimension and the shape of the write. Code: `memory.fields.predicted_law`, `hartree_mass`.
- **Original publications:**
  - P. C. Hohenberg and B. I. Halperin, Rev. Mod. Phys. 49, 435 (1977), [doi:10.1103/RevModPhys.49.435](https://doi.org/10.1103/RevModPhys.49.435): relaxational dynamics without (Model A) and with (Model B) conservation.
  - S. M. Allen and J. W. Cahn, Acta Metall. 27, 1085 (1979), [doi:10.1016/0001-6160(79)90196-2](https://doi.org/10.1016/0001-6160(79)90196-2): the relaxation of a non-conserved order parameter.

### Power-law loss

Canonical form: SNR ∝ t<sup>−d/2−n</sup>. A conserved density relaxes slowly at long wavelengths: the trace of a write decays as a power of time.

- **Defined:** [Glossary, Structure and transfer](tutorial/memory_glossary.md#structure-and-transfer), conserved field.
- **Derived:** [Module 8, §3](tutorial/22_memory_time.md#3-fields-conservation-dimension-and-the-shape-of-the-write), conservation, dimension and the shape of the write. Code: `memory.fields.predicted_law`, `spectral_snr`, `exponent_check`.
- **Original publications:**
  - A. Fick, Ann. Phys. 170, 59 (1855), [doi:10.1002/andp.18551700105](https://doi.org/10.1002/andp.18551700105): the diffusion of a conserved density.
  - J. W. Cahn, Acta Metall. 9, 795 (1961), [doi:10.1016/0001-6160(61)90182-1](https://doi.org/10.1016/0001-6160(61)90182-1): the relaxation of a conserved composition.
  - P. C. Hohenberg and B. I. Halperin, Rev. Mod. Phys. 49, 435 (1977), [doi:10.1103/RevModPhys.49.435](https://doi.org/10.1103/RevModPhys.49.435): relaxational dynamics without (Model A) and with (Model B) conservation.

### Stochastic calculus convention

Canonical form: Itô μ − σ<sup>2</sup>/2, Stratonovich μ. The reading of the noise term is part of the closure: it changes the measured growth rate of log X by σ<sup>2</sup>/2.

- **Defined:** [Chapter 10](tutorial/10_stochastic_construction.md), why a change of stochastic coordinate adds a drift.
- **Derived:** [Chapter 10](tutorial/10_stochastic_construction.md#the-convention-changes-a-measurable-consequence), the convention changes a measurable consequence. Code: `verification.verify_construction`.
- **Original publications:**
  - K. Itô, Proc. Imp. Acad. Tokyo 20, 519 (1944), [doi:10.3792/pia/1195572786](https://doi.org/10.3792/pia/1195572786): the Itô integral.
  - E. Wong and M. Zakai, Ann. Math. Stat. 36, 1560 (1965), [doi:10.1214/aoms/1177699916](https://doi.org/10.1214/aoms/1177699916): smooth noise converges to the Stratonovich reading.
  - R. L. Stratonovich, SIAM J. Control 4, 362 (1966), [doi:10.1137/0304028](https://doi.org/10.1137/0304028): the Stratonovich integral.
  - N. G. van Kampen, J. Stat. Phys. 24, 175 (1981), [doi:10.1007/BF01007642](https://doi.org/10.1007/BF01007642): the convention as part of the physical model.
- **Reviews and textbooks:** H. Risken, The Fokker-Planck Equation (Springer, Berlin, 1984), [doi:10.1007/978-3-642-96807-5](https://doi.org/10.1007/978-3-642-96807-5).

## Laws of retention

A stored state is lost by one of three laws, set by the form of the landscape at the state ([Module 4, §1.3](tutorial/18_memory_writing_and_retention.md#13-retention)).

### Law 1

Relaxation in a curved minimum (κ &gt; 0): the information about the write decreases as e<sup>−2κt</sup>.

- **Defined:** [Module 4, §1.3](tutorial/18_memory_writing_and_retention.md#13-retention), the three laws of retention.
- **Derived:** [Module 8, §2](tutorial/22_memory_time.md#2-one-expression-for-the-three-regimes-of-kappa), Laws 1 and 2 from one expression. Code: `memory.regimes.gaussian_information`.
- **Original publications:**
  - G. E. Uhlenbeck and L. S. Ornstein, Phys. Rev. 36, 823 (1930), [doi:10.1103/PhysRev.36.823](https://doi.org/10.1103/PhysRev.36.823): the relaxation of a Brownian particle bound by a linear force.
  - C. E. Shannon, Proc. IRE 37, 10 (1949), [doi:10.1109/JRPROC.1949.232969](https://doi.org/10.1109/JRPROC.1949.232969): the information carried through Gaussian noise, ½ ln(1 + SNR).

### Law 2

Diffusion along a direction with κ = 0, such as the phase of a limit cycle: the information decreases as 1/t.

- **Defined:** [Module 4, §1.3](tutorial/18_memory_writing_and_retention.md#13-retention), the three laws of retention.
- **Derived:** [Module 8, §2](tutorial/22_memory_time.md#2-one-expression-for-the-three-regimes-of-kappa), Laws 1 and 2 from one expression; [Module 6](tutorial/20_memory_phase.md#1-concepts), the phase of a limit cycle. Code: `memory.regimes.gaussian_information`, `memory.phase.hold`.
- **Original publications:**
  - A. Einstein, Ann. Phys. 322, 549 (1905), [doi:10.1002/andp.19053220806](https://doi.org/10.1002/andp.19053220806): free diffusion, ⟨δx<sup>2</sup>⟩ = 2Dt.
  - J. Goldstone, Nuovo Cimento 19, 154 (1961), [doi:10.1007/BF02812722](https://doi.org/10.1007/BF02812722): a zero mode from a broken continuous symmetry.
  - M. Lax, Phys. Rev. 160, 290 (1967), [doi:10.1103/PhysRev.160.290](https://doi.org/10.1103/PhysRev.160.290): the diffusion of the phase of a self-sustained oscillator under noise.

### Law 3

Activation over a barrier ΔV between stored states, at the Kramers rate, proportional to e<sup>−ΔV/k<sub>B</sub>T</sup>.

- **Defined:** [Module 4, §1.3](tutorial/18_memory_writing_and_retention.md#13-retention), the three laws of retention.
- **Derived:** [Module 4, §1.3](tutorial/18_memory_writing_and_retention.md#13-retention), the law named on the memory card; the retention time is calculated. Code: `memory.analysis.hold_time`.
- **Original publications:**
  - S. Arrhenius, Z. Phys. Chem. 4, 226 (1889), [doi:10.1515/zpch-1889-0416](https://doi.org/10.1515/zpch-1889-0416): the rate of a reaction proportional to e<sup>−E/k<sub>B</sub>T</sup>.
  - H. A. Kramers, Physica 7, 284 (1940), [doi:10.1016/S0031-8914(40)90098-2](https://doi.org/10.1016/S0031-8914(40)90098-2): the rate of escape over a barrier by Brownian motion.
  - L. Néel, Ann. Géophys. 5, 99 (1949): the activated loss of the magnetization of fine grains.
- **Reviews and textbooks:** P. Hänggi, P. Talkner and M. Borkovec, Rev. Mod. Phys. 62, 251 (1990), [doi:10.1103/RevModPhys.62.251](https://doi.org/10.1103/RevModPhys.62.251); H. Risken, The Fokker-Planck Equation (Springer, Berlin, 1984), [doi:10.1007/978-3-642-96807-5](https://doi.org/10.1007/978-3-642-96807-5).

## Certified laws

The laws whose constant FieldBridge calculates in every realization that reaches the target, shown on the page under *One mechanism in different fields*.

### Rabi law of the Bloch rotation

Law: f(t) = cos<sup>2</sup>θ + sin<sup>2</sup>θ cos(|Ω|t).

- **Defined:** [Chapter 24, §1](tutorial/24_spin_language.md#1-the-words-of-the-language), Eq. (2).
- **Derived:** [Chapter 24, §5](tutorial/24_spin_language.md#5-co-discovery-one-rotation-in-five-fields), one rotation on carriers from five fields. Code: `quantum.language.rabi_law`.
- **Original publications:**
  - I. I. Rabi, Phys. Rev. 51, 652 (1937), [doi:10.1103/PhysRev.51.652](https://doi.org/10.1103/PhysRev.51.652): the probability of a transition of a moment in a rotating field (the Rabi law).
  - F. Bloch, Phys. Rev. 70, 460 (1946), [doi:10.1103/PhysRev.70.460](https://doi.org/10.1103/PhysRev.70.460): the precession of a nuclear magnetization in a field.

### Swept-write law, constant π<sup>1/4</sup>

Law: P = Φ(π<sup>1/4</sup> h<sub>s</sub> / (D<sub>s</sub><sup>1/2</sup> r<sup>1/4</sup>)).

- **Defined:** [Glossary, Derivations across fields](tutorial/memory_glossary.md#derivations-across-fields), constant of the write law.
- **Derived:** [Module 4, §1.2](tutorial/18_memory_writing_and_retention.md#12-swept-writes), derivation in the linear stage of the sweep; [Module 9, §6](tutorial/23_memory_codiscovery.md#6-two-invariants-of-the-end-point), the constant in every realization. Code: `memory.construct.swept_write_check`, `codiscovery.law_constant`.
- **Original publications:**
  - D. K. Kondepudi and G. W. Nelson, Nature 314, 438 (1985), [doi:10.1038/314438a0](https://doi.org/10.1038/314438a0): the probability that a weak bias selects the state when the bifurcation is crossed at a finite rate.
  - D. K. Kondepudi, F. Moss and P. V. E. McClintock, Physica D 21, 296 (1986), [doi:10.1016/0167-2789(86)90006-0](https://doi.org/10.1016/0167-2789(86)90006-0): this selection measured in a noisy electronic circuit.

### Delay of a switch, constant |a<sub>1</sub>′| = 1.0188

Law: μ<sub>switch</sub> = |a<sub>1</sub>′| r<sup>2/3</sup>.

- **Defined:** [Glossary, Derivations across fields](tutorial/memory_glossary.md#derivations-across-fields), delay of a switch.
- **Derived:** [Module 10, §3](tutorial/25_memory_threshold_write.md#3-the-delay-law), the delay law. Code: `memory.codiscovery.fold_delay_law`, `delay_constant`.
- **Original publications:**
  - R. Haberman, SIAM J. Appl. Math. 37, 69 (1979), [doi:10.1137/0137006](https://doi.org/10.1137/0137006): the delay of a slowly swept fold, given by the Airy function.
  - P. Jung, G. Gray, R. Roy and P. Mandel, Phys. Rev. Lett. 65, 1873 (1990), [doi:10.1103/PhysRevLett.65.1873](https://doi.org/10.1103/PhysRevLett.65.1873): the scaling of hysteresis with the rate of the sweep.

### Locking range of the Adler equation

Law: ψ̇ = Δω − K sin ψ: locked for |Δω| &lt; K.

- **Defined:** [Glossary, Derivations across fields](tutorial/memory_glossary.md#derivations-across-fields), phase locking, locking range, locking ratio.
- **Derived:** [Module 11, §1](tutorial/26_memory_phase_locking.md#1-the-target), the target and its law; [Module 11, §4](tutorial/26_memory_phase_locking.md#4-a-symmetry-sets-the-ratio), a symmetry sets the ratio. Code: `memory.phase_locking.adler_law`, `locked_rate`, `slip_frequency`.
- **Original publications:**
  - R. Adler, Proc. IRE 34, 351 (1946), [doi:10.1109/JRPROC.1946.229930](https://doi.org/10.1109/JRPROC.1946.229930): the locking of an oscillator to an injected signal (the Adler equation).
- **Reviews and textbooks:** A. Pikovsky, M. Rosenblum and J. Kurths, Synchronization (Cambridge University Press, 2001), [doi:10.1017/CBO9780511755743](https://doi.org/10.1017/CBO9780511755743).

## Mechanisms in preparation

Mechanisms that are studied in the tutorial or in the literature on memory and are not yet derivation targets with a certified law ([roadmap](ROADMAP.md)).

### Frustrated loops

Law: excess energy per bond 1 − cos(Φ/N) for N equal rotor bonds with mismatch Φ.

- **Defined:** [Glossary, Structure and transfer](tutorial/memory_glossary.md#structure-and-transfer), holonomy, frustration.
- **Derived:** [Module 3, §4](tutorial/17_memory_predictions.md#4-tests-on-many-loops-and-networks), tests on many loops and networks. Code: `memory.predict`, `memory.networks`, `memory.compose`.
- **Original publications:**
  - F. Harary, Michigan Math. J. 2, 143 (1953), [doi:10.1307/mmj/1028989917](https://doi.org/10.1307/mmj/1028989917): the balance of a signed graph.
  - G. Toulouse, Commun. Phys. 2, 115 (1977): the frustration of a plaquette in spin glasses.
  - J. Villain, J. Phys. C 10, 1717 (1977), [doi:10.1088/0022-3719/10/10/014](https://doi.org/10.1088/0022-3719/10/10/014): frustration without disorder.
  - R. Thomas, in Numerical Methods in the Study of Critical Phenomena, Springer Series in Synergetics 9 (Springer, Berlin, 1981), p. 180, [doi:10.1007/978-3-642-81703-8_24](https://doi.org/10.1007/978-3-642-81703-8_24): positive loops for several steady states, negative loops for oscillation.
  - S. Teitel and C. Jayaprakash, Phys. Rev. B 27, 598 (1983), [doi:10.1103/PhysRevB.27.598](https://doi.org/10.1103/PhysRevB.27.598): rotors with a mismatch around each plaquette.
  - C. Soulé, ComPlexUs 1, 123 (2003), [doi:10.1159/000076100](https://doi.org/10.1159/000076100): the proof of the condition on positive loops.

### Retention against rewriting

Law: the ratio of the retention time to the writing time depends only on the work E<sub>w</sub> of the write: it grows as ln E<sub>w</sub> in a curved minimum, linearly along a zero mode and as e<sup>E<sub>w</sub>/k<sub>B</sub>T</sup> behind a barrier.

- **Defined:** [Glossary, Retention](tutorial/memory_glossary.md#retention), relation between retention and writing times.
- **Derived:** [Module 4, §1.4](tutorial/18_memory_writing_and_retention.md#14-relation-between-retention-and-writing-times), the relation and the three kinds of protocol not subject to it. Code: `memory.analysis.hold_time`, `write_test`.
- **Original publications:**
  - H. A. Kramers, Physica 7, 284 (1940), [doi:10.1016/S0031-8914(40)90098-2](https://doi.org/10.1016/S0031-8914(40)90098-2): the rate of escape over a barrier by Brownian motion.
  - G. I. Bell, Science 200, 618 (1978), [doi:10.1126/science.347575](https://doi.org/10.1126/science.347575): a force lowers the barrier and accelerates its crossing.
  - S. Fusi, P. J. Drew and L. F. Abbott, Neuron 45, 599 (2005), [doi:10.1016/j.neuron.2005.02.001](https://doi.org/10.1016/j.neuron.2005.02.001): the conflict between the rate of learning and the time of retention in synapses.
  - M. H. Kryder, E. C. Gage, T. W. McDaniel, W. A. Challener, R. E. Rottmayer, G. Ju, Y.-T. Hsia and M. F. Erden, Proc. IEEE 96, 1810 (2008), [doi:10.1109/JPROC.2008.2004315](https://doi.org/10.1109/JPROC.2008.2004315): heating during the write: a landscape changed between writing and retention.
  - M. K. Benna and S. Fusi, Nat. Neurosci. 19, 1697 (2016), [doi:10.1038/nn.4401](https://doi.org/10.1038/nn.4401): memory consolidation in synapses.

### Onset of oscillation (Hopf)

Law: ż = (μ + iω)z − |z|<sup>2</sup>z: an amplitude proportional to √μ; under a slow sweep of μ the onset is delayed (Neishtadt, 1987).

- **Defined:** [Glossary, Writing](tutorial/memory_glossary.md#writing), Hopf bifurcation.
- **Derived:** [Module 4, §1.1](tutorial/18_memory_writing_and_retention.md#11-write-points-and-normal-forms), the kinds of write points; a Hopf bifurcation stops a write. Code: `memory.analysis.locate_writes`.
- **Original publications:**
  - E. Hopf, Ber. Math.-Phys. Kl. Sächs. Akad. Wiss. Leipzig 94, 1 (1942); English translation in J. E. Marsden and M. McCracken, The Hopf Bifurcation and Its Applications (Springer, New York, 1976), [doi:10.1007/978-1-4612-6374-6](https://doi.org/10.1007/978-1-4612-6374-6): a periodic solution that branches from an equilibrium when a complex pair of eigenvalues crosses the imaginary axis.
  - A. I. Neishtadt, Differ. Uravn. 23, 2060 (1987); English translation in Differential Equations 23, 1385 (1987): the delay of the loss of stability under a slow sweep.
  - S. M. Baer, T. Erneux and J. Rinzel, SIAM J. Appl. Math. 49, 55 (1989), [doi:10.1137/0149003](https://doi.org/10.1137/0149003): the delay of a slow passage through a Hopf bifurcation and its dependence on the initial state.
- **Reviews and textbooks:** Y. A. Kuznetsov, Elements of Applied Bifurcation Theory, 3rd edn (Springer, New York, 2004), [doi:10.1007/978-1-4757-3978-7](https://doi.org/10.1007/978-1-4757-3978-7).

### Synchronization of a population (Kuramoto)

Law: oscillators with a symmetric, unimodal density g of natural frequencies synchronize above the coupling K<sub>c</sub> = 2/(π g(0)) (Kuramoto, 1984).

- **Defined:** [Roadmap](ROADMAP.md#mechanisms-in-preparation), mechanisms in preparation.
- **Derived:** not yet in FieldBridge.
- **Original publications:**
  - A. T. Winfree, J. Theor. Biol. 16, 15 (1967), [doi:10.1016/0022-5193(67)90051-3](https://doi.org/10.1016/0022-5193(67)90051-3): the phase response and the synchronization of populations of oscillators.
  - Y. Kuramoto, in International Symposium on Mathematical Problems in Theoretical Physics, Lecture Notes in Physics 39 (Springer, Berlin, 1975), p. 420, [doi:10.1007/BFb0013365](https://doi.org/10.1007/BFb0013365): the synchronization threshold of a population of coupled phase oscillators.
- **Reviews and textbooks:** Y. Kuramoto, Chemical Oscillations, Waves, and Turbulence (Springer, Berlin, 1984), [doi:10.1007/978-3-642-69689-3](https://doi.org/10.1007/978-3-642-69689-3); S. H. Strogatz, Physica D 143, 1 (2000), [doi:10.1016/S0167-2789(00)00094-4](https://doi.org/10.1016/S0167-2789(00)00094-4).
