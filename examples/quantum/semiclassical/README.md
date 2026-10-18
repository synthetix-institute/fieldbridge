# Quantum-computing elements as semiclassical drift equations

Thirty-one elements of quantum-computing hardware, each written as the drift of its mean field (the quadratures of
the mode, the phase of a junction, the Bloch components of a spin ensemble) with the vacuum noise of its loss channel
as the noise term. Every specification follows the schema `fieldbridge-memory/1` of the memory module, kind
`equations`: it names the control that an experiment sweeps, the parameters that a design can change (`material`),
the closure that produced the equations from the Lindblad description (`closure`), and the observable that is read
out. The sources are the published models; the truncation to a mean field with vacuum noise is the assumption that
each file states.

The set comes from a census of the hardware reviews (superconducting circuits, trapped ions, cavity optomechanics,
Rydberg ensembles, cavity magnonics, spin ensembles) and of the quantum-computing literature: every element with a
classical write, that is a mean-field drift that has a threshold or a bistability along an experimental control,
is here. Elements whose classical image is linear or relaxes to a single state are not, since the write rules say
nothing about them: couplers kept monostable by design, travelling-wave amplifiers, bosonic-code cavities whose
codewords have zero mean field, spin and ion qubits whose operations are rotations and whose pumping has one fixed
point, and linear optics. The exact Lindblad evolution of
the Kerr parametric oscillator is in `../open/`; the equations here are its classical image and the images of its
relatives. The Kerr terms are written as the literature writes the mean field, with the rotation rate 2K r²; the
truncated-Wigner image of the Hamiltonian term −K a†a†aa rotates at 2K(r² − 1), a shift that is negligible for
K ≪ κ and comparable to the other rates for the Kerr-cat elements (K ≳ κ), where a classical image is in any case
a poor guide to the write: below about ten stored photons the exact evolution keeps far more of a bias than either
classical image does.

## Running the rules

The write rules and the completion rung over the whole folder:

```bash
python3 -m fieldbridge memory predictions examples/quantum/semiclassical --out-dir out/qc --workers 2
```

For each element the report gives the three write targets (symmetric write, threshold write, phase locking) as
`memory codiscover` derives them, with the law constant where the law was checked, and, where the symmetric write is
obstructed, what the completion rung keeps and the law of the completed model with its certificate. One element at a
time:

```bash
python3 -m fieldbridge memory codiscover examples/quantum/semiclassical/kpo_triangle.json --target symmetric-write
python3 -m fieldbridge memory complete examples/quantum/semiclassical/kpo_triangle.json
```

The report of the whole folder is kept under [docs/predictions](../../../docs/predictions/quantum_semiclassical.md).

## The elements

| File | Element | Variables | Control | Source |
| --- | --- | --- | --- | --- |
| `bright_state_resonator.json` | driven Jaynes-Cummings resonator at high power (bright state, punch-out) | x, y, sx, sy, sz | xi | Reed, M. D., DiCarlo, L., Johnson, B. R., Sun, L., Schuster, D. I., Frunzio, L., Schoelkopf, R. J., 'High-f... |
| `cim_dopo.json` | degenerate optical parametric oscillator with pump depletion (coherent Ising machine node) | s, q | Ep | Drummond, McNeil & Walls, Opt. Acta 28, 211 (1981); Marandi et al., Nat. Photon. 8, 937 (2014) |
| `cim_pair.json` | coherent Ising machine, two coupled DOPOs (ferromagnetic) | x1, x2 | a | Yamamoto, Y., Aihara, K., Leleu, T., Kawarabayashi, K., Kako, S., Fejer, M., Inoue, K., Takesue, H., 'Coher... |
| `cim_triangle.json` | coherent Ising machine, three DOPOs on an antiferromagnetic triangle | x1, x2, x3 | a | Yamamoto, Y., Aihara, K., Leleu, T., Kawarabayashi, K., Kako, S., Fejer, M., Inoue, K., Takesue, H., 'Coher... |
| `dissipative_cat.json` | dissipative cat qubit (two-photon loss) | x, y | A | Mirrahimi et al., New J. Phys. 16, 045014 (2014); Lescanne et al., Nat. Phys. 16, 509 (2020) |
| `flux_qubit_tilt.json` | persistent-current flux qubit, double well tilted by the loop flux | th, v | eps | Krantz, P., Kjaergaard, M., Yan, F., Orlando, T. P., Gustavsson, S., Oliver, W. D., 'A quantum engineer's g... |
| `fluxonium.json` | fluxonium (superinductance and junction), classical phase particle | ph, v | phe | Krantz, P., Kjaergaard, M., Yan, F., Orlando, T. P., Gustavsson, S., Oliver, W. D., 'A quantum engineer's g... |
| `ion_duffing.json` | driven single-ion Duffing oscillator (octupole anharmonicity) | x, y | Delta | Akerman, N., Kotler, S., Glickman, Y., Dallal, Y., Keselman, A., Ozeri, R., 'A single-ion nonlinear mechani... |
| `ion_phonon_laser.json` | trapped-ion phonon laser (Van der Pol image) with injection | x, y | F | Knuenz, S., 'A trapped-ion phonon laser and the detection of ultra-weak forces', PhD thesis, LMU Munich (20... |
| `jba_resonator.json` | driven Kerr resonator of a Josephson bifurcation amplifier | x, y | F | Siddiqi et al., PRL 93, 207002 (2004); Vijay, Devoret & Siddiqi, RSI 80, 111101 (2009) |
| `jpc_nondegenerate.json` | non-degenerate parametric oscillator (Josephson parametric converter above threshold) | ax, ay, bx, by, cx, cy | Ein | Bergeal, N., Schackert, F., Metcalfe, M., Vijay, R., Manucharyan, V. E., Frunzio, L., Prober, D. E., Schoel... |
| `kerr_cat_kpo.json` | Kerr parametric oscillator (Kerr-cat qubit) | x, y | p | Goto, Sci. Rep. 6, 21686 (2016); Puri et al., npj Quantum Inf. 3, 18 (2017); Grimm et al., Nature 584, 205... |
| `kpo_pair.json` | two coupled Kerr parametric oscillators (two cat qubits or a two-spin annealer) | x1, y1, x2, y2 | p | Goto, Sci. Rep. 6, 21686 (2016); Puri et al., Nat. Commun. 8, 15785 (2017) |
| `kpo_triangle.json` | three Kerr parametric oscillators on an antiferromagnetic triangle | x1, y1, x2, y2, x3, y3 | p | Goto, Sci. Rep. 6, 21686 (2016); frustrated Ising triangle |
| `magnon_kerr.json` | cavity magnon polariton with magnon Kerr nonlinearity (YIG sphere in a microwave cavity) | ax, ay, bx, by | eps | Wang, Zhang, Shen, Zhang, Tang, You, Phys. Rev. Lett. 120, 057202 (2018) |
| `nuclear_bath.json` | dynamically polarised nuclear-spin bath (Overhauser feedback, rate equation) | x | d | Vink, I. T., Nowack, K. C., Koppens, F. H. L., Danon, J., Nazarov, Yu. V., Vandersypen, L. M. K., 'Locking... |
| `optomech_bistable.json` | cavity optomechanical element, static bistability | ar, ai, X, V | Ein | Aspelmeyer, M., Kippenberg, T. J., Marquardt, F., 'Cavity optomechanics', Rev. Mod. Phys. 86, 1391 (2014)... |
| `optomech_selfosc.json` | cavity optomechanical element, self-sustained oscillation | ar, ai, X, V | Ein | Aspelmeyer, M., Kippenberg, T. J., Marquardt, F., 'Cavity optomechanics', Rev. Mod. Phys. 86, 1391 (2014)... |
| `parametric_amplifier.json` | parametric amplifier below threshold | x, y | p | Castellanos-Beltran et al., Nat. Phys. 4, 929 (2008) |
| `phase_qubit_rfsquid.json` | flux-biased rf-SQUID phase qubit (and Josephson photomultiplier) | ph, v | phe | Wendin, G., 'Quantum information processing with superconducting circuits: a review', Rep. Prog. Phys. 80... |
| `po_amplitude_oscillator.json` | Josephson parametric oscillator for qubit measurement, amplitude encoded (Krantz 2016) | x, y | eps | Krantz, P., Bengtsson, A., Simoen, M., Gustavsson, S., Shumeiko, V., Oliver, W. D., Wilson, C. M., Delsing... |
| `pplo_oscillator.json` | Josephson parametric phase-locked oscillator for qubit measurement (Lin 2014) | qx, qy | p | Lin, Z. R., Inomata, K., Koshino, K., Oliver, W. D., Nakamura, Y., Tsai, J. S., Yamamoto, T., 'Josephson pa... |
| `qfp_latch.json` | quantum flux parametron (QFP) latch | th, v | f | Berkley, A. J., Johnson, M. W., Bunyk, P., Harris, R., Johansson, J., Lanting, T., Ladizinsky, E., Tolkache... |
| `quantum_vdp.json` | self-sustained oscillator with one-quantum gain and two-quantum loss (maser, phonon laser) | x, y | F | Lee & Sadeghpour, PRL 111, 234101 (2013); Walter, Nunnenkamp & Bruder, PRL 112, 094102 (2014) |
| `rfsquid_annealer.json` | rf-SQUID annealer qubit (CCJJ), classical phase particle | th, v | beta | Harris, R. et al., 'Experimental demonstration of a robust and scalable flux qubit', Phys. Rev. B 81, 13451... |
| `rfsquid_pair.json` | two inductively coupled rf-SQUID annealer qubits (antiferromagnetic) | th1, v1, th2, v2 | beta | Harris, R. et al., 'Experimental demonstration of a robust and scalable flux qubit', Phys. Rev. B 81, 13451... |
| `rfsquid_triangle.json` | three rf-SQUID annealer qubits on an antiferromagnetic triangle | th1, v1, th2, v2, th3, v3 | beta | Harris, R. et al., 'Experimental demonstration of a robust and scalable flux qubit', Phys. Rev. B 81, 13451... |
| `rydberg_bistability.json` | driven-dissipative Rydberg ensemble (mean-field optical bistability) | u, v, w | Delta | Carr, C., Ritter, R., Wade, C. G., Adams, C. S., Weatherill, K. J., 'Nonequilibrium phase transition in a d... |
| `snail_po.json` | SNAIL parametric oscillator (three-wave mixing, small Kerr) | x, y | g | Frattini, N. E., Vool, U., Shankar, S., Narla, A., Sliwa, K. M., Devoret, M. H., '3-wave mixing Josephson d... |
| `spin_ensemble_cavity.json` | spin ensemble strongly coupled to a microwave cavity (driven Tavis-Cummings, absorptive bistability) | x, y, sx, sy, sz | eta | Angerer et al., Sci. Adv. 3, e1701626 (2017) |
| `transmon_driven.json` | driven transmon (Duffing image) | x, y | F | Krantz, P., Kjaergaard, M., Yan, F., Orlando, T. P., Gustavsson, S., Oliver, W. D., 'A quantum engineer's g... |
