# Roadmap

## Mechanisms in preparation

A mechanism enters the map of the web page and `memory codiscover` when FieldBridge derives it in models from
several fields and certifies a law whose constant does not depend on the field: a canonical form, a law checked by
simulating the full model, and letters that verify each step. Four laws meet this rule: the swept write (constant
π<sup>1/4</sup>), the delayed fold (|a₁′| = 1.0188) and the half-width of phase locking (1 in units of K), which are
the targets of `fieldbridge/memory/codiscovery.py`, and the Rabi law of the Bloch rotation
(`fieldbridge/quantum/language.py`). Return-point memory is on the map with a structural condition in place of a law constant:
the return to a turning point of a slow drive is exact when no loop of the network that counts the drive as an
element is frustrated ([Module 12](tutorial/27_memory_return_point.md), `memory hysterons`), and six realizations
from three fields are compared with it. The return to a set point after a step of an input is a second family with a
structural condition: an integrator together with a stable steady state makes the return exact
([Regulation, Module 1](tutorial/28_regulation_set_point.md), `regulation card`), and published models from five
fields are compared with it, in the column "regulation: set point" of the map on the web page. Computation by a
driven body is a third family whose classes follow from structure: fading memory over the modes that the input
reaches and the observables see, an affine body that holds only degree 1, and an odd body that holds only odd
degrees ([Computation, Module 1](tutorial/29_computation_capacity.md), `computation card`); published models from six
fields are compared with it, in the column "computation: capacity". The quantum write is a fourth addition: a
parametric oscillator swept through its threshold with a bias follows the classical write law with the noise fixed
by the loss and the temperature, exact for a quadratic generator, and leaves it for the selection of its biased
steady state when the stored states hold a few photons ([Quantum write](tutorial/29_quantum_write.md),
`quantum write`), in the column "quantum: open evolution". The mechanisms below are studied in the tutorial or in the literature on
memory and are not yet derivation targets. The web page lists them under "Mechanisms in preparation".

| Mechanism | Canonical form or law | What exists | What is missing |
| --- | --- | --- | --- |
| Frustrated loops | excess energy per bond 1 − cos(Φ/N) for N equal rotor bonds with mismatch Φ | the structural prediction and `memory loops` ([Module 3](tutorial/17_memory_predictions.md#4-tests-on-many-loops-and-networks)): 91 loops of genes, spins and rotors behave as predicted; rotor frustration agrees with 1 − cos(Φ/N) to 3 × 10⁻¹² | a derivation target with its letters; a browser engine for networks |
| Retention against rewriting | the ratio of the retention time to the writing time depends only on the work E<sub>w</sub> of the write: it grows as ln E<sub>w</sub> in a curved minimum, linearly along a zero mode and as e<sup>E<sub>w</sub>/k<sub>B</sub>T</sup> behind a barrier ([Module 4, Section 1.4](tutorial/18_memory_writing_and_retention.md#14-relation-between-retention-and-writing-times)) | the memory card computes the ratio for each writing protocol | a target in which the write of bounded work is a letter; its comparison across fields |
| Onset of oscillation (Hopf) | ż = (μ + iω)z − \|z\|²z: an amplitude proportional to √μ; under a slow sweep of μ the onset is delayed (Neishtadt, 1987) | Hopf bifurcations are located and reported where they stop a write, as in the ring of three repressors | the reduction to the complex amplitude; a certified law |
| Synchronization of a population (Kuramoto) | oscillators with a symmetric, unimodal density g of natural frequencies synchronize above K<sub>c</sub> = 2/(π g(0)) (Kuramoto, 1984) | the phase reduction and the law of locking of one oscillator ([Module 11](tutorial/26_memory_phase_locking.md)) | a carrier for populations of oscillators |

The memory specifications (`fieldbridge/memory/spec.py`) describe drift equations on a Euclidean, torus or orthant
carrier and networks of rotors, spins or genes; `fieldbridge/memory/fields.py` reads linear lattice fields, and
`fieldbridge/memory/hysterons.py` interacting hysterons under a slow drive (kind `hysterons`).

Sources: A. I. Neishtadt, Differential Equations 23, 1385 (1987); Y. Kuramoto, *Chemical Oscillations, Waves, and
Turbulence* (Springer, 1984). The original publications of these and of
the mechanisms on the map are listed, with their DOIs, in [mechanisms.md](mechanisms.md#mechanisms-in-preparation).

## Version 0.1

- deterministic route/fiber fingerprint;
- three field packs;
- PDF-folder field adapter builder;
- field-native constructor role and substrate profiles;
- mechanism extraction from one paper or fragment;
- paper-to-paper mechanism comparison;
- analog search over existing mechanism records;
- cross-field mechanism translation sheet;
- clear evidence boundary.

## Version 0.2

- import larger mechanism indexes from JSON or Parquet;
- attach public mechanism-page links;
- attach arXiv and DOI references;
- support patent/invention claim language as a fourth field pack;
- export Markdown and JSON reports.
- corpus-level adapter quality checks.

## Version 0.3

- optional vector database backend;
- high-dimensional fingerprints;
- k-nearest-neighbor search over millions of equation records;
- confidence decomposition by route, fiber, equation evidence, and reference
  support.

## Version 0.4

- web demo;
- community field-pack registry;
- reviewer mode for checking whether a proposed analogy is only word-level or
  genuinely mechanism-level.

## Scientific Boundary

FieldBridge should never claim that an analogy is already experimentally true.
It should produce a mechanism, equations, variables, controls, and references
that make the claim testable.
