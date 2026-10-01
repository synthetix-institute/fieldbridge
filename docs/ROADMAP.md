# Roadmap

## Mechanisms in preparation

A mechanism enters the map of the web page and `memory codiscover` when FieldBridge derives it in models from
several fields and certifies a law whose constant does not depend on the field: a canonical form, a law checked by
simulating the full model, and letters that verify each step. Four laws meet this rule: the swept write (constant
π<sup>1/4</sup>), the delayed fold (|a₁′| = 1.0188) and the half-width of phase locking (1 in units of K), which are
the targets of `fieldbridge/memory/codiscovery.py`, and the Rabi law of the Bloch rotation
(`fieldbridge/quantum/language.py`). The mechanisms below are studied in the tutorial or in the literature on
memory and are not yet derivation targets. The web page lists them under "Mechanisms in preparation".

| Mechanism | Canonical form or law | What exists | What is missing |
| --- | --- | --- | --- |
| Frustrated loops | excess energy per bond 1 − cos(Φ/N) for N equal rotor bonds with mismatch Φ | the structural prediction and `memory loops` ([Module 3](tutorial/17_memory_predictions.md#4-tests-on-many-loops-and-networks)): 91 loops of genes, spins and rotors behave as predicted; rotor frustration agrees with 1 − cos(Φ/N) to 3 × 10⁻¹² | a derivation target with its letters; a browser engine for networks |
| Retention against rewriting | the ratio of the retention time to the writing time depends only on the work E<sub>w</sub> of the write: it grows as ln E<sub>w</sub> in a curved minimum, linearly along a zero mode and as e<sup>E<sub>w</sub>/k<sub>B</sub>T</sup> behind a barrier ([Module 4, Section 1.4](tutorial/18_memory_writing_and_retention.md#14-relation-between-retention-and-writing-times)) | the memory card computes the ratio for each writing protocol | a target in which the write of bounded work is a letter; its comparison across fields |
| Onset of oscillation (Hopf) | ż = (μ + iω)z − \|z\|²z: an amplitude proportional to √μ; under a slow sweep of μ the onset is delayed (Neishtadt, 1987) | Hopf bifurcations are located and reported where they stop a write, as in the ring of three repressors | the reduction to the complex amplitude; a certified law |
| Synchronization of a population (Kuramoto) | oscillators with a symmetric, unimodal density g of natural frequencies synchronize above K<sub>c</sub> = 2/(π g(0)) (Kuramoto, 1984) | the phase reduction and the law of locking of one oscillator ([Module 11](tutorial/26_memory_phase_locking.md)) | a carrier for populations of oscillators |
| Return-point memory | rate-independent hysteresis: the state returns to the same configuration when the field returns to an earlier extremum, and it is set by the sequence of extrema (Sethna et al., 1993; Keim et al., 2019) | nothing yet | a specification of rate-independent elements (hysterons) driven without thermal noise |

The memory specifications (`fieldbridge/memory/spec.py`) describe drift equations on a Euclidean, torus or orthant
carrier and networks of rotors, spins or genes; `fieldbridge/memory/fields.py` reads linear lattice fields. Return-point
memory, and other rate-independent or athermal materials, need a new kind of specification.

Sources: A. I. Neishtadt, Differential Equations 23, 1385 (1987); Y. Kuramoto, *Chemical Oscillations, Waves, and
Turbulence* (Springer, 1984); J. P. Sethna et al., Phys. Rev. Lett. 70, 3347 (1993); N. C. Keim, J. D. Paulsen,
Z. Zeravcic, S. Sastry and S. R. Nagel, Rev. Mod. Phys. 91, 035002 (2019). The original publications of these and of
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
