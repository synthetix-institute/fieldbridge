# Learning objectives as material specifications

Two specifications written from the objective of Balestriero and LeCun, *Proving and Improving the Platonic
Representation Hypothesis with JEPAs*, arXiv:2610.39302v1. The paper states a population objective on a finite
world and its minimizers. The specifications are the gradient flow of that objective (Proposition 1 in its spectral
contrastive form, centred embeddings) written on the eigenfunctions of the positive-pair operator; the batch noise is
taken as additive. No network is trained.

| File | Setting of the paper | Variables | Control |
| --- | --- | --- | --- |
| `jepa_views.json` | a fraction q of good views, the rest product corruptions (Proposition 4 ii, Theorem 5), one embedding direction | w on the scene, u on the nuisance | q |
| `jepa_rotation.json` | shared positives (Theorem 1), embedding dimension two, a world of four states | six coefficients | third eigenvalue m3 |

```bash
python3 -B -m fieldbridge memory card examples/memory/learning_objectives/jepa_views.json --out-dir /tmp/jepa_views --quick
python3 -B -m fieldbridge memory card examples/memory/learning_objectives/jepa_rotation.json --out-dir /tmp/jepa_rotation --quick
```

The first card places the exchange of stability between the nuisance and the scene at
q* = (beta - alpha)/(beta - alpha + mu) = 0.429, where the drift vanishes along the critical direction, and gives the
eigenvalue gap as the slowest rate at a stable state. The second returns a continuum of neutral states, the rotation
of the theorem, with diffusion along it as the law of loss, and the gap m2 - m3 as the slowest rate of a distortion.
The structural pre-screen tests permutations and sign reversals and does not report the rotation symmetry; the
calculation finds the zero mode.

These files are in a subfolder and are not part of the catalog of materials.
