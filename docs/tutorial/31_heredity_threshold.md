# Heredity, Module 1. Inheritance through a threshold

**Learning objectives.** After this module you can

1. describe a body that grows and divides, and whose order exists only above a critical size, in a specification of
   schema `fieldbridge-heredity/1` (equations, a one-dimensional field on cells, or an elastic chain);
2. state from the equations whether a daughter can inherit its parent's order: the threshold and its pitchfork, the
   dip that division creates, the gain of the order over a generation, and whether division moves the threshold;
3. compute the probability that a daughter keeps its parent's sign, P = Φ(φ<sub>c</sub>/σ<sub>c</sub>), and compare it
   with lineages of the full body;
4. run `heredity predict`, `heredity card` and `heredity survey`, and read the class and the conditions of the law.

**Prerequisites.** Modules [1](15_memory_first_card.md), [2](16_memory_specification.md) and
[4](18_memory_writing_and_retention.md) of the memory part, for the realization, its specification and the swept
write. **Time.** About 45 minutes.

## 1. The question

Many bodies keep an order only above a critical size. A reaction volume with autocatalysis keeps an excess of one
handedness only when it holds enough molecules; a domain of a Turing system holds its peak at one end only when it is
long enough; a stripe of an active nematic flows only when it is wide enough; a ferroelectric film is polar only above
a critical thickness; a filament under a fixed load buckles only above a critical length. When such a body grows and
divides in the middle, each daughter can be smaller than the critical size. Its order then decays, and the daughter
has to grow through the threshold again before it can be ordered. The question of this module is whether, and with
what probability, the daughter's new order has the sign of its parent's.

## 2. The law

Near the threshold L<sub>c</sub> of the symmetric state, the order s, the projection of the state on the critical
mode, obeys the normal form of a pitchfork with the size as its control,

$$\dot s = a\,(L - L_c)\,s - b\,s^3 + \sqrt{2D_s}\,\xi(t),$$

with b > 0 for a supercritical pitchfork. A daughter is born at L<sub>div</sub>/2 < L<sub>c</sub> with the order of its
parent. In the dip its order decays; growth carries it back to L<sub>c</sub> with the order φ<sub>c</sub>, and noise
accumulated in the slow passage through the threshold adds a Gaussian deviation of standard deviation σ<sub>c</sub>.
The daughter keeps its parent's sign with probability

$$P = \Phi\!\left(\frac{\varphi_c}{\sigma_c}\right),\qquad
\sigma_c^2 = 2D_s\sqrt{\frac{\pi}{a r}}\;\Phi\!\left(a\mu_0\sqrt{\frac{2}{a r}}\right),$$

where Φ is the normal distribution function, r = dL/dt at L<sub>c</sub>, and μ<sub>0</sub> = L<sub>c</sub> −
L<sub>div</sub>/2 the depth of the dip. For a sustained bias h in place of an initial order the same calculation gives
the law of the swept write of Module 4, P = Φ(π<sup>1/4</sup>h/(D<sup>1/2</sup>(ar)<sup>1/4</sup>)) (Kondepudi and
Nelson, 1983); the Gaussian statistics of a slow passage from an offset are those of Lythe (1996) and of the delayed
bifurcations with noise (van den Broeck and Mandel, 1987).

φ<sub>c</sub> must be the order that the body's own equations carry to the crossing. The normal form with a linear
ramp gives a decay e<sup>−aμ<sub>0</sub><sup>2</sup>/2r</sup> and, with its cubic term, a faster decay of a large order;
where the dip is deep, the critical eigenvalue of the body is not linear in the size (diffusion enters as 1/L², bending
as 1/L⁴), and that decay is wrong. The card therefore computes φ<sub>c</sub> by integrating the body without noise
from the daughter's state at birth to L<sub>c</sub> (the version "body"), and reports the normal form beside it.

The law needs three conditions, which the card checks:

1. **Division leaves the threshold in place.** The daughter starts near the symmetric state of its own size. When the
   order is a redistribution of a conserved amount that sets the threshold, division gives the two halves different
   amounts and moves their thresholds instead.
2. **The order survives a generation without noise.** The critical eigenvalue λ integrated over a generation,
   ln G = ∫λ dt from L<sub>div</sub>/2 to L<sub>div</sub>, must be positive: the regrowth above L<sub>c</sub> must
   restore what the dip removes.
3. **The noise is small in the crossing window.** Λ = a r/(b D<sub>s</sub>) ≫ 1. When it is not, the new order is
   flipped by the noise just past the threshold, and fewer daughters keep the sign than the law gives.

## 3. The prediction from structure

`heredity predict` finds the threshold along the size (the largest eigenvalue of the Jacobian at the symmetric state
crosses zero), reduces the body there to a, b and the critical vectors, and follows one lineage without noise: a parent
that leaves the threshold with the noise amplitude σ<sub>c</sub>, as a noisy parent does, grows to L<sub>div</sub> and
divides, and its daughters follow. The class:

| Class | Condition |
| --- | --- |
| kept above the threshold | L<sub>div</sub>/2 ≥ L<sub>c</sub>: no dip; the order passes behind its barrier |
| lost in the dip | ln G < 0: without noise the order dies out over the generations |
| threshold moved by division | ln G > 0 along the symmetric state, but the lineage without noise loses the order |
| inherited through a threshold | the daughter crosses again and keeps the sign with P = Φ(φ<sub>c</sub>/σ<sub>c</sub>) |

```bash
python -m fieldbridge heredity predict examples/heredity/benchmarks/normal_form.json
```

For the normal form with a = b = L<sub>c</sub> = 1, exponential growth at the rate 0.05 and division at 1.5 the
daughter is born at 0.75, ln G = (0.75 − ln 2)/0.05 = 1.14, Λ = 16.7 at the noise 3·10⁻³, and the law gives
P = 0.787 for the daughter of the lineage without noise.

## 4. Running the command

A specification has the header of every example, a key `body` and a block `lineage`:

| Key | Meaning |
| --- | --- |
| `size` | the declared parameter that growth raises and division halves |
| `threshold` | the interval in which the threshold is sought |
| `growth` | `exponential` or `linear`, its rate (a number or a parameter), and for a field the mode `uniform` or `faces` |
| `division` | `at` a size or `at_relative` to L<sub>c</sub>; for equations, `keep` or `halve` per variable and a `binomial` partition |
| `noise` | a number, or for equations a variance rate per variable (for example the chemical Langevin noise) |
| `symmetric_state`, `dt`, `start_spread`, `expect` | a starting point for the symmetric state, the time step, the initial spread of the order, the expected class |

The three bodies: `equations` (the drift of the memory module, read by the same restricted parser), `field` (variables
on cells with local reactions, a mobility that may depend on the local variables, and faces with zero flux or an
extrapolation length) and `chain` (an inextensible elastic filament under a dead load, with anisotropic drag).

```bash
python -m fieldbridge heredity card examples/heredity/chiral_autocatalysis_saito2007.json --out-dir out/chiral
```

The reaction volume of Saito, Sugimori and Hyuga (2007) keeps an excess of R or S only above a critical content of
39.5 molecules per source volume. Supplied with substrate and transferred at the content 64, it is diluted to 32: the
card gives L<sub>c</sub> = 39.53, a = 0.0316, b = 7.5·10⁻⁴, ln G = 6.52 and, in a vessel 16 times the source volume,
Λ = 42.7 and P = 0.886.

## 5. Controls

Each control changes one component and changes the class:

| Control | Change | Class |
| --- | --- | --- |
| `controls/normal_form_no_dip.json` | P: division at 1.5 → 2.4 L<sub>c</sub> | inherited → kept above the threshold |
| `controls/normal_form_lost.json` | P: division at 1.5 → 1.3 L<sub>c</sub> (ln G = −0.86) | inherited → lost in the dip |
| `controls/chiral_autocatalysis_no_dip.json` | P: transfer at the content 64 → 100 | inherited → kept above the threshold |
| `controls/filament_short_division.json` | P: severing at 1.6 → 1.4 L<sub>c</sub> (ln G = −4.1) | inherited → lost in the dip |
| `polarity_brauns2020.json` | an order made of a conserved protein | threshold moved by division |

The polarity of Brauns, Halatek and Frey (2020) exists only within a window of mean densities. A polarized cell
divided in the middle gives one half more protein than the other, both outside the window; neither daughter polarizes
again from the parent's order, although ln G = 3.5 along the symmetric state. The order at the successive divisions of
the lineage without noise falls from 1.3 to below 10⁻⁷: the parent's polarity becomes a difference between the
daughters, an asymmetric division.

## 6. Lineages

```bash
python -m fieldbridge heredity card examples/heredity/filament_baczynski2007.json --simulate --lineages 1000 --seed 2 --out-dir out/fil
```

With `--simulate` the card runs lineages of the full body with noise over three generations (the first not counted)
and compares the fraction of daughters that keep their parent's sign with the law averaged over the daughters' states
at birth, within |P − P<sub>law</sub>| < 4·stderr + 0.02. With 1000 lineages:

| Body | Field | Measured | Law (body) | Law (history) | Normal form | Λ |
| --- | --- | --- | --- | --- | --- | --- |
| Pitchfork normal form | — | 0.715 ± 0.014 | 0.752 | 0.751 | 0.781 | 17 |
| Chiral autocatalysis, vessel of 16 source volumes | chemistry | 0.877 ± 0.010 | 0.885 | 0.882 | 0.864 | 43 |
| Turing domain of one half-stripe | developmental pattern formation | 0.762 ± 0.013 | 0.773 | 0.769 | 0.822 | 70 |
| Active nematic stripe | active matter | 0.924 ± 0.008 | 0.926 | 0.925 | 0.935 | 64 |
| Ferroelectric film, δ = ξ | condensed-matter physics | 0.909 ± 0.009 | 0.898 | 0.898 | 0.909 | 61 |
| Filament under a dead load | polymer mechanics | 0.693 ± 0.015 | 0.692 | 0.690 | 0.823 | 447 |

All six agree within the gate, in the body and the history versions of the law (`--seed 2`). For the normal form,
with Λ = 17, two runs of 4000 lineages give 0.736 and 0.732 against 0.748 and 0.746: when Λ is not large the new
order is sometimes flipped just past the threshold, and the measured fraction falls slightly below the law.

The law holds when φ<sub>c</sub> comes from the body's own equations. The normal form with a linear ramp is off by
0.13 for the filament, whose critical eigenvalue goes as L⁻⁴ through the dip, and by 0.06 for the Turing domain.

Growth, division and, except for the chemical noise of the reaction volume, the noise are assumptions stated in each
specification; none of the sources has them. With the published time scales of the three physical bodies a daughter's
residual order relaxes much faster than the daughter regrows (for a microtubule severed under load, by a factor of
about 10⁴), so its order is lost and P is ½. The examples therefore take growth within a few relaxation times of the
critical mode, as in a viscous or crowded medium; the law gives P for any growth rate.

## 7. The published bodies

```bash
python -m fieldbridge heredity survey examples/heredity
```

The survey gives the class and the law for every specification in a folder:

| Body | Field | Class | L_c | L_div / L_c | ln G | Λ | P (body) | P (normal form) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Pitchfork normal form with the size as its control | — | inherited through a threshold | 1 | 1.5 | 1.14 | 16.7 | 0.787 | 0.819 |
| Chiral autocatalysis transferred at a content above twice the critical content | chemistry | kept above the threshold | 39.53 | 2.53 | — | 42.7 | — | — |
| Filament severed at 1.4 times its critical length | polymer mechanics | lost in the dip | 0.9959 | 1.4 | -4.11 | 447 | — | — |
| Pitchfork normal form divided at 1.3 times the threshold | — | lost in the dip | 1 | 1.3 | -0.863 | 16.7 | — | — |
| Pitchfork normal form divided above twice the threshold | — | kept above the threshold | 1 | 2.4 | — | 16.7 | — | — |
| Active nematic stripe that flows only above a critical width | active matter | inherited through a threshold | 3.183 | 1.6 | 0.895 | 64 | 0.949 | 0.961 |
| Chiral autocatalysis with recycling in a reaction volume diluted by serial transfer | chemistry | inherited through a threshold | 39.53 | 1.62 | 6.52 | 42.7 | 0.886 | 0.865 |
| Ferroelectric film polar only above a critical thickness (dimensionless) | condensed-matter physics | inherited through a threshold | 1.594 | 1.6 | 0.791 | 61 | 0.944 | 0.957 |
| Filament that buckles only above a critical length under a dead compressive load | polymer mechanics | inherited through a threshold | 0.9959 | 1.6 | 0.267 | 447 | 0.759 | 0.946 |
| Cell polarity of a mass-conserving reaction-diffusion system | cell biology | threshold moved by division | 8.855 | 1.58 | 3.52 | 9.2e+03 | — | — |
| Turing system on a growing domain of one half-stripe (Lengyel-Epstein kinetics) | developmental pattern formation | inherited through a threshold | 0.2955 | 1.69 | 2.46 | 69.9 | 0.896 | 0.903 |

P is the law for the daughter of the lineage without noise of Section 3; the lineages of Section 6 average it over the
daughters of noisy parents.

The ferroelectric film is dimensionless: the sources read give no material values of the gradient coefficient and the
extrapolation length δ, and δ/ξ = 1 is an assumption (δ > 0 suppresses the polarization at the faces and makes the
threshold).

## 8. A model from your field

Write the equations of your body, name the extensive size that division halves as a parameter, and give the interval
in which the threshold lies. Run `heredity predict` first: it tells you whether the threshold is supercritical,
whether division leaves it in place, and whether the order survives a generation without noise. Then run
`heredity card --simulate` at the noise of your system. A specification in `examples/heredity` appears on the web page
in the column of its class.

## 9. Exercises

1. In `benchmarks/normal_form.json`, find the smallest division size at which the order survives without noise.
   Show that for exponential growth ln G = (L<sub>div</sub>/2 − L<sub>c</sub> ln 2)/r, so that the bound is
   L<sub>div</sub> = 2 L<sub>c</sub> ln 2 = 1.386 L<sub>c</sub>, whatever the growth rate.
2. Double the growth rate of the Turing domain. Does P rise or fall? Which term of the law changes most?
3. Give the filament a smaller load. How does the critical length change, and does the class change at the same
   division size?
4. Run the chiral volume at a vessel of 1 instead of 16 source volumes. What is Λ, and how do the lineages compare
   with the law?

## Summary

- A body whose order exists only above a critical size passes its order to a daughter born below that size with
  P = Φ(φ<sub>c</sub>/σ<sub>c</sub>): the order that the body's own equations carry to the threshold against the noise
  of the slow passage.
- The law needs a supercritical threshold that division leaves in place, a positive gain of the order over a
  generation, and small noise in the crossing window (Λ ≫ 1).
- The controls change the class by one component: a division above twice the threshold keeps the order without a dip;
  a division too close to the threshold loses it; an order made of a conserved amount is converted into a difference
  between the daughters.

## Reference

`fieldbridge.heredity`: `spec` (the schema), `bodies` (equations, fields, chains), `reduce` (the threshold and its
reduction), `law`, `lineage` (lineages and the law in three versions), `predict`, `card`, `cli`. Examples in
`examples/heredity`: six published bodies, the normal form and four controls. Sources: Saito, Sugimori and Hyuga (2007),
Painter, Maini and Othmer (1999), Duclos et al. (2018), Voituriez, Joanny and Prost (2005), Kretschmer and Binder
(1979), Tilley and Zeks (1984), Baczynski, Lipowsky and Kierfeld (2007), Hallatschek, Frey and Kroy (2007), Grassia,
Hinch and Nitsche (1995), Brauns, Halatek and Frey (2020), Kondepudi and Nelson (1983), Lythe (1996), van den Broeck
and Mandel (1987); see docs/mechanisms.md.
