# Reviewing the claims of a manuscript

**Learning objectives.** After this chapter you can

1. write the claims of a manuscript in a review file, with the system on which each quantity was computed and the
   system for which it is claimed;
2. run the four checks of `fieldbridge review` and read what each failure means;
3. see why each check follows from the nested description of a model, and what the checks cannot find.

**Prerequisites.** The nested description of [the system map](00_system_map.md): an operation acting on a carrier
(the states), completed by closures (the constraints), an observable and a protocol (preparation and observation), and
realized by a material (the parameters). **Time.** About 20 minutes; the example runs in a second.

## 1. Why four checks

Every quantity in a model belongs to a level of the nested description, and a claim holds only at the level at which
the quantity was computed. Manuscripts most often go wrong where a claim moves between levels without saying so, or
leans on a result outside its range. Four checks test these places.

1. **Same system.** Was the quantity computed for the system it is claimed for? A count made while some parameters
   were free is not a count for the system in which they are fixed. The compression principle of the constructor
   states the general rule: a description may drop a distinction only while that distinction cannot change the
   predicted measurement. Where the review gives equation and unknown counts, the check also asks whether the claimed
   state exists for a prescribed system at all.
2. **Independence of the protocol.** Does the quantity change with the run length, the rate or the preparation? A
   property of the state does not; a property that does belongs to the protocol and the observable as well. In
   particular, a thermodynamic statement cannot rest on a kinetic observable alone: how long a prepared state survives
   depends on how it was prepared and on how long one watches. The check also tests a fitted divergence against the
   observation window: whether it sits at the edge of the window, carries an uncertainty smaller than the unmeasured
   interval, is determined by the fit at all, or competes with a form that has no divergence.
3. **Scope of the theorems.** Is every theorem applied inside its hypotheses, and does it support the conclusion
   drawn? The registry (`fieldbridge review theorems`) states the hypotheses of each theorem as conditions on the case,
   with the conclusions it supports and conclusions often drawn from it that it does not support.
4. **Prior art.** Does each claim of novelty, and each statement the review relies on, carry a recorded search with a
   verdict (published, partly published, not found), and does any published result contradict the claim?

These are necessary conditions. A failed check names a fault; a passed check does not make a claim right. Arithmetic
errors in a derivation, numerical artefacts and the physics that only a calculation reveals are outside them, and the
review file records the calculations and readings on which its entries rest.

## 2. The review file

A review file (schema `fieldbridge-review/1`) has three parts.

- `paper`: title, reference, and where the numbers in the file come from.
- `systems`: each system on which a quantity is computed or for which it is claimed, with its carrier (state
  variables, free and fixed parameters), closures, realization and, where useful, the numbers of equations and
  unknowns as expressions in `N`.
- `claims`: each claim with its `kind` (count, structural, mechanical, kinetic, thermodynamic), the systems
  `computed_on` and `claimed_for`, the `observable_type` when a thermodynamic claim rests on a measured time, the
  protocol-dependent estimates or the observation window of a fitted divergence, the `theorems` used with the case and
  the inference, and the `prior_art` searches.

## 3. A worked example

`examples/review/ideal_glass_2404.07492.json` records six claims of arXiv:2404.07492, which builds triangulated
packings of polydisperse disks by letting the radii vary during preparation and concludes that they are ideal glasses
with zero configurational entropy.

```bash
python3 -B -m fieldbridge review check examples/review/ideal_glass_2404.07492.json --out-dir build/review
```

The zero-entropy claim fails the first check: it is computed while the radii are free, and claimed for disks whose
radii are fixed. The counts show why the difference matters: a triangulation of the torus gives 3N contact equations,
which with the radii fixed act on only 2N + 1 unknowns, so a prescribed set of radii generically has no triangulated
packing, and the radii of each such packing are tuned to it. The same claim fails the third check twice: the Gaussian
count it uses is a theorem only inside 1 < m/n < 3, not at the triangulated end m/n = 3, and a vanishing fraction of
graphs can still be an exponential number of them.

The melting-temperature claim fails the second check. It is a thermodynamic statement resting on a fitted divergence
of a relaxation time, without a free energy or an entropy, and in an independent reproduction of the measurement the
fitted divergence sits at the edge of the observation window, on a bound of the fit, with a residual that changes by
less than a factor of two over most of the range, while a form without a divergence fits as well. The statement about
the Mermin–Wagner–Hohenberg theorem fails the third check: a positive melting temperature is not long-range order of a
continuous symmetry and cannot violate it. Published results contradict three claims (fourth check).

The three claims computed on the fixed, tuned system warn under the first check rather than fail: they are made for
the system on which they were computed, but that system exists only because its radii were tuned.

## 4. Exercises

1. Add to the example a claim that the bulk modulus is finite at zero pressure, computed and claimed for the tuned
   system, using `maxwell_constraint_counting` with the inference `finite_shear_modulus_at_zero_pressure_when_hyperstatic`.
   Which checks pass, and why does a pass not support the ideal-glass interpretation?
2. Give the melting-temperature claim three estimates of T_m from runs of increasing length that drift downward, and
   read the slope reported by the second check.
3. Write a review file for a model of your own field with one claim that moves between systems, and see which part of
   the first check's explanation names the difference.
