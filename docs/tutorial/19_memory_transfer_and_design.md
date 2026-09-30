# Module 5. Transfer to another carrier and removal of an obstruction

**Learning objectives.** After this module you can

1. transfer a known memory mechanism to another carrier and read which of its properties are kept and which
   change;
2. identify the obstruction that prevents a symmetric write on the new carrier;
3. find the parameter setting at which a one-sided write becomes symmetric, and the direction of a sweep through
   it.

**Prerequisites.** Modules 1–4. **Time.** About 40 minutes.

## 1. Concepts

A memory mechanism can be separated from the carrier on which it was found and attached to another carrier:

$$
I^{s}_{\mathrm{real}} \;\to\; I^{s}_{\mathrm{op}} \;\to\; (\Omega_s, 0) \;\to\; (\Omega_s, \Xi_t) \;\to\; I^{t}_{\mathrm{real}} .
$$

From the source realization $s$ the program keeps its memory signature: the properties that follow from the
structure of the operation $\Omega_s$ (Module 3). It then attaches the operation to the target carrier $\Xi_t$ and
closes it with the bath, observable and protocol of the target $t$. The report lists three groups.

| Group | Content | Key in `attach.json` |
| --- | --- | --- |
| Kept | Entries of the signature shared by source and target: reciprocity, kind of write, retention law, oscillation, multistability, role of the control, zero modes, satisfiability of loops | `kept` |
| Changed | Entries that differ; each is a prediction about the target that a calculation or an experiment must confirm | `changed` |
| Carrier and material | Differences that are expected by construction, such as the period of the angles or the ratio of two couplings | `carrier_and_material` |

The role of the control takes one of three values. A **scale** control multiplies the whole drift and acts as an
inverse temperature (Module 4). A **bias** control adds a constant tilt to the landscape, so a sweep of it favours
one state. A **shape** control changes the form of the landscape and can create states at a write point.

If a changed entry is a one-sided write, the program can search for the setting of a second parameter at which
the quadratic term $a_2$ of the normal form vanishes. At such a point, a cusp, the fold becomes a pitchfork. The
program also returns the direction of a sweep through the cusp along which the bias on the unstable direction
does not change.

## 2. Worked example: from capillary rods to point dipoles

Point dipoles in a plane, at the positions of the colloids of Module 4, interact through

$$
E_{ij} = -\frac{\mu^2}{r_{ij}^3}\Bigl[\tfrac12\cos(\psi_i-\psi_j) + \tfrac32\cos(\psi_i+\psi_j-2\varphi_{ij})\Bigr],
$$

where $\psi_i$ is the direction of dipole $i$ and $\varphi_{ij}$ the direction of the bond. The first term aligns
the two dipoles; the second reflects the direction of one dipole across the bond axis. The capillary rods have
couplings of the same two kinds, with a different ratio and range, and with angles of period $\pi$ instead of
$2\pi$.

```bash
python3 -B -m fieldbridge memory attach --from examples/memory/colloid_patch.json \
  --to examples/memory/dipole_patch.json --out-dir build/attach_dipoles
```

| Entry (`attach.json`) | Value | Interpretation |
| --- | --- | --- |
| `changed` | empty | The signature is kept: reciprocal couplings, a scale control, discrete states written by a field, activated retention ([Law 3](18_memory_writing_and_retention.md#13-retention)), no oscillation |
| `carrier_and_material.carrier` | period $\pi \to 2\pi$ | Apolar rods are replaced by polar dipoles |
| `carrier_and_material.reflect_to_align` | median 1.25 (range 0.50–7.5) $\to$ 3 | The capillary ratio $5(r_0/r)^4 : 4$ depends on the distance; the dipolar ratio is fixed at $3:1$ |

The transfer predicts that an array of dipoles at these positions stores discrete orientation patterns written
by a field and retained by activation. It does not predict the numbers: the ratio of the two couplings sets the
frustration of each loop and therefore the energies and barriers, which must be calculated for the target
(`memory card examples/memory/dipole_patch.json`). The gallery of Module 7 contains this calculation: sampling finds
ten stable states of the dipole network, written by a uniform field and retained by activation (Law 3), as the
transfer predicted.

## 3. Worked example: from the toggle switch to a chemical reactor

```bash
python3 -B -m fieldbridge memory attach --from examples/memory/toggle.json \
  --to examples/memory/schlogl.json --out-dir build/attach_schlogl
```

| Changed entry | Toggle switch | Schlögl reactor |
| --- | --- | --- |
| Writing | a pitchfork where the symmetric state loses stability along a reversed mode, selected by a weak bias | one-sided writes at folds; a symmetric write needs one tuned parameter (a cusp) |
| Role of the control | shape (the promoter strength creates the second state) | bias (the feed $b$ tilts the landscape) |
| Reciprocal couplings | no | yes (one variable: a gradient drift) |
| Oscillation | not expected (no negative loop) | impossible |
| Multistability | possible (a positive loop) | not predicted: the rule on loop signs is not applied to a gradient drift |

The reactor has no exchange symmetry, so the symmetric write of the toggle switch is not kept. Both the fold and
the bias of the control are obstructions to a symmetric write.

## 4. Removing the obstruction

```bash
python3 -B -m fieldbridge memory design examples/memory/schlogl.json --free a --out-dir build/design_schlogl
python3 -B -m fieldbridge memory design examples/memory/tubes_unequal.json --free L2 --out-dir build/design_tubes
python3 -B -m fieldbridge memory design examples/memory/toggle_unequal.json --free gamma --out-dir build/design_toggle
```

The program follows the write point as the free parameter changes until $a_2$ vanishes, refines the point, and
reports the kind of write along the sweep through it (`design.md`).

| Realization | Free parameter | Fold (obstruction) | Cusp | Sweep direction | Kind along the sweep |
| --- | --- | --- | --- | --- | --- |
| Schlögl reactor | autocatalytic feed $a$ | $b = 1.49$ and $2.26$; $a_2 = \mp 1.73$ | $a = 4.153$, $b = 2.654$, $x = 1.384$ | $db/da = -1.917$ | supercritical pitchfork |
| Two tubes of unequal length | length $L_2$ | $\mu = 1.213$; $a_2 = 0.59$ | $L_2 = 1.000 = L_1$, $\mu = 1.383$ | along $\mu$ | subcritical pitchfork |
| Toggle with unequal promoters | promoter ratio $\gamma$ | $\alpha = 3.273$; $a_2 = -0.19$ | $\gamma = 1.000$, $\alpha = 2.000$ | along $\alpha$ | supercritical pitchfork |

Two results can be checked by hand. For the reactor, the drift $-x^3 + a x^2 - k_3 x + b$ has a triple root when
$x = a/3$, $k_3 = a^2/3$ and $b = a^3/27$; with $k_3 = 5.75$ this gives $a = 4.153$, $b = 2.654$ and
$x = 1.384$. Keeping the drift at $x = a/3$ unchanged to first order requires $db = -x^2\, da$, the reported
direction. For the toggle switch and the tubes, the cusp lies at the parameter value that restores the exchange
symmetry, $\gamma = 1$ and $L_2 = L_1$.

The flow network differs from the other two. With unequal lengths it writes one-sidedly: the flow selects one
tube irrespective of a weak bias. At equal lengths the write is symmetric but subcritical: beyond the write point at
$\mu = 1.383$ the state moves to a distant branch, and at the operating value $\mu = 2$ each stable state carries
the whole flow in one tube. A graded write would require a saturating term in the adaptation law.

## 5. Exercises

1. Why is the signature kept in the transfer from rods to dipoles although the ratio of the two couplings changes
   from about 1.25 to 3?
2. A second reactor has the drift $-x^3 + a x^2 - k_3 x + b$ with $k_3 = 12$. Where is its cusp?
3. The design for the toggle with unequal promoters returns $\gamma = 1$ to ten digits. Why is this value exact?

<details><summary>Answers</summary>

1. The signature depends on the kinds of bond (an alignment and a reflection across the bond axis), on the
   reciprocity of the couplings and on the role of the control; rods and dipoles share all three. The ratio of
   the couplings sets energies and barriers, which are not part of the signature.
2. $a = (3 k_3)^{1/2} = 6$, $b = a^3/27 = 8$, $x = a/3 = 2$.
3. At $\gamma = 1$ the exchange of the two genes is a symmetry of the equations. A symmetric write point then has
   $a_2 = 0$ by symmetry, not by a numerical coincidence.

</details>

## Summary

- A transfer keeps the memory signature determined by the structure of the operation and lists what differs on
  the new carrier; each difference is a prediction to test.
- A one-sided write (a fold) becomes symmetric at a cusp, which requires one additional tuned parameter.
- The design command locates the cusp and the direction of a sweep through it without a bias.

## Reference

| Result | Function | Test |
| --- | --- | --- |
| Transfer | [`cli.cmd_attach`](../../fieldbridge/memory/cli.py) | `test_attaching_the_colloid_memory_to_dipoles_keeps_its_class` |
| Cusp and sweep direction | [`construct.cancel_asymmetry`](../../fieldbridge/memory/construct.py), [`cli.cmd_design`](../../fieldbridge/memory/cli.py) | `test_schlogl_cusp_removes_the_one_sided_write`, `test_unequal_tubes_are_cancelled_only_by_equal_lengths`, `test_equal_promoters_remove_the_toggles_one_sided_write` |

[Previous: Module 4](18_memory_writing_and_retention.md) · [Next: Module 6, memory in the phase of an oscillator](20_memory_phase.md) · [Tutorial index](index.md) · [Glossary](memory_glossary.md)
