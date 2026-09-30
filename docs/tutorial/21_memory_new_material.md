# Module 7. A material from your own field

**Learning objectives.** After this module you can

1. describe a material from your field as a realization, naming each component in the terms of that field;
2. apply the analyses of Modules 2–6 to it in a fixed order and keep a research record;
3. read the gallery of realizations and compare memory mechanisms across fields;
4. state what the calculations establish and what they do not.

**Prerequisites.** Modules 1–6. **Time.** About 60 minutes, and the time to write your own model.

## 1. The components in the terms of six fields

| Component | Colloids at an interface | Magnetic dipole arrays | Gene circuits | Chemical reactors | Flow networks | Chemical oscillators |
| --- | --- | --- | --- | --- | --- | --- |
| Carrier $\Xi$ | particle orientations (period $\pi$) | moment directions (period $2\pi$) | protein concentrations | concentrations | tube conductances | concentrations along a cycle |
| Operation $\Omega$ | capillary and alignment couplings | dipolar energy | activation and repression | mass-action kinetics | adaptation of the conductance to the flux | rate laws |
| Control | coupling strength relative to $k_BT$ | moment or temperature | promoter strength | feed rate | reinforcement exponent | a rate constant |
| Writing protocol $P$ | a uniform or rotating field | a magnetic field | a pulse of inducer | a pulse of feed | a transient pressure difference | a pulse at a chosen phase |
| Observable $R$ | overlap with a written pattern | local moments | reporter fluorescence | concentration | flow in each tube | timing of the oscillation |
| Closure $C$ | rotational Brownian motion | thermal fluctuations | intrinsic noise of expression | concentration fluctuations | fluctuations of the flow | phase noise |

## 2. The sequence of analyses

1. Write the specification (Module 2), starting from `memory new NAME --carrier euclid|orthant|torus`, which
   writes a template that loads and runs. List every uncertain modelling choice in `assumptions`, and check the
   file with `memory check`, which reports the fields still left as placeholders.
2. Run `memory predict` and record the structural predictions before any calculation (Module 3).
3. Run `memory card`: stable states, write points, retention law and writing protocols, and the agreement with the
   predictions (Module 4).
4. Run `memory attach --from <a known memory>` to see which properties of a known mechanism your carrier keeps
   (Module 5).
5. If the write is one-sided, run `memory design --free <parameter>` to find the setting that makes it symmetric
   (Module 5).
6. If the realization oscillates, run `memory phase` (Module 6); if it is a continuous field, run `memory field`
   (Module 8).
7. Perform a control calculation: change one structural feature that must change the answer, such as a symmetry,
   the sign of a loop or only the mobility, and confirm that the answer changes as predicted.
8. To add the material to the repository, regenerate the [catalog of materials](../materials.md) with
   `memory catalog` and open a pull request ([CONTRIBUTING.md](../../CONTRIBUTING.md)).

## 3. The gallery of realizations

```bash
python3 -B -m fieldbridge memory gallery --quick --out-dir build/gallery
```

The command computes one card for each example in `examples/memory` and repeats the loop tests of Module 3. Open
`build/gallery/index.html` in a browser; `gallery.md` lists the same results as text.

| Realization | Writing mechanism | Stable states | Retention |
| --- | --- | --- | --- |
| caged capillary rotors | write by a uniform field | 4 | Law 3, activation |
| three compartments with an unobserved imbalance | single state: relaxation | 1 | Law 1, relaxation |
| caged in-plane dipoles | write by a uniform field | 10 (found by sampling) | Law 3, activation |
| pitchfork normal form | symmetric write (supercritical pitchfork) | 2 | Law 3, activation |
| ring of 3 repressors | limit cycle: phase memory | none | Law 2 along the cycle |
| ring of 4 repressors | symmetric write (supercritical pitchfork) | 2 | Law 3, activation |
| Schlögl reactor | one-sided write (fold); symmetric at a cusp | 2 | Law 3, activation |
| genetic toggle switch | symmetric write (supercritical pitchfork) | 2 | Law 3, activation |
| toggle switch, unequal promoters | one-sided write (fold); symmetric at a cusp | 2 | Law 3, activation |
| adaptive network, two equal tubes | symmetric write to a distant state (subcritical pitchfork) | 2 | Law 3, activation |
| adaptive network, unequal tubes | one-sided write (fold); symmetric at a cusp | 2 | Law 3, activation |

The retention laws are defined in [Module 4, Section 1.3](18_memory_writing_and_retention.md#13-retention). The classes describe the models; none of them is a failed calculation. A fold, a subcritical write or a single
state is a property of the material as modelled, and Modules 4–6 show what each implies for writing and
retention.

## 4. A research record

Keep one record per material. The example below collects the results of Modules 4, 5 and 7 for the flow
network.

| Entry | Flow network with two tubes |
| --- | --- |
| Question | Can two tubes retain which of them carried the flow? |
| Structural prediction (`memory predict`) | The exchange of the two tubes is a symmetry: if the state of equal flows loses stability along the mode that the exchange reverses, it does so at a pitchfork (a symmetric write); a positive loop allows two states |
| Calculation | Subcritical pitchfork at $\mu = 1.38$; two stable states, each carrying the whole flow in one tube; Law 3 |
| Control calculation | Unequal lengths ($L_2 = 1.15$): a fold at $\mu = 1.21$; the write becomes one-sided |
| Design | Equal lengths restore the symmetric write; it remains subcritical |
| Open question | Does a real network have the exchange symmetry, and at which reinforcement exponent? |

## 5. What the calculations establish

Every report carries `novelty_established: false` and states its evidence boundary. A structural prediction and a
calculation describe the supplied model. Whether the model describes a material is decided by experiment, and
whether a result is new is decided by the literature. Many of the rules used here are established in their own
fields: loop signs in gene regulation (Thomas), frustration in spin systems (Toulouse), symmetry in bifurcation
theory, swept bifurcations (Kondepudi). The program applies them together, to any carrier, and names the rule
behind each prediction.

## 6. Exercises

1. Choose a material from your field and fill the table of Section 1 for it.
2. Your card predicts a symmetric write. Which control calculation tests this prediction?
3. A card reports "write by a uniform field" and Law 3. Which protocol writes a state with a field below the
   threshold?

<details><summary>Answers</summary>

1. The answer depends on the material; every entry should be a measurable quantity or a stated assumption.
2. Break the symmetry on which the prediction rests, for example by making two couplings unequal. The write must
   then become a fold, as for the toggle switch in Module 1.
3. A protocol that changes the landscape: the field acts while the barriers are low, and the barriers are raised
   after the write, for example by raising the coupling or lowering the temperature (Module 4).

</details>

## Summary

- A material from any field is described by the same six components.
- The analyses are applied in a fixed order, and every structural prediction is tested by a control calculation.
- The gallery compares the mechanisms of eleven realizations; the classes are descriptive.
- Calculations describe models; experiments and the literature decide whether a model describes a material and
  whether a result is new.

## Reference

Code: [`cli.py`](../../fieldbridge/memory/cli.py) (`cmd_gallery`), [`discovery.py`](../../fieldbridge/memory/discovery.py)
(`evaluate`, `verdict`, `mechanism_class`), [`visual.py`](../../fieldbridge/memory/visual.py) (`html_report`).
Test: `test_command_line_writes_reports_with_provenance`.

[Previous: Module 6](20_memory_phase.md) · [Next: Module 8, time](22_memory_time.md) · [Tutorial index](index.md) · [Glossary](memory_glossary.md)
