# Repair campaign on the law departures of the quantum-computing elements

`scripts/close_repair_loop.py docs/records/quantum_semiclassical.jsonl` on the cluster, 2026-10-30, under the rule of commit
d2e1c9a, which precedes every campaign run. Machine-readable results: `repair_loop_first_pass.json`, `repair_loop_second_pass.json`,
`repair_loop_certified_control.json` in this folder. The libraries with all attempts are on the cluster in
`~/Desktop/Hyperion/repair_loop_20261030/`.

## Two critical directions, 2026-11-02 (provisional)

Commit dd530a8 closes the loop at k real, semisimple critical directions. Before it the four failed laws of `kpo_triangle`
and `rfsquid_triangle` were refused. The run took the merged records, the cluster, fresh libraries and 3600 s per law, under
the rule of d2e1c9a. Results: `repair_loop_two_directions_campaign.json` and `repair_loop_two_directions_certified_control.json`.
The records still carry the one-direction outcome-map defect of 87427d6 (`pplo_oscillator`, `kerr_cat_kpo`, `kpo_pair`,
`phase_qubit_rfsquid`), so these counts are provisional until a rerun on the corrected records.

**7 of the 20 failed laws get a verified smaller description, and 0 of the 20 at the scale of the written states.** The 16
laws with one critical direction reproduce their earlier status, order, radius and frozen candidate, and so do the 18
such certified laws. The new verified law is `rfsquid_triangle` 0, six source coordinates replaced by two, at radius 0.0625 of
the separation of the two nearest written states. Its withheld error is 0.0026, the failed-law drift is wrong by 0.32, the power
control is rejected, and the plain critical coordinates miss by 0.039.

| Law (two critical directions) | Outcome |
|---|---|
| `rfsquid_triangle` 0 | verified smaller description, 6 -> 2 |
| `kpo_triangle` 0 | insufficient excitation (failed-law drift within 0.0095 at radius 0.0625) |
| `kpo_triangle` 1 | timed out at 3600 s |
| `kpo_triangle` 2 | no order passes development |

At the KPO triangle the reduced drift on the critical plane has quadratic terms of about 6.7. A displacement of 0.036 therefore
reaches the written states, about 1 away, within some 10 time units, and the truncated series stops at a spurious equilibrium
near 0.1. This is the reach of a local series, not a defect of the closure.

| Control | Result |
|---|---|
| exact identity | 0 of 20; the triangles have six source coordinates, beyond its four |
| certified laws | 9 of 20 (45 %, against 35 % on failed laws), median radius 0.125 against 0.0625; `cim_triangle` 0, the one certified law with two directions, verified at 0.125 |
| memory hint against the direct ladder | 636 against 707 development checks on 19 laws |

Tests at dd530a8: 933 passed, 7 skipped, 1 failed. The failure is `test_no_locale_dependent_file_reads`, because the encoding commit
of this branch was reverted (889e73a) in favour of the same fix on origin/main (0f86208). This branch does not merge it.

## Rerun on the merged records, 2026-11-02

Commit f9b6850, main at 0e19b29 merged in, run on the cluster with fresh libraries and a budget of 3600 s per law. Thresholds,
interventions and selection rule are those of d2e1c9a. Results: `repair_loop_merged_campaign.json`,
`repair_loop_merged_certified_control.json`.

**6 of the 20 failed laws (17 nonlinear-stage, 3 equilibrium-write) get a verified smaller description, and 0 of the 20 at the
scale of the written states.** The six are the five of the first campaign and `rfsquid_pair` 0, all nonlinear-stage. They hold at
radii 0.03 to 0.125 of the written-state separation, and each replaces two to four source coordinates by one. None of the three
equilibrium writes is verified. The 14 laws that already had a row get the same outcome as before.

| Outcome on the 20 | Count | Laws |
|---|---:|---|
| verified smaller description | 6 | `ion_duffing` 0, 1; `jba_resonator` 0, 1; `rydberg_bistability`; `rfsquid_pair` 0 |
| no order passes development | 5 | `cim_pair` 2; `phase_qubit_rfsquid` 0, 1; `ion_phonon_laser`; `optomech_bistable` |
| refused: a second critical eigenvalue | 4 | `kpo_triangle` 0, 1, 2; `rfsquid_triangle` 0 |
| no test power | 2 | `transmon_driven` 0, 1 |
| insufficient excitation | 2 | `pplo_oscillator` 0; `rfsquid_pair` 1 |
| withheld test failed | 1 | `magnon_kerr` |

| Control | Result |
|---|---|
| exact identity (search and symbolic, 120 s each) | 0 of 20 verified: budget exhausted, numerical failure, or non-polynomial refused |
| same closure on certified linear-stage laws | 8 of 19 verified (42 %, against 30 % on failed laws); 2 timed out (`kpo_pair` 1, 2) |
| median verified radius | 0.1875 on certified laws, 0.0625 on failed laws |
| plain critical coordinate within the tolerance | 2 of 6 (`jba_resonator` 0, 1: 0.004 and 0.006) |
| strongest nonlinear term scaled by 1.25 | rejected on all 6 (errors 0.013 to 0.060) |
| memory hint against the direct ladder | 591 against 602 development checks on 16 laws; one hint saves 22 and selects another pair, 11 cost one, 4 have none |

The closure does not separate failed from certified laws, and the margin of `rfsquid_pair` 0 is narrow, with a withheld error of
0.0095 against 0.01. The four refused laws lie at a write point with two critical directions, the limit already stated. Of the 230
retained outcomes, 20 get a task and 210 are `no_failed_law`. The campaign below, on the stored records, is kept as it was run.

## Result on the stored records, 2026-10-30

**5 of the 13 failed laws (9 nonlinear-stage, 4 equilibrium-write) get a verified smaller description, one coordinate in place
of two to four, and 0 of the 13 get one at the scale of the written states.** The five hold at 3 % to 6 % of the separation of
the written states. All five are nonlinear-stage laws; none of the four equilibrium writes is verified. The second pass did
not change the count.

This does not repair the failed laws. The test concerns the deterministic drift near the write point, with the noise of each
specification set to zero. The selection probabilities of the failed law were not re-derived, and the departure of a law arises
at the sweep's full excursion. There no order up to 7 passes development, and the selected order misses the coordinate on withheld
interventions by 0.6 to 1.2 of its amplitude, or diverges.

| Outcome on the 13 | First pass | After the second pass |
|---|---:|---:|
| verified smaller description | 5 | 5 |
| reaches the written states (radius 1) | 0 | 0 |
| no order passes development | 0 | 2 |
| falsely reported as resonant (defect below) | 2 | 0 |
| insufficient excitation (the failed-law drift is itself adequate at the radius) | 3 | 3 |
| no test power (scaling the strongest nonlinear term by 1.25 is not detected) | 2 | 2 |
| withheld test failed | 1 | 1 |

Each decision rests on three withheld interventions at the selected radius.

## Controls

**Exact identity on the same tasks.** The search of `learning.learn` and the direct-symbolic proposer of `repair_families.construct`
verify 0 of 13. Their budgets of 20 candidates end in `no_polynomial_closure`, with two exceptions: `rfsquid_pair`, whose drift is not
polynomial, is refused, and on `ion_phonon_laser` the symbolic path finds one candidate that leaves the numerical range during
development and is not tested further. None was verified, so the size rule had no exact extension to reclassify in this campaign.

**Laws that did not fail.** The same closure on the certified one-direction linear-stage laws verifies 9 of 20, against 5 of 13.
Being verified therefore does not separate failed from certified laws, and it cannot be read as an account of the failure. The
reach differs, on small numbers: the verified radii have a median of 0.125 of the written-state separation for the certified laws
(0.03 to 0.5) and 0.0625 for the failed ones (0.03 to 0.06). `cim_dopo` was also the pilot used to set the development margin, and
two `kpo_pair` laws ran out of the 900 s budget of this control. The control excludes the 5 certified laws with more than one
critical direction.

**The correction of the coordinate.** With the same reduced drift on the plain coordinate w.(x - x*), 2 of the 5 verified laws
(both of `jba_resonator`) also stay within the tolerance of 0.01, so the fibres were not needed there at that radius. The other 3
need them (plain error 0.016 to 0.044).

**The failed-law drift.** The linear-stage drift diverges over the 56.6 time units of the sweep on all five verified laws, so
the improvement requirement against it is met trivially there; it distinguishes little beyond the plain-coordinate control.

**Memory against the direct ladder.** On the 11 laws with a development table: 416 development checks without a hint and 348 with
an order and radius from an earlier verified law of the class. Two hints saved everything (35 to 1 and 40 to 1), five cost one
check more, and four laws had no hint. One saving comes from a sibling write point of the same model (`ion_duffing`), one from a
related oscillator (`jba_resonator` to `transmon_driven`), and the saving depends on the order of the campaign, which is the
sorted hash of each departure. This is a saving over the ladder of this closure, not over the exact symbolic baseline, which
closes nothing here. A hint below the radius reached counts as a failed check, because the ladder did not evaluate that pair.

## Table of the 13 (after the second pass)

Error, plain and altered are maximum errors at the selected radius, as a fraction of the coordinate's amplitude; the tolerance is
0.01. "Altered" is the power control.

| Model | Law | Class | Status | N | Radius | Error | Plain | Altered | Exact search |
|---|---|---|---|---:|---:|---:|---:|---:|---|
| ion_phonon_laser | 0/unfolding/0 | equilibrium write | no_order_passes_development | – | – | – | – | – | search_budget_exhausted |
| kerr_cat_kpo | 0/law | equilibrium write | insufficient_excitation | 3 | 0.03125 | 0.00053 | 0.00068 | 0.00053 | search_budget_exhausted |
| kpo_pair | 0/law | equilibrium write | insufficient_excitation | 2 | 0.03125 | 0.002 | 0.002 | 0.002 | search_budget_exhausted |
| rfsquid_pair | 1/law | equilibrium write | insufficient_excitation | 5 | 0.0625 | 0.0015 | 0.31 | 0.044 | unsupported_or_invalid |
| ion_duffing | 0/unfolding/0 | nonlinear stage | verified_smaller_description | 7 | 0.0625 | 0.0078 | 0.044 | 0.021 | search_budget_exhausted |
| ion_duffing | 1/unfolding/0 | nonlinear stage | verified_smaller_description | 7 | 0.0625 | 2.3e-06 | 0.018 | 0.026 | search_budget_exhausted |
| jba_resonator | 0/unfolding/0 | nonlinear stage | verified_smaller_description | 7 | 0.03125 | 2.6e-05 | 0.0045 | 0.013 | search_budget_exhausted |
| jba_resonator | 1/unfolding/0 | nonlinear stage | verified_smaller_description | 5 | 0.03125 | 0.00069 | 0.0061 | 0.013 | search_budget_exhausted |
| magnon_kerr | 0/unfolding/0 | nonlinear stage | withheld_failed | 5 | 0.03125 | 0.011 | 0.027 | 0.017 | search_budget_exhausted |
| optomech_bistable | 0/unfolding/0 | nonlinear stage | no_order_passes_development | – | – | – | – | – | search_budget_exhausted |
| rydberg_bistability | 0/unfolding/0 | nonlinear stage | verified_smaller_description | 6 | 0.0625 | 0.00069 | 0.016 | 0.025 | search_budget_exhausted |
| transmon_driven | 0/unfolding/0 | nonlinear stage | no_test_power | 7 | 0.03125 | 0.0013 | 0.029 | 0.0074 | search_budget_exhausted |
| transmon_driven | 1/unfolding/0 | nonlinear stage | no_test_power | 5 | 0.03125 | 0.00039 | 0.0092 | 0.0068 | search_budget_exhausted |

## Certified control

| Model | Law | Status | N | Radius | Error | Plain |
|---|---|---|---:|---:|---:|---:|
| cim_dopo | 0/law | verified_smaller_description | 7 | 0.25 | 0.004 | 0.06 |
| cim_pair | 0/law | verified_smaller_description | 3 | 0.25 | 0.00019 | 0.0011 |
| cim_pair | 1/law | no_order_passes_development | – | – | – | – |
| cim_pair | 2/law | no_order_passes_development | – | – | – | – |
| dissipative_cat | 0/law | verified_smaller_description | 3 | 0.0625 | 0.0021 | 0.0046 |
| flux_qubit_tilt | 0/unfolding/0 | insufficient_excitation | 4 | 0.03125 | 0.0008 | 0.014 |
| flux_qubit_tilt | 1/unfolding/0 | no_test_power | 7 | 0.0625 | 0.004 | 0.011 |
| fluxonium | 0/unfolding/0 | no_test_power | 5 | 0.03125 | 0.0036 | 0.022 |
| fluxonium | 1/unfolding/0 | verified_smaller_description | 5 | 0.03125 | 0.0076 | 0.019 |
| kpo_pair | 1/law | timed_out | – | – | – | – |
| kpo_pair | 2/law | timed_out | – | – | – | – |
| nuclear_bath | 0/unfolding/0 | withheld_failed | 7 | 0.5 | 0.02 | 0.02 |
| nuclear_bath | 1/unfolding/0 | verified_smaller_description | 7 | 0.5 | 0.00059 | 0.00059 |
| po_amplitude_oscillator | 0/law | verified_smaller_description | 7 | 0.125 | 0.0054 | 0.066 |
| pplo_oscillator | 0/law | insufficient_excitation | 2 | 0.125 | 9.8e-05 | 9.7e-05 |
| qfp_latch | 0/law | verified_smaller_description | 7 | 0.25 | 0.002 | 0.0056 |
| rfsquid_annealer | 0/law | verified_smaller_description | 7 | 0.0625 | 0.0051 | 0.0092 |
| rfsquid_pair | 0/law | verified_smaller_description | 7 | 0.125 | 0.0095 | 0.025 |
| rfsquid_pair | 2/law | withheld_failed | 5 | 0.0625 | 0.012 | 0.35 |
| snail_po | 0/law | insufficient_excitation | 3 | 0.125 | 0.0026 | 0.012 |

## The loop

`memory departures` on this file retains 227 outcomes; `repair_departures` returned `needs_equation_task` for all 227. Now every one
has a row: 13 receive a task and one campaign, and 214 are `no_failed_law`, with their retained kind (47 target mismatches, 42
reached targets, 5 numerical failures, 120 unclassified). No failure was discarded: each attempt stays in the library, and a rerun
keeps its earlier report.

## Rule

Fixed in d2e1c9a, on a synthetic system and a certified law, before any campaign run. A ladder of six radii (1 to 1/32 of the
written-state separation, with the parameter excursion scaling alike) carries three development and three withheld interventions
each, over the sweep's 56.6 time units, including preparations off the centre manifold. Selection on development takes the
largest radius at which an order 1 to 7 stays below half the tolerance, then the smallest such order. The closure must pass the
identity modulo order N before its withheld set is consumed, once, through `MoveLibrary.reserve_withheld`. Verified means: every
withheld error at that radius below 0.01, the failed-law drift wrong by at least 0.1 and improved on by at least 90 %, and the
altered coefficient detected. The thresholds were chosen, not derived.

## Defects found by this campaign

Corrected in bf37dff, with thresholds, interventions and selection unchanged. The first-pass result is retained as it was run.

- The resonance test read singular values, which a large bias makes tiny without a resonance. `ion_phonon_laser` (other rate -0.667)
  and `optomech_bistable` (-0.079 ± 2.068i, -0.941) were reported as resonant. Neither had reached its withheld test, so both were
  run once more, and neither passes development at any order or radius.
- The generated coordinate name `u` collided with a variable of the Rydberg model, so its exact baseline did not run. Rerun: no closure.
- The verifier scaled its limit to the first omitted coefficient, which is lenient where coefficients grow. The stricter one (Taylor
  coefficients of both sides at 60 digits, each order against its own size) passes all 11 selected candidates of the first pass.

The certified control ran on the corrected code. Its first run stalled on one law under the earlier code and was stopped.

## Limits

- **Euclidean carrier: lifted** in the series backend. A torus is treated in its covering space, an orthant at interior points; a
  boundary point is reported. Tested on synthetic torus and orthant checks, not in this campaign, whose 13 laws are all Euclidean.
- **Four source coordinates: lifted** in the series backend, to a budget of 4000 monomials; tested on six coordinates. None of the 13
  exceeds four, so this campaign does not exercise it. `learning._parse_task` keeps its limit of four, and the exact path still refuses more.
- **Non-polynomial drifts: lifted** in the series backend, which is what admits `rfsquid_pair` to the closure.
- **Remaining:** one critical direction with a real eigenvalue (a second critical eigenvalue or a complex pair is refused, which is
  why Hopf points and critical spaces of higher dimension are outside this rule); kinds other than equations; a resonance stops the series below it.
- **Scope:** deterministic drift only, as above. A verified description here makes no claim of novelty.
- **Dependency:** `fieldbridge/memory/learning.py` and `repair_families.py` were uncommitted in the main checkout. The campaign used
  a read-only copy (sha256 bba095b7…, 51700433…), not committed here, so `repair_loop.py` cannot be imported from this branch alone
  and its tests skip until those modules land.
- **Environment:** cluster Python 3.9, numpy 1.24, scipy 1.13, sympy 1.14; the laptop stays canonical for the law record, and no
  stored law value changes.
- `memory repair` is not wired to this loop, because `cli.py` also has uncommitted edits in the main checkout.
- The records file predates the fix of `cancel_asymmetry`; each cusp is recomputed from its write point. On the 13, none moved from
  the recorded cusp by more than 1e-3.
- **Tests:** the full suite on bf37dff gave 743 passed, 8 skipped, 1 failed in 3 h 23 min. The failure was
  `test_no_locale_dependent_file_reads`, on two reads without an encoding in `scripts/close_repair_loop.py`, corrected here. The full
  suite was not repeated after this correction, which touches only that script and `tests/test_repair_loop.py`; the locale test passes
  and `tests/test_repair_loop.py` skips without the two modules, and gives 31 passed with the read-only copy in place.
- **When the two modules land:** the locale test will fail on their own reads, at `learning.py:530`, `:535` and
  `repair_families.py:165`, `:182` of the copy, until these take an explicit encoding.

## The campaign set moves with the records

The campaign read the stored records of 2026-10-23, which predate commits before a06823f. The regenerated records on the series branch
(b4f8c15, not on main on 2026-10-31) give 20 failed laws in the two classes, 17 nonlinear-stage and 3 equilibrium-write, in place of
9 and 4. All five verified laws stay in the set. This is a comparison of the two files and the rows already run, not a new campaign.

| Law | Stored records | Series records | Row already run, at the frozen rule |
|---|---|---|---|
| `kerr_cat_kpo` 0, `kpo_pair` 0 | failed, equilibrium write | certified, linear stage | insufficient excitation, both |
| `optomech_bistable` 0 | failed, nonlinear stage | failed, equilibrium write | no order passes development |
| `rfsquid_pair` 0 | certified, linear stage | failed, nonlinear stage | verified at radius 0.125 |
| `pplo_oscillator` 0 | certified, linear stage | failed, nonlinear stage | insufficient excitation |
| `cim_pair` 2 | certified, linear stage | failed, nonlinear stage | no order passes development |
| `phase_qubit_rfsquid` 0, 1 | certified, linear stage | failed, nonlinear stage | none |
| `rfsquid_triangle` 0 | certified, linear stage | failed, nonlinear stage | none |
| `kpo_triangle` 0, 1, 2 | not derived | failed, nonlinear stage | none |

On the rows already run, 6 of the 14 laws of the new set that have a row are verified, and 6 laws have no row. The certified control
loses `rfsquid_pair` 0, `pplo_oscillator` 0 and `cim_pair` 2. The count of 5 of 13 therefore holds for the stored records only, and a
count on the new set needs the series branch on main, followed by a campaign with fresh libraries. That campaign is the rerun above.
- **Tests on the merged state (f9b6850):** the full suite gave 928 passed, 7 skipped, 2 failed in 2 h 34 min, at a laptop
  load average near 290. `test_no_locale_dependent_file_reads` fails on main itself, on 63 reads and writes in 16 files that
  came with the learning modules, none of them in the repair loop. `test_a_child_that_outruns_its_time_limit_keeps_the_steps_it_had_written`
  found the symmetric write not derived within its time limit, and it passes alone (38.6 s, load 35). This branch does not
  change that test or the code it exercises.
