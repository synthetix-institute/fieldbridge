# Predictions fixed before the ensembles of entries 31 and 32 were run

Written on 31 October 2026, after the starting conditions and their energies existed and before any trajectory of the ensembles had run. The boundary comes from
the 16 trajectories in `results/displaced.json`, each started from a chosen torsion with the other modes at rest.

Claim under test: the torsion of the CHR end at TS1 has a boundary near -87 degrees in both reactions. A trajectory started on the far side of it (more negative) reaches
the Z radical and one started on the near side reaches the E radical, and the Z share follows from the thermal spread of the torsion across that boundary.

The torsion at the starts has mean and spread -109.8 and 29.8 degrees (entry 31) and -100.6 and 28.4 degrees (entry 32), with harmonic classical positions. Reweighting
1000 candidates by their DFT energies gives -106.7 and -95.8 degrees, with effective sample sizes of only 30 and 145.

| Quantity | Entry 31 (CH3) | Entry 32 (CH2F) |
|---|---|---|
| Boundary in the torsion at TS1 | -87 +- 10 degrees | -87 +- 10 degrees |
| Z share from the positions alone, harmonic | 0.83 | 0.71 |
| Z share from the positions alone, reweighted by energy | 0.80 | 0.65 |
| Estimate reported before these starts existed | 0.81-0.91 | 0.67-0.81 |
| Trajectory counts of Wu, Zhang and Hong (2024), abstract | 0.956 +- 0.018 | 0.823 +- 0.036 |

What would refute the claim:
1. the boundaries fitted in the two reactions differ by more than their errors, or lie outside the stated range;
2. the probability of Z does not change monotonically across the boundary, or the width of that change exceeds the spread of the torsion, so that the torsion does not decide;
3. the ensemble's share falls outside the range given by the boundary and spread, once the width of the change is included as P = Phi((a* - <a>) / sqrt(sigma^2 + w^2)).

These positions already sit below the published shares, so the claim and the published counts cannot both be right unless the velocity or the sampling changes the share.
