# Rule R8 against the published hysteron platforms

The rule: a network of hysterons under a scalar drive returns to every turning point when the signed network of
couplings with the drive as a node is balanced; the frustration index F counts the edges that must change sign. The
consequence measured in the sequence-memory study on random networks: the number of distinguishable states exceeds
that of the ferromagnetic reference by about 0.11 to 0.17 bit per unsatisfiable edge above a threshold, almost all of
the added states need a drive with three or more reversals, and the first switching fields of a descending sweep
identify them. The table lists what each platform's coupling structure implies and what its paper measured, from the
abstracts and the full-text summaries of the prior-art record (2026-10-11); the sign structure of each platform should
be confirmed from the paper before any claim is made.

| Platform | Coupling signs as reported | Drive | F with the drive | Measured in the paper | Not measured: the prediction to test |
| --- | --- | --- | --- | --- | --- |
| Serially coupled bistable beams, Liu et al., PNAS 121, e2308414121 (2024); n <= 3 experiments, <= 5 enumerated | serial coupling: negative (an element's switch hinders its neighbour) | common compression | > 0 for every negative coupling (unsatisfiable triangle with the drive) | pathways that violate return-point memory, sequential processing | the count of distinguishable states against F for chains of 8 to 30 beams, its linear growth and threshold, the share of states needing three reversals |
| Hysteretic springs in series and parallel, Shohat and van Hecke, PRL 134, 188201 (2025); <= 3 springs | series: negative; parallel: positive | common displacement | > 0 with any series coupling; 0 for all-parallel | geometric control of the pathways | the gain over the all-parallel reference for the same springs, the return test on mixed networks |
| Hysteron latch, Lindeman, Jalowiec and Keim, Sci. Adv. 11, eadr5933 (2025); pairs | a frustrated pair (one negative coupling) | one-sided drive from a reference state | 1 per latch | 2^M sequences stored by M latches under descending sequences | the identification of the stored sequence from the switching fields of one sweep; the capacity of randomly coupled latches against F |
| Mechanical counters, Kwakernaak and van Hecke, PRL 130, 268204 (2023); up to 11 cells | contacts between neighbouring beams, designed | cyclic compression | > 0 (counting requires return-point memory to fail) | counting of cycles, order of three magnitudes | the number of states reachable by non-periodic histories against F, which the counter design does not use |
| Frustrated square lattice, Sirote-Katz et al., Nat. Commun. 15, 4008 (2024) | antiferromagnetic couplings on an even lattice with stabilizing bars | local flips and a global drive | > 0 under a uniform drive (every negative coupling frustrates with the drive); 0 under a staggered drive after relabeling | order-dependent final states on one lattice | the gain under a uniform against a staggered drive on the same lattice: the rule predicts zero gain and exact return under the staggered drive |
| Compressed corrugated sheets, Bense and van Hecke, PNAS 118, e2111436118 (2021) | positive and negative, from the sheet's geometry | common compression | > 0 | complex pathways and memory | the count of pathways against the number of unsatisfiable contacts |
| Soft spots of sheared amorphous solids, Keim and Paulsen, Sci. Adv. 7, eabg7685 (2021); Mungan et al., PRL 123, 178002 (2019) | both signs, from elasticity | cyclic shear | > 0 | multiperiodic orbits, transition graphs | the capacity of a sample against the frustration of its soft-spot network, which requires the sign of every interaction |
| One-dimensional random antiferromagnet, Deutsch, Dhar and Narayan, PRL 92, 227203 (2004) | all negative, chain | uniform field from saturation | N/2 (relabel one sublattice: a ferromagnetic chain with a staggered drive) | return-point memory holds in their model | the exception to the converse, checked here on chains of 12 elements with random fields and half-widths (`counterpart.construct`, seeds 0 to 3): three of four chains return in every excursion with F = 6, one fails 230 of 1,750; all four hold 0.9 to 2.8 bits more than their ferromagnetic reference with no deep memory. The gain of a chain is the absence of locking, since an antiferromagnetic chain has no avalanches; the return test, not the gain, decides whether the memory is beyond return-point. The rule is sufficient, not necessary. |
| Three-dimensional spin glass at T = 0, Deutsch and Narayan, PRL 91, 200601 (2003) | random signs | uniform field | extensive | subharmonic and aperiodic loops; the response period grows with size | the capacity per spin against the frustration per spin |

Reading of the table. Every platform with a negative coupling under a common drive has F > 0 and is predicted to hold
states beyond return-point memory; in each case the papers report the qualitative fact and none reports the count of
states against F, the threshold, the share of deep states, or the identification from a partial sweep. Those are the
measurements the rule asks for. The staggered-drive prediction for the frustrated lattice and the chain exception are
the two cases where the rule can fail or be sharpened, and they cost one enumeration each.
