"""The realizations shown on the web page, the single-component edits that connect them, and the stepped sequences.

A node is a realization I_real = ((Omega, Xi); C, R, P; A): a specification file, or a specification derived from
another node by one change (a parameter value, a term, the observable, an attachment to another carrier). An edge is
one change of one component. site_data.py computes every node with the FieldBridge functions, checks that every
edge changes only the component it names, and fills the numbers of the texts from those calculations.

Texts are HTML. A placeholder {key} is a fact of the node the edge leads to (the node itself for a first step);
{from.key} is a fact of the node it starts from. No number is written here by hand.
"""
from __future__ import annotations

from math import pi

# ------------------------------------------------------------------------------------------------ the language
SLOTS = {
    "Omega": {"symbol": "Ω", "name": "operation",
              "definition": "the generator of the dynamics: a Hamiltonian, or a drift"},
    "Xi": {"symbol": "Ξ", "name": "carrier",
           "definition": "the states and the operators or variables that act on them"},
    "C": {"symbol": "C", "name": "closure",
          "definition": "what is specified or discarded to obtain closed equations: constitutive relations, "
                        "admissible states and operator domains, boundaries, imposed conservation laws, "
                        "eliminated degrees of freedom"},
    "R": {"symbol": "R", "name": "observable", "definition": "what is measured"},
    "P": {"symbol": "P", "name": "protocol", "definition": "how the system is prepared and driven"},
    "A": {"symbol": "A", "name": "parameters",
          "definition": "the material or apparatus that implements the model, entered as parameter values"},
}

# mechanism classes: the columns of the map, in reading order
CLASSES = {
    "rotation": "Bloch rotation",
    "conserved": "conserved observable",
    "obstructed": "several frequencies",
    "single-state": "single stable state",
    "symmetric-write": "symmetric write (supercritical pitchfork)",
    "threshold-write": "one-sided write (fold)",
    "subcritical-write": "write to a distant state (subcritical pitchfork)",
    "field-write": "write by a uniform field",
    "return-point": "return-point memory",
    "no-return": "return not exact",
    "perfect-adaptation": "perfect adaptation",
    "fine-tuned-adaptation": "fine-tuned adaptation",
    "partial-adaptation": "partial adaptation",
    "no-adaptation": "no adaptation",
    "linear-memory": "linear memory",
    "odd-capacity": "odd degrees only",
    "nonlinear-capacity": "nonlinear capacity",
    "integrating": "no fading memory",
    "inherited-through-threshold": "inherited through a threshold",
    "kept-above-threshold": "kept above the threshold",
    "lost-in-the-dip": "lost in the dip",
    "threshold-moved": "threshold moved by division",
    "synchronizes": "synchronizes above a critical coupling",
    "no-onset": "no onset of synchrony",
    "follows-the-bias": "follows the bias (swept write)",
    "set-by-the-sample": "set by the sample",
    "reflection-seed": "passage from a one-component seed",
    "rotation-seed": "passage from a two-component seed",
    "linear-stage-write": "linear-stage write (vacuum seed)",
    "equilibrium-write": "equilibrium write (balance of the wells)",
    "oscillation": "limit cycle",
    "neutral-cycles": "neutral cycles",
    "exponential-loss": "exponential loss",
    "power-loss": "power-law loss",
    "convention": "stochastic calculus convention",
}

# the mechanism of each class written without a field: its canonical form, what it does, and the realization that the
# map of mechanisms opens for it (a canonical form where the page has one)
MECHANISMS = {
    # three operators rotate under the generator, i[H, J_a] = eps_abc Omega_b J_c: they generate su(2) together with
    # H, or they are the closure of the observable within a larger algebra
    "rotation": {"canonical": "<b>ṁ</b> = <b>Ω</b> × <b>m</b>",
                 "text": "A three-vector of expectations rotates about a fixed axis, and the measured signal follows "
                         "the Rabi law.", "node": "rotation_canonical"},
    "conserved": {"canonical": "[H, R] = 0",
                  "text": "The observable commutes with the generator: its measured value does not change.",
                  "node": "rotation_axis"},
    # the closure of the observable, span{R, [H, R], [H, [H, R]], ...}, decides the class: one operator is a
    # conserved observable, two or three move with one frequency, more move with several
    "obstructed": {"canonical": "closure of R larger than three operators",
                   "text": "The commutators of the generator with the observable do not close on three operators: "
                           "the observable moves with several frequencies, and the measured signal is not that of "
                           "one rotation.", "node": "spin1_easy_axis"},
    "single-state": {"canonical": "ẋ = −κx, κ &gt; 0",
                     "text": "Every preparation relaxes to one state, and nothing of the preparation is kept.",
                     "node": "pitchfork_below"},
    "symmetric-write": {"canonical": "ẋ = εx − x<sup>3</sup> + h",
                        "text": "Two stable states appear together at a supercritical pitchfork; a weak bias during "
                                "the crossing selects the state that is written.", "node": "pitchfork"},
    "threshold-write": {"canonical": "ẋ = μ + x<sup>2</sup>",
                        "text": "A stored state disappears at a fold; a field past the threshold switches the state "
                                "after a delay set by the Airy law.", "node": "pitchfork_bias"},
    "subcritical-write": {"canonical": "ẋ = εx + ax<sup>3</sup> − x<sup>5</sup>, a &gt; 0",
                          "text": "The state loses stability at a subcritical pitchfork and jumps to a distant state.",
                          "node": "pitchfork_subcritical"},
    "field-write": {"canonical": "g U(q) − h·m(q)",
                    "text": "The control only rescales the energy; a uniform field writes a state, and the barriers "
                            "between states keep it.", "node": "colloid_patch"},
    # interacting hysterons under a slow drive: the drive counted as an element of the network; without a frustrated
    # loop through it a relabeling gives cooperative couplings and a uniform drive (no passing)
    "return-point": {"canonical": "f<sub>i</sub> = Σ<sub>j</sub> K<sub>ij</sub>σ<sub>j</sub> + h<sub>i</sub> + H, "
                                  "K<sub>ij</sub> ≥ 0",
                     "text": "Elements switch at thresholds of a slow drive. With cooperative couplings and a drive that "
                             "pushes every element the same way, the state at a turning point is recovered exactly "
                             "after any excursion inside it.", "node": "rfim_ferromagnet"},
    "no-return": {"canonical": "a subloop that does not return; it needs η<sub>i</sub>J<sub>ij</sub>η<sub>j</sub> "
                               "&lt; 0 on some coupling",
                  "text": "After an excursion the state at a turning point differs in a measured subloop. This needs a "
                          "frustrated loop through the drive, which allows a failure without forcing one.",
                  "node": "rfim_antiferromagnet"},
    # regulation: the return of an output to its set point after a step of an input; an integrator of the error with a
    # stable steady state makes the return exact at every parameter value (the internal model principle)
    "perfect-adaptation": {"canonical": "dφ/dt = g(q) (y − y<sub>0</sub>), g of one sign",
                           "text": "An integrator of the error and a stable steady state: after a step of the input the "
                                   "output returns exactly to its set point, whatever the values of the rates.",
                           "node": "reg_pi_loop"},
    "fine-tuned-adaptation": {"canonical": "G = 0 only where k<sub>1</sub>k<sub>4</sub> = k<sub>2</sub>k<sub>3</sub>",
                              "text": "The output returns exactly only at tuned parameters: the integrator exists only "
                                      "on the surface where two paths of the input cancel.",
                              "node": "reg_feedforward_subtractive"},
    "partial-adaptation": {"canonical": "dφ/dt = k (y − y<sub>0</sub>) − δφ: a fraction 1/(1 + kg/δ) remains",
                           "text": "A leaky integrator: after a step the output returns part of the way, and the "
                                   "fraction that remains falls as the leak falls.", "node": "reg_leaky_integrator"},
    "no-adaptation": {"canonical": "no integrator: G = dy/du ≠ 0 and no return after the peak",
                      "text": "The output moves to its new value and stays there; feedback that is only proportional "
                              "reduces the change without returning it.", "node": "cruise_control_p"},
    # computation: which functions of the input history a linear combination of the measured observables of a driven
    # body represents (the information processing capacity of Dambre et al. 2012), decided by its equations
    "linear-memory": {"canonical": "C(k) = (1 − a²) a<sup>2k</sup>; Σ<sub>k</sub> C(k) = n<sub>lin</sub>",
                      "text": "A linear body with linear observables represents past inputs, not their products: the "
                              "capacity lies at degree 1, falls with the delay as the modes decay, and sums to the "
                              "rank of the linear response.", "node": "comp_one_mode"},
    "odd-capacity": {"canonical": "F(−x, −u) = −F(x, u): C = 0 at every even degree",
                     "text": "A body odd about its steady state, measured by odd observables and driven by a symmetric "
                             "input, represents products of an odd number of past inputs only.",
                     "node": "comp_odd_oscillator"},
    "nonlinear-capacity": {"canonical": "C > 0 at even degrees; Σ<sub>degrees, delays</sub> C = number of independent "
                                        "measured signals",
                           "text": "Nonlinear terms move capacity to products of past inputs; with fading memory the "
                                   "capacities of all degrees together equal the number of independent measured "
                                   "signals.", "node": "comp_square_cascade"},
    "integrating": {"canonical": "κ = 0: x<sub>t</sub> = Σ<sub>k</sub> u<sub>t−k</sub>",
                    "text": "A mode without decay keeps the running sum of the input; its correlation with the input at "
                            "any fixed delay vanishes as the sum grows.", "node": "comp_running_sum"},
    # heredity: a body whose order exists only above a critical size divides below it and regrows through the pitchfork
    "inherited-through-threshold": {"canonical": "P = Φ(φ<sub>c</sub>/σ<sub>c</sub>), σ<sub>c</sub><sup>2</sup> = "
                                                 "2D<sub>s</sub>(π/ar)<sup>1/2</sup>Φ(aμ<sub>0</sub>(2/ar)<sup>1/2</sup>)",
                                    "text": "A daughter born below the threshold loses part of its order in the dip; "
                                            "growth carries it back through the pitchfork, where it keeps its "
                                            "parent's sign when the order that remains outweighs the noise of the "
                                            "crossing.", "node": "her_normal_form"},
    "kept-above-threshold": {"canonical": "L<sub>div</sub>/2 &gt; L<sub>c</sub>: no dip",
                             "text": "The daughter is born above the threshold, and the order passes to it, kept "
                                     "behind its barrier.", "node": "her_normal_form_no_dip"},
    "lost-in-the-dip": {"canonical": "ln G = ∫λ dt over a generation &lt; 0",
                        "text": "The dip removes more order than the regrowth restores: without noise the order dies "
                                "out over the generations.", "node": "her_normal_form_lost"},
    "threshold-moved": {"canonical": "the halves differ in a conserved amount that sets L<sub>c</sub>",
                        "text": "The order is a redistribution of a conserved amount: division gives the daughters "
                                "different amounts, which moves their thresholds, and the parent's order is not passed "
                                "on.", "node": "polarity_brauns2020"},
    # decision: a population acts as one; it synchronizes, follows a bias through its collective threshold, or leaves an
    # unstable state from the seed of its own fluctuations
    "synchronizes": {"canonical": "2/K<sub>c</sub> = (b<sub>1</sub> − i a<sub>1</sub>) [π g(Ω) − i PV∫g(ω)/(ω − Ω) dω]",
                     "text": "Limit-cycle units with a spread of frequencies, coupled through an observable. Reduced to "
                             "phases, they lose incoherence at the coupling K<sub>c</sub> set by the density g(ω) and the "
                             "first harmonic of the coupling function H, and synchrony grows at the rate of the unstable "
                             "root.", "node": "dec_van_der_pol"},
    "no-onset": {"canonical": "b<sub>1</sub> = 0: the dispersion relation has no root inside the band",
                 "text": "The first harmonic of the coupling function has no sine part: the coupling only shifts the "
                         "frequencies, and no coupling makes the incoherent state unstable.", "node": "dec_josephson"},
    "follows-the-bias": {"canonical": "P = Φ(π<sup>1/4</sup> h<sub>s</sub> h / (D<sub>s</sub><sup>1/2</sup>"
                                      "(ar)<sup>1/4</sup>)), D<sub>s</sub> ∝ 1/N",
                         "text": "A population swept through its collective pitchfork selects the state favoured by a "
                                 "weak bias with a probability set by the bias, the sweep rate and the noise of the "
                                 "collective mode, which falls as 1/N.", "node": "dec_ising"},
    "set-by-the-sample": {"canonical": "P = Φ(h / (σ<sub>th</sub><sup>2</sup> + s<sub>q</sub><sup>2</sup>)<sup>1/2</sup>),"
                                       " s<sub>q</sub> ∝ N<sup>−1/2</sup>",
                          "text": "A finite sample of diverse units is not symmetric: its mean carries a frozen bias. "
                                  "When that bias exceeds the thermal spread, the outcome is set by the composition of "
                                  "the sample, and a slower sweep does not average it out.",
                          "node": "dec_diverse_random"},
    "reflection-seed": {"canonical": "Λ(τ<sub>90</sub>) − Λ(τ<sub>10</sub>) = ln 13.09",
                        "text": "After a step one real eigenvalue leads. The unstable mode grows from a one-component "
                                "Gaussian seed, and the 10–90% window of passage times is ln 13.09 in the growth "
                                "exponent Λ = ∫λ dt, whatever the noise; the median moves by 1/2 per factor e in N.",
                        "node": "dec_cim_step"},
    "rotation-seed": {"canonical": "Λ(τ<sub>90</sub>) − Λ(τ<sub>10</sub>) = ln 4.675",
                      "text": "A complex pair leads: the mode rotates while it grows from a two-component seed, and the "
                              "window of passage times is ln 4.675 in the growth exponent.", "node": "dec_macrospin"},
    # the quantum write: a parametric oscillator swept through its threshold with a weak bias chooses one of its two
    # states; in the linear stage the noise of the seed is fixed by the loss and the temperature (the vacuum)
    "linear-stage-write": {"canonical": "ẋ = (ε<sub>2</sub>(t) − κ/2) x + h + √(2D) ξ, 2D = κ(2n̄ + 1)/4",
                           "text": "A parametric oscillator swept through its threshold with a weak bias chooses one of "
                                   "its two states in the linear stage: the probability of the favoured state is "
                                   "Φ(h I<sub>1</sub>/√(σ<sub>0</sub><sup>2</sup> + 2D I<sub>2</sub>)), the classical write "
                                   "law with the noise fixed by the loss and the temperature and no free parameter.",
                           "node": "kpo_27"},
    "equilibrium-write": {"canonical": "dp/dt = Γ(ε<sub>2</sub>) [P<sub>eq</sub>(ε<sub>2</sub>) − p]",
                          "text": "At a few stored photons the two wells are shallow and switch faster than the sweep "
                                  "passes: the probability relaxes to the selection of the biased steady state at the "
                                  "rate of the Lindbladian's gap, below the linear-stage law.",
                          "node": "kpo_few"},
    "oscillation": {"canonical": "driven: φ̇ = ν − K sin φ",
                    "text": "The preparations settle on a limit cycle whose phase is neutral; a periodic drive locks "
                            "the phase inside the Adler range.", "node": "van_der_pol"},
    "neutral-cycles": {"canonical": "dI/dt = 0",
                       "text": "A conserved quantity fills the plane with closed orbits, none of which attracts its "
                               "neighbours.", "node": "lotka_volterra"},
    "exponential-loss": {"canonical": "SNR ∝ e<sup>−2Mκ<sub>0</sub>t</sup>",
                         "text": "Each mode of the field relaxes at a finite rate κ(k) (Law 1): the trace of a write is "
                                 "lost exponentially.",
                         "node": "field_nonconserved"},
    "power-loss": {"canonical": "SNR ∝ t<sup>−d/2−n</sup>",
                   "text": "Each mode of the field relaxes at its own rate (Law 1), and for a conserved density the "
                           "rate κ(k) = Mk² vanishes at long wavelengths: the trace of a write decays as a power of "
                           "time.", "node": "field_charge_1d"},
    "convention": {"canonical": "Itô μ − σ<sup>2</sup>/2, Stratonovich μ",
                   "text": "The reading of the noise term is part of the closure: it changes the measured growth rate "
                           "of log X by σ<sup>2</sup>/2.", "node": "log_ito"},
}
# the realization the page opens on: a mechanism written without a field
START = "pitchfork"

# specifications in these folders that no node names are added to the page automatically (site_data.auto_nodes):
# a contributed material appears in the column of its mechanism, joined to a realization of the same mechanism
AUTO_DIRS = {"examples/quantum": "unitary", "examples/memory": "dissipative", "examples/memory/oscillators": "dissipative",
             "examples/memory/fields": "field", "examples/memory/hysterons": "hysterons",
             "examples/regulation": "regulation", "examples/computation": "computation", "examples/quantum/open": "open",
             "examples/heredity": "heredity", "examples/decision": "decision"}

Q_CH24 = "docs/tutorial/24_spin_language.md"
Q_CH11 = "docs/tutorial/11_quantum_closure.md"
M9 = "docs/tutorial/23_memory_codiscovery.md"
M10 = "docs/tutorial/25_memory_threshold_write.md"
M11 = "docs/tutorial/26_memory_phase_locking.md"
M12 = "docs/tutorial/27_memory_return_point.md"
R1 = "docs/tutorial/28_regulation_set_point.md"
C1 = "docs/tutorial/29_computation_capacity.md"
QW = "docs/tutorial/30_quantum_write.md"
H1 = "docs/tutorial/31_heredity_threshold.md"
D1 = "docs/tutorial/32_decision_population.md"

# ------------------------------------------------------------------------------------------------ nodes
# family: unitary (fieldbridge.quantum), dissipative (fieldbridge.memory, equations and networks), field
# (fieldbridge.memory.fields), hysterons (fieldbridge.memory.hysterons), stochastic (fieldbridge.verification). A node
# without "spec" is built from "base".
NODES = [
    # -- the spins of Chapter 11 and their edits
    {"id": "two_spins", "family": "unitary", "spec": "examples/quantum/two_spins.json",
     "name": "two coupled spins (Chapter 11)", "tutorial": Q_CH24 + "#3-detaching-the-rotation-from-the-spins-of-chapter-11"},
    {"id": "two_spins_h0", "base": "two_spins", "params": {"h": 0.0},
     "name": "two coupled spins without a transverse field (Chapter 11)", "tutorial": Q_CH11,
     "preparations": [
         {"id": "top", "label": "top eigenstate of X<sub>0</sub>"},
         {"id": "plus_y", "label": "first spin along +y, second along +z", "product": ["+y", "+z"]},
         {"id": "minus_y", "label": "first spin along −y, second along +z", "product": ["-y", "+z"]}]},
    {"id": "two_spins_z0", "base": "two_spins_h0", "observable": [{"coefficient": 1, "operator": "Z0"}],
     "name": "two coupled spins measured through Z<sub>0</sub>", "parent": "two_spins_h0", "tutorial": Q_CH11},
    {"id": "two_spins_x1", "base": "two_spins", "replace_terms": {1: {"coefficient": "h", "operator": "X1"}},
     "name": "two coupled spins with the field on the second spin", "parent": "two_spins",
     "tutorial": Q_CH24 + "#6-a-larger-algebra-the-closure-of-the-observable"},
    {"id": "two_spins_both", "base": "two_spins", "add_terms": [{"coefficient": "h", "operator": "X1"}],
     "name": "two coupled spins with the field on both spins", "parent": "two_spins",
     "tutorial": Q_CH24 + "#6-a-larger-algebra-the-closure-of-the-observable"},
    # -- the bridge to the magnet
    {"id": "spin1_transverse", "attach": ("two_spins_h0", "spin", 2),
     "name": "spin 1 in a transverse field", "tutorial": Q_CH24 + "#4-attaching-the-rotation-to-the-exchange-chain-of-chapter-14"},
    {"id": "spin1_easy_axis", "base": "spin1_transverse", "params": {"D": -2.0},
     "add_terms": [{"coefficient": "D", "operator": "Jz Jz"}],
     "name": "spin 1 with easy-axis anisotropy", "parent": "spin1_transverse", "tutorial": Q_CH24 + "#7-obstructions"},
    {"id": "nv_centre", "family": "unitary", "spec": "examples/quantum/nv_centre.json", "parent": "spin1_atom",
     "tutorial": Q_CH24 + "#7-obstructions"},
    {"id": "spin1_atom", "family": "unitary", "spec": "examples/quantum/spin1_atom.json",
     "tutorial": Q_CH24 + "#5-co-discovery-one-rotation-in-five-fields"},
    {"id": "stoner_wohlfarth", "family": "dissipative", "spec": "examples/memory/stoner_wohlfarth.json",
     "tutorial": M10 + "#2-running-the-constructor"},
    {"id": "sw_isotropic", "base": "stoner_wohlfarth", "drift": {"phi": "h*cos(phi + psi)"},
     "potential": "-h*sin(phi + psi)", "name": "damped classical spin without anisotropy",
     "question": "Does a damped magnetic moment without anisotropy keep its direction?",
     "tutorial": M10 + "#2-running-the-constructor"},
    {"id": "sw_oblique", "base": "stoner_wohlfarth", "params": {"psi": pi / 9},
     "name": "Stoner-Wohlfarth particle, field at 20° to the easy axis",
     "tutorial": M10 + "#2-running-the-constructor"},
    {"id": "sw_easy", "base": "stoner_wohlfarth", "params": {"psi": 0.0},
     "name": "Stoner-Wohlfarth particle, field along the easy axis",
     "tutorial": M9 + "#5-where-a-derivation-stops"},
    # -- one rotation on other carriers (Chapter 24)
    {"id": "nmr_spin", "family": "unitary", "spec": "examples/quantum/nmr_spin.json",
     "tutorial": Q_CH24 + "#4-attaching-the-rotation-to-the-exchange-chain-of-chapter-14"},
    {"id": "nmr_chain", "attach": ("nmr_spin", "chain", 4), "name": "exchange chain of four spins, one flipped",
     "tutorial": Q_CH24 + "#4-attaching-the-rotation-to-the-exchange-chain-of-chapter-14"},
    {"id": "nmr_bosons", "attach": ("nmr_spin", "bosons", 4), "name": "four bosons in two wells",
     "tutorial": Q_CH24 + "#4-attaching-the-rotation-to-the-exchange-chain-of-chapter-14"},
    {"id": "nmr_pair", "attach": ("nmr_spin", "fermion-pair", None), "name": "Cooper-pair level",
     "tutorial": Q_CH24 + "#4-attaching-the-rotation-to-the-exchange-chain-of-chapter-14"},
    {"id": "nmr_spin32", "attach": ("nmr_spin", "spin", 3), "name": "one spin 3/2",
     "tutorial": Q_CH24 + "#4-attaching-the-rotation-to-the-exchange-chain-of-chapter-14"},
    {"id": "nmr_collective", "attach": ("nmr_spin", "collective", 3), "name": "three spins-½ driven together",
     "tutorial": Q_CH24 + "#7-obstructions"},
    {"id": "collective_ising", "base": "nmr_collective", "params": {"lam": 0.5},
     "add_terms": [{"coefficient": "lam", "operator": "Z0 Z1 + Z0 Z2 + Z1 Z2"}],
     "name": "three driven spins with a collective Ising coupling", "parent": "nmr_collective",
     "tutorial": Q_CH24 + "#7-obstructions"},
    {"id": "chain_ising", "base": "nmr_chain", "params": {"lam": 0.5},
     "add_terms": [{"coefficient": "lam", "operator": " + ".join(f"Z{j} Z{k}" for j in range(4)
                                                                 for k in range(j + 1, 4))}],
     "name": "exchange chain with a collective Ising coupling", "parent": "nmr_chain",
     "tutorial": Q_CH24 + "#7-obstructions"},
    {"id": "bose_josephson", "family": "unitary", "spec": "examples/quantum/bose_josephson.json",
     "tutorial": Q_CH24 + "#5-co-discovery-one-rotation-in-five-fields"},
    {"id": "interacting_bosons", "family": "unitary", "spec": "examples/quantum/interacting_bosons.json",
     "parent": "bose_josephson", "tutorial": Q_CH24 + "#7-obstructions"},
    {"id": "cooper_pair", "family": "unitary", "spec": "examples/quantum/cooper_pair.json",
     "tutorial": Q_CH24 + "#5-co-discovery-one-rotation-in-five-fields"},
    {"id": "state_transfer_chain", "family": "unitary", "spec": "examples/quantum/state_transfer_chain.json",
     "tutorial": Q_CH24 + "#4-attaching-the-rotation-to-the-exchange-chain-of-chapter-14"},
    {"id": "module14_chain", "family": "unitary", "spec": "examples/quantum/module14_chain.json",
     "name": "exchange chain of Chapter 14 (bonds 3 and 4)", "parent": "module14_equal", "tutorial": Q_CH24 + "#7-obstructions"},
    {"id": "module14_equal", "base": "module14_chain", "coefficients": {1: 3},
     "name": "exchange chain of three spins with equal bonds", "tutorial": Q_CH24 + "#7-obstructions"},
    {"id": "heteronuclear_spins", "family": "unitary", "spec": "examples/quantum/heteronuclear_spins.json",
     "parent": "homonuclear_spins", "tutorial": Q_CH24 + "#7-obstructions"},
    {"id": "homonuclear_spins", "base": "heteronuclear_spins", "params": {"w2": 1.0},
     "name": "two nuclear spins of one species", "tutorial": Q_CH24 + "#7-obstructions"},
    # -- canonical forms: mechanisms written without a field ("universal"); the page starts from the pitchfork
    {"id": "pitchfork", "family": "dissipative", "spec": "examples/memory/pitchfork.json", "universal": True,
     "tutorial": "docs/tutorial/18_memory_writing_and_retention.md#the-normal-form-on-the-web-page"},
    {"id": "pitchfork_below", "base": "pitchfork", "params": {"eps": -0.5}, "universal": True,
     "name": "pitchfork normal form below the transition",
     "question": "What does the normal form keep when ε is below the transition?",
     "tutorial": "docs/tutorial/18_memory_writing_and_retention.md#13-retention"},
    {"id": "pitchfork_bias", "base": "pitchfork", "drift": {"x": "eps*x - x**3 + h"},
     "potential": "-eps*x**2/2 + x**4/4 - h*x", "params": {"h": 0.2}, "universal": True,
     "name": "pitchfork normal form with a constant bias",
     "question": "How does a constant bias change the write at the pitchfork?",
     "tutorial": "docs/tutorial/18_memory_writing_and_retention.md#the-normal-form-on-the-web-page"},
    {"id": "pitchfork_subcritical", "base": "pitchfork", "drift": {"x": "eps*x + 2*x**3 - x**5"},
     "potential": "-eps*x**2/2 - x**4/2 + x**6/6", "universal": True,
     "name": "subcritical pitchfork normal form",
     "question": "Where does the state go when x = 0 loses stability at a subcritical pitchfork?",
     "tutorial": "docs/tutorial/18_memory_writing_and_retention.md#11-write-points-and-normal-forms"},
    {"id": "rotation_canonical", "attach": ("two_spins", "qubit", None), "universal": True,
     "name": "rotation on its smallest carrier", "tutorial": Q_CH24 + "#3-detaching-the-rotation-from-the-spins-of-chapter-11"},
    {"id": "rotation_axis", "base": "rotation_canonical", "universal": True, "parent": "rotation_canonical",
     "observable": [{"coefficient": 1.0, "operator": "X0"}, {"coefficient": 0.5, "operator": "Z0"}],
     "name": "rotation measured along its axis", "tutorial": Q_CH11},
    # -- memory in model materials
    {"id": "laser", "family": "dissipative", "spec": "examples/memory/laser.json",
     "tutorial": M9 + "#3-running-the-constructor"},
    {"id": "toggle", "family": "dissipative", "spec": "examples/memory/toggle.json",
     "tutorial": "docs/tutorial/15_memory_first_card.md"},
    {"id": "toggle_unequal", "family": "dissipative", "spec": "examples/memory/toggle_unequal.json",
     "tutorial": "docs/tutorial/19_memory_transfer_and_design.md#4-removing-the-obstruction"},
    {"id": "repressor_ring4", "family": "dissipative", "spec": "examples/memory/repressor_ring4.json",
     "tutorial": "docs/tutorial/17_memory_predictions.md#2-worked-example-two-rings-of-repressors"},
    {"id": "ring4_activation", "base": "repressor_ring4", "drift": {"u0": "alpha*u3**n/(1 + u3**n) - u0"},
     "name": "ring of 4 genes with one activation",
     "question": "What does the ring of four repressors do when one repression is replaced by an activation?",
     "tutorial": "docs/tutorial/17_memory_predictions.md#3-control-calculation-one-activation-in-the-ring-of-four"},
    {"id": "repressilator", "family": "dissipative", "spec": "examples/memory/repressilator.json",
     "tutorial": "docs/tutorial/20_memory_phase.md#2-worked-example-the-repressilator"},
    {"id": "repressilator_activation", "base": "repressilator", "drift": {"u0": "alpha*u2**n/(1 + u2**n) - u0"},
     "name": "ring of 3 genes with one activation",
     "question": "What does the ring of three repressors do when one repression is replaced by an activation?",
     "tutorial": "docs/tutorial/17_memory_predictions.md#1-concepts"},
    {"id": "schlogl", "family": "dissipative", "spec": "examples/memory/schlogl.json",
     "tutorial": "docs/tutorial/19_memory_transfer_and_design.md#3-worked-example-from-the-toggle-switch-to-a-chemical-reactor"},
    {"id": "tubes", "family": "dissipative", "spec": "examples/memory/tubes.json",
     "tutorial": M9 + "#5-where-a-derivation-stops"},
    {"id": "tubes_unequal", "family": "dissipative", "spec": "examples/memory/tubes_unequal.json",
     "tutorial": "docs/tutorial/19_memory_transfer_and_design.md#4-removing-the-obstruction"},
    {"id": "colloid_patch", "family": "dissipative", "spec": "examples/memory/colloid_patch.json",
     "tutorial": "docs/tutorial/18_memory_writing_and_retention.md#2-worked-example-anisotropic-colloids-at-a-fluid-interface"},
    {"id": "dipole_patch", "family": "dissipative", "spec": "examples/memory/dipole_patch.json",
     "tutorial": "docs/tutorial/19_memory_transfer_and_design.md#2-worked-example-from-capillary-rods-to-point-dipoles"},
    # -- oscillators and phase locking
    {"id": "van_der_pol", "family": "dissipative", "spec": "examples/memory/oscillators/van_der_pol.json",
     "tutorial": M11 + "#3-running-the-constructor"},
    {"id": "vdp_stiffness", "base": "van_der_pol", "drift": {"y": "mu*(1 - x**2)*y - k*x + bias"},
     "params": {"k": 1.0}, "control": {"name": "k", "range": [0.5, 1.5]},
     "name": "van der Pol oscillator with a modulated stiffness",
     "question": "At which ratio does a modulation of the stiffness lock the van der Pol oscillator?",
     "tutorial": M11 + "#3-running-the-constructor"},
    {"id": "parametron", "family": "dissipative", "spec": "examples/memory/oscillators/parametron.json",
     "tutorial": M11 + "#3-running-the-constructor"},
    {"id": "brusselator", "family": "dissipative", "spec": "examples/memory/oscillators/brusselator.json",
     "tutorial": M11 + "#3-running-the-constructor"},
    {"id": "goodwin", "family": "dissipative", "spec": "examples/memory/oscillators/goodwin.json",
     "tutorial": M11 + "#3-running-the-constructor"},
    {"id": "fitzhugh_nagumo", "family": "dissipative", "spec": "examples/memory/oscillators/fitzhugh_nagumo.json",
     "tutorial": M11 + "#3-running-the-constructor"},
    {"id": "predator_prey", "family": "dissipative", "spec": "examples/memory/oscillators/predator_prey.json",
     "tutorial": M11 + "#3-running-the-constructor"},
    {"id": "lotka_volterra", "family": "dissipative", "spec": "examples/memory/oscillators/lotka_volterra.json",
     "tutorial": M11 + "#3-running-the-constructor"},
    {"id": "josephson", "family": "dissipative", "spec": "examples/memory/oscillators/josephson.json",
     "tutorial": M11 + "#3-running-the-constructor"},
    # -- fields: conservation, dimension and the shape of the write
    {"id": "field_nonconserved", "family": "field", "spec": "examples/memory/fields/nonconserved_1d.json",
     "tutorial": "docs/tutorial/22_memory_time.md#3-fields-conservation-dimension-and-the-shape-of-the-write"},
    {"id": "field_charge_1d", "family": "field", "spec": "examples/memory/fields/conserved_1d_charge.json",
     "tutorial": "docs/tutorial/22_memory_time.md#3-fields-conservation-dimension-and-the-shape-of-the-write"},
    {"id": "field_dipole_1d", "family": "field", "spec": "examples/memory/fields/conserved_1d_dipole.json",
     "tutorial": "docs/tutorial/22_memory_time.md#3-fields-conservation-dimension-and-the-shape-of-the-write"},
    {"id": "field_charge_2d", "family": "field", "spec": "examples/memory/fields/conserved_2d_charge.json",
     "tutorial": "docs/tutorial/22_memory_time.md#3-fields-conservation-dimension-and-the-shape-of-the-write"},
    {"id": "field_dipole_2d", "family": "field", "spec": "examples/memory/fields/conserved_2d_dipole.json",
     "tutorial": "docs/tutorial/22_memory_time.md#3-fields-conservation-dimension-and-the-shape-of-the-write"},
    # -- hysterons: the return to a turning point of a slow drive
    {"id": "rfim_ferromagnet", "family": "hysterons", "spec": "examples/memory/hysterons/rfim_ferromagnet.json",
     "tutorial": M12 + "#3-the-prediction-from-structure"},
    {"id": "rfim_antiferromagnet", "family": "hysterons", "spec": "examples/memory/hysterons/rfim_antiferromagnet.json",
     "tutorial": M12 + "#5-frustration-of-the-couplings-is-not-the-criterion"},
    {"id": "rfim_antiferromagnet_staggered", "family": "hysterons",
     "spec": "examples/memory/hysterons/rfim_antiferromagnet_staggered.json",
     "tutorial": M12 + "#5-frustration-of-the-couplings-is-not-the-criterion"},
    {"id": "antiferromagnetic_chain", "family": "hysterons",
     "spec": "examples/memory/hysterons/antiferromagnetic_chain.json", "tutorial": M12 + "#6-sufficient-not-necessary"},
    {"id": "adsorption_pores", "family": "hysterons", "spec": "examples/memory/hysterons/adsorption_pores.json",
     "tutorial": M12 + "#4-running-the-command"},
    {"id": "soft_spots", "family": "hysterons", "spec": "examples/memory/hysterons/soft_spots.json",
     "tutorial": M12 + "#4-running-the-command"},
    # -- regulation: the return of an output to its set point after a step of an input
    {"id": "reg_pi_loop", "family": "regulation", "spec": "examples/regulation/benchmarks/pi_loop.json",
     "tutorial": R1 + "#3-the-prediction-from-structure", "universal": True},
    {"id": "reg_leaky_integrator", "family": "regulation", "spec": "examples/regulation/benchmarks/leaky_integrator.json",
     "tutorial": R1 + "#3-the-prediction-from-structure", "universal": True},
    {"id": "reg_feedforward_subtractive", "family": "regulation",
     "spec": "examples/regulation/benchmarks/feedforward_subtractive.json", "tutorial": R1 + "#9-exercises"},
    {"id": "chemotaxis_tu2008", "family": "regulation", "spec": "examples/regulation/chemotaxis_tu2008.json",
     "tutorial": R1 + "#4-running-the-command"},
    {"id": "chemotaxis_turnover", "family": "regulation",
     "spec": "examples/regulation/controls/chemotaxis_tu2008_control1.json", "tutorial": R1 + "#5-controls"},
    {"id": "antithetic_briat2016", "family": "regulation", "spec": "examples/regulation/antithetic_briat2016.json",
     "tutorial": R1 + "#4-running-the-command"},
    {"id": "antithetic_hill", "family": "regulation",
     "spec": "examples/regulation/controls/antithetic_briat2016_control1.json", "tutorial": R1 + "#5-controls"},
    {"id": "cruise_control_pi", "family": "regulation", "spec": "examples/regulation/cruise_control_pi.json",
     "tutorial": R1 + "#4-running-the-command"},
    {"id": "cruise_control_p", "family": "regulation",
     "spec": "examples/regulation/controls/cruise_control_pi_control1.json", "tutorial": R1 + "#5-controls"},
    {"id": "envz_ompr", "family": "regulation", "spec": "examples/regulation/envz_ompr.json",
     "tutorial": R1 + "#7-certificates"},
    {"id": "envz_ompr_phosphatase", "family": "regulation",
     "spec": "examples/regulation/controls/envz_ompr_control1.json", "tutorial": R1 + "#5-controls"},
    {"id": "qian2018_quasi", "family": "regulation", "spec": "examples/regulation/qian2018_quasi.json",
     "tutorial": R1 + "#6-leaks"},
    {"id": "qian2018_leaky", "family": "regulation",
     "spec": "examples/regulation/controls/qian2018_quasi_control1.json", "tutorial": R1 + "#6-leaks"},
    {"id": "qian2018_ideal", "family": "regulation",
     "spec": "examples/regulation/controls/qian2018_quasi_control2.json", "tutorial": R1 + "#6-leaks"},
    {"id": "ma2009_nfblb", "family": "regulation", "spec": "examples/regulation/ma2009_nfblb.json",
     "tutorial": R1 + "#6-leaks"},
    {"id": "ma2009_nfblb_saturated", "family": "regulation", "spec": "examples/regulation/ma2009_nfblb_saturated.json",
     "tutorial": R1 + "#6-leaks"},
    {"id": "ma2009_ifflp", "family": "regulation", "spec": "examples/regulation/ma2009_ifflp.json",
     "tutorial": R1 + "#4-running-the-command"},
    # -- computation: which functions of the input history the measured observables of a driven body represent
    {"id": "comp_one_mode", "family": "computation", "spec": "examples/computation/benchmarks/one_mode_map.json",
     "tutorial": C1 + "#3-the-prediction-from-structure", "universal": True},
    {"id": "comp_linear_chain", "family": "computation", "spec": "examples/computation/benchmarks/linear_chain.json",
     "tutorial": C1 + "#3-the-prediction-from-structure", "universal": True},
    {"id": "comp_linear_chain_no_decay", "family": "computation",
     "spec": "examples/computation/controls/linear_chain_no_decay.json", "tutorial": C1 + "#5-controls",
     "universal": True},
    {"id": "comp_running_sum", "family": "computation", "spec": "examples/computation/benchmarks/running_sum.json",
     "tutorial": C1 + "#5-controls", "universal": True},
    {"id": "comp_odd_oscillator", "family": "computation", "spec": "examples/computation/benchmarks/odd_oscillator.json",
     "tutorial": C1 + "#3-the-prediction-from-structure", "universal": True},
    {"id": "comp_odd_oscillator_bias", "family": "computation",
     "spec": "examples/computation/controls/odd_oscillator_bias.json", "tutorial": C1 + "#5-controls",
     "universal": True},
    {"id": "comp_square_cascade", "family": "computation", "spec": "examples/computation/benchmarks/square_cascade.json",
     "tutorial": C1 + "#3-the-prediction-from-structure", "universal": True},
    {"id": "chemotaxis_methylation_tu2008", "family": "computation",
     "spec": "examples/computation/chemotaxis_methylation_tu2008.json", "tutorial": C1 + "#4-running-the-command"},
    {"id": "spin_torque_furuta2018", "family": "computation", "spec": "examples/computation/spin_torque_furuta2018.json",
     "tutorial": C1 + "#6-measurement-noise"},
    {"id": "hodgkin_huxley1952", "family": "computation", "spec": "examples/computation/hodgkin_huxley1952.json",
     "tutorial": C1 + "#7-the-published-bodies"},
    {"id": "mapk_huang_ferrell1996", "family": "computation", "spec": "examples/computation/mapk_huang_ferrell1996.json",
     "tutorial": C1 + "#7-the-published-bodies"},
    {"id": "echo_state_dambre2012", "family": "computation", "spec": "examples/computation/echo_state_dambre2012.json",
     "tutorial": C1 + "#7-the-published-bodies"},
    {"id": "mass_spring_hauser2011", "family": "computation", "spec": "examples/computation/mass_spring_hauser2011.json",
     "tutorial": C1 + "#7-the-published-bodies"},
    # -- heredity: inheritance through a threshold of the size
    {"id": "her_normal_form", "family": "heredity", "spec": "examples/heredity/benchmarks/normal_form.json",
     "tutorial": H1 + "#3-the-prediction-from-structure", "universal": True},
    {"id": "her_normal_form_no_dip", "family": "heredity", "spec": "examples/heredity/controls/normal_form_no_dip.json",
     "tutorial": H1 + "#5-controls", "universal": True},
    {"id": "her_normal_form_lost", "family": "heredity", "spec": "examples/heredity/controls/normal_form_lost.json",
     "tutorial": H1 + "#5-controls", "universal": True},
    {"id": "chiral_autocatalysis_saito2007", "family": "heredity",
     "spec": "examples/heredity/chiral_autocatalysis_saito2007.json", "tutorial": H1 + "#4-running-the-command"},
    {"id": "chiral_autocatalysis_no_dip", "family": "heredity",
     "spec": "examples/heredity/controls/chiral_autocatalysis_no_dip.json", "tutorial": H1 + "#5-controls"},
    {"id": "turing_painter1999", "family": "heredity", "spec": "examples/heredity/turing_painter1999.json",
     "tutorial": H1 + "#7-the-published-bodies"},
    {"id": "active_nematic_duclos2018", "family": "heredity", "spec": "examples/heredity/active_nematic_duclos2018.json",
     "tutorial": H1 + "#7-the-published-bodies"},
    {"id": "ferroelectric_film_lgd", "family": "heredity", "spec": "examples/heredity/ferroelectric_film_lgd.json",
     "tutorial": H1 + "#7-the-published-bodies"},
    {"id": "filament_baczynski2007", "family": "heredity", "spec": "examples/heredity/filament_baczynski2007.json",
     "tutorial": H1 + "#6-lineages"},
    {"id": "filament_short_division", "family": "heredity",
     "spec": "examples/heredity/controls/filament_short_division.json", "tutorial": H1 + "#5-controls"},
    {"id": "polarity_brauns2020", "family": "heredity", "spec": "examples/heredity/polarity_brauns2020.json",
     "tutorial": H1 + "#5-controls"},
    # -- decision: populations that synchronize, follow a bias through a threshold, or leave an unstable state
    {"id": "dec_van_der_pol", "family": "decision", "spec": "examples/decision/synchronization/van_der_pol.json",
     "tutorial": D1 + "#2-synchronization", "universal": True},
    {"id": "dec_fitzhugh_nagumo", "family": "decision", "spec": "examples/decision/synchronization/fitzhugh_nagumo.json",
     "tutorial": D1 + "#7-the-published-populations"},
    {"id": "dec_brusselator", "family": "decision", "spec": "examples/decision/synchronization/brusselator.json",
     "tutorial": D1 + "#7-the-published-populations"},
    {"id": "dec_goodwin", "family": "decision", "spec": "examples/decision/synchronization/goodwin.json",
     "tutorial": D1 + "#7-the-published-populations"},
    {"id": "dec_predator_prey", "family": "decision", "spec": "examples/decision/synchronization/predator_prey.json",
     "tutorial": D1 + "#7-the-published-populations"},
    {"id": "dec_josephson", "family": "decision", "spec": "examples/decision/controls/josephson_no_onset.json",
     "tutorial": D1 + "#6-controls"},
    {"id": "dec_ising", "family": "decision", "spec": "examples/decision/write/ising_glauber1963.json",
     "tutorial": D1 + "#3-the-swept-collective-write", "universal": True},
    {"id": "dec_honeybees", "family": "decision", "spec": "examples/decision/write/honeybees_pais2013.json",
     "tutorial": D1 + "#7-the-published-populations"},
    {"id": "dec_chiral_write", "family": "decision", "spec": "examples/decision/write/chiral_autocatalysis_saito2007.json",
     "tutorial": D1 + "#7-the-published-populations"},
    {"id": "dec_wong_wang_write", "family": "decision", "spec": "examples/decision/write/decision_network_wong_wang2006.json",
     "tutorial": D1 + "#3-the-swept-collective-write"},
    {"id": "dec_cim_write", "family": "decision", "spec": "examples/decision/write/coherent_ising_machine_wang2013.json",
     "tutorial": D1 + "#7-the-published-populations"},
    {"id": "dec_diverse", "family": "decision", "spec": "examples/decision/write/diverse_units_tessone2006.json",
     "tutorial": D1 + "#6-controls"},
    {"id": "dec_diverse_random", "family": "decision", "spec": "examples/decision/controls/diverse_units_random_sample.json",
     "tutorial": D1 + "#6-controls"},
    {"id": "dec_cim_step", "family": "decision", "spec": "examples/decision/passage/coherent_ising_machine_step.json",
     "tutorial": D1 + "#4-passage-from-a-seed", "universal": True},
    {"id": "dec_chiral_step", "family": "decision", "spec": "examples/decision/passage/chiral_autocatalysis_step.json",
     "tutorial": D1 + "#7-the-published-populations"},
    {"id": "dec_wong_wang_step", "family": "decision", "spec": "examples/decision/passage/decision_network_step.json",
     "tutorial": D1 + "#4-passage-from-a-seed"},
    {"id": "dec_macrospin", "family": "decision", "spec": "examples/decision/passage/macrospin_stoner_wohlfarth.json",
     "tutorial": D1 + "#4-passage-from-a-seed"},
    # -- the quantum write: a parametric oscillator swept through its threshold with a bias
    {"id": "kpo_27", "family": "open", "spec": "examples/quantum/open/kerr_parametric_oscillator.json",
     "tutorial": QW + "#3-the-worked-example", "universal": True},
    {"id": "kpo_few", "family": "open", "spec": "examples/quantum/open/few_photon_oscillator.json",
     "tutorial": QW + "#5-the-control-calculation"},
    {"id": "kpo_thermal", "family": "open", "spec": "examples/quantum/open/thermal_control.json",
     "tutorial": QW + "#4-the-results"},
    # -- the convention of a stochastic calculation
    {"id": "log_ito", "family": "stochastic", "spec": "examples/construction/log_signal_ito.json",
     "tutorial": "docs/tutorial/10_stochastic_construction.md#the-convention-changes-a-measurable-consequence"},
    {"id": "log_stratonovich", "family": "stochastic", "spec": "examples/construction/log_signal_stratonovich.json",
     "tutorial": "docs/tutorial/10_stochastic_construction.md#the-convention-changes-a-measurable-consequence"},
]

# ------------------------------------------------------------------------------------------------ edges
# kind: param (A), term (Omega), observable (R), attach (Xi: the detached mechanism written on another carrier),
# codiscovery (Xi: the same target derived in another realization), closure (C), protocol (P)
EDGES = [
    # the canonical forms: one component of the pitchfork normal form changed
    {"id": "pf_below", "from": "pitchfork", "to": "pitchfork_below", "slot": "A", "kind": "param",
     "change": "ε: 1 → −0.5",
     "text": "Below ε = 0 the potential −εx<sup>2</sup>/2 + x<sup>4</sup>/4 has a single minimum. Every "
             "preparation relaxes to x = 0, and nothing of the preparation is kept: {states_text} stable state."},
    {"id": "pf_bias", "from": "pitchfork", "to": "pitchfork_bias", "slot": "A", "kind": "param",
     "change": "h: 0 → 0.2", "reduces": {"h": 0.0},
     "text": "A constant bias h unfolds the pitchfork. Along ε the state of the unfavoured sign now appears at a "
             "fold, at ε = {write_point}, and a field past the threshold switches the state: the write becomes "
             "one-sided (derivation {thr_word}). The same unfolding turns the write of the toggle switch with "
             "unequal promoters and of the magnet in an oblique field into a fold."},
    {"id": "pf_subcritical", "from": "pitchfork", "to": "pitchfork_subcritical", "slot": "Omega", "kind": "term",
     "change": "−x<sup>3</sup> → +2x<sup>3</sup> − x<sup>5</sup>",
     "text": "With a destabilizing cubic term and a quintic term that bounds the motion, x = 0 loses stability at "
             "a subcritical pitchfork, at ε = {write_point}. The state does not grow continuously from zero: it "
             "jumps to one of {states_text} distant states, as the magnetization of a particle switched along its "
             "easy axis."},
    {"id": "pf_to_magnet", "from": "pitchfork", "to": "stoner_wohlfarth", "slot": "Xi", "kind": "codiscovery",
     "target": "symmetric-write", "change": "order parameter → magnetization angle",
     "text": "The normal form is written on the angle of a magnetization: a single-domain particle in a field along "
             "its hard axis. Its derivation ({sym_word}) ends in the same canonical form, "
             "εx − x<sup>3</sup> + h, at h = {write_point}: the two stored states are the two directions of the "
             "magnetization along the easy axis."},
    {"id": "pf_to_toggle", "from": "pitchfork", "to": "toggle", "slot": "Xi", "kind": "codiscovery",
     "target": "symmetric-write", "change": "order parameter → two repressor concentrations",
     "text": "The normal form is written on two genes that repress each other. The derivation ({sym_word}) ends in "
             "the same canonical form at α = {write_point}: the two stored states are the two genes, one "
             "expressed and the other repressed."},
    # the canonical rotation: the rotation detached from the two spins, written on a single spin-1/2
    {"id": "rotation_detached", "from": "two_spins", "to": "rotation_canonical", "slot": "Xi", "kind": "attach",
     "change": "two spins-½ → one spin-½",
     "text": "The rotation detached from the two spins is written on the smallest carrier of su(2), one spin-½: "
             "H = {rate} (sin θ X + cos θ Z)/2 with θ = {theta}°. The rate and the angle are those of the two "
             "spins (derivation {word})."},
    {"id": "rotation_on_axis", "from": "rotation_canonical", "to": "rotation_axis", "slot": "R", "kind": "observable",
     "change": "Z → the rotation axis, X + Z/2",
     "text": "Measured along the axis of the rotation, the observable commutes with H. Its closure has dimension "
             "{closure}: the measured value is conserved."},
    {"id": "rotation_nmr_values", "from": "rotation_canonical", "to": "nmr_spin", "slot": "A", "kind": "param",
     "change": "the values of a nuclear spin: Ω<sub>R</sub> = 2, δ = 0.5",
     "text": "With the fields of a nuclear spin in a radio-frequency field, Ω<sub>R</sub> = 2 and δ = 0.5 in the "
             "rotating frame, the canonical rotation is the spin of magnetic resonance: rate {rate}, angle "
             "{theta}° (derivation {word}). The parameters name the apparatus; the operators and the mechanism "
             "are unchanged."},
    {"id": "rotation_to_chain3", "from": "rotation_canonical", "to": "module14_equal", "slot": "Xi",
     "kind": "codiscovery", "target": "rotation", "change": "one spin-½ → exchange chain of three spins",
     "text": "An exchange chain of three spins with equal bonds carries a rotation as well: rate {rate}, angle "
             "{theta}° (derivation {word})."},
    {"id": "rotation_to_species", "from": "rotation_canonical", "to": "homonuclear_spins", "slot": "Xi",
     "kind": "codiscovery", "target": "rotation", "change": "one spin-½ → two nuclear spins driven together",
     "text": "Two nuclear spins of one species driven by the same field rotate as one collective spin: rate "
             "{rate}, angle {theta}° (derivation {word})."},
    # the spins of Chapter 11
    {"id": "spins_field_off", "from": "two_spins", "to": "two_spins_h0", "slot": "A", "kind": "param",
     "change": "h: 0.5 → 0",
     "text": "Without the transverse field the rotation axis is the correlation Z<sub>0</sub>Z<sub>1</sub> and the "
             "angle is {theta}°. The measured X<sub>0</sub> oscillates as cos({rate:.3g}t): the result of Chapter 11, "
             "x(t) = x(0) cos 2gt."},
    {"id": "spins_measure_z", "from": "two_spins_h0", "to": "two_spins_z0", "slot": "R", "kind": "observable",
     "change": "X<sub>0</sub> → Z<sub>0</sub>",
     "text": "Z<sub>0</sub> commutes with g Z<sub>0</sub>Z<sub>1</sub>. Its closure has dimension {closure}: the "
             "measured value is conserved, and nothing else is needed to predict it. For X<sub>0</sub> the closure "
             "had dimension {from.closure}, X<sub>0</sub> and the correlation Y<sub>0</sub>Z<sub>1</sub>: the same "
             "material requires a different description for a different measurement."},
    {"id": "spins_field_both", "from": "two_spins", "to": "two_spins_both", "slot": "Omega", "kind": "term",
     "change": "+ h X<sub>1</sub>",
     "text": "With a field on the second spin as well, H and X<sub>0</sub> generate an algebra of dimension {dim}; "
             "without the term {cause} it is su(2). The closure of X<sub>0</sub> grows from {from.closure} to "
             "{closure} operators, which move with {frequencies_text} frequencies ({frequency_list}): the measured "
             "signal leaves the Rabi law."},
    {"id": "spins_field_moved", "from": "two_spins", "to": "two_spins_x1", "slot": "Omega", "kind": "term",
     "change": "h X<sub>0</sub> → h X<sub>1</sub>",
     "text": "With the field on the second spin instead, H and X<sub>0</sub> generate an algebra of dimension "
             "{dim}, not su(2); without the term {cause} it is su(2). The rotation is reached through the closure "
             "of the observable (derivation {word}): the commutators of H with X<sub>0</sub> close on {closure} "
             "operators, X<sub>0</sub>, Y<sub>0</sub>Z<sub>1</sub> and Y<sub>0</sub>Y<sub>1</sub>, which move with "
             "{frequencies_text} frequency ({frequency_list}). The measured signal is that of the two spins, with "
             "rate {rate} and angle {theta}°. That H and the observable generate su(2) is sufficient for the Rabi "
             "law, and not necessary."},
    # the bridge: from the rotation of a spin to a stored magnetization
    {"id": "spins_to_spin1", "from": "two_spins_h0", "to": "spin1_transverse", "slot": "Xi", "kind": "attach",
     "change": "two spins-½ → one spin 1",
     "text": "The rotation detached from the two spins is written on one spin 1: a field along x, "
             "H = {rate} J<sub>x</sub>, and the observable 2J<sub>z</sub>. The rate ({rate}) and the angle "
             "({theta}°) are those of the two spins; the representation changes from {from.rep} to {rep}."},
    {"id": "spin1_anisotropy", "from": "spin1_transverse", "to": "spin1_easy_axis", "slot": "Omega", "kind": "term",
     "change": "+ D J<sub>z</sub><sup>2</sup>, D = −2",
     "text": "The term D J<sub>z</sub><sup>2</sup> is quadratic in <b>J</b>. With it the algebra is su(3), of "
             "dimension {dim}, and the constructor names {cause} as the obstruction. The expectation of <b>J</b> "
             "leaves the sphere: the state spreads over the three levels of the spin instead of rotating."},
    {"id": "spin1_to_nv", "from": "spin1_easy_axis", "to": "nv_centre", "slot": "A", "kind": "param",
     "change": "D = −2 → 3, weak fields",
     "text": "The nitrogen-vacancy centre has the same operators with D &gt; 0, its zero-field splitting, and "
             "weak fields. Its algebra has dimension {dim}: the same term obstructs the rotation."},
    {"id": "spin1_classical_aniso", "from": "spin1_easy_axis", "to": "stoner_wohlfarth", "slot": "C",
     "kind": "closure", "change": "closed unitary → classical direction, damped, in the plane of the field",
     "text": "The closure is changed. The spin is replaced by its classical direction, the limit of a large spin "
             "in which the anisotropy energy becomes D j<sup>2</sup> cos<sup>2</sup>θ with a correction of order "
             "1/j; a bath adds damping and noise; the motion is restricted to the plane of the easy axis and the "
             "field. The term that prevented the rotation now gives the energy two minima: the Stoner–Wohlfarth "
             "particle keeps one of {states_text} directions of its magnetization."},
    {"id": "spin1_classical", "from": "spin1_transverse", "to": "sw_isotropic", "slot": "C", "kind": "closure",
     "change": "closed unitary → classical direction, damped, in the plane of the field",
     "text": "The same change of closure without the anisotropy leaves {states_text} stable direction, along the field. "
             "The damped moment relaxes to it and keeps nothing of its preparation."},
    {"id": "sw_anisotropy", "from": "sw_isotropic", "to": "stoner_wohlfarth", "slot": "Omega", "kind": "term",
     "change": "+ anisotropy ½ cos<sup>2</sup>φ",
     "text": "Adding the anisotropy to the damped moment gives {states_text} stored directions. Applied to the closed "
             "spin, the same term prevented the rotation. The two routes from the spin in a field, the term first "
             "or the closure first, end at the same realization."},
    {"id": "sw_oblique_field", "from": "stoner_wohlfarth", "to": "sw_oblique", "slot": "A", "kind": "param",
     "change": "ψ: 90° → 20°",
     "text": "A field at 20° to the easy axis removes one of the two states at a fold, at h = {write_point}, a "
             "point of the Stoner–Wohlfarth astroid. The write becomes one-sided (derivation {thr_word}): a field "
             "past the threshold switches the magnetization, with the delay of the Airy law."},
    {"id": "sw_easy_field", "from": "sw_oblique", "to": "sw_easy", "slot": "A", "kind": "param",
     "change": "ψ: 20° → 0°",
     "text": "A field along the easy axis leaves the reversed state an equilibrium up to h = {write_point}. There "
             "it loses stability at a subcritical pitchfork, and the magnetization jumps to the distant state "
             "along the field."},
    # the same write in other materials
    {"id": "sw_to_toggle", "from": "stoner_wohlfarth", "to": "toggle", "slot": "Xi", "kind": "codiscovery", "target": "symmetric-write",
     "change": "magnetization angle → two repressor concentrations",
     "text": "The symmetric write of the magnet (derivation {from.sym_word}) is derived again in the genetic toggle "
             "switch (derivation {sym_word}). At the write point both reduce to the canonical form "
             "εx − x<sup>3</sup> + h: one bias field chooses which of the two states is written."},
    {"id": "toggle_to_laser", "from": "toggle", "to": "laser", "slot": "Xi", "kind": "codiscovery", "target": "symmetric-write",
     "change": "repressor concentrations → field amplitude and inversion",
     "text": "The single-mode laser reaches the same write at the lasing threshold, P = {write_point} (derivation "
             "{sym_word}): the two signs of the field amplitude are the two stored states."},
    {"id": "laser_to_normal_form", "from": "laser", "to": "pitchfork", "slot": "Xi", "kind": "codiscovery", "target": "symmetric-write",
     "change": "field amplitude and inversion → one order parameter",
     "text": "The pitchfork normal form is the end point of every derivation of the symmetric write (derivation "
             "{sym_word})."},
    {"id": "toggle_to_ring4", "from": "toggle", "to": "repressor_ring4", "slot": "Xi", "kind": "codiscovery", "target": "symmetric-write",
     "change": "two genes → four genes in a ring",
     "text": "A ring of four repressors is also a positive loop and has {states_text} stable states, written at "
             "α = {write_point} (derivation {sym_word})."},
    {"id": "toggle_promoters", "from": "toggle", "to": "toggle_unequal", "slot": "A", "kind": "param",
     "change": "γ: 1 → 1.25", "reduces": {"gamma": 1.0},
     "text": "With one promoter stronger the pitchfork unfolds into a fold at α = {write_point}: the write becomes "
             "one-sided (derivation {thr_word}). Returning γ to 1 restores the symmetric write; the constructor "
             "finds this cusp itself (derivation {sym_word})."},
    {"id": "ring4_activation_edit", "from": "repressor_ring4", "to": "ring4_activation", "slot": "Omega",
     "kind": "term", "change": "one repression → activation",
     "text": "With one repression replaced by an activation the loop is negative. The {from.states_text} stored states "
             "disappear and the concentrations oscillate: a feedback loop can hold two stable states only if it "
             "is positive (Thomas's rule)."},
    {"id": "ring3_activation_edit", "from": "repressilator", "to": "repressilator_activation", "slot": "Omega",
     "kind": "term", "change": "one repression → activation",
     "text": "The ring of three repressors is a negative loop and oscillates. Replacing one repression by an "
             "activation makes the loop positive, and {states_text} stable states appear."},
    {"id": "toggle_unequal_to_schlogl", "from": "toggle_unequal", "to": "schlogl", "slot": "Xi",
     "kind": "codiscovery", "target": "threshold-write", "change": "two repressor concentrations → one concentration",
     "text": "The Schlögl reactor writes one-sidedly at two folds of its feed rate b, the first at "
             "b = {write_point} (derivation {thr_word}). Tuning the autocatalytic feed to the cusp and sweeping through it "
             "gives the symmetric write (derivation {sym_word})."},
    {"id": "sw_oblique_to_toggle_unequal", "from": "sw_oblique", "to": "toggle_unequal", "slot": "Xi",
     "kind": "codiscovery", "target": "threshold-write", "change": "magnetization angle → two repressor concentrations",
     "text": "The one-sided write of the magnet in an oblique field and that of the toggle with unequal promoters "
             "are one mechanism: both reduce at the threshold to the fold μ + x<sup>2</sup> (derivations "
             "{from.thr_word} and {thr_word}), and a swept switch is delayed by the same Airy law."},
    {"id": "sw_easy_to_tubes", "from": "sw_easy", "to": "tubes", "slot": "Xi", "kind": "codiscovery",
     "change": "magnetization angle → two tube conductances",
     "text": "Two equal tubes of an adaptive transport network also lose a state at a subcritical pitchfork, at "
             "μ = {write_point}: the written state is a distant branch, as for a magnet switched along its easy "
             "axis."},
    {"id": "tubes_unequal_edit", "from": "tubes", "to": "tubes_unequal", "slot": "A", "kind": "param",
     "change": "L<sub>2</sub>: L<sub>1</sub> → longer",
     "text": "Two equal tubes of an adaptive network write at a subcritical pitchfork: the written state is a "
             "distant branch. Unequal lengths turn it into a fold at μ = {write_point}."},
    {"id": "colloids_to_dipoles", "from": "colloid_patch", "to": "dipole_patch", "slot": "Xi", "kind": "codiscovery",
     "change": "capillary rotors (period π) → point dipoles (period 2π)",
     "text": "The memory of the caged rods is carried over to in-plane dipoles: the control only rescales the "
             "energy in both, and a state is written by a uniform field and lost by activation."},
    # oscillators
    {"id": "ring3_is_clock", "from": "ring4_activation", "to": "repressilator", "slot": "Xi", "kind": "codiscovery",
     "change": "four genes → three genes",
     "text": "Both negative rings oscillate. The ring of three repressors locks to a drive at {lock_ratio}:1, set "
             "by its cyclic symmetry (derivation {lock_word})."},
    {"id": "ring3_to_vdp", "from": "repressilator", "to": "van_der_pol", "slot": "Xi", "kind": "codiscovery", "target": "phase-locking",
     "change": "gene concentrations → voltage and current",
     "text": "The van der Pol oscillator is locked by a modulation of its bias at {lock_ratio}:1 (derivation "
             "{lock_word}); the locking range has half-width {lock_width} in units of the coupling K."},
    {"id": "vdp_parametron", "from": "van_der_pol", "to": "vdp_stiffness", "slot": "P", "kind": "protocol",
     "change": "bias modulated → stiffness modulated",
     "text": "The same oscillator driven through its stiffness instead of its bias locks at {lock_ratio}:1 "
             "(derivation {lock_word}): the reflection (x, y) → (−x, −y) leaves the drift unchanged at every "
             "stiffness and maps the cycle onto itself half a period later, so the drive acts through the second "
             "harmonic. The locked phase takes two values half a period apart: the "
             "parametron keeps a state in its phase."},
    {"id": "parametron_damping", "from": "vdp_stiffness", "to": "parametron", "slot": "A", "kind": "param",
     "change": "μ: 1 → 0.5",
     "text": "The pumped oscillator of the parametron examples, with weaker nonlinear damping, locks at "
             "{lock_ratio}:1 with half-width {lock_width} K (derivation {lock_word})."},
    {"id": "vdp_to_brusselator", "from": "van_der_pol", "to": "brusselator", "slot": "Xi", "kind": "codiscovery", "target": "phase-locking",
     "change": "electronics → chemical kinetics",
     "text": "The Brusselator locks at {lock_ratio}:1 (derivation {lock_word}), with half-width {lock_width} K."},
    {"id": "brusselator_to_goodwin", "from": "brusselator", "to": "goodwin", "slot": "Xi", "kind": "codiscovery", "target": "phase-locking",
     "change": "chemical kinetics → gene expression",
     "text": "The Goodwin clock locks at {lock_ratio}:1 (derivation {lock_word}), with half-width {lock_width} K."},
    {"id": "goodwin_to_fhn", "from": "goodwin", "to": "fitzhugh_nagumo", "slot": "Xi", "kind": "codiscovery", "target": "phase-locking",
     "change": "gene expression → membrane potential",
     "text": "The FitzHugh–Nagumo neuron is first moved into its oscillating range (letter C) and then locks at "
             "{lock_ratio}:1 (derivation {lock_word}), with half-width {lock_width} K."},
    {"id": "fhn_to_predator_prey", "from": "fitzhugh_nagumo", "to": "predator_prey", "slot": "Xi",
     "kind": "codiscovery", "target": "phase-locking", "change": "membrane potential → prey and predator",
     "text": "The Rosenzweig–MacArthur populations lock at {lock_ratio}:1 (derivation {lock_word}), with "
             "half-width {lock_width} K."},
    {"id": "predator_prey_to_lv", "from": "predator_prey", "to": "lotka_volterra", "slot": "Omega", "kind": "term",
     "change": "logistic growth and saturating predation → linear terms",
     "text": "Without logistic growth and saturating predation the model conserves a quantity. Every closed orbit "
             "is neutral, there is no isolated cycle, and the derivation of phase locking stops at R."},
    {"id": "predator_prey_to_josephson", "from": "predator_prey", "to": "josephson", "slot": "Xi",
     "kind": "codiscovery", "target": "phase-locking", "change": "prey and predator → superconducting phase",
     "text": "The overdamped Josephson junction, driven past its critical current, locks to a microwave current at "
             "{lock_ratio}:1 (derivation {lock_word}), with half-width {lock_width} K: the Shapiro step."},
    # fields
    {"id": "field_conservation", "from": "field_nonconserved", "to": "field_charge_1d", "slot": "C",
     "kind": "closure", "change": "no conservation → conserved density",
     "text": "With the density conserved, modes of long wavelength relax slowly, κ ∝ k<sup>2</sup>: the trace of a "
             "write decays as t<sup>{exponent}</sup> instead of exponentially, and the total it added is kept."},
    {"id": "field_dipole_write", "from": "field_charge_1d", "to": "field_dipole_1d", "slot": "P", "kind": "protocol",
     "change": "a write that adds material → a write that moves material",
     "text": "A write that only moves material carries no conserved charge. Its trace decays faster, as "
             "t<sup>{exponent}</sup>."},
    {"id": "field_dimension", "from": "field_charge_1d", "to": "field_charge_2d", "slot": "Xi", "kind": "attach",
     "change": "d = 1 → 2",
     "text": "In two dimensions the write spreads over an area: t<sup>{exponent}</sup>. The exponent is set by the "
             "conservation law, the dimension and the shape of the write, before any simulation."},
    {"id": "field_dipole_2d_write", "from": "field_charge_2d", "to": "field_dipole_2d", "slot": "P",
     "kind": "protocol", "change": "a write that adds material → a write that moves material",
     "text": "Moving material in two dimensions: t<sup>{exponent}</sup>."},
    # hysterons: the return to a turning point
    {"id": "hyst_couplings_reversed", "from": "rfim_ferromagnet", "to": "rfim_antiferromagnet", "slot": "Omega",
     "kind": "term", "change": "J → −J",
     "text": "Reversing every coupling turns the ferromagnet into an antiferromagnet. The square lattice is bipartite, "
             "so no plaquette becomes frustrated, but each of the {bonds} couplings now closes a frustrated loop through "
             "the uniform field, and {failed} of {total} subloops do not return to their turning point."},
    {"id": "hyst_drive_staggered", "from": "rfim_antiferromagnet", "to": "rfim_antiferromagnet_staggered", "slot": "P",
     "kind": "protocol", "change": "uniform field → staggered field",
     "text": "A field of opposite signs on the two sublattices. Relabeling one sublattice turns these couplings into "
             "the ferromagnet in a uniform field: no loop through the drive is frustrated, and {failed} of {total} "
             "subloops do not return."},
    {"id": "hyst_chain", "from": "rfim_antiferromagnet", "to": "antiferromagnetic_chain", "slot": "Xi",
     "kind": "attach", "change": "square lattice → chain",
     "text": "The antiferromagnet on a chain. Every loop through the field is still frustrated, yet {failed} of {total} "
             "subloops do not return: random antiferromagnetic chains return exactly (Deutsch, Dhar and Narayan, "
             "2004). A balanced network is sufficient for the return, and not necessary."},
    # regulation: the leak, the turnover, a second phosphatase, the speed of the controller, dilution, saturation
    {"id": "reg_leak", "from": "reg_pi_loop", "to": "reg_leaky_integrator", "slot": "Omega", "kind": "term",
     "change": "+ −δz in the integral",
     "text": "A leak of the integral. The integral of the error is lost at rate δ, and a fraction {gain_ratio} of the "
             "step remains: the clamp of the integral gives δλ/(δλ + k<sub>I</sub>)."},
    {"id": "reg_methylation_turnover", "from": "chemotaxis_tu2008", "to": "chemotaxis_turnover", "slot": "Omega",
     "kind": "term", "change": "+ −k<sub>d</sub>(m − m<sub>0</sub>)",
     "text": "A turnover of the methylation level. Its rate no longer depends on the receptor activity alone, the "
             "integrator is gone, and a fraction {gain_ratio} of the step remains."},
    {"id": "reg_second_phosphatase", "from": "envz_ompr", "to": "envz_ompr_phosphatase", "slot": "Omega",
     "kind": "term", "change": "+ OmpR-P → OmpR",
     "text": "A second phosphatase of OmpR-P that does not pass through EnvZ-ADP. The structure of absolute "
             "concentration robustness is broken, and a fraction {gain_ratio} of a step of the total OmpR remains."},
    {"id": "reg_slow_controller", "from": "qian2018_quasi", "to": "qian2018_leaky", "slot": "A", "kind": "param",
     "change": "ε = 0.02 → 1",
     "text": "Controller reactions as slow as dilution. The memory of the error leaks, and {gain_ratio} of the step "
             "remains, against {from.gain_ratio} with fast controller reactions."},
    {"id": "reg_no_dilution", "from": "qian2018_quasi", "to": "qian2018_ideal", "slot": "A", "kind": "param",
     "change": "γ = 1 → 0",
     "text": "No dilution of the controller species. The difference z₂ − z₁ integrates the error exactly, and the "
             "output returns to its set point."},
    {"id": "reg_saturation", "from": "ma2009_nfblb", "to": "ma2009_nfblb_saturated", "slot": "A", "kind": "param",
     "change": "K<sub>CB</sub> = K′<sub>FBB</sub> = 0.1 → 0.01",
     "text": "The enzymes on the buffer node closer to saturation. The buffer node integrates the output more nearly, "
             "and the fraction that remains falls from {from.gain_ratio} to {gain_ratio}."},
    # computation: a bias of the input, the decay of a stage
    {"id": "comp_bias", "from": "comp_odd_oscillator", "to": "comp_odd_oscillator_bias", "slot": "P", "kind": "protocol",
     "change": "offset of the input 0 → 0.5",
     "text": "A constant bias of the input. The body is no longer odd about its steady state, and degree 2 takes "
             "{degree_2} of the capacity, against {from.degree_2} without the bias."},
    {"id": "comp_no_decay", "from": "comp_linear_chain", "to": "comp_linear_chain_no_decay", "slot": "A", "kind": "param",
     "change": "k<sub>1</sub> = 0.5 → 0",
     "text": "The first stage loses its decay. It keeps the running sum of the input, and the chain has no fading "
             "memory: the capacity at every fixed delay vanishes as the sum grows."},
    # heredity: the size at which the body divides
    {"id": "her_no_dip", "from": "her_normal_form", "to": "her_normal_form_no_dip", "slot": "P", "kind": "protocol",
     "change": "division at 1.5 → 2.4 L<sub>c</sub>",
     "text": "The body divides above twice its threshold. The daughters are born at {L_birth}, above L<sub>c</sub> = "
             "{L_c}, and the order passes to them without a dip."},
    {"id": "her_lost", "from": "her_normal_form", "to": "her_normal_form_lost", "slot": "P", "kind": "protocol",
     "change": "division at 1.5 → 1.3 L<sub>c</sub>",
     "text": "The body divides closer to its threshold. The dip is long against the regrowth, ln G = {lnG} over a "
             "generation, and the order dies out without noise."},
    {"id": "her_chiral_no_dip", "from": "chiral_autocatalysis_saito2007", "to": "chiral_autocatalysis_no_dip",
     "slot": "P", "kind": "protocol", "change": "transfer at the content 64 → 100",
     "text": "The volume is transferred at a larger content. The daughter starts at {L_birth} molecules per source "
             "volume, above the critical content {L_c}, and keeps the handedness without a dip."},
    {"id": "her_filament_short", "from": "filament_baczynski2007", "to": "filament_short_division", "slot": "P",
     "kind": "protocol", "change": "severing at 1.6 → 1.4 L<sub>c</sub>",
     "text": "The filament is severed closer to its critical length. The regrowth above L<sub>c</sub> no longer "
             "restores what the dip removes (ln G = {lnG}), and the buckle is lost without noise."},
    # decision: the protocol (a sweep or a step) and the sample of the diversity
    {"id": "dec_sample", "from": "dec_diverse", "to": "dec_diverse_random", "slot": "P", "kind": "protocol",
     "change": "quantile sample → random sample",
     "text": "Every population draws its own random sample of the diversity. Its mean carries a frozen bias of spread "
             "{s_frozen}, against the thermal spread {sigma_thermal}: the outcome is set by the sample, and P = {P} at "
             "a bias of {z} thermal spreads."},
    {"id": "dec_chiral_to_step", "from": "dec_chiral_write", "to": "dec_chiral_step", "slot": "P", "kind": "protocol",
     "change": "sweep → step of the content",
     "text": "The substrate is added at once instead of supplied slowly. The racemic state is left from the seed of the "
             "molecular noise, with d = {d} and the window ln 13.09 in the growth exponent (rate {rate_end} at the "
             "end)."},
    {"id": "dec_ww_to_step", "from": "dec_wong_wang_write", "to": "dec_wong_wang_step", "slot": "P",
     "kind": "protocol", "change": "ramp → step of the stimulus to 30 Hz",
     "text": "The stimulus is switched on at once, as in the source. The spontaneous state becomes a saddle, and the "
             "pools separate from the seed of the gate noise at the rate {rate_end}/s; the window is ln 13.09 in the "
             "growth exponent along the relaxing mean state."},
    {"id": "dec_cim_to_step", "from": "dec_cim_write", "to": "dec_cim_step", "slot": "P", "kind": "protocol",
     "change": "pump ramp → pump step",
     "text": "The pump is stepped above the collective threshold. The in-phase mode grows at {rate_end} from the vacuum "
             "seed, and the window of passage times is {window_law} in the growth exponent."},
    # the convention of a stochastic calculation
    {"id": "stochastic_convention", "from": "log_ito", "to": "log_stratonovich", "slot": "C", "kind": "closure",
     "change": "Itô → Stratonovich",
     "text": "The same equation read in the Stratonovich convention gives log X the mean growth {growth} instead of "
             "{from.growth}: the convention is part of the closure, and it changes a measured rate."},
]

# ------------------------------------------------------------------------------------------------ sequences
# steps: {"node": id} starts; {"edge": id} applies an edge; {"edge": id, "reverse": True} returns along it;
# {"prepare": id} changes the preparation of a unitary node (P, a live edit)
SEQUENCES = [
    # memory in model materials first: the writes and the retention of a stored state
    {"id": "memory-writes", "title": "Writing and retention in magnets, genes and lasers",
     "steps": [
         {"node": "stoner_wohlfarth",
          "text": "A single-domain particle with uniaxial anisotropy in a field along its hard axis keeps one of "
                  "{states_text} directions of its magnetization. Lowering the field through h = {write_point} writes "
                  "one of them (derivation {sym_word})."},
         {"edge": "sw_to_toggle"},
         {"edge": "toggle_promoters"},
         {"edge": "toggle_promoters", "reverse": True, "text": "Equal promoters again."},
         {"edge": "toggle_to_laser"},
         {"edge": "laser_to_normal_form"},
         {"edge": "laser_to_normal_form", "reverse": True, "text": "Back to the laser."},
         {"edge": "toggle_to_laser", "reverse": True, "text": "Back to the toggle switch."},
         {"edge": "toggle_to_ring4"},
         {"edge": "ring4_activation_edit"},
     ]},
    {"id": "canonical", "title": "From one normal form to other mechanisms",
     "steps": [
         {"node": "pitchfork",
          "text": "The pitchfork normal form, dx/dt = εx − x<sup>3</sup>, written without any field. For ε &gt; 0 it "
                  "has {states_text} stable states, and sweeping ε upward through 0 with a weak bias writes one of "
                  "them (derivation {sym_word}). Every symmetric write on this page reduces to this form."},
         {"edge": "pf_below"},
         {"edge": "pf_below", "reverse": True, "text": "Above the transition again: {states_text} stable states."},
         {"edge": "pf_bias"},
         {"edge": "pf_bias", "reverse": True, "text": "Without the bias the two states appear together again."},
         {"edge": "pf_subcritical"},
         {"edge": "pf_subcritical", "reverse": True,
          "text": "The stabilizing cubic term again: the supercritical pitchfork."},
         {"edge": "pf_to_magnet"},
         {"edge": "spin1_classical_aniso", "reverse": True,
          "text": "The closure is changed back: the damped classical direction becomes a spin 1 under closed "
                  "evolution, with the same anisotropy term. The term that gave two stored directions now prevents "
                  "a rotation: the algebra is su(3), of dimension {dim}."},
         {"edge": "spin1_anisotropy", "reverse": True,
          "text": "Without the anisotropy term the spin 1 in a transverse field rotates: rate {rate}, angle "
                  "{theta}° (derivation {word})."},
         {"edge": "spins_to_spin1", "reverse": True,
          "text": "The same rotation written on two coupled spins, H = g Z<sub>0</sub>Z<sub>1</sub> measured "
                  "through X<sub>0</sub>: the correlation dynamics of Chapter 11, rate {rate} (derivation {word})."},
         {"edge": "spins_measure_z"},
     ]},
    {"id": "spin-to-magnet", "title": "From a spin rotation to a stored magnetization",
     "steps": [
         {"node": "two_spins",
          "text": "Two spins with H = g Z<sub>0</sub>Z<sub>1</sub> + h X<sub>0</sub>, measured through X<sub>0</sub>. "
                  "The commutators of H and X<sub>0</sub> close on X<sub>0</sub>, Y<sub>0</sub>Z<sub>1</sub> and "
                  "Z<sub>0</sub>Z<sub>1</sub>: an su(2) made of one magnetization and two correlations. The "
                  "expectation of this three-vector rotates at {rate} about an axis at {theta}° to the observable "
                  "(derivation {word})."},
         {"edge": "spins_field_both"},
         {"edge": "spins_field_both", "reverse": True,
          "text": "Without the field on the second spin the algebra is su(2) again, and the closure of "
                  "X<sub>0</sub> has {closure} operators."},
         {"edge": "spins_field_off"},
         {"prepare": "plus_y",
          "text": "Two product states with the same x(0) = 0, the first spin along +y or −y, differ in the "
                  "correlation Y<sub>0</sub>Z<sub>1</sub>, and their signals have opposite signs. X<sub>0</sub> at "
                  "t = 0 does not determine its later values."},
         {"prepare": "top", "text": "Back to the top eigenstate of X<sub>0</sub>."},
         {"edge": "spins_measure_z"},
         {"edge": "spins_measure_z", "reverse": True,
          "text": "Measuring X<sub>0</sub> again: the closure has dimension {closure}."},
         {"edge": "spins_to_spin1"},
         {"edge": "spin1_anisotropy"},
         {"edge": "spin1_classical_aniso"},
         {"edge": "sw_oblique_field"},
         {"edge": "sw_easy_field"},
     ]},
    {"id": "one-rotation", "title": "One rotation on six carriers",
     "steps": [
         {"node": "nmr_spin",
          "text": "A nuclear spin in a rotating field: H = (Ω<sub>R</sub>/2) X + (δ/2) Z. It rotates at {rate} about "
                  "an axis at {theta}° to the measured Z (derivation {word})."},
         {"edge": "nmr_to_chain"},
         {"edge": "nmr_to_bosons"},
         {"edge": "nmr_to_pair"},
         {"edge": "nmr_to_spin32"},
         {"edge": "nmr_to_collective"},
         {"edge": "collective_ising_edit"},
         {"edge": "collective_ising_edit", "reverse": True, "text": "Removing the collective coupling."},
         {"edge": "nmr_to_collective", "reverse": True, "text": "Back to the nuclear spin."},
         {"edge": "nmr_to_chain"},
         {"edge": "chain_ising_edit"},
     ]},
    {"id": "phase-locking", "title": "Phase locking in eight oscillators",
     "steps": [
         {"node": "repressilator",
          "text": "Three genes in a negative ring have no stable state: the concentrations follow a limit cycle, "
                  "whose phase is kept only as long as noise allows. A modulation of the synthesis rate fixes "
                  "that phase at {lock_ratio}:1 (derivation {lock_word})."},
         {"edge": "ring3_to_vdp"},
         {"edge": "vdp_parametron"},
         {"edge": "parametron_damping"},
         {"edge": "parametron_damping", "reverse": True, "text": "Stronger damping again."},
         {"edge": "vdp_parametron", "reverse": True, "text": "The bias is modulated again."},
         {"edge": "vdp_to_brusselator"},
         {"edge": "brusselator_to_goodwin"},
         {"edge": "goodwin_to_fhn"},
         {"edge": "fhn_to_predator_prey"},
         {"edge": "predator_prey_to_lv"},
         {"edge": "predator_prey_to_lv", "reverse": True, "text": "Back to the Rosenzweig–MacArthur model."},
         {"edge": "predator_prey_to_josephson"},
     ]},
    {"id": "law-of-loss", "title": "Conservation and the law of loss",
     "steps": [
         {"node": "field_nonconserved",
          "text": "A local write into a field without a conservation law relaxes at the rate of its slowest mode: "
                  "the trace is lost exponentially."},
         {"edge": "field_conservation"},
         {"edge": "field_dipole_write"},
         {"edge": "field_dipole_write", "reverse": True, "text": "A write that adds material again."},
         {"edge": "field_dimension"},
         {"edge": "field_dipole_2d_write"},
     ]},
    {"id": "set-point", "title": "Return to a set point",
     "steps": [
         {"node": "reg_pi_loop",
          "text": "An output y held at a set point by the integral z of its error. After a step of the input the output "
                  "moves and returns exactly to its set point, {set_point}."},
         {"edge": "reg_leak"},
         {"node": "chemotaxis_tu2008",
          "text": "Bacterial chemotaxis. The methylation level integrates a function of the receptor activity, and "
                  "the activity returns to {set_point} after a step of the attractant."},
         {"edge": "reg_methylation_turnover"},
         {"node": "envz_ompr",
          "text": "EnvZ–OmpR. The integrator has a gain proportional to EnvZ-ADP; OmpR-P returns to {set_point} after "
                  "a step of the total OmpR."},
         {"edge": "reg_second_phosphatase"},
         {"edge": "reg_slow_controller"},
     ]},
    {"id": "computation", "title": "Computation by a body",
     "steps": [
         {"node": "comp_one_mode",
          "text": "One linear mode. A linear combination of its state reproduces the input k intervals back with "
                  "capacity (1 − a²)a<sup>2k</sup>; the capacities sum to {n_lin}."},
         {"node": "comp_linear_chain",
          "text": "Three linear stages: {n_lin} independent functions of the input history, all of degree 1."},
         {"edge": "comp_no_decay"},
         {"node": "comp_odd_oscillator",
          "text": "An oscillator with a cubic restoring force. It is odd about rest, so only odd degrees carry capacity: "
                  "{exact}."},
         {"edge": "comp_bias"},
         {"node": "chemotaxis_methylation_tu2008",
          "text": "Bacterial chemotaxis. The kinase activity follows the attractant at once and the methylation keeps "
                  "about one interval of it: {exact}."},
         {"node": "spin_torque_furuta2018",
          "text": "A magnetic tunnel junction measured 50 times per pulse. Its linear response has {n_lin} directions; "
                  "the nonlinear signals add many weak ones, so the capacity depends on how precisely the "
                  "resistance is measured."},
     ]},
    {"id": "heredity", "title": "Inheritance through growth and division",
     "steps": [
         {"node": "her_normal_form",
          "text": "A pitchfork with the size as its control. A daughter is born at {L_birth}, below L<sub>c</sub> = {L_c}; "
                  "growth carries it back through the threshold, and it keeps its parent's sign with P = {P_body}."},
         {"edge": "her_lost"},
         {"node": "chiral_autocatalysis_saito2007",
          "text": "Chiral autocatalysis in a volume diluted by serial transfer. The handedness passes through the "
                  "critical content with P = {P_body}."},
         {"edge": "her_chiral_no_dip"},
         {"node": "filament_baczynski2007",
          "text": "A filament severed below its Euler length under a fixed load. The side of its buckle passes with "
                  "P = {P_body}; the normal form with a linear ramp would give {P_normal_form}, since the decay in the "
                  "dip goes as L<sup>−4</sup>."},
         {"node": "polarity_brauns2020",
          "text": "A polarity made of a redistributed conserved protein. Division gives the halves different amounts, "
                  "and neither daughter polarizes again from the parent's order."},
     ]},
    {"id": "decision", "title": "A population acts as one",
     "steps": [
         {"node": "dec_van_der_pol",
          "text": "Van der Pol units with a 1% spread of frequencies, coupled through their rate of change. The "
                  "reduction gives a<sub>1</sub> = {a1}, b<sub>1</sub> = {b1}, and the incoherent state loses stability "
                  "at K<sub>c</sub> = {K_c}."},
         {"node": "dec_josephson",
          "text": "Overdamped junctions coupled through the supercurrent: b<sub>1</sub> = {b1}, and no coupling "
                  "synchronizes them."},
         {"node": "dec_ising",
          "text": "Spins swept through the Curie point with a weak field choose its direction with P = {P}; the noise "
                  "of the magnetization falls as 1/N."},
         {"node": "dec_wong_wang_write",
          "text": "Two neural pools with a stimulus ramped through the loss of the spontaneous state. The mean state "
                  "lags the ramp, so the law along the passage, P = {P}, exceeds the closed form {P_linear}."},
         {"edge": "dec_ww_to_step"},
         {"node": "dec_diverse",
          "text": "Bistable units with a Gaussian diversity, a symmetric sample: the swept write holds, P = {P}."},
         {"edge": "dec_sample"},
         {"node": "dec_macrospin",
          "text": "A macrospin whose field is reversed: the transverse magnetization rotates while it grows from the "
                  "thermal seed, d = {d}, and the window is {window_law} in the growth exponent."},
     ]},
    {"id": "turning-points", "title": "Return to a turning point",
     "steps": [
         {"node": "rfim_ferromagnet",
          "text": "Domains of a disordered ferromagnet switch at thresholds of a slow field. After an excursion inside "
                  "a turning point the configuration returns: {failed} of {total} subloops do not."},
         {"edge": "hyst_couplings_reversed"},
         {"edge": "hyst_drive_staggered"},
         {"edge": "hyst_drive_staggered", "reverse": True, "text": "The uniform field again."},
         {"edge": "hyst_chain"},
     ]},
]

# the attachments and obstructions of Chapter 24, referred to by the sequence "one-rotation"
EDGES += [
    {"id": "nmr_to_chain", "from": "nmr_spin", "to": "nmr_chain", "slot": "Xi", "kind": "attach",
     "change": "one spin-½ → exchange chain of four spins",
     "text": "On the chain the rotating vector is one flipped spin moving along the chain. The attachment writes "
             "the exchange couplings {couplings}, in the ratio √3 : 2 : √3, and site fields for the detuning. Rate "
             "{rate} and angle {theta}° are kept (derivation {word})."},
    {"id": "nmr_to_bosons", "from": "nmr_spin", "to": "nmr_bosons", "slot": "Xi", "kind": "attach",
     "change": "one spin-½ → four bosons in two wells",
     "text": "Four bosons in two wells carry the rotation as a tunnelling of {tunnelling} with an energy "
             "difference {energy_difference}; the representation is {rep}. Rate {rate}, angle {theta}°."},
    {"id": "nmr_to_pair", "from": "nmr_spin", "to": "nmr_pair", "slot": "Xi", "kind": "attach",
     "change": "one spin-½ → Cooper-pair level",
     "text": "On a pair level k, −k the pairing field rotates Anderson's pseudospin, the empty and the doubly "
             "occupied states; the two singly occupied states do not move. Rate {rate}, angle {theta}°."},
    {"id": "nmr_to_spin32", "from": "nmr_spin", "to": "nmr_spin32", "slot": "Xi", "kind": "attach",
     "change": "one spin-½ → one spin 3/2",
     "text": "A spin 3/2 in the same fields: four states, representation {rep}, the same rate {rate} and angle "
             "{theta}°."},
    {"id": "nmr_to_collective", "from": "nmr_spin", "to": "nmr_collective", "slot": "Xi", "kind": "attach",
     "change": "one spin-½ → three spins-½ driven together",
     "text": "Three spins driven together carry the rotation in the representation {rep}: a collective spin 3/2 "
             "and two spins ½ that the drive does not distinguish. Rate {rate}, angle {theta}°."},
    {"id": "collective_ising_edit", "from": "nmr_collective", "to": "collective_ising", "slot": "Omega",
     "kind": "term", "change": "+ collective Ising coupling λ Σ Z<sub>j</sub>Z<sub>k</sub>",
     "text": "The collective Ising coupling is quadratic in the collective spin. On three freely driven spins it "
             "enlarges the algebra to dimension {dim}, and the constructor names {cause}."},
    {"id": "chain_ising_edit", "from": "nmr_chain", "to": "chain_ising", "slot": "Omega", "kind": "term",
     "change": "+ collective Ising coupling λ Σ Z<sub>j</sub>Z<sub>k</sub>",
     "text": "The same coupling leaves the rotation of the chain unchanged: the sector of one flipped spin fixes "
             "ΣZ<sub>j</sub>, on which the coupling is a constant. The algebra stays su(2) (derivation {word}). The closure, which "
             "declares the sector, removes the obstruction."},
    {"id": "bosons_interaction", "from": "bose_josephson", "to": "interacting_bosons", "slot": "Omega",
     "kind": "term", "change": "+ U (n<sub>a</sub><sup>2</sup> + n<sub>b</sub><sup>2</sup>)/2",
     "text": "The on-site interaction is U J<sub>z</sub><sup>2</sup> within the sector of fixed atom number: "
             "one-axis twisting. The algebra becomes su(5), of dimension {dim}, and the constructor names {cause}."},
    {"id": "chain_bond_edit", "from": "module14_equal", "to": "module14_chain", "slot": "A", "kind": "param",
     "change": "bonds 3, 3 → 3, 4",
     "text": "A chain of three spins carries the rotation only if its two bonds are equal. With bonds 3 and 4 the "
             "algebra is su(3), of dimension {dim}: a parameter value, not a term, removes the mechanism."},
    {"id": "spins_species", "from": "homonuclear_spins", "to": "heteronuclear_spins", "slot": "A", "kind": "param",
     "change": "equal Larmor frequencies → w<sub>2</sub> = 0.4 w<sub>1</sub>",
     "text": "Two spins of one species driven together rotate as one. With different Larmor frequencies they are "
             "two rotations, not one: the algebra has dimension {dim}."},
    {"id": "nmr_to_spin1_atom", "from": "nmr_spin", "to": "spin1_atom", "slot": "Xi", "kind": "codiscovery", "target": "rotation",
     "change": "nuclear spin-½ → atomic spin 1",
     "text": "A spin-1 hyperfine level in a weak field rotates at {rate} (derivation {word})."},
    {"id": "nmr_to_bose_josephson", "from": "nmr_spin", "to": "bose_josephson", "slot": "Xi",
     "kind": "codiscovery", "target": "rotation", "change": "nuclear spin → atoms in two wells",
     "text": "The Bose–Josephson junction without interaction, with its own tunnelling and detuning, rotates at "
             "{rate} about an axis at {theta}° (derivation {word}): a sector of fixed atom number first selects "
             "the states on which the rotation acts."},
    {"id": "nmr_to_cooper", "from": "nmr_spin", "to": "cooper_pair", "slot": "Xi", "kind": "codiscovery",
     "target": "rotation", "change": "nuclear spin → Cooper-pair level",
     "text": "The Cooper-pair level of the co-discovery table rotates at {rate}, angle {theta}° (derivation "
             "{word})."},
    {"id": "chain3_to_chain4", "from": "module14_equal", "to": "state_transfer_chain", "slot": "Xi",
     "kind": "codiscovery", "target": "rotation", "change": "three spins → four spins",
     "text": "A chain of four spins with bonds in the ratio √3 : 2 : √3 transfers a flipped spin from one end to "
             "the other at t = π/{rate} (derivation {word})."},
]

# ------------------------------------------------------------------------------------------------ labels of the map
EDGES += [
    # the quantum write: the Kerr coefficient (the photon number of the stored states) and the thermal bath
    {"id": "kpo_kerr", "from": "kpo_27", "to": "kpo_few", "slot": "A", "kind": "param",
     "change": "K: 0.02 → 1",
     "text": "The Kerr coefficient rises fifty-fold and the stored states hold half a photon instead of 27. The two "
             "shallow wells switch faster than the sweep passes, and the probability of the favoured state leaves the "
             "linear-stage law for the selection of the biased steady state: 0.684 against the law's 0.750."},
    {"id": "kpo_bath", "from": "kpo_27", "to": "kpo_thermal", "slot": "A", "kind": "param",
     "change": "n̄: 0 → 0.6",
     "text": "A thermal bath at occupation 0.6 raises the noise of the seed by 2n̄ + 1 = 2.2. The probability of the "
             "favoured state falls from 0.750 to 0.676, as the law predicts, still with no free parameter."},
]

CLASS_SHORT = {"rotation": "rotation", "conserved": "conserved", "obstructed": "several frequencies",
               "single-state": "one state", "symmetric-write": "symmetric write", "threshold-write": "one-sided write",
               "subcritical-write": "distant write", "field-write": "field write", "return-point": "return point",
               "no-return": "return not exact", "perfect-adaptation": "perfect adaptation",
               "fine-tuned-adaptation": "fine-tuned", "partial-adaptation": "partial adaptation",
               "no-adaptation": "no adaptation", "linear-memory": "linear memory", "odd-capacity": "odd degrees",
               "nonlinear-capacity": "nonlinear capacity", "integrating": "no fading memory",
               "inherited-through-threshold": "inherited", "kept-above-threshold": "kept above",
               "lost-in-the-dip": "lost in the dip", "threshold-moved": "threshold moved",
               "synchronizes": "synchronizes", "no-onset": "no onset", "follows-the-bias": "follows the bias",
               "set-by-the-sample": "set by the sample", "reflection-seed": "one-component seed",
               "rotation-seed": "two-component seed",
               "linear-stage-write": "linear-stage write", "equilibrium-write": "equilibrium write", "oscillation": "limit cycle",
               "neutral-cycles": "neutral cycles", "exponential-loss": "exponential loss", "power-loss": "power-law loss",
               "convention": "convention"}

# the outcomes in which a mechanism of its group is absent: nothing is stored (memory), or nothing rotates as one
# three-vector (quantum); the map marks them with this label and reason
CLASS_ABSENT = {
    "single-state": ["no memory", "every preparation relaxes to one state: nothing of it is kept"],
    "neutral-cycles": ["no memory", "no isolated phase: nothing restores a written phase"],
    "conserved": ["no rotation", "the observable is conserved: nothing rotates"],
    "obstructed": ["no rotation", "the observable moves with several frequencies, not as one rotation"],
    "no-adaptation": ["no regulation", "the output moves to its new value and stays there: nothing returns it"],
    "integrating": ["no capacity", "a mode without decay keeps a running sum: no input at a fixed delay can be "
                                        "recovered from it"],
    "lost-in-the-dip": ["not inherited", "the order dies out over the generations without noise"],
    "threshold-moved": ["not inherited", "division moves the daughters' thresholds; the parent's order becomes a "
                                         "difference between the daughters"],
    "no-onset": ["no synchrony", "the coupling only shifts the frequencies: no coupling makes the incoherent state "
                                 "unstable"],
    "no-return": ["no return point", "a measured subloop does not return to the state at its turning point; a "
                                     "frustrated loop through the drive allows this without forcing it"],
}

# the three laws by which a stored state is lost, defined in Module 4, Section 1.3; the page links every mention
M4_RETENTION = "docs/tutorial/18_memory_writing_and_retention.md#13-retention"
RETENTION_LAWS = {
    "1": "relaxation in a curved minimum (κ &gt; 0): the information about the write decreases as e<sup>−2κt</sup>",
    "2": "diffusion along a direction with κ = 0, such as the phase of a limit cycle: the information decreases as 1/t",
    "3": "activation over a barrier ΔV between stored states, at the Kramers rate, proportional to "
         "e<sup>−ΔV/k<sub>B</sub>T</sup>",
}

# mechanisms in preparation: studied in the tutorial or in the literature on memory, not yet derivation targets with a
# certified law. docs/ROADMAP.md lists the same items.
PLANNED = [
    {"id": "frustrated-loops", "name": "frustrated loops", "group": "memory: writing a state",
     "law": "excess energy per bond 1 − cos(Φ/N) for N equal rotor bonds with mismatch Φ",
     "exists": "the structural prediction and <code>memory loops</code>: 91 loops of genes, spins and rotors behave as "
               "predicted, and rotor frustration agrees with 1 − cos(Φ/N) to 3 × 10<sup>−12</sup>",
     "missing": "a derivation target with its letters, and a browser engine for networks",
     "tutorial": "docs/tutorial/17_memory_predictions.md#4-tests-on-many-loops-and-networks"},
    {"id": "retention-rewriting", "name": "retention against rewriting", "group": "memory: retention",
     "law": "the ratio of the retention time to the writing time depends only on the work E<sub>w</sub> of the write: "
            "it grows as ln E<sub>w</sub> in a curved minimum (Law 1), linearly along a zero mode (Law 2) and as "
            "e<sup>E<sub>w</sub>/k<sub>B</sub>T</sup> behind a barrier (Law 3)",
     "exists": "the memory card computes the ratio for each writing protocol",
     "missing": "a target in which the write of bounded work is a letter, and its comparison across fields",
     "tutorial": "docs/tutorial/18_memory_writing_and_retention.md#14-relation-between-retention-and-writing-times"},
    {"id": "hopf-onset", "name": "onset of oscillation (Hopf)", "group": "memory: phase",
     "law": "ż = (μ + iω)z − |z|<sup>2</sup>z: an amplitude proportional to √μ; under a slow sweep of μ the onset is "
            "delayed (Neishtadt, 1987)",
     "exists": "FieldBridge locates Hopf bifurcations and reports them where they stop a write, as in the ring of "
               "three repressors",
     "missing": "the reduction to the complex amplitude and a certified law"},
]

SHORT = {
    "kpo_27": "27 photons", "kpo_few": "½ photon", "kpo_thermal": "bath, n̄ = 0.6",
    "reg_pi_loop": "integral control", "reg_leaky_integrator": "leaky integral", "reg_feedforward_subtractive":
    "subtractive feedforward", "chemotaxis_tu2008": "E. coli chemotaxis", "chemotaxis_turnover": "chemotaxis + turnover",
    "antithetic_briat2016": "antithetic controller", "antithetic_hill": "Hill controller",
    "cruise_control_pi": "cruise control, PI", "cruise_control_p": "cruise control, P", "envz_ompr": "EnvZ–OmpR",
    "envz_ompr_phosphatase": "EnvZ–OmpR + phosphatase", "qian2018_quasi": "quasi-integral, ε = 0.02",
    "qian2018_leaky": "leaky, ε = 1", "qian2018_ideal": "without dilution", "ma2009_nfblb": "buffer node, K = 0.1",
    "ma2009_nfblb_saturated": "buffer node, K = 0.01", "ma2009_ifflp": "proportioner node",
    "comp_one_mode": "one linear mode", "comp_linear_chain": "three linear stages",
    "comp_linear_chain_no_decay": "stages, first without decay", "comp_running_sum": "running sum",
    "comp_odd_oscillator": "cubic oscillator", "comp_odd_oscillator_bias": "cubic oscillator + bias",
    "comp_square_cascade": "square cascade", "chemotaxis_methylation_tu2008": "chemotaxis, methylation",
    "spin_torque_furuta2018": "tunnel junction", "hodgkin_huxley1952": "Hodgkin–Huxley membrane",
    "mapk_huang_ferrell1996": "MAPK cascade", "echo_state_dambre2012": "echo state network",
    "mass_spring_hauser2011": "mass–spring network",
    "her_normal_form": "normal form", "her_normal_form_no_dip": "divided at 2.4 L_c", "her_normal_form_lost":
    "divided at 1.3 L_c", "chiral_autocatalysis_saito2007": "chiral autocatalysis", "chiral_autocatalysis_no_dip":
    "transfer at 100", "turing_painter1999": "Turing domain", "active_nematic_duclos2018": "active nematic stripe",
    "ferroelectric_film_lgd": "ferroelectric film", "filament_baczynski2007": "filament under load",
    "filament_short_division": "severed at 1.4 L_c", "polarity_brauns2020": "mass-conserving polarity",
    "dec_van_der_pol": "van der Pol units", "dec_fitzhugh_nagumo": "FitzHugh–Nagumo units",
    "dec_brusselator": "Brusselator units", "dec_goodwin": "Goodwin clocks", "dec_predator_prey": "predator–prey patches",
    "dec_josephson": "Josephson junctions", "dec_ising": "mean-field Ising", "dec_honeybees": "honeybee scouts",
    "dec_chiral_write": "chiral autocatalysis, swept", "dec_wong_wang_write": "two neural pools, ramp",
    "dec_cim_write": "parametric oscillators, ramp", "dec_diverse": "diverse units", "dec_diverse_random":
    "diverse units, random sample", "dec_cim_step": "parametric oscillators, step", "dec_chiral_step":
    "chiral autocatalysis, step", "dec_wong_wang_step": "two neural pools, step", "dec_macrospin": "macrospin",
    "two_spins": "two spins", "two_spins_h0": "two spins, h = 0", "two_spins_z0": "measured Z₀",
    "two_spins_x1": "field on spin 1", "two_spins_both": "field on both spins", "spin1_transverse": "spin 1", "spin1_easy_axis": "spin 1 + DJz²",
    "nv_centre": "NV centre", "spin1_atom": "spin-1 atom", "stoner_wohlfarth": "magnet, 90°",
    "sw_isotropic": "damped moment", "sw_oblique": "magnet, 20°", "sw_easy": "magnet, 0°",
    "nmr_spin": "nuclear spin", "nmr_chain": "chain of 4", "nmr_bosons": "bosons, 2 wells", "nmr_pair": "pair level",
    "nmr_spin32": "spin 3/2", "nmr_collective": "3 spins driven", "collective_ising": "3 spins + Ising",
    "chain_ising": "chain + Ising", "bose_josephson": "Bose–Josephson", "interacting_bosons": "bosons + U",
    "cooper_pair": "Cooper pair", "state_transfer_chain": "transfer chain", "module14_chain": "chain 3, 4",
    "module14_equal": "chain 3, 3", "heteronuclear_spins": "two species", "homonuclear_spins": "one species",
    "pitchfork": "pitchfork", "pitchfork_below": "below ε = 0", "pitchfork_bias": "with bias h",
    "pitchfork_subcritical": "subcritical form", "rotation_canonical": "canonical rotation",
    "rotation_axis": "measured on the axis", "laser": "laser", "toggle": "toggle switch", "toggle_unequal": "unequal toggle",
    "repressor_ring4": "ring of 4", "ring4_activation": "ring 4, activation", "repressilator": "repressilator",
    "repressilator_activation": "ring 3, activation", "schlogl": "Schlögl", "tubes": "two tubes",
    "tubes_unequal": "unequal tubes", "colloid_patch": "capillary rotors", "dipole_patch": "dipoles",
    "van_der_pol": "van der Pol", "vdp_stiffness": "stiffness drive", "parametron": "parametron",
    "brusselator": "Brusselator", "goodwin": "Goodwin", "fitzhugh_nagumo": "FitzHugh–Nagumo",
    "predator_prey": "predator–prey", "lotka_volterra": "Lotka–Volterra", "josephson": "Josephson",
    "field_nonconserved": "field, free", "field_charge_1d": "conserved, d = 1", "field_dipole_1d": "dipole, d = 1",
    "field_charge_2d": "conserved, d = 2", "field_dipole_2d": "dipole, d = 2", "log_ito": "Itô", "log_stratonovich": "Stratonovich",
    "rfim_ferromagnet": "ferromagnet", "rfim_antiferromagnet": "antiferromagnet",
    "rfim_antiferromagnet_staggered": "staggered field", "antiferromagnetic_chain": "chain, J < 0",
    "adsorption_pores": "pores", "soft_spots": "soft spots",
}
