"""Theorems with their hypotheses and the conclusions they do and do not support (check 3 of ``fieldbridge review``).

Each entry states the hypotheses as conditions on named attributes of the case to which a manuscript applies the
theorem, the conclusions the theorem supports, and conclusions that are often drawn from it but that it does not
support. A case satisfies a hypothesis only when the attribute is given and meets the condition; an attribute that is
not given leaves the hypothesis unchecked, which the report states.

Condition forms: ``("equals", v)``, ``("at_most", v)``, ``("at_least", v)``, ``("in_open", a, b)``,
``("in_closed_subinterval_of_open", a, b)`` (strictly inside the open interval (a, b)), ``("is", True/False)``,
``("one_of", [...])``.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

THEOREMS: Dict[str, Dict] = {
    "mermin_wagner_hohenberg": {
        "statement": "For short-range interactions in one or two dimensions at positive temperature, a continuous "
                     "symmetry is not broken: there is no conventional long-range order of its order parameter.",
        "hypotheses": [
            ("dimension", ("at_most", 2), "the dimension is one or two"),
            ("short_range_interactions", ("is", True), "the interactions are short-ranged"),
            ("temperature_positive", ("is", True), "the temperature is positive"),
            ("continuous_symmetry", ("is", True), "the order in question breaks a continuous symmetry"),
        ],
        "supports": ["no_long_range_order_of_continuous_symmetry", "logarithmic_displacement_fluctuations"],
        "does_not_support": {
            "zero_melting_temperature": "Two-dimensional crystals melt at a positive temperature through dislocation "
                                        "and disclination unbinding (Halperin and Nelson 1978; Young 1979), in "
                                        "agreement with the theorem.",
            "theorem_violated_by_finite_melting_temperature": "A positive melting temperature is not a long-range "
                                                              "order of a continuous symmetry; it cannot violate the "
                                                              "theorem.",
            "no_bond_orientational_order": "Bond-orientational order is not excluded for crystals (Mermin 1968).",
        },
        "references": ["N. D. Mermin and H. Wagner, Phys. Rev. Lett. 17, 1133 (1966)",
                       "P. C. Hohenberg, Phys. Rev. 158, 383 (1967)",
                       "N. D. Mermin, Phys. Rev. 176, 250 (1968)"],
    },
    "bender_gao_wormald_local_limit": {
        "statement": "The number of labelled 2-connected planar graphs with n vertices and m edges has a Gaussian "
                     "local limit law in m/n.",
        "hypotheses": [
            ("graph_class", ("one_of", ["labelled 2-connected planar"]), "the graphs are labelled, 2-connected and "
                                                                          "planar"),
            ("edges_per_vertex", ("in_closed_subinterval_of_open", 1.0, 3.0),
             "the ratio m/n lies in a closed subinterval of 1 < m/n < 3"),
        ],
        "supports": ["gaussian_count_inside_interval"],
        "does_not_support": {
            "gaussian_count_at_triangulations": "The endpoint m/n = 3 (triangulations) is outside the theorem; "
                                                "triangulations grow with rate 256/27 (Tutte 1962).",
            "rarity_implies_zero_entropy": "A vanishing fraction of graphs can still be an exponential number of "
                                           "graphs.",
        },
        "references": ["E. A. Bender, Z. Gao and N. C. Wormald, Electron. J. Combin. 9, R43 (2002)",
                       "W. T. Tutte, Canad. J. Math. 14, 21 (1962)"],
    },
    "circle_packing_uniqueness_torus": {
        "statement": "A triangulation of the torus has a circle packing with that contact graph, unique up to "
                     "similarity: the flat lattice and every radius are fixed up to one overall scale.",
        "hypotheses": [
            ("contact_graph_triangulation", ("is", True), "the contact graph is a triangulation"),
            ("surface", ("one_of", ["torus", "sphere"]), "the surface is the torus (or the sphere, up to Moebius "
                                                         "maps)"),
        ],
        "supports": ["packing_exists_for_triangulation", "radii_fixed_by_triangulation_up_to_scale"],
        "does_not_support": {
            "radii_can_be_prescribed": "The radii are outputs of the triangulation; a prescribed generic set of radii "
                                       "has no triangulated packing (Connelly et al. 2020, Thm. 8.1).",
            "entropy_of_fixed_system_from_graph_count": "Counting triangulations counts packings with free radii, "
                                                        "not states of a fixed set of disks.",
        },
        "references": ["W. P. Thurston, The Geometry and Topology of Three-Manifolds, Ch. 13 (1980)",
                       "A. I. Bobenko and B. A. Springborn, Trans. AMS 356, 659 (2004)",
                       "R. Connelly, S. Gortler, E. Solomonides and M. Yampolskaya, Discrete Comput. Geom. 64, 734 "
                       "(2020)"],
    },
    "kramers_escape": {
        "statement": "Escape over a barrier E much larger than the temperature occurs at a rate proportional to "
                     "exp(-E/T).",
        "hypotheses": [
            ("barrier_over_temperature", ("at_least", 5.0), "the barrier is large compared with the temperature"),
            ("temperature_positive", ("is", True), "the temperature is positive"),
        ],
        "supports": ["arrhenius_escape_rate"],
        "does_not_support": {
            "divergence_at_positive_temperature": "Activated escape over a finite barrier is slow but finite at "
                                                  "every positive temperature; a divergence at T > 0 needs a "
                                                  "thermodynamically stable state or a barrier that diverges.",
        },
        "references": ["H. A. Kramers, Physica 7, 284 (1940)"],
    },
    "maxwell_constraint_counting": {
        "statement": "A network with z contacts per particle in d dimensions has at least d - z/2 floppy modes per "
                     "particle; for z > 2d the network is hyperstatic and rigid at zero stress.",
        "hypotheses": [
            ("coordination", ("at_least", 0.0), "the mean contact number is known"),
            ("dimension", ("at_least", 1), "the dimension is known"),
        ],
        "supports": ["finite_shear_modulus_at_zero_pressure_when_hyperstatic"],
        "does_not_support": {
            "zero_configurational_entropy": "Rigidity of one network says nothing about how many networks the "
                                            "system has.",
            "ideal_glass": "Hyperstatic stiffness follows from coordination alone, ordered or not.",
        },
        "references": ["J. C. Maxwell, Phil. Mag. 27, 294 (1864)",
                       "M. Wyart, S. R. Nagel and T. A. Witten, Europhys. Lett. 72, 486 (2005)"],
    },
    "centre_manifold_reduction": {
        "statement": "Near a loss of stability with k critical directions and the other modes decaying faster, the "
                     "dynamics reduces to the k critical amplitudes.",
        "hypotheses": [
            ("critical_directions", ("at_least", 1), "the number of critical directions is known"),
            ("spectral_gap", ("is", True), "the other modes decay faster than the critical ones near onset"),
        ],
        "supports": ["normal_form_near_onset"],
        "does_not_support": {
            "one_coordinate_law_with_degenerate_modes": "With two or more critical directions, a one-coordinate "
                                                        "reduction discards part of the dynamics.",
            "global_barriers_from_local_form": "The reduction is local; barriers between distant states lie outside "
                                               "it.",
        },
        "references": ["J. Guckenheimer and P. Holmes, Nonlinear Oscillations, Dynamical Systems, and Bifurcations "
                       "of Vector Fields (1983)"],
    },
}


def _meets(value, cond) -> bool:
    kind = cond[0]
    if kind == "equals":
        return value == cond[1]
    if kind == "at_most":
        return float(value) <= cond[1]
    if kind == "at_least":
        return float(value) >= cond[1]
    if kind == "in_open":
        return cond[1] < float(value) < cond[2]
    if kind == "in_closed_subinterval_of_open":
        v = value if isinstance(value, (list, tuple)) else [value, value]
        return all(cond[1] < float(x) < cond[2] for x in v)
    if kind == "is":
        return bool(value) is cond[1]
    if kind == "one_of":
        return value in cond[1]
    raise ValueError(f"unknown condition {kind}")


def check_use(name: str, case: Dict, inference: str = "") -> Dict:
    """Whether a case satisfies the hypotheses of a theorem and whether the inference drawn is one it supports."""
    if name not in THEOREMS:
        return {"theorem": name, "status": "unknown", "explanation": "The theorem is not in the registry; its "
                "hypotheses must be checked by hand."}
    th = THEOREMS[name]
    failed: List[str] = []
    unchecked: List[str] = []
    for attr, cond, text in th["hypotheses"]:
        if attr not in case:
            unchecked.append(text)
        elif not _meets(case[attr], cond):
            failed.append(f"{text} (given: {case[attr]})")
    reasons = []
    status = "pass"
    if failed:
        status = "fail"
        reasons.append("The case is outside the theorem: " + "; ".join(failed) + ".")
    if inference and inference in th["does_not_support"]:
        status = "fail"
        reasons.append(f"The inference '{inference}' does not follow: {th['does_not_support'][inference]}")
    elif inference and inference not in th["supports"]:
        status = "warn" if status == "pass" else status
        reasons.append(f"The inference '{inference}' is not among the conclusions the registry lists for this "
                       "theorem; check it by hand.")
    if unchecked and status == "pass":
        status = "warn"
    if unchecked:
        reasons.append("Not given, so not checked: " + "; ".join(unchecked) + ".")
    if status == "pass":
        reasons.append("The case satisfies the hypotheses and the inference is a supported conclusion.")
    return {"theorem": name, "status": status, "explanation": " ".join(reasons), "references": th["references"]}


def table() -> List[Tuple[str, str]]:
    return [(k, v["statement"]) for k, v in THEOREMS.items()]
