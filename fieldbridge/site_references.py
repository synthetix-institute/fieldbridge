"""Where each mechanism is defined, where it is derived, and the publications in which it was first found, derived or
measured.

One entry for every mechanism of the map (the classes of site_registry.CLASSES), every law of retention, every law
whose constant FieldBridge certifies, and every mechanism in preparation (site_registry.PLANNED). The web page shows
the entry of the selected mechanism; docs/mechanisms.md is written from these entries by

    python3 -B -m fieldbridge mechanisms --out docs/mechanisms.md

and a test checks that the file is current, that every link reaches a section and that every DOI is well formed.
Every reference was checked against its DOI record (Crossref) or, for works without a DOI, against the citations of
later papers.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List

from . import site_registry as reg

ROOT = Path(__file__).resolve().parents[1]
DOC = "docs/mechanisms.md"
T = "docs/tutorial/"
GLOSSARY = T + "memory_glossary.md"

# ------------------------------------------------------------------------------------------------ the publications
# cite: the reference as printed; short: author and year, as the page shows it; doi: empty for works without one
BIB = {
    # rotation and the quantum closure
    "rabi1937": {"cite": "I. I. Rabi, Phys. Rev. 51, 652 (1937)", "short": "Rabi, 1937", "doi": "10.1103/PhysRev.51.652"},
    "bloch1946": {"cite": "F. Bloch, Phys. Rev. 70, 460 (1946)", "short": "Bloch, 1946", "doi": "10.1103/PhysRev.70.460"},
    "feynman1957": {"cite": "R. P. Feynman, F. L. Vernon and R. W. Hellwarth, J. Appl. Phys. 28, 49 (1957)",
                    "short": "Feynman, Vernon and Hellwarth, 1957", "doi": "10.1063/1.1722572"},
    "bornjordan1925": {"cite": "M. Born and P. Jordan, Z. Phys. 34, 858 (1925)", "short": "Born and Jordan, 1925",
                       "doi": "10.1007/BF01328531"},
    "bhj1926": {"cite": "M. Born, W. Heisenberg and P. Jordan, Z. Phys. 35, 557 (1926)",
                "short": "Born, Heisenberg and Jordan, 1926", "doi": "10.1007/BF01379806"},
    "mori1965": {"cite": "H. Mori, Prog. Theor. Phys. 34, 399 (1965)", "short": "Mori, 1965", "doi": "10.1143/PTP.34.399"},
    "kitagawa1993": {"cite": "M. Kitagawa and M. Ueda, Phys. Rev. A 47, 5138 (1993)", "short": "Kitagawa and Ueda, 1993",
                     "doi": "10.1103/PhysRevA.47.5138"},
    "viswanath1994": {"cite": "V. S. Viswanath and G. Müller, The Recursion Method (Springer, Berlin, 1994)",
                      "short": "Viswanath and Müller, 1994", "doi": "10.1007/978-3-540-48651-0"},
    # stability and bifurcations
    "lyapunov1892": {"cite": "A. M. Lyapunov, The general problem of the stability of motion (Kharkov, 1892); "
                             "English translation in Int. J. Control 55, 531 (1992)",
                     "short": "Lyapunov, 1892", "doi": "10.1080/00207179208934253"},
    "poincare1885": {"cite": "H. Poincaré, Acta Math. 7, 259 (1885)", "short": "Poincaré, 1885", "doi": "10.1007/BF02402204"},
    "landau1937": {"cite": "L. D. Landau, Zh. Eksp. Teor. Fiz. 7, 19 (1937); English translation in Collected Papers of "
                           "L. D. Landau (Pergamon, Oxford, 1965), p. 193",
                   "short": "Landau, 1937", "doi": "10.1016/B978-0-08-010586-4.50034-1"},
    "kondepudi1983": {"cite": "D. K. Kondepudi and G. W. Nelson, Phys. Rev. Lett. 50, 1023 (1983)",
                      "short": "Kondepudi and Nelson, 1983", "doi": "10.1103/PhysRevLett.50.1023"},
    "kondepudi1985": {"cite": "D. K. Kondepudi and G. W. Nelson, Nature 314, 438 (1985)",
                      "short": "Kondepudi and Nelson, 1985", "doi": "10.1038/314438a0"},
    "kondepudi1986": {"cite": "D. K. Kondepudi, F. Moss and P. V. E. McClintock, Physica D 21, 296 (1986)",
                      "short": "Kondepudi, Moss and McClintock, 1986", "doi": "10.1016/0167-2789(86)90006-0"},
    "griffiths1970": {"cite": "R. B. Griffiths, Phys. Rev. Lett. 24, 715 (1970)", "short": "Griffiths, 1970",
                      "doi": "10.1103/PhysRevLett.24.715"},
    "haberman1979": {"cite": "R. Haberman, SIAM J. Appl. Math. 37, 69 (1979)", "short": "Haberman, 1979",
                     "doi": "10.1137/0137006"},
    "jung1990": {"cite": "P. Jung, G. Gray, R. Roy and P. Mandel, Phys. Rev. Lett. 65, 1873 (1990)",
                 "short": "Jung, Gray, Roy and Mandel, 1990", "doi": "10.1103/PhysRevLett.65.1873"},
    "hopf1942": {"cite": "E. Hopf, Ber. Math.-Phys. Kl. Sächs. Akad. Wiss. Leipzig 94, 1 (1942); English translation "
                         "in J. E. Marsden and M. McCracken, The Hopf Bifurcation and Its Applications (Springer, New "
                         "York, 1976)",
                 "short": "Hopf, 1942", "doi": "10.1007/978-1-4612-6374-6"},
    "neishtadt1987": {"cite": "A. I. Neishtadt, Differ. Uravn. 23, 2060 (1987); English translation in Differential "
                              "Equations 23, 1385 (1987)", "short": "Neishtadt, 1987", "doi": ""},
    "baer1989": {"cite": "S. M. Baer, T. Erneux and J. Rinzel, SIAM J. Appl. Math. 49, 55 (1989)",
                 "short": "Baer, Erneux and Rinzel, 1989", "doi": "10.1137/0149003"},
    # magnets
    "stoner1948": {"cite": "E. C. Stoner and E. P. Wohlfarth, Phil. Trans. R. Soc. A 240, 599 (1948)",
                   "short": "Stoner and Wohlfarth, 1948", "doi": "10.1098/rsta.1948.0007"},
    "neel1949": {"cite": "L. Néel, Ann. Géophys. 5, 99 (1949)", "short": "Néel, 1949", "doi": ""},
    "brown1963": {"cite": "W. F. Brown, Phys. Rev. 130, 1677 (1963)", "short": "Brown, 1963", "doi": "10.1103/PhysRev.130.1677"},
    # oscillators
    "poincare1881": {"cite": "H. Poincaré, J. Math. Pures Appl. (3) 7, 375 (1881) and 8, 251 (1882)",
                     "short": "Poincaré, 1881", "doi": ""},
    "vanderpol1926": {"cite": "B. van der Pol, Phil. Mag. 2, 978 (1926)", "short": "van der Pol, 1926",
                      "doi": "10.1080/14786442608564127"},
    "andronov1929": {"cite": "A. A. Andronov, C. R. Acad. Sci. Paris 189, 559 (1929)", "short": "Andronov, 1929", "doi": ""},
    "adler1946": {"cite": "R. Adler, Proc. IRE 34, 351 (1946)", "short": "Adler, 1946", "doi": "10.1109/JRPROC.1946.229930"},
    "winfree1967": {"cite": "A. T. Winfree, J. Theor. Biol. 16, 15 (1967)", "short": "Winfree, 1967",
                    "doi": "10.1016/0022-5193(67)90051-3"},
    "guckenheimer1975": {"cite": "J. Guckenheimer, J. Math. Biol. 1, 259 (1975)", "short": "Guckenheimer, 1975",
                         "doi": "10.1007/BF01273747"},
    "lax1967": {"cite": "M. Lax, Phys. Rev. 160, 290 (1967)", "short": "Lax, 1967", "doi": "10.1103/PhysRev.160.290"},
    "kuramoto1975": {"cite": "Y. Kuramoto, in International Symposium on Mathematical Problems in Theoretical Physics, "
                             "Lecture Notes in Physics 39 (Springer, Berlin, 1975), p. 420",
                     "short": "Kuramoto, 1975", "doi": "10.1007/BFb0013365"},
    "strogatz2000": {"cite": "S. H. Strogatz, Physica D 143, 1 (2000)", "short": "Strogatz, 2000",
                     "doi": "10.1016/S0167-2789(00)00094-4"},
    "lotka1920": {"cite": "A. J. Lotka, Proc. Natl. Acad. Sci. USA 6, 410 (1920)", "short": "Lotka, 1920",
                  "doi": "10.1073/pnas.6.7.410"},
    "volterra1926": {"cite": "V. Volterra, Nature 118, 558 (1926)", "short": "Volterra, 1926", "doi": "10.1038/118558a0"},
    # noise, diffusion and fields
    "uhlenbeck1930": {"cite": "G. E. Uhlenbeck and L. S. Ornstein, Phys. Rev. 36, 823 (1930)",
                      "short": "Uhlenbeck and Ornstein, 1930", "doi": "10.1103/PhysRev.36.823"},
    "shannon1949": {"cite": "C. E. Shannon, Proc. IRE 37, 10 (1949)", "short": "Shannon, 1949",
                    "doi": "10.1109/JRPROC.1949.232969"},
    "einstein1905": {"cite": "A. Einstein, Ann. Phys. 322, 549 (1905)", "short": "Einstein, 1905",
                     "doi": "10.1002/andp.19053220806"},
    "goldstone1961": {"cite": "J. Goldstone, Nuovo Cimento 19, 154 (1961)", "short": "Goldstone, 1961",
                      "doi": "10.1007/BF02812722"},
    "arrhenius1889": {"cite": "S. Arrhenius, Z. Phys. Chem. 4, 226 (1889)", "short": "Arrhenius, 1889",
                      "doi": "10.1515/zpch-1889-0416"},
    "kramers1940": {"cite": "H. A. Kramers, Physica 7, 284 (1940)", "short": "Kramers, 1940",
                    "doi": "10.1016/S0031-8914(40)90098-2"},
    "hanggi1990": {"cite": "P. Hänggi, P. Talkner and M. Borkovec, Rev. Mod. Phys. 62, 251 (1990)",
                   "short": "Hänggi, Talkner and Borkovec, 1990", "doi": "10.1103/RevModPhys.62.251"},
    "fick1855": {"cite": "A. Fick, Ann. Phys. 170, 59 (1855)", "short": "Fick, 1855", "doi": "10.1002/andp.18551700105"},
    "cahn1961": {"cite": "J. W. Cahn, Acta Metall. 9, 795 (1961)", "short": "Cahn, 1961",
                 "doi": "10.1016/0001-6160(61)90182-1"},
    "allen1979": {"cite": "S. M. Allen and J. W. Cahn, Acta Metall. 27, 1085 (1979)", "short": "Allen and Cahn, 1979",
                  "doi": "10.1016/0001-6160(79)90196-2"},
    "hohenberg1977": {"cite": "P. C. Hohenberg and B. I. Halperin, Rev. Mod. Phys. 49, 435 (1977)",
                      "short": "Hohenberg and Halperin, 1977", "doi": "10.1103/RevModPhys.49.435"},
    "ito1944": {"cite": "K. Itô, Proc. Imp. Acad. Tokyo 20, 519 (1944)", "short": "Itô, 1944",
                "doi": "10.3792/pia/1195572786"},
    "stratonovich1966": {"cite": "R. L. Stratonovich, SIAM J. Control 4, 362 (1966)", "short": "Stratonovich, 1966",
                         "doi": "10.1137/0304028"},
    "wong1965": {"cite": "E. Wong and M. Zakai, Ann. Math. Stat. 36, 1560 (1965)", "short": "Wong and Zakai, 1965",
                 "doi": "10.1214/aoms/1177699916"},
    "vankampen1981": {"cite": "N. G. van Kampen, J. Stat. Phys. 24, 175 (1981)", "short": "van Kampen, 1981",
                      "doi": "10.1007/BF01007642"},
    # loops, hysteresis and rewriting
    "toulouse1977": {"cite": "G. Toulouse, Commun. Phys. 2, 115 (1977)", "short": "Toulouse, 1977", "doi": ""},
    "villain1977": {"cite": "J. Villain, J. Phys. C 10, 1717 (1977)", "short": "Villain, 1977",
                    "doi": "10.1088/0022-3719/10/10/014"},
    "harary1953": {"cite": "F. Harary, Michigan Math. J. 2, 143 (1953)", "short": "Harary, 1953",
                   "doi": "10.1307/mmj/1028989917"},
    "teitel1983": {"cite": "S. Teitel and C. Jayaprakash, Phys. Rev. B 27, 598 (1983)", "short": "Teitel and Jayaprakash, 1983",
                   "doi": "10.1103/PhysRevB.27.598"},
    "thomas1981": {"cite": "R. Thomas, in Numerical Methods in the Study of Critical Phenomena, Springer Series in "
                           "Synergetics 9 (Springer, Berlin, 1981), p. 180",
                   "short": "Thomas, 1981", "doi": "10.1007/978-3-642-81703-8_24"},
    "soule2003": {"cite": "C. Soulé, ComPlexUs 1, 123 (2003)", "short": "Soulé, 2003", "doi": "10.1159/000076100"},
    "bell1978": {"cite": "G. I. Bell, Science 200, 618 (1978)", "short": "Bell, 1978", "doi": "10.1126/science.347575"},
    "kryder2008": {"cite": "M. H. Kryder, E. C. Gage, T. W. McDaniel, W. A. Challener, R. E. Rottmayer, G. Ju, Y.-T. "
                           "Hsia and M. F. Erden, Proc. IEEE 96, 1810 (2008)",
                   "short": "Kryder et al., 2008", "doi": "10.1109/JPROC.2008.2004315"},
    "fusi2005": {"cite": "S. Fusi, P. J. Drew and L. F. Abbott, Neuron 45, 599 (2005)", "short": "Fusi, Drew and Abbott, 2005",
                 "doi": "10.1016/j.neuron.2005.02.001"},
    "benna2016": {"cite": "M. K. Benna and S. Fusi, Nat. Neurosci. 19, 1697 (2016)", "short": "Benna and Fusi, 2016",
                  "doi": "10.1038/nn.4401"},
    "preisach1935": {"cite": "F. Preisach, Z. Phys. 94, 277 (1935)", "short": "Preisach, 1935", "doi": "10.1007/BF01349418"},
    "barker1983": {"cite": "J. A. Barker, D. E. Schreiber, B. G. Huth and D. H. Everett, Proc. R. Soc. Lond. A 386, 251 "
                           "(1983)", "short": "Barker et al., 1983", "doi": "10.1098/rspa.1983.0035"},
    "middleton1992": {"cite": "A. A. Middleton, Phys. Rev. Lett. 68, 670 (1992)", "short": "Middleton, 1992",
                      "doi": "10.1103/PhysRevLett.68.670"},
    "sethna1993": {"cite": "J. P. Sethna, K. Dahmen, S. Kartha, J. A. Krumhansl, B. W. Roberts and J. D. Shore, Phys. "
                           "Rev. Lett. 70, 3347 (1993)", "short": "Sethna et al., 1993", "doi": "10.1103/PhysRevLett.70.3347"},
    "keim2019": {"cite": "N. C. Keim, J. D. Paulsen, Z. Zeravcic, S. Sastry and S. R. Nagel, Rev. Mod. Phys. 91, 035002 "
                         "(2019)", "short": "Keim et al., 2019", "doi": "10.1103/RevModPhys.91.035002"},
    "angeli2003": {"cite": "D. Angeli and E. D. Sontag, IEEE Trans. Automat. Control 48, 1684 (2003)",
                   "short": "Angeli and Sontag, 2003", "doi": "10.1109/TAC.2003.817920"},
    "deutsch2004": {"cite": "J. M. Deutsch, A. Dhar and O. Narayan, Phys. Rev. Lett. 92, 227203 (2004)",
                    "short": "Deutsch, Dhar and Narayan, 2004", "doi": "10.1103/PhysRevLett.92.227203"},
    "vanhecke2021": {"cite": "M. van Hecke, Phys. Rev. E 104, 054608 (2021)", "short": "van Hecke, 2021",
                     "doi": "10.1103/PhysRevE.104.054608"},
    # textbooks
    "strogatz1994": {"cite": "S. H. Strogatz, Nonlinear Dynamics and Chaos (Addison-Wesley, Reading, 1994)",
                     "short": "Strogatz, 1994", "doi": ""},
    "kuznetsov2004": {"cite": "Y. A. Kuznetsov, Elements of Applied Bifurcation Theory, 3rd edn (Springer, New York, 2004)",
                      "short": "Kuznetsov, 2004", "doi": "10.1007/978-1-4757-3978-7"},
    "golubitsky1988": {"cite": "M. Golubitsky, I. Stewart and D. G. Schaeffer, Singularities and Groups in Bifurcation "
                               "Theory, Vol. II (Springer, New York, 1988)",
                       "short": "Golubitsky, Stewart and Schaeffer, 1988", "doi": "10.1007/978-1-4612-4574-2"},
    "kuramoto1984": {"cite": "Y. Kuramoto, Chemical Oscillations, Waves, and Turbulence (Springer, Berlin, 1984)",
                     "short": "Kuramoto, 1984", "doi": "10.1007/978-3-642-69689-3"},
    "pikovsky2001": {"cite": "A. Pikovsky, M. Rosenblum and J. Kurths, Synchronization (Cambridge University Press, 2001)",
                     "short": "Pikovsky, Rosenblum and Kurths, 2001", "doi": "10.1017/CBO9780511755743"},
    "risken1984": {"cite": "H. Risken, The Fokker-Planck Equation (Springer, Berlin, 1984)", "short": "Risken, 1984",
                   "doi": "10.1007/978-3-642-96807-5"},
}

# ------------------------------------------------------------------------------------------------ the entries
# defined and derived: (label, what is there, link); sources: the original publications, (key of BIB, what the
# publication found); textbooks: reviews and textbooks; code: the functions named in the tutorial's tables of reference
S24 = T + "24_spin_language.md"
M4 = T + "18_memory_writing_and_retention.md"
M8 = T + "22_memory_time.md"
M12 = T + "27_memory_return_point.md"

READING: Dict[str, Dict] = {
    # ---------------------------------------------------------------- mechanisms of the map (site_registry.CLASSES)
    "rotation": {
        "defined": ("Chapter 24, §1", "Eqs. (1) and (2): su(2) and the Rabi law", S24 + "#1-the-words-of-the-language"),
        "derived": [("Chapter 24, §3", "the rotation detached from two spins and derived on six carriers",
                     S24 + "#3-detaching-the-rotation-from-the-spins-of-chapter-11"),
                    ("Chapter 24, §6", "the rotation reached through the closure of the observable",
                     S24 + "#6-a-larger-algebra-the-closure-of-the-observable")],
        "code": "quantum.language.derive_bloch_rotation, rabi_law, closure_basis",
        "sources": [("rabi1937", "the probability of a transition of a moment in a rotating field (the Rabi law)"),
                    ("bloch1946", "the precession of a nuclear magnetization in a field"),
                    ("feynman1957", "every two-level system as a vector that rotates in three dimensions")],
    },
    "conserved": {
        "defined": ("Chapter 24, §6", "a closure of the observable that contains one operator",
                    S24 + "#6-a-larger-algebra-the-closure-of-the-observable"),
        "derived": [("Chapter 11", "the closure of an observable under the Hamiltonian", T + "11_quantum_closure.md"
                     + "#how-the-program-finds-the-additional-observable")],
        "code": "quantum.language.closure_basis, observable_frequencies",
        "sources": [("bornjordan1925", "the equation of motion of an observable in matrix mechanics"),
                    ("bhj1926", "quantities that commute with the Hamiltonian are constants of the motion")],
    },
    "obstructed": {
        "defined": ("Chapter 24, §6", "a closure of more than three operators moves with several frequencies",
                    S24 + "#6-a-larger-algebra-the-closure-of-the-observable"),
        "derived": [("Chapter 24, §7", "the obstructions, and the term whose removal restores one rotation",
                     S24 + "#7-obstructions")],
        "code": "quantum.language.observable_frequencies, _single_term_cause",
        "sources": [("mori1965", "the time dependence of an observable from its successive commutators with the "
                                 "generator"),
                    ("kitagawa1993", "one-axis twisting, the term J<sub>z</sub><sup>2</sup> named for the junction "
                                     "with interaction")],
        "textbooks": ["viswanath1994"],
    },
    "single-state": {
        "defined": ("Glossary, Writing", "stable state and local relaxation rate κ", GLOSSARY + "#writing"),
        "derived": [("Module 4, §1.1", "the state below the write point", M4 + "#11-write-points-and-normal-forms"),
                    ("Module 8, §2", "the loss of a write at κ &gt; 0",
                     M8 + "#2-one-expression-for-the-three-regimes-of-kappa")],
        "code": "memory.analysis.locate_writes",
        "sources": [("lyapunov1892", "the stability of an equilibrium decided by the linearized motion")],
        "textbooks": ["strogatz1994"],
    },
    "symmetric-write": {
        "defined": ("Glossary, Writing", "supercritical pitchfork, swept write", GLOSSARY + "#writing"),
        "derived": [("Module 4, §1.1–1.2", "the normal form and the law of the swept write, with its derivation",
                     M4 + "#12-swept-writes"),
                    ("Module 9", "the same canonical form and law constant in models from different fields",
                     T + "23_memory_codiscovery.md#4-four-classes-of-derivation")],
        "code": "memory.analysis.normal_form, construct.swept_write_check, codiscovery.derive_symmetric_write",
        "sources": [("poincare1885", "the bifurcation of a family of equilibria as a parameter changes"),
                    ("landau1937", "the order parameter and the symmetric expansion of the free energy at a "
                                   "continuous transition"),
                    ("kondepudi1983", "symmetry breaking at a pitchfork far from equilibrium and its sensitivity to "
                                      "a small bias"),
                    ("kondepudi1985", "the probability that a weak bias selects the state when the bifurcation is "
                                      "crossed at a finite rate"),
                    ("kondepudi1986", "this selection measured in a noisy electronic circuit")],
        "textbooks": ["golubitsky1988", "strogatz1994"],
    },
    "threshold-write": {
        "defined": ("Glossary, Derivations across fields", "threshold write, write field, delay of a switch",
                    GLOSSARY + "#derivations-across-fields"),
        "derived": [("Module 10, §1", "the target and the derivation in seven fields",
                     T + "25_memory_threshold_write.md#1-the-target"),
                    ("Module 10, §3", "the delay law", T + "25_memory_threshold_write.md#3-the-delay-law")],
        "code": "memory.codiscovery.derive_threshold_write, fold_delay_law, delay_constant",
        "sources": [("poincare1885", "the bifurcation of a family of equilibria as a parameter changes"),
                    ("stoner1948", "switching of a single-domain particle when the field removes the occupied "
                                   "minimum (the coercive field)"),
                    ("haberman1979", "the delay of a slowly swept fold, given by the Airy function"),
                    ("jung1990", "the scaling of hysteresis with the rate of the sweep")],
        "textbooks": ["strogatz1994", "kuznetsov2004"],
    },
    "subcritical-write": {
        "defined": ("Glossary, Writing", "subcritical pitchfork", GLOSSARY + "#writing"),
        "derived": [("Module 4, §1.1", "the normal forms and the kind of each write",
                     M4 + "#11-write-points-and-normal-forms")],
        "code": "memory.analysis.normal_form",
        "sources": [("griffiths1970", "the tricritical point, at which the quartic coefficient changes sign and a "
                                      "continuous transition becomes discontinuous")],
        "textbooks": ["strogatz1994"],
    },
    "field-write": {
        "defined": ("Glossary, Description of a material", "role of the control: scale, bias or shape",
                    GLOSSARY + "#description-of-a-material"),
        "derived": [("Module 4, §2", "anisotropic colloids at a fluid interface, written by a field",
                     M4 + "#2-worked-example-anisotropic-colloids-at-a-fluid-interface"),
                    ("Module 4, §1.4", "the relation between retention and writing times",
                     M4 + "#14-relation-between-retention-and-writing-times")],
        "code": "memory.analysis.hold_time, write_test",
        "sources": [("stoner1948", "hysteresis of single-domain particles written by a uniform field"),
                    ("neel1949", "thermally activated loss of the magnetization of fine grains"),
                    ("brown1963", "the thermal fluctuations of a single-domain particle")],
    },
    "return-point": {
        "defined": ("Module 12, §1", "return-point memory", M12 + "#1-the-question"),
        "derived": [("Module 12, §3", "the drive counted as an element: the relabeling and no passing",
                     M12 + "#3-the-prediction-from-structure"),
                    ("Module 12, §4", "six realizations from three fields", M12 + "#4-running-the-command")],
        "code": "memory.hysterons.predict, relabel, check",
        "sources": [("preisach1935", "hysteresis as a population of elementary loops (hysterons)"),
                    ("barker1983", "return-point memory and minor loops in magnets"),
                    ("middleton1992", "the order of states under a monotonic drive (no passing)"),
                    ("sethna1993", "return-point memory in a disordered model of first-order transitions"),
                    ("harary1953", "the balance of a signed graph, here with the drive as one of its vertices"),
                    ("angeli2003", "the same sign condition for a monotone system with an input")],
        "textbooks": ["keim2019"],
    },
    "no-return": {
        "defined": ("Module 12, §5", "the frustration of the couplings is not the criterion",
                    M12 + "#5-frustration-of-the-couplings-is-not-the-criterion"),
        "derived": [("Module 12, §6", "a frustrated loop through the drive allows a failure and does not force one",
                     M12 + "#6-sufficient-not-necessary")],
        "code": "memory.hysterons.predict, check",
        "sources": [("deutsch2004", "random antiferromagnetic chains return exactly, although no passing fails"),
                    ("vanhecke2021", "most transition graphs of three interacting hysterons violate return-point "
                                     "memory")],
    },
    "oscillation": {
        "defined": ("Glossary, Retention", "phase memory, phase response curve, phase diffusion",
                    GLOSSARY + "#retention"),
        "derived": [("Module 6", "memory in the phase of an oscillator", T + "20_memory_phase.md#1-concepts"),
                    ("Module 11", "phase locking in eight oscillators and the locking ratio",
                     T + "26_memory_phase_locking.md#1-the-target")],
        "code": "memory.phase.response_curve, hold, lock; memory.phase_locking.derive, adler_law",
        "sources": [("poincare1881", "the limit cycle"),
                    ("vanderpol1926", "relaxation oscillations of a triode circuit"),
                    ("andronov1929", "self-sustained oscillations as Poincaré's limit cycles"),
                    ("adler1946", "the locking of an oscillator to an injected signal (the Adler equation)"),
                    ("winfree1967", "the phase response and the synchronization of populations of oscillators"),
                    ("guckenheimer1975", "the isochrons, which give the phase of a state near a limit cycle"),
                    ("lax1967", "the diffusion of the phase of a self-sustained oscillator under noise")],
        "textbooks": ["kuramoto1984", "pikovsky2001"],
    },
    "neutral-cycles": {
        "defined": ("Glossary, Derivations across fields", "neutral cycles", GLOSSARY + "#derivations-across-fields"),
        "derived": [("Module 11, §6", "why a drive cannot fix an isolated phase",
                     T + "26_memory_phase_locking.md#6-where-the-derivation-stops")],
        "code": "memory.phase_locking.floquet",
        "sources": [("lotka1920", "undamped oscillations in a model of interacting species"),
                    ("volterra1926", "the predator-prey cycles and their conserved quantity")],
    },
    "exponential-loss": {
        "defined": ("Glossary, Structure and transfer", "non-conserved field", GLOSSARY + "#structure-and-transfer"),
        "derived": [("Module 8, §3", "conservation, dimension and the shape of the write",
                     M8 + "#3-fields-conservation-dimension-and-the-shape-of-the-write")],
        "code": "memory.fields.predicted_law, hartree_mass",
        "sources": [("allen1979", "the relaxation of a non-conserved order parameter"),
                    ("hohenberg1977", "relaxational dynamics without (Model A) and with (Model B) conservation")],
    },
    "power-loss": {
        "defined": ("Glossary, Structure and transfer", "conserved field", GLOSSARY + "#structure-and-transfer"),
        "derived": [("Module 8, §3", "conservation, dimension and the shape of the write",
                     M8 + "#3-fields-conservation-dimension-and-the-shape-of-the-write")],
        "code": "memory.fields.predicted_law, spectral_snr, exponent_check",
        "sources": [("fick1855", "the diffusion of a conserved density"),
                    ("cahn1961", "the relaxation of a conserved composition"),
                    ("hohenberg1977", "relaxational dynamics without (Model A) and with (Model B) conservation")],
    },
    "convention": {
        "defined": ("Chapter 10", "why a change of stochastic coordinate adds a drift",
                    T + "10_stochastic_construction.md"),
        "derived": [("Chapter 10", "the convention changes a measurable consequence",
                     T + "10_stochastic_construction.md#the-convention-changes-a-measurable-consequence")],
        "code": "verification.verify_construction",
        "sources": [("ito1944", "the Itô integral"),
                    ("stratonovich1966", "the Stratonovich integral"),
                    ("wong1965", "smooth noise converges to the Stratonovich reading"),
                    ("vankampen1981", "the convention as part of the physical model")],
        "textbooks": ["risken1984"],
    },
    # ---------------------------------------------------------------- the laws of retention (site_registry.RETENTION_LAWS)
    "law-1": {
        "defined": ("Module 4, §1.3", "the three laws of retention", reg.M4_RETENTION),
        "derived": [("Module 8, §2", "Laws 1 and 2 from one expression",
                     M8 + "#2-one-expression-for-the-three-regimes-of-kappa")],
        "code": "memory.regimes.gaussian_information",
        "sources": [("uhlenbeck1930", "the relaxation of a Brownian particle bound by a linear force"),
                    ("shannon1949", "the information carried through Gaussian noise, ½ ln(1 + SNR)")],
    },
    "law-2": {
        "defined": ("Module 4, §1.3", "the three laws of retention", reg.M4_RETENTION),
        "derived": [("Module 8, §2", "Laws 1 and 2 from one expression",
                     M8 + "#2-one-expression-for-the-three-regimes-of-kappa"),
                    ("Module 6", "the phase of a limit cycle", T + "20_memory_phase.md#1-concepts")],
        "code": "memory.regimes.gaussian_information, memory.phase.hold",
        "sources": [("einstein1905", "free diffusion, ⟨δx<sup>2</sup>⟩ = 2Dt"),
                    ("goldstone1961", "a zero mode from a broken continuous symmetry"),
                    ("lax1967", "the diffusion of the phase of a self-sustained oscillator under noise")],
    },
    "law-3": {
        "defined": ("Module 4, §1.3", "the three laws of retention", reg.M4_RETENTION),
        "derived": [("Module 4, §1.3", "the law named on the memory card; the retention time is calculated",
                     reg.M4_RETENTION)],
        "code": "memory.analysis.hold_time",
        "sources": [("arrhenius1889", "the rate of a reaction proportional to e<sup>−E/k<sub>B</sub>T</sup>"),
                    ("kramers1940", "the rate of escape over a barrier by Brownian motion"),
                    ("neel1949", "the activated loss of the magnetization of fine grains")],
        "textbooks": ["hanggi1990", "risken1984"],
    },
    # ---------------------------------------------------------------- the certified laws (One mechanism in different fields)
    "law-rotation": {
        "name": "Rabi law of the Bloch rotation",
        "law": "f(t) = cos<sup>2</sup>θ + sin<sup>2</sup>θ cos(|Ω|t)",
        "defined": ("Chapter 24, §1", "Eq. (2)", S24 + "#1-the-words-of-the-language"),
        "derived": [("Chapter 24, §5", "one rotation on carriers from five fields",
                     S24 + "#5-co-discovery-one-rotation-in-five-fields")],
        "code": "quantum.language.rabi_law",
        "sources": [("rabi1937", "the probability of a transition of a moment in a rotating field (the Rabi law)"),
                    ("bloch1946", "the precession of a nuclear magnetization in a field")],
    },
    "law-symmetric-write": {
        "name": "swept-write law, constant π<sup>1/4</sup>",
        "law": "P = Φ(π<sup>1/4</sup> h<sub>s</sub> / (D<sub>s</sub><sup>1/2</sup> r<sup>1/4</sup>))",
        "defined": ("Glossary, Derivations across fields", "constant of the write law",
                    GLOSSARY + "#derivations-across-fields"),
        "derived": [("Module 4, §1.2", "derivation in the linear stage of the sweep", M4 + "#12-swept-writes"),
                    ("Module 9, §6", "the constant in every realization",
                     T + "23_memory_codiscovery.md#6-two-invariants-of-the-end-point")],
        "code": "memory.construct.swept_write_check, codiscovery.law_constant",
        "sources": [("kondepudi1985", "the probability that a weak bias selects the state when the bifurcation is "
                                      "crossed at a finite rate"),
                    ("kondepudi1986", "this selection measured in a noisy electronic circuit")],
    },
    "law-threshold-write": {
        "name": "delay of a switch, constant |a<sub>1</sub>′| = 1.0188",
        "law": "μ<sub>switch</sub> = |a<sub>1</sub>′| r<sup>2/3</sup>",
        "defined": ("Glossary, Derivations across fields", "delay of a switch", GLOSSARY + "#derivations-across-fields"),
        "derived": [("Module 10, §3", "the delay law", T + "25_memory_threshold_write.md#3-the-delay-law")],
        "code": "memory.codiscovery.fold_delay_law, delay_constant",
        "sources": [("haberman1979", "the delay of a slowly swept fold, given by the Airy function"),
                    ("jung1990", "the scaling of hysteresis with the rate of the sweep")],
    },
    "law-phase-locking": {
        "name": "locking range of the Adler equation",
        "law": "ψ̇ = Δω − K sin ψ: locked for |Δω| &lt; K",
        "defined": ("Glossary, Derivations across fields", "phase locking, locking range, locking ratio",
                    GLOSSARY + "#derivations-across-fields"),
        "derived": [("Module 11, §1", "the target and its law", T + "26_memory_phase_locking.md#1-the-target"),
                    ("Module 11, §4", "a symmetry sets the ratio", T + "26_memory_phase_locking.md#4-a-symmetry-sets-the-ratio")],
        "code": "memory.phase_locking.adler_law, locked_rate, slip_frequency",
        "sources": [("adler1946", "the locking of an oscillator to an injected signal (the Adler equation)")],
        "textbooks": ["pikovsky2001"],
    },
    # ---------------------------------------------------------------- mechanisms in preparation (site_registry.PLANNED)
    "frustrated-loops": {
        "defined": ("Glossary, Structure and transfer", "holonomy, frustration", GLOSSARY + "#structure-and-transfer"),
        "derived": [("Module 3, §4", "tests on many loops and networks",
                     T + "17_memory_predictions.md#4-tests-on-many-loops-and-networks")],
        "code": "memory.predict, memory.networks, memory.compose",
        "sources": [("toulouse1977", "the frustration of a plaquette in spin glasses"),
                    ("villain1977", "frustration without disorder"),
                    ("harary1953", "the balance of a signed graph"),
                    ("teitel1983", "rotors with a mismatch around each plaquette"),
                    ("thomas1981", "positive loops for several steady states, negative loops for oscillation"),
                    ("soule2003", "the proof of the condition on positive loops")],
    },
    "retention-rewriting": {
        "defined": ("Glossary, Retention", "relation between retention and writing times", GLOSSARY + "#retention"),
        "derived": [("Module 4, §1.4", "the relation and the three kinds of protocol not subject to it",
                     M4 + "#14-relation-between-retention-and-writing-times")],
        "code": "memory.analysis.hold_time, write_test",
        "sources": [("kramers1940", "the rate of escape over a barrier by Brownian motion"),
                    ("bell1978", "a force lowers the barrier and accelerates its crossing"),
                    ("kryder2008", "heating during the write: a landscape changed between writing and retention"),
                    ("fusi2005", "the conflict between the rate of learning and the time of retention in synapses"),
                    ("benna2016", "memory consolidation in synapses")],
    },
    "hopf-onset": {
        "defined": ("Glossary, Writing", "Hopf bifurcation", GLOSSARY + "#writing"),
        "derived": [("Module 4, §1.1", "the kinds of write points; a Hopf bifurcation stops a write",
                     M4 + "#11-write-points-and-normal-forms")],
        "code": "memory.analysis.locate_writes",
        "sources": [("hopf1942", "a periodic solution that branches from an equilibrium when a complex pair of "
                                 "eigenvalues crosses the imaginary axis"),
                    ("neishtadt1987", "the delay of the loss of stability under a slow sweep"),
                    ("baer1989", "the delay of a slow passage through a Hopf bifurcation and its dependence on "
                                 "the initial state")],
        "textbooks": ["kuznetsov2004"],
    },
    "kuramoto": {
        "defined": ("Roadmap", "mechanisms in preparation", "docs/ROADMAP.md#mechanisms-in-preparation"),
        "derived": [],
        "code": "",
        "sources": [("winfree1967", "the phase response and the synchronization of populations of oscillators"),
                    ("kuramoto1975", "the synchronization threshold of a population of coupled phase oscillators")],
        "textbooks": ["kuramoto1984", "strogatz2000"],
    },
}

CERTIFIED = ["law-rotation", "law-symmetric-write", "law-threshold-write", "law-phase-locking"]
RETENTION = ["law-1", "law-2", "law-3"]


# ------------------------------------------------------------------------------------------------ what the page uses
def slug(heading: str) -> str:
    """The anchor that GitHub gives a heading."""
    h = re.sub(r"<[^>]+>", "", heading).strip().lower()
    h = re.sub(r"[^\w\- ]", "", h)
    return h.replace(" ", "-")


def _year(r: Dict) -> int:
    return int(r["short"].rsplit(", ", 1)[1])


def _ref(key: str, what: str = "") -> Dict:
    b = BIB[key]
    out = {"key": key, "cite": b["cite"], "short": b["short"], "url": f"https://doi.org/{b['doi']}" if b["doi"] else ""}
    if what:
        out["what"] = what
    return out


def heading(key: str) -> str:
    """The heading of an entry in docs/mechanisms.md."""
    if key in reg.CLASSES:
        return reg.CLASSES[key][0].upper() + reg.CLASSES[key][1:]
    if key in RETENTION:
        return f"Law {key[-1]}"
    if key in READING and "name" in READING[key]:
        n = READING[key]["name"]
        return n[0].upper() + n[1:]
    p = next(p for p in reg.PLANNED if p["id"] == key)
    return p["name"][0].upper() + p["name"][1:]


def reading(key: str) -> Dict:
    """The entry as the page uses it: links, publications with their DOI addresses, and the section of
    docs/mechanisms.md. Empty for a key without an entry (a class added to the registry before its entry here; the test
    of the references names it)."""
    if key not in READING:
        return {}
    e = READING[key]
    link = lambda t: {"label": t[0], "what": t[1], "link": t[2]}
    return {"defined": link(e["defined"]), "derived": [link(d) for d in e["derived"]], "code": e.get("code", ""),
            "sources": sorted((_ref(k, w) for k, w in e["sources"]), key=_year),
            "textbooks": [_ref(k) for k in e.get("textbooks", [])],
            "doc": DOC + "#" + slug(heading(key))}


# ------------------------------------------------------------------------------------------------ docs/mechanisms.md
def _md_link(label: str, link: str) -> str:
    return f"[{label}]({link[len('docs/'):] if link.startswith('docs/') else '../' + link})"


def _md_ref(r: Dict) -> str:
    text = r["cite"] + (f", [doi:{r['url'][len('https://doi.org/'):]}]({r['url']})" if r["url"] else "")
    return text + (f": {r['what']}" if r.get("what") else "")


def _md_entry(key: str, lead: str) -> List[str]:
    r = reading(key)
    if not r:
        return []
    lines = [f"### {heading(key)}", ""]
    if lead:
        lines += [lead, ""]
    d = r["defined"]
    lines.append(f"- **Defined:** {_md_link(d['label'], d['link'])}, {d['what']}.")
    if r["derived"]:
        lines.append("- **Derived:** " + "; ".join(f"{_md_link(x['label'], x['link'])}, {x['what']}" for x in r["derived"])
                     + (". Code: " + ", ".join(f"`{f}`" for f in re.split(r"[,;] ", r["code"])) + "." if r["code"] else "."))
    else:
        lines.append("- **Derived:** not yet in FieldBridge.")
    lines.append("- **Original publications:**")
    lines += [f"  - {_md_ref(s)}." for s in r["sources"]]
    if r["textbooks"]:
        lines.append("- **Reviews and textbooks:** " + "; ".join(_md_ref(t) for t in r["textbooks"]) + ".")
    return lines + [""]


def mechanisms_doc() -> str:
    """The text of docs/mechanisms.md."""
    out = ["# Mechanisms: definitions, derivations and original publications", "",
           "For each mechanism of the [web page](https://synthetix-institute.github.io/fieldbridge/), each law of "
           "retention, each law whose constant FieldBridge certifies and each mechanism in preparation: where it is "
           "defined, where the tutorial derives it (with the functions that calculate it), and the publications in "
           "which it was first found, derived or measured. The terms of the memory modules are defined in the "
           "[glossary](tutorial/memory_glossary.md), those of the quantum chapters in "
           "[Chapter 24, §1](tutorial/24_spin_language.md#1-the-words-of-the-language). The realizations of each "
           "mechanism cite their own sources in their specification files and in the "
           "[catalog of materials](materials.md).", "",
           "This page is written from [`fieldbridge/site_references.py`](../fieldbridge/site_references.py) by "
           "`python3 -B -m fieldbridge mechanisms --out docs/mechanisms.md`; a test checks that it is current.", "",
           "| Mechanism | Canonical form | Defined | Derived | Original publications |", "| --- | --- | --- | --- | --- |"]
    for c in (c for c in reg.CLASSES if c in READING):
        r = reading(c)
        der = ", ".join(_md_link(x["label"], x["link"]) for x in r["derived"])
        pubs = "; ".join(s["short"] for s in r["sources"])
        out.append(f"| [{heading(c)}](#{slug(heading(c))}) | {reg.MECHANISMS[c]['canonical']} | "
                   f"{_md_link(r['defined']['label'], r['defined']['link'])} | {der} | {pubs} |")
    out += ["", "## Mechanisms of the map", ""]
    for c in reg.CLASSES:
        m = reg.MECHANISMS.get(c, {})
        out += _md_entry(c, f"Canonical form: {m.get('canonical', '')}. {m.get('text', '')}")
    out += ["## Laws of retention", "",
            f"A stored state is lost by one of three laws, set by the form of the landscape at the state "
            f"({_md_link('Module 4, §1.3', reg.M4_RETENTION)}). Laws 1 and 2 are the loss for the two signs κ > 0 and "
            f"κ = 0 of the curvature along the written direction; where κ < 0 the information stops decreasing once the "
            f"expansion outruns the noise, at ½ ln(1 + W), which is how a state is written. Law 3 needs wells separated by a barrier. A field loses a written pattern mode by mode, "
            f"each mode at its own rate κ(k) "
            f"({_md_link('Module 8, §3', M8 + '#3-fields-conservation-dimension-and-the-shape-of-the-write')}; the "
            f"mechanisms [exponential loss](#exponential-loss) and [power-law loss](#power-law-loss)).", ""]
    for k in RETENTION:
        out += _md_entry(k, f"{reg.RETENTION_LAWS[k[-1]][0].upper()}{reg.RETENTION_LAWS[k[-1]][1:]}.")
    out += ["## Certified laws", "",
            "The laws whose constant FieldBridge calculates in every realization that reaches the target, shown on "
            "the page under *One mechanism in different fields*. They describe how a state is written, or how a "
            "quantum observable rotates, and not how a state is lost.", ""]
    for k in CERTIFIED:
        out += _md_entry(k, f"Law: {READING[k]['law']}.")
    out += ["## Mechanisms in preparation", "",
            "Mechanisms that are studied in the tutorial or in the literature on memory and are not yet derivation "
            "targets with a certified law ([roadmap](ROADMAP.md)).", ""]
    for p in reg.PLANNED:
        if p["id"] not in reg.CLASSES:  # a planned mechanism that has entered the map is listed with the map
            out += _md_entry(p["id"], f"Law: {p['law']}.")
    return "\n".join(out).rstrip() + "\n"


def write_doc(out: str | Path = ROOT / DOC) -> Path:
    path = Path(out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(mechanisms_doc(), encoding="utf-8")
    return path
