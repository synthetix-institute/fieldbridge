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
    # regulation: the return of an output to its set point
    "francis1976": {"cite": "B. A. Francis and W. M. Wonham, Automatica 12, 457 (1976)", "short": "Francis and Wonham, 1976",
                    "doi": "10.1016/0005-1098(76)90006-6"},
    "barkai1997": {"cite": "N. Barkai and S. Leibler, Nature 387, 913 (1997)", "short": "Barkai and Leibler, 1997",
                   "doi": "10.1038/43199"},
    "yi2000": {"cite": "T.-M. Yi, Y. Huang, M. I. Simon and J. Doyle, Proc. Natl. Acad. Sci. USA 97, 4649 (2000)",
               "short": "Yi et al., 2000", "doi": "10.1073/pnas.97.9.4649"},
    "sontag2003": {"cite": "E. D. Sontag, Syst. Control Lett. 50, 119 (2003)", "short": "Sontag, 2003",
                   "doi": "10.1016/S0167-6911(03)00136-1"},
    "tu2008": {"cite": "Y. Tu, T. S. Shimizu and H. C. Berg, Proc. Natl. Acad. Sci. USA 105, 14855 (2008)",
               "short": "Tu, Shimizu and Berg, 2008", "doi": "10.1073/pnas.0807569105"},
    "ma2009": {"cite": "W. Ma, A. Trusina, H. El-Samad, W. A. Lim and C. Tang, Cell 138, 760 (2009)",
               "short": "Ma et al., 2009", "doi": "10.1016/j.cell.2009.06.013"},
    "shoval2010": {"cite": "O. Shoval, L. Goentoro, Y. Hart, A. Mayo, E. Sontag and U. Alon, Proc. Natl. Acad. Sci. USA "
                           "107, 15995 (2010)", "short": "Shoval et al., 2010", "doi": "10.1073/pnas.1002352107"},
    "shinar2010": {"cite": "G. Shinar and M. Feinberg, Science 327, 1389 (2010)", "short": "Shinar and Feinberg, 2010",
                   "doi": "10.1126/science.1183372"},
    "briat2016": {"cite": "C. Briat, A. Gupta and M. Khammash, Cell Syst. 2, 15 (2016)",
                  "short": "Briat, Gupta and Khammash, 2016", "doi": "10.1016/j.cels.2016.01.004"},
    "qian2018": {"cite": "Y. Qian and D. Del Vecchio, J. R. Soc. Interface 15, 20170902 (2018)",
                 "short": "Qian and Del Vecchio, 2018", "doi": "10.1098/rsif.2017.0902"},
    "cappelletti2020": {"cite": "D. Cappelletti, A. Gupta and M. Khammash, J. R. Soc. Interface 17, 20200437 (2020)",
                        "short": "Cappelletti, Gupta and Khammash, 2020", "doi": "10.1098/rsif.2020.0437"},
    "ferrell2016": {"cite": "J. E. Ferrell Jr., Cell Syst. 2, 62 (2016)", "short": "Ferrell, 2016",
                    "doi": "10.1016/j.cels.2016.02.006"},
    "khammash2021": {"cite": "M. H. Khammash, Cell Syst. 12, 509 (2021)", "short": "Khammash, 2021",
                     "doi": "10.1016/j.cels.2021.05.020"},
    "astrom2021": {"cite": "K. J. Åström and R. M. Murray, Feedback Systems, 2nd edn (Princeton University Press, 2021)",
                   "short": "Åström and Murray, 2021", "doi": ""},
    # computation: which functions of an input history a driven body represents
    "jaeger2001": {"cite": "H. Jaeger, The \"echo state\" approach to analysing and training recurrent neural networks, "
                           "GMD Report 148 (German National Research Center for Information Technology, 2001)",
                   "short": "Jaeger, 2001", "doi": ""},
    "jaeger2002": {"cite": "H. Jaeger, Short term memory in echo state networks, GMD Report 152 (German National Research "
                           "Center for Information Technology, 2002)", "short": "Jaeger, 2002", "doi": ""},
    "boyd1985": {"cite": "S. Boyd and L. O. Chua, IEEE Trans. Circuits Syst. 32, 1150 (1985)", "short": "Boyd and Chua, 1985",
                 "doi": "10.1109/TCS.1985.1085649"},
    "dambre2012": {"cite": "J. Dambre, D. Verstraeten, B. Schrauwen and S. Massar, Sci. Rep. 2, 514 (2012)",
                   "short": "Dambre et al., 2012", "doi": "10.1038/srep00514"},
    "grigoryeva2015": {"cite": "L. Grigoryeva, J. Henriques, L. Larger and J.-P. Ortega, Sci. Rep. 5, 12858 (2015)",
                       "short": "Grigoryeva et al., 2015", "doi": "10.1038/srep12858"},
    "gonon2020": {"cite": "L. Gonon, L. Grigoryeva and J.-P. Ortega, Physica D 414, 132721 (2020)",
                  "short": "Gonon, Grigoryeva and Ortega, 2020", "doi": "10.1016/j.physd.2020.132721"},
    "herteux2020": {"cite": "J. Herteux and C. Räth, Chaos 30, 123142 (2020)", "short": "Herteux and Räth, 2020",
                    "doi": "10.1063/5.0028993"},
    "kubota2021": {"cite": "T. Kubota, H. Takahashi and K. Nakajima, Phys. Rev. Research 3, 043135 (2021)",
                   "short": "Kubota, Takahashi and Nakajima, 2021", "doi": "10.1103/PhysRevResearch.3.043135"},
    "hulser2023": {"cite": "T. Hülser, F. Köster, L. Jaurigue and K. Lüdge, Nanophotonics 12, 937 (2023)",
                   "short": "Hülser et al., 2023", "doi": "10.1515/nanoph-2022-0415"},
    "koster2024": {"cite": "F. Köster, S. Yanchuk and K. Lüdge, IEEE Trans. Neural Netw. Learn. Syst. 35, 7712 (2024)",
                   "short": "Köster, Yanchuk and Lüdge, 2024", "doi": "10.1109/TNNLS.2022.3220532"},
    "ballarin2024": {"cite": "G. Ballarin, L. Grigoryeva and J.-P. Ortega, J. Mach. Learn. Res. 25(243), 1 (2024)",
                     "short": "Ballarin, Grigoryeva and Ortega, 2024", "doi": ""},
    "hodgkin1952": {"cite": "A. L. Hodgkin and A. F. Huxley, J. Physiol. 117, 500 (1952)", "short": "Hodgkin and Huxley, 1952",
                    "doi": "10.1113/jphysiol.1952.sp004764"},
    "huang1996": {"cite": "C.-Y. F. Huang and J. E. Ferrell, Proc. Natl. Acad. Sci. USA 93, 10078 (1996)",
                  "short": "Huang and Ferrell, 1996", "doi": "10.1073/pnas.93.19.10078"},
    "saito2007": {"cite": "Y. Saito, T. Sugimori and H. Hyuga, J. Phys. Soc. Jpn. 76, 044802 (2007)",
                  "short": "Saito, Sugimori and Hyuga, 2007", "doi": "10.1143/JPSJ.76.044802"},
    "painter1999": {"cite": "K. J. Painter, P. K. Maini and H. G. Othmer, Proc. Natl. Acad. Sci. USA 96, 5549 (1999)",
                    "short": "Painter, Maini and Othmer, 1999", "doi": "10.1073/pnas.96.10.5549"},
    "voituriez2005": {"cite": "R. Voituriez, J.-F. Joanny and J. Prost, Europhys. Lett. 70, 404 (2005)",
                      "short": "Voituriez, Joanny and Prost, 2005", "doi": "10.1209/epl/i2004-10501-2"},
    "duclos2018": {"cite": "G. Duclos et al., Nat. Phys. 14, 728 (2018)", "short": "Duclos et al., 2018",
                   "doi": "10.1038/s41567-018-0099-7"},
    "kretschmer1979": {"cite": "R. Kretschmer and K. Binder, Phys. Rev. B 20, 1065 (1979)",
                       "short": "Kretschmer and Binder, 1979", "doi": "10.1103/PhysRevB.20.1065"},
    "tilley1984": {"cite": "D. R. Tilley and B. Zeks, Solid State Commun. 49, 823 (1984)", "short": "Tilley and Zeks, 1984",
                   "doi": "10.1016/0038-1098(84)90089-9"},
    "baczynski2007": {"cite": "K. Baczynski, R. Lipowsky and J. Kierfeld, Phys. Rev. E 76, 061914 (2007)",
                      "short": "Baczynski, Lipowsky and Kierfeld, 2007", "doi": "10.1103/PhysRevE.76.061914"},
    "hallatschek2007": {"cite": "O. Hallatschek, E. Frey and K. Kroy, Phys. Rev. E 75, 031905 (2007)",
                        "short": "Hallatschek, Frey and Kroy, 2007", "doi": "10.1103/PhysRevE.75.031905"},
    "brauns2020": {"cite": "F. Brauns, J. Halatek and E. Frey, Phys. Rev. X 10, 041036 (2020)",
                   "short": "Brauns, Halatek and Frey, 2020", "doi": "10.1103/PhysRevX.10.041036"},
    "kondepudi1983": {"cite": "D. K. Kondepudi and G. W. Nelson, Phys. Rev. Lett. 50, 1023 (1983)",
                      "short": "Kondepudi and Nelson, 1983", "doi": "10.1103/PhysRevLett.50.1023"},
    "lythe1996": {"cite": "G. D. Lythe, Phys. Rev. E 53, R4271 (1996)", "short": "Lythe, 1996",
                  "doi": "10.1103/PhysRevE.53.R4271"},
    "vandenbroeck1987": {"cite": "C. van den Broeck and P. Mandel, Phys. Lett. A 122, 36 (1987)",
                         "short": "van den Broeck and Mandel, 1987", "doi": "10.1016/0375-9601(87)90771-7"},
    "hauser2011": {"cite": "H. Hauser, A. J. Ijspeert, R. M. Füchslin, R. Pfeifer and W. Maass, Biol. Cybern. 105, 355 (2011)",
                   "short": "Hauser et al., 2011", "doi": "10.1007/s00422-012-0471-0"},
    "torrejon2017": {"cite": "J. Torrejon et al., Nature 547, 428 (2017)", "short": "Torrejon et al., 2017",
                     "doi": "10.1038/nature23011"},
    "furuta2018": {"cite": "T. Furuta, K. Fujii, K. Nakajima, S. Tsunegi, H. Kubota, Y. Suzuki and S. Miwa, Phys. Rev. "
                           "Applied 10, 034063 (2018)", "short": "Furuta et al., 2018", "doi": "10.1103/PhysRevApplied.10.034063"},
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
R1 = T + "28_regulation_set_point.md"
C1 = T + "29_computation_capacity.md"
H1 = T + "31_heredity_threshold.md"

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
    "perfect-adaptation": {
        "defined": ("Regulation 1, §3", "an integrator of the error and a stable steady state",
                    R1 + "#3-the-prediction-from-structure"),
        "derived": [("Regulation 1, §4", "published models from five fields, each with a control",
                     R1 + "#4-running-the-command"),
                    ("Regulation 1, §7", "certificates checked without repeating the search", R1 + "#7-certificates")],
        "code": "regulation.integrator.find_integrator, card.card, certificate.check",
        "sources": [("francis1976", "the internal model principle: a robust regulator contains a model of the signals "
                                    "it rejects"),
                    ("barkai1997", "adaptation in bacterial chemotaxis follows from the structure of the network, not "
                                   "from tuned rates"),
                    ("yi2000", "the chemotaxis network contains integral feedback"),
                    ("sontag2003", "adaptation to a class of inputs implies an internal model of that class"),
                    ("briat2016", "the antithetic pair of controller species integrates the error, also with noise"),
                    ("cappelletti2020", "an integrator, logarithmic or with a gain that depends on the state, behind "
                                        "absolute concentration robustness"),
                    ("tu2008", "the methylation rate of chemoreceptors as a function of their activity alone"),
                    ("shinar2010", "the structural source of absolute concentration robustness, EnvZ–OmpR")],
        "textbooks": ["astrom2021", "khammash2021"],
    },
    "fine-tuned-adaptation": {
        "defined": ("Regulation 1, §3", "a vanishing gain at the stated parameters only",
                    R1 + "#3-the-prediction-from-structure"),
        "derived": [("Regulation 1, §9", "the subtractive feedforward loop: an integrator only on the surface where "
                                         "two paths cancel", R1 + "#9-exercises")],
        "code": "regulation.card.card, sample_parameters",
        "sources": [("ma2009", "three-node topologies that adapt in regions of parameter space"),
                    ("shoval2010", "incoherent feedforward loops that adapt exactly and detect fold changes")],
        "textbooks": ["ferrell2016"],
    },
    "partial-adaptation": {
        "defined": ("Regulation 1, §3", "a leaky integrator leaves the fraction 1/(1 + G) of the step",
                    R1 + "#3-the-prediction-from-structure"),
        "derived": [("Regulation 1, §6", "leaks: dilution in growing cells, Michaelis constants of a buffer node",
                     R1 + "#6-leaks")],
        "code": "regulation.gains.clamped_gain, step.step_response",
        "sources": [("qian2018", "dilution makes integral control in growing cells leaky; fast controller reactions "
                                 "make the error small"),
                    ("ma2009", "near-perfect adaptation when the enzymes on a buffer node are saturated")],
        "textbooks": ["astrom2021"],
    },
    "no-adaptation": {
        "defined": ("Regulation 1, §5", "controls that remove the integrator", R1 + "#5-controls"),
        "derived": [("Regulation 1, §5", "proportional feedback reduces the change and does not return it",
                     R1 + "#5-controls")],
        "code": "regulation.card.card, gains.attenuation",
        "sources": [("briat2016", "a static controller of Hill type does not adapt to a change of the process")],
        "textbooks": ["astrom2021"],
    },
    "linear-memory": {
        "defined": ("Computation 1, §3", "a linear body with linear observables: capacity at degree 1 only",
                    C1 + "#3-the-prediction-from-structure"),
        "derived": [("Computation 1, §3", "the profile of one mode and the rank of the linear response, from the "
                                          "equations", C1 + "#3-the-prediction-from-structure"),
                    ("Computation 1, §6", "the rank counts directions far below any measurement", C1 + "#6-measurement-noise")],
        "code": "computation.predict.linear, ipc.capacities",
        "sources": [("jaeger2002", "the short-term memory of an echo state network; at most the number of its variables"),
                    ("dambre2012", "the capacity of a dynamical system for functions of its input history"),
                    ("gonon2020", "the memory capacity of a linear network equals the rank of its controllability matrix"),
                    ("ballarin2024", "the numerical rank of that matrix and the bias of the estimator")],
        "textbooks": [],
    },
    "odd-capacity": {
        "defined": ("Computation 1, §3", "an odd body, odd observables and a symmetric input: odd degrees only",
                    C1 + "#3-the-prediction-from-structure"),
        "derived": [("Computation 1, §5", "a bias of the input brings the even degrees", C1 + "#5-controls")],
        "code": "computation.predict.odd, ipc.capacities",
        "sources": [("dambre2012", "the capacities of an echo state network vanish at even degrees"),
                    ("herteux2020", "the symmetry of the equations of an echo state network and how a bias breaks it"),
                    ("kubota2021", "the capacities as coefficients of an expansion of the measured signals")],
        "textbooks": [],
    },
    "nonlinear-capacity": {
        "defined": ("Computation 1, §2", "with fading memory the capacities of all degrees sum to the number of "
                                         "independent measured signals", C1 + "#2-the-measure"),
        "derived": [("Computation 1, §4", "exact capacities of every degree and delay for one or two state variables",
                     C1 + "#4-running-the-command"),
                    ("Computation 1, §7", "published bodies from six fields", C1 + "#7-the-published-bodies")],
        "code": "computation.exact.ExactCapacities, card.card, ipc.capacities",
        "sources": [("boyd1985", "fading memory and the approximation of time-invariant filters by Volterra series"),
                    ("dambre2012", "the total capacity and its distribution over degrees"),
                    ("grigoryeva2015", "closed-form capacities of linear and quadratic tasks for delay systems"),
                    ("koster2024", "the linear capacity at every delay from the linearization of a delay system"),
                    ("hulser2023", "the error on a named task from the capacity profile"),
                    ("hodgkin1952", "the membrane of the squid giant axon"),
                    ("huang1996", "ultrasensitivity of the MAPK cascade"),
                    ("hauser2011", "networks of nonlinear springs and masses as computing bodies"),
                    ("furuta2018", "a macrospin tunnel junction measured at many times per input"),
                    ("torrejon2017", "computing with a spin-torque oscillator in experiment")],
        "textbooks": [],
    },
    "integrating": {
        "defined": ("Computation 1, §5", "a mode without decay: no capacity at a fixed delay", C1 + "#5-controls"),
        "derived": [("Computation 1, §5", "a rate set to zero turns the class", C1 + "#5-controls")],
        "code": "computation.predict.linear",
        "sources": [("jaeger2001", "the echo state property: the state forgets its initial condition"),
                    ("boyd1985", "fading memory")],
        "textbooks": [],
    },
    "inherited-through-threshold": {
        "defined": ("Heredity 1, §2", "a daughter below the threshold keeps its parent's sign with P = Φ(φ_c/σ_c)",
                    H1 + "#2-the-law"),
        "derived": [("Heredity 1, §4", "the law for a daughter, from the body's own equations",
                     H1 + "#4-running-the-command"),
                    ("Heredity 1, §6", "lineages of the full body in five fields against the law", H1 + "#6-lineages"),
                    ("Heredity 1, §7", "the published bodies", H1 + "#7-the-published-bodies")],
        "code": "heredity.predict.predict, lineage.condition, law.probability",
        "sources": [("lythe1996", "a slow passage from an initial offset: the Gaussian statistics of the order"),
                    ("vandenbroeck1987", "delayed bifurcations in the presence of noise"),
                    ("kondepudi1983", "the sign selected by a swept pitchfork with a bias"),
                    ("saito2007", "chiral autocatalysis with a critical number of molecules"),
                    ("painter1999", "Turing patterns on a growing domain"),
                    ("voituriez2005", "the spontaneous flow transition of an active polar gel above a critical thickness"),
                    ("duclos2018", "the spontaneous shear flow of a confined cellular nematic above a critical width"),
                    ("kretschmer1979", "the critical thickness of a ferroelectric film with surface terms"),
                    ("tilley1984", "the Landau theory of a ferroelectric film of finite thickness"),
                    ("baczynski2007", "the critical length of a filament under a fixed compressive load"),
                    ("hallatschek2007", "the overdamped dynamics of a semiflexible filament with thermal noise")],
        "textbooks": [],
    },
    "kept-above-threshold": {
        "defined": ("Heredity 1, §3", "the daughter is born above the threshold: no dip",
                    H1 + "#3-the-prediction-from-structure"),
        "derived": [("Heredity 1, §5", "division above twice the threshold", H1 + "#5-controls")],
        "code": "heredity.predict.predict",
        "sources": [("saito2007", "chiral autocatalysis above its critical number of molecules")],
        "textbooks": [],
    },
    "lost-in-the-dip": {
        "defined": ("Heredity 1, §3", "ln G = ∫λ dt over a generation < 0", H1 + "#3-the-prediction-from-structure"),
        "derived": [("Heredity 1, §5", "division closer to the threshold", H1 + "#5-controls")],
        "code": "heredity.lineage.ln_gain",
        "sources": [("baczynski2007", "the relaxation of a filament below its critical length")],
        "textbooks": [],
    },
    "threshold-moved": {
        "defined": ("Heredity 1, §3", "division moves the daughters away from the symmetric state of their size",
                    H1 + "#3-the-prediction-from-structure"),
        "derived": [("Heredity 1, §5", "a polarity of a conserved protein", H1 + "#5-controls")],
        "code": "heredity.predict.predict",
        "sources": [("brauns2020", "polarity by mass-conserving reaction and diffusion: the window of densities")],
        "textbooks": [],
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


# ------------------------------------------------------------------------------------------------ the quantum write
BIB.update({
    "marthaler2006": {"cite": "M. Marthaler and M. I. Dykman, Phys. Rev. A 73, 042108 (2006)", "short": "Marthaler and Dykman, 2006",
                      "doi": "10.1103/PhysRevA.73.042108"},
    "roquescarmes2023": {"cite": "C. Roques-Carmes et al., Science 381, 205 (2023)", "short": "Roques-Carmes et al., 2023",
                         "doi": "10.1126/science.adh4920"},
    "frattini2024": {"cite": "N. E. Frattini et al., Phys. Rev. X 14, 031040 (2024)", "short": "Frattini et al., 2024",
                     "doi": "10.1103/PhysRevX.14.031040"},
    "gu2025": {"cite": "A. Gu et al., Phys. Rev. Research 7, L022056 (2025)", "short": "Gu et al., 2025",
               "doi": "10.1103/PhysRevResearch.7.L022056"},
    "yamaji2025": {"cite": "T. Yamaji et al., Phys. Rev. Applied 24, 064052 (2025)", "short": "Yamaji et al., 2025",
                   "doi": "10.1103/PhysRevApplied.24.064052"},
    "grimm2020": {"cite": "A. Grimm et al., Nature 584, 205 (2020)", "short": "Grimm et al., 2020",
                  "doi": "10.1038/s41586-020-2587-z"},
})
QW = T + "30_quantum_write.md"
READING.update({
    "linear-stage-write": {
        "defined": ("Quantum write, §2", "the Gaussian stage of a swept parametric oscillator: the amplified quadrature "
                                         "obeys the classical write equation with the noise fixed by the loss and the "
                                         "temperature", QW + "#2-the-mechanism"),
        "derived": [("Quantum write, §3", "the 27-photon oscillator: the exact evolution of the density operator against "
                                          "the law along the protocol", QW + "#3-the-worked-example"),
                    ("Quantum write, §4", "the thermal control and the crossover between the dissipative and the closed "
                                          "limits", QW + "#4-the-results"),
                    ("Quantum write, §6", "two devices with published parameters, the Josephson oscillator of Yamaji et "
                                          "al. and the Kerr-cat resonator of Grimm et al.: both leave the law by 0.04 "
                                          "with frozen wells, the choice made in the nonlinear stage",
                     QW + "#6-measured-oscillators")],
        "code": "quantum.open.derive_quantum_write, law_along_protocol, evolve, prob_positive",
        "sources": [("kondepudi1985", "the selection law of a pitchfork swept with a bias and noise, "
                                      "P = Φ(π^{1/4} h / (D^{1/2} r^{1/4}))"),
                    ("gu2025", "the bias-probability relation of a parametric oscillator for a step of the gain, with the "
                               "initial quantum state and the gain noise"),
                    ("roquescarmes2023", "bias pulses below one photon set the probabilities of the two states of an "
                                         "optical parametric oscillator"),
                    ("yamaji2025", "a Josephson oscillator with twelve photons, its pump raised with a one-photon bias, "
                                   "and the probability of the favoured state against the bias"),
                    ("grimm2020", "a Kerr-cat resonator with 2.6 photons per well, its squeezing drive raised along a "
                                  "tanh ramp; the parameters of the second device of §6")],
    },
    "equilibrium-write": {
        "defined": ("Quantum write, §5", "at a few stored photons the switching between the wells outruns the sweep and "
                                         "the probability relaxes to the selection of the biased steady state",
                    QW + "#5-the-control-calculation"),
        "derived": [("Quantum write, §5", "the half-photon oscillator: the exact probability against the biased steady "
                                          "state and the gap of the Lindbladian", QW + "#5-the-control-calculation")],
        "code": "quantum.open.equilibrium, derive_quantum_write",
        "sources": [("marthaler2006", "switching between the two states of a parametrically modulated oscillator by "
                                      "quantum activation"),
                    ("frattini2024", "a few-photon Kerr parametric oscillator and the quantum regime of the Arrhenius law")],
    },
})
