# Wanted materials

Each row is a model, well documented in its field, that the [catalog](materials.md) does not yet contain. Any of
them is a self-contained first contribution: write its equations as a specification
([CONTRIBUTING.md](../CONTRIBUTING.md)), check it, and open a pull request. The last column says what the model is
known for in its field; whether and how it stores a state is what the memory card calculates. Check the source
before writing the specification: the parameter values must be those of the paper, or be listed under
`assumptions`.

| Field | Model | Source | Known for |
| --- | --- | --- | --- |
| ferroelectrics | Landau–Khalatnikov dynamics of the polarization in a double-well free energy | L. D. Landau and I. M. Khalatnikov, Dokl. Akad. Nauk SSSR 96, 469 (1954) | polarization switching under an applied field |
| optics | Optical bistability of atoms in a ring cavity, mean-field equations | R. Bonifacio and L. A. Lugiato, Phys. Rev. A 18, 1129 (1978) | two transmitted intensities for one input intensity |
| mechanics | Driven Duffing oscillator, written as the averaged (slow-flow) equations of amplitude and phase | G. Duffing, *Erzwungene Schwingungen bei veränderlicher Eigenfrequenz* (Vieweg, 1918) | jumps between two response amplitudes when the drive frequency is swept |
| climate | Stommel's two-box model of the thermohaline circulation | H. Stommel, Tellus 13, 224 (1961) | two circulation states and hysteresis under freshwater forcing |
| climate | Budyko energy-balance model with ice–albedo feedback | M. I. Budyko, Tellus 21, 611 (1969) | an ice-covered and an ice-free climate for the same forcing |
| ecology | Spruce budworm outbreaks | D. Ludwig, D. D. Jones and C. S. Holling, J. Anim. Ecol. 47, 315 (1978) | an outbreak threshold and hysteresis of the population |
| cell biology | Lac operon induction | E. M. Ozbudak, M. Thattai, H. N. Lim, B. I. Shraiman and A. van Oudenaarden, Nature 427, 737 (2004) | bistable induction and memory of the inducer history |
| cell biology | Cdc2–cyclin B mitotic switch | J. R. Pomerening, E. D. Sontag and J. E. Ferrell, Nat. Cell Biol. 5, 346 (2003) | hysteresis of mitotic entry |
| neuroscience | Wilson–Cowan excitatory and inhibitory populations | H. R. Wilson and J. D. Cowan, Biophys. J. 12, 1 (1972) | coexisting activity states and oscillations |
| neuroscience | Morris–Lecar neuron | C. Morris and H. Lecar, Biophys. J. 35, 193 (1981) | onset of firing through different bifurcations depending on parameters |
| chemistry | Oregonator model of the Belousov–Zhabotinsky reaction | R. J. Field and R. M. Noyes, J. Chem. Phys. 60, 1877 (1974) | chemical oscillations and their entrainment |
| neural networks | Hopfield network of two graded-response units | J. J. Hopfield, Proc. Natl. Acad. Sci. USA 81, 3088 (1984) | stored patterns as stable states |
| statistical physics | Mean-field Glauber dynamics of an Ising magnet (Curie–Weiss) | R. J. Glauber, J. Math. Phys. 4, 294 (1963) | spontaneous magnetization below the critical temperature and hysteresis in a field |

To propose a model that is not listed, use the
[material proposal form](https://github.com/synthetix-institute/fieldbridge/issues/new?template=material.yml).
