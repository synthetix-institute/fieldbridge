"""Memory constructor: build, predict and transfer memory on any carrier.

A material is described as a realization I_real = ((Omega, Xi); C, R, P; A) in a JSON specification
(``fieldbridge-memory/1``, see spec.py). The constructor then answers, for that material:

  predict   what its structure alone implies: whether it can store, frustrate or oscillate, how a state is
            written (pitchfork, fold, Hopf), which law loses a stored state, and which protocols are not limited by
            the relation between retention and writing times (predict.py, networks.py)
  card      the same questions calculated: stored states, kappa spectra, write points and their normal forms,
            barriers, retention and rewriting times and their ratio (analysis.py, construct.py, discovery.py)
  attach    which properties are kept when a memory is detached from one carrier and attached to another
  design    the material change that removes an obstruction (a cusp and the sweep through it)

Numerical routines need the memory extra: pip install -e '.[memory]'.
"""
from .carriers import Carrier, euclid, orthant, torus
from .identity import Realization

__all__ = ["Carrier", "Realization", "euclid", "orthant", "torus"]
