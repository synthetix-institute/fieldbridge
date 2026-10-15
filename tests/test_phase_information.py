"""Diffusion on a circle and on its unwrapped covering line have different information laws."""
import numpy as np


def test_wrapped_phase_information_has_exponential_tail():
    angle = np.linspace(0, 2*np.pi, 32768, endpoint=False)
    information = []
    for dt in (4.0, 8.0):
        # Heat kernel divided by the uniform density, for a uniformly written phase.
        ratio = np.ones_like(angle)
        for n in range(1, 20):
            ratio += 2*np.exp(-n*n*dt)*np.cos(n*angle)
        assert abs(ratio.mean() - 1) < 1e-14
        information.append(np.mean(ratio*np.log(ratio)))
    assert np.isclose(information[-1], np.exp(-16), rtol=1e-6)
    assert np.isclose(information[-1]/information[0], np.exp(-8), rtol=2e-4)
    unwrapped = 0.5*np.log1p(1/(2*8.0))
    assert unwrapped > 10000*information[-1]
