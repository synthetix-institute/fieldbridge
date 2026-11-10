"""The sweep of the symmetric-write law check must stay on the passage the law describes.

Two corpus models of the kinetic Ising antiferromagnet (staggered magnetization m_s or s, uniform magnetization m, the
temperature T_b in 1/T_b) gave constants far below pi^(1/4) in the records: 0.30 +- 0.10 (arXiv:2502.01758, swept in
T_b) and 0.62 +- 0.10 (arXiv:2204.03995, swept in h through a cusp). Cooled through the Neel point T_b = 1 in an
independent simulation, the first gives 1.31 +- 0.02.

The check sweeps for 4/sqrt(rate) on each side of the crossing, which at the growth slope 1/2 of these models takes
T_b from 2.13 down to -0.13: through T_b = 0, where 1/T_b changes sign, the disordered state regains its stability and
the choice is erased. sweep_window holds the control where the rates at the critical state are no longer resolved by
the time step. In the second model the refined cusp lies 1e-6 from the symmetric point h = 0, the bias gradient left
over from that distance passed the test for a vanishing gradient, and its normal, along h, is tangent to the critical
line: the critical mode never grew along the sweep (growth slope 8e-7). cancel_asymmetry now sweeps along the growth
gradient when the zero-bias direction does not cross the critical line."""
from __future__ import annotations

import numpy as np

from fieldbridge.memory import codiscovery as cd
from fieldbridge.memory import construct as cs
from fieldbridge.memory import spec as mspec


def _ising(name, drift, control, value, rng):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": name, "field": "cond-mat",
            "question": "test", "assumptions": [], "provenance": {"source": name, "contributor": "test"},
            "carrier": {"kind": "euclid", "name": "magnetizations", "variables": list(drift), "scale": 1.0},
            "parameters": value, "drift": drift, "noise": 0.01, "dt": 0.01,
            "control": {"name": control, "range": rng}, "material": [k for k in value if k != control]}


AFM_TEMPERATURE = _ising(
    "arXiv:2502.01758 model 0 (test)",
    {"m_s": "-m_s/2 - tanh((h - m - m_s)/T_b)/4 + tanh((h - m + m_s)/T_b)/4",
     "m": "-1*m/2 + tanh((h - m - m_s)/T_b)/4 + tanh((h - m + m_s)/T_b)/4"},
    "T_b", {"T_b": 1.0, "h": 0.0}, [0.5, 1.5])
AFM_FIELD = _ising(
    "arXiv:2204.03995 model 0 (test)",
    {"s": "-1*s/2 - tanh((-s + (h - m))/T_b)/4 + tanh((s + (h - m))/T_b)/4",
     "m": "-1*m/2 + tanh((-s + (h - m))/T_b)/4 - tanh((s + (h - m))/T_b)/4"},
    "h", {"T_b": 0.2, "h": 0.5}, [0.25, 0.75])


def _law_check_input(real):
    seen = {}

    def spy(r, nf, rng, **kw):
        seen["real"], seen["nf"] = r, nf
        return {}

    original = cd.swept_write_check
    cd.swept_write_check = spy
    try:
        res = cd.derive_symmetric_write(real, np.random.default_rng(0), check_law=True, n_traj=10)
    finally:
        cd.swept_write_check = original
    return res, seen["real"], seen["nf"]


def _window(real, nf, rate=0.02):
    t_half = cs.sweep_half_time(real, nf, rate)
    dv = np.sign(-nf["kappa_slope"]) * rate / abs(nf["kappa_slope"])
    v0 = nf["value"] - dv * t_half
    return cs.sweep_window(real, nf, lambda t: v0 + dv * t, t_half), (lambda t: v0 + dv * t), t_half


def test_a_temperature_is_not_swept_through_zero():
    real = mspec.load(AFM_TEMPERATURE)
    _, r, nf = _law_check_input(real)
    window, path, t_half = _window(r, nf)
    assert window is not None
    assert path(2 * t_half) < 0.0                  # the whole sweep would end below T_b = 0
    assert 0.0 < path(window[1]) < 0.05            # it is held just above it, where 1/(2 T_b) reaches 0.3/dt
    assert window[0] == 0.0                        # the start, T_b = 2.13, is regular


def test_the_law_constant_of_the_antiferromagnet_cooled_through_its_neel_point():
    res = cd.derive_symmetric_write(mspec.load(AFM_TEMPERATURE), np.random.default_rng(0), check_law=True,
                                    n_traj=400)
    lc = res["law_constant"]
    assert res["class"] == "S(reflection) C R(pitchfork) K L"
    assert abs(lc["constant"] - np.pi ** 0.25) < 3 * lc["stderr"]      # 0.30 +- 0.10 before the window


def test_the_sweep_through_a_symmetric_cusp_crosses_the_critical_line():
    res, r, nf = _law_check_input(mspec.load(AFM_FIELD))
    U = next(s for s in res["steps"] if s["letter"] == "U")
    d = U["sweep_direction"]
    assert abs(d["T_b"]) > 0.99 and abs(d["h"]) < 1e-3        # across the Neel line, not along it
    assert abs(nf["kappa_slope"]) > 0.4                        # 7.8e-7 along h
    window, path, t_half = _window(r, nf)
    assert window is not None and window[1] < 2 * t_half       # cooled along T_b, held above T_b = 0


def test_the_law_constant_of_the_antiferromagnet_through_its_cusp():
    res = cd.derive_symmetric_write(mspec.load(AFM_FIELD), np.random.default_rng(0), check_law=True, n_traj=400)
    lc = res["law_constant"]
    assert res["class"] == "C R(transcritical) U R(pitchfork) K L"
    assert abs(lc["constant"] - np.pi ** 0.25) < 3 * lc["stderr"]      # 0.62 +- 0.10 before the two changes


def test_a_regular_sweep_is_not_held():
    from fieldbridge import site_data as sd
    spec = sd.resolve(only=["pitchfork"])["pitchfork"]["spec"]
    _, r, nf = _law_check_input(mspec.load(spec))
    assert _window(r, nf)[0] is None
