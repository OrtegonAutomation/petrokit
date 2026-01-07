import numpy as np
import pytest

from petrokit.ipr import jones_ipr, standing_ipr, ipr_curve_jones, ipr_curve_standing


def test_jones_endpoints_and_monotonicity():
    p_res = 3000.0
    C = 2.0
    D = 0.001

    # Endpoints
    assert jones_ipr(p_res, C, D, p_res) == 0.0
    q0 = jones_ipr(p_res, C, D, 0.0)
    assert q0 > 0.0

    # Monotonic: as pwf increases, q should not increase
    pwf = np.linspace(0.0, p_res, 50)
    q = jones_ipr(p_res, C, D, pwf)
    # q should be non-increasing with pwf
    assert np.all(np.diff(q) <= 1e-9)


def test_jones_negative_D_raises():
    with pytest.raises(ValueError):
        jones_ipr(3000.0, 2.0, -0.1, 1000.0)


def test_standing_is_continuous_at_pb_and_linear_above():
    p_res = 3000.0
    p_b = 2200.0
    J = 1.5

    # Continuity at bubble point
    q_at_pb = standing_ipr(p_res, p_b, J, p_b)
    q_linear_at_pb = J * (p_res - p_b)
    assert np.isclose(q_at_pb, q_linear_at_pb, rtol=1e-12, atol=1e-9)

    # Above Pb -> linear
    pwf_hi = 2500.0
    q_hi = standing_ipr(p_res, p_b, J, pwf_hi)
    assert np.isclose(q_hi, J * (p_res - pwf_hi), rtol=1e-12, atol=1e-9)

    # Below Pb -> should be >= the linear extrapolation at same pwf (typical behavior)
    # Below Pb -> debe aumentar respecto a q(Pb) y estar entre q(Pb) y q(0)
    pwf_lo = 1500.0
    q_lo = standing_ipr(p_res, p_b, J, pwf_lo)
    q_at_0 = standing_ipr(p_res, p_b, J, 0.0)

    assert q_lo > q_at_pb
    assert q_lo < q_at_0



def test_ipr_curves_return_shapes():
    pwf_j, q_j = ipr_curve_jones(3000.0, C=2.0, D=0.001, npts=31)
    assert pwf_j.shape == q_j.shape == (31,)
    assert np.isclose(pwf_j[0], 0.0)
    assert np.isclose(pwf_j[-1], 3000.0)

    pwf_s, q_s = ipr_curve_standing(3000.0, p_b=2200.0, J=1.5, npts=31)
    assert pwf_s.shape == q_s.shape == (31,)
    assert np.isclose(pwf_s[0], 0.0)
    assert np.isclose(pwf_s[-1], 3000.0)
