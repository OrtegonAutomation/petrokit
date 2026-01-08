import numpy as np
import pytest

from petrokit.multiphase_hb import hagedorn_brown_pressure_gradient


def test_hb_core_outputs_are_finite_and_bounded():
    res = hagedorn_brown_pressure_gradient(
        v_sl_ft_s=1.0,
        v_sg_ft_s=2.0,
        d_in=2.992,
        theta_deg=90.0,
        rho_l_lbm_ft3=60.0,
        rho_g_lbm_ft3=2.0,
        mu_l_cp=1.0,
        mu_g_cp=0.02,
    )
    assert np.isfinite(res.dp_dz_psi_per_ft)
    assert res.dp_dz_psi_per_ft > 0.0
    assert 0.0 < res.holdup < 1.0


def test_hb_angle_effect_small_angles():
    base = dict(
        v_sl_ft_s=1.0,
        v_sg_ft_s=2.0,
        d_in=2.992,
        rho_l_lbm_ft3=60.0,
        rho_g_lbm_ft3=2.0,
        mu_l_cp=1.0,
        mu_g_cp=0.02,
    )
    down = hagedorn_brown_pressure_gradient(theta_deg=-5.0, **base).dp_dz_psi_per_ft
    horz = hagedorn_brown_pressure_gradient(theta_deg=0.0, **base).dp_dz_psi_per_ft
    up = hagedorn_brown_pressure_gradient(theta_deg=5.0, **base).dp_dz_psi_per_ft

    assert up > horz
    assert horz > down

    # En horizontal no hay término hidrostático -> queda fricción pura -> debe ser positiva
    assert horz > 0.0



@pytest.mark.parametrize("bad", [
    dict(d_in=0.0),
    dict(rho_l_lbm_ft3=-1.0),
    dict(mu_l_cp=0.0),
    dict(theta_deg=120.0),
    dict(v_sl_ft_s=-1.0),
])
def test_hb_invalid_inputs_raise(bad):
    base = dict(
        v_sl_ft_s=1.0,
        v_sg_ft_s=2.0,
        d_in=2.992,
        theta_deg=90.0,
        rho_l_lbm_ft3=60.0,
        rho_g_lbm_ft3=2.0,
        mu_l_cp=1.0,
        mu_g_cp=0.02,
    )
    base.update(bad)
    with pytest.raises(ValueError):
        hagedorn_brown_pressure_gradient(**base)
