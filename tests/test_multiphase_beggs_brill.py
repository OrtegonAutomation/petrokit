import math
import pytest

from petrokit.multiphase import (
    beggs_brill_boundaries,
    beggs_brill_flow_regime,
    beggs_brill_liquid_holdup,
    beggs_brill_pressure_gradient,
)


def _velocities_for_fr(c_l: float, fr: float, d_ft: float) -> tuple[float, float]:
    g = 32.174
    v_m = math.sqrt(fr * g * d_ft)
    v_sl = c_l * v_m
    v_sg = max(v_m - v_sl, 0.0)
    return v_sl, v_sg


def test_flow_regime_map_basic_regions():
    c_l = 0.2
    L1, L2, L3, L4 = beggs_brill_boundaries(c_l)
    d_ft = 3.0 / 12.0

    # Segregado (Fr < L2)
    fr = 0.5 * L2
    v_sl, v_sg = _velocities_for_fr(c_l, fr, d_ft)
    regime, _, _ = beggs_brill_flow_regime(v_sl, v_sg, d_ft)
    assert regime == "segregated"

    # Transicion (L2 < Fr < L3)
    fr = 0.5 * (L2 + L3)
    v_sl, v_sg = _velocities_for_fr(c_l, fr, d_ft)
    regime, _, _ = beggs_brill_flow_regime(v_sl, v_sg, d_ft)
    assert regime == "transition"

    # Intermintente (Fr > L3 and <= L1) para CL<0.4
    fr = max(L3 * 1.1, min(L1 * 0.5, L1 - 1e-6))
    v_sl, v_sg = _velocities_for_fr(c_l, fr, d_ft)
    regime, _, _ = beggs_brill_flow_regime(v_sl, v_sg, d_ft)
    assert regime in ("intermittent", "distributed")  # depending on map edges

    # Muy Alto Fr deberia ser distribuido
    fr = 2.0 * L4
    v_sl, v_sg = _velocities_for_fr(c_l, fr, d_ft)
    regime, _, _ = beggs_brill_flow_regime(v_sl, v_sg, d_ft)
    assert regime == "distributed"


def test_holdup_bounds():
    d_ft = 3.0 / 12.0
    v_sl = 2.0
    v_sg = 3.0
    c_l = v_sl / (v_sl + v_sg)

    holdup, regime = beggs_brill_liquid_holdup(
        v_sl_ft_s=v_sl,
        v_sg_ft_s=v_sg,
        d_ft=d_ft,
        theta_deg=30.0,
        rho_l_lbm_ft3=55.0,
        sigma_dyn_cm=30.0,
    )
    assert holdup >= c_l
    assert holdup <= 1.0


def test_dp_dz_increases_with_inclination():
   
   #Algunos rates, mas inclinación hacia arriba -> mayor componente hidrostática -> mayor dp/dz
    v_sl, v_sg = 2.0, 3.0

    r0 = beggs_brill_pressure_gradient(
        v_sl_ft_s=v_sl,
        v_sg_ft_s=v_sg,
        d_in=2.992,
        theta_deg=0.0,
        rho_l_lbm_ft3=55.0,
        rho_g_lbm_ft3=3.0,
        mu_l_cp=1.0,
        mu_g_cp=0.02,
    )
    r60 = beggs_brill_pressure_gradient(
        v_sl_ft_s=v_sl,
        v_sg_ft_s=v_sg,
        d_in=2.992,
        theta_deg=60.0,
        rho_l_lbm_ft3=55.0,
        rho_g_lbm_ft3=3.0,
        mu_l_cp=1.0,
        mu_g_cp=0.02,
    )

    assert r60.dp_dz_psi_per_ft > r0.dp_dz_psi_per_ft


def test_single_phase_liquid_limits_to_holdup_1():
    # v_sg=0 => c_l ~ 1 y holdup debe clampearse a 1
    r = beggs_brill_pressure_gradient(
        v_sl_ft_s=2.0,
        v_sg_ft_s=0.0,
        d_in=2.992,
        theta_deg=30.0,
        rho_l_lbm_ft3=55.0,
        rho_g_lbm_ft3=0.1,  
        mu_l_cp=1.0,
        mu_g_cp=0.02,
    )
    assert r.holdup == pytest.approx(1.0, abs=1e-6)
    assert r.dp_dz_psi_per_ft > 0.0
