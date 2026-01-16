import numpy as np
from petrokit.vlp_advanced import beggs_brill_dp, hagedorn_brown_dp, vlp_curve_beggs, vlp_curve_hagedorn

def test_beggs_brill_positive_and_holdup():
    dp, H = beggs_brill_dp(q_liq_stb_d=500, q_gas_mscf_d=200, d_in=3.5, L_ft=8000, theta_deg=90,
                          rho_l=60, rho_g=0.08, mu_l=1.0, mu_g=0.02)
    assert dp >= 0
    assert 0 < H < 1

def test_hagedorn_brown_positive_and_holdup():
    dp, H = hagedorn_brown_dp(q_liq_stb_d=500, q_gas_mscf_d=200, d_in=3.5, L_ft=8000, theta_deg=90,
                              rho_l=60, rho_g=0.08, mu_l=1.0, mu_g=0.02)
    assert dp >= 0
    assert 0 < H < 1

def test_vlp_curves_lengths():
    q_range = np.linspace(100, 1500, 10)
    q, pwf_bb, H_bb = vlp_curve_beggs(q_range, 200, 3.5, 8000, 90, 60, 0.08, 1.0, 0.02)
    assert pwf_bb.shape == q.shape
    assert np.all(pwf_bb >= 0)
    q2, pwf_hb, H_hb = vlp_curve_hagedorn(q_range, 200, 3.5, 8000, 90, 60, 0.08, 1.0, 0.02)
    assert pwf_hb.shape == q2.shape
    assert np.all(pwf_hb >= 0)
