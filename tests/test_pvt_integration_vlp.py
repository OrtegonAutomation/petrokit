import numpy as np
from petrokit.vlp import vlp_curve_beggs_brill_blackoil

def test_vlp_blackoil_wrapper_runs_and_shapes():
    q = np.array([0, 200, 500], dtype=float)
    pwf = vlp_curve_beggs_brill_blackoil(
        q_range=q,
        well_depth=8000,
        d_in=2.875,
        t_f=180,
        api=35,
        gamma_g=0.8,
        mu_l_cp=1.5,
        p_ref_psia=2000,
        q_gas_mscf_d=500,
        theta_deg=90,
        p_wh=100,
    )
    assert pwf.shape == q.shape
    assert np.all(np.isfinite(pwf))
