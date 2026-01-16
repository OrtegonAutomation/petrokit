import numpy as np
from petrokit.vlp import vlp_curve_model

def test_vlp_beggs_brill_blackoil_runs():
    q = np.array([0.0, 200.0, 500.0])

    pwf = vlp_curve_model(
        "beggs_brill_blackoil",
        q_range=q,
        well_depth=8000,
        rho=50.0,     # se ignora
        mu=1.5,       # mu_l
        d=2.875,
        # PVT inputs:
        t_f=180,
        api=35,
        gamma_g=0.8,
        p_ref_psia=2000,
        # VLP inputs:
        q_gas_mscf_d=500,
        theta_deg=90,
        p_wh=100,
    )

    assert pwf.shape == q.shape
    assert np.all(np.isfinite(pwf))
