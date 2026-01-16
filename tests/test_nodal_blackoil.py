import math
from petrokit.nodal import nodal_analysis

def test_nodal_runs_with_beggs_brill_blackoil():
    q_star, pwf_star = nodal_analysis(
        p_res=3000,
        q_max=800,
        well_depth=8000,
        rho=50.0,
        mu=1.5,
        d=2.875,
        npts=10,  # >= 10 según nodal_analysis_detail
        ipr_model="vogel",
        vlp_model="beggs_brill_blackoil",
        ipr_kwargs={},
        vlp_kwargs={
            "t_f": 180,
            "api": 35,
            "gamma_g": 0.8,
            "p_ref_psia": 2000,
            "q_gas_mscf_d": 500,
            "theta_deg": 90,
            "p_wh": 100,
        },
    )

    assert isinstance(q_star, (int, float))
    assert isinstance(pwf_star, (int, float))
    assert math.isfinite(q_star)
    assert math.isfinite(pwf_star)
    assert 0 <= q_star <= 800
    assert 0 <= pwf_star <= 3000
