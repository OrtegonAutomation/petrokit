import numpy as np
from petrokit.vlp import vlp_curve_model


def test_vlp_curve_model_hagedorn_brown_runs_and_monotonic():
    q = np.linspace(0, 1200, 25)
    pwf = vlp_curve_model(
        "hagedorn_brown",
        q_range=q,
        well_depth=8000,
        rho=60,
        mu=1.0,
        d=2.992,
        q_gas_mscf_d=200.0,
    )
    assert len(pwf) == len(q)
    assert np.all(np.isfinite(pwf))

    dif = np.diff(pwf)
    assert np.all(dif >= -1e-8)
