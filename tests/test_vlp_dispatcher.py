import numpy as np

from petrokit.vlp import available_vlp_models, vlp_curve_model


def test_available_vlp_models_contains_expected():
    models = available_vlp_models()
    assert "darcy" in models
    assert "beggs_brill" in models
    assert "hagedorn_brown" in models


def test_vlp_curve_model_darcy_runs():
    q = np.linspace(0, 1200, 30)
    pwf = vlp_curve_model("darcy", q, well_depth=8000, rho=60, mu=1, d=2.992, f=0.02)
    assert len(pwf) == len(q)
    assert np.all(np.isfinite(pwf))


def test_vlp_curve_model_beggs_brill_runs_and_monotonic():
    q = np.linspace(0, 1200, 30)
    pwf = vlp_curve_model("beggs_brill", q, well_depth=8000, rho=60, mu=1, d=2.992)
    assert len(pwf) == len(q)
    assert np.all(np.isfinite(pwf))

    # Monótono no-decreciente (por tolerancia numérica)
    dif = np.diff(pwf)
    assert np.all(dif >= -1e-8)
