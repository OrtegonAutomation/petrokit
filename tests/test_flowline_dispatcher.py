import pytest
from petrokit.flowline import (
    available_flowline_models,
    flowline_pressure_drop,
    flowline_pressure_drop_model,
)

def test_available_flowline_models():
    models = available_flowline_models()
    assert "darcy" in models
    assert "beggs_brill" in models

def test_flowline_dispatcher_darcy_matches_legacy():
    dp_legacy = flowline_pressure_drop(q=1000, L=5000, d=6, rho=55, mu=1, f=0.02, elev=100)
    dp_disp = flowline_pressure_drop_model("darcy", q=1000, L=5000, d=6, rho=55, mu=1, f=0.02, elev=100)
    assert dp_disp == pytest.approx(dp_legacy, rel=1e-12, abs=1e-12)

def test_flowline_beggs_brill_increasing_with_rate():
    dp_low = flowline_pressure_drop_model("beggs_brill", q=500, L=5000, d=6, rho=55, mu=1, elev=0, q_gas_mscf_d=0.0)
    dp_high = flowline_pressure_drop_model("beggs_brill", q=5000, L=5000, d=6, rho=55, mu=1, elev=0, q_gas_mscf_d=0.0)
    assert dp_high > dp_low

def test_flowline_beggs_brill_elevation_effect_uphill():
    dp_flat = flowline_pressure_drop_model("beggs_brill", q=1000, L=5000, d=6, rho=55, mu=1, elev=0, q_gas_mscf_d=0.0)
    dp_up = flowline_pressure_drop_model("beggs_brill", q=1000, L=5000, d=6, rho=55, mu=1, elev=500, q_gas_mscf_d=0.0)
    assert dp_up > dp_flat

def test_flowline_beggs_brill_downhill_reduces_drop():
    dp_flat = flowline_pressure_drop_model("beggs_brill", q=1000, L=5000, d=6, rho=55, mu=1, elev=0, q_gas_mscf_d=0.0)
    dp_down = flowline_pressure_drop_model("beggs_brill", q=1000, L=5000, d=6, rho=55, mu=1, elev=-500, q_gas_mscf_d=0.0)
    assert dp_down < dp_flat  # puede incluso ser negativo, y está bien
