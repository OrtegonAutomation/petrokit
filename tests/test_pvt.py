import numpy as np
import pytest

from petrokit.pvt import (
    oil_rs_standing,
    oil_pb_standing,
    oil_bo_standing,
    gas_z_factor_papay,
    gas_bg_rb_per_scf,
    build_pvt_table,
)

def test_rs_increases_with_pressure_typical_range():
    p = np.array([500, 1000, 2000, 3000], dtype=float)
    rs = oil_rs_standing(p, t_f=180, api=35, gamma_g=0.8)
    assert np.all(np.diff(rs) >= 0)

def test_pb_and_rs_consistency_at_bubblepoint():
    rsb = 600.0
    pb = oil_pb_standing(rsb_scf_stb=rsb, t_f=180, api=35, gamma_g=0.8)
    rs_at_pb = float(oil_rs_standing(pb, t_f=180, api=35, gamma_g=0.8))
    # tras el fix, debe quedar bastante cercano
    assert rs_at_pb == pytest.approx(rsb, rel=0.05)

def test_bo_positive_and_increases_with_rs():
    rs = np.array([0, 200, 600, 1000], dtype=float)
    bo = oil_bo_standing(rs, t_f=180, api=35, gamma_g=0.8)
    assert np.all(bo > 0)
    assert np.all(np.diff(bo) >= 0)

def test_z_factor_bounds_due_to_clip():
    p = np.linspace(14.7, 8000, 50)
    z = gas_z_factor_papay(p, t_f=180, gamma_g=0.8)
    assert np.all(z >= 0.2)
    assert np.all(z <= 2.0)

def test_bg_positive_and_decreases_with_pressure_if_z_constant():
    p = np.array([500, 1000, 2000, 3000], dtype=float)
    z = np.ones_like(p)
    bg = gas_bg_rb_per_scf(p, t_f=180, z=z)
    assert np.all(bg > 0)
    assert np.all(np.diff(bg) < 0)

def test_build_pvt_table_shapes_and_keys():
    p = np.linspace(500, 3000, 10)
    out = build_pvt_table(p, t_f=180, api=35, gamma_g=0.8, z_method="papay")

    for k in ["P_psia", "Rs_scf_stb", "Bo_rb_stb", "Z", "Bg_rb_scf", "pb_psia", "rsb_scf_stb"]:
        assert k in out

    assert out["P_psia"].shape == (10,)
    assert out["Rs_scf_stb"].shape == (10,)
    assert out["Bo_rb_stb"].shape == (10,)
    assert out["Z"].shape == (10,)
    assert out["Bg_rb_scf"].shape == (10,)
    assert out["pb_psia"].shape == (1,)
    assert out["rsb_scf_stb"].shape == (1,)
