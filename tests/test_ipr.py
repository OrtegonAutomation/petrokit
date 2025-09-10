import pytest
from petrokit.ipr import vogel_ipr, fetkovich_ipr

def test_vogel_ipr_zero_pwf():
    # Si pwf = 0 → Vogel debe devolver q_max
    q = vogel_ipr(p_res=3000, q_max=1000, pwf=0)
    assert pytest.approx(q, rel=1e-3) == 1000

def test_vogel_ipr_reservoir_pressure():
    # Si pwf = p_res → q = 0
    q = vogel_ipr(p_res=3000, q_max=1000, pwf=3000)
    assert pytest.approx(q, rel=1e-3) == 0

def test_fetkovich_ipr_linear():
    # Fetkovich es lineal
    q = fetkovich_ipr(p_res=2500, J=2, pwf=2500)
    assert q == 0
    q = fetkovich_ipr(p_res=2500, J=2, pwf=2000)
    assert q == 1000
