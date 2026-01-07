import numpy as np
from petrokit.ipr_advanced import standing_ipr, jones_ipr

def test_standing_ipr_zero():
    assert standing_ipr(p_res=3000, q_max=1500, pwf=3000) == 0

def test_standing_ipr_max():
    assert standing_ipr(p_res=3000, q_max=1500, pwf=0) == 1500

def test_jones_ipr_linear():
    q = jones_ipr(p_res=3000, J=2, pwf=2000, s=0)
    assert np.isclose(q, 2000)  # (3000-2000)*2

def test_jones_ipr_skin_effect():
    q_no_skin = jones_ipr(p_res=3000, J=2, pwf=2000, s=0)
    q_stimulated = jones_ipr(p_res=3000, J=2, pwf=2000, s=-0.5)
    assert q_stimulated > q_no_skin
