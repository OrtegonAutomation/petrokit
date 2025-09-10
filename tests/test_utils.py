import math
from petrokit import utils

def test_pressure_conversions():
    psi = 1000
    pa = utils.psi_to_pa(psi)
    assert math.isclose(utils.pa_to_psi(pa), psi, rel_tol=1e-6)


def test_volume_conversions():
    stb = 100
    m3 = utils.stb_to_m3(stb)
    assert math.isclose(utils.m3_to_stb(m3), stb, rel_tol=1e-6)


def test_length_conversions():
    ft = 3280.84
    m = utils.ft_to_m(ft)
    assert math.isclose(utils.m_to_ft(m), ft, rel_tol=1e-6)


def test_reynolds_number_laminar():
    """
    Con valores bajos de velocidad y diámetro,
    el número de Reynolds debe ser menor a 2000 (flujo laminar).
    """
    Re = utils.reynolds_number(q=0.001, d=0.1, mu=10, rho=62.4)
    assert Re < 2000


def test_reynolds_number_turbulent():
    """
    Con valores altos de caudal y diámetro,
    el número de Reynolds debe ser mayor a 4000 (flujo turbulento).
    """
    Re = utils.reynolds_number(q=5, d=0.5, mu=1, rho=62.4)
    assert Re > 4000
