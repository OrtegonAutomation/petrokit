import numpy as np
from petrokit.flowline import flowline_pressure_drop

def test_flowline_positive():
    """
    La pérdida de presión en un flowline siempre debe ser positiva.
    """
    dp = flowline_pressure_drop(q=1000, L=5000, d=6, rho=55, mu=1)
    assert dp > 0


def test_flowline_increasing_with_rate():
    """
    A mayor caudal, la pérdida de presión debe ser mayor.
    """
    dp_low = flowline_pressure_drop(q=500, L=5000, d=6, rho=55, mu=1)
    dp_high = flowline_pressure_drop(q=5000, L=5000, d=6, rho=55, mu=1)
    assert dp_high > dp_low


def test_flowline_elevation_effect():
    """
    Si hay elevación positiva, la pérdida debe aumentar.
    """
    dp_no_elev = flowline_pressure_drop(q=1000, L=5000, d=6, rho=55, mu=1, elev=0)
    dp_with_elev = flowline_pressure_drop(q=1000, L=5000, d=6, rho=55, mu=1, elev=500)
    assert dp_with_elev > dp_no_elev
