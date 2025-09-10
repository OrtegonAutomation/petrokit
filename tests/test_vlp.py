import numpy as np
from petrokit.vlp import vlp_curve

def test_vlp_monotonic():
    """
    La curva VLP debe ser creciente con el caudal,
    ya que mayor caudal implica mayor caída de presión.
    """
    q_range = np.linspace(100, 2000, 10)
    pwf = vlp_curve(q_range, well_depth=8000, rho=60, mu=1, d=2.992)

    # La presión debe crecer con el caudal
    assert all(np.diff(pwf) > 0)


def test_vlp_positive_pressures():
    """
    Todas las presiones calculadas deben ser positivas.
    """
    q_range = np.linspace(100, 2000, 10)
    pwf = vlp_curve(q_range, well_depth=8000, rho=60, mu=1, d=2.992)

    assert all(p > 0 for p in pwf)


def test_vlp_low_rate_limit():
    """
    A caudal muy bajo, la presión debe ser cercana a la hidrostática.
    """
    q_range = np.array([0.0001])  # casi cero caudal
    pwf = vlp_curve(q_range, well_depth=8000, rho=60, mu=1, d=2.992)

    # Presión hidrostática esperada
    p_hydro = (60/144) * 8000  # psi
    assert abs(pwf[0] - p_hydro) < 5  # tolerancia de ±5 psi
