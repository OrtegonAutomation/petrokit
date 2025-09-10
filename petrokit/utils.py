import math

def psi_to_pa(psi):
    """
    Convierte presión de psi a Pascal.
    1 psi = 6894.76 Pa
    """
    return psi * 6894.76


def pa_to_psi(pa):
    """
    Convierte presión de Pascal a psi.
    """
    return pa / 6894.76


def stb_to_m3(stb):
    """
    Convierte barriles estándar (STB) a metros cúbicos.
    1 STB = 0.158987 m³
    """
    return stb * 0.158987


def m3_to_stb(m3):
    """
    Convierte metros cúbicos a barriles estándar (STB).
    """
    return m3 / 0.158987


def ft_to_m(ft):
    """
    Convierte pies a metros.
    """
    return ft * 0.3048


def m_to_ft(m):
    """
    Convierte metros a pies.
    """
    return m / 0.3048


def reynolds_number(q, d, mu, rho):
    """
    Calcula el número de Reynolds en tuberías.
    
    Parámetros:
    -----------
    q : float
        Caudal volumétrico [ft³/s]
    d : float
        Diámetro interno [ft]
    mu : float
        Viscosidad dinámica [cP]
    rho : float
        Densidad del fluido [lb/ft³]
    
    Retorna:
    --------
    Re : float
        Número de Reynolds (adimensional)
    """
    # Convertir viscosidad a lb/ft-s (1 cP = 0.000672 lb/ft-s)
    mu_lbfts = mu * 0.000672
    area = math.pi * (d / 2) ** 2
    v = q / area  # velocidad [ft/s]
    Re = (rho * v * d) / mu_lbfts
    return Re
