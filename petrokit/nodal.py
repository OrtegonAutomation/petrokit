import numpy as np
import matplotlib.pyplot as plt
from petrokit.ipr import vogel_ipr
from petrokit.vlp import vlp_curve

def nodal_analysis(p_res, q_max, well_depth, rho, mu, d, npts=50):
    """
    Realiza análisis nodal cruzando IPR (Vogel) y VLP (simplificada).
    
    Parámetros:
    -----------
    p_res : float
        Presión de yacimiento [psi]
    q_max : float
        Caudal máximo a pwf=0 (STB/d)
    well_depth : float
        Profundidad del pozo [ft]
    rho : float
        Densidad de fluido [lb/ft³]
    mu : float
        Viscosidad del fluido [cP]
    d : float
        Diámetro del tubing [in]
    npts : int
        Número de puntos para el cruce de curvas
    
    Retorna:
    --------
    (q_op, pwf_op) : tuple
        Punto de operación (caudal y presión de fondo)
    """

    # Rango de presiones de fondo
    pwf_range = np.linspace(0, p_res, npts)

    # Curva IPR (Vogel)
    q_ipr = [vogel_ipr(p_res, q_max, pwf) for pwf in pwf_range]

    # Curva VLP (simplificada con función existente)
    q_range = np.linspace(100, q_max * 1.2, npts)
    pwf_vlp = vlp_curve(q_range, well_depth, rho, mu, d)

    # Interpolación para encontrar cruce
    from numpy import interp
    pwf_common = np.linspace(0, p_res, npts)
    q_ipr_interp = interp(pwf_common, pwf_range, q_ipr)
    q_vlp_interp = interp(pwf_common, pwf_vlp, q_range)

    diff = np.abs(q_ipr_interp - q_vlp_interp)
    idx = np.argmin(diff)

    q_op = q_ipr_interp[idx]
    pwf_op = pwf_common[idx]

    return q_op, pwf_op


def plot_nodal(p_res, q_max, well_depth, rho, mu, d, npts=50):
    """
    Dibuja el gráfico del análisis nodal mostrando el punto de operación.
    """
    pwf_range = np.linspace(0, p_res, npts)
    q_ipr = [vogel_ipr(p_res, q_max, pwf) for pwf in pwf_range]

    q_range = np.linspace(100, q_max * 1.2, npts)
    pwf_vlp = vlp_curve(q_range, well_depth, rho, mu, d)

    q_op, pwf_op = nodal_analysis(p_res, q_max, well_depth, rho, mu, d, npts)

    plt.figure(figsize=(7,6))
    plt.plot(q_ipr, pwf_range, label="IPR (Vogel)", color="blue")
    plt.plot(q_range, pwf_vlp, label="VLP", color="red")
    plt.scatter(q_op, pwf_op, color="green", s=80, label=f"Operación\nQ={q_op:.1f}, pwf={pwf_op:.1f}")
    plt.xlabel("Caudal [STB/d]")
    plt.ylabel("Presión Fondo Fluyente pwf [psi]")
    plt.title("Análisis Nodal")
    plt.legend()
    plt.grid(True)
    plt.show()

    return q_op, pwf_op
