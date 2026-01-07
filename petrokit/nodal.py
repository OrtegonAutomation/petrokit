import numpy as np
import matplotlib.pyplot as plt

from petrokit.ipr import vogel_ipr
from petrokit.vlp import vlp_curve
from petrokit.ipr_advanced import standing_ipr, jones_ipr
from petrokit.vlp_advanced import vlp_curve_beggs, vlp_curve_hagedorn


def find_operating_point(pwf_ipr, q_ipr, pwf_vlp, q_vlp, npts=100):
    """
    Encuentra el punto de operación (intersección entre curvas IPR y VLP)
    mediante interpolación lineal.

    Parámetros:
    -----------
    pwf_ipr : array
        Presiones correspondientes a la curva IPR [psi].
    q_ipr : array
        Caudales correspondientes a la curva IPR [STB/d].
    pwf_vlp : array
        Presiones correspondientes a la curva VLP [psi].
    q_vlp : array
        Caudales correspondientes a la curva VLP [STB/d].
    npts : int
        Número de puntos de interpolación.

    Retorna:
    --------
    (q_op, pwf_op) : tuple
        Punto de operación (caudal y presión de fondo fluyente).
    """
    from numpy import interp

    pwf_common = np.linspace(0, min(pwf_ipr[-1], pwf_vlp[-1]), npts)
    q_ipr_interp = interp(pwf_common, pwf_ipr, q_ipr)
    q_vlp_interp = interp(pwf_common, pwf_vlp, q_vlp)

    diff = np.abs(q_ipr_interp - q_vlp_interp)
    idx = np.argmin(diff)

    q_op = q_ipr_interp[idx]
    pwf_op = pwf_common[idx]

    return q_op, pwf_op


def nodal_analysis(p_res, q_max, model_ipr='vogel', model_vlp='simple',
                   well_depth=8000, rho=55, mu=2, d=2.5, J=2, n=1.8, s=0,
                   npts=50):
    """
    Realiza análisis nodal genérico cruzando IPR y VLP según modelos seleccionados.

    Parámetros:
    -----------
    p_res : float
        Presión de yacimiento [psi]
    q_max : float
        Caudal máximo a pwf = 0 [STB/d] (para modelos tipo Vogel / Standing)
    model_ipr : str
        Modelo IPR a usar ('vogel', 'standing', 'jones')
    model_vlp : str
        Modelo VLP a usar ('simple', 'beggs', 'hagedorn')
    well_depth : float
        Profundidad del pozo [ft]
    rho : float
        Densidad de fluido [lb/ft³]
    mu : float
        Viscosidad [cP]
    d : float
        Diámetro interno del tubing [in]
    J : float
        Índice de productividad (para modelo Jones)
    n : float
        Exponente empírico (para modelo Standing)
    s : float
        Skin factor (para modelo Jones)
    npts : int
        Número de puntos para discretización

    Retorna:
    --------
    (q_op, pwf_op) : tuple
        Punto de operación (caudal y presión de fondo fluyente).
    """

    # --- IPR Curve ---
    pwf_range = np.linspace(0, p_res, npts)

    if model_ipr == 'vogel':
        q_ipr = [vogel_ipr(p_res, q_max, pwf) for pwf in pwf_range]
    elif model_ipr == 'standing':
        q_ipr = [standing_ipr(p_res, q_max, pwf, n=n) for pwf in pwf_range]
    elif model_ipr == 'jones':
        q_ipr = [jones_ipr(p_res, J, pwf, s=s) for pwf in pwf_range]
    else:
        raise ValueError("Modelo IPR no reconocido. Usa 'vogel', 'standing' o 'jones'.")

    # --- VLP Curve ---
    q_range = np.linspace(100, q_max * 1.5, npts)

    if model_vlp == 'simple':
        pwf_vlp = vlp_curve(q_range, well_depth, rho, mu, d)
    elif model_vlp == 'beggs':
        pwf_vlp = vlp_curve_beggs(q_range, well_depth, rho, mu, d)
    elif model_vlp == 'hagedorn':
        pwf_vlp = vlp_curve_hagedorn(q_range, well_depth, rho, mu, d)
    else:
        raise ValueError("Modelo VLP no reconocido. Usa 'simple', 'beggs' o 'hagedorn'.")

    # --- Operating Point ---
    q_op, pwf_op = find_operating_point(pwf_range, q_ipr, pwf_vlp, q_range, npts)

    return q_op, pwf_op


def plot_nodal(p_res, q_max, model_ipr='vogel', model_vlp='simple',
               well_depth=8000, rho=55, mu=2, d=2.5, J=2, n=1.8, s=0,
               npts=50, ax=None):
    """
    Grafica las curvas IPR y VLP seleccionadas mostrando el punto de operación.
    """

    pwf_range = np.linspace(0, p_res, npts)
    q_range = np.linspace(100, q_max * 1.5, npts)

    # IPR
    if model_ipr == 'vogel':
        q_ipr = [vogel_ipr(p_res, q_max, pwf) for pwf in pwf_range]
        ipr_label = "Vogel"
    elif model_ipr == 'standing':
        q_ipr = [standing_ipr(p_res, q_max, pwf, n=n) for pwf in pwf_range]
        ipr_label = f"Standing (n={n})"
    elif model_ipr == 'jones':
        q_ipr = [jones_ipr(p_res, J, pwf, s=s) for pwf in pwf_range]
        ipr_label = f"Jones (s={s})"
    else:
        raise ValueError("Modelo IPR inválido")

    # VLP
    if model_vlp == 'simple':
        pwf_vlp = vlp_curve(q_range, well_depth, rho, mu, d)
        vlp_label = "VLP simple"
    elif model_vlp == 'beggs':
        pwf_vlp = vlp_curve_beggs(q_range, well_depth, rho, mu, d)
        vlp_label = "Beggs & Brill"
    elif model_vlp == 'hagedorn':
        pwf_vlp = vlp_curve_hagedorn(q_range, well_depth, rho, mu, d)
        vlp_label = "Hagedorn & Brown"
    else:
        raise ValueError("Modelo VLP inválido")

    q_op, pwf_op = find_operating_point(pwf_range, q_ipr, pwf_vlp, q_range, npts)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))

    ax.plot(q_ipr, pwf_range, label=f"IPR ({ipr_label})", color="blue", linewidth=2)
    ax.plot(q_range, pwf_vlp, label=f"VLP ({vlp_label})", color="red", linewidth=2)
    ax.scatter(q_op, pwf_op, color="green", s=80,
               label=f"Operación\nQ={q_op:.1f} STB/d\npwf={pwf_op:.1f} psi")
    ax.set_xlabel("Caudal [STB/d]")
    ax.set_ylabel("Presión Fondo Fluyente pwf [psi]")
    ax.set_title(f"Análisis Nodal ({ipr_label} + {vlp_label})")
    ax.legend()
    ax.grid(True)

    if ax is None:
        plt.show()

    return q_op, pwf_op
