import numpy as np
import matplotlib.pyplot as plt

def flowline_pressure_drop(q, L, d, rho, mu, f=0.02, elev=0):
    """
    Calcula la pérdida de presión en una línea de flujo usando Darcy-Weisbach.

    Parámetros
    ----------
    q : float
        Caudal [STB/d]
    L : float
        Longitud de la línea [ft]
    d : float
        Diámetro interno de la tubería [in]
    rho : float
        Densidad del fluido [lb/ft3]
    mu : float
        Viscosidad dinámica [cp] (no usado por ahora, se deja para mejorar el modelo)
    f : float, opcional
        Factor de fricción (default=0.02)
    elev : float, opcional
        Diferencia de elevación entre extremos de la línea [ft]

    Retorna
    -------
    delta_p : float
        Pérdida de presión total [psi]
    """

    # Conversión de caudal a ft3/s
    stb_to_ft3 = 5.615 / (24*3600)  # 1 STB/d = 5.615 ft3/d
    q_ft3_s = q * stb_to_ft3

    # Área de la tubería [ft2]
    d_ft = d / 12
    A = np.pi * (d_ft/2)**2

    # Velocidad [ft/s]
    v = q_ft3_s / A

    # Constante
    g = 32.174  # ft/s2

    # Pérdida por fricción (psi)
    delta_p_fric = f * (L/d_ft) * (rho/144) * (v**2 / (2*g))

    # Pérdida/ganancia por elevación (psi)
    delta_p_elev = (rho/144) * elev

    # Total
    delta_p = delta_p_fric + delta_p_elev

    return delta_p


def plot_flowline(q_range, L, d, rho, mu, f=0.02, elev=0):
    """
    Grafica la caída de presión en la línea de flujo en función del caudal.
    """
    delta_p = [flowline_pressure_drop(q, L, d, rho, mu, f, elev) for q in q_range]

    import matplotlib.pyplot as plt
    plt.figure(figsize=(6,5))
    plt.plot(q_range, delta_p, color="green", label="Flowline ΔP")
    plt.xlabel("Caudal [STB/d]")
    plt.ylabel("ΔP [psi]")
    plt.title("Pérdida de Presión en Flowline")
    plt.grid(True)
    plt.legend()
    plt.show()

    return delta_p
