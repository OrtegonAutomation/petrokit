import numpy as np
import matplotlib.pyplot as plt

def vogel_ipr(p_res, q_max, pwf):
    """
    Vogel IPR equation (for solution gas drive reservoirs)
    p_res: presión de yacimiento [psi]
    q_max: caudal máximo a pwf=0 [STB/d]
    pwf: presión de fondo fluyente [psi]
    """
    return q_max * (1 - 0.2 * (pwf / p_res) - 0.8 * (pwf / p_res) ** 2)

def fetkovich_ipr(p_res, J, pwf):
    """
    Fetkovich IPR equation (para reservorios lineales)
    J: índice de productividad [STB/d/psi]
    """
    return J * (p_res - pwf)

def ipr_curve_vogel(p_res, q_max, npts=30):
    pwf_range = np.linspace(0, p_res, npts)
    q_range = [vogel_ipr(p_res, q_max, pwf) for pwf in pwf_range]
    return pwf_range, q_range

def plot_ipr_vogel(p_res, q_max):
    pwf, q = ipr_curve_vogel(p_res, q_max)
    plt.figure(figsize=(6,5))
    plt.plot(q, pwf, label="Vogel IPR", color="blue")
    plt.xlabel("Caudal [STB/d]")
    plt.ylabel("Presión Fondo Fluyente pwf [psi]")
    plt.title("Curva IPR (Vogel)")
    plt.grid(True)
    plt.legend()
    plt.show()
