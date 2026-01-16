# ================================================
# IPR Avanzado - PetroKit 
# Modelos: Standing y Jones + sensibilidad
# ================================================

import numpy as np
import matplotlib.pyplot as plt


# ================================================
# MODELOS DE IPR
# ================================================

def standing_ipr(p_res, q_max, pwf, n=1.8):
    """
    Modelo IPR de Standing (reservorios con drive de gas en solución).
    
    Ecuación:
        q = q_max * (1 - (pwf / p_res)^n)

    Parámetros:
    -----------
    p_res : float
        Presión de yacimiento [psi]
    q_max : float
        Caudal máximo a pwf = 0 [STB/d]
    pwf : float o array
        Presión de fondo fluyente [psi]
    n : float
        Exponente empírico (1.6 – 2.0, default = 1.8)
    
    Retorna:
    --------
    q : float o array
        Caudal [STB/d]
    """
    if p_res <= 0 or q_max <= 0:
        raise ValueError("p_res y q_max deben ser positivos")
    if n <= 0:
        raise ValueError("El exponente n debe ser positivo")

    # Asegurar que siempre sea al menos un arreglo 1D
    pwf = np.atleast_1d(pwf).astype(float)

    q = q_max * (1 - (pwf / p_res) ** n)
    q = np.where(q < 0, 0, q)  # evita caudales negativos

    # Si el usuario pasó un escalar, devolver un escalar
    return float(q[0]) if q.size == 1 else q



def jones_ipr(p_res, J, pwf, s=0):
    """
    Modelo IPR de Jones (pozos fracturados, efecto de skin).
    
    Ecuación:
        q = J_eff * (p_res - pwf)
        con J_eff = J / (1 + s)
    
    Parámetros:
    -----------
    p_res : float
        Presión de yacimiento [psi]
    J : float
        Índice de productividad [STB/d/psi]
    pwf : float o array
        Presión de fondo fluyente [psi]
    s : float
        Factor de daño/estimulación 
        (s=0 pozo ideal, s>0 dañado, s<0 estimulado)
    
    Retorna:
    --------
    q : float o array
        Caudal [STB/d]
    """
    if p_res <= 0 or J <= 0:
        raise ValueError("p_res y J deben ser positivos")
    # Jones simplificado: no dejamos que 1+s sea 0 para evitar división por cero
    # s < -1 representa estimulación teórica infinita en este modelo simplificado
    denom = (1 + s)
    if denom <= 0:
        denom = 1e-3 # Valor mínimo para evitar indeterminación


    pwf = np.atleast_1d(pwf).astype(float)
    delta_p = p_res - pwf
    delta_p = np.where(delta_p < 0, 0, delta_p)

    q_eff = (J / denom) * delta_p
    return q_eff[0] if q_eff.size == 1 else q_eff


# ================================================
# FUNCIONES GENERALES PARA CURVAS
# ================================================

def ipr_curve(model, p_res, **kwargs):
    """
    Genera curvas IPR para cualquier modelo (Standing, Jones, etc.)
    
    Parámetros:
    -----------
    model : función
        Modelo de IPR a usar (ej. standing_ipr, jones_ipr)
    p_res : float
        Presión de yacimiento [psi]
    kwargs : parámetros adicionales para el modelo
    
    Retorna:
    --------
    pwf_range : np.array
        Rango de presiones pwf [psi]
    q_range : np.array
        Rango de caudales [STB/d]
    """
    pwf_range = np.linspace(0, p_res, 30)
    q_range = [model(p_res, pwf=p, **kwargs) for p in pwf_range]
    return pwf_range, np.array(q_range)


def plot_ipr(models, p_res, labels=None, kwargs_list=None):
    """
    Graficar múltiples curvas IPR en un mismo gráfico.

    Parámetros:
    -----------
    models : list of funciones
        Lista de modelos a graficar
    p_res : float
        Presión de yacimiento [psi]
    labels : list of str
        Etiquetas para cada modelo
    kwargs_list : list of dict
        Lista de parámetros para cada modelo
    """
    plt.figure(figsize=(7, 5))
    for i, model in enumerate(models):
        args = kwargs_list[i] if kwargs_list else {}
        pwf, q = ipr_curve(model, p_res, **args)
        plt.plot(q, pwf, linewidth=2, label=labels[i] if labels else f"Modelo {i+1}")
    
    plt.xlabel("Caudal q [STB/d]")
    plt.ylabel("Presión de fondo fluyente pwf [psi]")
    plt.gca().invert_yaxis()
    plt.legend()
    plt.title("Comparación de curvas IPR avanzadas")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.show()


# ================================================
# SENSIBILIDAD DE PARÁMETROS
# ================================================

def plot_standing_sensitivity(p_res, q_max, pwf_points=None, n_values=[1.6, 1.8, 2.0]):
    """
    Graficar sensibilidad del exponente n en el modelo de Standing.
    """
    plt.figure(figsize=(7, 5))
    for n in n_values:
        pwf, q = ipr_curve(standing_ipr, p_res, q_max=q_max, n=n)
        plt.plot(q, pwf, label=f"Standing (n={n})")

    if pwf_points:
        plt.scatter([0], [p_res], marker="o", color="black", label="p_res")

    plt.xlabel("Caudal q [STB/d]")
    plt.ylabel("Presión pwf [psi]")
    plt.title("Sensibilidad modelo Standing")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.legend()
    plt.gca().invert_yaxis()
    plt.show()


def plot_jones_sensitivity(p_res, J, s_values=[-0.5, 0, 2]):
    """
    Graficar sensibilidad del skin en el modelo de Jones.
    """
    plt.figure(figsize=(7, 5))
    for s in s_values:
        pwf, q = ipr_curve(jones_ipr, p_res, J=J, s=s)
        plt.plot(q, pwf, label=f"Jones (s={s})")

    plt.xlabel("Caudal q [STB/d]")
    plt.ylabel("Presión pwf [psi]")
    plt.title("Sensibilidad modelo Jones (skin)")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.legend()
    plt.gca().invert_yaxis()
    plt.show()


# ================================================
# EJEMPLO DE USO
# ================================================
if __name__ == "__main__":
    p_res = 3000  # psi
    q_max = 1500  # STB/d
    J = 2.0       # STB/d/psi

    # Comparación Standing vs Jones
    plot_ipr(
        models=[standing_ipr, jones_ipr],
        p_res=p_res,
        labels=["Standing (n=1.8)", "Jones (s=0)"],
        kwargs_list=[{"q_max": q_max, "n": 1.8}, {"J": J, "s": 0}]
    )

    # Sensibilidad Standing
    plot_standing_sensitivity(p_res=p_res, q_max=q_max)

    # Sensibilidad Jones
    plot_jones_sensitivity(p_res=p_res, J=J)


__all__ = [
    "standing_ipr",
    "jones_ipr",
    "ipr_curve",
    "plot_ipr",
    "plot_standing_sensitivity",
    "plot_jones_sensitivity",
]
