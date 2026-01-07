"""
petrokit.ipr
=============

Modelos IPR (Inflow Performance Relationship).

Unidades esperadas (consistentes con el README del proyecto):
- Presiones: psi
- Caudales (oil): STB/d

Incluye (Fase 1):
- Vogel (solution-gas drive, saturado)
- Fetkovich (lineal por PI)

Incluye (Fase 2 / Producción):
- Jones (turbulencia / no-Darcy): Δp = C*q + D*q^2
- Standing (undersaturated / compuesto): lineal sobre Pb + Vogel bajo Pb
"""

from __future__ import annotations

from typing import Tuple, Union
import numpy as np
import matplotlib.pyplot as plt

ArrayLike = Union[float, int, np.ndarray]


# ---------------------------
# Helpers
# ---------------------------
def _as_array(x: ArrayLike) -> np.ndarray:
    return np.asarray(x, dtype=float)


def _maybe_scalar(x_in: ArrayLike, x_out: np.ndarray) -> ArrayLike:
    # Si el input fue escalar, devuelve escalar. Si no, ndarray.
    if np.isscalar(x_in):
        return float(np.asarray(x_out).reshape(-1)[0])
    return x_out


def _require_positive(name: str, value: float, allow_zero: bool = False) -> None:
    if allow_zero:
        if not (value >= 0):
            raise ValueError(f"{name} debe ser >= 0. Recibido: {value}")
    else:
        if not (value > 0):
            raise ValueError(f"{name} debe ser > 0. Recibido: {value}")


def _clip_pwf(pwf: np.ndarray, p_res: float) -> np.ndarray:
    # Para robustez numérica y evitar valores no físicos
    return np.clip(pwf, 0.0, p_res)


# ---------------------------
# Modelos IPR (rates q = f(pwf))
# ---------------------------
def vogel_ipr(p_res: float, q_max: float, pwf: ArrayLike) -> ArrayLike:
    """
    Vogel IPR (solution-gas drive, saturado):

        q = q_max * (1 - 0.2*(pwf/p_res) - 0.8*(pwf/p_res)^2)

    Donde:
      - p_res: presión promedio de yacimiento (psi)
      - q_max: caudal a pwf=0 (STB/d)
      - pwf: presión de fondo fluyente (psi), escalar o vector

    Devuelve:
      - q (STB/d), mismo shape que pwf

    Nota:
      - pwf se “clipea” a [0, p_res] para evitar valores no físicos.
    """
    _require_positive("p_res", float(p_res))
    _require_positive("q_max", float(q_max), allow_zero=True)

    pwf_arr = _as_array(pwf)
    pwf_arr = _clip_pwf(pwf_arr, float(p_res))

    x = pwf_arr / float(p_res)
    q = float(q_max) * (1.0 - 0.2 * x - 0.8 * x**2)
    q = np.maximum(q, 0.0)

    return _maybe_scalar(pwf, q)


def fetkovich_ipr(p_res: float, J: float, pwf: ArrayLike) -> ArrayLike:
    """
    Fetkovich (forma lineal por PI):

        q = J * (p_res - pwf)

    Donde:
      - J: Productivity Index (STB/d/psi)

    Nota:
      - pwf se clipea a [0, p_res] y q se limita a >= 0.
    """
    _require_positive("p_res", float(p_res))
    _require_positive("J", float(J), allow_zero=False)

    pwf_arr = _as_array(pwf)
    pwf_arr = _clip_pwf(pwf_arr, float(p_res))

    q = float(J) * (float(p_res) - pwf_arr)
    q = np.maximum(q, 0.0)

    return _maybe_scalar(pwf, q)


def jones_ipr(p_res: float, C: float, D: float, pwf: ArrayLike) -> ArrayLike:
    """
    Jones (turbulencia / no-Darcy) usando la forma cuadrática:

        Δp = (p_res - pwf) = C*q + D*q^2

    => D*q^2 + C*q - Δp = 0

    Solución física (q >= 0):
      - si D > 0: q = (-C + sqrt(C^2 + 4*D*Δp)) / (2*D)
      - si D = 0: q = Δp / C

    Parámetros:
      - C: coef. laminar (psi / (STB/d))
      - D: coef. turbulento (psi / (STB/d)^2)

    Nota:
      - pwf se clipea a [0, p_res]
      - si D < 0 se considera no físico y se lanza error
    """
    _require_positive("p_res", float(p_res))
    _require_positive("C", float(C), allow_zero=False)

    if float(D) < 0:
        raise ValueError("D (coeficiente turbulento) no puede ser negativo.")

    pwf_arr = _as_array(pwf)
    pwf_arr = _clip_pwf(pwf_arr, float(p_res))
    dp = float(p_res) - pwf_arr  # Δp >= 0

    if float(D) == 0.0:
        q = dp / float(C)
        q = np.maximum(q, 0.0)
        return _maybe_scalar(pwf, q)

    disc = float(C) ** 2 + 4.0 * float(D) * dp
    disc = np.maximum(disc, 0.0)
    q = (-float(C) + np.sqrt(disc)) / (2.0 * float(D))
    q = np.maximum(q, 0.0)

    return _maybe_scalar(pwf, q)


def standing_ipr(p_res: float, p_b: float, J: float, pwf: ArrayLike) -> ArrayLike:
    """
    Standing (undersaturated / compuesto "gas solution"):
    - Para pwf >= Pb (flujo monofásico): q = J*(p_res - pwf)
    - Para pwf <  Pb (dos fases): Vogel, pero usando el mismo p_res
      y calibrando q_max para que sea continuo en Pb.

    Pasos:
      1) q_b = J*(p_res - p_b)      (rate en bubble point)
      2) q_max = q_b / [1 - 0.2*(p_b/p_res) - 0.8*(p_b/p_res)^2]
      3) si pwf < p_b: q = q_max*(1 - 0.2*(pwf/p_res) - 0.8*(pwf/p_res)^2)

    Requisitos:
      - p_res > p_b > 0
      - J > 0
    """
    _require_positive("p_res", float(p_res))
    _require_positive("p_b", float(p_b))
    _require_positive("J", float(J), allow_zero=False)

    if not (0.0 < float(p_b) < float(p_res)):
        raise ValueError(f"Se requiere 0 < p_b < p_res. Recibido p_b={p_b}, p_res={p_res}")

    pwf_arr = _as_array(pwf)
    pwf_arr = _clip_pwf(pwf_arr, float(p_res))

    q_linear = float(J) * (float(p_res) - pwf_arr)

    # Calibración en Pb
    q_b = float(J) * (float(p_res) - float(p_b))
    x_b = float(p_b) / float(p_res)
    denom = 1.0 - 0.2 * x_b - 0.8 * x_b**2
    if denom <= 0:
        raise ValueError("Parámetros no físicos: denom <= 0 al calibrar Standing/Vogel en Pb.")

    q_max = q_b / denom

    x = pwf_arr / float(p_res)
    q_vogel = q_max * (1.0 - 0.2 * x - 0.8 * x**2)

    q = np.where(pwf_arr >= float(p_b), q_linear, q_vogel)
    q = np.maximum(q, 0.0)

    return _maybe_scalar(pwf, q)


# ---------------------------
# Curvas + plots
# ---------------------------
def ipr_curve_vogel(p_res: float, q_max: float, npts: int = 50) -> Tuple[np.ndarray, np.ndarray]:
    _require_positive("p_res", float(p_res))
    _require_positive("q_max", float(q_max), allow_zero=True)
    if int(npts) < 2:
        raise ValueError("npts debe ser >= 2")

    pwf = np.linspace(0.0, float(p_res), int(npts))
    q = _as_array(vogel_ipr(p_res, q_max, pwf))
    return pwf, q


def ipr_curve_fetkovich(p_res: float, J: float, npts: int = 50) -> Tuple[np.ndarray, np.ndarray]:
    _require_positive("p_res", float(p_res))
    _require_positive("J", float(J), allow_zero=False)
    if int(npts) < 2:
        raise ValueError("npts debe ser >= 2")

    pwf = np.linspace(0.0, float(p_res), int(npts))
    q = _as_array(fetkovich_ipr(p_res, J, pwf))
    return pwf, q


def ipr_curve_jones(p_res: float, C: float, D: float, npts: int = 50) -> Tuple[np.ndarray, np.ndarray]:
    _require_positive("p_res", float(p_res))
    _require_positive("C", float(C), allow_zero=False)
    if float(D) < 0:
        raise ValueError("D no puede ser negativo.")
    if int(npts) < 2:
        raise ValueError("npts debe ser >= 2")

    pwf = np.linspace(0.0, float(p_res), int(npts))
    q = _as_array(jones_ipr(p_res, C, D, pwf))
    return pwf, q


def ipr_curve_standing(p_res: float, p_b: float, J: float, npts: int = 50) -> Tuple[np.ndarray, np.ndarray]:
    _require_positive("p_res", float(p_res))
    _require_positive("p_b", float(p_b))
    _require_positive("J", float(J), allow_zero=False)
    if int(npts) < 2:
        raise ValueError("npts debe ser >= 2")

    pwf = np.linspace(0.0, float(p_res), int(npts))
    q = _as_array(standing_ipr(p_res, p_b, J, pwf))
    return pwf, q


def plot_ipr_vogel(p_res: float, q_max: float, npts: int = 50) -> None:
    pwf, q = ipr_curve_vogel(p_res, q_max, npts=npts)
    plt.figure()
    plt.plot(q, pwf)
    plt.xlabel("q (STB/d)")
    plt.ylabel("pwf (psi)")
    plt.title("IPR - Vogel")
    plt.grid(True)
    plt.show()


def plot_ipr_fetkovich(p_res: float, J: float, npts: int = 50) -> None:
    pwf, q = ipr_curve_fetkovich(p_res, J, npts=npts)
    plt.figure()
    plt.plot(q, pwf)
    plt.xlabel("q (STB/d)")
    plt.ylabel("pwf (psi)")
    plt.title("IPR - Fetkovich")
    plt.grid(True)
    plt.show()


def plot_ipr_jones(p_res: float, C: float, D: float, npts: int = 50) -> None:
    pwf, q = ipr_curve_jones(p_res, C, D, npts=npts)
    plt.figure()
    plt.plot(q, pwf)
    plt.xlabel("q (STB/d)")
    plt.ylabel("pwf (psi)")
    plt.title("IPR - Jones (C*q + D*q²)")
    plt.grid(True)
    plt.show()


def plot_ipr_standing(p_res: float, p_b: float, J: float, npts: int = 50) -> None:
    pwf, q = ipr_curve_standing(p_res, p_b, J, npts=npts)
    plt.figure()
    plt.plot(q, pwf)
    plt.xlabel("q (STB/d)")
    plt.ylabel("pwf (psi)")
    plt.title("IPR - Standing (compuesto sobre/baixo Pb)")
    plt.grid(True)
    plt.show()


__all__ = [
    "vogel_ipr",
    "fetkovich_ipr",
    "jones_ipr",
    "standing_ipr",
    "ipr_curve_vogel",
    "ipr_curve_fetkovich",
    "ipr_curve_jones",
    "ipr_curve_standing",
    "plot_ipr_vogel",
    "plot_ipr_fetkovich",
    "plot_ipr_jones",
    "plot_ipr_standing",
]
