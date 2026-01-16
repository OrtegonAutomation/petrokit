"""
petrokit.pvt

Utilidades PVT *simplificadas* (black-oil) para generar tablas rápidas.

Estas correlaciones NO pretenden reemplazar un PVT de laboratorio ni un simulador:
son aproximaciones típicas para ejercicios de ingeniería, screening y notebooks.

Unidades (convención):
- Presión: psia
- Temperatura: °F
- API: grados API
- Gravedad específica gas: gamma_g (aire = 1)
- Rs: scf/STB
- Bo: RB/STB
- Bg: RB/scf

Correlaciones incluidas:
- Standing (1947): Rs(P), Pb(Rsb), Bo(Rs) (aceite saturado)
- Sutton (1985): pseudo-críticas Tpc/Ppc a partir de gamma_g (gas dulce)
- Papay: z-factor explícito (con Ppr, Tpr)
- Bg: factor volumétrico de gas (RB/scf)

Notas:
- Rango de validez limitado (depende de cada correlación).
- Para P > Pb (aceite sub-saturado), esta fase usa aproximación simple:
  Rs = Rsb y Bo = Bob (constante). En fases posteriores se puede agregar co(P).
"""

from __future__ import annotations

from typing import Dict, Union, Optional
import numpy as np

ArrayLike = Union[float, int, np.ndarray]


def _to_array(x: ArrayLike) -> np.ndarray:
    """Convierte escalar/lista/ndarray a ndarray(float)."""
    return np.asarray(x, dtype=float)


def _to_rankine(t_f: float) -> float:
    """°F -> °R."""
    return float(t_f) + 459.67


def oil_gamma_o(api: float) -> float:
    """Gravedad específica del aceite (agua=1) a partir de API."""
    api = float(api)
    if api <= 0:
        raise ValueError("API debe ser > 0")
    return 141.5 / (api + 131.5)


def oil_rs_standing(p_psia: ArrayLike, t_f: float, api: float, gamma_g: float) -> np.ndarray:
    """
    Standing (1947) - Rs (gas en solución) para aceite saturado.

    Rs = γg * [ (18.2 P + 1.4) * 10^X ]^1.2048
    X = 0.0125*API - 0.00091*(T - 460)

    Inputs:
      p_psia: presión (psia)
      t_f: temperatura (°F)
      api: °API
      gamma_g: gravedad específica del gas (aire=1)

    Returns:
      Rs [scf/STB] (ndarray)
    """
    p = _to_array(p_psia)
    t = float(t_f)
    api = float(api)
    gg = float(gamma_g)

    if gg <= 0:
        raise ValueError("gamma_g debe ser > 0")

    x = 0.0125 * api - 0.00091 * (t - 460.0)
    rs = gg * np.power(((p / 18.2 + 1.4) * np.power(10.0, x)), 1.2048)
    return np.maximum(rs, 0.0)


def oil_pb_standing(rsb_scf_stb: float, t_f: float, api: float, gamma_g: float) -> float:
    """
    Standing (1947) - presión de burbuja Pb a partir de Rsb.

    Pb = 18.2 * [ (Rsb/γg)^0.83 * 10^a - 1.4 ]
    a  = 0.00091*(T - 460) - 0.0125*API

    Returns:
      Pb [psia]
    """
    rsb = float(rsb_scf_stb)
    t = float(t_f)
    api = float(api)
    gg = float(gamma_g)

    if rsb < 0:
        raise ValueError("rsb_scf_stb debe ser >= 0")
    if gg <= 0:
        raise ValueError("gamma_g debe ser > 0")

    a = 0.00091 * (t - 460.0) - 0.0125 * api
    pb = 18.2 * (np.power(rsb / gg, 0.83) * np.power(10.0, a) - 1.4)
    return float(max(pb, 0.0))


def oil_bo_standing(rs_scf_stb: ArrayLike, t_f: float, api: float, gamma_g: float) -> np.ndarray:
    """
    Standing (1947) - Bo (oil FVF) para aceite saturado.

    Bo = 0.9759 + 0.000120 * F^1.2
    F  = Rs*(γo/γg) + 1.25*T

    Returns:
      Bo [RB/STB] (ndarray)
    """
    rs = _to_array(rs_scf_stb)
    t = float(t_f)
    api = float(api)
    gg = float(gamma_g)

    if gg <= 0:
        raise ValueError("gamma_g debe ser > 0")

    go = oil_gamma_o(api)
    f = rs * (go / gg) + 1.25 * t
    bo = 0.9759 + 0.000120 * np.power(f, 1.2)
    return np.maximum(bo, 1e-9)


def gas_pseudocritical_sutton(gamma_g: float) -> Dict[str, float]:
    """
    Sutton (1985) - pseudo-críticas para gas dulce (sweet gas).

    Ppc = 787 - 147*γg - 7.5*γg^2  [psia]
    Tpc = 169 + 314*γg - 71.3*γg^2 [°R]
    """
    gg = float(gamma_g)
    if gg <= 0:
        raise ValueError("gamma_g debe ser > 0")

    ppc = 787.0 - 147.0 * gg - 7.5 * gg**2
    tpc = 169.0 + 314.0 * gg - 71.3 * gg**2
    return {"Ppc_psia": float(ppc), "Tpc_R": float(tpc)}


def gas_z_factor_papay(p_psia: ArrayLike, t_f: float, gamma_g: float) -> np.ndarray:
    """
    Papay - z-factor explícito (aprox) en función de Ppr y Tpr.

    Z = 1 - (3.53*Ppr)/(10^(0.9813*Tpr)) + (0.274*Ppr^2)/(10^(0.815*Tpr))

    Inputs:
      p_psia: presión (psia)
      t_f: temperatura (°F)
      gamma_g: gravedad específica gas (aire=1)

    Returns:
      Z (ndarray)
    """
    p = _to_array(p_psia)
    t_r = _to_rankine(t_f)

    pc = gas_pseudocritical_sutton(gamma_g)
    ppc = pc["Ppc_psia"]
    tpc = pc["Tpc_R"]

    if ppc <= 0 or tpc <= 0:
        raise ValueError("Pseudo-críticas inválidas (revisa gamma_g)")

    ppr = p / ppc
    tpr = t_r / tpc

    z = 1.0 - (3.53 * ppr) / np.power(10.0, 0.9813 * tpr) + (0.274 * ppr**2) / np.power(10.0, 0.815 * tpr)
    # Evitar valores no físicos por fuera de rango:
    return np.clip(z, 0.2, 2.0)


def gas_bg_rb_per_scf(p_psia: ArrayLike, t_f: float, z: ArrayLike) -> np.ndarray:
    """
    Factor volumétrico de gas Bg (RB/scf).

    Bg = 0.00504 * Z * T(°R) / P(psia)
    """
    p = _to_array(p_psia)
    z = _to_array(z)
    t_r = _to_rankine(t_f)

    if np.any(p <= 0):
        raise ValueError("La presión debe ser > 0")

    bg = 0.00504 * z * t_r / p
    return np.maximum(bg, 1e-12)


def build_pvt_table(
    p_psia: ArrayLike,
    t_f: float,
    api: float,
    gamma_g: float,
    pb_psia: Optional[float] = None,
    rsb_scf_stb: Optional[float] = None,
    z_method: str = "papay",
) -> Dict[str, np.ndarray]:
    """
    Construye una tabla PVT mínima (arrays) para usar en notebooks o cálculos rápidos.

    Parámetros:
      p_psia: vector de presiones (psia)
      t_f: temperatura (°F)
      api: °API
      gamma_g: gas gravity (aire=1)
      pb_psia: (opcional) bubble point (psia)
      rsb_scf_stb: (opcional) Rs a Pb (scf/STB)
      z_method: 'papay' o 'ideal'

    Reglas:
      - Si pb_psia se da, calcula rsb con Standing a Pb.
      - Si rsb_scf_stb se da, calcula pb con Standing.
      - Si no se da ninguno, asume pb = max(P) (solo para generar una tabla).
      - Para P > Pb: Rs = Rsb y Bo = Bob (aprox simple).

    Returns dict con:
      P_psia, Rs_scf_stb, Bo_rb_stb, Z, Bg_rb_scf, pb_psia, rsb_scf_stb
    """
    p = _to_array(p_psia)

    if pb_psia is None and rsb_scf_stb is None:
        pb = float(np.nanmax(p))
        rsb = float(oil_rs_standing(pb, t_f, api, gamma_g))
    elif pb_psia is not None and rsb_scf_stb is None:
        pb = float(pb_psia)
        rsb = float(oil_rs_standing(pb, t_f, api, gamma_g))
    elif pb_psia is None and rsb_scf_stb is not None:
        rsb = float(rsb_scf_stb)
        pb = float(oil_pb_standing(rsb, t_f, api, gamma_g))
    else:
        pb = float(pb_psia)
        rsb = float(rsb_scf_stb)

    # Rs: saturado (clamp a Rsb)
    rs_sat = oil_rs_standing(p, t_f, api, gamma_g)
    rs = np.minimum(rs_sat, rsb)

    # Bo: saturado a partir de Rs; por encima de Pb => Bob (constante)
    bo_sat = oil_bo_standing(rs, t_f, api, gamma_g)
    bob = float(oil_bo_standing(rsb, t_f, api, gamma_g))
    bo = np.where(p > pb, bob, bo_sat)

    # Z
    z_method = str(z_method).lower().strip()
    if z_method == "ideal":
        z = np.ones_like(p, dtype=float)
    elif z_method == "papay":
        z = gas_z_factor_papay(p, t_f, gamma_g)
    else:
        raise ValueError("z_method debe ser 'papay' o 'ideal'")

    bg = gas_bg_rb_per_scf(p, t_f, z)

    return {
        "P_psia": p,
        "Rs_scf_stb": rs,
        "Bo_rb_stb": bo,
        "Z": z,
        "Bg_rb_scf": bg,
        "pb_psia": float(pb),
        "rsb_scf_stb": float(rsb),
    }


__all__ = [
    "oil_gamma_o",
    "oil_rs_standing",
    "oil_pb_standing",
    "oil_bo_standing",
    "gas_pseudocritical_sutton",
    "gas_z_factor_papay",
    "gas_bg_rb_per_scf",
    "build_pvt_table",
]

