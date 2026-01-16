import numpy as np
import matplotlib.pyplot as plt
import math


def flowline_pressure_drop(q, L, d, rho, mu, f=0.02, elev=0):
    """
    Calcula la pérdida de presión en una línea de flujo usando Darcy-Weisbach (monofásico).

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
        Viscosidad [cP] (no se usa en esta versión simplificada)
    f : float
        Factor de fricción (adimensional)
    elev : float
        Cambio de elevación [ft] (positivo = subida)

    Retorna
    -------
    delta_p : float
        Pérdida de presión total [psi]
    """
    # Conversión de caudal a ft3/s
    stb_to_ft3 = 5.615 / (24 * 3600)  # 1 STB/d = 5.615 ft3/d
    q_ft3_s = q * stb_to_ft3

    # Área de la tubería [ft2]
    d_ft = d / 12
    A = np.pi * (d_ft / 2) ** 2

    # Velocidad [ft/s]
    v = q_ft3_s / A

    # Constante
    g = 32.174  # ft/s2

    # Pérdida por fricción (psi)
    delta_p_fric = f * (L / d_ft) * (rho / 144) * (v ** 2 / (2 * g))

    # Pérdida/ganancia por elevación (psi)
    delta_p_elev = (rho / 144) * elev

    # Total
    delta_p = delta_p_fric + delta_p_elev
    return delta_p


def plot_flowline(q_range, L, d, rho, mu, f=0.02, elev=0):
    """
    Grafica la caída de presión en la línea de flujo en función del caudal (Darcy monofásico).
    """
    delta_p = [flowline_pressure_drop(q, L, d, rho, mu, f, elev) for q in q_range]

    plt.figure(figsize=(6, 5))
    plt.plot(q_range, delta_p, color="green", label="Flowline ΔP")
    plt.xlabel("Caudal [STB/d]")
    plt.ylabel("ΔP [psi]")
    plt.title("Pérdida de Presión en Flowline")
    plt.grid(True)
    plt.legend()
    plt.show()

    return delta_p


# ---------------------------
# Phase 2: Dispatcher + multifásico
# ---------------------------

def available_flowline_models():
    """
    Modelos disponibles para flowline ΔP:
      - "darcy" (legacy monofásico)
      - "beggs_brill" (multifásico, incluye inclinación y régimen)
    """
    return ["darcy", "beggs_brill"]


def _norm_model(name: str) -> str:
    return (name or "").strip().lower().replace("-", "_").replace(" ", "_")


def _theta_from_elev(L_ft: float, elev_ft: float) -> float:
    """
    Convierte elevación (ft) en ángulo desde horizontal (deg).
    elev>0 => uphill (theta>0). elev<0 => downhill (theta<0).
    """
    L = float(L_ft)
    if L <= 0:
        raise ValueError("L must be > 0")
    x = float(elev_ft) / L
    x = max(min(x, 1.0), -1.0)
    return math.degrees(math.asin(x))


def flowline_pressure_drop_beggs_brill(
    q_liq_stb_d: float,
    L_ft: float,
    d_in: float,
    rho_l_lbm_ft3: float,
    mu_l_cp: float,
    elev_ft: float = 0.0,
    q_gas_mscf_d: float = 0.0,
    rho_g_lbm_ft3: float = 2.0,
    mu_g_cp: float = 0.02,
    sigma_dyn_cm: float = 30.0,
    eps_in: float = 0.0006,
    theta_deg: float | None = None,
    p_psia: float | None = None,
) -> float:
    """
    Flowline ΔP (psi) usando Beggs & Brill core (dp/dz * L).

    Convenciones:
      - theta_deg: ángulo desde horizontal (+subida, -bajada).
        Si no se provee, se calcula desde elev_ft/L_ft.
      - ΔP puede ser negativo si hay bajada suficientemente grande (ganancia de presión).
    """
    from .multiphase import beggs_brill_pressure_gradient

    L = float(L_ft)
    if L <= 0:
        raise ValueError("L_ft must be > 0")
    d = float(d_in)
    if d <= 0:
        raise ValueError("d_in must be > 0")

    if theta_deg is None:
        theta = _theta_from_elev(L, elev_ft)
    else:
        theta = float(theta_deg)

    # conversiones a ft3/s
    STB_TO_FT3 = 5.615
    MSCF_TO_FT3 = 1000.0
    DAY_TO_S = 86400.0

    ql_ft3_s = (float(q_liq_stb_d) * STB_TO_FT3) / DAY_TO_S
    qg_ft3_s = (float(q_gas_mscf_d) * MSCF_TO_FT3) / DAY_TO_S

    d_ft = d / 12.0
    area = math.pi * (d_ft / 2.0) ** 2

    # caso no-flow: solo columna estática (líquido)
    if ql_ft3_s <= 0.0 and qg_ft3_s <= 0.0:
        # sin(theta) = elev/L cuando theta viene de elev, pero soporta theta explícito
        dp_dz_hydro = (float(rho_l_lbm_ft3) * math.sin(math.radians(theta))) / 144.0
        return float(dp_dz_hydro * L)

    v_sl = ql_ft3_s / area
    v_sg = qg_ft3_s / area

    res = beggs_brill_pressure_gradient(
        v_sl_ft_s=float(v_sl),
        v_sg_ft_s=float(v_sg),
        d_in=d,
        theta_deg=theta,
        rho_l_lbm_ft3=float(rho_l_lbm_ft3),
        rho_g_lbm_ft3=float(rho_g_lbm_ft3),
        mu_l_cp=float(mu_l_cp),
        mu_g_cp=float(mu_g_cp),
        sigma_dyn_cm=float(sigma_dyn_cm),
        eps_in=float(eps_in),
        p_psia=p_psia,
    )
    return float(res.dp_dz_psi_per_ft) * L


def flowline_pressure_drop_model(model, q, L, d, rho, mu, f=0.02, elev=0, **kwargs):
    """
    Dispatcher de Flowline ΔP por modelo.

    Mantiene la firma legacy para "darcy":
      q [STB/d], L [ft], d [in], rho [lb/ft3], mu [cP], elev [ft]

    Para "beggs_brill" acepta kwargs:
      q_gas_mscf_d, rho_g_lbm_ft3 (o rho_g), mu_g_cp (o mu_g),
      sigma_dyn_cm, eps_in, theta_deg (si no, se calcula desde elev/L),
      p_psia (opcional).
    """
    m = _norm_model(model)

    if m in ("darcy", "darcy_weisbach", "darcy-weisbach"):
        return flowline_pressure_drop(q=q, L=L, d=d, rho=rho, mu=mu, f=f, elev=elev)

    if m in ("beggs_brill", "beggs-brill", "bb"):
        q_gas = float(kwargs.get("q_gas_mscf_d", 0.0))
        rho_g = float(kwargs.get("rho_g_lbm_ft3", kwargs.get("rho_g", 2.0)))
        mu_g = float(kwargs.get("mu_g_cp", kwargs.get("mu_g", 0.02)))
        sigma = float(kwargs.get("sigma_dyn_cm", 30.0))
        eps_in = float(kwargs.get("eps_in", 0.0006))
        theta_deg = kwargs.get("theta_deg", None)
        p_psia = kwargs.get("p_psia", None)
        if p_psia is not None:
            p_psia = float(p_psia)

        return flowline_pressure_drop_beggs_brill(
            q_liq_stb_d=float(q),
            L_ft=float(L),
            d_in=float(d),
            rho_l_lbm_ft3=float(rho),
            mu_l_cp=float(mu),
            elev_ft=float(elev),
            q_gas_mscf_d=q_gas,
            rho_g_lbm_ft3=rho_g,
            mu_g_cp=mu_g,
            sigma_dyn_cm=sigma,
            eps_in=eps_in,
            theta_deg=None if theta_deg is None else float(theta_deg),
            p_psia=p_psia,
        )

    raise ValueError(f"Modelo flowline inválido: {model}. Usa: {available_flowline_models()}")


def plot_flowline_model(model, q_range, L, d, rho, mu, f=0.02, elev=0, **kwargs):
    """
    Grafica ΔP vs q para el modelo elegido.
    """
    q_range = list(q_range)
    dp = [flowline_pressure_drop_model(model, q, L, d, rho, mu, f=f, elev=elev, **kwargs) for q in q_range]

    plt.figure(figsize=(6, 5))
    plt.plot(q_range, dp, label=f"Flowline ΔP ({model})")
    plt.xlabel("Caudal [STB/d]")
    plt.ylabel("ΔP [psi]")
    plt.title("Pérdida de Presión en Flowline (Model)")
    plt.grid(True)
    plt.legend()
    plt.show()

    return dp


__all__ = [
    "flowline_pressure_drop",
    "plot_flowline",
    "available_flowline_models",
    "flowline_pressure_drop_model",
    "flowline_pressure_drop_beggs_brill",
    "plot_flowline_model",
]
