import numpy as np
import matplotlib.pyplot as plt


# ---------------------------
# VLP monofásico (legacy)
# ---------------------------

def vlp_curve(q_range, well_depth, rho, mu, d, f=0.02):
    """
    Calcula la curva VLP usando una fórmula simplificada Darcy-Weisbach.

    q_range : array de caudales [STB/d]
    well_depth : profundidad medida del pozo [ft]
    rho : densidad del fluido [lb/ft3]
    mu : viscosidad dinámica [cp] (no se usa aquí pero se deja para ampliar)
    d : diámetro interior de la tubería [in]
    f : factor de fricción (adimensional)

    Retorna:
    pwf_range : presión en fondo fluyente [psi]
    """
    q_range = np.asarray(q_range, dtype=float)

    # Conversión de caudal a ft3/s
    stb_to_ft3 = 5.615 / (24 * 3600)  # 1 STB/d = 5.615 ft3/d
    q_ft3_s = q_range * stb_to_ft3

    # Área de la tubería [ft2]
    d_ft = d / 12
    A = np.pi * (d_ft / 2) ** 2

    # Velocidad [ft/s]
    v = q_ft3_s / A

    # Head loss por fricción (nota: fórmula legacy simplificada)
    g = 32.174  # ft/s2
    delta_p_fric = f * (well_depth / d_ft) * (rho / 144) * (v ** 2 / 2 / g)

    # Head estático (columna fluido)
    delta_p_hydro = (rho / 144) * well_depth

    # Presión en fondo fluyente [psi]
    pwf_range = delta_p_fric + delta_p_hydro

    return pwf_range


def plot_vlp(q_range, well_depth, rho, mu, d, f=0.02):
    """
    Genera y grafica la curva VLP (Darcy simplificado).
    """
    pwf_range = vlp_curve(q_range, well_depth, rho, mu, d, f)
    plt.figure(figsize=(6, 5))
    plt.plot(q_range, pwf_range, color="red", label="VLP")
    plt.xlabel("Caudal [STB/d]")
    plt.ylabel("Presión Fondo Fluyente pwf [psi]")
    plt.title("Curva VLP")
    plt.grid(True)
    plt.legend()
    plt.show()


# ---------------------------
# API extendida (Phase 2): dispatcher de modelos
# ---------------------------

def available_vlp_models():
    """
    Modelos disponibles para VLP (según README):
      - "darcy"
      - "beggs_brill"
      - "hagedorn_brown"
    """
    return ["darcy", "beggs_brill", "hagedorn_brown"]


def _norm_model(name: str) -> str:
    return (name or "").strip().lower().replace("-", "_").replace(" ", "_")


def vlp_curve_beggs_brill(q_range, well_depth, rho, mu, d, **kwargs):
    """
    VLP usando Beggs & Brill (core en petrokit.multiphase).

    Mantiene la firma base:
      q_range [STB/d], well_depth [ft], rho [lbm/ft3], mu [cP], d [in]

    kwargs opcionales:
      theta_deg: ángulo desde horizontal (deg). Default 90 (vertical)
      p_wh: presión en cabeza (psi). Default 0
      q_gas_mscf_d: gas constante (MSCF/d). Default 0
      q_gas_range_mscf_d: array gas (MSCF/d) mismo largo que q_range
      rho_g: densidad gas (lbm/ft3). Default 2.0
      mu_g: viscosidad gas (cP). Default 0.02
      sigma_dyn_cm: tensión superficial (dyn/cm). Default 30
      eps_in: rugosidad (in). Default 0.0006
      p_psia: presión para corrección aceleración (psia). Default None
    """
    from .multiphase import beggs_brill_pressure_gradient

    q = np.asarray(q_range, dtype=float)
    L = float(well_depth)
    rho_l = float(rho)
    mu_l = float(mu)
    d_in = float(d)

    theta_deg = float(kwargs.get("theta_deg", 90.0))
    p_wh = float(kwargs.get("p_wh", 0.0))

    qg_const = kwargs.get("q_gas_mscf_d", None)
    qg_range = kwargs.get("q_gas_range_mscf_d", None)
    if (qg_const is not None) and (qg_range is not None):
        raise ValueError("Usa solo uno: q_gas_mscf_d o q_gas_range_mscf_d")

    if qg_range is not None:
        qg_mscf_d = np.asarray(qg_range, dtype=float)
        if qg_mscf_d.shape != q.shape:
            raise ValueError("q_gas_range_mscf_d debe tener el mismo tamaño que q_range")
    else:
        qg_mscf_d = np.full_like(q, float(qg_const) if qg_const is not None else 0.0)

    rho_g = float(kwargs.get("rho_g", 2.0))
    mu_g = float(kwargs.get("mu_g", 0.02))
    sigma = float(kwargs.get("sigma_dyn_cm", 30.0))
    eps_in = float(kwargs.get("eps_in", 0.0006))
    p_psia = kwargs.get("p_psia", None)
    if p_psia is not None:
        p_psia = float(p_psia)

    # conversions
    STB_TO_FT3 = 5.615
    MSCF_TO_FT3 = 1000.0
    DAY_TO_S = 86400.0

    d_ft = d_in / 12.0
    area = np.pi * (d_ft / 2.0) ** 2

    pwf = np.zeros_like(q, dtype=float)

    sin_theta = np.sin(np.deg2rad(theta_deg))

    for i, (ql_stb_d, qg_m) in enumerate(zip(q, qg_mscf_d)):
        ql_ft3_s = (ql_stb_d * STB_TO_FT3) / DAY_TO_S
        qg_ft3_s = (qg_m * MSCF_TO_FT3) / DAY_TO_S

        # Caso q=0: solo columna estática (evita v_m=0 en el core)
        if ql_ft3_s <= 0.0 and qg_ft3_s <= 0.0:
            dp_dz = (rho_l * sin_theta) / 144.0  # psi/ft
            pwf[i] = p_wh + dp_dz * L
            continue

        v_sl = ql_ft3_s / area
        v_sg = qg_ft3_s / area

        res = beggs_brill_pressure_gradient(
            v_sl_ft_s=float(v_sl),
            v_sg_ft_s=float(v_sg),
            d_in=d_in,
            theta_deg=theta_deg,
            rho_l_lbm_ft3=rho_l,
            rho_g_lbm_ft3=rho_g,
            mu_l_cp=mu_l,
            mu_g_cp=mu_g,
            sigma_dyn_cm=sigma,
            eps_in=eps_in,
            p_psia=p_psia,
        )
        pwf[i] = p_wh + float(res.dp_dz_psi_per_ft) * L

    return pwf


def vlp_curve_hagedorn_brown(q_range, well_depth, rho, mu, d, **kwargs):
    """
    VLP usando Hagedorn & Brown (prototipo en petrokit.vlp_advanced).

    kwargs opcionales:
      theta_deg: default 90
      p_wh: default 0
      q_gas_mscf_d: default 0
      rho_g, mu_g: defaults razonables
      roughness_ft: default 0.0005  (ojo: este prototipo asume roughness en ft)
      f_guess: opcional
    """
    from .vlp_advanced import hagedorn_brown_dp

    q = np.asarray(q_range, dtype=float)
    L = float(well_depth)
    rho_l = float(rho)
    mu_l = float(mu)
    d_in = float(d)

    theta_deg = float(kwargs.get("theta_deg", 90.0))
    p_wh = float(kwargs.get("p_wh", 0.0))

    qg = float(kwargs.get("q_gas_mscf_d", 0.0))
    rho_g = float(kwargs.get("rho_g", 2.0))
    mu_g = float(kwargs.get("mu_g", 0.02))
    roughness_ft = float(kwargs.get("roughness_ft", 0.0005))
    f_guess = kwargs.get("f_guess", None)

    pwf = np.zeros_like(q, dtype=float)
    for i, ql in enumerate(q):
        dp_total_psi, _holdup = hagedorn_brown_dp(
            q_liq_stb_d=float(ql),
            q_gas_mscf_d=qg,
            d_in=d_in,
            L_ft=L,
            rho_l=rho_l,
            rho_g=rho_g,
            mu_l=mu_l,
            mu_g=mu_g,
            theta_deg=theta_deg,
            roughness=roughness_ft,
            f_guess=f_guess,
        )
        pwf[i] = p_wh + float(dp_total_psi)

    return pwf


def vlp_curve_model(model, q_range, well_depth, rho, mu, d, **kwargs):
    """
    Dispatcher de VLP por modelo.
    """
    m = _norm_model(model)

    if m in ("darcy", "darcy_weisbach", "darcy-weisbach"):
        f = float(kwargs.get("f", 0.02))
        return vlp_curve(q_range, well_depth=well_depth, rho=rho, mu=mu, d=d, f=f)

    if m in ("beggs_brill", "beggs-brill", "bb"):
        return vlp_curve_beggs_brill(q_range, well_depth=well_depth, rho=rho, mu=mu, d=d, **kwargs)

    if m in ("hagedorn_brown", "hagedorn-brown", "hb"):
        return vlp_curve_hagedorn_brown(q_range, well_depth=well_depth, rho=rho, mu=mu, d=d, **kwargs)

    raise ValueError(f"Modelo VLP inválido: {model}. Usa: {available_vlp_models()}")


__all__ = [
    "vlp_curve",
    "plot_vlp",
    "available_vlp_models",
    "vlp_curve_model",
    "vlp_curve_beggs_brill",
    "vlp_curve_hagedorn_brown",
]
