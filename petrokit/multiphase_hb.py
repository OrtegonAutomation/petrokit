from __future__ import annotations

import math
from dataclasses import dataclass


G_FT_S2 = 32.174
LBFT3_TO_PSI_PER_FT = 1.0 / 144.0  # rho(lb/ft3) * ft /144 = psi


def _cp_to_lbm_ft_s(mu_cp: float) -> float:
    # 1 cP = 0.000671969 lbm/(ft*s) (aprox)
    return float(mu_cp) * 0.000671969


def _darcy_friction_factor_swamee_jain(re: float, rel_roughness: float) -> float:
    """
    Darcy friction factor using Swamee-Jain.
    - Laminar: f = 64/Re
    """
    re = float(re)
    if re <= 0:
        raise ValueError("Reynolds must be > 0")

    if re < 2000.0:
        return 64.0 / re

    rr = max(float(rel_roughness), 0.0)
    return 0.25 / (math.log10(rr / 3.7 + 5.74 / (re ** 0.9)) ** 2)


@dataclass(frozen=True)
class HagedornBrownResult:
    dp_dz_psi_per_ft: float
    holdup: float
    rho_m_lbm_ft3: float
    re_m: float
    f_darcy: float


def hagedorn_brown_holdup(
    v_sl_ft_s: float,
    v_sg_ft_s: float,
    k: float = 0.15,
) -> float:
    """
    Holdup empírico (versión "core" alineada con tu prototipo actual).
    Basado en fracción no-slip y una corrección suave por gas.
    """
    v_sl = float(v_sl_ft_s)
    v_sg = float(v_sg_ft_s)
    if v_sl < 0 or v_sg < 0:
        raise ValueError("Superficial velocities must be >= 0")

    v_m = v_sl + v_sg
    if v_m <= 0:
        return 1.0  # no-flow -> todo líquido

    lam_l = v_sl / v_m  # no-slip liquid fraction

    # En tu prototipo: H = lam_l * (1 + k*sqrt(v_sg/v_m))
    h = lam_l * (1.0 + float(k) * math.sqrt(max(v_sg / (v_m + 1e-12), 0.0)))
    # clamps físicos
    h = max(h, lam_l)
    return min(max(h, 1e-6), 0.999999)


def hagedorn_brown_pressure_gradient(
    v_sl_ft_s: float,
    v_sg_ft_s: float,
    d_in: float,
    theta_deg: float,
    rho_l_lbm_ft3: float,
    rho_g_lbm_ft3: float,
    mu_l_cp: float,
    mu_g_cp: float,
    eps_in: float = 0.0006,
    k_holdup: float = 0.15,
) -> HagedornBrownResult:
    """
    Hagedorn & Brown (core simplificado):
    - Devuelve gradiente total dp/dz [psi/ft] = fricción + elevación
    - Ignora aceleración (se puede agregar después)

    Inputs:
      v_sl_ft_s, v_sg_ft_s : ft/s
      d_in                : in
      theta_deg           : deg (desde horizontal; vertical = 90)
      rho_*               : lbm/ft3
      mu_*                : cP
      eps_in              : in (rugosidad)
    """
    if d_in <= 0:
        raise ValueError("d_in must be > 0")
    if rho_l_lbm_ft3 <= 0 or rho_g_lbm_ft3 <= 0:
        raise ValueError("densities must be > 0")
    if mu_l_cp <= 0 or mu_g_cp <= 0:
        raise ValueError("viscosities must be > 0")
    if theta_deg < -90 or theta_deg > 90:
        raise ValueError("theta_deg must be within [-90, 90]")

    v_sl = float(v_sl_ft_s)
    v_sg = float(v_sg_ft_s)
    if v_sl < 0 or v_sg < 0:
        raise ValueError("superficial velocities must be >= 0")

    d_ft = float(d_in) / 12.0
    v_m = v_sl + v_sg

    # no-flow: solo hidrostática líquida
    if v_m <= 0:
        dpdz = (float(rho_l_lbm_ft3) * math.sin(math.radians(theta_deg))) * LBFT3_TO_PSI_PER_FT
        return HagedornBrownResult(
            dp_dz_psi_per_ft=float(dpdz),
            holdup=1.0,
            rho_m_lbm_ft3=float(rho_l_lbm_ft3),
            re_m=0.0,
            f_darcy=0.0,
        )

    # Holdup
    h = hagedorn_brown_holdup(v_sl, v_sg, k=float(k_holdup))

    rho_m = h * float(rho_l_lbm_ft3) + (1.0 - h) * float(rho_g_lbm_ft3)

    # Viscosidad mezcla (simple)
    mu_mix_cp = h * float(mu_l_cp) + (1.0 - h) * float(mu_g_cp)
    mu_mix = _cp_to_lbm_ft_s(mu_mix_cp)

    # Reynolds mezcla
    re_m = (rho_m * v_m * d_ft) / max(mu_mix, 1e-12)

    # Fricción (Darcy)
    rel_rough = (float(eps_in) / 12.0) / d_ft
    f = _darcy_friction_factor_swamee_jain(re_m, rel_rough)

    # dp/dz fricción: f*(rho*v^2)/(2*g*D)  -> psi/ft
    dpdz_fric = f * (rho_m * (v_m ** 2)) / (2.0 * G_FT_S2 * d_ft) * LBFT3_TO_PSI_PER_FT

    # dp/dz elevación
    dpdz_hydro = rho_m * math.sin(math.radians(theta_deg)) * LBFT3_TO_PSI_PER_FT

    dpdz_total = float(dpdz_fric + dpdz_hydro)

    return HagedornBrownResult(
        dp_dz_psi_per_ft=dpdz_total,
        holdup=h,
        rho_m_lbm_ft3=rho_m,
        re_m=float(re_m),
        f_darcy=float(f),
    )
