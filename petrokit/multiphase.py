# petrokit/multiphase.py
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, Optional, Tuple

FlowRegime = Literal["segregated", "transition", "intermittent", "distributed"]


@dataclass(frozen=True)
class BeggsBrillResult:
    dp_dz_psi_per_ft: float
    holdup: float
    regime: FlowRegime
    friction_factor_tp: float
    rho_m_lbm_ft3: float
    rho_ns_lbm_ft3: float
    re_ns: float


def _cp_to_lbm_ft_s(mu_cp: float) -> float:
    # 1 cP = 1e-3 Pa*s
    # 1 Pa*s = 0.671969 lbm/(ft*s)  -> 1 cP = 0.000671969 lbm/(ft*s)
    return mu_cp * 0.000671969


def _darcy_friction_factor_swamee_jain(re: float, rel_roughness: float) -> float:
    if re <= 0:
        raise ValueError("Reynolds must be positive.")
    if re < 2000.0:
        return 64.0 / re
    # Swamee-Jain explicit approximation for Darcy friction factor
    return 0.25 / (math.log10(rel_roughness / 3.7 + 5.74 / (re ** 0.9)) ** 2)


def beggs_brill_boundaries(c_l: float) -> Tuple[float, float, float, float]:
    """
    Beggs & Brill flow map boundaries.
    Returns (L1, L2, L3, L4) as described in common Beggs & Brill summaries.
    """
    if not (0.0 < c_l <= 1.0):
        raise ValueError("c_l must be in (0, 1].")

    L1 = 316.0 * (c_l ** 0.302)
    L2 = 0.0009252 * (c_l ** (-2.4684))
    L3 = 0.1 * (c_l ** (-1.4516))
    L4 = 0.5 * (c_l ** (-6.738))
    return L1, L2, L3, L4


def beggs_brill_flow_regime(v_sl_ft_s: float, v_sg_ft_s: float, d_ft: float) -> Tuple[FlowRegime, float, float]:
    """
    Identify Beggs & Brill flow regime based on mixture Froude number and no-slip liquid fraction.
    Returns (regime, c_l, fr).
    """
    if d_ft <= 0:
        raise ValueError("d_ft must be > 0.")
    if v_sl_ft_s < 0 or v_sg_ft_s < 0:
        raise ValueError("Superficial velocities must be >= 0.")

    v_m = v_sl_ft_s + v_sg_ft_s
    if v_m <= 0:
        raise ValueError("Mixture velocity must be > 0.")

    c_l = v_sl_ft_s / v_m
    # Guard to avoid 0**negative in boundaries; correlation is two-phase anyway.
    c_l_eff = max(min(c_l, 1.0), 1e-12)

    g = 32.174  # ft/s^2
    fr = (v_m ** 2) / (g * d_ft)

    L1, L2, L3, L4 = beggs_brill_boundaries(c_l_eff)

    # Classification per commonly used summary (includes "transition")
    if (c_l_eff < 0.01 and fr < L1) or (c_l_eff >= 0.01 and fr < L2):
        regime: FlowRegime = "segregated"
    elif L2 < fr < L3:
        regime = "transition"
    elif (0.01 <= c_l_eff < 0.4 and L3 < fr <= L1) or (c_l_eff >= 0.4 and L3 < fr <= L4):
        regime = "intermittent"
    else:
        # Distributed as remaining (high Fr region per map usage)
        regime = "distributed"

    return regime, c_l_eff, fr


def _holdup_horizontal(c_l: float, fr: float, regime: FlowRegime) -> float:
    # E_L(0) = a * C_L^b / Fr^c
    if fr <= 0:
        raise ValueError("Froude must be > 0.")

    if regime == "segregated":
        a, b, c = 0.98, 0.4846, 0.0868
    elif regime == "intermittent":
        a, b, c = 0.845, 0.5351, 0.0173
    else:
        # distributed (and a fallback for transition handled elsewhere)
        a, b, c = 1.065, 0.5824, 0.0609

    el0 = a * (c_l ** b) / (fr ** c)
    # Constraint: E_L(0) >= C_L
    return max(el0, c_l)


def _beta_inclination(
    c_l: float,
    fr: float,
    n_lv: float,
    theta_deg: float,
    regime: FlowRegime,
) -> float:
    """
    beta = (1 - C_L)*ln(d*C_L^e*N_LV^f*Fr^g)
    Coefficients depend on uphill/downhill and regime.
    """
    # If horizontal, correction term contributes ~0 anyway; keep stable.
    if abs(theta_deg) < 1e-12:
        return 0.0

    uphill = theta_deg > 0.0

    # Coeffs as commonly tabulated in summaries
    if uphill:
        if regime == "segregated":
            d, e, f, g = 0.011, -3.768, 3.539, -1.614
        elif regime == "intermittent":
            d, e, f, g = 2.96, 0.305, -0.4473, 0.0978
        else:
            # Distributed: beta=0
            return 0.0
    else:
        # Downhill: same coeffs for all regimes in the summary table
        d, e, f, g = 4.7, -0.3692, 0.1244, -0.5056

    arg = d * (c_l ** e) * (n_lv ** f) * (fr ** g)
    # Avoid log of non-positive (shouldn't happen for physical inputs)
    arg = max(arg, 1e-30)
    return (1.0 - c_l) * math.log(arg)


def beggs_brill_liquid_holdup(
    v_sl_ft_s: float,
    v_sg_ft_s: float,
    d_ft: float,
    theta_deg: float,
    rho_l_lbm_ft3: float,
    sigma_dyn_cm: float = 30.0,
) -> Tuple[float, FlowRegime]:
    """
    Returns (E_L(theta), regime).

    Units:
      - v_sl_ft_s, v_sg_ft_s: ft/s
      - d_ft: ft
      - theta_deg: degrees from horizontal (uphill +, downhill -)
      - rho_l_lbm_ft3: lbm/ft^3
      - sigma_dyn_cm: dyn/cm (used only for N_LV)
    """
    if rho_l_lbm_ft3 <= 0:
        raise ValueError("rho_l_lbm_ft3 must be > 0.")
    if sigma_dyn_cm <= 0:
        raise ValueError("sigma_dyn_cm must be > 0.")

    regime, c_l, fr = beggs_brill_flow_regime(v_sl_ft_s, v_sg_ft_s, d_ft)

    v_m = v_sl_ft_s + v_sg_ft_s
    el0 = _holdup_horizontal(c_l, fr, regime)

    # Liquid velocity number (dimensionless) as per common summary
    g = 32.174
    n_lv = 1.938 * v_sl_ft_s * ((rho_l_lbm_ft3 / (g * sigma_dyn_cm)) ** 0.25)

    beta = _beta_inclination(c_l=c_l, fr=fr, n_lv=n_lv, theta_deg=theta_deg, regime=regime)

    # B(theta) = 1 + beta*(sin(1.8θ) - (1/3)sin^3(1.8θ))
    # θ is in degrees; 1.8θ treated in degrees in most summaries -> convert to radians after multiplying.
    s = math.sin(math.radians(1.8 * theta_deg))
    b_theta = 1.0 + beta * (s - (1.0 / 3.0) * (s ** 3))

    el_theta_main = el0 * b_theta
    el_theta_main = min(max(el_theta_main, c_l), 1.0)

    if regime != "transition":
        return el_theta_main, regime

    # Transition: weighted between segregated and intermittent solutions
    L1, L2, L3, L4 = beggs_brill_boundaries(c_l)
    A = (L3 - fr) / (L3 - L2)
    A = min(max(A, 0.0), 1.0)
    B = 1.0 - A

    # Compute segregated and intermittent holdups (at same inputs) and blend
    el0_seg = _holdup_horizontal(c_l, fr, "segregated")
    beta_seg = _beta_inclination(c_l=c_l, fr=fr, n_lv=n_lv, theta_deg=theta_deg, regime="segregated")
    b_seg = 1.0 + beta_seg * (s - (1.0 / 3.0) * (s ** 3))
    el_seg = min(max(el0_seg * b_seg, c_l), 1.0)

    el0_int = _holdup_horizontal(c_l, fr, "intermittent")
    beta_int = _beta_inclination(c_l=c_l, fr=fr, n_lv=n_lv, theta_deg=theta_deg, regime="intermittent")
    b_int = 1.0 + beta_int * (s - (1.0 / 3.0) * (s ** 3))
    el_int = min(max(el0_int * b_int, c_l), 1.0)

    el_tr = A * el_seg + B * el_int
    el_tr = min(max(el_tr, c_l), 1.0)
    return el_tr, "transition"


def beggs_brill_pressure_gradient(
    v_sl_ft_s: float,
    v_sg_ft_s: float,
    d_in: float,
    theta_deg: float,
    rho_l_lbm_ft3: float,
    rho_g_lbm_ft3: float,
    mu_l_cp: float,
    mu_g_cp: float = 0.02,
    sigma_dyn_cm: float = 30.0,
    eps_in: float = 0.0006,
    p_psia: Optional[float] = None,
) -> BeggsBrillResult:
    """
    Beggs & Brill pressure gradient (psi/ft).

    Inputs:
      - v_sl_ft_s, v_sg_ft_s: superficial velocities [ft/s]
      - d_in: pipe ID [in]
      - theta_deg: angle from horizontal (uphill +, downhill -) [deg]
      - rho_l_lbm_ft3, rho_g_lbm_ft3: densities [lbm/ft^3]
      - mu_l_cp, mu_g_cp: viscosities [cP]
      - sigma_dyn_cm: surface tension [dyn/cm] (for N_LV in holdup)
      - eps_in: roughness [in]
      - p_psia: if provided, applies acceleration correction via Ek (optional)

    Returns:
      BeggsBrillResult with dp/dz in psi/ft and intermediate values.
    """
    if d_in <= 0:
        raise ValueError("d_in must be > 0.")
    if rho_l_lbm_ft3 <= 0 or rho_g_lbm_ft3 <= 0:
        raise ValueError("Densities must be > 0.")
    if mu_l_cp <= 0 or mu_g_cp <= 0:
        raise ValueError("Viscosities must be > 0.")
    if v_sl_ft_s < 0 or v_sg_ft_s < 0:
        raise ValueError("Superficial velocities must be >= 0.")
    if p_psia is not None and p_psia <= 0:
        raise ValueError("p_psia must be > 0 if provided.")

    d_ft = d_in / 12.0
    v_m = v_sl_ft_s + v_sg_ft_s
    if v_m <= 0:
        raise ValueError("Mixture velocity must be > 0.")

    # Flow regime + holdup
    holdup, regime = beggs_brill_liquid_holdup(
        v_sl_ft_s=v_sl_ft_s,
        v_sg_ft_s=v_sg_ft_s,
        d_ft=d_ft,
        theta_deg=theta_deg,
        rho_l_lbm_ft3=rho_l_lbm_ft3,
        sigma_dyn_cm=sigma_dyn_cm,
    )

    # No-slip liquid fraction for no-slip properties
    c_l = v_sl_ft_s / v_m

    rho_ns = rho_l_lbm_ft3 * c_l + rho_g_lbm_ft3 * (1.0 - c_l)
    rho_m = rho_l_lbm_ft3 * holdup + rho_g_lbm_ft3 * (1.0 - holdup)

    mu_l = _cp_to_lbm_ft_s(mu_l_cp)
    mu_g = _cp_to_lbm_ft_s(mu_g_cp)
    mu_ns = mu_l * c_l + mu_g * (1.0 - c_l)

    re_ns = rho_ns * v_m * d_ft / mu_ns

    rel_roughness = (eps_in / 12.0) / d_ft
    f_ns = _darcy_friction_factor_swamee_jain(re_ns, rel_roughness)

    # Two-phase friction factor multiplier
    y = c_l / (holdup ** 2)
    y = max(y, 1e-12)
    ln_y = math.log(y)

    if 1.0 < y < 1.2:
        S = math.log(2.2 * y - 1.2)
    else:
        denom = (-0.0523 + 3.182 * ln_y - 0.8725 * (ln_y ** 2) + 0.01853 * (ln_y ** 4))
        denom = denom if abs(denom) > 1e-12 else (1e-12 if denom >= 0 else -1e-12)
        S = ln_y / denom

    f_tp = f_ns * math.exp(S)

    # Pressure gradients (psi/ft)
    # Hydrostatic component (English units simplification g/gc cancels)
    dp_dz_ele = (rho_m * math.sin(math.radians(theta_deg))) / 144.0

    # Friction component
    g_c = 32.174
    dp_dz_fric = (2.0 * f_tp * (v_m ** 2) * rho_ns) / (144.0 * g_c * d_ft)

    # Acceleration correction factor Ek (optional)
    if p_psia is None:
        ek = 0.0
    else:
        ek = (rho_m * v_m * v_sg_ft_s) / (g_c * p_psia)

    ek = min(max(ek, 0.0), 0.95)  # keep stable

    dp_dz_total = (dp_dz_fric + dp_dz_ele) / (1.0 - ek)

    return BeggsBrillResult(
        dp_dz_psi_per_ft=dp_dz_total,
        holdup=holdup,
        regime=regime,
        friction_factor_tp=f_tp,
        rho_m_lbm_ft3=rho_m,
        rho_ns_lbm_ft3=rho_ns,
        re_ns=re_ns,
    )
