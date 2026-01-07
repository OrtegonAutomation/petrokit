"""
petrokit.vlp_advanced
Versión prototipo/robusta de modelos VLP multifásicos:
- beggs_brill_dp(...) : implementación inspirada en Beggs & Brill (simplificada)
- hagedorn_brown_dp(...): implementación simplificada de Hagedorn & Brown

Salida:
- presión diferencial total estimada (psi) sobre la longitud L_ft
- holdup estimado (fracción volumétrica líquida)

NOTA:
- Implementación educativa/prototipo. Validar con datos de campo antes de tomar decisiones.
"""

from typing import Tuple, Iterable
import numpy as np

# Conversiones y constantes
STB_TO_FT3 = 5.615      # 1 STB ≈ 5.615 ft3
MSCF_TO_FT3 = 1000.0    # 1 MSCF = 1000 ft3
G_FT_S2 = 32.174        # gravedad [ft/s^2]
LBFT3_TO_PSI = 1 / 144.0  # lb/ft3 -> psi/ft factor (rho/144 -> psi/ft when multiplied by ft)


# ---------------------------
# Helpers de conversión y fricción
# ---------------------------

def _stb_day_to_ft3_s(q_stb_d: float) -> float:
    return float(q_stb_d) * STB_TO_FT3 / 86400.0

def _mscf_day_to_ft3_s(q_mscf_d: float) -> float:
    return float(q_mscf_d) * MSCF_TO_FT3 / 86400.0

def _in_to_ft(x_in: float) -> float:
    return float(x_in) / 12.0

def swamee_jain_f(re, rel_roughness):
    """
    Swamee-Jain explicit formula for friction factor f (valid for turbulent flow, Re>3000 typically).
    re: Reynolds number
    rel_roughness: epsilon/D
    """
    re = float(re)
    if re <= 0:
        return 0.02  # fallback
    # avoid log of zero
    A = (rel_roughness / 3.7) ** 1.11
    f = 0.25 / (np.log10(A + 5.74 / re ** 0.9) ** 2)
    return float(max(f, 1e-6))


# ---------------------------
# Beggs & Brill (prototipo)
# ---------------------------

def beggs_brill_dp(q_liq_stb_d: float,
                   q_gas_mscf_d: float,
                   d_in: float,
                   L_ft: float,
                   theta_deg: float,
                   rho_l: float,
                   rho_g: float,
                   mu_l: float,
                   mu_g: float,
                   roughness: float = 0.0005,
                   f_guess: float = None) -> Tuple[float, float]:
    """
    Estima ΔP total (psi) sobre una tubería de longitud L_ft usando
    una versión simplificada inspirada en Beggs & Brill.

    Retorna:
      (dp_total_psi, holdup)
    """
    # validaciones
    if L_ft <= 0 or d_in <= 0:
        raise ValueError("L_ft y d_in deben ser positivos")
    if q_liq_stb_d < 0 or q_gas_mscf_d < 0:
        raise ValueError("Los caudales no pueden ser negativos")

    # conversión de caudales a ft3/s
    ql = _stb_day_to_ft3_s(q_liq_stb_d)
    qg = _mscf_day_to_ft3_s(q_gas_mscf_d)
    d_ft = _in_to_ft(d_in)
    area = np.pi * (d_ft / 2.0) ** 2

    v_l = ql / area if area > 0 else 0.0
    v_g = qg / area if area > 0 else 0.0
    v_m = v_l + v_g

    # fracciones superficiales
    Vf_l = ql / (ql + qg) if (ql + qg) > 0 else 1.0
    Vf_l = np.clip(Vf_l, 0.0, 1.0)

    # holdup prototipo
    a = 0.25
    H = Vf_l * (1.0 + a * (v_g / (v_m + 1e-12)))
    H = float(np.clip(H, 1e-6, 0.9999))

    # densidad mixta
    rho_m = H * rho_l + (1.0 - H) * rho_g  # lb/ft3

    # Reynolds mixto aproximado (usar propiedades mixtas simplificadas)
    # viscosidad mixta en cP -> convert to lb/ft*s? keep relative for Re estimate
    mu_mix_cp = H * mu_l + (1.0 - H) * mu_g
    # conv: cP = 0.001 Pa*s ; but for rough estimate use: nu ≈ mu/(rho) in ft2/s
    # approximate kinematic viscosity ft^2/s:
    # 1 cP = 0.000672 lb/(ft·s) approximated? we'll approximate Re using v_m*d/nu using a rough nu:
    nu_ft2_s = (mu_mix_cp * 1e-3) / (rho_m * 0.03108095) if rho_m > 0 else 1e-5  # crude conversion
    Re = abs(v_m) * d_ft / (nu_ft2_s + 1e-12)

    # friction factor
    rel_rough = roughness / d_ft
    if f_guess is None:
        f = swamee_jain_f(Re, rel_rough)
    else:
        f = float(f_guess)

    # pérdida por fricción (Darcy-Weisbach -> psi)
    dp_fric_psi = f * (L_ft / d_ft) * (rho_m * (v_m ** 2) / (2.0 * G_FT_S2)) * LBFT3_TO_PSI

    # hidrostática (elevation)
    elevation_ft = L_ft * np.sin(np.deg2rad(theta_deg))
    dp_hydro_psi = (rho_m) * elevation_ft * LBFT3_TO_PSI

    dp_total_psi = float(dp_fric_psi + dp_hydro_psi)
    return dp_total_psi, H


# ---------------------------
# Hagedorn & Brown (simplificado)
# ---------------------------

def hagedorn_brown_dp(q_liq_stb_d: float,
                      q_gas_mscf_d: float,
                      d_in: float,
                      L_ft: float,
                      rho_l: float,
                      rho_g: float,
                      mu_l: float,
                      mu_g: float,
                      theta_deg: float = 90.0,
                      roughness: float = 0.0005,
                      f_guess: float = None) -> Tuple[float, float]:
    """
    Estimación simplificada de ΔP basada en Hagedorn & Brown (orientada a columnas verticales).
    """
    if L_ft <= 0 or d_in <= 0:
        raise ValueError("L_ft y d_in deben ser positivos")

    ql = _stb_day_to_ft3_s(q_liq_stb_d)
    qg = _mscf_day_to_ft3_s(q_gas_mscf_d)
    d_ft = _in_to_ft(d_in)
    area = np.pi * (d_ft / 2.0) ** 2

    v_l = ql / area if area > 0 else 0.0
    v_g = qg / area if area > 0 else 0.0
    v_m = v_l + v_g

    # Holdup empírico (prototipo)
    k = 0.15
    denom = (ql + qg) if (ql + qg) > 0 else 1.0
    H = (ql / denom) * (1.0 + k * np.sqrt(qg / (denom + 1e-12)))
    H = float(np.clip(H, 1e-6, 0.9999))

    rho_m = H * rho_l + (1.0 - H) * rho_g

    # Re y f
    mu_mix_cp = H * mu_l + (1.0 - H) * mu_g
    nu_ft2_s = (mu_mix_cp * 1e-3) / (rho_m * 0.03108095) if rho_m > 0 else 1e-5
    Re = abs(v_m) * d_ft / (nu_ft2_s + 1e-12)
    rel_rough = roughness / d_ft
    if f_guess is None:
        f = swamee_jain_f(Re, rel_rough)
    else:
        f = float(f_guess)

    dp_fric_psi = f * (L_ft / d_ft) * (rho_m * (v_m ** 2) / (2.0 * G_FT_S2)) * LBFT3_TO_PSI
    elevation_ft = L_ft * np.sin(np.deg2rad(theta_deg))
    dp_hydro_psi = (rho_m) * elevation_ft * LBFT3_TO_PSI

    dp_total_psi = float(dp_fric_psi + dp_hydro_psi)
    return dp_total_psi, H


# ---------------------------
# Vectorizados y plotting
# ---------------------------

def vlp_curve_beggs(q_liq_range_stb_d: Iterable[float],
                    q_gas_mscf_d: float,
                    d_in: float,
                    well_length_ft: float,
                    theta_deg: float,
                    rho_l: float,
                    rho_g: float,
                    mu_l: float,
                    mu_g: float,
                    roughness: float = 0.0005,
                    f_guess: float = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    q_range = np.asarray(q_liq_range_stb_d, dtype=float)
    pwf_arr = []
    H_arr = []
    for q in q_range:
        dp, H = beggs_brill_dp(q, q_gas_mscf_d, d_in, well_length_ft, theta_deg,
                                rho_l, rho_g, mu_l, mu_g, roughness, f_guess)
        pwf_arr.append(dp)
        H_arr.append(H)
    return q_range, np.array(pwf_arr), np.array(H_arr)


def vlp_curve_hagedorn(q_liq_range_stb_d: Iterable[float],
                       q_gas_mscf_d: float,
                       d_in: float,
                       well_length_ft: float,
                       rho_l: float,
                       rho_g: float,
                       mu_l: float,
                       mu_g: float,
                       theta_deg: float = 90.0,
                       roughness: float = 0.0005,
                       f_guess: float = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    q_range = np.asarray(q_liq_range_stb_d, dtype=float)
    pwf_arr = []
    H_arr = []
    for q in q_range:
        dp, H = hagedorn_brown_dp(q, q_gas_mscf_d, d_in, well_length_ft,
                                  rho_l, rho_g, mu_l, mu_g, theta_deg, roughness, f_guess)
        pwf_arr.append(dp)
        H_arr.append(H)
    return q_range, np.array(pwf_arr), np.array(H_arr)


# plotting
import matplotlib.pyplot as _plt

def plot_vlp_compare(q_liq_range, q_gas, d_in, L_ft, theta_deg, rho_l, rho_g, mu_l, mu_g, f_guess=None):
    q_arr, pwf_bb, H_bb = vlp_curve_beggs(q_liq_range, q_gas, d_in, L_ft, theta_deg, rho_l, rho_g, mu_l, mu_g, f_guess=f_guess)
    _, pwf_hb, H_hb = vlp_curve_hagedorn(q_liq_range, q_gas, d_in, L_ft, rho_l, rho_g, mu_l, mu_g, theta_deg, f_guess)

    _plt.figure(figsize=(8,6))
    _plt.plot(q_arr, pwf_bb, label='Beggs & Brill (prot)', linewidth=2)
    _plt.plot(q_arr, pwf_hb, label='Hagedorn & Brown (prot)', linewidth=2, linestyle='--')
    _plt.xlabel('Caudal líquido Q [STB/d]')
    _plt.ylabel('ΔP estimada en tubería [psi]')
    _plt.title(f'VLP: comparativa (L={L_ft} ft, d={d_in} in, angle={theta_deg}°)')
    _plt.grid(True, linestyle='--', alpha=0.6)
    _plt.legend()
    _plt.show()


# ejemplo
if __name__ == "__main__":
    q_liq = 500.0     # STB/d
    q_gas = 200.0     # MSCF/d
    d = 3.5           # in
    L = 8000.0        # ft
    theta = 90.0      # vertical
    rho_l = 60.0      # lb/ft3
    rho_g = 0.08      # lb/ft3
    mu_l = 1.0
    mu_g = 0.02

    dp_bb, H_bb = beggs_brill_dp(q_liq, q_gas, d, L, theta, rho_l, rho_g, mu_l, mu_g)
    dp_hb, H_hb = hagedorn_brown_dp(q_liq, q_gas, d, L, rho_l, rho_g, mu_l, mu_g, theta)

    print(f"Beggs&Brill (prot): ΔP = {dp_bb:.2f} psi, Holdup = {H_bb:.3f}")
    print(f"Hagedorn&Brown (prot): ΔP = {dp_hb:.2f} psi, Holdup = {H_hb:.3f}")
"""
petrokit.vlp_advanced
Versión prototipo/robusta de modelos VLP multifásicos:
- beggs_brill_dp(...) : implementación inspirada en Beggs & Brill (simplificada)
- hagedorn_brown_dp(...): implementación simplificada de Hagedorn & Brown

Salida:
- presión diferencial total estimada (psi) sobre la longitud L_ft
- holdup estimado (fracción volumétrica líquida)

NOTA:
- Implementación educativa/prototipo. Validar con datos de campo antes de tomar decisiones.
"""

from typing import Tuple, Iterable
import numpy as np

# Conversiones y constantes
STB_TO_FT3 = 5.615      # 1 STB ≈ 5.615 ft3
MSCF_TO_FT3 = 1000.0    # 1 MSCF = 1000 ft3
G_FT_S2 = 32.174        # gravedad [ft/s^2]
LBFT3_TO_PSI = 1 / 144.0  # lb/ft3 -> psi/ft factor (rho/144 -> psi/ft when multiplied by ft)


# ---------------------------
# Helpers de conversión y fricción
# ---------------------------

def _stb_day_to_ft3_s(q_stb_d: float) -> float:
    return float(q_stb_d) * STB_TO_FT3 / 86400.0

def _mscf_day_to_ft3_s(q_mscf_d: float) -> float:
    return float(q_mscf_d) * MSCF_TO_FT3 / 86400.0

def _in_to_ft(x_in: float) -> float:
    return float(x_in) / 12.0

def swamee_jain_f(re, rel_roughness):
    """
    Swamee-Jain explicit formula for friction factor f (valid for turbulent flow, Re>3000 typically).
    re: Reynolds number
    rel_roughness: epsilon/D
    """
    re = float(re)
    if re <= 0:
        return 0.02  # fallback
    # avoid log of zero
    A = (rel_roughness / 3.7) ** 1.11
    f = 0.25 / (np.log10(A + 5.74 / re ** 0.9) ** 2)
    return float(max(f, 1e-6))


# ---------------------------
# Beggs & Brill (prototipo)
# ---------------------------

def beggs_brill_dp(q_liq_stb_d: float,
                   q_gas_mscf_d: float,
                   d_in: float,
                   L_ft: float,
                   theta_deg: float,
                   rho_l: float,
                   rho_g: float,
                   mu_l: float,
                   mu_g: float,
                   roughness: float = 0.0005,
                   f_guess: float = None) -> Tuple[float, float]:
    """
    Estima ΔP total (psi) sobre una tubería de longitud L_ft usando
    una versión simplificada inspirada en Beggs & Brill.

    Retorna:
      (dp_total_psi, holdup)
    """
    # validaciones
    if L_ft <= 0 or d_in <= 0:
        raise ValueError("L_ft y d_in deben ser positivos")
    if q_liq_stb_d < 0 or q_gas_mscf_d < 0:
        raise ValueError("Los caudales no pueden ser negativos")

    # conversión de caudales a ft3/s
    ql = _stb_day_to_ft3_s(q_liq_stb_d)
    qg = _mscf_day_to_ft3_s(q_gas_mscf_d)
    d_ft = _in_to_ft(d_in)
    area = np.pi * (d_ft / 2.0) ** 2

    v_l = ql / area if area > 0 else 0.0
    v_g = qg / area if area > 0 else 0.0
    v_m = v_l + v_g

    # fracciones superficiales
    Vf_l = ql / (ql + qg) if (ql + qg) > 0 else 1.0
    Vf_l = np.clip(Vf_l, 0.0, 1.0)

    # holdup prototipo
    a = 0.25
    H = Vf_l * (1.0 + a * (v_g / (v_m + 1e-12)))
    H = float(np.clip(H, 1e-6, 0.9999))

    # densidad mixta
    rho_m = H * rho_l + (1.0 - H) * rho_g  # lb/ft3

    # Reynolds mixto aproximado (usar propiedades mixtas simplificadas)
    # viscosidad mixta en cP -> convert to lb/ft*s? keep relative for Re estimate
    mu_mix_cp = H * mu_l + (1.0 - H) * mu_g
    # conv: cP = 0.001 Pa*s ; but for rough estimate use: nu ≈ mu/(rho) in ft2/s
    # approximate kinematic viscosity ft^2/s:
    # 1 cP = 0.000672 lb/(ft·s) approximated? we'll approximate Re using v_m*d/nu using a rough nu:
    nu_ft2_s = (mu_mix_cp * 1e-3) / (rho_m * 0.03108095) if rho_m > 0 else 1e-5  # crude conversion
    Re = abs(v_m) * d_ft / (nu_ft2_s + 1e-12)

    # friction factor
    rel_rough = roughness / d_ft
    if f_guess is None:
        f = swamee_jain_f(Re, rel_rough)
    else:
        f = float(f_guess)

    # pérdida por fricción (Darcy-Weisbach -> psi)
    dp_fric_psi = f * (L_ft / d_ft) * (rho_m * (v_m ** 2) / (2.0 * G_FT_S2)) * LBFT3_TO_PSI

    # hidrostática (elevation)
    elevation_ft = L_ft * np.sin(np.deg2rad(theta_deg))
    dp_hydro_psi = (rho_m) * elevation_ft * LBFT3_TO_PSI

    dp_total_psi = float(dp_fric_psi + dp_hydro_psi)
    return dp_total_psi, H


# ---------------------------
# Hagedorn & Brown (simplificado)
# ---------------------------

def hagedorn_brown_dp(q_liq_stb_d: float,
                      q_gas_mscf_d: float,
                      d_in: float,
                      L_ft: float,
                      rho_l: float,
                      rho_g: float,
                      mu_l: float,
                      mu_g: float,
                      theta_deg: float = 90.0,
                      roughness: float = 0.0005,
                      f_guess: float = None) -> Tuple[float, float]:
    """
    Estimación simplificada de ΔP basada en Hagedorn & Brown (orientada a columnas verticales).
    """
    if L_ft <= 0 or d_in <= 0:
        raise ValueError("L_ft y d_in deben ser positivos")

    ql = _stb_day_to_ft3_s(q_liq_stb_d)
    qg = _mscf_day_to_ft3_s(q_gas_mscf_d)
    d_ft = _in_to_ft(d_in)
    area = np.pi * (d_ft / 2.0) ** 2

    v_l = ql / area if area > 0 else 0.0
    v_g = qg / area if area > 0 else 0.0
    v_m = v_l + v_g

    # Holdup empírico (prototipo)
    k = 0.15
    denom = (ql + qg) if (ql + qg) > 0 else 1.0
    H = (ql / denom) * (1.0 + k * np.sqrt(qg / (denom + 1e-12)))
    H = float(np.clip(H, 1e-6, 0.9999))

    rho_m = H * rho_l + (1.0 - H) * rho_g

    # Re y f
    mu_mix_cp = H * mu_l + (1.0 - H) * mu_g
    nu_ft2_s = (mu_mix_cp * 1e-3) / (rho_m * 0.03108095) if rho_m > 0 else 1e-5
    Re = abs(v_m) * d_ft / (nu_ft2_s + 1e-12)
    rel_rough = roughness / d_ft
    if f_guess is None:
        f = swamee_jain_f(Re, rel_rough)
    else:
        f = float(f_guess)

    dp_fric_psi = f * (L_ft / d_ft) * (rho_m * (v_m ** 2) / (2.0 * G_FT_S2)) * LBFT3_TO_PSI
    elevation_ft = L_ft * np.sin(np.deg2rad(theta_deg))
    dp_hydro_psi = (rho_m) * elevation_ft * LBFT3_TO_PSI

    dp_total_psi = float(dp_fric_psi + dp_hydro_psi)
    return dp_total_psi, H


# ---------------------------
# Vectorizados y plotting
# ---------------------------

def vlp_curve_beggs(q_liq_range_stb_d: Iterable[float],
                    q_gas_mscf_d: float,
                    d_in: float,
                    well_length_ft: float,
                    theta_deg: float,
                    rho_l: float,
                    rho_g: float,
                    mu_l: float,
                    mu_g: float,
                    roughness: float = 0.0005,
                    f_guess: float = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    q_range = np.asarray(q_liq_range_stb_d, dtype=float)
    pwf_arr = []
    H_arr = []
    for q in q_range:
        dp, H = beggs_brill_dp(q, q_gas_mscf_d, d_in, well_length_ft, theta_deg,
                                rho_l, rho_g, mu_l, mu_g, roughness, f_guess)
        pwf_arr.append(dp)
        H_arr.append(H)
    return q_range, np.array(pwf_arr), np.array(H_arr)


def vlp_curve_hagedorn(q_liq_range_stb_d: Iterable[float],
                       q_gas_mscf_d: float,
                       d_in: float,
                       well_length_ft: float,
                       rho_l: float,
                       rho_g: float,
                       mu_l: float,
                       mu_g: float,
                       theta_deg: float = 90.0,
                       roughness: float = 0.0005,
                       f_guess: float = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    q_range = np.asarray(q_liq_range_stb_d, dtype=float)
    pwf_arr = []
    H_arr = []
    for q in q_range:
        dp, H = hagedorn_brown_dp(q, q_gas_mscf_d, d_in, well_length_ft,
                                  rho_l, rho_g, mu_l, mu_g, theta_deg, roughness, f_guess)
        pwf_arr.append(dp)
        H_arr.append(H)
    return q_range, np.array(pwf_arr), np.array(H_arr)


# plotting
import matplotlib.pyplot as _plt

def plot_vlp_compare(q_liq_range, q_gas, d_in, L_ft, theta_deg, rho_l, rho_g, mu_l, mu_g, f_guess=None):
    q_arr, pwf_bb, H_bb = vlp_curve_beggs(q_liq_range, q_gas, d_in, L_ft, theta_deg, rho_l, rho_g, mu_l, mu_g, f_guess=f_guess)
    _, pwf_hb, H_hb = vlp_curve_hagedorn(q_liq_range, q_gas, d_in, L_ft, rho_l, rho_g, mu_l, mu_g, theta_deg, f_guess)

    _plt.figure(figsize=(8,6))
    _plt.plot(q_arr, pwf_bb, label='Beggs & Brill (prot)', linewidth=2)
    _plt.plot(q_arr, pwf_hb, label='Hagedorn & Brown (prot)', linewidth=2, linestyle='--')
    _plt.xlabel('Caudal líquido Q [STB/d]')
    _plt.ylabel('ΔP estimada en tubería [psi]')
    _plt.title(f'VLP: comparativa (L={L_ft} ft, d={d_in} in, angle={theta_deg}°)')
    _plt.grid(True, linestyle='--', alpha=0.6)
    _plt.legend()
    _plt.show()


# ejemplo
if __name__ == "__main__":
    q_liq = 500.0     # STB/d
    q_gas = 200.0     # MSCF/d
    d = 3.5           # in
    L = 8000.0        # ft
    theta = 90.0      # vertical
    rho_l = 60.0      # lb/ft3
    rho_g = 0.08      # lb/ft3
    mu_l = 1.0
    mu_g = 0.02

    dp_bb, H_bb = beggs_brill_dp(q_liq, q_gas, d, L, theta, rho_l, rho_g, mu_l, mu_g)
    dp_hb, H_hb = hagedorn_brown_dp(q_liq, q_gas, d, L, rho_l, rho_g, mu_l, mu_g, theta)

    print(f"Beggs&Brill (prot): ΔP = {dp_bb:.2f} psi, Holdup = {H_bb:.3f}")
    print(f"Hagedorn&Brown (prot): ΔP = {dp_hb:.2f} psi, Holdup = {H_hb:.3f}")
