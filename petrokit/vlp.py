import numpy as np
import matplotlib.pyplot as plt

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
    # Conversión de caudal a ft3/s
    stb_to_ft3 = 5.615 / (24*3600)  # 1 STB/d = 5.615 ft3/d
    q_ft3_s = q_range * stb_to_ft3
    
    # Área de la tubería [ft2]
    d_ft = d / 12
    A = np.pi * (d_ft/2)**2
    
    # Velocidad [ft/s]
    v = q_ft3_s / A
    
    # Head loss por fricción
    g = 32.174  # ft/s2
    delta_p_fric = f * (well_depth/d_ft) * (rho/144) * (v**2/2/g)
    
    # Head estático (columna fluido)
    delta_p_hydro = (rho/144) * well_depth
    
    # Presión en fondo fluyente [psi]
    pwf_range = delta_p_fric + delta_p_hydro
    
    return pwf_range

def plot_vlp(q_range, well_depth, rho, mu, d, f=0.02):
    """
    Genera y grafica la curva VLP.
    """
    pwf_range = vlp_curve(q_range, well_depth, rho, mu, d, f)
    plt.figure(figsize=(6,5))
    plt.plot(q_range, pwf_range, color="red", label="VLP")
    plt.xlabel("Caudal [STB/d]")
    plt.ylabel("Presión Fondo Fluyente pwf [psi]")
    plt.title("Curva VLP")
    plt.grid(True)
    plt.legend()
    plt.show()
