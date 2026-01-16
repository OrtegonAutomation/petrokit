"""
PetroKit: Herramientas para Ingeniería de Petróleos
Módulos disponibles:
- ipr: curvas de IPR (Vogel, Fetkovich, etc.)
- vlp: curvas de VLP
- flowline: pérdidas de presión en tuberías
- nodal: análisis nodal
- utils: utilidades generales
"""

# Importar funciones principales directamente
from .ipr import vogel_ipr, fetkovich_ipr, ipr_curve_vogel, plot_ipr_vogel
from .vlp import *       # en caso que ya tengas funciones ahí
from .flowline import *  # idem
from .nodal import *     # idem
from .utils import *     # idem

from .ipr import jones_ipr, standing_ipr, vogel_ipr, fetkovich_ipr
from .vlp import vlp_curve_model, available_vlp_models
from .nodal import nodal_analysis, nodal_analysis_detail

__all__ = [
    "jones_ipr", "standing_ipr", "vogel_ipr", "fetkovich_ipr",
    "vlp_curve_model", "available_vlp_models",
    "nodal_analysis", "nodal_analysis_detail",
]
