Inicio Rápido
=============

Este tutorial te guiará a través de los conceptos básicos de PetroKit.

Importar módulos
----------------

.. code-block:: python

   from petrokit.ipr import plot_ipr_vogel, ipr_curve_fetkovich
   from petrokit.vlp import vlp_curve
   from petrokit.nodal import nodal_analysis, plot_nodal
   from petrokit.utils import psi_to_pa, stb_to_m3

Ejemplo básico: IPR
-------------------

Calcular y graficar una curva IPR usando el modelo de Vogel:

.. code-block:: python

   from petrokit.ipr import plot_ipr_vogel
   
   p_res = 3000  # psi - Presión de reservorio
   q_max = 1200  # STB/d - Caudal máximo
   
   # Graficar IPR
   plot_ipr_vogel(p_res, q_max)

Ejemplo básico: VLP
-------------------

Calcular la presión de fondo para diferentes caudales:

.. code-block:: python

   import numpy as np
   from petrokit.vlp import vlp_curve
   
   q_range = np.linspace(0, 1200, 50)  # Rango de caudales [STB/d]
   well_depth = 8000  # ft - Profundidad del pozo
   rho = 60           # lb/ft³ - Densidad del fluido
   mu = 1             # cP - Viscosidad
   d = 2.992          # in - Diámetro interno del tubing
   
   pwf = vlp_curve(q_range, well_depth, rho, mu, d)

Análisis Nodal completo
-----------------------

Combinar IPR y VLP para encontrar el punto de operación:

.. code-block:: python

   from petrokit.nodal import nodal_analysis, plot_nodal
   
   p_res = 3000       # psi
   q_max = 1200       # STB/d
   well_depth = 8000  # ft
   rho = 60           # lb/ft³
   mu = 1             # cP
   d = 2.992          # in
   
   # Calcular punto de operación
   q_op, pwf_op = nodal_analysis(p_res, q_max, well_depth, rho, mu, d)
   print(f"Punto de operación: Q = {q_op:.2f} STB/d, pwf = {pwf_op:.2f} psi")
   
   # Graficar nodal completo
   plot_nodal(p_res, q_max, well_depth, rho, mu, d)

Modelos avanzados
-----------------

PetroKit también soporta modelos más avanzados:

.. code-block:: python

   from petrokit.ipr import ipr_curve_jones, ipr_curve_standing
   from petrokit.vlp import vlp_curve_model, available_vlp_models
   import numpy as np
   
   # Ver modelos VLP disponibles
   print(available_vlp_models())  # ['darcy', 'beggs_brill', 'hagedorn_brown', ...]
   
   # IPR Jones (pozos fracturados - no-Darcy flow)
   pwf_jones, q_jones = ipr_curve_jones(p_res=3000, C=0.5, D=0.001)
   
   # IPR Standing (gas solution - compuesto sobre/bajo bubble point)
   pwf_standing, q_standing = ipr_curve_standing(
       p_res=3000, p_b=2200, J=1.5
   )
   
   # VLP con Beggs & Brill (flujo multifásico)
   q_range = np.linspace(0, 1200, 50)
   pwf_bb = vlp_curve_model(
       "beggs_brill", q_range, 
       well_depth=8000, rho=60, mu=1, d=2.992,
       ql_bpd=800,  # Caudal de líquido [bbl/d]
       qg_scfd=50000,  # Caudal de gas [scf/d]
       angle_deg=90  # Ángulo vertical [grados]
   )

Próximos pasos
--------------

* Explora la :doc:`api` completa
* Revisa los :doc:`examples` para casos más complejos
* Lee sobre cómo :doc:`contributing` al proyecto
