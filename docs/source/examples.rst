Ejemplos
========

PetroKit incluye varios notebooks de ejemplo que demuestran diferentes aspectos de la librería.

Notebooks disponibles
---------------------

Análisis Nodal
~~~~~~~~~~~~~~

* **analisis_nodal_español.ipynb**: Ejemplo completo de análisis nodal en español
* **nodal_analysis_english.ipynb**: Complete nodal analysis example in English

IPR Avanzado
~~~~~~~~~~~~

* **ipr_avanzado.ipynb**: Modelos IPR avanzados (Jones, Standing)
* **ipr_advanced_english.ipynb**: Advanced IPR models (Jones, Standing)

VLP Avanzado
~~~~~~~~~~~~

* **vlp_avanzado.ipynb**: Correlaciones VLP avanzadas (Beggs & Brill, Hagedorn & Brown)
* **vlp_advanced_english.ipynb**: Advanced VLP correlations

PVT y Flowlines
~~~~~~~~~~~~~~~

* **ingenieria_pvt_flowline.ipynb**: Propiedades PVT y análisis de flowlines
* **engineering_pvt_flowline_english.ipynb**: PVT properties and flowline analysis

Acceder a los ejemplos
-----------------------

Los notebooks se encuentran en el directorio ``examples/`` del repositorio:

.. code-block:: bash

   cd examples
   jupyter notebook

Ejemplos interactivos con Plotly
---------------------------------

Los ejemplos más recientes incluyen visualizaciones interactivas usando Plotly:

.. code-block:: python

   import plotly.graph_objects as go
   from petrokit.nodal import nodal_analysis
   from petrokit.ipr import ipr_curve_vogel
   from petrokit.vlp import vlp_curve
   import numpy as np
   
   # Parámetros
   p_res = 3000
   q_max = 1200
   well_depth = 8000
   rho = 60
   mu = 1
   d = 2.992
   
   # Calcular curvas
   pwf_ipr, q_ipr = ipr_curve_vogel(p_res, q_max)
   q_vlp = np.linspace(0, 1200, 50)
   pwf_vlp = vlp_curve(q_vlp, well_depth, rho, mu, d)
   
   # Crear gráfico interactivo
   fig = go.Figure()
   fig.add_trace(go.Scatter(x=q_ipr, y=pwf_ipr, name='IPR', mode='lines'))
   fig.add_trace(go.Scatter(x=q_vlp, y=pwf_vlp, name='VLP', mode='lines'))
   
   fig.update_layout(
       title='Análisis Nodal Interactivo',
       xaxis_title='Caudal (STB/d)',
       yaxis_title='Presión (psi)',
       hovermode='x unified'
   )
   
   fig.show()

Análisis de sensibilidad con ipywidgets
----------------------------------------

Ejemplos con widgets interactivos para análisis de sensibilidad paramétrica:

.. code-block:: python

   from ipywidgets import interact, FloatSlider
   import matplotlib.pyplot as plt
   from petrokit.ipr import ipr_curve_vogel
   
   def plot_sensitivity(q_max):
       p_res = 3000
       pwf, q = ipr_curve_vogel(p_res, q_max)
       
       plt.figure(figsize=(10, 6))
       plt.plot(q, pwf, linewidth=2)
       plt.xlabel('Caudal (STB/d)')
       plt.ylabel('Presión de fondo (psi)')
       plt.title(f'IPR - Vogel (q_max = {q_max} STB/d)')
       plt.grid(True, alpha=0.3)
       plt.show()
   
   # Crear widget interactivo
   interact(plot_sensitivity, 
            q_max=FloatSlider(min=500, max=2000, step=100, value=1200,
                             description='q_max (STB/d):'))
