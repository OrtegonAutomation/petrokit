Instalación
===========

Requisitos
----------

* Python >= 3.10
* pip

Instalación desde PyPI
----------------------

Una vez publicado en PyPI, podrás instalar PetroKit con:

.. code-block:: bash

   pip install petrokit

Instalación para desarrollo
----------------------------

Para instalar en modo editable con dependencias de desarrollo:

.. code-block:: bash

   git clone https://github.com/OrtegonAutomation/petrokit.git
   cd petrokit
   pip install -e ".[dev]"

Verificar la instalación
-------------------------

Para verificar que PetroKit está correctamente instalado:

.. code-block:: python

   import petrokit
   from petrokit.ipr import ipr_curve_fetkovich
   
   # Si no hay errores, la instalación fue exitosa
   print("PetroKit instalado correctamente")

Dependencias
------------

Dependencias principales:

* ``numpy``: Cálculos numéricos
* ``matplotlib``: Visualización básica
* ``plotly``: Gráficos interactivos
* ``ipywidgets``: Widgets interactivos para Jupyter

Dependencias de desarrollo:

* ``pytest``: Testing
* ``pytest-cov``: Cobertura de tests
* ``sphinx``: Documentación
* ``sphinx-rtd-theme``: Tema ReadTheDocs
