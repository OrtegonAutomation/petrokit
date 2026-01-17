Contribuir a PetroKit
=====================

¡Gracias por tu interés en contribuir a PetroKit! Este documento te guiará a través del proceso.

Cómo contribuir
---------------

1. **Fork del repositorio**

   Haz un fork del repositorio en GitHub:
   
   https://github.com/OrtegonAutomation/petrokit

2. **Clonar tu fork**

   .. code-block:: bash
   
      git clone https://github.com/tu-usuario/petrokit.git
      cd petrokit

3. **Crear una rama**

   .. code-block:: bash
   
      git checkout -b feature/nueva-funcionalidad

4. **Instalar dependencias de desarrollo**

   .. code-block:: bash
   
      pip install -e ".[dev]"

5. **Hacer cambios**

   * Añade tu código
   * Añade tests para tu código
   * Actualiza la documentación si es necesario

6. **Ejecutar tests**

   .. code-block:: bash
   
      pytest -v

7. **Commit y push**

   .. code-block:: bash
   
      git add .
      git commit -m "Descripción clara de los cambios"
      git push origin feature/nueva-funcionalidad

8. **Abrir Pull Request**

   Abre un Pull Request en GitHub describiendo tus cambios.

Guía de estilo
--------------

* Sigue PEP 8 para el estilo de código Python
* Usa docstrings para documentar funciones y clases
* Añade type hints cuando sea posible
* Mantén las líneas de código bajo 100 caracteres

Escribir tests
--------------

Todos los cambios deben incluir tests:

.. code-block:: python

   # tests/test_mi_modulo.py
   
   def test_mi_funcion():
       resultado = mi_funcion(parametro=10)
       assert resultado > 0
       assert resultado < 100

Ejecutar tests localmente:

.. code-block:: bash

   # Todos los tests
   pytest -v
   
   # Un solo archivo
   pytest tests/test_mi_modulo.py
   
   # Con cobertura
   pytest --cov=petrokit --cov-report=html

Documentación
-------------

Si añades nuevas funcionalidades, actualiza la documentación:

1. Añade docstrings en el código
2. Actualiza los archivos .rst en ``docs/source/``
3. Añade ejemplos si es necesario

Para generar la documentación localmente:

.. code-block:: bash

   cd docs
   make html
   # Abre docs/build/html/index.html en tu navegador

Reportar bugs
-------------

Si encuentras un bug:

1. Verifica que no esté ya reportado en Issues
2. Crea un nuevo Issue con:
   
   * Descripción clara del problema
   * Código para reproducir el bug
   * Comportamiento esperado vs. actual
   * Versión de Python y PetroKit

Sugerir mejoras
---------------

Para sugerir nuevas funcionalidades:

1. Abre un Issue describiendo:
   
   * La funcionalidad propuesta
   * Casos de uso
   * Ejemplos de cómo se usaría

2. Espera feedback antes de implementar

Código de conducta
------------------

* Sé respetuoso con otros colaboradores
* Acepta críticas constructivas
* Enfócate en lo mejor para el proyecto
* Usa un lenguaje inclusivo

Licencia
--------

Al contribuir a PetroKit, aceptas que tus contribuciones se licenciarán bajo la licencia MIT del proyecto.
