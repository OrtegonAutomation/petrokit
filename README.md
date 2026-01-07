
# PetroKit — Biblioteca Python para Ingeniería de Producción y Transporte de Hidrocarburos

![PyPI - Placeholder](https://img.shields.io/badge/PyPI-pending-lightgrey)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

**PetroKit** es una librería en Python diseñada para ingenieros de petróleos y científicos de datos que trabajan en **producción y transporte de hidrocarburos**.

Su objetivo es proporcionar un **entorno abierto, reproducible y extensible** que permita:

* Modelar curvas **IPR** (Inflow Performance Relationship).
* Calcular curvas **VLP** (Vertical Lift Performance).
* Estimar **pérdidas de presión en flowlines**.
* Integrar ambos en un **Análisis Nodal**.
* Usar funciones auxiliares de conversión de unidades y números adimensionales.

---

## Tabla de contenidos

1. Visión general
2. Instalación
3. Conceptos físicos y fórmulas
4. API principal
5. Ejemplo práctico
6. Validación y pruebas
7. Buenas prácticas de uso
8. Roadmap técnico
9. Integración y CI/CD
10. Contribuciones y citación
11. Licencia

---

## Visión general

PetroKit nace como un **toolkit de producción** para cubrir la brecha entre:

* **Scripts aislados en Excel/Python** (no reproducibles ni escalables).
* **Software comercial cerrado** (costo alto, poca transparencia).

### En la Fase 1 (MVP):

* Modelos básicos de IPR (Vogel, Fetkovich).
* VLP simplificada con Darcy–Weisbach.
* Nodal analysis básico.
* Flowline monofásico con fricción constante.
* Pytest y ejemplos en Jupyter.

### Fases futuras:

* Correlaciones multifásicas (Beggs & Brill, Hagedorn & Brown).
* Modelos de levantamiento artificial (gas lift, ESP, PCP).
* Dashboards interactivos y API REST.

---

## Instalación

Clonar e instalar en modo editable:

```bash
git clone https://github.com/OrtegonAutomation/petrokit.git
cd petrokit
pip install -e .
```

Dependencias principales:

* `numpy`
* `matplotlib`
* `pytest` (solo desarrollo)

---

## Conceptos físicos y fórmulas

**IPR — Vogel**

$$
q = q_{max} \left(1 - 0.2 \frac{p_{wf}}{p_{res}} - 0.8 \left(\frac{p_{wf}}{p_{res}}\right)^2\right)
$$

**IPR — Fetkovich**

$$
q = J \cdot (p_{res} - p_{wf})
$$

**IPR — Jones (no-Darcy / fracturados)**

$$ \Delta p = p_{res} - p_{wf} = C q + D q^2 $$

(IPR se obtiene resolviendo la raíz positiva de la cuadrática en q).

**IPR — Standing (gas solution)**

Para $p_{wf} \ge p_b$:

$$ q = J \cdot (p_{res} - p_{wf}) $$

Para $p_{wf} < p_b$ (dos fases), se usa Vogel calibrando $q_{max}$ para continuidad en $p_b$:

$$ q = q_{max}\left(1 - 0.2\frac{p_{wf}}{p_{res}} - 0.8\left(\frac{p_{wf}}{p_{res}}\right)^2\right) $$


**VLP — Darcy–Weisbach (simplificado, monofásico)**

$$
\Delta p_{fric} = f \cdot \frac{L}{D} \cdot \frac{\rho v^2}{2}
$$

**Columna estática**

$$
\Delta p_{hydro} = \rho g H
$$

---

## API principal

### petrokit.ipr

  * `jones_ipr(p_res, C, D, pwf)` → caudal con no-Darcy (C·q + D·q²).
  * `standing_ipr(p_res, p_b, J, pwf)` → IPR compuesto (sobre/bajo bubble point).
  * `ipr_curve_fetkovich(p_res, J, npts=50)` → arrays pwf, q.
  * `ipr_curve_jones(p_res, C, D, npts=50)` → arrays pwf, q.
  * `ipr_curve_standing(p_res, p_b, J, npts=50)` → arrays pwf, q.
  * `plot_ipr_fetkovich(...)`, `plot_ipr_jones(...)`, `plot_ipr_standing(...)` → gráficos.


### petrokit.vlp

* `vlp_curve(q_range, well_depth, rho, mu, d, f=0.02)` → pwf \[psi].
* `plot_vlp(...)` → gráfico.

### petrokit.flowline

* `flowline_pressure_drop(q, L, d, rho, mu, f=0.02, elev=0)` → ΔP \[psi].
* `plot_flowline(...)` → curva ΔP vs q.

### petrokit.nodal

* `nodal_analysis(p_res, q_max, well_depth, rho, mu, d, npts=50)` → (q\_op, pwf\_op).
* `plot_nodal(...)` → intersección IPR–VLP.

### petrokit.utils

* Conversión de unidades: `psi_to_pa`, `stb_to_m3`, etc.
* `reynolds_number(q, d, mu, rho)` → número de Reynolds.

---

## Ejemplo práctico

```python
from petrokit.ipr import plot_ipr_vogel
from petrokit.nodal import nodal_analysis, plot_nodal

p_res = 3000       # psi
q_max = 1200       # STB/d
well_depth = 8000  # ft
rho = 60           # lb/ft³
mu = 1             # cP
d = 2.992          # in

# Graficar IPR
plot_ipr_vogel(p_res, q_max)

# Calcular punto de operación (Nodal Analysis)
q_op, pwf_op = nodal_analysis(p_res, q_max, well_depth, rho, mu, d)
print(f"Punto de operación: Q = {q_op:.2f} STB/d, pwf = {pwf_op:.2f} psi")

# Graficar nodal completo
plot_nodal(p_res, q_max, well_depth, rho, mu, d)
```

**Salida esperada:**

```
Punto de operación: Q ≈ 850 STB/d, pwf ≈ 1200 psi
```

Revisa el notebook `examples/analisis_nodal_español.ipynb` para un estudio completo con gráficas.

---

## Validación y pruebas

Ejecutar pruebas unitarias:

```bash
pytest -v
```

Cobertura en Fase 1:

* IPR: condiciones límite (pwf=0, pwf=p_res) + casos Jones/Standing.
* VLP: monotonicidad, valores positivos.
* Flowline: efecto de elevación y longitud.
* Nodal: punto válido dentro de rango.

---

## Buenas prácticas de uso

* Mantener **coherencia de unidades**: psi, STB/d, ft, in, lb/ft³.
* Evitar parámetros no físicos (`densidad ≤ 0`, `diámetro ≤ 0`).
* Para sensibilidades usa vectores `numpy`.
* Verificar siempre contra datos de campo o software comercial antes de usar en proyectos reales.

---


# 🛣️ Roadmap del Proyecto **PetroKit**

---

## 🔹 Fase 1: **MVP (Producto Mínimo Viable)**

🎯 Objetivo: Tener un paquete funcional, probado y con ejemplos básicos.

1. **Estructura del paquete**

   * `petrokit/ipr.py` → Modelos IPR (Vogel, Fetkovich).
   * `petrokit/vlp.py` → Modelos VLP simplificados (Darcy–Weisbach).
   * `petrokit/flowline.py` → Caída de presión en líneas de flujo.
   * `petrokit/nodal.py` → Intersección IPR-VLP.
   * `petrokit/utils.py` → Conversión de unidades y helpers.
   * `tests/` → Pruebas unitarias con `pytest`.
   * `examples/` → Notebooks ilustrativos en español.
   * `README.md` + `requirements.txt` + `LICENSE` + `.gitignore`.

2. **Tests unitarios**

   * Validar ecuaciones básicas.
   * Sensibilidad de parámetros (ejemplo: q\_max, diámetro).
   * Integración básica (análisis nodal completo).

3. **Documentación inicial**

   * README con instalación y ejemplos.
   * Un `examples/nodal_analysis.ipynb` en español.


---

## 🔹 Fase 2: **Ampliación de Modelos**

🎯 Objetivo: Pasar de un demo académico a un toolkit más completo.

1. **Producción**

   * Extender modelos IPR: Jones (fracturados) ✅, Standing (gas solution) ✅.
   * VLP con correlaciones: Beggs & Brill, Hagedorn & Brown.

2. **Transporte**

   * Caída de presión multifásica en tuberías y flowlines.
   * Consideración de ángulo de inclinación y régimen de flujo.

3. **Utilidades**

   * Tablas PVT simplificadas.
   * Conversión entre unidades (STB ↔ m³, psi ↔ Pa).

4. **Más ejemplos**

   * Casos de pozos de petróleo y gas.
   * Comparación de correlaciones en notebooks.

---

## 🔹 Fase 3: **Profesionalización**

🎯 Objetivo: Convertirlo en un paquete distribuible e instalable con `pip`.

1. **Infraestructura**

   * Publicación en **PyPI** (`pip install petrokit`).
   * Configuración de `pyproject.toml` o `setup.py`.
   * GitHub Actions para correr `pytest` en cada commit.

2. **Documentación profesional**

   * Página web con **Sphinx + ReadTheDocs**.
   * Ejemplos con gráficos interactivos (Plotly).

3. **Ejemplos avanzados**

   * Sensibilidad paramétrica con `ipywidgets`.
   * Comparación entre varios escenarios (e.g., tubing vs flowline).

---

## 🔹 Fase 4: **Integración Industrial**

🎯 Objetivo: Acercarse a la práctica profesional de ingeniería de petróleo.

1. **Modelos de superficie**

   * Redes de flujo (múltiples pozos + manifold).
   * Simulación de facilities (presión en separadores).

2. **PVT más robusto**

   * Ecuaciones de estado (Peng–Robinson, SRK).
   * Matching de laboratorio.

3. **Optimización**

   * Maximización de producción (ajuste de choke, diámetro de tubing).
   * Algoritmos de sensibilidad y escenarios.

4. **Exportación de resultados**

   * Generar reportes en PDF/Excel.
   * Compatibilidad con simuladores comerciales (exportación de datos).

---

## 🔹 Fase 5: **Escalamiento & Comunidad**

🎯 Objetivo: Convertir a **PetroKit** en un proyecto colaborativo y escalable.

1. **Colaboración**

   * Abrir Issues y PRs en GitHub.
   * Crear documentación clara para contribuir.

2. **Integración Data Science**

   * Machine Learning para predicción de IPR.
   * Modelos de regresión para curvas PVT.

3. **Visualización avanzada**

   * Dashboards con **Streamlit** o **Dash**.
   * Gráficas interactivas 3D.

4. **Casos reales**

   * Datasets sintéticos y de la literatura.
   * Benchmarking contra simuladores industriales.

---

# 🚀 Visión Final

PetroKit podría evolucionar en:
✅ Un **toolbox académico** para estudiantes de Ingeniería de Petróleos.
✅ Un **prototipo industrial** para ingenieros de producción.
✅ Una **plataforma open-source** de referencia en simulaciones de producción y transporte.


---

## Integración y CI/CD

* **Tests automáticos**: GitHub Actions (`.github/workflows/tests.yml`).
* **Empaquetado PyPI**: `pyproject.toml` + `twine upload`.
* **Documentación**: recomendación `mkdocs` o `sphinx` para docs técnicas.

---

## Contribuciones y citación

### Cómo contribuir

1. Haz un fork.
2. Crea una rama: `git checkout -b feature/nueva-funcionalidad`.
3. Haz cambios y añade tests.
4. Abre un Pull Request.

### Cita PetroKit

Si usas PetroKit en investigación:

```bibtex
@software{petrokit2025,
  author = {Camilo Andrés Ortegon Cuevas},
  title = {Librería para análisis de transporte y producción en ingeniería de petróleos},
  year = {2025},
  url = {https://github.com/OrtegonAutomation/petrokit},
}
```

---

## Licencia

Este proyecto está bajo licencia MIT.
Eres libre de usarlo, modificarlo y distribuirlo citando la fuente.

---
