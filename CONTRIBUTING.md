# Guía de Contribución a PetroKit

¡Gracias por tu interés en contribuir a PetroKit! Este documento proporciona pautas para contribuir al proyecto.

## Tabla de Contenidos

1. [Código de Conducta](#código-de-conducta)
2. [Cómo Empezar](#cómo-empezar)
3. [Proceso de Desarrollo](#proceso-de-desarrollo)
4. [Estándares de Código](#estándares-de-código)
5. [Tests](#tests)
6. [Documentación](#documentación)
7. [Envío de Pull Requests](#envío-de-pull-requests)

## Código de Conducta

Este proyecto se adhiere a estándares de comportamiento profesional y respetuoso. Se espera que todos los contribuyentes:

- Sean respetuosos con otros colaboradores
- Acepten críticas constructivas
- Se enfoquen en lo que es mejor para la comunidad
- Muestren empatía hacia otros miembros de la comunidad

## Cómo Empezar

### Configuración del Entorno de Desarrollo

1. **Fork el repositorio**
   ```bash
   # En GitHub, haz clic en "Fork"
   ```

2. **Clona tu fork**
   ```bash
   git clone https://github.com/TU-USUARIO/petrokit.git
   cd petrokit
   ```

3. **Instala las dependencias de desarrollo**
   ```bash
   pip install -e ".[dev]"
   ```

4. **Verifica que todo funciona**
   ```bash
   pytest -v
   ```

### Estructura del Proyecto

```
petrokit/
├── petrokit/           # Código fuente principal
│   ├── ipr.py         # Modelos IPR
│   ├── vlp.py         # Modelos VLP
│   ├── flowline.py    # Cálculos de flowline
│   ├── nodal.py       # Análisis nodal
│   ├── pvt.py         # Propiedades PVT
│   └── utils.py       # Utilidades
├── tests/             # Tests unitarios
├── examples/          # Notebooks de ejemplo
├── docs/              # Documentación Sphinx
└── README.md          # Documentación principal
```

## Proceso de Desarrollo

### 1. Crear una Rama

```bash
git checkout -b feature/nueva-funcionalidad
# o
git checkout -b bugfix/correccion-error
```

**Convención de nombres de ramas:**
- `feature/` - Nueva funcionalidad
- `bugfix/` - Corrección de errores
- `docs/` - Cambios en documentación
- `test/` - Añadir o mejorar tests

### 2. Hacer Cambios

- Haz commits pequeños y atómicos
- Escribe mensajes de commit descriptivos
- Sigue las convenciones de código del proyecto

### 3. Ejecutar Tests

```bash
# Ejecutar todos los tests
pytest -v

# Ejecutar tests con cobertura
pytest -v --cov=petrokit --cov-report=term

# Ejecutar tests específicos
pytest tests/test_ipr.py -v
```

### 4. Actualizar Documentación

Si añades nueva funcionalidad:
- Actualiza docstrings en el código
- Añade ejemplos en `examples/` si es apropiado
- Actualiza `docs/source/api.rst` si es necesario

## Estándares de Código

### Estilo Python

- Seguir [PEP 8](https://peps.python.org/pep-0008/)
- Usar nombres descriptivos para variables y funciones
- Mantener funciones pequeñas y enfocadas
- Añadir docstrings a todas las funciones públicas

### Formato de Docstrings

Usar el formato Google/NumPy:

```python
def funcion_ejemplo(param1: float, param2: str) -> float:
    """
    Breve descripción de la función.

    Descripción más detallada si es necesario.

    Args:
        param1: Descripción del primer parámetro
        param2: Descripción del segundo parámetro

    Returns:
        Descripción del valor retornado

    Raises:
        ValueError: Cuando ocurre este error

    Example:
        >>> resultado = funcion_ejemplo(1.0, "test")
        >>> print(resultado)
        1.0
    """
    return param1
```

### Convenciones Específicas de PetroKit

- **Unidades**: Usar unidades del campo petrolero (psi, STB/d, ft, °F) por defecto
- **Nomenclatura**: 
  - `p_res` para presión de reservorio
  - `pwf` para presión de fondo fluyente
  - `q` para caudal
  - `rho` para densidad
  - `mu` para viscosidad

## Tests

### Escribir Tests

- Todo código nuevo debe incluir tests
- Los tests deben ser claros y descriptivos
- Probar casos límite y condiciones de error

```python
def test_funcion_ejemplo():
    """Test que verifica comportamiento básico."""
    resultado = funcion_ejemplo(1.0)
    assert resultado > 0
    assert isinstance(resultado, float)

def test_funcion_ejemplo_error():
    """Test que verifica manejo de errores."""
    with pytest.raises(ValueError):
        funcion_ejemplo(-1.0)
```

### Cobertura de Tests

- Apuntar a >80% de cobertura de código
- Verificar cobertura con: `pytest --cov=petrokit --cov-report=html`
- Ver reporte en `htmlcov/index.html`

## Documentación

### Construir la Documentación Localmente

```bash
cd docs
make html
# Abrir docs/build/html/index.html en tu navegador
```

### Actualizar la Documentación

1. Actualizar docstrings en el código
2. Modificar archivos `.rst` en `docs/source/` si es necesario
3. Reconstruir la documentación
4. Verificar que no hay warnings

## Envío de Pull Requests

### Antes de Enviar

Verifica que:
- [ ] Todos los tests pasan: `pytest -v`
- [ ] El código sigue los estándares de estilo
- [ ] Añadiste tests para nueva funcionalidad
- [ ] Actualizaste la documentación
- [ ] Tu rama está actualizada con `main`

### Proceso de PR

1. **Push a tu fork**
   ```bash
   git push origin feature/nueva-funcionalidad
   ```

2. **Abrir Pull Request en GitHub**
   - Título descriptivo
   - Descripción clara de los cambios
   - Referencia a issues relacionados (#123)

3. **Descripción del PR debe incluir:**
   - ¿Qué cambia este PR?
   - ¿Por qué es necesario?
   - ¿Cómo se probó?
   - Screenshots si aplica (para cambios visuales)

### Template de PR

```markdown
## Descripción
Breve descripción de los cambios

## Tipo de cambio
- [ ] Bug fix (cambio que corrige un issue)
- [ ] Nueva funcionalidad (cambio que añade funcionalidad)
- [ ] Breaking change (fix o feature que causa que funcionalidad existente no funcione como antes)
- [ ] Documentación

## ¿Cómo se ha probado?
Describe las pruebas realizadas

## Checklist
- [ ] Mi código sigue el estilo del proyecto
- [ ] He realizado una auto-revisión de mi código
- [ ] He comentado mi código en áreas difíciles de entender
- [ ] He actualizado la documentación
- [ ] Mis cambios no generan nuevos warnings
- [ ] He añadido tests que prueban que mi fix es efectivo o que mi feature funciona
- [ ] Los tests nuevos y existentes pasan localmente
```

## Preguntas

Si tienes preguntas, puedes:
- Abrir un [Issue](https://github.com/OrtegonAutomation/petrokit/issues)
- Revisar la [documentación](https://petrokit.readthedocs.io)
- Contactar a los mantenedores

## Agradecimientos

¡Gracias por contribuir a PetroKit! Tu ayuda hace que este proyecto sea mejor para toda la comunidad de ingeniería de petróleos.

---

**Nota**: Esta guía está en constante evolución. Si tienes sugerencias para mejorarla, ¡no dudes en proponer cambios!
