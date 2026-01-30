# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Phase 3 - Profesionalización ✅

#### Added
- GitHub Actions CI/CD workflow for automated testing
  - Multi-platform support (Ubuntu, Windows, macOS)
  - Multiple Python versions (3.10, 3.11, 3.12)
  - Code coverage integration with Codecov
- Professional Sphinx documentation
  - Complete API reference with autodoc
  - Installation guide
  - Quick start tutorial
  - Examples gallery
  - Contributing guidelines
- ReadTheDocs configuration (`.readthedocs.yml`)
- Interactive visualization examples
  - Plotly interactive charts (`examples/interactive_plotly_visualizations.ipynb`)
  - ipywidgets parametric sensitivity analysis (`examples/parametric_sensitivity_interactive.ipynb`)
- Professional package configuration
  - Updated `pyproject.toml` with proper metadata and URLs
  - Updated `setup.py` for PyPI distribution
  - Development dependencies properly organized

#### Fixed
- `build_pvt_table()` now returns consistent array types for all outputs
- Sphinx documentation build warnings (missing _static directory)

#### Documentation
- Complete Phase 3 summary (`docs/FASE3_RESUMEN.md`)
- Enhanced README with professional badges
- Added this CHANGELOG

### Phase 2 - Ampliación de Modelos ✅

#### Added
- Extended IPR models: Jones (fractured), Standing (gas solution)
- VLP correlations: Beggs & Brill, Hagedorn & Brown
- Multiphase pressure drop in pipes and flowlines
- Inclination angle and flow regime considerations
- Simplified PVT tables (integrated via black-oil VLP)
- Unit conversion utilities (STB ↔ m³, psi ↔ Pa)
- Additional examples for oil and gas wells
- Correlation comparison notebooks

### Phase 1 - MVP ✅

#### Added
- Core package structure
  - `petrokit/ipr.py` - IPR models (Vogel, Fetkovich)
  - `petrokit/vlp.py` - Simplified VLP models (Darcy-Weisbach)
  - `petrokit/flowline.py` - Pressure drop in flowlines
  - `petrokit/nodal.py` - IPR-VLP intersection
  - `petrokit/utils.py` - Unit conversions and helpers
- Unit tests with pytest
- Basic examples in Spanish (`examples/analisis_nodal_español.ipynb`)
- Initial README with installation instructions
- MIT License
- `.gitignore` configuration

## [0.1.0] - 2025-01

### Added
- Initial release
- Basic IPR and VLP calculations
- Nodal analysis functionality
- Unit conversion utilities
- Comprehensive test suite
- Documentation and examples

[Unreleased]: https://github.com/OrtegonAutomation/petrokit/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/OrtegonAutomation/petrokit/releases/tag/v0.1.0
