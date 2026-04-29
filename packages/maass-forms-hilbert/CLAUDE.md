# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

### Installation and Building
- `make install` - Clean build, create source distribution, and install with sage
- `make build` - Build Cython extensions in-place using `sage -python setup.py build_ext --inplace`
- `make sdist` - Create source distribution
- `make clean` - Remove build artifacts, compiled extensions, and temporary files

### Testing
- `make test` - Run Sage doctests on src/ directory (`sage -t src/*`)
- `pytest` - Run pytest tests (requires `sage -i pytest pytest_xdist` first)
- `tox -e doctest` - Run Sage doctester with tox
- `tox -e coverage` - Get doctest coverage information
- `tox -e pytest` - Run pytest via tox

### Code Quality
- `ruff check .` - Run linting checks
- `ruff format .` - Format code
- `tox -e ruff` - Run both ruff check and format via tox
- `tox -e codespell` - Check for misspelled words
- `tox -e relint` - Check for forbidden patterns

### Full Test Suite
- `make tox` - Install tox and run full test suite including doctests, coverage, and linting

## Architecture Overview

This is a Python package for algorithms related to Hilbert Maass forms, built on top of SageMath/PassageMath.

### Core Structure
- **`src/maass_forms_hilbert/`** - Main package directory
  - **`functions/`** - Mathematical functions and computational routines
    - Contains Cython extensions (`functions_cy.pyx`, `bessel/besselk_dp.pyx`) for performance-critical code
  - **`modform/`** - Hilbert modular forms implementation
    - `hilbert_maass_element.py` - Individual form elements
    - `hilbert_maass_space.py` - Spaces of forms
    - `coefficients.py` - Coefficient computation
    - `eisenstein_series.py` - Eisenstein series
    - `hecke_operator.py` - Hecke operators
  - **`database/`** - Database models and persistence
  - **`search/`** - Search functionality

### Build System
- Uses setuptools with setuptools_scm for versioning
- Cython extensions require SageMath/PassageMath environment
- Automatic detection of Homebrew paths on macOS for library dependencies
- Extensions include Bessel function computations and optimized mathematical routines

### Dependencies
- **Runtime**: `hilbert_modular_group`, `comp_manager`, `jupyter`
- **Build**: `cython>=3.0.8`, `passagemath-environment`, `passagemath-flint`, `passagemath-modules`
- **Development**: `tox`, `ruff`, `pytest`, `pytest-cov`, `mongomock`, `freezegun`

### SageMath Integration
This package is designed to work within the SageMath ecosystem:
- Uses `sage` command for building, testing, and running
- Cython extensions link against SageMath libraries
- Doctests follow SageMath conventions
- Installation via `sage -pip install`

### Development Notes
- Ruff configuration allows some PEP8 violations for mathematical variable names
- Line length set to 100 characters
- Cython files generate corresponding .c files during build
- Docker support available for containerized development

### Development Guidelines
- Always follow SageMath development guidelines and add doctests to all new functions and make sure that they pass

### Testing Goals
- Doctest coverage should be 100%
- Always run tests and code checks in a virtual environment venv_test, where all dependencies are installed