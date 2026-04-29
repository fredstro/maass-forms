# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a mathematical software ecosystem for Maass forms computation consisting of two domain-specific packages and planned shared infrastructure. The project focuses on automorphic forms over different mathematical domains.

## Repository Structure (Monorepo)

- **`packages/maass-form-core/`** - Shared infrastructure (database, coefficients, utilities)
- **`packages/maass-forms-hilbert/`** - Hilbert Maass forms over number fields (algebraic number theory)
- **`packages/maass-forms-klein/`** - Kleinian Maass forms on hyperbolic 3-manifolds (geometric topology)
- **`docs/`** - Architecture documentation and migration plans

## Development Commands

### Top-Level Commands
- `make install` - Install all packages (core first, then domain packages)
- `make test` - Run tests across all packages
- `make lint` - Lint all packages with ruff
- `make format` - Format all packages with ruff
- `make clean` - Clean build artifacts across all packages

### Package-Level Commands (run in packages/maass-forms-hilbert/ or packages/maass-forms-klein/)

#### Installation and Building
- `make install` - Clean build, create source distribution, and install with sage
- `make build` - Build Cython extensions in-place using `sage -python setup.py build_ext --inplace`
- `make sdist` - Create source distribution
- `make clean` - Remove build artifacts, compiled extensions, and temporary files

#### Testing
- `make test` - Run Sage doctests on src/ directory (`sage -t src/*`)
- `pytest` - Run pytest tests (requires `sage -i pytest pytest_xdist` first)
- `make tox` - Install tox and run full test suite including doctests, coverage, and linting

#### Code Quality
- `ruff check .` - Run linting checks
- `ruff format .` - Format code
- `tox -e ruff` - Run both ruff check and format via tox
- `tox -e codespell` - Check for misspelled words
- `tox -e relint` - Check for forbidden patterns

#### Development Environment
- `make examples` - Start local Jupyter with SageMath kernel (port 8888)
- `make docker` - Build Docker container for isolated development
- `make docker-test` - Run tests in Docker container
- `make docker-examples` - Start Jupyter notebook server in Docker (port 8888)

## Architecture and Dependencies

### Core Framework
- **SageMath/PassageMath**: Primary mathematical computing framework
- **Python 3.11+**: Programming language
- **Cython**: Performance-critical mathematical computations
- **MongoDB**: Database persistence via MongoEngine

### Build Requirements
- `sage` command must be available in PATH
- Cython extensions require SageMath environment
- Automatic detection of Homebrew paths on macOS for library dependencies

### Mathematical Domains

#### Hilbert Maass Forms (maass-forms-hilbert/)
- **Groups**: Hilbert modular groups over number fields
- **Geometry**: Products of upper half-planes
- **Spectral Parameters**: Multi-dimensional
- **Dependencies**: `hilbert_modular_group`, `comp_manager`

#### Kleinian Maass Forms (maass-forms-klein/)
- **Groups**: Kleinian groups (discrete subgroups of PSL(2,ℂ))
- **Geometry**: Hyperbolic 3-space quotients
- **Spectral Parameters**: 1-dimensional complex
- **Dependencies**: `snappy`, `snappy-manifolds`, `maass_forms_hilbert`

## Development Guidelines

### SageMath Integration
- Use `sage` command for building, testing, and running
- Follow SageMath docstring format with `EXAMPLES::` sections
- Doctests must pass with 100% coverage goal
- Mathematical classes inherit from SageMath base classes (`Parent`, `Element`)
- Use SageMath's caching decorators (`@cached_method`) for expensive computations

### Code Quality Standards
- Line length: 100 characters
- Ruff configuration allows mathematical variable naming conventions
- Max McCabe complexity: 10 for maass-forms-klein, slightly higher for maass-forms-hilbert
- Always run tests in virtual environment (`venv_test`) with all dependencies

### Testing Requirements
- Primary testing via SageMath doctests (not pytest)
- All public methods require doctest examples
- Examples should demonstrate actual usage and verify correctness
- Test coverage tracking via `tox -e coverage`

### Mathematical Documentation
- Comprehensive docstrings following SageMath conventions
- Mathematical context and background in docstrings
- INPUT/OUTPUT/EXAMPLES sections required
- Bridge pure mathematics with computational implementation

## Architectural Migration Plan

### Current State (v1.x)
Both packages operate independently with ~70% code duplication in database patterns, coefficient management, and mathematical space abstractions.

### Future Architecture (v2.0+)
**maass-forms-core/** shared library will contain:
- Common database patterns and abstractions
- Coefficient management and storage
- Mathematical space base classes  
- Shared utilities, validation, logging
- Common test fixtures and assertions

This will eliminate code duplication while maintaining mathematical domain separation.

## Quality Assurance

### Pre-commit Requirements
1. All Sage doctests must pass: `make test`
2. Linting must pass: `ruff check .`
3. Code formatting: `ruff format .`
4. Full test suite: `make tox`

### Development Environment
- Always work in `venv_test` virtual environment with all dependencies installed
- Use Docker for isolated testing when needed
- Jupyter notebooks available for interactive development and examples

## Project Philosophy

This ecosystem maintains strict separation between mathematical infrastructure (shared) and domain-specific algorithms (separate packages), enabling both mathematical rigor and software engineering best practices.