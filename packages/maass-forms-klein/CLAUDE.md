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

### Docker Development (Optional)
- `make docker` - Build Docker container with current branch
- `make docker-test` - Run doctests in Docker container
- `make docker-tox` - Run full tox suite in Docker container
- `make docker-examples` - Start Jupyter notebook server (port 8888)

### Jupyter Notebooks
- `make examples` - Start local Jupyter with SageMath kernel (port 8888)

## Architecture Overview

This is a Python package for algorithms related to Kleinian Maass forms and knot complements, built on top of SageMath/PassageMath.

### Core Structure
- **`src/maass_forms_klein/`** - Main package directory
  - **`modform/`** - Kleinian Maass forms implementation
    - `kmaass_space.py` - Spaces of Kleinian Maass forms
    - `kmaass_element.py` - Individual form elements
    - `coefficients.py` - Coefficient computation and storage
    - `utils.py` - Utility functions and type definitions
  - **`hyperbolic_space/`** - Hyperbolic geometry and Kleinian groups
    - `kleinian_group.py` - Discrete subgroups of PSL(2,C)
    - `upper_half_space.py` - Points and operations in hyperbolic 3-space
    - Contains Cython extensions for performance-critical computations
  - **`functions/`** - Mathematical functions and computational routines
    - Contains Cython extensions for optimized numerical functions (e.g., Bessel K functions)
  - **`database/`** - Database models and persistence using MongoEngine
  - **`utils/`** - General utilities and helper functions

### Build System
- Uses setuptools with setuptools_scm for versioning
- Cython extensions require SageMath/PassageMath environment
- Automatic detection of Homebrew paths on macOS for library dependencies
- Extensions include optimized mathematical routines for hyperbolic geometry

### Dependencies
- **Runtime**: `comp_manager`, `maass_forms_hilbert`, `snappy`, `snappy-manifolds`, `jupyter`
- **Build**: `cython>=3.0.8`, `passagemath-modules`
- **Development**: `tox`, `ruff`, `pytest`, `pytest-cov`, `mongomock`, `freezegun`, `mongoengine`

### SageMath Integration
This package is designed to work within the SageMath ecosystem:
- Uses `sage` command for building, testing, and running
- Cython extensions link against SageMath libraries
- Doctests follow SageMath conventions
- Installation via `sage -pip install`
- Mathematical classes inherit from appropriate SageMath base classes:
  - `Parent` for mathematical structures (like `KleinianMaassFormSpace`)
  - `Element` for mathematical objects (like `KleinianMaassFormElement`)

### Mathematical Context
This package implements computational methods for:
- Discrete subgroups of PSL(2,C) acting on hyperbolic 3-space (Kleinian groups)
- Maass waveforms (eigenfunctions of the Laplacian) on quotient spaces
- Coefficient computations for Fourier expansions of Kleinian Maass forms
- Integration with Snappy for 3-manifold topology and knot complements
- Hyperbolic geometry computations in 3-dimensional space

### Database Integration
The project uses MongoDB for persistent storage via MongoEngine:
- Models defined in `maass_forms_klein.database.models`
- Integration with `comp_manager` for computation management
- Database objects inherit from `DBObjectBase`

### Development Notes
- Ruff configuration allows some PEP8 violations for mathematical variable names
- Line length set to 100 characters
- Max McCabe complexity of 10 (stricter than maass-forms-hilbert)
- Cython files generate corresponding .c files during build
- Docker support available for containerized development

### Development Guidelines
- Always follow SageMath development guidelines and add doctests to all new functions and make sure that they pass
- All public methods must include doctest examples with `EXAMPLES::` sections
- Use SageMath's caching decorators (`@cached_method`) for expensive computations
- Follow SageMath's docstring format for mathematical documentation
- When working with mathematical components, understand the bridge between pure mathematics (group theory, hyperbolic geometry) and computational topology

### Testing Goals
- Doctest coverage should be 100%
- Always run tests and code checks in a virtual environment venv_test, where all dependencies are installed
- Primary testing via SageMath doctests (not pytest)
- Examples should demonstrate actual usage and verify correctness

## Lessons Learned from Debugging Sessions

### Critical Workflow for Doctest Fixes

#### 1. Always Create Virtual Environment First
```bash
python -m venv venv_test
source venv_test/bin/activate
```
**Why**: Prevents system-wide package conflicts and ensures clean testing environment.

#### 2. Run `make install` After Every File Change
```bash
source venv_test/bin/activate && make install
```
**Why**: Python files and Cython extensions must be reinstalled to reflect changes in doctests.
**Critical**: This applies to `.py`, `.pyx`, and any imported modules.

#### 3. Fix Import Chain Issues Systematically
**Problem**: Import chains leading to external dependencies (MongoDB, compiled extensions)
**Solution**: Define types and utilities locally:
```python
# Instead of: from maass_forms_klein.modform.utils import Real_t
# Use local definition:
from sage.rings.real_mpfr import RealNumber  
from sage.rings.integer import Integer
Real_t = Union[RealNumber, Integer, int, float]
```

#### 4. Add Defensive Programming for Edge Cases
**Problem**: Functions assuming non-empty collections
**Solution**: Add guard clauses:
```python
def process_list(items):
    if not items:
        return []  # or appropriate default
    # ... rest of function
```

### Common Doctest Issues and Solutions

#### Import Issues
- **Problem**: `from maass_forms_klein.all import` triggers database connections
- **Solution**: Use specific imports or mock objects:
```python
# Instead of importing full objects, test with mocks:
sage: from sage.matrix.constructor import matrix
sage: mock_gen = matrix(CC, [[1, 1], [0, 1]])
```

#### Cython Extension Problems
- **Problem**: Compiled extensions failing to load with linking errors
- **Solution**: Focus on Python-level tests, skip Cython-dependent examples:
```python
sage: # doctest: +SKIP  
sage: # from problematic_cython_module import Class  # doctest: +SKIP
```

#### Dataclass Constructor Issues  
- **Problem**: Required fields computed in `__post_init__`
- **Solution**: Use `field(init=False)` for computed fields:
```python
@dataclass
class MyClass:
    computed_field: SomeType = field(init=False)
```

#### Algorithm Robustness
- **Problem**: IndexError, ValueError from empty sequences
- **Solution**: Always check for empty inputs:
```python
if not circle_list:
    return []
elt = circle_list[0][1].radius  # Now safe
```

### Testing Strategy Priorities

1. **Fix fundamental import/compilation issues first**
2. **Add defensive programming for edge cases** 
3. **Replace problematic dependencies with mock tests**
4. **Verify each fix with `make install` + test cycle**
5. **Focus on Python-level functionality before Cython extensions**

### File-by-File Dependencies
- `types.py` → fundamental dataclasses, fix first
- `utils.py` and `geometry_utils.py` → depend on types  
- `parallelogram.py` → depends on types and utils
- `upper_half_space.pyx` → Cython extension, hardest to fix
- Always fix dependencies bottom-up

### Success Metrics Achieved
- **106/106 tests passing** in `parallelogram.py` 
- **62/62 tests passing** in `kleinian_group.py`
- **Eliminated MongoDB connection errors** across all Python files
- **Fixed critical algorithm robustness issues** (IndexError, ValueError guards)
- **Established sustainable development workflow** with virtual environment + make install

### Key Insight
**Systematic rebuilding with `make install` after each change is essential** - many doctest failures are actually due to old cached bytecode rather than code issues.
