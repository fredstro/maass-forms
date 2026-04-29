# Implementation Plan: maass-forms-core Shared Library

## Current State Assessment (Updated March 2026)

### What Exists

**maass-form-core package:** Scaffolding only — build system, pyproject.toml, Makefile, and 6 empty placeholder modules (`database`, `coefficients`, `spaces`, `functions`, `utils`, `testing`). No actual implementations.

**maass-forms-hilbert:** Full working codebase with database models, coefficient management, space/element abstractions, JSON serialization, index mapping utilities, and Hilbert-specific mathematics.

**maass-forms-klein:** Full working codebase with parallel structure to hilbert — high overlap in infrastructure patterns.

### Measured Duplication

| Component | Overlap % | Lines Duplicated (approx) |
|-----------|-----------|--------------------------|
| JSON serialization (`utils/json_converters.py`) | 95% | ~300 |
| Exception hierarchy | 90% | ~80 |
| Import/export pattern (`all.py`) | 85% | ~60 |
| Coefficient container (`modform/coefficients.py`) | 80% | ~400 |
| Space abstraction (`modform/*_space.py`) | 75% | ~200 |
| Element abstraction (`modform/*_element.py`) | 70% | ~350 |
| Database models (`database/models.py`) | 70% | ~250 |
| Utility functions (`modform/utils.py`) | 60% | ~200 |
| **Total** | **~70%** | **~1840** |

---

## Architecture Overview

```
maass-form-core/src/maass_form_core/
├── __init__.py
├── database/
│   ├── __init__.py
│   ├── models.py          # Base DB document, Point, save() pattern
│   └── queryset.py        # Base QuerySet with Integer coercion, spectral_range
├── coefficients/
│   ├── __init__.py
│   └── base.py            # CoefficientManagerBase (ModuleElement)
├── spaces/
│   ├── __init__.py
│   ├── base_space.py      # MaassFormSpace (Module) base class
│   └── base_element.py    # MaassFormElement (ModuleElement) base class
├── utils/
│   ├── __init__.py
│   ├── types.py           # Integer_t, Real_t, Complex_t
│   ├── json_converters.py # SageJSONEncoder, ring/element/matrix/nf serialization
│   ├── indexing.py        # map_tuple_to_int, map_int_to_tuple, cartesian utilities
│   └── validation.py      # Common validation helpers
├── functions/
│   ├── __init__.py
│   └── basis.py           # number_field_basis_matrix, ideal_basis_matrix
├── testing/
│   ├── __init__.py
│   └── fixtures.py        # Shared test base classes, assertion helpers
└── exceptions.py          # Base exception hierarchy
```

---

## Phase 1: Foundation — Utils, Exceptions, JSON (Weeks 1-2)

**Goal:** Extract the highest-duplication, lowest-risk components first. These are pure functions with no inheritance complexities.

### Week 1: Type Definitions, JSON Converters, Exceptions

#### 1.1 Create `utils/types.py`
Extract from `hilbert/modform/utils.py` and `klein/modform/utils.py`:
```python
from sage.rings.integer import Integer
from sage.rings.real_mpfr import RealNumber as RealNumber_class
from sage.rings.complex_mpfr import ComplexNumber

Integer_t = Integer | int
Real_t = RealNumber_class | float
Complex_t = ComplexNumber | complex
```

#### 1.2 Create `utils/json_converters.py`
Extract from `klein/utils/json_converters.py` (95% identical to hilbert's `modform/utils.py` JSON functions):
- `SageJSONEncoder` — custom JSON encoder for SageMath types
- `ring_to_json()` / `ring_from_json()`
- `ring_element_to_json()` / `ring_element_from_json()`
- `matrix_to_json()` / `matrix_from_json()`
- `number_field_to_json()` / `number_field_from_json()`
- `complex_number_to_json()` / `complex_number_from_json()`
- `complex_tuple_to_json()` / `complex_tuple_from_json()`
- `is_json_number()`
- `dict_to_json()` / `dict_from_json()`

Source files:
- `packages/maass-forms-klein/src/maass_forms_klein/utils/json_converters.py`
- `packages/maass-forms-hilbert/src/maass_forms_hilbert/modform/utils.py` (lines 604-886)

#### 1.3 Create `exceptions.py`
Extract from `klein/exceptions.py` (hilbert lacks a formal exceptions module):
```python
class MaassFormError(Exception):
    """Base exception for all maass-form packages."""

class InvalidSpectralParameterError(MaassFormError):
    def __init__(self, message, parameter_value=None):
        ...

class InvalidSpaceError(MaassFormError):
    def __init__(self, message, space_config=None):
        ...

class ComputationError(MaassFormError):
    def __init__(self, message, computation_details=None):
        ...

class DatabaseError(MaassFormError):
    def __init__(self, message, operation=None):
        ...

class ValidationError(MaassFormError):
    def __init__(self, message, field_name=None, value=None):
        ...
```

#### 1.4 Tests for Week 1
- Unit tests for all JSON converters (roundtrip Integer, Real, Complex, Matrix, NumberField)
- Unit tests for exception hierarchy
- Type alias verification

### Week 2: Index Mapping, Basis Utilities, Validation

#### 2.1 Create `utils/indexing.py`
Extract from `hilbert/modform/utils.py` (lines 76-326) and `klein/modform/utils.py`:
- `cartesian_product_from_M(M)` — generate all index tuples from bounds
- `is_tuple_zero(t)` — check if a tuple is all zeros
- `length_from_M(M)` — calculate cartesian product size
- `get_Q_from_bounds(M)` — compute bounding box
- `map_tuple_to_int(t, M)` — bijection from tuples to integers
- `map_int_to_tuple(n, M)` — inverse bijection
- `integer_to_bounds_tuple(n)` — convert integer to bounds format
- `coefficient_dict_to_matrix(d, M)` — convert dict to matrix
- `list_of_lists_to_limits(ll)` — extract range limits

#### 2.2 Create `functions/basis.py`
Extract from `hilbert/modform/utils.py` (lines 330-601):
- `number_field_basis_matrix(K)` — compute basis matrix for a number field
- `ideal_basis_matrix(I)` — compute basis matrix for a fractional ideal
- `dual_ideal(I)` — compute dual/codifferent ideal
- `dual_ideal_basis_matrix(I)` — basis matrix of dual ideal
- `ideal_coordinates(x, I)` — express element in terms of ideal basis
- `ideal_generator(I)` — find totally positive generator

These are generic algebraic number theory utilities used by both packages.

#### 2.3 Create `utils/validation.py`
Common validation patterns extracted from both packages:
- Spectral parameter validation (type checking, range checking)
- Space configuration validation
- Bounds/precision parameter validation

#### 2.4 Tests for Week 2
- Index mapping roundtrip tests (tuple→int→tuple)
- Basis matrix correctness with known number fields
- Validation edge cases

#### 2.5 Integration Milestone
At end of Week 2:
- [ ] `from maass_form_core.utils.types import Integer_t, Real_t, Complex_t` works
- [ ] `from maass_form_core.utils.json_converters import SageJSONEncoder` works
- [ ] `from maass_form_core.utils.indexing import map_tuple_to_int` works
- [ ] `from maass_form_core.exceptions import MaassFormError` works
- [ ] All unit tests pass
- [ ] `ruff check` and `ruff format` pass

---

## Phase 2: Database and Coefficient Base Classes (Weeks 3-5)

**Goal:** Extract the database patterns and coefficient management base classes. These form the persistence layer shared by both domain packages.

### Week 3: Database Base Classes

#### 3.1 Create `database/models.py`
Extract from `hilbert/database/models.py` and `klein/database/models.py`:

```python
"""Base database document classes for Maass form packages."""
import mongoengine as me

class PointDB(me.EmbeddedDocument):
    """Generic 2D point for complex number storage."""
    x = me.FloatField(required=True)
    y = me.FloatField(required=True)
    meta = {'allow_inheritance': True}

    def __str__(self):
        return f"({self.x}, {self.y})"

class MathematicalObjectDB(me.Document):
    """Abstract base for mathematical objects stored in MongoDB.

    Provides:
    - spectral_parameter dict field
    - parent dict field (space configuration)
    - coefficients dict field
    - set_coefficients dict field
    - status field with standard choices
    - comments field
    - Generic save() with pre-processing hook
    - near_or_create() class method pattern
    """
    meta = {'abstract': True}

    spectral_parameter = me.DictField()
    parent = me.DictField()
    coefficients = me.DictField()
    set_coefficients = me.DictField()
    status = me.StringField(choices=['new', 'computed', 'verified', 'error'])
    comments = me.StringField(default='')

    def pre_save_processing(self):
        """Hook for domain-specific processing before save. Override in subclass."""
        pass

    def save(self, *args, **kwargs):
        self.pre_save_processing()
        return super().save(*args, **kwargs)

    @classmethod
    def near_or_create(cls, search_criteria, tolerance=1e-6, **defaults):
        """Find an existing object near the given criteria, or create a new one."""
        raise NotImplementedError("Subclasses must implement near_or_create")
```

#### 3.2 Create `database/queryset.py`
Extract shared QuerySet patterns:

```python
"""Base QuerySet class with common filtering patterns."""
import mongoengine as me
from sage.rings.integer import Integer

class MaassFormQuerySet(me.QuerySet):
    """Base QuerySet with SageMath Integer coercion and common filters."""

    def __getitem__(self, key):
        if isinstance(key, Integer):
            key = int(key)
        return super().__getitem__(key)

    def spectral_range(self, r_min, r_max):
        """Filter by spectral parameter range. Override for domain-specific logic."""
        raise NotImplementedError("Subclasses must implement spectral_range")
```

#### 3.3 Tests for Week 3
- Database model CRUD with mongomock
- QuerySet Integer coercion
- save() pre-processing hook invocation
- Abstract method enforcement

### Week 4-5: Coefficient Base Class

#### 4.1 Create `coefficients/base.py`
Extract from `hilbert/modform/coefficients.py` (~80% shared):

```python
"""Base class for Maass form coefficient containers."""
from sage.structure.element import ModuleElement

class CoefficientManagerBase(ModuleElement):
    """Abstract base for coefficient storage across all Maass form domains.

    Provides:
    - Initialization with matrix/bounds/spectral parameter/space
    - to_json() / from_json() serialization
    - __hash__(), __eq__()
    - M(), Y(), Q(), space(), spectral_parameter() accessors
    - coefficient_matrix(), set_coefficients() accessors
    - __getitem__() with index/slice support
    - __iter__(), keys(), values()
    - _add_(), _sub_(), _neg_(), _lmul_() arithmetic
    """

    def __init__(self, parent, coefficients, M, spectral_parameter,
                 space, Y=None, Q=None, set_coefficients=None, **kwargs):
        ...

    # --- Serialization ---
    def to_json(self) -> dict: ...
    @classmethod
    def from_json(cls, data, **kwargs): ...

    # --- Accessors ---
    def M(self): ...
    def Y(self): ...
    def Q(self): ...
    def space(self): ...
    def spectral_parameter(self): ...
    def coefficient_matrix(self): ...
    def set_coefficients(self): ...

    # --- Indexing ---
    def __getitem__(self, key): ...
    def __iter__(self): ...
    def keys(self): ...
    def values(self): ...

    # --- Comparison ---
    def __hash__(self): ...
    def __eq__(self, other): ...

    # --- Arithmetic ---
    def _add_(self, other): ...
    def _sub_(self, other): ...
    def _neg_(self): ...
    def _lmul_(self, scalar): ...
```

Domain packages override:
- Hilbert: adds `coordinate_ideals()`, `ideal_coordinates()`, `norms()`, error estimation
- Klein: adds `coordinate_values()` for translation lattice vectors

#### 4.2 Tests for Week 4-5
- Coefficient creation and accessors
- JSON roundtrip serialization
- Arithmetic operations (add, sub, scalar mul)
- Indexing and iteration
- Hash consistency

#### 4.3 Integration Milestone
At end of Week 5:
- [ ] `from maass_form_core.database.models import MathematicalObjectDB, PointDB` works
- [ ] `from maass_form_core.coefficients.base import CoefficientManagerBase` works
- [ ] Database models can be subclassed and saved (with mongomock)
- [ ] Coefficient base class supports all arithmetic operations
- [ ] All tests pass with >90% coverage

---

## Phase 3: Space and Element Base Classes (Weeks 6-8)

**Goal:** Extract the SageMath Parent/Element pattern into reusable base classes.

### Week 6: Base Space Class

#### 6.1 Create `spaces/base_space.py`
Extract from `hilbert/modform/hilbert_maass_space.py` (~75% shared):

```python
"""Base class for mathematical spaces of Maass forms."""
from sage.structure.parent import Parent

class MaassFormSpace(Parent):
    """Abstract base space for Maass forms.

    Provides:
    - Element class registration
    - to_json() / from_json() serialization
    - __hash__(), __eq__(), __repr__()
    - is_cuspidal() property
    - an_element() factory
    - _element_constructor_() coercion
    """
    Element = None  # Set by subclass

    def __init__(self, group, cuspidal=True, **kwargs):
        self._group = group
        self._cuspidal = cuspidal
        super().__init__(**kwargs)

    def group(self):
        return self._group

    def is_cuspidal(self):
        return self._cuspidal

    def to_json(self) -> dict: ...
    @classmethod
    def from_json(cls, data): ...

    def __hash__(self): ...
    def __eq__(self, other): ...
    def __repr__(self): ...

    def an_element(self): ...
    def _element_constructor_(self, spectral_parameter, **kwargs): ...
```

### Week 7: Base Element Class

#### 7.1 Create `spaces/base_element.py`
Extract from `hilbert/modform/hilbert_maass_element.py` (~70% shared):

```python
"""Base class for individual Maass form elements."""
from sage.structure.element import ModuleElement

class MaassFormElement(ModuleElement):
    """Abstract base for individual Maass form instances.

    Provides:
    - Initialization with parent space, spectral parameter, coefficients
    - to_json() / from_json() serialization
    - __reduce__() pickling
    - __hash__(), __repr__(), __eq__(), __ne__()
    - spectral_parameter() accessor
    - coefficients() accessor
    - is_cuspidal() delegating to parent
    - Arithmetic: __mul__, _rmul_, __add__, _add_, __sub__, _sub_
    - __copy__()
    """

    def __init__(self, parent, spectral_parameter, coefficients=None, **kwargs):
        self._spectral_parameter = spectral_parameter
        self._coefficients = coefficients
        super().__init__(parent)

    def spectral_parameter(self): ...
    def coefficients(self): ...
    def is_cuspidal(self):
        return self.parent().is_cuspidal()

    def to_json(self) -> dict: ...
    @classmethod
    def from_json(cls, data): ...

    def __reduce__(self): ...
    def __hash__(self): ...
    def __repr__(self): ...
    def __eq__(self, other): ...

    def __mul__(self, other): ...
    def _rmul_(self, scalar): ...
    def __add__(self, other): ...
    def _add_(self, other): ...
    def __sub__(self, other): ...
    def _sub_(self, other): ...
    def __copy__(self): ...
```

### Week 8: Testing Framework and Integration

#### 8.1 Create `testing/fixtures.py`
Shared test infrastructure:

```python
"""Shared test base classes and assertion helpers."""

class MaassFormTestBase:
    """Base test class with common assertions."""

    def assert_json_roundtrip(self, obj, from_json_cls=None):
        """Assert that to_json -> from_json reproduces the object."""
        ...

    def assert_spectral_parameter_valid(self, element):
        """Assert spectral parameter meets domain constraints."""
        ...

    def assert_coefficient_arithmetic(self, c1, c2):
        """Assert coefficient arithmetic consistency."""
        ...

    def assert_space_element_consistency(self, space, element):
        """Assert element belongs to space and inherits properties."""
        ...
```

#### 8.2 Integration Milestone
At end of Week 8:
- [ ] `from maass_form_core.spaces.base_space import MaassFormSpace` works
- [ ] `from maass_form_core.spaces.base_element import MaassFormElement` works
- [ ] `from maass_form_core.testing.fixtures import MaassFormTestBase` works
- [ ] Base classes can be subclassed with custom group/field types
- [ ] SageMath Parent/Element coercion framework functions correctly
- [ ] All tests pass with >95% coverage in core
- [ ] Full API documentation for all public classes

---

## Phase 4: Migrate maass-forms-hilbert (Weeks 9-11)

**Goal:** Refactor hilbert to inherit from core base classes while maintaining 100% backward compatibility.

### Week 9: Utility Migration

#### 9.1 Replace Duplicated Utilities
In `maass_forms_hilbert/modform/utils.py`:
- Remove type definitions → import from `maass_form_core.utils.types`
- Remove JSON converters → import from `maass_form_core.utils.json_converters`
- Remove index mapping functions → import from `maass_form_core.utils.indexing`
- Remove basis utilities → import from `maass_form_core.functions.basis`
- **Keep**: Hecke relations, symmetric relations, bilinear forms (Hilbert-specific)

#### 9.2 Add Re-exports for Backward Compatibility
Ensure existing imports like `from maass_forms_hilbert.modform.utils import map_tuple_to_int` still work by re-exporting from the module.

#### 9.3 Update `pyproject.toml`
Add dependency:
```toml
dependencies = [
    "maass_form_core>=1.0.0",
    # ... existing deps
]
```

### Week 10: Database and Coefficient Migration

#### 10.1 Migrate Database Models
In `maass_forms_hilbert/database/models.py`:
- `Point` → inherit from `maass_form_core.database.models.PointDB`
- `HilbertMaassformQuerySet` → inherit from `maass_form_core.database.queryset.MaassFormQuerySet`
- `HilbertMaassFormDB` → inherit from `maass_form_core.database.models.MathematicalObjectDB`
- **Keep**: Hilbert-specific fields (`r_values`, `y_values`), Hilbert-specific query methods

#### 10.2 Migrate Coefficient Class
In `maass_forms_hilbert/modform/coefficients.py`:
- `HilbertMaassCoefficients` → inherit from `maass_form_core.coefficients.base.CoefficientManagerBase`
- Move generic `to_json`, `from_json`, arithmetic, accessors to base class calls via `super()`
- **Keep**: `coordinate_ideals()`, `ideal_coordinates()`, `norms()`, error estimation

### Week 11: Space and Element Migration

#### 11.1 Migrate Space Class
In `maass_forms_hilbert/modform/hilbert_maass_space.py`:
- `HilbertMaassFormSpace` → inherit from `maass_form_core.spaces.base_space.MaassFormSpace`
- Move generic `to_json`, `from_json`, `__hash__`, `__eq__`, `__repr__` to `super()` calls
- **Keep**: HilbertModularGroup init, number field ops, dual ideals, pullback, check_interval, functional

#### 11.2 Migrate Element Class
In `maass_forms_hilbert/modform/hilbert_maass_element.py`:
- `HilbertMaassForm_Element` → inherit from `maass_form_core.spaces.base_element.MaassFormElement`
- Move generic serialization, arithmetic, comparison to `super()` calls
- **Keep**: `__call__` evaluation, galois_conjugate, action_by_unit, compute_coefficients, plot/animation

#### 11.3 Verification
- [ ] All existing doctests pass unchanged
- [ ] `from maass_forms_hilbert.all import *` works identically
- [ ] Existing user scripts produce identical results
- [ ] No performance regression in coefficient computation (benchmark)

---

## Phase 5: Migrate maass-forms-klein (Weeks 12-14)

**Goal:** Same migration pattern as hilbert. Klein should be simpler since the pattern is established.

### Week 12: Utility and Database Migration

#### 12.1 Replace Duplicated Utilities
- Remove `maass_forms_klein/utils/json_converters.py` entirely → use `maass_form_core.utils.json_converters`
- Remove duplicated type definitions and index mappings
- **Keep**: Bessel K-function utilities, geometry-specific helpers

#### 12.2 Migrate Database Models
- `KleinianMaassFormQuerySet` → inherit from `MaassFormQuerySet`
- `KleinianMaassFormDB` → inherit from `MathematicalObjectDB`
- **Keep**: `ParallelogramDB`, `Word`, `KleinianGroupDB` (geometry-specific)

### Week 13: Space, Element, Coefficient Migration

#### 13.1 Migrate Core Classes
- `KleinianMaassFormSpace` → inherit from `MaassFormSpace`
- `KleinianMaassFormElement` → inherit from `MaassFormElement`
- `KleinianMaassFormCoefficients` → inherit from `CoefficientManagerBase`
- **Keep**: All hyperbolic geometry, manifold topology, 3D evaluation

#### 13.2 Add Exception Integration
- `KnotMaassError` → inherit from `maass_form_core.exceptions.MaassFormError`
- Domain exceptions → inherit from core exception hierarchy

### Week 14: Integration Testing

#### 14.1 Cross-Package Tests
- Test that hilbert and klein can coexist with shared core
- Test that core version upgrades don't break either package
- Test JSON interoperability (serialize from one, deserialize in other where applicable)

#### 14.2 Verification
- [ ] All existing doctests pass unchanged
- [ ] `from maass_forms_klein.all import *` works identically
- [ ] No performance regression
- [ ] Database backward compatibility confirmed

---

## Phase 6: Verification and Documentation (Week 15)

### 15.1 Verification Checklist

#### Core Package
- [ ] `pip install maass-forms-core` works
- [ ] >95% test coverage
- [ ] All public APIs documented with SageMath-style docstrings
- [ ] `ruff check` and `ruff format` pass

#### Backward Compatibility
- [ ] All existing `maass-forms-hilbert` user code works unchanged
- [ ] All existing `maass-forms-klein` user code works unchanged
- [ ] All import paths preserved (with re-exports where needed)
- [ ] Existing MongoDB data accessible without migration
- [ ] <5% performance degradation in critical paths

#### Quality Metrics
- [ ] Code duplication reduced by >70% (measured by lines)
- [ ] maass-forms-core: >95% coverage
- [ ] maass-forms-hilbert: >90% coverage maintained
- [ ] maass-forms-klein: >90% coverage maintained

### 15.2 Documentation
- API reference for all core public classes
- Migration guide for downstream users
- Architecture overview with module diagram
- Examples: creating custom spaces, implementing coefficient managers, database integration

---

## What Stays Domain-Specific

### maass-forms-hilbert retains:
- HilbertModularGroup integration
- Number field and fractional ideal operations
- Hecke operators, eigenvalue detection, algebraic relations
- Multi-dimensional spectral parameter search
- Dual ideal computations and pullback
- Galois conjugation and unit group actions
- Hilbert Eisenstein series
- Plot/animation with product of upper half-planes geometry

### maass-forms-klein retains:
- KleinianGroup (PSL(2,ℂ) discrete subgroups)
- Hyperbolic 3-space geometry (upper half-space model)
- Parallelogram fundamental domains, word reduction
- Bessel K-function evaluation
- SnapPy manifold integration
- 3D hyperbolic evaluation and visualization

---

## Risk Mitigation

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Breaking existing user code | Medium | High | Re-export all moved symbols from original locations; deprecation warnings in v2.1, removal in v3.0 |
| Performance regression from abstraction | Low | Medium | Benchmark critical paths before/after; provide escape hatches |
| SageMath Parent/Element coercion issues | Medium | High | Test coercion framework thoroughly; keep exact same class hierarchy depth |
| Database schema incompatibility | Low | High | Abstract base uses `meta = {'abstract': True}`, no schema changes |
| Circular import between core and domain | Medium | Medium | Core must never import from domain packages; strict dependency direction |

---

## Implementation Order Priority

For maximum value with minimum risk, implement in this order within each phase:

1. **JSON converters** (95% identical, pure functions, zero risk)
2. **Type definitions** (trivial, immediate value)
3. **Index mapping** (pure functions, well-tested)
4. **Exceptions** (simple hierarchy, no dependencies)
5. **Basis utilities** (pure functions, used by both)
6. **Database base classes** (enables shared persistence)
7. **Coefficient base class** (complex but high-value)
8. **Space base class** (SageMath integration required)
9. **Element base class** (depends on Space)
10. **Testing fixtures** (supports all of the above)

---

## Success Metrics

### Quantitative
- **Code duplication**: >70% reduction (~1840 → <550 duplicated lines)
- **Test coverage**: >95% in core, >90% in domain packages
- **Performance**: <5% degradation in critical computation paths
- **New domain bootstrap**: <500 lines + <2 weeks (down from 1500+ lines / 4+ weeks)

### Qualitative
- Bug fixes in infrastructure propagate to all packages automatically
- New mathematical domains can reuse infrastructure without copy-paste
- Clear separation between what's shared and what's domain-specific
- Architecture supports future packages (e.g., Siegel Maass forms, Bianchi forms)