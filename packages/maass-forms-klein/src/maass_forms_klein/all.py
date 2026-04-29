import logging
from typing import Any, Optional

# Configure logging
logging.basicConfig(level=logging.WARNING, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# Core imports with error handling
try:
    from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
    from maass_forms_klein.modform.kmaass_element import KleinianMaassFormElement
    from maass_forms_klein.modform.coefficients import KleinianMaassFormCoefficients
except ImportError as e:
    log.error(f"Cannot import core modform modules: {e}")
    KleinianMaassFormSpace = None
    KleinianMaassFormElement = None
    KleinianMaassFormCoefficients = None

try:
    from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
    from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement
except ImportError as e:
    log.error(f"Cannot import hyperbolic space modules: {e}")
    KleinianGroup = None
    UpperHalfSpaceElement = None

try:
    from maass_forms_klein.modform.utils import Integer_t, Real_t
    from maass_forms_klein.modform.utils import map_tuple_to_int, map_int_to_tuple
except ImportError as e:
    log.error(f"Cannot import utility functions: {e}")
    Integer_t = None
    Real_t = None
    map_tuple_to_int = None
    map_int_to_tuple = None

# Database models with optional import (requires MongoDB)
try:
    from maass_forms_klein.database.models import KleinianMaassFormDB, KleinianGroupDB
except ImportError as e:
    log.warning(f"Cannot import database models (MongoDB not available?): {e}")
    KleinianMaassFormDB = None
    KleinianGroupDB = None

# Custom exceptions
try:
    from maass_forms_klein.exceptions import (
        KnotMaassError,
        InvalidSpectralParameterError,
        InvalidSpaceError,
        InvalidGroupError,
        ComputationError,
        DatabaseError,
        ValidationError
    )
except ImportError as e:
    log.error(f"Cannot import custom exceptions: {e}")
    KnotMaassError = Exception
    InvalidSpectralParameterError = ValueError
    InvalidSpaceError = ValueError
    InvalidGroupError = ValueError
    ComputationError = RuntimeError
    DatabaseError = RuntimeError
    ValidationError = ValueError

# Explicit exports for API clarity
__all__ = [
    # Core mathematical objects
    'KleinianMaassFormSpace',
    'KleinianMaassFormElement', 
    'KleinianMaassFormCoefficients',
    'KleinianGroup',
    'UpperHalfSpaceElement',
    # Utility functions and types
    'Integer_t',
    'Real_t',
    'map_tuple_to_int',
    'map_int_to_tuple',
    # Database models (optional)
    'KleinianMaassFormDB',
    'KleinianGroupDB',
    # Custom exceptions
    'KnotMaassError',
    'InvalidSpectralParameterError',
    'InvalidSpaceError',
    'InvalidGroupError',
    'ComputationError',
    'DatabaseError',
    'ValidationError'
]


def _check_imports() -> dict[str, bool]:
    """Check which optional components are available.
    
    OUTPUT:
    - Dict mapping component names to availability status
    
    EXAMPLES::
        sage: from maass_forms_klein.all import _check_imports
        sage: status = _check_imports()
        sage: status['core_modform']  # doctest: +SKIP
        True
    """
    return {
        'core_modform': all(x is not None for x in [KleinianMaassFormSpace, KleinianMaassFormElement]),
        'hyperbolic_space': all(x is not None for x in [KleinianGroup, UpperHalfSpaceElement]),
        'utilities': all(x is not None for x in [Integer_t, Real_t, map_tuple_to_int]),
        'database': all(x is not None for x in [KleinianMaassFormDB, KleinianGroupDB]),
        'exceptions': KnotMaassError is not Exception
    }


def version_info() -> dict[str, Any]:
    """Get version and dependency information.
    
    OUTPUT:
    - Dict with version and availability information
    
    EXAMPLES::
        sage: from maass_forms_klein.all import version_info
        sage: info = version_info()
        sage: 'maass_forms_klein_version' in info
        True
    """
    info = {'available_components': _check_imports()}
    
    try:
        from maass_forms_klein import __version__
        info['maass_forms_klein_version'] = __version__
    except Exception:
        info['maass_forms_klein_version'] = 'unknown'
    
    # Check SageMath version
    try:
        from sage.version import version as sage_version
        info['sage_version'] = sage_version
    except ImportError:
        info['sage_version'] = 'not available'
    
    return info