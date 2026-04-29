"""
Common type aliases for SageMath numerical types.

These type aliases provide unified type annotations across all Maass form packages,
supporting both native Python types and SageMath types.

EXAMPLES::

    sage: from maass_form_core.utils.types import Integer_t, Real_t, Complex_t
    sage: isinstance(1, Integer_t)
    True
    sage: isinstance(ZZ(1), Integer_t)
    True
    sage: isinstance(1.0, Real_t)
    True
    sage: isinstance(RR(1.0), Real_t)
    True
"""

from sage.rings.complex_mpfr import ComplexNumber
from sage.rings.integer import Integer
from sage.rings.real_mpfr import RealNumber as RealNumber_class

Integer_t = Integer | int
Real_t = RealNumber_class | float
Complex_t = ComplexNumber | complex
