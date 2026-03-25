import logging
from typing import Any

try:
    from maass_form_core.functions.bessel.besselk_dp import besselk_dp
except ImportError:
    besselk_dp = None
from sage.functions.bessel import bessel_K
from sage.misc.functional import log
from sage.misc.misc_c import prod
from sage.rings.complex_mpfr import ComplexNumber
from sage.rings.integer import Integer
from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_gaussian
from sage.rings.rational import Rational
from sage.rings.real_mpfr import RealNumber as RealNumber_class, RR, RealNumber
from sage.misc.cachefunc import cached_function
from sage.functions.log import exp

this_log = logging.getLogger(__name__)
# User defined type for either Python or Sage type
Integer_t = Integer | int
Real_t = RealNumber_class | float
Complex_t = ComplexNumber | complex | NumberFieldElement_gaussian
Circle_t = tuple[tuple[Real_t, Real_t], Real_t]
Rectangle_t = tuple[tuple[Real_t], tuple[Real_t]]
class Rectangle:
    def __init__(self, x1, x2, y1, y2):
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2

@cached_function
def map_tuple_to_int(index_tuple: tuple, tuple_limits: tuple[tuple[Integer_t]],
                     tuple_len: int = None) -> int:
    r"""
    Map a tuple `(a0, a1)` with `min_i < ai < max_i` to an integer
     $\sum_i=0^(n-1) (max_0 - min_0 + 1)*(a0 - min_0) + (a1 - min_1)$

    NOTE: This function together with `map_int_to_tuple` provides an isomorphism between
            $[min_0,...,max_0] x [min_1,...,max_1] x ... x [min_{n-1},...,max_{n-1}]$
             and $[0,...,N]$ where $N = prod(max_i - min_i + 1)$.

    INPUT:

    - ``index_tuple`` -- tuple
    - ``tuple_limits`` -- tuple of tuples
    - ``tuple_len`` -- integer: number of tuples (default: None) if positive then the tuple_limits
                       are duplicated that number of times.
    EXAMPLES::

        sage: from maass_forms_klein.modform.utils import map_tuple_to_int
        sage: map_tuple_to_int((-1,), ((-1, 1),), 1)
        0
        sage: map_tuple_to_int((-1, -1),((-1, 1),), 2)
        0
        sage: map_tuple_to_int((-1, -1, -1), ((-1, 1),), 3)
        0
        sage: map_tuple_to_int((-1, -1),((-1, 1), (-1, 1)))
        0
        sage: map_tuple_to_int((0, 0),((-1, 1), (-1, 1)))
        4
        sage: map_tuple_to_int((-1, -3),((-1, 1), (-3, 1)))
        0
        sage: map_tuple_to_int((0, -3),((-1, 1), (-3, 1)))
        1
        sage: map_tuple_to_int((0, -3),((-1, 1), (-3, 1)))
        1
        sage: map_tuple_to_int((0, -4), ((0,5),(-5,5)))
        6
        sage: map_tuple_to_int((1, -4), ((0,5),(-5,5)))
        7

    TESTS::

        sage: map_tuple_to_int((-1, -1), ((-1, 1), (-1, -2)))
        Traceback (most recent call last):
        ...
        ValueError: tuple_limits ((-1, 1), (-1, -2)) do not give positive length intervals

        sage: map_tuple_to_int((-2, -1), ((-1, 1),(-1, 1)))
        Traceback (most recent call last):
        ...
        IndexError: Tuple element (-2, -1) is out of bounds!
        sage: map_tuple_to_int((-1, 2), ((-1,1),(-1,1)))
        Traceback (most recent call last):
        ...
        IndexError: Tuple element (-1, 2) is out of bounds!
    """
    if not isinstance(index_tuple, tuple) or not isinstance(tuple_limits, tuple):
        raise ValueError("Call with tuples!")
    if len(tuple_limits) == 1 and isinstance(tuple_len, (Integer, int)) and tuple_len > 1:
        tuple_limits = tuple_limits * tuple_len
    if len(index_tuple) != len(tuple_limits):
        raise ValueError(f"lengths differ: {len(index_tuple)} != {len(tuple_limits)}")
    if any(x[1] - x[0] + 1 <= 0 for x in tuple_limits):
        raise ValueError(f"tuple_limits {tuple_limits} do not give positive length intervals")
    if any(index_tuple[i] < min_tix or index_tuple[i] > max_tix
           for i, (min_tix, max_tix) in enumerate(tuple_limits)):
        raise IndexError(f"Tuple element {index_tuple} is out of bounds!")
    n = len(index_tuple)
    # Calculate the index of the tuple
    return int(sum((tuple_limits[i-1][1] - tuple_limits[i-1][0] + 1)**i*(index_tuple[i] - min_tix)
               for i, (min_tix, max_tix) in enumerate(tuple_limits)))


@cached_function
def map_int_to_tuple(index: Integer_t, tuple_limits: tuple[tuple[Integer_t]],
                     tuple_len: Integer_t = None) -> tuple:
    r"""
    Map integer to tuple (the inverse of map_tuple_to_int) by modding recursively
    modulo the lengths of the integer intervals.

    INPUT:

    - ``index`` -- integer
    - ``tuple_limits`` -- tuple of tuples of limits
    - ``tuple_len`` -- integer (number of tuples - duplicates the input tuple_limits)


    EXAMPLES::

    sage: from maass_forms_klein.modform.utils import map_int_to_tuple
    sage: map_int_to_tuple(0,((-1,1),), 1) == (-1, )
    True
    sage: map_int_to_tuple(0,((-1,1),), 2)
    (-1, -1)
    sage: map_int_to_tuple(0,((-1,1),), 3)
    (-1, -1, -1)
    sage: map_int_to_tuple(0,((-1,1),(-1,1)))
    (-1, -1)
    sage: map_int_to_tuple(0,((-1,1),(-3,1)))
    (-1, -3)
    sage: map_int_to_tuple(1,((-1,1),(-3,1)))
    (0, -3)

    TESTS::

    sage: map_int_to_tuple(0,((-1,1),(-1,-2)))
    Traceback (most recent call last):
    ...
    ValueError: tuple_limits ((-1, 1), (-1, -2)) do not give positive length intervals

    sage: map_int_to_tuple(9,((-1,1),(-1,1)))
    Traceback (most recent call last):
    ...
    IndexError: Index 9 is out of bounds!
    sage: map_int_to_tuple(-1,((-1,1),(-1,1)))
    Traceback (most recent call last):
    ...
    IndexError: Index -1 is out of bounds!
    """
    if not isinstance(tuple_limits, tuple):
        raise ValueError("Call with tuple!")
    if not isinstance(index, (int, Integer)):
        raise ValueError("Call with integer!")
    if len(tuple_limits) == 1 and isinstance(tuple_len, (Integer, int)) and tuple_len > 1:
        tuple_limits = tuple_limits * tuple_len
    if any(x[1] - x[0] + 1 <= 0 for x in tuple_limits):
        raise ValueError(f"tuple_limits {tuple_limits} do not give positive length intervals")
    if index < 0 or index >= prod(x[1] - x[0] + 1 for x in tuple_limits):
        raise IndexError(f"Index {index} is out of bounds!")
    ix_t = []
    for min_tix, max_tix in tuple_limits:
        range_tix = max_tix - min_tix + 1
        t = index % range_tix
        ix_t.append(t + min_tix)
        index = (index - t) / range_tix
    return tuple(ix_t)


@cached_function
def dual_ideal_element(index_tuple, ideal):
    b1, b2 = ideal.integral_basis()
    return b1 * index_tuple[0] + b2 * index_tuple[1]


@cached_function
def bessel_function(absv, y, s, pre_factor=1, sgn='+'):
    if absv == 0:
        if s == 1.0:
            return y if sgn == '+' else y * log(y)
        return y ** s if sgn == '+' else y ** (2 - s)
    elif absv < 0:
        raise ValueError("absv must be non-negative")
    # print("S in bes=",s,type(s))
    if isinstance(s, float) or (isinstance(s, RealNumber) and s.prec() == 53):
        return besselk_dp(s, absv * y, pref=pre_factor)
    else:
        if pre_factor:
            factor = s.parent().pi() * s / 2.0
            if isinstance(factor, RealNumber):
                factor = (s.parent().pi() * s / 2.0).exp()
            else:
                factor = exp(s.parent().pi() * s / 2.0)

        else:
            factor = 1
        return bessel_K(s, absv * y) * factor


@cached_function
def exp_trace(z):
    if isinstance(z, complex):
        z = 2 * RR.pi() * z.real
    else:
        z = 2 * z.parent().pi() * z.real()
    return z.exp()


def get_prec(x: Any) -> int:
    if isinstance(x, (RealNumber, ComplexNumber)):
        return int(x.prec())
    if isinstance(x, (Integer_t, Rational)):
        return int(0)
    return int(53)

def get_epsilon(x: Any) -> tuple[float,Any]:
    if isinstance(x, float):
        eps = 2 ** - 53
    elif hasattr(x, 'parent'):
        eps = x.parent().base_ring().epsilon()
    else:
        raise ValueError("Could not find base ring")
    return eps


