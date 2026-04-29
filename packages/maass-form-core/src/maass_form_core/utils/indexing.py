"""
Index mapping utilities for Cartesian product spaces.

Provides bijections between multi-dimensional integer tuples and linear indices,
used for coefficient storage and retrieval across all Maass form packages.

EXAMPLES::

    sage: from maass_form_core.utils.indexing import map_tuple_to_int, map_int_to_tuple
    sage: map_tuple_to_int((0, 0), ((-1, 1), (-1, 1)))
    4
    sage: map_int_to_tuple(4, ((-1, 1), (-1, 1)))
    (0, 0)
"""

from sage.categories.sets_cat import cartesian_product
from sage.misc.cachefunc import cached_function
from sage.misc.misc_c import prod
from sage.rings.integer import Integer


@cached_function
def cartesian_product_from_M(M: tuple) -> list:
    r"""
    Return a list of all vectors in the Cartesian product of integer ranges defined by M.

    INPUT:

    - ``M`` -- tuple of tuples; each inner tuple (min, max) defines an integer range

    OUTPUT:

    - List of tuples in the Cartesian product of the ranges

    EXAMPLES::

        sage: from maass_form_core.utils.indexing import cartesian_product_from_M
        sage: cartesian_product_from_M(((-1, 1), (-2, 2)))
        [(-1, -2), (-1, -1), (-1, 0), (-1, 1), (-1, 2),
         (0, -2), (0, -1), (0, 0), (0, 1), (0, 2),
         (1, -2), (1, -1), (1, 0), (1, 1), (1, 2)]
    """
    return list(cartesian_product([range(m0[0], m0[1] + 1) for m0 in M]))


def is_tuple_zero(t: tuple) -> bool:
    r"""
    Check if a tuple consists entirely of zeros.

    INPUT:

    - ``t`` -- tuple of integers

    OUTPUT:

    - Boolean; True if all elements are zero, False otherwise

    EXAMPLES::

        sage: from maass_form_core.utils.indexing import is_tuple_zero
        sage: is_tuple_zero((0, 0, 0))
        True
        sage: is_tuple_zero((0, 1, 0))
        False
    """
    return all(t0 == 0 for t0 in t)


def length_from_M(M: tuple) -> int:
    r"""
    Calculate the total number of points in the Cartesian product of integer ranges.

    INPUT:

    - ``M`` -- tuple of tuples of integers; each inner tuple (min, max) defines a range

    OUTPUT:

    - Integer representing the total number of points

    EXAMPLES::

        sage: from maass_form_core.utils.indexing import length_from_M
        sage: length_from_M(((-1, 1), (-2, 2)))
        15
        sage: length_from_M(((-5, 5), (-5, 5)))
        121
    """
    return prod([m0[1] - m0[0] + 1 for m0 in M])


def get_Q_from_bounds(M: tuple) -> tuple:
    r"""
    Find a bounding box for the integer coordinates based on the given bounds.

    Creates a bounding box for the cube [-M1,M1] x [-M2,M2] x...
    for the integer coordinates corresponding to the box [-b,b]^n in the lattice.

    INPUT:

    - ``M`` -- tuple of tuples of integers; each inner tuple (min, max) defines bounds

    OUTPUT:

    - Tuple of integers representing the bounding box

    EXAMPLES::

        sage: from maass_form_core.utils.indexing import get_Q_from_bounds
        sage: get_Q_from_bounds(((-1, 1), (-2, 2)))
        (4, 4)
        sage: get_Q_from_bounds(((-5, 5), (-10, 10)))
        (12, 12)
    """
    C = max(max(abs(b0), abs(b1)) for b0, b1 in M) + 2
    return (C,) * len(M)


@cached_function
def map_tuple_to_int(index_tuple: tuple, tuple_limits: tuple, tuple_len: int = None) -> int:
    r"""
    Map a tuple (a0,a1,...,a[n-1]) with min_i <= ai <= max_i to an integer.

    Computes $\sum_{i=0}^{n-1} (\prod_{j<i} (max_j - min_j + 1)) \cdot (a_i - min_i)$.

    Together with :func:`map_int_to_tuple`, this provides an isomorphism between
    $[min_0, max_0] \times \cdots \times [min_{n-1}, max_{n-1}]$ and $[0, N-1]$
    where $N = \prod (max_i - min_i + 1)$.

    INPUT:

    - ``index_tuple`` -- tuple of integers
    - ``tuple_limits`` -- tuple of (min, max) tuples
    - ``tuple_len`` -- integer (default: None); if positive, duplicate tuple_limits

    OUTPUT:

    - int; the linear index

    EXAMPLES::

        sage: from maass_form_core.utils.indexing import map_tuple_to_int
        sage: map_tuple_to_int((0, 1), ((-5, 5), (-5, 5)), 2)
        71
        sage: map_tuple_to_int((-1,), ((-1, 1),), 1)
        0
        sage: map_tuple_to_int((-1, -1), ((-1, 1),), 2)
        0
        sage: map_tuple_to_int((-1, -1, -1), ((-1, 1),), 3)
        0
        sage: map_tuple_to_int((-1, -1), ((-1, 1), (-1, 1)))
        0
        sage: map_tuple_to_int((0, 0), ((-1, 1), (-1, 1)))
        4
        sage: map_tuple_to_int((-1, -3), ((-1, 1), (-3, 1)))
        0
        sage: map_tuple_to_int((0, -3), ((-1, 1), (-3, 1)))
        1
        sage: map_tuple_to_int((0, -4), ((0, 5), (-5, 5)))
        6
        sage: map_tuple_to_int((1, -4), ((0, 5), (-5, 5)))
        7

    TESTS::

        sage: map_tuple_to_int((-1, -1), ((-1, 1), (-1, -2)))
        Traceback (most recent call last):
        ...
        ValueError: tuple_limits ((-1, 1), (-1, -2)) do not give positive length intervals

        sage: map_tuple_to_int((-2, -1), ((-1, 1), (-1, 1)))
        Traceback (most recent call last):
        ...
        IndexError: Tuple element (-2, -1) is out of bounds!
        sage: map_tuple_to_int((-1, 2), ((-1, 1), (-1, 1)))
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
    if any(
        index_tuple[i] < min_tix or index_tuple[i] > max_tix
        for i, (min_tix, max_tix) in enumerate(tuple_limits)
    ):
        raise IndexError(f"Tuple element {index_tuple} is out of bounds!")
    return int(
        sum(
            (tuple_limits[i - 1][1] - tuple_limits[i - 1][0] + 1) ** i * (index_tuple[i] - min_tix)
            for i, (min_tix, max_tix) in enumerate(tuple_limits)
        )
    )


@cached_function
def map_int_to_tuple(index, tuple_limits: tuple, tuple_len=None, order: str = "r_l") -> tuple:
    r"""
    Map integer to tuple (the inverse of :func:`map_tuple_to_int`).

    Recovers the tuple by successive modular arithmetic on the interval lengths.

    INPUT:

    - ``index`` -- integer
    - ``tuple_limits`` -- tuple of (min, max) tuples
    - ``tuple_len`` -- integer (default: None); if positive, duplicate tuple_limits
    - ``order`` -- string (default: 'r_l'); 'l_r' or 'r_l' ordering

    OUTPUT:

    - tuple of integers

    EXAMPLES::

        sage: from maass_form_core.utils.indexing import map_int_to_tuple
        sage: map_int_to_tuple(0, ((-1, 1),), 1)
        (-1,)
        sage: map_int_to_tuple(0, ((-1, 1),), 2)
        (-1, -1)
        sage: map_int_to_tuple(0, ((-1, 1),), 3)
        (-1, -1, -1)
        sage: map_int_to_tuple(0, ((-1, 1), (-1, 1)))
        (-1, -1)
        sage: map_int_to_tuple(0, ((-1, 1), (-3, 1)))
        (-1, -3)
        sage: map_int_to_tuple(1, ((-1, 1), (-3, 1)))
        (0, -3)

    TESTS::

        sage: map_int_to_tuple(0, ((-1, 1), (-1, -2)))
        Traceback (most recent call last):
        ...
        ValueError: tuple_limits ((-1, 1), (-1, -2)) do not give positive length intervals

        sage: map_int_to_tuple(9, ((-1, 1), (-1, 1)))
        Traceback (most recent call last):
        ...
        IndexError: Index 9 is out of bounds!
        sage: map_int_to_tuple(-1, ((-1, 1), (-1, 1)))
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
        if order == "r_l":
            ix_t = [*ix_t, t + min_tix]
        else:
            ix_t = [t + min_tix, *ix_t]
        index = (index - t) / range_tix
    return tuple(ix_t)


def integer_to_bounds_tuple(m, degree) -> tuple:
    r"""
    Convert a positive integer to a tuple of symmetric bounds.

    INPUT:

    - ``m`` -- positive integer
    - ``degree`` -- positive integer

    OUTPUT:

    - tuple of tuples; ``((-m, m),) * degree``

    EXAMPLES::

        sage: from maass_form_core.utils.indexing import integer_to_bounds_tuple
        sage: integer_to_bounds_tuple(1, 2)
        ((-1, 1), (-1, 1))
        sage: integer_to_bounds_tuple(2, 3)
        ((-2, 2), (-2, 2), (-2, 2))
    """
    if m <= 0 or not isinstance(m, (Integer, int)):
        raise ValueError("m must be positive")
    if degree <= 0 or not isinstance(degree, (Integer, int)):
        raise ValueError("degree must be positive")
    return ((-m, m),) * degree
