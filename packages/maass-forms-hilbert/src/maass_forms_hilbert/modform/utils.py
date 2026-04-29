import json
import logging
from collections.abc import Iterable
from functools import wraps
from typing import TYPE_CHECKING, Any, Literal

from sage.all import CC
from sage.arith.misc import factor
from sage.categories.sets_cat import cartesian_product
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_function
from sage.misc.functional import round
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import NumberField
from sage.rings.number_field.number_field_element import NumberFieldElement
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.rings.number_field.unit_group import UnitGroup
from sage.rings.real_mpfr import RealNumber as RealNumber_class
from sage.structure.element import Matrix, Vector

if TYPE_CHECKING:
    from maass_forms_hilbert.modform.hilbert_maass_space import HilbertMaassFormSpace

try:
    from comp_manager.core.decorators import mongo_cache
except (ModuleNotFoundError, ImportError, Exception):
    # Handle the case when comp_manager is not installed or MongoDB is not available
    logging.info("Not using MongoCache!")

    def mongo_cache_dummy(func):
        r"""
        Mock of decorator for caching function results in MongoDB.

        This is a wrapper for the comp_manager.core.decorators.mongo_cache decorator.
        If comp_manager is not installed, this provides a no-op decorator instead.

        INPUT:

        - ``*args`` -- Variable length argument list
        - ``**kwargs`` -- Arbitrary keyword arguments

        OUTPUT:

        - Decorator function

        EXAMPLES::

            sage: from maass_forms_hilbert.modform.utils import mongo_cache
            ...
            sage: @mongo_cache
            ....: def my_function(*args, **kwargs):
            ....:     return args[0]^2
            sage: my_function(5)
            25
        """

        @wraps(func)
        def wrapped_func(*args: Any, **kwargs: dict[str, Any]) -> Any:
            return func(*args, **kwargs)

        return wrapped_func

    mongo_cache = mongo_cache_dummy

log = logging.getLogger(__name__)
# User defined type for either Python int or Sage Integer
Integer_t = Integer | int
Real_t = RealNumber_class | float
Complex_t = ComplexNumber | complex


@cached_function
def cartesian_product_from_M(M: tuple[tuple[Integer_t]]) -> Iterable[Vector]:
    r"""
    Return a list of all vectors in the cartesian product of integer ranges defined by M.

    INPUT:

    - ``M`` -- tuple of tuples; each inner tuple (min, max) defines an integer range

    OUTPUT:

    - List of vectors in the cartesian product of the ranges

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import cartesian_product_from_M
        sage: cartesian_product_from_M(((-1, 1), (-2, 2)))
        [(-1, -2), (-1, -1), (-1, 0), (-1, 1), (-1, 2),
         (0, -2), (0, -1), (0, 0), (0, 1), (0, 2),
         (1, -2), (1, -1), (1, 0), (1, 1), (1, 2)]
    """
    return list(cartesian_product([range(m0[0], m0[1] + 1) for m0 in M]))


def is_tuple_zero(t: tuple[Integer_t]) -> bool:
    r"""
    Check if a tuple consists entirely of zeros.

    INPUT:

    - ``t`` -- tuple of integers

    OUTPUT:

    - Boolean; True if all elements are zero, False otherwise

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import is_tuple_zero
        sage: is_tuple_zero((0, 0, 0))
        True
        sage: is_tuple_zero((0, 1, 0))
        False
    """
    return all(t0 == 0 for t0 in t)


def length_from_M(M: tuple[tuple[Integer_t]]) -> int:
    r"""
    Calculate the total number of points in the cartesian product of integer ranges.

    INPUT:

    - ``M`` -- tuple of tuples of integers; each inner tuple (min, max) defines an integer range

    OUTPUT:

    - Integer representing the total number of points

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import length_from_M
        sage: length_from_M(((-1, 1), (-2, 2)))
        15
        sage: length_from_M(((-5, 5), (-5, 5)))
        121
    """
    return prod([m0[1] - m0[0] + 1 for m0 in M])


def get_Q_from_bounds(M: tuple[tuple[Integer_t]]) -> tuple:
    r"""
    Find a bounding box for the integer coordinates based on the given bounds.

    This function creates a bounding box for the cube [-M1,M1] x [-M2,M2] x...
    for the integer coordinates corresponding to the box [-b,b]^n in the lattice.

    INPUT:

    - ``M`` -- tuple of tuples of integers; each inner tuple (min, max) defines bounds

    OUTPUT:

    - Tuple of integers representing the bounding box

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import get_Q_from_bounds
        sage: get_Q_from_bounds(((-1, 1), (-2, 2)))
        (4, 4)
        sage: get_Q_from_bounds(((-5, 5), (-10, 10)))
        (12, 12)
    """
    C = max(max(abs(b0), abs(b1)) for b0, b1 in M) + 2
    return (C,) * len(M)


@cached_function
def map_tuple_to_int(
    index_tuple: tuple, tuple_limits: tuple[tuple[Integer_t]], tuple_len: int = None
) -> int:
    r"""
    Map a tuple (a0,a1,...,a[n-1]) with min_i < ai < max_i to an integer
     $\sum_i=0^(n-1) (max_i - min_i + 1)**(n - 1 - i)*(ai - min_i)$

    NOTE: This function together with `map_int_to_tuple` provides an isomorphism between
            [min_0,...,max_0] x [min_1,...,max_1] x ... x [min_{n-1},...,max_{n-1}]
             and [0,...,N-1] where N = prod(max_i - min_i + 1).

    INPUT:

    - ``index_tuple`` -- tuple
    - ``tuple_limits`` -- tuple of tuples
    - ``tuple_len`` -- integer: number of tuples (default: None) if positive then the tuple_limits
                       are duplicated that number of times.

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import map_tuple_to_int
        sage: map_tuple_to_int((0,1), ((-5, 5),(-5,5)), 2)
        71
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
    if any(
        index_tuple[i] < min_tix or index_tuple[i] > max_tix
        for i, (min_tix, max_tix) in enumerate(tuple_limits)
    ):
        raise IndexError(f"Tuple element {index_tuple} is out of bounds!")
    # Calculate the index of the tuple
    return int(
        sum(
            (tuple_limits[i - 1][1] - tuple_limits[i - 1][0] + 1) ** i * (index_tuple[i] - min_tix)
            for i, (min_tix, max_tix) in enumerate(tuple_limits)
        )
    )


@cached_function
def map_int_to_tuple(
    index: Integer_t,
    tuple_limits: tuple[tuple[Integer_t, Integer_t], ...],
    tuple_len: Integer_t = None,
    order: str = "r_l",
) -> tuple:
    r"""
    Map integer to tuple (the inverse of map_tuple_to_int) by modding recursively
    modulo the lengths of the integer intervals.

    INPUT:

    - ``index`` -- integer
    - ``tuple_limits`` -- tuple of tuples of limits
    - ``tuple_len`` -- integer (number of tuples - duplicates the input tuple_limits)
    - ``order`` -- string: 'l_r' or 'r_l' (default: 'l_r') - If 'l_r' then the tuples run through
          the first coordinate first

    EXAMPLES::

    sage: from maass_forms_hilbert.modform.utils import map_int_to_tuple
    sage: map_int_to_tuple(0,((-1,1),), 1)
    (-1,)
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
        if order == "r_l":
            ix_t = [*ix_t, t + min_tix]
        else:
            ix_t = [t + min_tix, *ix_t]
        index = (index - t) / range_tix
    return tuple(ix_t)


@cached_function()
def number_field_basis_matrix(number_field: NumberField, prec: int = 53) -> Matrix:
    r"""
    Return the basis matrix with respect to the '.vector' property of number field elements.

    This function computes a matrix B such that for any number field element e,
    B * e.vector() == e.complex_embeddings(prec).

    INPUT:

    - ``number_field`` -- a NumberField; the field whose basis matrix we want
    - ``prec`` -- integer (default: 53); precision in bits for complex embeddings

    OUTPUT:

    - A matrix whose rows are the complex embeddings of the number field basis elements

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import number_field_basis_matrix
        sage: B = number_field_basis_matrix(QuadraticField(-1)); B
        [   1.00000000000000 -1.00000000000000*I]
        [   1.00000000000000  1.00000000000000*I]
        sage: a = QuadraticField(-1).gen()
        sage: B * a.vector() == vector(a.complex_embeddings())
        True
        sage: K = QuadraticField(5)
        sage: B = number_field_basis_matrix(K)
        sage: B
        [ 1.00000000000000 -2.23606797749979]
        [ 1.00000000000000  2.23606797749979]
        sage: a = K.gen()
        sage: B * a.vector() == vector(a.complex_embeddings())
        True

    """
    V, f, _ = number_field.vector_space()
    return matrix([f(b).complex_embeddings(prec) for b in V.basis()]).transpose()


def ideal_basis_matrix(ideal: NumberFieldFractionalIdeal, prec: int = 53) -> Matrix:
    r"""
    Return the matrix whose rows are the complex embeddings of the integral basis elements of the ideal.

    INPUT:

    - ``ideal`` -- NumberFieldFractionalIdeal; the ideal whose basis matrix we want
    - ``prec`` -- integer (default: 53); precision in bits for complex embeddings

    OUTPUT:

    - A matrix whose columns contain the complex embeddings of the ideal's integral basis elements

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import ideal_basis_matrix
        sage: K = QuadraticField(5)
        sage: I = K.ideal(2)
        sage: B = ideal_basis_matrix(I)
        sage: B
        [ 2.00000000000000  -3.23606797749979]
        [ 2.00000000000000   1.23606797749979]
        sage: all(abs(B[i,j] - I.integral_basis()[j].complex_embeddings()[i]) < 1e-10
        ....:     for i in range(2) for j in range(2))
        True
    """
    return matrix([b.complex_embeddings(prec) for b in ideal.integral_basis()]).transpose()


@cached_function()
def dual_ideal(ideala: NumberFieldFractionalIdeal) -> NumberFieldFractionalIdeal:
    r"""
    Return the dual ideal of a given fractional ideal.

    This function computes the dual ideal of the given ideal, defined as
    $\mathfrak{a}^{-1}\mathcal{D}^{-1}$ where $\mathcal{D}$ is the different ideal.

    INPUT:

    - ``ideala`` -- NumberFieldFractionalIdeal; the ideal whose dual we want to compute

    OUTPUT:

    - NumberFieldFractionalIdeal; the dual ideal

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import dual_ideal
        sage: K = QuadraticField(5)
        sage: I = K.ideal(2)
        sage: D = dual_ideal(I)
        sage: D * I * K.different()
        Fractional ideal (1)
        sage: dual_ideal(K.ideal(1)) == K.different()**(-1)
        True
    """
    return ideala**-1 * ideala.number_field().different() ** -1


@cached_function()
def dual_ideal_basis_matrix(ideal: NumberFieldFractionalIdeal, prec: int = 53) -> Matrix:
    r"""
    Return the basis matrix of the dual ideal.

    This function computes the basis matrix of the dual ideal with respect to the
    complex embeddings. The dual ideal is defined as $\mathfrak{a}^{-1}\mathcal{D}^{-1}$
    where $\mathcal{D}$ is the different ideal.

    INPUT:

    - ``ideal`` -- NumberFieldFractionalIdeal; the ideal whose dual basis matrix we want
    - ``prec`` -- integer (default: 53); precision in bits for complex embeddings

    OUTPUT:

    - A matrix whose rows are the complex embeddings of the dual ideal's integral basis elements

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import dual_ideal_basis_matrix
        sage: K = QuadraticField(5)
        sage: I = K.ideal(1)
        sage: B = dual_ideal_basis_matrix(I)
        sage: B
        [ 1.00000000000000 0.276393202250021]
        [ 1.00000000000000 0.723606797749979]

    The dual ideal basis matrix times the ideal basis matrix should give the identity::

        sage: from maass_forms_hilbert.modform.utils import ideal_basis_matrix
        sage: D = dual_ideal_basis_matrix(I)
        sage: M = ideal_basis_matrix(I)
        sage: abs((D * M).det() - 1) < 1e-10
        True
    """
    dual = ideal**-1 * ideal.number_field().different() ** -1
    return ideal_basis_matrix(dual, prec)


def ideal_coordinates(
    ideala: NumberFieldFractionalIdeal, element: NumberFieldElement, check: bool = False
):
    """
    Find the coordinates of an ideal element with respect to an integral basis of that ideal.

    INPUT:
     - `ideala` -- an ideal
     - `element` -- an element of `ideala`
     - `check` -- if True, check that the coordinates are correct

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import (dual_ideal, ideal_generator,
        ....:   ideal_coordinates, number_field_basis_matrix, ideal_basis_matrix)
        sage: from sage.rings.number_field.number_field import QuadraticField
        sage: ideala_dual = dual_ideal(QuadraticField(3).ideal(1))
        sage: delta = ideal_generator(ideala_dual)
        sage: coords = ideal_coordinates(ideala_dual, delta, check=False)
        sage: coords
        (0, 1)
        sage: sum([c * ideala_dual.integral_basis()[i] for i, c in enumerate(coords)]) == delta
        True
        sage: ideala_dual = dual_ideal(QuadraticField(5).ideal(1))
        sage: delta = ideal_generator(ideala_dual)
        sage: ideal_coordinates(ideala_dual, delta, check=False)
        (1, -1)

    """
    nf = ideala.number_field()
    if not isinstance(element, NumberFieldElement):
        element = nf(element)
    if element not in ideala:
        raise ValueError(f"Element {element} not in ideal: {ideala}")
    basis_change_matrix = ideal_basis_matrix(ideala) ** -1 * number_field_basis_matrix(nf)
    coordinates = basis_change_matrix * element.vector()
    eps = basis_change_matrix.base_ring().epsilon() * 2**3
    coordinates_int = tuple(round(c.real()) for c in coordinates if abs(c.imag()) < eps)
    if len(coordinates_int) != len(coordinates):
        raise ArithmeticError(
            f"Can not find lattice coordinates for delta={element}."
            f" coordinates={coordinates}"
            f" coordinates_int={coordinates_int}"
        )
    if check:
        assert (
            sum([c * ideala.integral_basis()[i] for i, c in enumerate(coordinates_int)]) == element
        )
    return coordinates_int


@cached_function()
def dual_ideal_element(
    coordinates: tuple[Integer_t] or vector,
    ideal: NumberFieldFractionalIdeal,
    as_nf_element=False,
):
    r"""
    Create an element in the dual ideal from the given coordinates.

    INPUT:

    - ``coordinates`` -- tuple of integers or vector; coordinates with respect to the dual ideal basis
    - ``ideal`` -- NumberFieldFractionalIdeal; the ideal whose dual will contain the element
    - ``as_nf_element`` -- boolean (default: False); if True, return the element as a number field element

    OUTPUT:

    - Either a vector of complex embeddings of the element (if as_nf_element=False)
      or the number field element itself (if as_nf_element=True)

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import dual_ideal_element, dual_ideal
        sage: K = QuadraticField(5)
        sage: ideal = K.ideal(1)
        sage: dual_ideal_element((1, 2), ideal, as_nf_element=True)
        1/5*a + 2
        sage: v = dual_ideal_element((1, 2), ideal)
        sage: v  # complex embeddings of the element
        (1.55278640450004, 2.44721359549996)
    """
    if not as_nf_element:
        return dual_ideal_basis_matrix(ideal) * vector(coordinates)
    dual = ideal**-1 * ideal.number_field().different() ** -1
    return sum([c * dual.integral_basis()[i] for i, c in enumerate(coordinates)])


def ideal_generator(ideala: NumberFieldFractionalIdeal) -> NumberFieldFractionalIdeal:
    r"""
    Find a totally positive generator for an ideal if possible, else use a reduced generator.

    INPUT:

    - ``ideala`` -- NumberFieldFractionalIdeal

    OUTPUT:

    - A generator for the ideal, preferably totally positive if possible.

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import ideal_generator
        sage: K = QuadraticField(5)
        sage: I = K.ideal(1)
        sage: g = ideal_generator(I)
        sage: I == K.fractional_ideal(g)
        True
        sage: g.is_totally_positive()
        True
        sage: I2 = K.ideal(2)
        sage: g2 = ideal_generator(I2)
        sage: I2 == K.fractional_ideal(g2)
        True
    """
    narrow_class_number = ideala.number_field().narrow_class_group().order()
    if narrow_class_number == 1:
        u = UnitGroup(ideala.number_field()).gens_values()[1]
        x = ideala.gens_reduced()[0]
        # f = space.number_field().galois_group()
        test = [x, x * u, -x * u]
        delta = None
        for delta_test in test:
            if not (delta_test).is_totally_positive():
                continue
            if ideala != ideala.number_field().fractional_ideal(delta_test):
                continue
            delta = delta_test
            break
        if not delta:
            raise ArithmeticError(f"Cannot find a totally positive generator for {ideala}")
    else:
        delta = ideala.gens_reduced()[0]
    return delta


def complex_number_to_json(s: ComplexNumber) -> dict:
    r"""
    Convert a ComplexNumber to a JSON-serializable dictionary.

    INPUT:

    - ``s`` -- ComplexNumber or coercible value

    OUTPUT:

    - dict; JSON representation of the complex number

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import complex_number_to_json, complex_number_from_json
        sage: from sage.rings.complex_mpfr import ComplexField
        sage: CF = ComplexField(53)
        sage: s = CF(1 + 2j)
        sage: d = complex_number_to_json(s)
        sage: d['prec']
        53
        sage: 'val' in d
        True
        sage: complex_number_from_json(d) == s
        True
    """
    if is_json_number(s):
        return s
    if not isinstance(s, ComplexNumber):
        s = CC(s)
    return {"prec": s.parent().prec(), "val": str(s)}


def is_json_number(data: dict | str | Any) -> bool:
    r"""
    Check if the input is a JSON representation of a complex number.

    INPUT:

    - ``data`` -- dict, str, or any; the data to check

    OUTPUT:

    - Boolean; True if the data is a JSON complex number, False otherwise

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import is_json_number
        sage: is_json_number({'prec': 53, 'val': '1.0'})
        True
        sage: is_json_number({'foo': 1})
        False
        sage: is_json_number('not a dict')
        False
    """
    return isinstance(data, dict) and "prec" in data and "val" in data and len(data.keys()) == 2


def complex_number_from_json(json_complex: dict | str) -> ComplexNumber:
    r"""
    Create a ComplexNumber from JSON data.

    This function takes either a JSON string or a dictionary containing the precision
    and value representation of a complex number and returns a SageMath ComplexNumber.

    INPUT:

    - ``json_complex`` -- dict or str; either a dictionary with 'prec' and 'val' keys,
      or a JSON string representing such a dictionary

    OUTPUT:

    - ComplexNumber; the reconstructed complex number

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import complex_number_from_json
        sage: complex_number_from_json({'prec': 53, 'val': '1.0 + 2.0*I'})
        1.00000000000000 + 2.00000000000000*I
        sage: complex_number_from_json('{"prec": 53, "val": "1.0"}')
        1.00000000000000
        sage: complex_number_from_json({'prec': 100, 'val': '1.0'})  # higher precision
        1.0000000000000000000000000000

    TESTS::

        sage: complex_number_from_json('not a json')
        Traceback (most recent call last):
        ...
        JSONDecodeError: ...
    """
    if isinstance(json_complex, str):
        json_complex = json.loads(json_complex)
    try:
        return ComplexField(json_complex["prec"])(json_complex["val"])
    except TypeError:
        return ComplexField(json_complex)


def complex_tuple_to_json(complex_tuple: tuple[ComplexNumber]) -> list[dict]:
    r"""
    Convert a tuple of complex numbers to a JSON-serializable list.

    INPUT:

    - ``complex_tuple`` -- tuple of ComplexNumber; the complex numbers to convert

    OUTPUT:

    - list of dicts; JSON representation of the complex numbers

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import complex_tuple_to_json
        sage: from sage.rings.complex_mpfr import ComplexField
        sage: CF = ComplexField(53)
        sage: t = (CF(1), CF(2 + 1j))
        sage: j = complex_tuple_to_json(t)
        sage: len(j)
        2
        sage: j[0]['prec']
        53
        sage: j[1]['val']
        '2.00000000000000 + 1.00000000000000*I'
    """
    return [complex_number_to_json(s) for s in complex_tuple]


def complex_tuple_from_json(json_list: list[dict] | str) -> tuple[ComplexNumber]:
    r"""
    Convert a JSON list to a tuple of complex numbers.

    INPUT:

    - ``json_list`` -- list of dicts or a string; JSON representation of complex numbers

    OUTPUT:

    - Tuple of ComplexNumber objects

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import complex_tuple_to_json, complex_tuple_from_json
        sage: from sage.rings.complex_mpfr import ComplexField
        sage: CF = ComplexField(53)
        sage: t = (CF(1), CF(2))
        sage: j = complex_tuple_to_json(t)
        sage: complex_tuple_from_json(j)
        (1.00000000000000, 2.00000000000000)
        sage: complex_tuple_from_json('[{"prec": 53, "val": "1.0"}]')
        (1.00000000000000,)
    """
    if isinstance(json_list, str):
        json_list = json.loads(json_list)
    return tuple(complex_number_from_json(s) for s in json_list)


def number_field_to_json(nf: NumberField) -> dict:
    r"""
    Return a JSON representation of a number field.

    NOTE: Any information about embeddings is ignored.

    INPUT:

    - ``nf`` -- NumberField; the number field to serialize

    OUTPUT:

    - dict; JSON-serializable dictionary

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import number_field_to_json, number_field_from_json
        sage: K = QuadraticField(5, 'a')
        sage: d = number_field_to_json(K)
        sage: d['polynomial']
        'x^2 - 5'
        sage: K2 = number_field_from_json(d)
        sage: K2.polynomial()
        x^2 - 5
    """
    return {"polynomial": str(nf.polynomial()), "names": list(nf._names)}


def number_field_from_json(data: dict | str) -> NumberField:
    r"""
    Create a number field from JSON data.

    INPUT:

    - ``data`` -- dict or str; JSON representation of a number field

    OUTPUT:

    - NumberField; the reconstructed number field

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import number_field_to_json, number_field_from_json
        sage: K = QuadraticField(5, 'a')
        sage: d = number_field_to_json(K)
        sage: K2 = number_field_from_json(d)
        sage: K2.polynomial()
        x^2 - 5
    """
    if isinstance(data, str):
        data = json.loads(data)
    return NumberField(ZZ["x"](data["polynomial"]), names=tuple(data["names"]))


def coefficient_dict_to_json(coeff_dict: dict) -> dict:
    r"""
    Convert a coefficient dict to JSON format.

    INPUT:

    - ``coeff_dict`` -- dict of coefficients with keys as integer tuples

    OUTPUT:

    - dict; JSON-serializable dictionary with stringified keys and complex number values

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import coefficient_dict_to_json
        sage: data = { (0, 0): 1, (1, 1): 1, (2, 2): 1 }
        sage: coefficient_dict_to_json(data)
         {'[0, 0]': {'prec': 53, 'val': '1.00000000000000'},
         '[1, 1]': {'prec': 53, 'val': '1.00000000000000'},
         '[2, 2]': {'prec': 53, 'val': '1.00000000000000'}}
        sage: from sage.rings.complex_mpfr import ComplexField
        sage: CF = ComplexField(103)
        sage: data = { (0, 0): 1, (1, 1): CF(1), (2, 2): CF(0,1) }
        sage: coefficient_dict_to_json(data)
         {'[0, 0]': {'prec': 53, 'val': '1.00000000000000'},
          '[1, 1]': {'prec': 103, 'val': '1.00000000000000000000000000000'},
          '[2, 2]': {'prec': 103, 'val': '1.00000000000000000000000000000*I'}}

    """
    if not coeff_dict:
        return {}
    if all(isinstance(x, str) and is_json_number(x) for x in coeff_dict.items()):
        return coeff_dict
    return {
        json.dumps([int(ki) for ki in k]): complex_number_to_json(v) for k, v in coeff_dict.items()
    }


def coefficient_dict_from_json(data: dict | str) -> dict:
    r"""
    Convert a coefficient dict from JSON format.

    INPUT:

    - ``data`` -- dict or str; JSON data as dict or string

    OUTPUT:

    - dict; coefficient dictionary with tuple keys and ComplexNumber values

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import (coefficient_dict_from_json,
        ....:                                       coefficient_dict_to_json)
        sage: data = { (0, 0): 1, (1, 1): 1, (2, 2): 1 }
        sage: json_string = coefficient_dict_to_json(data)
        sage: coefficient_dict_from_json(json_string)
        {(0, 0): 1.00000000000000, (1, 1): 1.00000000000000, (2, 2): 1.00000000000000}
        sage: from sage.rings.complex_mpfr import ComplexField
        sage: CF = ComplexField(103)
        sage: data = { (0, 0): 1, (1, 1): CF(1), (2, 2): CF(0,1) }
        sage: json_string = coefficient_dict_to_json(data)
        sage: coefficient_dict_from_json(json_string)
         {(0, 0): 1.00000000000000,
         (1, 1): 1.00000000000000000000000000000,
         (2, 2): 1.00000000000000000000000000000*I}

    """
    if isinstance(data, str):
        data = json.loads(data)
    return {tuple(json.loads(k)): complex_number_from_json(v) for k, v in data.items()}


def integer_to_bounds_tuple(m: Integer_t, degree: Integer_t) -> tuple[tuple[Integer_t]]:
    """
    Convert a positive integer to a tuple of bounds.

    INPUT:

    -``m`` -- positive integer
    -``degree`` -- positive integer

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import integer_to_bounds_tuple
        sage: integer_to_bounds_tuple(1, 2)
        ((-1, 1), (-1, 1))
        sage: integer_to_bounds_tuple(2, 3)
        ((-2, 2), (-2, 2), (-2, 2))
    """
    if m <= 0 or not isinstance(m, Integer_t):
        raise ValueError("m must be positive")
    if degree <= 0 or not isinstance(m, Integer_t):
        raise ValueError("degree must be positive")
    return ((-m, m),) * degree


def unit_relations(space: "HilbertMaassFormSpace", m: Integer_t = 6):
    """
    produce the sets of integer lattice points which are related by the
    automorphy a(unit^2 v)=a(v) for the quadratic real fields.

    INPUT:

    -``m`` -- positive integer
    -``space`` -- HilbertMaassFormSpace

    EXAMPLES::

        sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
        sage: from maass_forms_hilbert.modform.utils import unit_relations
        sage: space = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
        sage: unit_relations(space, 6)  # doctest: +ELLIPSIS
        [[(-6, -6), (-6, 6)],
         [(-5, -6), (-3, 2)],
         [(-5, -5), (-5, 5)],
         [(-5, 6), (-3, -2)],
         [(-4, -6), (0, -2), (4, -6)],
         ...]
    """
    ideala = space.number_field().ideal(1)
    dual_ideala = dual_ideal(ideala)
    store1 = cartesian_product_from_M(((-m, m), (-m, m)))
    store = []
    for r in store1:
        if abs(r[0]) <= m and abs(r[1]) <= m:
            store.append(r)
    store.remove((0, 0))
    u = UnitGroup(space.number_field()).gens_values()[1]
    unit = u**2
    temp = []
    while store != []:
        r = store[0]
        r_element = dual_ideal_element(r, ideala, as_nf_element=True)
        d = r
        x = 0
        use = r_element
        kemp = []
        while abs(d[0]) <= m and abs(d[1]) <= m:
            kemp.append(tuple(d))
            if d in store:
                store.remove(tuple(d))
            use = use * unit
            d = ideal_coordinates(dual_ideala, use)
            x = x + 1
        use = r_element * unit ** (-1)
        d = ideal_coordinates(dual_ideala, use)
        while abs(d[0]) <= m and abs(d[1]) <= m:
            kemp.append(d)
            use = use * unit ** (-1)
            if d in store:
                store.remove((d[0], d[1]))
            d = ideal_coordinates(dual_ideala, use)
            x = x + 1
        if x <= 1:
            kemp.pop()
        else:
            temp.append(kemp)
    return temp


def hecke_relations_coprime(space: "HilbertMaassFormSpace", m: Integer_t = 6):
    """
    produce the sets of integer lattice points which are related by the
    the Hecke relation of the form a(delta m)a(delta n)=a(delta mn), where delta is a generator of
    dual ideal

    INPUT:

    -``m`` -- positive integer
    -``space`` -- HilbertMaassFormSpace

    EXAMPLES::

        sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
        sage: from maass_forms_hilbert.modform.utils import hecke_relations_coprime
        sage: space = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
        sage: hecke_relations_coprime(space, 6) # long time
          [[(-5, -5), (-1, 0), (5, 0)],
           [(-5, -5), (1, 0), (-5, 0)],
           [(-5, 5), (-3, -4), (5, 0)],
           [(-5, 5), (3, 4), (-5, 0)],
        ...
           [(5, 5), (1, 0), (5, 0)]]
    """
    ideala = space.number_field().ideal(1)
    dual_ideala = dual_ideal(ideala)
    t = ideal_generator(dual_ideala)
    store1 = cartesian_product_from_M(((-m, m), (-m, m)))
    store = []
    for r in store1:
        if abs(r[0]) <= m and abs(r[1]) <= m:
            store.append(r)
    store.remove((0, 0))
    temp = []
    for r in store:
        for s in store:
            r_element = dual_ideal_element(r, ideala, as_nf_element=True) / t
            r_ideal = space.number_field().ideal(r_element)
            s_element = dual_ideal_element(s, ideala, as_nf_element=True) / t
            s_ideal = space.number_field().ideal(s_element)
            if r != s and r_ideal != ideala and s_ideal != ideala and r_ideal + s_ideal == ideala:
                y = t * r_element * s_element
                d = ideal_coordinates(dual_ideala, y)
                if (d[0], d[1]) in store:
                    temp.append([r, s, d])
    return temp


def hecke_relations_prime_power(space: "HilbertMaassFormSpace", m: Integer_t = 6):
    """
    produce the sets of integer lattice points which are related by the
    the Hecke relation of the form a(delta p^n)=a(delta p)a(delta p^{n-1})-a(delta p^{n-2}), where delta is a generator of
    dual ideal

    INPUT:

    -``m`` -- positive integer
    -``space`` -- HilbertMaassFormSpace

    EXAMPLES::

        sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
        sage: from maass_forms_hilbert.modform.utils import hecke_relations_prime_power
        sage: space = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
        sage: hecke_relations_prime_power(space, 6)
        [[(-2, -1), (5, -1)],
         [(-1, -3), (5, 1)],
         [(-1, -2), (2, 2), (-2, -4), (4, 4)],
         [(-1, 0), (2, -2)],
         [(1, 0), (2, -2)],
         [(1, 2), (2, 2), (2, 4), (4, 4)],
         [(1, 3), (5, 1)],
         [(2, 1), (5, -1)]]
    """
    ideala = space.number_field().ideal(1)
    dual_ideala = dual_ideal(ideala)
    t = ideal_generator(dual_ideala)
    # f = space.number_field().galois_group()
    store1 = cartesian_product_from_M(((-m, m), (-m, m)))
    store = []
    for r in store1:
        if abs(r[0]) <= m and abs(r[1]) <= m:
            store.append(r)
    store.remove((0, 0))
    temp = []
    for r in store:
        r_element = dual_ideal_element(r, ideala, as_nf_element=True) / t
        r_ideal = space.number_field().ideal(r_element)
        if r_ideal.is_prime():
            a = r_element
            x = 1
            d = r
            kemp = []
            while (d[0], d[1]) in store:
                kemp.append(d)
                a = a * r_element
                y = a * t
                d = ideal_coordinates(dual_ideala, y)
                x = x + 1
            if x <= 2:
                kemp.pop()
            else:
                temp.append(kemp)
    return temp


def symmetric_relations(space: "HilbertMaassFormSpace"):
    """
    produces a set of integer lattice points which are related by
    reflection relation (flipping operators) a(u_p v)=ta(v), where u is a unit which is negative
    at the place (prime) p and positive at other places.

    INPUT:

    -``space`` -- HilbertMaassFormSpace

    EXAMPLES::

        sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
        sage: from maass_forms_hilbert.modform.utils import symmetric_relations
        sage: space = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
        sage: symmetric_relations(space)
        [(1, 1), (0, -1), (0, 1)]

        Comment: The first element in the output represents the positive generator of the dual ideal.
        These three elements are non-zero for a hilbert Maass form.

        sage: space = HilbertMaassFormSpace(QuadraticField(5), cuspidal=True)
        sage: symmetric_relations(space)
        [(1, -1), (1, -2), (-1, 2)]

        sage: space = HilbertMaassFormSpace(QuadraticField(17), cuspidal=True)
        sage: symmetric_relations(space)
        [(5, -8), (1, -2), (-1, 2)]

        sage: space = HilbertMaassFormSpace(QuadraticField(41), cuspidal=True)
        sage: symmetric_relations(space)
        [(37, -64), (1, -2), (-1, 2)]

        Comments: For QuadraticField with bigger discriminant we have to choose the smallest of these
        3 and assign the value 1 for constructing object using hejhal' algorithm.

    """
    ideala = space.number_field().ideal(1)
    dual_ideala = dual_ideal(ideala)
    t = ideal_generator(dual_ideala)
    u = UnitGroup(space.number_field()).gens_values()[1]
    ideala = space.number_field().ideal(1)
    if u > 0:
        u = -u
    set_check = [t, t * u**-1, -t * u**-1]
    kemp = []
    for use in set_check:
        d = ideal_coordinates(dual_ideala, use)
        kemp.append(d)
    return kemp


def bilinear_form(space: "HilbertMaassFormSpace", x: "NumberFieldElement", y: "NumberFieldElement"):
    """
    produce the value of bilinear forms f(x, y)=x_1y_1+...+x_ny_n

    INPUT:
    -``space`` -- HilbertMaassFormSpace
    -``x, y``  ---NumberFieldElement, NumberFieldElement

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import bilinear_form
        sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
        sage: space = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
        sage: a=space.number_field().gen()
        sage: bilinear_form(space, a, a)
        4

    """
    f = space.number_field().galois_group()
    t = len(f)
    kep = 0
    for r in range(0, t):
        s = f[r](x) * f[r](y)
        kep = kep + s
    return kep


def best_hecke_relation(
    space: "HilbertMaassFormSpace",
    m: Integer_t = 6,
    check: Literal["coprime", "prime_power", "unit", "symmetric"] = "unit",
    epsilon: Integer_t = 25,
    same_norm: bool = False,
):
    """
    produce the best hecke relation in the sense those with smallest bilinear norm value
    less than epsilon.

    INPUT:
    -``space`` -- HilbertMaassFormSpace
    -``m ``  --Integer_t  ( bound)
    -``check`` --'{coprime, prime_power, unit}' Enter one value out of these three
    -``epsilon`` --Integer_t=25
    -`` same_norm`` bool (produces the hecke relation of same norm if it set to True)

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import best_hecke_relation
        sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
        sage: space = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
        sage: best_hecke_relation(space = space, check = 'coprime', same_norm = True) # slow
         [[(-1, -3), (-1, 0), (-1, 4)],
          [(-1, -3), (1, 0), (1, -4)],
          [(-1, 0), (-1, -3), (-1, 4)],
         ...
          [(3, 5), (2, -2), (-2, 6)]]
        sage: best_hecke_relation(space = space, check = 'coprime')  # doctest: +ELLIPSIS
         [[(-1, -3), (-1, 0), (-1, 4)],
          [(-2, -1), (-1, -2), (1, 4)],
          [(-1, -3), (-1, -2), (3, 2)],
         ...
          [(-3, -5), (-2, 2), (-2, 6)]]
    """

    if check == "coprime":
        ak = hecke_relations_coprime(space, m)
    elif check == "prime_power":
        ak = hecke_relations_prime_power(space, m)
    elif check == "unit":
        ak = unit_relations(space, m)
    elif check == "symmetric":
        ak = symmetric_relations(space)
    else:
        raise ValueError("Enter one value out of 'coprime', 'prime_power', 'unit', 'symmetric'")

    def custom_function(x, space):
        total = 0
        ideala = space.number_field().ideal(1)
        for r in x:
            t = dual_ideal_element(r, ideala, as_nf_element=True)
            total = total + bilinear_form(space, t, t)
        return total

    paired_elements = [(x, custom_function(x, space)) for x in ak]
    sorted_paired_elements = sorted(paired_elements, key=lambda pair: pair[1])
    value_sorted_paired_elements = []
    if same_norm:
        for pair in sorted_paired_elements:
            x, value = pair
            if value < epsilon:
                value_sorted_paired_elements.append(pair)
    else:
        for pair in sorted_paired_elements:
            x, value = pair
            t = 1
            for pair1 in value_sorted_paired_elements:
                x1, value1 = pair1
                if value == value1:
                    t = 0
                    break
            if value < epsilon and t == 1:
                value_sorted_paired_elements.append(pair)
    sorted_elements = [pair[0] for pair in value_sorted_paired_elements]
    return sorted_elements


def ideal_factors(ida):
    r"""
    Compute all ideal factors of a given ideal.

    INPUT:

    - ``ida`` -- an ideal (NumberFieldFractionalIdeal)

    OUTPUT:

    - list; all ideal factors of ``ida``

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import QuadraticField
        sage: from maass_forms_hilbert.modform.utils import ideal_factors
        sage: K = QuadraticField(5)
        sage: I = K.ideal(6)
        sage: F = ideal_factors(I)
        sage: all((I / f).is_integral() for f in F)
        True
        sage: K.ideal(1) in F
        True
        sage: K.ideal(6) in F
        True
    """
    prime_factors = [x[0] for x in factor(ida)]
    exponents = [range(x[1] + 1) for x in factor(ida)]
    new_exponents = list(cartesian_product(exponents))
    factors = []
    for ex in new_exponents:
        idb = prod([p ** ex[i] for i, p in enumerate(prime_factors)])
        factors.append(idb)
    return factors


def list_of_lists_to_limits(list_of_lists: list[list[int]]) -> tuple[tuple[int, int], ...]:
    """
    Convert list of lists to tuple of tuples

    INPUT:

    - ``list_of_lists`` -- list of lists with integer entries

    OUTPUT:

    - tuple of tuples consisting of max and min of each list

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import list_of_lists_to_limits
        sage: list_of_lists_to_limits([[1,2,3],[4,5,6]])
        ((1, 3), (4, 6))
        sage: list_of_lists_to_limits([[-2, -1, 0, 1, 2], [-1, 0, 1, 2, 3]])
        ((-2, 2), (-1, 3))
    """
    min_max = [(int(min(x)), int(max(x))) for x in list_of_lists]
    return tuple(min_max)


def coefficient_dict_to_matrix(coefficient_dict: dict[tuple[int, ...], Any]):
    """
    Create a matrix of coefficients from dictionary

    INPUT:

    - ``coefficient_dict`` -- dictionary mapping tuples to coefficients

    OUTPUT:

    - matrix; matrix of coefficients

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.utils import coefficient_dict_to_matrix
        sage: coefficient_dict_to_matrix({(1, 1, 1): 1, (2, 2, 2): 2})
        Traceback (most recent call last):
        ...
        ValueError: Coefficient dict does not have a contiguous set of keys.
        sage: dict = {(1, 1, 1): 1, (2, 1, 1): 2, (1, 2, 1): 3, (2, 2, 1): 4, (1, 1, 2): 5,
        ....:  (2, 1, 2): 6,  (1, 2, 2): 7, (2, 2, 2): 8}
        sage: coefficient_dict_to_matrix(dict).transpose()
        [1 2 3 4 5 6 7 8]
    """
    degree = len(next(iter(coefficient_dict)))
    index_slices = [[x[i] for x in coefficient_dict] for i in range(degree)]
    M = list_of_lists_to_limits(index_slices)
    nmax = prod([t[1] + 1 - t[0] for t in M])
    try:
        coefficients = [[coefficient_dict[map_int_to_tuple(n, M)]] for n in range(nmax)]
    except KeyError as err:
        raise ValueError("Coefficient dict does not have a contiguous set of keys.") from err
    return matrix(coefficients)
