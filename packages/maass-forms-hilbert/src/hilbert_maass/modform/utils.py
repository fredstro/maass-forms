import json
from typing import Iterable

from hilbert_modgroup.pullback import HilbertPullback
from sage.all import ZZ, CC
from sage.categories.sets_cat import cartesian_product
from sage.functions.other import ceil
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_function
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.complex_mpfr import ComplexNumber, ComplexField
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.number_field.number_field import NumberField
from sage.rings.number_field.number_field_element import NumberFieldElement
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.rings.real_lazy import RLF
from sage.rings.real_mpfr import RealNumber as RealNumber_class
from sage.structure.element import Matrix, Vector

# User defined type for either Python int or Sage Integer
Integer_t = Integer | int
Real_t = RealNumber_class | float
Complex_t = ComplexNumber | complex

@cached_function
def cartesian_product_from_M(M: tuple[tuple[Integer_t]]) -> Iterable[Vector]:
    return [v for v in
            cartesian_product([range(m0[0], m0[1] + 1) for m0 in M])]


def is_tuple_zero(t: tuple[Integer_t]) -> bool:
    return all(t0 == 0 for t0 in t)

def length_from_M(M: tuple[tuple[Integer_t]]) -> int:
    """
    
    INPUT:

    - ``M`` -- tuple of tuples of integers

    Examples::

    """
    return prod([m0[1] - m0[0] + 1 for m0 in M])


def get_Q_from_bounds(P: HilbertPullback, M: tuple[tuple[Integer_t]]) -> tuple:
    """
    Find a bounding box for the cube [-M1,M1]x[-M2,M2],... for the integer coordinates correponding to the box [-b,b]^n in the lattice.
    """
    C = 1
    for ida in P.group().ideal_cusp_representatives():
        t = matrix(P.basis_matrix_ideal(ida)).transpose().norm(Infinity)
        if t > C:
            C = t
    C = ceil(C)
    C = C * max(max(abs(b0), abs(b1)) for b0, b1 in M)
    return (C,) * len(M)


@cached_function
def map_tuple_to_int(index_tuple: tuple, tuple_limits: tuple[tuple[Integer_t]],
                     tuple_len: int = None) -> int:
    r"""
    Map a tuple (a0,a1,...,a[n-1]) with min_i < ai < max_i to an integer
     $\sum_i=0^(n-1) (max_i - min_i + 1)**(n - 1 - i)*(ai - min_i)$

    NOTE: This function together with `map_int_to_tuple` provides an isomorphism between
            [min_0,...,max_0] x [min_1,...,max_1] x ... x [min_{n-1},...,max_{n-1}]
             and [0,...,N] where N = prod(max_i - min_i + 1).

    INPUT:

    - ``index_tuple`` -- tuple
    - ``tuple_limits`` -- tuple of tuples
    - ``tuple_len`` -- integer: number of tuples (default: None) if positive then the tuple_limits
                       are duplicated that number of times.
    EXAMPLES::

        sage: from hilbert_maass.modform.utils import map_tuple_to_int
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

    sage: from hilbert_maass.modform.utils import map_int_to_tuple
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
        ix_t.append(t + min_tix)
        index = (index - t) / range_tix
    return tuple(ix_t)

@cached_function()
def number_field_basis_matrix(number_field: NumberField, prec: int = 53) -> Matrix:
    """
    Basis matrix with respect to the '.vector' property of number field elements.

    Note: number_field_basis_matrix(nf, prec)*elt.vector() == elt.complex_embeddings(prec)
    """
    V, f, _ = number_field.vector_space()
    return matrix([
        f(b).complex_embeddings(prec) for b in V.basis()
    ]).transpose()


def ideal_basis_matrix(ideal: NumberFieldFractionalIdeal, prec: int = 53) -> Matrix:
    """
    :param number_field:
    :param prec:
    :return:
    """
    return matrix([
                   b.complex_embeddings(prec)
                   for b in ideal.integral_basis()
                   ]).transpose()

@cached_function()
def dual_ideal(ideala: NumberFieldFractionalIdeal) -> NumberFieldFractionalIdeal:
    return ideala ** -1 * ideala.number_field().different() ** -1

@cached_function()
def dual_ideal_basis_matrix(ideal: NumberFieldFractionalIdeal, prec: int = 53) -> Matrix:
    """
    :param number_field:
    :param prec:
    :return:
    """
    dual = ideal ** -1 * ideal.number_field().different() ** -1
    return ideal_basis_matrix(dual, prec)


def ideal_coordinates(ideala: NumberFieldFractionalIdeal,
                      element: NumberFieldElement,
                      check: bool = False):
    """
    Find the coordinates of an ideal element with respect to an integral basis of that element.

    :param ideala:
    :param element:
    :return:
    """
    nf = ideala.number_field()
    if not isinstance(element, NumberFieldElement):
        element = nf(element)
    if element not in ideala:
        raise ValueError(f"Element {element} not in ideal: {ideala}")
    basis_change_matrix = ideal_basis_matrix(ideala) ** -1 * number_field_basis_matrix(nf)
    ideal_coordinates = basis_change_matrix * element.vector()
    if check:
        assert sum([c * ideala.integral_basis[i]
                    for i, c in enumerate(ideal_coordinates)]) == element
    coordinates = basis_change_matrix * element.vector()
    coordinates_int = tuple(int(c.real()) for c in coordinates if c.imag() == 0)
    if len(coordinates_int) != len(coordinates):
        raise ArithmeticError(f"Can not find lattice coordinates for delta={element}."
                              f" coordinates={coordinates}"
                              f" coordinates_int={coordinates_int}")
    return coordinates_int


@cached_function()
def dual_ideal_element(coordinates: tuple[Integer_t] or vector,
                       ideal: NumberFieldFractionalIdeal,
                       as_nf_element=False):
    """
    Element in dual ideal given by coordinates.

    INPUT:

    - `
    """
    if not as_nf_element:
        return dual_ideal_basis_matrix(ideal)*vector(coordinates)
    dual = ideal ** -1 * ideal.number_field().different() ** -1
    return sum([c * dual.integral_basis()[i] for i, c in enumerate(coordinates)])

def totally_positive_generator(ideala: NumberFieldFractionalIdeal) -> NumberFieldFractionalIdeal:
    """
    Find a totally positive generator for an ideal.

    INPUT:

    - ``ideala`` -- NumberFieldFractionalIdeal


    """
    x, y = ideala.gens_two()
    delta = None
    for delta_test in [x + y, x - y, -x - y, -x + y]:
        if not delta_test.is_totally_positive():
            continue
        if ideala != ideala.number_field().fractional_ideal(delta_test):
            continue
        delta = delta_test
        break
    if not delta:
        raise ArithmeticError(f"Cannot find a totally positive generator for {ideala_dual}")
    return delta

def complex_number_to_json(s: ComplexNumber) -> dict:
    """
    Create json data from a ComplexNumber.

    INPUT:

    - ``s`` -- complex number
    """
    if not isinstance(s, ComplexNumber):
        s = CC(s)
    return {'prec': s.parent().prec(), 'val': str(s)}


def complex_number_from_json(json_complex: dict | str) -> ComplexNumber:
    """
    Create a ComplexNumber from json data.

    INPUT:

    - ``s`` -- string or dict representing a complex number
    """
    if isinstance(json_complex, str):
        json_complex = json.loads(json_complex)
    try:
        return ComplexField(json_complex['prec'])(json_complex['val'])
    except TypeError:
        return ComplexField(json_complex)


def complex_tuple_to_json(complex_tuple: tuple[ComplexNumber]) -> list[dict]:
    return [complex_number_to_json(s) for s in complex_tuple]


def complex_tuple_from_json(json_list: list[dict] | str) -> tuple[ComplexNumber]:
    """
    Tuple of complex numbers from a json list

    :param json_list:
    :return:
    """
    if isinstance(json_list, str):
        json_list = json.loads(json_list)
    return tuple(complex_number_from_json(s) for s in json_list)


def number_field_to_json(nf: NumberField) -> dict:
    """
    Json representation of number field.

    NOTE: Any information about embeddings is ignored.

    """
    return {'polynomial': str(nf.polynomial()), 'names': list(nf._names)}


def number_field_from_json(data: dict | str) -> NumberField:
    """
    Create number field from json.

    """
    if isinstance(data, str):
        data = json.loads(data)
    return NumberField(ZZ['x'](data['polynomial']), names=tuple(data['names']))


def coefficient_dict_to_json(coeff_dict: dict) -> dict:
    """
    Convert a coefficient dict to JSON format.

    INPUT:

    - ``coeff_dict`` -- dict of coefficients with keys as integer tuples

    EXAMPLES::

        sage: from hilbert_maass.modform.utils import coefficient_dict_to_json
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
    return {json.dumps([int(ki) for ki in k]): complex_number_to_json(v) for k,v in coeff_dict.items() }


def coefficient_dict_from_json(data: dict | str) -> dict:
    """
    Convert a coefficient dict from JSON format.

    INPUT:

    - ``data`` -- json data as dict or string

    EXAMPLES::

        sage: from hilbert_maass.modform.utils import (coefficient_dict_from_json,
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
    return { tuple(json.loads(k)): complex_number_from_json(v) for k,v in data.items() }