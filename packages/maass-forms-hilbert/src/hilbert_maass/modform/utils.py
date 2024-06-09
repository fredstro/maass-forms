import json
from typing import Iterable
import logging

from hilbert_modgroup.pullback import HilbertPullback
from sage.all import ZZ, CC
from sage.arith.misc import factor
from sage.categories.sets_cat import cartesian_product
from sage.functions.other import ceil
from sage.misc.functional import round
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
from sage.rings.real_mpfr import RealNumber as RealNumber_class
from sage.structure.element import Matrix, Vector
from sage.rings.number_field.unit_group import UnitGroup

try:
    from comp_manager.decorators import mongo_cache
except ModuleNotFoundError:
    # Handle the case when comp_manager is not installed
    def mongo_cache(*args, **kwargs):
        def decorator(func):
            return func
        return decorator

log = logging.getLogger(__name__)
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
    C = C * (max(max(abs(b0), abs(b1)) for b0, b1 in M)+1)
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

    INPUT:
     - `ideala` -- an ideal
     - `element` -- an element of `ideala`
     - `check` -- if True, check that the coordinates are correct

    EXAMPLES::

        sage: from hilbert_maass.modform.utils import (dual_ideal, ideal_generator,
        ....:   ideal_coordinates, number_field_basis_matrix, ideal_basis_matrix)
        sage: ideala_dual = dual_ideal(QuadraticField(3).ideal(1))
        sage: delta = ideal_generator(ideala_dual)
        sage: ideal_coordinates(ideala_dual, delta, check=True)
        (0, -1)
        sage: ideala_dual = dual_ideal(QuadraticField(5).ideal(1))
        sage: delta = ideal_generator(ideala_dual)
        sage: ideal_coordinates(ideala_dual, delta, check=True)
        (1, -1)

    """
    nf = ideala.number_field()
    if not isinstance(element, NumberFieldElement):
        element = nf(element)
    if element not in ideala:
        raise ValueError(f"Element {element} not in ideal: {ideala}")
    basis_change_matrix = ideal_basis_matrix(ideala) ** -1 * number_field_basis_matrix(nf)
    coordinates = basis_change_matrix * element.vector()
    eps = basis_change_matrix.base_ring().epsilon() * 2 ** 3
    coordinates_int = tuple(round(c.real()) for c in coordinates if abs(c.imag()) < eps)
    if len(coordinates_int) != len(coordinates):
        raise ArithmeticError(f"Can not find lattice coordinates for delta={element}."
                              f" coordinates={coordinates}"
                              f" coordinates_int={coordinates_int}")
    if check:
        assert sum([c * ideala.integral_basis()[i]
                    for i, c in enumerate(coordinates_int)]) == element
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

def ideal_generator(ideala: NumberFieldFractionalIdeal) -> NumberFieldFractionalIdeal:
    """
    Find a totally positive generator for an ideal if possible, else use a reduced generator.

    INPUT:

    - ``ideala`` -- NumberFieldFractionalIdeal


    """
    narrow_class_number=ideala.number_field().narrow_class_group().order()
    if (narrow_class_number==1):
        u = UnitGroup(space.number_field()).gens_values()[1]
        x = ideala.gens_reduced()[0]
        #f = space.number_field().galois_group()
        test =[x,  x*u, -x*u]
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
        delta= ideala.gens_reduced()[0]
    return delta


def complex_number_to_json(s: ComplexNumber) -> dict:
    """
    Create json data from a ComplexNumber.

    INPUT:

    - ``s`` -- complex number
    """
    if is_json_number(s):
        return s
    if not isinstance(s, ComplexNumber):
        s = CC(s)
    return {'prec': s.parent().prec(), 'val': str(s)}


def is_json_number(data: dict | str) -> bool:
    return isinstance(data, dict) and 'prec' in data and 'val' in data and len(dict.keys())==2

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
    if all(isinstance(x, str) and is_json_number(x) for x in coeff_dict.items()):
        return coeff_dict
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


def integer_to_bounds_tuple(m: Integer_t, degree: Integer_t) -> tuple[tuple[Integer_t]]:
    """
    Convert a positive integer to a tuple of bounds.

    INPUT:

    -``m`` -- positive integer
    -``degree`` -- positive integer

    EXAMPLES::

        sage: from hilbert_maass.modform.utils import integer_to_bounds_tuple
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


def unit_relations(m: Integer_t, space: 'HilbertMaassFormSpace'):
    """
        produce the sets of integer lattice points which are related by the
        automorphy a(unit^2 v)=a(v) for the quadratic real fields.

        INPUT:

        -``m`` -- positive integer
        -``space`` -- HilbertMaassFormSpace

        EXAMPLES::

            sage: from hilbert_maass.modform.utils import unit_relations
            sage: space = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
            sage: unit_relations(6, space)
            [[(-6, -6), (-6, 6)],
             [(-5, -6), (-3, 2)],
             [(-5, -5), (-5, 5)],
             [(-5, 6), (-3, -2)],
            [(-4, -6), (0, -2), (4, -6)],......]
        """
    ideala = space.number_field().ideal(1)
    dual_ideala = dual_ideal(ideala)
    store1 = cartesian_product_from_M(((-m, m), (-m, m)))
    store = []
    for r in store1:
        if (abs(r[0]) <=m and abs(r[1]) <=m):
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
        while (abs(d[0])<=m and abs(d[1])<=m):
            kemp.append(tuple(d))
            if (d in store):
                store.remove(tuple(d))
            use = use * unit
            d = ideal_coordinates(dual_ideala, use)
            x = x + 1
        use = r_element * unit ** (-1)
        d = ideal_coordinates(dual_ideala, use)
        while abs(d[0]) <=m and abs(d[1])<= m:
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


def hecke_relations_coprime(m: Integer_t, space: 'HilbertMaassFormSpace'):
    """
            produce the sets of integer lattice points which are related by the
            the Hecke relation of the form a(\delta m)a(\delta n)=a(\delta mn), where \delta is a generator of
            dual ideal

            INPUT:

            -``m`` -- positive integer
            -``space`` -- HilbertMaassFormSpace

            EXAMPLES::

                sage: from hilbert_maass.modform.utils import hecke_relations_coprime
                sage: space = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
                sage: hecke_relations_coprime(6, space)
                [[(-6, 6), (-1, -1), (6, 6)],
                 [(-6, 6), (-1, 0), (6, 0)],
                 [(-6, 6), (1, 0), (-6, 0)],
                 [(-6, 6), (1, 1), (-6, -6)],......]
            """
    ideala = space.number_field().ideal(1)
    dual_ideala = dual_ideal(ideala)
    t = ideal_generator(dual_ideala)
    store1 = cartesian_product_from_M(((-m, m), (-m, m)))
    store = []
    for r in store1:
        if abs(r[0])<=m and abs(r[1])<=m:
            store.append(r)
    store.remove((0, 0))
    temp = []
    for r in store:
        for s in store:
            r_element = (dual_ideal_element(r, ideala, as_nf_element=True) / t)
            r_ideal = space.number_field().ideal(r_element)
            s_element = (dual_ideal_element(s, ideala, as_nf_element=True) / t)
            s_ideal = space.number_field().ideal(s_element)
            if r != s and r_ideal != ideala and s_ideal != ideala and r_ideal + s_ideal == ideala:
                y = t * r_element * s_element
                d = ideal_coordinates(dual_ideala, y)
                if (d[0], d[1]) in store:
                    temp.append([r, s, d])
    return temp



def hecke_relations_prime_power(m: Integer_t, space: 'HilbertMaassFormSpace'):
    """
        produce the sets of integer lattice points which are related by the
        the Hecke relation of the form a(\delta p^n)=a(\delta p)a(\delta p^{n-1})-a(\delta p^{n-2}), where \delta is a generator of
        dual ideal

        INPUT:

        -``m`` -- positive integer
        -``space`` -- HilbertMaassFormSpace

        EXAMPLES::

            sage: from hilbert_maass.modform.utils import hecke_relations_prime_power
            sage: space = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
            sage: hecke_relations_prime_power(6, space)
            [[(1, (-3, 5)), (2, (5, -5))],
            [(1, (-2, 1)), (2, (3, 4))],
            [(1, (-2, 2)), (2, (4, -4))],
            [(1, (-2, 4)), (2, (0, 4))],
            [(1, (-1, 0)), (2, (0, 5))],
            [(1, (1, 0)), (2, (0, 5))],..........]
        """
    ideala = space.number_field().ideal(1)
    dual_ideala = dual_ideal(ideala)
    t = ideal_generator(dual_ideala)
    f = space.number_field().galois_group()
    store1 = cartesian_product_from_M(((-m, m), (-m, m)))
    store = []
    for r in store1:
        if (abs(r[0])<=6 and abs(r[1])<=6):
            store.append(r)
    store.remove((0, 0))
    temp = []
    for r in store:
        r_element = (dual_ideal_element(r, ideala, as_nf_element=True) / t)
        r_ideal = space.number_field().ideal(r_element)
        if r_ideal.is_prime():
            a = r_element
            x = 1
            d = r
            kemp = []
            while (d[0], d[1]) in store:
                kemp.append((x, d))
                a = a * r_element
                y = a * t
                d = ideal_coordinates(dual_ideala, y)
                x = x + 1
            if (x <= 2):
                kemp.pop()
            else:
                temp.append(kemp)
    return temp


def symmetric_relations(space: 'HilbertMaassFormSpace'):
    """
            produces a set of integer lattice points which are related by
            reflection relation (flipping operators) a(u_p v)=ta(v), where u is a unit which is neagtive
            at the place (prime) p and positive at other places.

            INPUT:

            -``space`` -- HilbertMaassFormSpace

            EXAMPLES::

                sage: from hilbert_maass.modform.utils import symmetric_relations
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
    set_check = [t, t * u ** -1, -t * u ** -1]
    kemp = []
    for use in set_check:
        d = ideal_coordinates(dual_ideala, use)
        kemp.append(d)
    return (kemp)

def ideal_factors(ida):
    prime_factors = [x[0] for x in factor(ida)]
    exponents = [range(x[1]+1) for x in factor(ida)]
    new_exponents = list(cartesian_product(exponents))
    factors = []
    for ex in new_exponents:
        idb = prod([p**ex[i] for i, p in enumerate(prime_factors)])
        factors.append(idb)
    return factors
