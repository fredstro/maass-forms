"""
Number field basis and ideal utilities.

Provides functions for computing basis matrices, dual ideals, and ideal
coordinates used across all Maass form packages that work with number fields.

EXAMPLES::

    sage: from maass_form_core.functions.basis import number_field_basis_matrix
    sage: K = QuadraticField(5)
    sage: B = number_field_basis_matrix(K)
    sage: a = K.gen()
    sage: B * a.vector() == vector(a.complex_embeddings())
    True
"""

from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_function
from sage.modules.free_module_element import vector
from sage.rings.number_field.number_field_element import NumberFieldElement
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.rings.number_field.unit_group import UnitGroup
from sage.structure.element import Matrix


@cached_function()
def number_field_basis_matrix(number_field, prec: int = 53) -> Matrix:
    r"""
    Return the basis matrix with respect to the '.vector' property of number field elements.

    Computes a matrix B such that for any number field element e,
    ``B * e.vector() == e.complex_embeddings(prec)``.

    INPUT:

    - ``number_field`` -- a NumberField
    - ``prec`` -- integer (default: 53); precision in bits for complex embeddings

    OUTPUT:

    - A matrix whose rows are the complex embeddings of the number field basis elements

    EXAMPLES::

        sage: from maass_form_core.functions.basis import number_field_basis_matrix
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
    Return the matrix whose columns contain the complex embeddings of the
    ideal's integral basis elements.

    INPUT:

    - ``ideal`` -- NumberFieldFractionalIdeal
    - ``prec`` -- integer (default: 53); precision in bits

    OUTPUT:

    - A matrix of complex embeddings

    EXAMPLES::

        sage: from maass_form_core.functions.basis import ideal_basis_matrix
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
    Return the dual ideal (codifferent) of a given fractional ideal.

    Computes `\mathfrak{a}^{-1}\mathcal{D}^{-1}` where `\mathcal{D}`
    is the different ideal of the number field.

    INPUT:

    - ``ideala`` -- NumberFieldFractionalIdeal

    OUTPUT:

    - NumberFieldFractionalIdeal; the dual ideal

    EXAMPLES::

        sage: from maass_form_core.functions.basis import dual_ideal
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

    INPUT:

    - ``ideal`` -- NumberFieldFractionalIdeal
    - ``prec`` -- integer (default: 53); precision in bits

    OUTPUT:

    - A matrix whose rows are the complex embeddings of the dual ideal's
      integral basis elements

    EXAMPLES::

        sage: from maass_form_core.functions.basis import dual_ideal_basis_matrix
        sage: K = QuadraticField(5)
        sage: I = K.ideal(1)
        sage: B = dual_ideal_basis_matrix(I)
        sage: B
        [ 1.00000000000000 0.276393202250021]
        [ 1.00000000000000 0.723606797749979]

    The dual ideal basis matrix times the ideal basis matrix should give
    a matrix with determinant close to 1::

        sage: from maass_form_core.functions.basis import ideal_basis_matrix
        sage: D = dual_ideal_basis_matrix(I)
        sage: M = ideal_basis_matrix(I)
        sage: abs((D * M).det() - 1) < 1e-10
        True
    """
    dual = ideal**-1 * ideal.number_field().different() ** -1
    return ideal_basis_matrix(dual, prec)


def ideal_coordinates(ideala: NumberFieldFractionalIdeal, element, check: bool = False) -> tuple:
    r"""
    Find the coordinates of an element with respect to an integral basis of an ideal.

    INPUT:

    - ``ideala`` -- NumberFieldFractionalIdeal
    - ``element`` -- a number field element in the ideal
    - ``check`` -- boolean (default: False); verify coordinates are correct

    OUTPUT:

    - Tuple of integers giving the coordinates

    EXAMPLES::

        sage: from maass_form_core.functions.basis import (
        ....:     dual_ideal, ideal_generator, ideal_coordinates)
        sage: ideala_dual = dual_ideal(QuadraticField(3).ideal(1))
        sage: delta = ideal_generator(ideala_dual)
        sage: coords = ideal_coordinates(ideala_dual, delta)
        sage: coords
        (0, 1)
        sage: sum([c * ideala_dual.integral_basis()[i] for i, c in enumerate(coords)]) == delta
        True
        sage: ideala_dual = dual_ideal(QuadraticField(5).ideal(1))
        sage: delta = ideal_generator(ideala_dual)
        sage: ideal_coordinates(ideala_dual, delta)
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
    from sage.misc.functional import round

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
def dual_ideal_element(coordinates, ideal, as_nf_element=False):
    r"""
    Create an element in the dual ideal from the given coordinates.

    INPUT:

    - ``coordinates`` -- tuple of integers or vector
    - ``ideal`` -- NumberFieldFractionalIdeal
    - ``as_nf_element`` -- boolean (default: False); return as number field element

    OUTPUT:

    - Vector of complex embeddings, or a number field element if as_nf_element=True

    EXAMPLES::

        sage: from maass_form_core.functions.basis import dual_ideal_element, dual_ideal
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


def ideal_generator(ideala: NumberFieldFractionalIdeal):
    r"""
    Find a totally positive generator for an ideal if possible.

    For number fields with narrow class number 1, searches for a totally
    positive generator. Otherwise uses a reduced generator.

    INPUT:

    - ``ideala`` -- NumberFieldFractionalIdeal

    OUTPUT:

    - A generator for the ideal, preferably totally positive

    EXAMPLES::

        sage: from maass_form_core.functions.basis import ideal_generator
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
        test = [x, x * u, -x * u]
        delta = None
        for delta_test in test:
            if not delta_test.is_totally_positive():
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
