r"""
Classes to represent point in the Hyperbolic upper half-space.
In Cython to make it faster if needed.
"""
from typing import ParamSpec, Any, Iterable

from maass_forms_klein.modform.utils import Real_t, Integer_t, Complex_t
from sage.all import I
from sage.functions.other import real, imag
from sage.matrix.constructor import matrix
from sage.modules.free_module import FreeModule_generic
from sage.modules.free_module_element import vector
from sage.rings.complex_mpfr import ComplexField, ComplexField_class
from sage.rings.integer_ring import IntegerRing_class
from sage.rings.number_field.number_field import NumberField_quadratic
from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_gaussian
from sage.rings.rational import Rational
from sage.rings.rational_field import RationalField
from sage.rings.real_mpfr import RealField_class, RealField
from sage.structure.element import Vector
from sage.structure.sequence import Sequence
from sage.symbolic.expression import Expression

P = ParamSpec('P')


class UpperHalfSpace(FreeModule_generic):
    r"""
        Class of upper half-space.
    """

    Element = UpperHalfSpaceElement__class


    def __init__(self, base_ring, degree = 3, sparse=False, category=None):
        # base_ring should be a real field
        if base_ring.base_ring() != base_ring:
            raise ValueError("base_ring should be the ring of integers, rationals or reals")
        super().__init__(base_ring, 3, 3)

    def __repr__(self):
        return 'UpperHalfSpace({})'.format(self.base_ring())

    def _an_element_(self):
        return self.element_class([0, 0, 1])

    def _element_constructor_(self, e: Any, *args: P.args, **kwargs: P.kwargs):
        r"""

        EXAMPLES:

            sage: from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpace
            sage: U = UpperHalfSpace(RR, 3)
            sage: U(vector((CC(1,1),0)))
            1.00000000000000 + 1.00000000000000i + 0.000000000000000j
            sage: U((CC(1,1),0))
            1.00000000000000 + 1.00000000000000i + 0.000000000000000j
            sage: U((1,1))
            1.00000000000000 + 1.00000000000000i + 0.000000000000000j
            sage: U(vector((1,1)))
            1.00000000000000 + 1.00000000000000i + 0.000000000000000j
            sage: U([1,1])
            1.00000000000000 + 1.00000000000000i + 0.000000000000000j
            sage: U = UpperHalfSpace(QQ, 3)
            sage: U((1+I,0))
            1 + 1i + 0j

        """
        if isinstance(e, (Complex_t, Real_t, Integer_t, Rational, Expression,
                                NumberFieldElement_gaussian)):
            e = [real(e), imag(e), 0]
        if isinstance(e, UpperHalfSpaceElement__class):
            return self.element_class(self, [e[0], e[1], e[2]])
        base_ring = Sequence(e).universe()
        if isinstance(base_ring, (ComplexField_class, NumberField_quadratic)):
            base_ring = base_ring.base_ring()
        orig_len = len(e)
        if len(e) == 2 and isinstance(e[0], Complex_t):
            e = [real(e[0]), imag(e[0]), base_ring(e[1])]
        elif len(e) == 2:
            e = [e[0], e[1], base_ring(0)]
        if len(e) != 3 or not all(isinstance(e, (Real_t, Integer_t, Rational, Expression))
                                        for e in e):
            raise ValueError(
                f"Arguments must be two or three real numbers, or one complex and one real: {e} "
                f"{type(e[0]), type(e[1]), type(e[2])} len={len(e)} base={base_ring}"
                f"")
        if e[2] < 0:
            raise ValueError("Argument 3 must be non-negative")
        if not isinstance(e, tuple):
            e = tuple(e)
        return self.element_class(self, e)

    def _coerce_map_from_(self, S):
        """
        Can Coerce to lists, tuples, upper half space elements, integers, rationals,
        real and complex numbers.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpace
            sage: from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement
            sage: UpperHalfSpace(RR, 2).coerce([1,1])
            1.00000000000000 + 1.00000000000000i + 0.000000000000000j
            sage: UpperHalfSpace(RR, 2).coerce([1+I,1])
            1.00000000000000 + 1.00000000000000i + 1.00000000000000j
            sage: z=UpperHalfSpaceElement([1+I*1,0])
            sage: UpperHalfSpace(RR, 2).coerce(z)
            1.00000000000000 + 1.00000000000000i + 0.000000000000000j
            sage: w=UpperHalfSpaceElement([1+I*1,0.0])
            sage: z - w
            0 + 0i + 0j
        """
        if S in [list, tuple, Vector]:
            return True
        if isinstance(S, (UpperHalfSpace, IntegerRing_class, RationalField,
                          RealField_class, ComplexField_class)):
            return True
        return False

    def coerce(self, x):
        r"""
        Coerce x to an element of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProduct
            sage: ComplexPlaneProduct(degree=2).coerce([1,1])
            [1.00000000000000, 1.00000000000000]
            sage: ComplexPlaneProduct(degree=2).coerce([1,1+I])
            [1.00000000000000, 1.00000000000000 + 1.00000000000000*I]

        """
        return self._element_constructor_(x)

cdef class UpperHalfSpaceElement__class(FreeModuleElement_generic_dense):
    r"""
        Class of elements in complex upper half-space.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement
            sage: z =  UpperHalfSpaceElement([1.0,2.0,3.0])
            sage: z
            1.00000000000000 + 2.00000000000000i + 3.00000000000000j
            sage: z.x0()
            1.00000000000000
            sage: z.x1()
            2.00000000000000
            sage: z.y()
            3.00000000000000
            sage: z.z()
            1.00000000000000 + 2.00000000000000*I
            sage: UpperHalfSpaceElement([1+2*I,3])
            1 + 2i + 3j
            sage: UpperHalfSpaceElement([CC(1,2),3])
            1.00000000000000 + 2.00000000000000i + 3.00000000000000j
    """

    # Parent = UpperHalfSpace__class
    def __repr__(self):
        r"""
        Representation of self.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement
            sage: z =  UpperHalfSpaceElement([1+2*I,3])
            sage: z
            1 + 2i + 3j
        """
        return (f"{self[0]} {'+' if self[1] >= 0 else '-'} {abs(self[1])}i "
                f"+ {self[2]}j")

    cpdef x0(self):
        return self[0]

    cpdef x1(self):
        return self[1]

    cpdef y(self):
        return self[2]

    cdef norm(self):
        return self[0] ** 2 + self[1] ** 2 + self[2] ** 2
    cpdef z(self):
        r"""
        Return x0 + i x1 component of self.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement
            sage: z =  UpperHalfSpaceElement([1+2*I,3])
            sage: z.z()
            2*I + 1
            sage: z =  UpperHalfSpaceElement([1.0,2.0,3.0])
            sage: z.z()            
            1.00000000000000 + 2.00000000000000*I
            
        """
        if isinstance(self[0], Real_t):
            return ComplexField(self.base_ring().prec())(self[0], self[1])
        else:
            return self[0] + I * self[1]

    def action(self, A, check=True):
        r"""
        Action of the matrix A in SL(2,C) on self.

        INPUT:
        - ``A`` -- matrix
        - ``check`` -- (default True) check that A is in SL(2,C)

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement
            sage: z =  UpperHalfSpaceElement([1.0, 2.0, 3.0])
            sage: A = matrix([[1, 0], [0, 1]])
            sage: z.action(A)
            1.00000000000000 + 2.00000000000000i + 3.00000000000000j
            sage: z =  UpperHalfSpaceElement([1, 2, 3])
            sage: z.action(matrix([[1, 0], [0, 1]]))
            1 + 2i + 3j
            sage: z.action(matrix([[1, 1], [0, 1]]))
            2 + 2i + 3j
            sage: z.action(matrix([[1, I], [0, 1]]))
            1 + 3i + 3j
            sage: z =  UpperHalfSpaceElement([1, 2, 3])
            sage: z.action(matrix([[0, I], [I, 0]]))
            1/14 - 1/7i + 3/14j
            sage: w = z.action(matrix([[0, -1], [1, 0]])); w
            -1/14 + 1/7i + 3/14j
            sage: w1 = w.action(matrix([[1, 1], [0, 1]])); w1
            13/14 + 1/7i + 3/14j
            sage: w2 = z.action(matrix([[1, -1], [1, 0]])); w2
            13/14 + 1/7i + 3/14j

        """
        A = matrix(A)
        if A.nrows() == 2 and A.ncols() == 2:
            a, b, c, d = A.list()
        else:
            raise ValueError("Input must be coercable to a matrix in SL(2,C).")
        eps = a.parent().epsilon()
        if check and abs(a * d - b * c  - 1) > 16 * eps:
            raise ValueError(f"Matrix must be in SL(2,C): det(A)-1 = "
                             f"{abs(a * d - b * c -1)} > {16 * eps}")
        z = self.z()
        y = self.y()
        denominator = abs(c * z + d) ** 2 + abs(c * y) ** 2
        numerator = (a * z + b) * (c * z + d).conjugate() + a * c.conjugate() * y ** 2
        x0 = real(numerator) / denominator
        x1 = imag(numerator) / denominator
        y = y / denominator
        return UpperHalfSpaceElement([x0, x1, y])

    def __add__(self, other: Any) -> UpperHalfSpaceElement__class:
        """
        Add other to self.

        EXAMPLES:

            sage: from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement
            sage: z =  UpperHalfSpaceElement([1,2,3])
            sage: z + z
            2 + 4i + 6j
            sage: z + [1,2,3]
            2 + 4i + 6j
            sage: z + vector([1,2,3])
            2 + 4i + 6j
            sage: z + vector([1.0,2.0])
            2.00000000000000 + 4.00000000000000i + 3.00000000000000j

        """
        if not isinstance(other, UpperHalfSpaceElement__class):
            other = UpperHalfSpaceElement(other)
        return UpperHalfSpaceElement([self[i]+other[i] for i in range(3)])

    def __sub__(self, other: Any) -> UpperHalfSpaceElement__class:
        """
            Subtract other from self.

            EXAMPLES:

                sage: from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement
                sage: z =  UpperHalfSpaceElement([1,2,3])
                sage: z - z
                0 + 0i + 0j
                sage: z - [1,2,3]
                0 + 0i + 0j
                sage: z - vector([1,2,3])
                0 + 0i + 0j
                sage: z - vector([1.0,2.0])
                0.000000000000000 + 0.000000000000000i + 3.00000000000000j

            """
        if not isinstance(other, UpperHalfSpaceElement__class):
            try:
                other = UpperHalfSpaceElement(other)
            except ValueError as e:
                raise e
                raise ValueError(f"Cannot subtract {other} from {self}. Type={type(other)}, Equal={other == vector((0.0,0.0))}")
        return UpperHalfSpaceElement([self[i] - other[i] for i in range(3)])

    def translate(self, z: Complex_t | Real_t):
        """
        Translate self by complex number z.
        
        INPUT:
        - ``z`` -- complex number
        
        EXAMPLES:
        
            sage: from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement
            sage: z =  UpperHalfSpaceElement([1.0,2.0,3.0])
            sage: z
            1.00000000000000 + 2.00000000000000i + 3.00000000000000j
            sage: z.translate(1+2*I)
            2.00000000000000 + 4.00000000000000i + 3.00000000000000j
            sage: z.translate(1.0)
            2.00000000000000 + 2.00000000000000i + 3.00000000000000j
        """
        elements = [self[0] + real(z), self[1] + imag(z), self[2]]
        return UpperHalfSpaceElement(elements)

    # cpdef add(self, UpperHalfSpaceElement__class other):
    #     r"""
    #     Add other to self.
    #
    #     INPUT:
    #     - ``other`` -- UpperHalfSpaceElement
    #
    #     EXAMPLES:
    #
    #         sage: from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement
    #         sage: z =  UpperHalfSpaceElement([1,2,3])
    #         sage: z
    #         1.00000000000000+2.00000000000000i+3.00000000000000j
    #         sage: z.add(z)
    #         2.00000000000000+4.00000000000000i+6.00000000000000j
    #         sage: z.add([1,2,3)
    #         2.00000000000000+4.00000000000000i+6.00000000000000j
    #         sage: z.add([1,2)
    #
    #     """
    #     if isinstance(other, (UpperHalfSpaceElement__class, list)):
    #         return UpperHalfSpaceElement__class([self[0] + other[0], self[1] + other[1],
    #                                              self[1] + other[2]])
    #     raise NotImplemented(f"Cannot add {other} to UpperHalfSpaceElement")
def UpperHalfSpaceElement(entries: Iterable | Any) -> UpperHalfSpaceElement__class:
    r"""
    Create upper half space element.

    EXAMPLES:
        sage: from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement
        sage: UpperHalfSpaceElement([1,2,3])
        1 + 2i + 3j
        sage: UpperHalfSpaceElement([1.0,2.0,3.0])
        1.00000000000000 + 2.00000000000000i + 3.00000000000000j
        sage: UpperHalfSpaceElement(vector([1.0,2.0,3.0]))
        1.00000000000000 + 2.00000000000000i + 3.00000000000000j
        sage: UpperHalfSpaceElement([1.0+2.0*I,3.0])
        1.00000000000000 + 2.00000000000000i + 3.00000000000000j
        sage: UpperHalfSpaceElement(1.0+2.0*I)
        1.00000000000000 + 2.00000000000000i + 0.000000000000000j
        sage: UpperHalfSpaceElement((1.0, 2.0))
        1.00000000000000 + 2.00000000000000i + 0.000000000000000j
        sage: UpperHalfSpaceElement(1.0)
        1.00000000000000 + 0.000000000000000i + 0.000000000000000j
    """
    if isinstance(entries, (UpperHalfSpaceElement__class, list, tuple, Vector)):
        x = entries[0]
        if hasattr(x, 'base_ring'):
            base_ring = x.base_ring()
        else:
            base_ring = RealField(53)
    else:
        base_ring = entries.base_ring()
    return UpperHalfSpace(base_ring, 3)(entries)