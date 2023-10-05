r"""
Classes to represent point in the Hyperbolic upper half-space.
In Cython to make it faster if needed.
"""
from typing import ParamSpec
from sage.structure.element cimport Element
from sage.rings.all import Integer, CC
from sage.rings.infinity import Infinity
from sage.structure.element import Matrix
from sage.structure.parent cimport Parent

P = ParamSpec('P')

cdef class UpperHalfSpace__class(Parent):

    Element = UpperHalfSpaceElement__class

    def __init__(self):
        pass

    def _an_element_(self):
        return UpperHalfSpaceElement__class([0, 0, 1])

cdef class UpperHalfSpaceElement__class(Element):
    r"""
        Class of elements in complex upper half-space.
        EXAMPLES::

            sage: from knot_maass.hyperbolic_space.upper_half_space import UpperHalfSpaceElement__class
            sage: z =  UpperHalfSpaceElement__class([1,2,3])
            sage: z
            1.0+2.0i+3.0j
            sage: z.x0()
            1.0
            sage: z.x1()
            2.0
            sage: z.y()
            3.0
    """

    Parent = UpperHalfSpace__class

    def __init__(self, x, *args: P.args, **kwargs: P.kwargs):
        if isinstance(x, list) and len(x) == 3:
            if x[2] <= 0:
                raise ValueError("Value of 'y' most be positive.")
            self._x0 = <double>x[0]
            self._x1 = <double>x[1]
            self._y = <double>x[2]
        else:
            raise NotImplemented()

    def __repr__(self):
        return f"{self._x0}+{self._x1}i+{self._y}j"

    cpdef double x0(self):
        return self._x0

    cpdef double x1(self):
        return self._x1

    cpdef double y(self):
        return self._y

    cpdef action(self, A):
        r"""
        TODO: Implement this
        """
        raise NotImplemented