from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.rings.complex_mpc cimport MPComplexField_class, MPComplexNumber
from sage.rings.real_mpfr cimport RealNumber
from sage.categories.map cimport Map

cdef class UpperHalfSpace__class(Parent):
    pass

cdef class UpperHalfSpaceElement__class(Element):
    """
    TODO: think about best representation.
    """
    # Double precision representation
    cdef double _x0
    cdef double _x1
    cdef double _y

    cpdef double x0(self)
    cpdef double x1(self)
    cpdef double y(self)
    cpdef action(self, A)