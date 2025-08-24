from sage.modules.free_module_element cimport FreeModuleElement_generic_dense
from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.rings.complex_mpc cimport MPComplexField_class, MPComplexNumber
from sage.rings.real_mpfr cimport RealNumber, RealField_class
from sage.categories.map cimport Map


cdef class UpperHalfSpaceElement__class(FreeModuleElement_generic_dense):
    """
    TODO: think about the best representation.
    """
    # Double precision representation
    # cdef MPComplexNumber _z
    # cdef int _prec
    # cdef list _vector
    # cdef RealField_class _base_ring
    cpdef x0(self)
    cpdef x1(self)
    cpdef y(self)
    cpdef z(self)
    cdef norm(self)
    # cpdef action(self, A)
    # cpdef add(self, UpperHalfSpaceElement__class other)