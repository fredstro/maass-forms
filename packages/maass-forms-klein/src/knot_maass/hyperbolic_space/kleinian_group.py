r"""
Class representing a Kleinian group, i.e. a discrete subgroup of PSL(2,C)
acting discretely on the upper half-space.

"""
from knot_maass.modform.upper_half_space.pyx import UpperHalfSpaceElement__class
from sage.all import ZZ
from sage.groups.matrix_gps.linear import LinearMatrixGroup_generic
from sage.rings.integer import Integer
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.number_field.order import Order


class KleinianGroup_class(LinearMatrixGroup_generic):
    r"""
    Class representing a Kleinian group, i.e. a discrete subgroup of PSL(2,C)
    acting discretely on the upper half-space.

    TODO: write examples.
    """
    def __init__(self, x, **kwargs):
        if isinstance(x, (int, Integer)) and x < 0 and ZZ(x).is_fundamental_discriminant():
            x = QuadraticField(x).ring_of_integers()
        if isinstance(x, Order):
            LinearMatrixGroup_generic.__init__(self, x, **kwargs)
        else:
            raise NotImplementedError

    def pullback(self, z: UpperHalfSpaceElement__class):
        r"""
        Pull back a point in the upper half-space to an element of the fundamental
        domain.
        """
        if self.base_ring().discriminant() == -4:
            return self._pullback_gaussian_integers(z)
        raise NotImplementedError

    def _pullback_gaussian_integers(self, z: UpperHalfSpaceElement__class):
        r"""
        Special case of Gaussian integers

        TODO: Implement this and write examples.
        NOTE: Could be made faster by e.g. calling a Cython function.
        """
        raise NotImplementedError
