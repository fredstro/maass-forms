r"""
Class representing a Kleinian group, i.e. a discrete subgroup of PSL(2,C)
acting discretely on the upper half-space.

"""
from knot_maass.modform.upper_half_space.pyx import UpperHalfSpaceElement__class
from sage.all import ZZ
from sage.groups.matrix_gps.matrix_group import MatrixGroup_generic
from sage.rings.integer import Integer
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.number_field.order import Order


class KleinianGroup_class(MatrixGroup_generic):
    r"""
    Class representing a Kleinian group, i.e. a discrete subgroup of PSL(2,C)
    acting discretely on the upper half-space.

    TODO: write examples.

    """
    def __init__(self, x, **kwargs):
        r"""

        INPUT:
        - ``x`` -- integer, order in imaginary quadratic number field, list of generators

        """
        MS = MatrixSpace(QQbar, 2, 2)
        if isinstance(x, (int, Integer)) and x < 0 and ZZ(x).is_fundamental_discriminant():
            K = QuadraticField(x)
            if K.class_number() > 1:
                raise NotImplementedError("Only class number 1 is supported.")
            x = QuadraticField(x).ring_of_integers().gens()
            for b in K.basis():
                gens.append()
        if isinstance(x, Order):
            LinearMatrixGroup_generic.__init__(self, x, **kwargs)
        elif isinstance(x, list):
            # x is a list of generators

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
