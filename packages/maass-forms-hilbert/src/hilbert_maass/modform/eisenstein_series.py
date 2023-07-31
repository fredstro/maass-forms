from hilbert_maass.modform.coefficients import HilbertMaassCoefficients
from hilbert_maass.utils import map_int_to_tuple
from sage.arith.misc import divisors
from sage.categories.sets_cat import cartesian_product
from sage.rings.complex_mpfr import ComplexNumber
from sage.rings.integer import Integer
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.number_field.number_field_element import NumberFieldElement
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.structure.sage_object import SageObject
from typing import NoReturn, ParamSpec

from sage.structure.element import ModuleElement

from ..utils import Integer_t

P = ParamSpec('P')


class HilbertEisensteinSeries(SageObject):
    """
    Class for non-holomorphic Hilbert Eisenstein series.

    """
    def __init__(self, number_field, ideal, spectral_parameter, **kwargs: P.kwargs) -> NoReturn:
        r"""

        INPUT:

        - ``number_field`` -- number field
        - ``ideal`` -- ideal
        - ``spectral_parameter`` -- tuple of complex numbers
        """
        super(HilbertEisensteinSeries, self).__init__()
        self.has_coefficients = False
        self._number_field = number_field
        self._ideal = ideal
        self._dual_ideal = ideal ** -1 * number_field.different() ** -1
        self._dual_basis = self._dual_ideal.integral_basis()
        self._spectral_parameter = spectral_parameter
        self._complex_field = spectral_parameter[0].parent()
        if self._number_field != QuadraticField(2):
            raise NotImplementedError(f'not implemented for {self.base_ring}')

    def number_field(self):
        return self._number_field

    def CF(self):
        return self._complex_field

    def _dual_ideal_element(self, v):
        if isinstance(v, NumberFieldElement) and v in self._dual_ideal:
            return v
        if isinstance(v, tuple):
            return sum(self._dual_basis[i] * v[i] for i in range(len(v)))
        raise ValueError(f'Dual element not implemented for {v}')

    def av(self, v: tuple) -> ComplexNumber:
        if all(vi == 0 for vi in v):
            return self.a0_minus()
        v = self._dual_ideal_element(v)
        half = self.CF()(1) / self.CF()(2)
        s = self._spectral_parameter[0]
        return 1 / self._number_field.discriminant().sqrt() * \
            (2 * self.CF().pi() ** s / s.gamma()) ** 2 * \
            abs(v.norm()) ** (s - half) * \
            self._sigma_K(v, 1 - 2 * s) / \
            self._number_field.zeta_function(self.CF().prec())(2 * s)

    def a0_minus(self):
        r"""
        The 0-th Fourier coefficient of this Eisenstein Series

        :return:
        """
        s = self._spectral_parameter[0]
        CF = self._complex_field
        s_minus_half = s - CF(1) / CF(2)
        result = CF.pi() / self.number_field().discriminant().sqrt() * (
                    s_minus_half.gamma() / s.gamma()) ** 2 * \
                    self.number_field().zeta_function(CF.prec())(
            2 * s - 1) / self.number_field().zeta_function(CF.prec())(2 * s)
        return CF(result)

    def _ideals_dividing_vD(self, v: tuple | NumberFieldElement) -> \
            list[NumberFieldFractionalIdeal]:
        """
        Returns divisors of the ideal v*[different ideal] for v in self._dual_ideal

        INPUT::

        - ``ideal`` -- ideal in number field
        - ``v`` -- element in number field

        """
        v = self._dual_ideal_element(v)
        idealv = self.number_field().fractional_ideal(v*self._ideal.number_field().different())
        return divisors(idealv)

    def _sigma_K(self, v: tuple | NumberFieldElement, s: ComplexNumber) -> ComplexNumber:
        r"""
        Divisor function in number field $\sum_{d | (v)*ideal} N(d)**(-s)$

        INPUT:

        - ``v`` -- element in number field
        - ``s`` -- complex number

        EXAMPLES::

            sage:

        """
        return sum(ideal_c.norm() ** s for ideal_c in self._ideals_dividing_vD(v))

    def coefficients(self, M: tuple[tuple[Integer_t]]) -> HilbertMaassCoefficients:
        different = self.number_field().different()
        dual_ideals = [ideal ** -1 * different ** -1 for ideal in self._dual_ideals]
        for V in cartesian_product([range(-m0[0], m0[1] + 1) for m0 in M]):
            for ideala_dual in dual_ideals:
                v = self.dual_ideal_element(V, self._ideal)

        C = HilbertMaassCoefficients(X, M, self._dual_ideals)