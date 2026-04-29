from typing import NoReturn, ParamSpec

from sage.arith.misc import divisors
from sage.rings.complex_mpfr import ComplexNumber
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.number_field.number_field_element import NumberFieldElement
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.structure.sage_object import SageObject

from maass_forms_hilbert.modform.coefficients import HilbertMaassCoefficients

from .utils import Integer_t

P = ParamSpec("P")


class HilbertEisensteinSeries(SageObject):
    """
    Class for non-holomorphic Hilbert Eisenstein series.

    """

    def __init__(self, number_field, ideal, spectral_parameter, **kwargs: P.kwargs) -> NoReturn:
        r"""
        Initialize the Hilbert Eisenstein series.

        INPUT:

        - ``number_field`` -- number field
        - ``ideal`` -- ideal
        - ``spectral_parameter`` -- tuple of complex numbers

        EXAMPLES::

            sage: from maass_forms_hilbert.modform.eisenstein_series import HilbertEisensteinSeries
            sage: K = QuadraticField(2)
            sage: ideal = K.ideal(1)
            sage: spectral_parameter = (CC(0.5,1), CC(0.5,1))
            sage: E = HilbertEisensteinSeries(K, ideal, spectral_parameter)
            sage: E.has_coefficients
            False
        """
        super().__init__()
        self.has_coefficients = False
        self._number_field = number_field
        self._ideal = ideal
        self._dual_ideal = ideal**-1 * number_field.different() ** -1
        self._dual_basis = self._dual_ideal.integral_basis()
        self._spectral_parameter = spectral_parameter
        self._complex_field = spectral_parameter[0].parent()
        if self._number_field != QuadraticField(2):
            raise NotImplementedError(f"not implemented for {self.base_ring}")

    def number_field(self):
        r"""
        Return the number field associated to this Eisenstein series.

        OUTPUT:

        - Number field

        EXAMPLES::

            sage: from maass_forms_hilbert.modform.eisenstein_series import HilbertEisensteinSeries
            sage: K = QuadraticField(2)
            sage: ideal = K.ideal(1)
            sage: spectral_parameter = (CC(0.5,1), CC(0.5,1))
            sage: E = HilbertEisensteinSeries(K, ideal, spectral_parameter)
            sage: E.number_field() == K
            True
        """
        return self._number_field

    def CF(self):
        r"""
        Return the complex field used for computations.

        OUTPUT:

        - Complex field

        EXAMPLES::

            sage: from maass_forms_hilbert.modform.eisenstein_series import HilbertEisensteinSeries
            sage: K = QuadraticField(2)
            sage: ideal = K.ideal(1)
            sage: spectral_parameter = (CC(0.5,1), CC(0.5,1))
            sage: E = HilbertEisensteinSeries(K, ideal, spectral_parameter)
            sage: E.CF().precision()
            53
        """
        return self._complex_field

    def _dual_ideal_element(self, v):
        r"""
        Return the dual ideal element for v.

        INPUT:

        - ``v`` -- tuple or NumberFieldElement

        OUTPUT:

        - NumberFieldElement in the dual ideal

        EXAMPLES::

            sage: from maass_forms_hilbert.modform.eisenstein_series import HilbertEisensteinSeries
            sage: K = QuadraticField(2)
            sage: ideal = K.ideal(1)
            sage: spectral_parameter = (CC(0.5,1), CC(0.5,1))
            sage: E = HilbertEisensteinSeries(K, ideal, spectral_parameter)
            sage: result = E._dual_ideal_element((1, 0))
            sage: result in K
            True
        """
        if isinstance(v, NumberFieldElement) and v in self._dual_ideal:
            return v
        if isinstance(v, tuple):
            return sum(self._dual_basis[i] * v[i] for i in range(len(v)))
        raise ValueError(f"Dual element not implemented for {v}")

    def av(self, v: tuple) -> ComplexNumber:
        r"""
        Return the v-th Fourier coefficient of this Eisenstein series.

        INPUT:

        - ``v`` -- tuple

        OUTPUT:

        - ComplexNumber, the v-th Fourier coefficient

        EXAMPLES::

            sage: from maass_forms_hilbert.modform.eisenstein_series import HilbertEisensteinSeries
            sage: K = QuadraticField(2)
            sage: ideal = K.ideal(1)
            sage: spectral_parameter = (CC(0.5,1), CC(0.5,1))
            sage: E = HilbertEisensteinSeries(K, ideal, spectral_parameter)
            sage: result = E.av((1, 0))
            sage: result is not None
            True
        """
        if all(vi == 0 for vi in v):
            return self.a0_minus()
        v = self._dual_ideal_element(v)
        half = self.CF()(1) / self.CF()(2)
        s = self._spectral_parameter[0]
        return (
            1
            / self._number_field.discriminant().sqrt()
            * (2 * self.CF().pi() ** s / s.gamma()) ** 2
            * abs(v.norm()) ** (s - half)
            * self._sigma_K(v, 1 - 2 * s)
            / self._number_field.zeta_function(self.CF().prec())(2 * s)
        )

    def a0_minus(self):
        r"""
        The 0-th Fourier coefficient of this Eisenstein Series.

        OUTPUT:

        - ComplexNumber, the constant term

        EXAMPLES::

            sage: from maass_forms_hilbert.modform.eisenstein_series import HilbertEisensteinSeries
            sage: K = QuadraticField(2)
            sage: ideal = K.ideal(1)
            sage: spectral_parameter = (CC(0.5,1), CC(0.5,1))
            sage: E = HilbertEisensteinSeries(K, ideal, spectral_parameter)
            sage: result = E.a0_minus()
            sage: result.parent() == CC
            True
        """
        s = self._spectral_parameter[0]
        CF = self._complex_field
        s_minus_half = s - CF(1) / CF(2)
        result = (
            CF.pi()
            / self.number_field().discriminant().sqrt()
            * (s_minus_half.gamma() / s.gamma()) ** 2
            * self.number_field().zeta_function(CF.prec())(2 * s - 1)
            / self.number_field().zeta_function(CF.prec())(2 * s)
        )
        return CF(result)

    def _ideals_dividing_vD(
        self, v: tuple | NumberFieldElement
    ) -> list[NumberFieldFractionalIdeal]:
        r"""
        Returns divisors of the ideal v*[different ideal] for v in self._dual_ideal.

        INPUT:

        - ``v`` -- element in number field

        OUTPUT:

        - list of NumberFieldFractionalIdeal

        EXAMPLES::

            sage: from maass_forms_hilbert.modform.eisenstein_series import HilbertEisensteinSeries
            sage: K = QuadraticField(2)
            sage: ideal = K.ideal(1)
            sage: spectral_parameter = (CC(0.5,1), CC(0.5,1))
            sage: E = HilbertEisensteinSeries(K, ideal, spectral_parameter)
            sage: result = E._ideals_dividing_vD((1, 0))
            sage: len(result) > 0
            True
        """
        v = self._dual_ideal_element(v)
        idealv = self.number_field().fractional_ideal(v * self._ideal.number_field().different())
        return divisors(idealv)

    def _sigma_K(self, v: tuple | NumberFieldElement, s: ComplexNumber) -> ComplexNumber:
        r"""
        Divisor function in number field $\sum_{d | (v)*ideal} N(d)**(-s)$

        INPUT:

        - ``v`` -- element in number field
        - ``s`` -- complex number

        OUTPUT:

        - ComplexNumber, the divisor sum

        EXAMPLES::

            sage: from maass_forms_hilbert.modform.eisenstein_series import HilbertEisensteinSeries
            sage: K = QuadraticField(2)
            sage: ideal = K.ideal(1)
            sage: spectral_parameter = (CC(0.5,1), CC(0.5,1))
            sage: E = HilbertEisensteinSeries(K, ideal, spectral_parameter)
            sage: result = E._sigma_K((1, 0), 1)
            sage: result is not None
            True
        """
        return sum(ideal_c.norm() ** s for ideal_c in self._ideals_dividing_vD(v))

    def coefficients(self, M: tuple[tuple[Integer_t]]) -> HilbertMaassCoefficients:
        r"""
        Return the coefficients of the Eisenstein series (not implemented).

        INPUT:

        - ``M`` -- tuple of integer bounds

        OUTPUT:

        - HilbertMaassCoefficients

        EXAMPLES::

            sage: from maass_forms_hilbert.modform.eisenstein_series import HilbertEisensteinSeries
            sage: K = QuadraticField(2)
            sage: ideal = K.ideal(1)
            sage: spectral_parameter = (CC(0.5,1), CC(0.5,1))
            sage: E = HilbertEisensteinSeries(K, ideal, spectral_parameter)
            sage: E.coefficients(((-1, 1), (-1, 1)))  # doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
        # different = self.number_field().different()
        # dual_ideals = [ideal**-1 * different**-1 for ideal in self._dual_ideals]
        # for V in cartesian_product([range(-m0[0], m0[1] + 1) for m0 in M]):
        #     for ideala_dual in dual_ideals:
        #         c = self.dual_ideal_element(V, self._ideal)
        #
        # HilbertMaassCoefficients(X, M, self._dual_ideals)
