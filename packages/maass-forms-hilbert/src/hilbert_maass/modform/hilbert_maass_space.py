"""
Class for spaces of Hilbert Maass forms.
"""
from hilbert_modgroup.all import HilbertModularGroup
from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup_class
from hilbert_modgroup.pullback import HilbertPullback
from numpy import linspace
from sage.modules.module import Module
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.number_field.number_field_base import NumberField
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.rings.real_mpfr import RealNumber
from typing import ParamSpec
from .hilbert_maass_element import HilbertMaassForm_Element

P = ParamSpec('P')


class HilbertMaassFormSpace(Module):

    Element = HilbertMaassForm_Element

    def __init__(self, group: HilbertModularGroup_class | NumberField, **kwargs: P.kwargs) -> None:
        if not isinstance(group, HilbertModularGroup_class):
            group = HilbertModularGroup(group)
        self._group = group
        self._numerical_precision = kwargs.pop('numerical_precision', 53)
        self._complex_field = ComplexField(self._numerical_precision)
        self._cuspidal = kwargs.pop('cuspidal', False)
        self._number_field = group.base_ring().number_field()
        different = self._number_field.different()
        representatives = self.group().ideal_cusp_representatives()
        self._dual_ideals = [ideal**-1*different**-1 for ideal in representatives]
        self._pullback = None
        super(HilbertMaassFormSpace, self).__init__(group.base_ring(), **kwargs)

    def group(self):
        """
        The group of self.

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: H.group()
            Hilbert Modular Group PSL(2) over Maximal Order in Number Field in a with defining...

        """
        return self._group

    def number_field(self):
        """
        The group of self.

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: H.number_field()
            Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?
            sage: H.number_field() == H.group().base_ring().number_field()
            True

        """

        return self._number_field

    def dual_ideals(self):
        return self._dual_ideals
    def is_cuspidal(self):
        """
        Is self cuspidal

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: H.is_cuspidal()
            False
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
            sage: H.is_cuspidal()
            True

        """
        return self._cuspidal

    def pullback(self):
        """
        An instance of HilbertPullback on the group of self.

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: H.pullback()
            Pullback class for Hilbert Modular Group PSL(2) over Maximal Order in Number Field...

        """
        if not self._pullback:
            self._pullback = HilbertPullback(self.group())
        return self._pullback

    def _element_constructor_(self, s: tuple | HilbertMaassForm_Element = None,
                              check: bool = False, **kwargs: P.kwargs):
        r"""
        Construct an element of this finite quadratic module.

        EXAMPLES::

            sage: from fqm_weil.all import FiniteQuadraticModule
        """
        if isinstance(s, self.element_class):
            if check and s.parent() != self:
                raise ValueError(f"Can not construct an element of {self} from {s}")
            # Only return x if the parent *is* self to avoid coercion problems.
            if s.parent() is self:
                return s
            s = s.list()
        if check:
            if isinstance(s, (RealNumber, ComplexNumber, int, float)):
                s = [s] * self._number_field.degree()
        return self.element_class(self, s, **kwargs)

    def check_interval(self, r_start: float | RealNumber,
                             r_stop: float | RealNumber,
                             fixed_params: tuple[float | RealNumber] = None,
                       nsteps: int = 10,
                       Y1: float | RealNumber = None,
                       Y2: float | RealNumber = None,
                       M: tuple[int] = None,
                       ideala: NumberFieldFractionalIdeal = None,
                       coeff: tuple = None
                       ) -> list:
        G = self.an_element()
        CF = self._complex_field
        ideala = ideala or self._number_field.fractional_ideal(1)
        Y1v = Y1 or (0.75, 0.75)
        Y2v = Y2 or (0.73, 0.73)
        M = M or (4, 4)
        result = []
        f0 = 0
        f1 = 0
        r_previous = r_start
        for r in linspace(r_start, r_stop, nsteps):
            if fixed_params is None:
                s = (CF(0.5, r),) * len(M)
            else:
                # Replace fixed params = None by the variable r
                rs = [r0 or r for r0 in fixed_params]
                s = tuple(CF(0.5, r0) for r0 in rs)
            print(s)
            G.spectral_parameter = s
            C1 = G.compute_coefficients(s, ideala, ideala, M=M, Y=Y1v)
            C2 = G.compute_coefficients(s, ideala, ideala, M=M, Y=Y2v)
            result.append((s, C1, C2))
            f0 = self.functional(C1, C2, coeff=coeff)
            if f0 * f1 < 0 and max(abs(f0), abs(f1)) < 1e-2:
                print(f"Sign change {coeff} in interval: [{r_previous}, {r}]")
            f0 = self.functional(C1, C2)
            r_previous = r
            f1 = f0
        return result

    def functional(self, C1, C2, coeff=(0, 1)):
        return C1[coeff].real() - C2[coeff].real()

