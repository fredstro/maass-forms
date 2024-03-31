"""
Class for spaces of Hilbert Maass forms.
"""
from hilbert_modgroup.all import HilbertModularGroup
from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup_class
from hilbert_modgroup.pullback import HilbertPullback
from numpy import linspace
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.modules.module import Module
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.number_field.number_field_base import NumberField as NumberFieldBase
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.rings.real_mpfr import RealNumber as RealNumber_class
from typing import ParamSpec, Any

from sage.structure.element import ModuleElement

from .hilbert_maass_element import HilbertMaassForm_Element
from .utils import number_field_from_json, number_field_to_json, Real_t, Integer_t

P = ParamSpec('P')


class HilbertMaassFormSpace(Module):

    Element = HilbertMaassForm_Element

    def __init__(self, group: HilbertModularGroup_class | NumberFieldBase, **kwargs: P.kwargs) -> None:
        if not isinstance(group, HilbertModularGroup_class):
            group = HilbertModularGroup(group)
        self._group = group
        self._cuspidal = kwargs.pop('cuspidal', False)
        self._number_field = group.base_ring().number_field()
        different = self._number_field.different()
        representatives = self.group().ideal_cusp_representatives()
        self._dual_ideals = [ideal**-1*different**-1 for ideal in representatives]
        self._pullback = None
        self._dual_ideal_basis_matrix = []
        super(HilbertMaassFormSpace, self).__init__(group.base_ring(), **kwargs)

    def to_json(self):
        """
        JSON representation of self.

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: H.to_json()
            {'cuspidal': False, 'number_field': {'names': ['a'], 'polynomial': 'x^2 - 2'}}

        """
        return {
            'number_field': number_field_to_json(self.number_field()),
            'cuspidal': self._cuspidal
        }

    @classmethod
    def from_json(cls, data):
        """
        Create an instance of self from JSON data.

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: H.from_json(H.to_json()) == H
            True

        """
        nf = number_field_from_json(data['number_field'])
        group = HilbertModularGroup(nf)
        return cls(group, cuspidal=data['cuspidal'])

    def __repr__(self):
        return f"HilbertMaassFormSpace({self.group()})"

    def __eq__(self, other: Any) -> bool:
        """
        Is self equal to other.

        Note: Isomorphic number fields (e.g. QuadraticField(2) and NumberField(x^2-2)
              have equal level, so we use the level for comparison of groups.

        INPUT:

        - `other` -- object to compare self with

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace
            sage: H1 = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: H1 == H1
            True
            sage: H1 != H1
            False
            sage: H2 = HilbertMaassFormSpace(NumberField(x^2-2, names='a'), cuspidal=False)
            sage: H1 == H2
            True
            sage: H1 != H2
            False
            sage: H3 = HilbertMaassFormSpace(NumberField(x^2-2, names='a'), cuspidal=True)
            sage: H1 == H3
            False
            sage: H1 != H3
            True


        """
        if not isinstance(other, HilbertMaassFormSpace):
            return False
        if self.is_cuspidal() != other.is_cuspidal():
            return False
        # Need to check isomorphic fields
        if not self.number_field().is_isomorphic(other.number_field()):
            return False
        this_level = self.group().level()
        other_level_gens = other.group().level().gens_reduced()
        return self.number_field().fractional_ideal(other_level_gens) == this_level

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

    @cached_method
    def dual_ideal_basis_matrix(self, ideal: NumberFieldFractionalIdeal | Integer_t = None):
        if ideal is None:
            ideal = 0
        if isinstance(ideal, NumberFieldFractionalIdeal):
            ideal = self._dual_ideals.index(ideal)
        if not self._dual_ideal_basis_matrix:
            self._dual_ideal_basis_matrix = [
                matrix(self.pullback()._get_lattice_and_ideal_basis(dual_ideal)[0])
                for dual_ideal in self.dual_ideals()
            ]
        return self._dual_ideal_basis_matrix[ideal]

    @cached_method
    def dual_ideal_element(self, coordinates: tuple[Integer_t] | ModuleElement,
                           ideal: NumberFieldFractionalIdeal | Integer_t = None):
        return self.dual_ideal_basis_matrix(ideal) * vector(coordinates)

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


        """
        if isinstance(s, self.element_class):
            if check and s.parent() != self:
                raise ValueError(f"Can not construct an element of {self} from {s}")
            # Only return x if the parent *is* self to avoid coercion problems.
            if s.parent() is self:
                return s
            s = s.list()
        if check:
            if isinstance(s, (RealNumber_class, ComplexNumber, int, float)):
                s = [s] * self._number_field.degree()
        return self.element_class(self, s, **kwargs)

    def check_interval(self, r_start: Real_t,
                             r_stop: Real_t,
                             fixed_params: tuple[Real_t] = None,
                       nsteps: int = 10,
                       Y1: Real_t = None,
                       Y2: Real_t = None,
                       M: tuple[int] = None,
                       ideala: NumberFieldFractionalIdeal = None,
                       coeff: tuple = None
                       ) -> list:
        G = self.an_element()
        CF = ComplexField(r_start.parent().prec())
        ideala = ideala or self._number_field.fractional_ideal(1)
        Y1v = Y1 or (0.75, 0.75)
        Y2v = Y2 or (0.73, 0.73)
        M = M or ((-4, 4),) * 2
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
            G._spectral_parameter = s
            C1 = G.compute_coefficients(s, ideala, ideala, M=M, Y=Y1v)
            C2 = G.compute_coefficients(s, ideala, ideala, M=M, Y=Y2v)
            result.append((s, C1, C2))
            coeff = coeff or (0, 1)
            f0 = self.functional(C1, C2, coeff=coeff)
            if f0 * f1 < 0 and max(abs(f0), abs(f1)) < 1e-2:
                print(f"Sign change {coeff} in interval: [{r_previous}, {r}]")
            f0 = self.functional(C1, C2)
            r_previous = r
            f1 = f0
        return result

    def functional(self, C1, C2, coeff=(0, 1)):
        return C1[coeff].real() - C2[coeff].real()

