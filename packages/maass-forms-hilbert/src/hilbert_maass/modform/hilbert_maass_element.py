"""
Classes For Hilbert-Maass forms

"""
import logging
from hilbert_maass.coefficients import get_pb_pts
from hilbert_maass.functions import bessel_prod
# from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
from hilbert_maass.utils import get_Q_from_bounds, map_tuple_to_int
from hilbert_modgroup.pullback import HilbertPullback
from sage.all import RR
from sage.categories.sets_cat import cartesian_product
from sage.functions.other import ceil
from sage.matrix.constructor import matrix
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.integer import Integer
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.rings.real_mpfr import RealNumber
from sage.structure.sage_object import SageObject
from typing import ParamSpec, NoReturn
from sage.structure.element import ModuleElement, Matrix

from hilbert_maass.coefficients import matrix_element
from hilbert_maass.utils import map_int_to_tuple

P = ParamSpec('P')
log = logging.getLogger(__name__)


class HilbertMaassForm_Element(ModuleElement):

    def __init__(self, parent: 'HilbertMaassFormSpace',
                 spectral_parameter: tuple[ComplexNumber | RealNumber],
                 **kwargs: P.kwargs) -> None:
        super(HilbertMaassForm_Element, self).__init__(parent, **kwargs)
        self.cuspidal = parent.is_cuspidal()
        self.spectral_parameter = spectral_parameter
        self._number_field = parent.number_field()
        self._complex_field = ComplexField(self.spectral_parameter[0].parent().prec())
        self.has_coefficients = False
        self._coefficients = {c: {} for c in range(self.parent().group().ncusps())}
        self._pullback = HilbertPullback(self.parent().group())
        different = self.parent().number_field().different()
        self._dual_ideal_matrix = {
            ideal:  matrix(self._pullback._get_lattice_and_ideal_basis(ideal**-1*different**-1)[0])
            for ideal in self.parent().group().ideal_cusp_representatives()
        }

    def dual_ideal_element(self, coordinates: tuple[int | Integer] or vector,
                           ideal: NumberFieldFractionalIdeal):
        return self._dual_ideal_matrix[ideal]*vector(coordinates)

    def pullback(self):
        return self._pullback

    def coefficients(self):
        return self._coefficients

    def get_pb_pts(self, M: tuple[int | Integer] = None,
                   Y: tuple = None,
                   ideala: NumberFieldFractionalIdeal = None) -> tuple:
        CF = self._complex_field
        if M is None:
            M0 = ceil((abs(self.spectral_parameter[0]) + 12) / (RR.pi() * 2) + 1)
            M = (M0, M0)
        if Y is None:
            Y = (CF(0.75), CF(0.75))
        Qs = get_Q_from_bounds(self.pullback(), M)
        Qs = tuple([ceil(q) + 5 for q in Qs])
        log.debug(f"M = {M}, Y = {Y}, Qs = {Qs}")
        if not ideala:
            ideala = self.pullback().number_field().ideal(1)
        zpb, zm = get_pb_pts(self.pullback(), Qs, ideala, Y)
        return zpb, zm, Qs, M, Y

    def compute_coefficients(self, s: tuple,
                             ideala: NumberFieldFractionalIdeal = None,
                             idealb: NumberFieldFractionalIdeal = None,
                             M: tuple[int | Integer] = None,
                             Y: tuple = None,
                             prec: int = 53,
                             sgn: str = '-', returnV: bool = False) -> Matrix:
        r"""

        INPUT:

        - ``ideala``  -- NumberField Fractional Ideal corresponding to cusp.
        - ``idealb``  --
        - ``s``       --
        - ``M``  --
        - ``Y``  --
        - ``prec``  --
        - ``cuspidal``  --
        - ``sgn``  --
        - ``returnV``  --

        EXAMPLES::

            sage: from hilbert_maass import *
            sage: M = (2,2)
            sage: spectral_parameter = [CC(1.5,1.5)]*2
            sage: X= compute_coefficients(H1, spectral_parameter, ideala, idealb,M=M)

        """
        ideala = ideala or self.parent().number_field().ideal(1)
        idealb = idealb or self.parent().number_field().ideal(1)
        zpb, zm, Qs, M, Y = self.get_pb_pts(M, Y, ideala)
        log.debug(f"M = {M}, Y = {Y}, Qs = {Qs}")
        matrixV = {}
        for V in cartesian_product([range(-M0, M0 + 1) for M0 in M]):
            v = self.dual_ideal_element(V, ideala)
            for W in cartesian_product([range(-M0, M0 + 1) for M0 in M]):
                # For cuspidal forms we don't need to compute the row corresponding to 0
                if self.cuspidal and (all(x == 0 for x in W) or all(x == 0 for x in V)):
                    matrixV[(V, W)] = 0
                else:
                    w = self.dual_ideal_element(W, idealb)
                    matrixV[(V, W)] = matrix_element(s, Qs, v, w, zpb, zm, sgn='-')
            matrixV[(V, V)] = matrixV[(V, V)] - bessel_prod(tuple(v), tuple(Y), s, sgn='-')
        RHS = {}
        # if is_cuspidal:
        # Set 0-th coefficient to 0 and 1st to 1
        t_0 = (0, ) * len(M)
        n_0 = map_tuple_to_int(t_0, -M[0], M[0])
        # tuple for 1
        t1 = self._dual_ideal_matrix[ideala] ** -1 * vector((1,) * len(M))
        t_1 = tuple(int(t11) for t11 in t1)
        n_1 = map_tuple_to_int(t_1, -M[0], M[0])

        # set_coefficient = (0,) * len(M)
        if not self.cuspidal:
            for V in cartesian_product([range(-M0, M0 + 1) for M0 in M]):
                v = self.dual_ideal_element(V, ideala)
                W = w = (0,) * len(v)
                RHS[(V, W)] = matrix_element(s, Qs, v, w, zpb, zm, sgn='+')
                if V == W:  # == 0,...,0
                    RHS[(V, W)] = RHS[(V, W)] - bessel_prod(v, tuple(Y), s, sgn='+')
        else:
            for V in cartesian_product([range(-M0, M0 + 1) for M0 in M]):
                RHS[(V, t_1)] = matrixV[(V, t_1)]
        n = prod(2 * M0 + 1 for M0 in M)
        Vmat = [[
                matrixV[map_int_to_tuple(r, -M[0], M[0], 2), map_int_to_tuple(k, -M[0], M[0], 2)]
                for k in range(n)
                ] for r in range(n)]
        Vmat = matrix(self._complex_field, n, n, Vmat)
        if not self.cuspidal:
            W = t_0
        else:
            W = t_1
        RHSmat = [
            RHS[(map_int_to_tuple(k, -M[0], M[0], 2), W)] for k in range(n)
        ]
        RHSmat = matrix(self._complex_field, n, 1, RHSmat)
        if returnV:
            return Vmat, RHSmat
        log.debug(f"n, n_0, n_1= {n, n_0, n_1}")
        if self.cuspidal:
            delete_rows = (n_0, n_1)
            Vmat = Vmat.delete_rows(delete_rows)
            RHSmat = RHSmat.delete_rows(delete_rows)
            Vmat = Vmat.delete_columns(delete_rows)
            X = Vmat.solve_right(-RHSmat)
            # Add back coefficients for 0 and 1
            rows = [r for i, r in enumerate(X.rows()) if i < n_0]
            rows += [(self._complex_field(0),)]
            rows += [r for i, r in enumerate(X.rows()) if i >= n_0 and i < n_1]
            rows += [(self._complex_field(1),)]
            rows += [r for i, r in enumerate(X.rows()) if i >= n_1]
            X = matrix(rows)
        else:
            X = Vmat.solve_right(-RHSmat)
        C = HilbertMaassCoefficients(X, M)
        self._coefficients = C
        return C


class HilbertMaassCoefficients(SageObject):


    def __init__(self, coefficients: Matrix, M: tuple[int | Integer], **kwargs: P.kwargs) -> None:
        super(HilbertMaassCoefficients, self).__init__(**kwargs)
        self._coefficients = coefficients
        self._M = M

    def __getitem__(self, key: tuple | int | Integer):
        if isinstance(key, tuple) and len(key) == 3:
            cusp, *v = key
        else:
            cusp = 0
            v = key
        if isinstance(v, tuple):
            v = map_tuple_to_int(v, -self._M[0], self._M[0])
        return self._coefficients[v][0]

    def keys(self, as_elements=False):
        keys = [map_int_to_tuple(v, -self._M[0], self._M[0], 2) for v in range(len(self._coefficients))]


    def __iter__(self):
        return iter(self._coefficients)


def HilbertMaassForm(group: 'HilbertModularGroup',
                     spectral_parameter: tuple[ComplexNumber | RealNumber],
                     **kwargs: P.kwargs) -> HilbertMaassForm_Element:
    from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
    space = HilbertMaassFormSpace(group, **kwargs)
    return HilbertMaassForm_Element(space, spectral_parameter=spectral_parameter)
