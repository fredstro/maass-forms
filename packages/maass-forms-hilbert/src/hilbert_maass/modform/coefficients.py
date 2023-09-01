import json

from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
from sage.categories.sets_cat import cartesian_product
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.complex_mpfr import ComplexNumber, ComplexField
from sage.rings.real_mpfr import RealField, RealNumber
from typing import ParamSpec

import logging
from hilbert_maass.modform.utils import ideal_coordinates, map_tuple_to_int, map_int_to_tuple, \
    get_Q_from_bounds
from sage.functions.other import ceil
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.rings.number_field.number_field_element import NumberFieldElement
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.structure.element import Matrix
from sage.structure.sage_object import SageObject

from ..functions.functions import bessel_prod, exp_trace_prod

from .utils import Integer_t, length_from_M, cartesian_product_from_M, dual_ideal_element, \
    complex_tuple_to_json, complex_tuple_from_json

P = ParamSpec('P')
log = logging.getLogger(__name__)


class HilbertMaassCoefficients(SageObject):
    """
    Coefficients of Hilbert Maass Forms

    """
    def __init__(self, coefficients: Matrix, M: tuple[Integer_t],
                 spectral_parameter: tuple[ComplexNumber | RealNumber],
                 space: 'HilbertMaassFormSpace',
                 coordinate_ideals: list[NumberFieldFractionalIdeal] = None,
                 check: bool = True,
                 **kwargs: P.kwargs) -> None:
        r"""
        Init self.

        INPUT:

        - ``coefficients`` -- matrix of coefficients
        - ``M`` -- tuple of integers giving the limits used to map the indices of the matrix to tuples
        - ``spectral_parameter`` -- tuple of complex or real numbers
        - ``space`` -- Hilbert Maass form space
        - ``coordinate_ideals`` -- list of fractional ideals indexing the different components
        - ``check`` -- boolean, if True, check that the coefficients are of the correct dimension.
        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = Matrix(RR, [[1],[2],[3],[4],[5],[6],[7],[8],[9]])
            sage: C = HilbertMaassCoefficients(Cmat, ((-1,1),(-1,1)), spectral_parameter, H); C
            Coefficients of a Hilbert Maass form with M=((-1, 1), (-1, 1)) and 1 cusp

        """
        super(HilbertMaassCoefficients, self).__init__(**kwargs)
        if not coordinate_ideals:
            representatives = space.group().ideal_cusp_representatives()
            different = space.group().base_ring().number_field().different()
            coordinate_ideals = [ideal ** -1 * different ** -1 for ideal in representatives]
        if check:
            if any(M0 == M1 == 0 for M0, M1 in M):
                raise ValueError('M must be non-zero')
            if coefficients.nrows() != length_from_M(M):
                raise ValueError(f'coefficients has wrong number of rows: {coefficients.nrows()}'
                                 f' != {length_from_M(M)}')
            if coefficients.ncols() != len(coordinate_ideals):
                raise ValueError('coefficients has number of cols: {coefficients.ncols()}'
                                 f' != {len(coordinate_ideals)}')
        from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
        if not coordinate_ideals and not space or not isinstance(space, HilbertMaassFormSpace):
            raise ValueError("Need to either specify coordinate ideals or a Maass form")

        self._coordinate_ideals = coordinate_ideals
        self._coefficients = coefficients
        self._M = M
        self._spectral_parameter = spectral_parameter
        self._space = space

    def to_json(self) -> dict:
        """
            JSON representation of self.

            EXAMPLES::

                sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
                sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
                sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
                sage: Cmat = Matrix(RR, [[1],[2],[3],[4],[5],[6],[7],[8],[9]])
                sage: C = HilbertMaassCoefficients(Cmat, ((-1,1),(-1,1)), spectral_parameter, H)
                sage: C.to_json()
                {'M': ((-1, 1), (-1, 1)),
                 'coefficients': [['1.00000000000000'],
                  ['2.00000000000000'],
                  ['3.00000000000000'],
                  ['4.00000000000000'],
                  ['5.00000000000000'],
                  ['6.00000000000000'],
                  ['7.00000000000000'],
                  ['8.00000000000000'],
                  ['9.00000000000000']],
                 'prec': 53,
                 'space': {'cuspidal': False,
                  'number_field': {'names': ('a',), 'polynomial': 'x^2 - 2'}},
                 'spectral_parameter': [{'prec': 53,
                         'val': '0.500000000000000 + 1.00000000000000*I'},
                        {'prec': 53, 'val': '0.500000000000000 + 1.00000000000000*I'}]}
        """
        return {
            'M': tuple((int(M0[0]), int(M0[1])) for M0 in self._M),
            'coefficients': [[str(x) for x in r] for r in self._coefficients],
            'prec': int(self._coefficients.base_ring().prec()),
            'spectral_parameter': complex_tuple_to_json(self._spectral_parameter),
            'space': self._space.to_json()
        }

    @classmethod
    def from_json(cls, data: str | dict) -> 'HilbertMaassCoefficients':
        """
            Create an instance of HilbertMaassCoefficients from json formatted dict or string.

            EXAMPLES::

                sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
                sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
                sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
                sage: Cmat = Matrix(RR, [[1],[2],[3],[4],[5],[6],[7],[8],[9]])
                sage: C = HilbertMaassCoefficients(Cmat, ((-1,1),(-1,1)), spectral_parameter, H)
                sage: HilbertMaassCoefficients.from_json(C.to_json()) == C
                True
                sage: import json
                sage: json_string = json.dumps(C.to_json())
                sage: HilbertMaassCoefficients.from_json(json_string) == C
                True

        """
        if isinstance(data, str):
            data = json.loads(data)
        CF = ComplexField(data['prec'])
        coefficients = matrix(CF, data['coefficients'])
        from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
        M = tuple(tuple(x) for x in data['M'])
        return cls(
            coefficients, M=M,
            spectral_parameter=complex_tuple_from_json(data['spectral_parameter']),
            space=HilbertMaassFormSpace.from_json(data['space'])
        )

    def space(self):
        return self._space

    def coefficient_matrix(self):
        return self._coefficients

    def spectral_parameter(self):
        return self._spectral_parameter

    def __eq__(self, other):
        if not isinstance(other, HilbertMaassCoefficients):
            return False
        return self.space() == other.space() and \
            self.coefficient_matrix() == other.coefficient_matrix() and \
            self.spectral_parameter() == other.spectral_parameter()

    def __repr__(self):
        """
        String representation of self.

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = Matrix(RR, [[1],[2],[3],[4],[5],[6],[7],[8],[9]])
            sage: C = HilbertMaassCoefficients(Cmat, ((-1,1),(-1,1)), spectral_parameter, H)
            sage: C.__repr__()
            'Coefficients of a Hilbert Maass form with M=((-1, 1), (-1, 1)) and 1 cusp'


        """
        num_cusps = len(self._coordinate_ideals)
        return f"Coefficients of a Hilbert Maass form with M={self._M} and" \
               f" {num_cusps} cusp{'s' if num_cusps > 1 else ''}"

    def __getitem__(self, key: tuple | Integer_t) -> ComplexNumber:
        """
        Get one coefficient of self.

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = Matrix(RR, [[1],[2],[3],[4],[5],[6],[7],[8],[9]])
            sage: C = HilbertMaassCoefficients(Cmat, ((-1,1),(-1,1)), spectral_parameter, H)
            sage: C[0]
            1.00000000000000
            sage: C[0, 1]
            8.00000000000000
            sage: C[0, (0, -1)]
            2.00000000000000
            sage: C[0, (0, 1)]
            8.00000000000000
            sage: C[0, (1, 1)]
            9.00000000000000

            TESTS::

            sage: C[0,(0, 2)]
            Traceback (most recent call last):
            ...
            IndexError: Tuple element (0, 2) is out of bounds!
            sage: C[1,(0, 0)]
            Traceback (most recent call last):
            ...
            ValueError: Cusp index 1 out of bounds.
        """
        cusp = 0
        v = None
        if isinstance(key, tuple) and len(key) == 2 and isinstance(key[1], tuple):
            cusp, v = key
        elif isinstance(key, (int, Integer)) or isinstance(key, tuple):
            v = key
        elif isinstance(key, NumberFieldElement):
            cusp = 0
            v = key
        if cusp < 0 or cusp >= len(self._coordinate_ideals):
            raise ValueError(f"Cusp index {cusp} out of bounds.")
        if isinstance(v, NumberFieldElement):
            v = ideal_coordinates(self._coordinate_ideals[cusp], v)
        if isinstance(v, tuple):
            v = map_tuple_to_int(v, self._M)
        try:
            return self._coefficients.column(cusp)[v]
        except IndexError:
            raise ValueError(f'Can not get coefficient for {v}')

    def keys(self, as_elements=False) -> list:
        keys_tuple = [map_int_to_tuple(v, self._M)
                      for v in range(self._coefficients.nrows())]
        if not as_elements:
            return keys_tuple
        else:
            basis = self._coordinate_ideals[0].integral_basis()
            return [sum(t[i] * basis[i] for i in range(len(basis))) for t in keys_tuple]

    def __iter__(self) -> iter:
        return iter(self._coefficients.column(0))


@cached_method
def get_pb_pts_set_params(space: 'HilbertMaassFormSpace',
                          spectral_parameter: tuple[ComplexNumber] = None,
                          M: tuple[tuple[Integer_t]] = None,
                          Y: tuple = None,
                          smax: float | RealNumber = None,
                          ideala: NumberFieldFractionalIdeal = None) -> tuple:
    spectral_parameter = spectral_parameter or (ComplexField(53)(0.5,10), )
    CF = spectral_parameter[0].parent()
    if not spectral_parameter and not smax:
        raise ValueError("Need either spectral parameter or smax set.")
    if spectral_parameter and not smax:
        smax = max(abs(s) for s in spectral_parameter)
    if not isinstance(M, tuple):
        M0 = M or ceil((smax + 12) / (6.28318530717959) + 1)
        M = ((-M0, M0),) * space.number_field().absolute_degree()
    if isinstance(M, tuple) and not isinstance(M[0], tuple):
        M = (M,) * space.number_field().absolute_degree()
    if Y is None:
        Y = (CF(0.75), ) * space.number_field().absolute_degree()
    Qs = get_Q_from_bounds(space.pullback(), M)
    Qs = tuple([ceil(q) + 5 for q in Qs])
    log.debug(f"M = {M}, Y = {Y}, Qs = {Qs}")
    if not ideala:
        ideala = space.pullback().number_field().ideal(1)
    zpb, zm = get_pb_pts(space.pullback(), Qs, ideala, Y)
    return zpb, zm, Qs, M, Y


def get_pb_pts(P, Q: tuple, ideala: NumberFieldFractionalIdeal, Y: tuple, prec: int = 53) -> tuple:
    """
    Get the list of points in the scaled lattice together with the corresponding pullbacks.
    """
    n = P.number_field().degree()
    CF = RealField(prec)
    ideala_matrix = matrix(P.basis_matrix_ideal(ideala))
    basis_matrix_m = ideala_matrix * diagonal_matrix([CF(1) / CF(2 * q) for q in Q])
    Q_combination = [range(1 - q, q + 1) for q in Q]
    zmpb = []
    zm = []
    for m in cartesian_product(Q_combination):
        xm = basis_matrix_m * vector(m)
        zm_elt = UpperHalfPlaneProductElement([(xm[i], Y[i]) for i in range(n)])
        zm.append(zm_elt)
        zmpb.append(P.reduce(zm_elt))
    return zmpb, zm


def compute_coefficients(space: 'HilbertMaassFormSpace',
                         spectral_parameter: tuple[complex | ComplexNumber],
                         ideala: NumberFieldFractionalIdeal = None,
                         idealb: NumberFieldFractionalIdeal = None,
                         M: tuple[Integer_t] = None,
                         Y: tuple[float | RealNumber] = None,
                         sgn: str = '-', returnV: bool = False) -> 'HilbertMaassCoefficients':
    r"""

    INPUT:

    - ``space`` -- Hilbert Maass form space
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

        sage: from hilbert_maass.all import HilbertMaassFormSpace
        sage: from hilbert_maass.modform.coefficients import compute_coefficients
        sage: M = (2,2)
        sage: s = CC(1.5,1.5), CC(1.5,1.5)
        sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
        sage: X= compute_coefficients(H, s, M = 2); X
        Coefficients of a Hilbert Maass form with M=((-2, 2), (-2, 2)) and 1 cusp
    """
    complex_field = spectral_parameter[0].parent()
    ideala = ideala or space.number_field().ideal(1)
    idealb = idealb or space.number_field().ideal(1)
    smax = max(ceil(abs(s0)) for s0 in spectral_parameter)
    zpb, zm, Qs, M, Y = get_pb_pts_set_params(space, smax=smax,
                                              M=M, Y=Y, ideala=ideala)
    log.debug(f"M = {M}, Y = {Y}, Qs = {Qs}")
    matrixV = {}
    for V in cartesian_product_from_M(M):
        v = dual_ideal_element(V, ideala)
        for W in cartesian_product_from_M(M):
            # For cuspidal forms we don't need to compute the row corresponding to 0
            if space.is_cuspidal() and (all(x == 0 for x in W) or all(x == 0 for x in V)):
                matrixV[(V, W)] = 0
            else:
                w = dual_ideal_element(W, idealb)
                matrixV[(V, W)] = matrix_element(spectral_parameter, Qs, v, w, zpb, zm, sgn='-')
        matrixV[(V, V)] = matrixV[(V, V)] - bessel_prod(tuple(v), tuple(Y), spectral_parameter,
                                                        sgn='-')
    RHS = {}
    # if is_cuspidal:
    # Set c(0)=0 and c(delta)=1 where delta >>0 is generator of the index ideal.
    t_0 = (0,) * space.number_field().absolute_degree()
    n_0 = map_tuple_to_int(t_0, M)
    # tuple for 1
    # Try this first:
    ideala_dual = ideala ** -1 * ideala.number_field().different() ** -1
    x, y = ideala_dual.gens_two()
    delta = None
    for delta_test in [x + y, x - y, -x - y, -x + y]:
        if not delta_test.is_totally_positive():
            continue
        if ideala_dual != ideala_dual.number_field().fractional_ideal(delta_test):
            continue
        delta = delta_test
        break
    if not delta:
        raise ArithmeticError(f"Cannot find a totally positive generator for {ideala_dual}")
    # Coordinate vector of delta
    t_1 = ideal_coordinates(ideala_dual, delta)
    n_1 = map_tuple_to_int(t_1, M)
    if not space.is_cuspidal():
        for V in cartesian_product_from_M(M):
            V = tuple(V)
            v = dual_ideal_element(V, ideala)
            W = w = (0,) * len(v)
            RHS[(V, W)] = matrix_element(spectral_parameter, Qs, v, w, zpb, zm, sgn='+')
            if V == W:  # == 0,...,0
                RHS[(V, W)] = RHS[(V, W)] - bessel_prod(v, tuple(Y), spectral_parameter, sgn='+')
    else:
        for V in cartesian_product_from_M(M):
            RHS[(V, t_1)] = matrixV[(V, t_1)]
    n = length_from_M(M)
    Vmat = [[
        matrixV[map_int_to_tuple(r, M), map_int_to_tuple(k, M)]
        for k in range(n)
    ] for r in range(n)]
    Vmat = matrix(complex_field, n, n, Vmat)
    if not space.is_cuspidal():
        W = t_0
    else:
        W = t_1
    RHSmat = [
        RHS[(map_int_to_tuple(k, M), W)] for k in range(n)
    ]
    RHSmat = matrix(complex_field, n, 1, RHSmat)
    if returnV:
        return Vmat, RHSmat
    log.debug(f"n, n_0, n_1= {n, n_0, n_1}")
    # Make sure that n0 < n1:
    c0 = complex_field(0)
    c1 = complex_field(1)
    if n_0 > n_1:
        n_0, n_1 = n_1, n_0
        c0, c1 = c1, c0
    if space.is_cuspidal():
        delete_rows = (n_0, n_1)
        Vmat = Vmat.delete_rows(delete_rows)
        RHSmat = RHSmat.delete_rows(delete_rows)
        Vmat = Vmat.delete_columns(delete_rows)
        X = Vmat.solve_right(-RHSmat)
        # Add back coefficients for 0 and 1
        rows = X.rows()[0:n_0] + [(c0,)] + X.rows()[n_0:n_1 - 1] \
                               + [(c1,)] + X.rows()[n_1 - 1:]
        X = matrix(rows)
    else:
        X = Vmat.solve_right(-RHSmat)
    return HilbertMaassCoefficients(X, M, spectral_parameter=spectral_parameter,
                                    space=space, coordinate_ideals=space.dual_ideals())


def matrix_element(s: tuple, Q: tuple, v: tuple, w: tuple, zpb_v: list, zm_v: list,
                   sgn: str = '+') -> ComplexNumber:
    factor = prod(2 * q for q in Q)
    summa = 0
    for m, zm in enumerate(zm_v):  # m in cartesian_product(Q_combination):
        xm = zm.real()
        zmpb = zpb_v[m]
        ympb = zmpb.imag()
        xmpb = zmpb.real()
        bes = bessel_prod(w, tuple(ympb), s, sgn=sgn)
        exp_arg = (xmpb[0] * w[0] - xm[0] * v[0], xmpb[1] * w[1] - xm[1] * v[1])
        exp_val = exp_trace_prod(exp_arg)
        term = bes * exp_val
        summa += term
    return summa / factor
