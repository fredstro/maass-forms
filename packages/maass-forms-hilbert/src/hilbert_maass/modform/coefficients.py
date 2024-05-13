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

from ..functions.functions import bessel_prod
from ..functions.functions_cy import exp_trace_prod_dp, bessel_prod_dp2

from .utils import Integer_t, length_from_M, cartesian_product_from_M, dual_ideal_element, \
    complex_tuple_to_json, complex_tuple_from_json, Real_t, dual_ideal, ideal_generator, \
    Complex_t, coefficient_dict_to_json, coefficient_dict_from_json, is_tuple_zero, mongo_cache

from comp_manager.decorators import mongo_cache

P = ParamSpec('P')
log = logging.getLogger(__name__)


class HilbertMaassCoefficients(SageObject):
    """
    Coefficients of Hilbert Maass Forms

    """
    def __init__(self, coefficients: Matrix,
                 M: tuple[Integer_t | tuple[Integer_t]],
                 spectral_parameter: tuple[Complex_t | Real_t],
                 space: 'HilbertMaassFormSpace', Y: tuple[Real_t] = None,
                 coordinate_ideals: tuple[NumberFieldFractionalIdeal] = None,
                 set_coefficients: dict = None,
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
        - ``set_coefficients`` -- dictionary with coefficients to set.

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
            coordinate_ideals = [dual_ideal(ideal) for ideal in representatives]
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
        self._set_coefficients = set_coefficients
        self._M = M
        self._spectral_parameter = spectral_parameter
        self._space = space
        self._Y = Y or tuple()

    def to_json(self) -> dict:
        """
            JSON representation of self.

            EXAMPLES::

                sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
                sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
                sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
                sage: Y = (1.0,1.0)
                sage: Cmat = Matrix(RR, [[1],[2],[3],[4],[5],[6],[7],[8],[9]])
                sage: set_coefficients = {(0,0): 0, (0,1): 1}
                sage: C = HilbertMaassCoefficients(Cmat, ((-1,1),(-1,1)), spectral_parameter, H, Y,
                ....:   set_coefficients=set_coefficients)
                sage: C.to_json()
                {'M': ((-1, 1), (-1, 1)),
                 'Y': (1.0, 1.0),
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
                  'set_coefficients': {'[0, 0]': {'prec': 53, 'val': '0.000000000000000'},
                                       '[0, 1]': {'prec': 53, 'val': '1.00000000000000'}},
                 'space': {'cuspidal': False,
                  'number_field': {'names': ['a'], 'polynomial': 'x^2 - 2'}},
                 'spectral_parameter': [{'prec': 53,
                         'val': '0.500000000000000 + 1.00000000000000*I'},
                        {'prec': 53, 'val': '0.500000000000000 + 1.00000000000000*I'}]}
        """
        return {
            'M': tuple((int(M0[0]), int(M0[1])) for M0 in self._M),
            'coefficients': [[str(x) for x in r] for r in self._coefficients],
            'prec': int(self._coefficients.base_ring().prec()),
            'spectral_parameter': complex_tuple_to_json(self._spectral_parameter),
            'set_coefficients': coefficient_dict_to_json(self._set_coefficients),
            'space': self._space.to_json(),
            'Y': tuple(float(y) for y in self._Y)
        }

    @classmethod
    def from_json(cls, data: str | dict) -> 'HilbertMaassCoefficients':
        """
            Create an instance of HilbertMaassCoefficients from json formatted dict or string.

            EXAMPLES::

                sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
                sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
                sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
                sage: Y = (0.5, 0.5)
                sage: Cmat = Matrix(RR, [[1],[2],[3],[4],[5],[6],[7],[8],[9]])
                sage: set_coefficients = {(0,0): 0, (0,1): 1}
                sage: C = HilbertMaassCoefficients(Cmat, ((-1,1),(-1,1)), spectral_parameter, H, Y,
                ....:   set_coefficients=set_coefficients)
                sage: HilbertMaassCoefficients.from_json(C.to_json()) == C
                True
                sage: import json
                sage: json_string = json.dumps(C.to_json())
                sage: HilbertMaassCoefficients.from_json(json_string) == C
                True

        """
        if isinstance(data, str):
            data = json.loads(data)
        if set(data.keys()) != \
           {'prec', 'coefficients', 'M', 'Y', 'spectral_parameter', 'space', 'set_coefficients'}:
            raise ValueError("Not a valid JSON representation of HilbertMaassCoefficients")
        CF = ComplexField(data['prec'])
        coefficients = matrix(CF, data['coefficients'])
        from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
        M = tuple(tuple(x) for x in data['M'])
        Y = tuple(RealField(data['prec'])(x) for x in data['Y'])
        return cls(
            coefficients, M=M,
            spectral_parameter=complex_tuple_from_json(data['spectral_parameter']),
            space=HilbertMaassFormSpace.from_json(data['space']),
            Y=Y,
            set_coefficients=coefficient_dict_from_json(data['set_coefficients']),
        )

    def M(self):
        return self._M

    def Y(self):
        return self._Y

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
            self.spectral_parameter() == other.spectral_parameter() and \
            self._Y == other._Y and \
            self._M == other._M and \
            self._set_coefficients == other._set_coefficients

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
            sage: C[(0, 1)]
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
                          ideala: NumberFieldFractionalIdeal = None,
                          use_symmetry: bool = False) -> tuple:
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
    if not ideala:
        ideala = space.pullback().number_field().ideal(1)
    Qs = get_Q_from_bounds(space.pullback(), M)
    # Try to find best value of Y
    Y = find_max_y(space, M, Qs=Qs, starting_Y=Y)
    # We try with given Y and if it doesn't work we keep decreasing Y until it does.
    try:
        zpb, zm = get_pb_pts(space, Qs, ideala, Y, use_symmetry=use_symmetry,
                             prec=CF.prec())
    except ArithmeticError as e:
        msg = f"Could not find good pullback points. Error: {e}"
        log.debug(msg)
        raise ArithmeticError(msg)
    return zpb, zm, Qs, M, Y


@cached_method()
def find_max_y(space: 'HilbertMaassFormSpace',
               M: tuple[tuple[Integer_t]],
               Qs: tuple[tuple[Integer_t]] = None,
               starting_Y: tuple[RealNumber] = None,
               max_iterations: int =100) -> tuple[Real_t]:
    """
    Find max allowed Y.

    EXAMPLES::


    """
    if not starting_Y:
        D = space.number_field().discriminant()
        if D == 5:
            starting_Y = 0.55
        elif D == 8:
            starting_Y = 0.31876213956380883
        else:
            starting_Y = 0.2
    if not isinstance(starting_Y, tuple):
        starting_Y = (starting_Y,) * space.number_field().absolute_degree()
    Y = starting_Y
    if not isinstance(M, tuple):
        M = ((-M, M),) * space.number_field().absolute_degree()
    if not Qs:
        Qs = get_Q_from_bounds(space.pullback(), M)
    log.debug(f"Trying: {Qs, Y}")
    for i in range(max_iterations):
        try:
            for id in space.group().ideal_cusp_representatives():
                get_pb_pts(space, Qs, id, Y)
        except ArithmeticError:
            log.debug("Arithmetic error, trying smaller Y")
            Y = tuple([y * 0.98 for y in Y])
        else:
            log.debug(f"Y={Y} is ok")
            break
    return Y

@mongo_cache()
def get_pb_pts(space: 'HilbertMaassFormSpace', Q: tuple, ideala: NumberFieldFractionalIdeal,
               Y: tuple, prec: int = 53, check: bool = True,
               use_symmetry: bool = False) -> tuple:
    """
    Get the list of points in the scaled lattice together with the corresponding pullbacks.
    """
    P = space.pullback()
    n = P.number_field().degree()
    CF = RealField(prec)
    ideala_matrix = matrix(P.basis_matrix_ideal(ideala))
    basis_matrix_m = ideala_matrix * diagonal_matrix([CF(1) / CF(2 * q) for q in Q])
    if use_symmetry:
        Q_combination = [range(1, Q[0] + 1), range(1-Q[1], Q[1]+1)]
    else:
        Q_combination = [range(1 - q, q + 1) for q in Q]
    zmpb = []
    zm = []
    log.info(f"Computing pullback for Q = {Q}, Y = {Y} idealamatrix={ideala_matrix}")
    for m in cartesian_product(Q_combination):
        xm = basis_matrix_m * vector(m)
        zm_elt = UpperHalfPlaneProductElement([(xm[i], Y[i]) for i in range(n)])
        zm.append(zm_elt)
        pbpt = P.reduce(zm_elt)
        if check and any(y <= Y[i] for i, y in enumerate(pbpt.imag())):
            raise ArithmeticError(f"Point {pbpt} has imaginary part smaller than {Y}. zm={zm_elt}")
        zmpb.append(pbpt)
    return zmpb, zm


def compute_coefficients(space: 'HilbertMaassFormSpace',
                         spectral_parameter: tuple[complex | ComplexNumber],
                         ideala: NumberFieldFractionalIdeal = None,
                         idealb: NumberFieldFractionalIdeal = None,
                         M: tuple[Integer_t] = None,
                         Y: tuple[float | RealNumber] = None,
                         returnV: bool = False,
                         set_coefficients: dict = None,
                         ncpus: int = 1) -> 'HilbertMaassCoefficients':
    r"""

    INPUT:

    - ``space`` -- Hilbert Maass form space
    - ``ideala``  -- NumberField Fractional Ideal corresponding to cusp.
    - ``idealb``  -- ? not used at the moment
    - ``s``       -- tuple of complex numbers - spectral parameter
    - ``M``       -- tuple of tuples of integers (or integer) - truncation bound
    - ``Y``       -- tuple of real numbers - height of sampling points
    - ``returnV`` -- boolean - return the matrix V - only used for debugging
    - ``ncpus``   -- number of cpus to use (need to be less than $SAGE_NUM_THREADS)
    - ``set_coefficients`` -- dictionary of coefficients to be set (default: None)

    EXAMPLES::

        sage: from hilbert_maass.all import HilbertMaassFormSpace
        sage: from hilbert_maass.modform.coefficients import compute_coefficients
        sage: M = (2,2)
        sage: s = CC(0.5,1.5), CC(0.5,1.5)
        sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
        sage: X = compute_coefficients(H, s, Y=(0.32, 0.32),M = 2); X
        Coefficients of a Hilbert Maass form with M=((-2, 2), (-2, 2)) and 1 cusp
        sage: X[(0,0)] == 0
        True
        sage: X[(1,1)] == 1
        True
        sage: X[(1,0)] # tol 1e-10
        0.0115350696785523 + 0.0000689324559423392*I
        sage: s = CC(0.5, 4.893781291438), CC(0.5, 4.893781291438)
        sage: H = HilbertMaassFormSpace(QuadraticField(5), cuspidal=True)
        sage: X = compute_coefficients(H, s, Y=(0.55, 0.55),M = 5); X
        Coefficients of a Hilbert Maass form with M=((-5, 5), (-5, 5)) and 1 cusp
        sage: X[(0,0)] == 0
        True
        sage: X[(1,-1)] == 1
        True
        sage: X[(0,-1)] # tol 1e-10
        1.00433823032911 + 0.000860350482435101*I
        sage: X[(0,1)] # tol 1e-10
        1.03339492015793 + 0.00102068814615568*I
    """
    complex_field = spectral_parameter[0].parent()
    ideala = ideala or space.number_field().ideal(1)
    idealb = idealb or space.number_field().ideal(1)
    smax = max(ceil(abs(s0)) for s0 in spectral_parameter)
    zpb, zm, Qs, M, Y = get_pb_pts_set_params(space, smax=smax,
                                              M=M, Y=Y, ideala=ideala)
    log.debug(f"M = {M}, Y = {Y}, Qs = {Qs}")
    use_iR = all((s-0.5).real() == 0 for s in spectral_parameter)
    matrixV = {}
    if ncpus > 1:
        matrix_arguments = []
        matrix_keys = []
    matrixV = setup_matrix(space, spectral_parameter,
                           ideala, idealb, Y, M, Qs, zpb, zm)
    RHS = {}
    normalisation = {}
    t_0 = (0,) * space.number_field().absolute_degree()
    n_0 = map_tuple_to_int(t_0, M)
    if not set_coefficients:
        set_coefficients = {}

    if space.is_cuspidal():
        # Set c(0)=0
        normalisation[n_0] = 0
        if not set_coefficients:
            # By default set c(delta)=1 where delta >>0 is generator of the index ideal.
            # tuple for delta
            ideala_dual = dual_ideal(ideala)
            delta = ideal_generator(ideala_dual)
            # Coordinate vector of delta
            t_1 = ideal_coordinates(ideala_dual, delta)
            n_1 = map_tuple_to_int(t_1, M)
            normalisation[n_1] = 1
    # Then update from set_coefficients
    for t, v in set_coefficients.items():
        normalisation[map_tuple_to_int(t, M)] = v
    if not space.is_cuspidal():
        W = w = (0,) * len(Y)
        for V in cartesian_product_from_M(M):
            V = tuple(V)
            v = dual_ideal_element(V, ideala)
            RHS[(V, W)] = matrix_element(spectral_parameter, Qs, v, w, zpb, zm, sgn='+')
            if V == W:
                tmp = bessel_prod(v, tuple(Y), spectral_parameter, sgn='+',
                                                use_iR=use_iR)
                RHS[(V, W)] = RHS[(V, W)] - tmp
    else:
        t_0 = (0,) * space.number_field().absolute_degree()
        for V in cartesian_product_from_M(M):
            RHS[(V, t_0)] = 0
            for n, v in normalisation.items():
                t = map_int_to_tuple(n, M)
                RHS[(V, t_0)] += matrixV[(V, t)] * v
    n = length_from_M(M)
    Vmat = [[
        matrixV[map_int_to_tuple(r, M), map_int_to_tuple(k, M)]
        for k in range(n)
    ] for r in range(n)]
    Vmat = matrix(complex_field, n, n, Vmat)
    RHSmat = [
        RHS[(map_int_to_tuple(k, M), t_0)] for k in range(n)
    ]
    RHSmat = matrix(complex_field, n, 1, RHSmat)
    if returnV and not space.is_cuspidal():
        return Vmat, RHSmat
    elif returnV:
        return Vmat
    # Find rows to delete
    delete_rows = list(normalisation.keys())
    delete_rows.sort()
    delete_rows = tuple(delete_rows)
    # Delete rows
    if delete_rows:
        # delete_rows = (n_0, n_1)
        Vmat = Vmat.delete_rows(delete_rows)
        RHSmat = RHSmat.delete_rows(delete_rows)
        Vmat = Vmat.delete_columns(delete_rows)
        X = Vmat.solve_right(-RHSmat)
        # Add back coefficients for 0 and 1
        skip_step = 0
        rows = []
        log.debug(f"normalisation = {normalisation}")
        for n in range(X.nrows() + len(normalisation)):
            if n in normalisation:
                rows.append((normalisation[n],))
                skip_step += 1
            else:
                rows.append(X.rows()[n-skip_step])
        # rows = X.rows()[0:n_0] + [(c0,)] + X.rows()[n_0:n_1 - 1] \
        #                        + [(c1,)] + X.rows()[n_1 - 1:]
        X = matrix(rows)
    else:
        X = Vmat.solve_right(-RHSmat)
    # Recreate the actual used set_coefficients dictionary
    set_coefficients_used = { map_int_to_tuple(k, M): v for k, v in normalisation.items() }
    return HilbertMaassCoefficients(X, M, spectral_parameter=spectral_parameter,
                                    space=space, coordinate_ideals=space.dual_ideals(),
                                    set_coefficients=set_coefficients_used,
                                    Y=Y)

def setup_matrix(space: 'HilbertMaassFormSpace',
                 spectral_parameter: tuple[complex | ComplexNumber],
                 ideala: NumberFieldFractionalIdeal,
                 idealb: NumberFieldFractionalIdeal,
                 Y,
                 M: tuple[tuple[Integer_t]],
                 Qs: tuple, zpb: list, zm: list) -> dict[tuple[tuple[Integer_t]]]:
    matrixV = {}
    use_iR = all((s-0.5).real() == 0 for s in spectral_parameter)
    bes_values = {}
    n = len(spectral_parameter)
    xms = [zmi.real() for zmi in zm]
    xpbs = [zpbi.real() for zpbi in zpb]
    ypbs = [zpbi.imag() for zpbi in zpb]
    # Pre-compute dual ideal elements
    dual_ideal_elements = {
        0: {},
        1: {}}
    for W in cartesian_product_from_M(M):
        dual_ideal_elements[0][W] = tuple(dual_ideal_element(W, ideala))
        dual_ideal_elements[1][W] = tuple(dual_ideal_element(W, idealb))
    # Pre-compute the Bessel product values
    for m, ympb in enumerate(ypbs):  # m in cartesian_product(Q_combination):
        bes_values[m] = {}
        for W in cartesian_product_from_M(M):
            # For cuspidal forms we don't need to compute the row corresponding to 0
            if space.is_cuspidal() and is_tuple_zero(W):
                bes_values[m][W] = 0
                continue
            w = dual_ideal_elements[1][W]
            if n == 2:
                bes = bessel_prod_dp2(w[0], w[1], ympb[0], ympb[1],
                                      spectral_parameter[0], spectral_parameter[1],
                                      sgn=0)
            else:
                bes = bessel_prod(tuple(w), tuple(ympb), spectral_parameter, sgn='-',
                                  use_iR=use_iR)
            exp_arg = (xpbs[m][0] * w[0], xpbs[m][1] * w[1])
            exp_val = exp_trace_prod_dp(exp_arg)
            bes_values[m][W] = bes * exp_val
    # Pre-compute the exponential values
    exp_values = {}
    for m, xm in enumerate(xms):
        exp_values[m] = {}
        for V in cartesian_product_from_M(M):
            v = dual_ideal_elements[0][V]
            exp_arg = (- xm[0] * v[0], - xm[1] * v[1])
            exp_values[m][V] = exp_trace_prod_dp(exp_arg)
    factor = prod(2 * q for q in Qs)
    for V in cartesian_product_from_M(M):
        for W in cartesian_product_from_M(M):
            # For cuspidal forms we don't need to compute the row corresponding to 0
            if space.is_cuspidal() and (is_tuple_zero(W) or is_tuple_zero(V)):
                matrixV[(V, W)] = 0
            else:
                summa = 0
                for m in range(len(xms)-1, -1, -1):
                    bes = bes_values[m][W]
                    term = bes * exp_values[m][V]
                    summa += term
                matrixV[(V, W)] = summa / factor
        v = dual_ideal_elements[0][V]
        if n == 2:
            bes = bessel_prod_dp2(v[0], v[1], Y[0], Y[1],
                                  spectral_parameter[0], spectral_parameter[1],
                                  sgn=0)
        else:
            bes = bessel_prod(v, tuple(Y), spectral_parameter, sgn='-', use_iR=use_iR)
        matrixV[(V, V)] = matrixV[(V, V)] - bes
    return matrixV
    RHS = {}
    normalisation = {}
    t_0 = (0,) * space.number_field().absolute_degree()
    n_0 = map_tuple_to_int(t_0, M)
    if space.is_cuspidal():
        # Set c(0)=0
        normalisation[n_0] = 0
    if not set_coefficients:
        set_coefficients = {}
        # By default set c(delta)=1 where delta >>0 is generator of the index ideal.
        # tuple for delta
        ideala_dual = dual_ideal(ideala)
        delta = ideal_generator(ideala_dual)
        # Coordinate vector of delta
        t_1 = ideal_coordinates(ideala_dual, delta)
        n_1 = map_tuple_to_int(t_1, M)
        normalisation[n_1] = 1
    # Then update from set_coefficients
    for t, v in set_coefficients.items():
        normalisation[map_tuple_to_int(t, M)] = v

    if not space.is_cuspidal():
        for V in cartesian_product_from_M(M):
            V = tuple(V)
            v = dual_ideal_element(V, ideala)
            W = w = (0,) * len(v)
            RHS[(V, W)] = matrix_element(spectral_parameter, Qs, v, w, zpb, zm, sgn='+')
            if V == W:
                RHS[(V, W)] = RHS[(V, W)] - bessel_prod(v, tuple(Y), spectral_parameter, sgn='+',
                                                        use_iR=use_iR)
    else:
        t_0 = (0,) * space.number_field().absolute_degree()
        for V in cartesian_product_from_M(M):
            RHS[(V, t_0)] = 0
            for n, v in normalisation.items():
                t = map_int_to_tuple(n, M)
                RHS[(V, t_0)] += matrixV[(V, t)] * v
    n = length_from_M(M)
    Vmat = [[
        matrixV[map_int_to_tuple(r, M), map_int_to_tuple(k, M)]
        for k in range(n)
    ] for r in range(n)]
    Vmat = matrix(complex_field, n, n, Vmat)
    if returnV:
        return Vmat
    RHSmat = [
        RHS[(map_int_to_tuple(k, M), t_0)] for k in range(n)
    ]
    RHSmat = matrix(complex_field, n, 1, RHSmat)
    # Find rows to delete
    delete_rows = list(normalisation.keys())
    delete_rows.sort()
    delete_rows = tuple(delete_rows)
    # Delete rows
    if delete_rows:
        # delete_rows = (n_0, n_1)
        Vmat = Vmat.delete_rows(delete_rows)
        RHSmat = RHSmat.delete_rows(delete_rows)
        Vmat = Vmat.delete_columns(delete_rows)
        X = Vmat.solve_right(-RHSmat)
        # Add back coefficients for 0 and 1
        skip_step = 0
        rows = []
        log.debug(f"normalisation = {normalisation}")
        for n in range(X.nrows() + len(normalisation)):
            if n in normalisation:
                rows.append((normalisation[n],))
                skip_step += 1
            else:
                rows.append(X.rows()[n-skip_step])
        # rows = X.rows()[0:n_0] + [(c0,)] + X.rows()[n_0:n_1 - 1] \
        #                        + [(c1,)] + X.rows()[n_1 - 1:]
        X = matrix(rows)
    else:
        X = Vmat.solve_right(-RHSmat)
    # Recreate the actual used set_coefficients dictionary
    set_coefficients_used = { map_int_to_tuple(k, M): v for k, v in normalisation.items() }
    return HilbertMaassCoefficients(X, M, spectral_parameter=spectral_parameter,
                                    space=space, coordinate_ideals=space.dual_ideals(),
                                    set_coefficients=set_coefficients_used,
                                    Y=Y)

def setup_matrix_symmetric(space: 'HilbertMaassFormSpace',
                 spectral_parameter: tuple[complex | ComplexNumber],
                 ideala: NumberFieldFractionalIdeal,
                 idealb: NumberFieldFractionalIdeal,
                 Y,
                 M: tuple[tuple[Integer_t]],
                 Qs: tuple, zpb: list, zm: list, sgn: str = '-',
                 eps: int = 1) -> dict[tuple[tuple[Integer_t]]]:
    """

    Args:
        space:
        spectral_parameter:
        ideala:
        idealb:
        Y:
        M:
        Qs:
        zpb:
        zm:
        sgn:
        eps: -- sign of symmetry

    Returns:

    """
    matrixV = {}
    use_iR = all((s-0.5).real() == 0 for s in spectral_parameter)
    bes_values = {}
    n = len(spectral_parameter)
    xms = [zmi.real() for zmi in zm]
    xpbs = [zpbi.real() for zpbi in zpb]
    ypbs = [zpbi.imag() for zpbi in zpb]
    # Pre-compute dual ideal elements
    dual_ideal_elements = {
        0: {},
        1: {}}
    list_of_coordinates = cartesian_product_from_M(M)
    for W in list_of_coordinates:
        dual_ideal_elements[0][W] = tuple(dual_ideal_element(W, ideala))
        dual_ideal_elements[1][W] = tuple(dual_ideal_element(W, idealb))
    # Pre-compute the Bessel product values
    for m, ympb in enumerate(ypbs):  # m in cartesian_product(Q_combination):
        bes_values[m] = {}
        for W in list_of_coordinates:
            # For cuspidal forms we don't need to compute the row corresponding to 0
            if space.is_cuspidal() and (all(x == 0 for x in W)):
                bes_values[m][W] = 0
                continue
            w = dual_ideal_elements[1][W]
            if n == 2:
                bes = bessel_prod_dp2(w[0], w[1], ympb[0], ympb[1],
                                      spectral_parameter[0], spectral_parameter[1],
                                      sgn == '+')
            else:
                bes = bessel_prod(tuple(w), tuple(ympb), spectral_parameter, sgn, use_iR=use_iR)
            exp_arg = (xpbs[m][0] * w[0], xpbs[m][1] * w[1])
            exp_val = exp_trace_prod_dp(exp_arg, symmetry=eps)
            bes_values[m][W] = bes * exp_val
    # Pre-compute the exponential values
    exp_values = {}
    for m, xm in enumerate(xms):
        exp_values[m] = {}
        for V in list_of_coordinates:
            v = dual_ideal_elements[0][V]
            exp_arg = (- xm[0] * v[0], - xm[1] * v[1])
            exp_values[m][V] = exp_trace_prod_dp(exp_arg, symmetry=eps)
    factor = prod(2 * q for q in Qs)
    for V in list_of_coordinates:
        for W in list_of_coordinates:
            # For cuspidal forms we don't need to compute the row corresponding to 0
            if space.is_cuspidal() and (all(x == 0 for x in W) or all(x == 0 for x in V)):
                matrixV[(V, W)] = 0
            else:
                summa = 0
                for m in range(len(xms)-1, -1, -1):
                    bes = bes_values[m][W]
                    term = bes * exp_values[m][V]
                    summa += term
                matrixV[(V, W)] = summa / factor
        v = dual_ideal_elements[0][V]
        if n == 2:
            bes = bessel_prod_dp2(v[0], v[1], Y[0], Y[1],
                                  spectral_parameter[0], spectral_parameter[1],
                                  sgn == '+')
        else:
            bes = bessel_prod(v, tuple(Y), spectral_parameter, sgn=sgn, use_iR=use_iR)
        matrixV[(V, V)] = matrixV[(V, V)] - bes
    return matrixV


def matrix_element(s: tuple, Q: tuple, v: tuple, w: tuple, zpb_v: list, zm_v: list,
                   sgn: str = '+') -> ComplexNumber:
    factor = prod(2 * q for q in Q)
    summa = 0
    sgn_bool = bool(sgn == '+')
    n = len(s)
    use_iR = all((si - 0.5).real() == 0 for si in s)
    for m, zm in enumerate(zm_v):  # m in cartesian_product(Q_combination):
        xm = zm.real()
        zmpb = zpb_v[m]
        ympb = zmpb.imag()
        xmpb = zmpb.real()
        if n == 2:
            bes = bessel_prod_dp2(w[0], w[1], ympb[0], ympb[1], s[0], s[1], sgn_bool)
        else:
            bes = bessel_prod(tuple(w), tuple(ympb), tuple(s), sgn, use_iR=use_iR)
        exp_arg = (xmpb[0] * w[0] - xm[0] * v[0], xmpb[1] * w[1] - xm[1] * v[1])
        exp_val = exp_trace_prod_dp(exp_arg)
        term = bes * exp_val
        summa += term
    return summa / factor
