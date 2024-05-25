from hilbert_maass.functions.functions_cy import bessel_prod_dp2, exp_trace_prod_dp
from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
from sage.categories.sets_cat import cartesian_product
from sage.functions.other import ceil, real
from sage.matrix.constructor import matrix, diagonal_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.complex_mpfr import ComplexNumber, ComplexField
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.rings.real_mpfr import RealNumber, RealField

from ..functions.functions import bessel_prod
from .coefficients import log, HilbertMaassCoefficients
from .utils import Integer_t, Real_t, mongo_cache, dual_ideal, \
    ideal_generator, cartesian_product_from_M, dual_ideal_element, length_from_M, is_tuple_zero, \
    get_Q_from_bounds, map_tuple_to_int, ideal_coordinates, map_int_to_tuple


@cached_method
def get_pb_pts_set_params(space: 'HilbertMaassFormSpace',
                          spectral_parameter: tuple[ComplexNumber] = None,
                          M: tuple[tuple[Integer_t]] = None,
                          Y: tuple = None,
                          smax: float | RealNumber = None,
                          ideala: NumberFieldFractionalIdeal = None,
                          use_symmetry: bool = False) -> tuple:
    spectral_parameter = spectral_parameter or (ComplexField(53)(0.5, 10), )
    CF = spectral_parameter[0].parent()
    if not spectral_parameter and not smax:
        raise ValueError("Need either spectral parameter or smax set.")
    if spectral_parameter and not smax:
        smax = max(abs(s) for s in spectral_parameter)
    if not isinstance(M, tuple):
        M0 = M or ceil((smax + 12) / (6.28318530717959) + 1)
        M = [(0 if use_symmetry else -M0, M0)]
        M = tuple(M + [(-M0, M0)] * (space.number_field().absolute_degree() - 1))
    if isinstance(M, tuple) and not isinstance(M[0], tuple):
        M = (M,) * space.number_field().absolute_degree()
    if use_symmetry and M[0][0] != 0:
        raise ValueError("If use_symmetry is True, M[0][0] must be 0")
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


@mongo_cache()
def get_pb_pts(space: 'HilbertMaassFormSpace', Q: tuple, ideala: NumberFieldFractionalIdeal,
               Y: tuple, prec: int = 53, check: bool = True,
               use_symmetry: bool = False) -> tuple:
    """
    Get the list of points in the scaled lattice together with the corresponding pullbacks.
    """
    P = space.pullback()
    n = P.number_field().degree()
    CF = RealField(prec+10)
    ideala_matrix = matrix(P.basis_matrix_ideal(ideala))
    basis_matrix_m = ideala_matrix * diagonal_matrix([CF(1) / CF(2 * q) for q in Q])
    if use_symmetry:
        Q_combination = [range(1, Q[0] + 1)] + [range(1 - q, q + 1) for q in Q[1:]]
    else:
        Q_combination = [range(1 - q, q + 1) for q in Q]
    zmpb = []
    zm = []
    log.info(f"Computing pullback for Q = {Q}, Y = {Y} idealamatrix={ideala_matrix}")
    half_vector = vector([CF(1)/CF(2)] * n)
    for m in cartesian_product(Q_combination):
        xm = basis_matrix_m * (vector(m) - half_vector)
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
        sage: from hilbert_maass.modform.compute_coefficients import compute_coefficients
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
        0.0108081473409434 + 3.77227326182247e-19*I
        sage: s = CC(0.5, 4.893781291438), CC(0.5, 4.893781291438)
        sage: H = HilbertMaassFormSpace(QuadraticField(5), cuspidal=True)
        sage: X = compute_coefficients(H, s, Y=(0.55, 0.55),M = 5); X
        Coefficients of a Hilbert Maass form with M=((-5, 5), (-5, 5)) and 1 cusp
        sage: X[(0,0)] == 0
        True
        sage: X[(1,-1)] == 1
        True
        sage: X[(0,-1)] # tol 1e-10
        1.00386877373754 + 2.59383671008929e-16*I
        sage: X[(0,1)] # tol 1e-10
        1.03204874572752 + 1.06392207755678e-15*I
    """
    if hasattr(spectral_parameter[0], 'parent'):
        complex_field = spectral_parameter[0].parent()
    else:
        complex_field = ComplexField(prec=53)
    ideala = ideala or space.number_field().ideal(1)
    idealb = idealb or space.number_field().ideal(1)
    smax = max(ceil(abs(s0)) for s0 in spectral_parameter)
    zpb, zm, Qs, M, Y = get_pb_pts_set_params(space, smax=smax,
                                              M=M, Y=Y, ideala=ideala)
    log.debug(f"M = {M}, Y = {Y}, Qs = {Qs}")
    use_iR = all(real(s - 0.5) == 0 for s in spectral_parameter)
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
        return Vmat, RHSmat
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
    use_iR = all(real(s-0.5) == 0 for s in spectral_parameter)
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
    factor = len(xms)
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


def compute_coefficients_symmetric(space: 'HilbertMaassFormSpace',
                         spectral_parameter: tuple[complex | ComplexNumber],
                         ideala: NumberFieldFractionalIdeal = None,
                         idealb: NumberFieldFractionalIdeal = None,
                         M: tuple[Integer_t] = None,
                         Y: tuple[float | RealNumber] = None,
                         returnV: bool = False,
                         set_coefficients: dict = None,
                         eps: int = 1,
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
        sage: from hilbert_maass.modform.compute_coefficients import compute_coefficients_symmetric
        sage: M = (2,2)
        sage: s = CC(0.5,1.5), CC(0.5,1.5)
        sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
        sage: X = compute_coefficients_symmetric(H, s, Y=(0.32, 0.32), eps=1, M = 2); X
        Coefficients of a Hilbert Maass form with M=((0, 2), (-2, 2)) and 1 cusp
        sage: X[(0,0)] == 0
        True
        sage: X[(1,1)] == 1
        True
        sage: X[(1,0)] # tol 1e-10
        0.0189237099737329
        sage: s = CC(0.5, 4.893781291438), CC(0.5, 4.893781291438)
        sage: H = HilbertMaassFormSpace(QuadraticField(5), cuspidal=True)
        sage: X = compute_coefficients_symmetric(H, s, Y=(0.55, 0.55),M = 5, eps=1); X
        Coefficients of a Hilbert Maass form with M=((0, 5), (-5, 5)) and 1 cusp
        sage: X[(0,0)] == 0
        True
        sage: X[(1,-1)] == 1
        True
        sage: X[(0,1)] # tol 1e-10
        0.989336294899374
        sage: X[(1,0)] # tol 1e-10
        0.797171647339938
    """
    complex_field = spectral_parameter[0].parent()
    ideala = ideala or space.number_field().ideal(1)
    idealb = idealb or space.number_field().ideal(1)
    smax = max(ceil(abs(s0)) for s0 in spectral_parameter)
    zpb, zm, Qs, M, Y = get_pb_pts_set_params(space, smax=smax,
                                              M=M, Y=Y, ideala=ideala,
                                              use_symmetry=True)
    log.debug(f"M = {M}, Y = {Y}, Qs = {Qs}")
    use_iR = all(real(s-0.5) == 0 for s in spectral_parameter)
    matrixV = {}
    if ncpus > 1:
        matrix_arguments = []
        matrix_keys = []
    list_of_coordinates_orig = cartesian_product_from_M(M)
    list_of_coordinates = []
    for W in list_of_coordinates_orig:
        if W[0] == 0 and W[1] < 0:
            continue
        list_of_coordinates.append(tuple(W))
    def map_int_to_tuple(index, limits):
        return list_of_coordinates[index]

    def map_tuple_to_int(tuple_index, limits):
        try:
            return list_of_coordinates.index(tuple_index)
        except IndexError:
            return None

    matrixV = setup_matrix_symmetric(space, spectral_parameter,
                           ideala, idealb, Y, M, Qs, zpb, zm, sgn='-',
                                     eps=eps)
    RHS = {}
    normalisation = {}
    t_0 = (0,) * space.number_field().absolute_degree()
    n_0 = map_tuple_to_int(t_0, M)
    log.debug(f"list={list_of_coordinates}")
    log.debug(f"t_0,n_0={t_0,n_0}")
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
        log.debug(f"t_0,n_0={t_1, n_1}")
        normalisation[n_1] = 1
    # Then update from set_coefficients
    for t, v in set_coefficients.items():
        normalisation[map_tuple_to_int(t, M)] = v
    if not space.is_cuspidal():
        for V in list_of_coordinates:
            V = tuple(V)
            v = dual_ideal_element(V, ideala)
            W = w = (0,) * len(v)
            RHS[(V, W)] = matrix_element(spectral_parameter, Qs, v, w, zpb, zm, sgn='+')
            if V == W:
                RHS[(V, W)] = RHS[(V, W)] - bessel_prod(v, tuple(Y), spectral_parameter, sgn='+',
                                                        use_iR=use_iR)
    else:
        t_0 = (0,) * space.number_field().absolute_degree()
        for V in list_of_coordinates:
            RHS[(V, t_0)] = 0
            for n, v in normalisation.items():
                t = map_int_to_tuple(n, M)
                RHS[(V, t_0)] += matrixV[(V, t)] * v
    n = len(list_of_coordinates) # length_from_M(M)
    Vmat = [[
        matrixV[map_int_to_tuple(r, M), map_int_to_tuple(k, M)]
        for k in range(n)
    ] for r in range(n)]
    Vmat = matrix(complex_field, n, n, Vmat)
    RHSmat = [
        RHS[(map_int_to_tuple(k, M), t_0)] for k in range(n)
    ]
    RHSmat = matrix(complex_field, n, 1, RHSmat)
    if returnV:
        return Vmat, RHSmat
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
    set_coefficients_used = {map_int_to_tuple(k, M): v for k, v in normalisation.items() }
    return HilbertMaassCoefficients(X, M, spectral_parameter=spectral_parameter,
                                    space=space, coordinate_ideals=space.dual_ideals(),
                                    set_coefficients=set_coefficients_used,
                                    index_tuples=[list_of_coordinates],
                                    check=False,
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
    if eps not in [-1, 1]:
        raise ValueError("eps must be -1 or 1")
    matrixV = {}
    use_iR = all(real(s-0.5) == 0 for s in spectral_parameter)
    bes_values = {}
    n = len(spectral_parameter)
    xms = [zmi.real() for zmi in zm]
    xpbs = [zpbi.real() for zpbi in zpb]
    ypbs = [zpbi.imag() for zpbi in zpb]
    # Pre-compute dual ideal elements
    dual_ideal_elements = {
        0: {},
        1: {}}
    list_of_coordinates_orig = cartesian_product_from_M(M)
    list_of_coordinates = []
    for W in list_of_coordinates_orig:
        if W[0] == 0 and W[1] < 0:
            continue
        list_of_coordinates.append(W)
    for W in list_of_coordinates:
        dual_ideal_elements[0][W] = tuple(dual_ideal_element(W, ideala))
        dual_ideal_elements[1][W] = tuple(dual_ideal_element(W, idealb))
    # Pre-compute the Bessel product values
    for m, ympb in enumerate(ypbs):  # m in cartesian_product(Q_combination):
        bes_values[m] = {}
        for W in list_of_coordinates:
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
                bes = bessel_prod(tuple(w), tuple(ympb), spectral_parameter, sgn='-', use_iR=use_iR)
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
    factor = len(xms) / 2
    for V in list_of_coordinates:
        for W in list_of_coordinates:
            # For cuspidal forms we don't need to compute the row corresponding to 0
            if space.is_cuspidal() and (is_tuple_zero(W) or is_tuple_zero(V)):
                matrixV[(V, W)] = 0
            else:
                summa = 0
                for m in range(len(xms)-1, -1, -1):
                # for m in range(len(xms)):
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


def matrix_element(s: tuple, Q: tuple, v: tuple, w: tuple, zpb_v: list, zm_v: list,
                   sgn: str = '+') -> ComplexNumber:
    factor = prod(2 * q for q in Q)
    summa = 0
    sgn_bool = bool(sgn == '+')
    n = len(s)
    use_iR = all(real(si - 0.5) == 0 for si in s)
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
