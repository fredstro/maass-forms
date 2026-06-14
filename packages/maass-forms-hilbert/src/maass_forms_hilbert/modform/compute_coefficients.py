from typing import TYPE_CHECKING

from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
from maass_form_core.functions.functions_cy import bessel_prod_dp2, exp_trace_prod_dp
from sage.categories.sets_cat import cartesian_product
from sage.functions.other import ceil, real
from sage.matrix.constructor import diagonal_matrix, matrix
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.infinity import Infinity
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.rings.real_mpfr import RR, RealField, RealNumber

from ..functions.functions import bessel_prod
from .coefficients import HilbertMaassCoefficients, log
from .utils import (
    Integer_t,
    Real_t,
    cartesian_product_from_M,
    dual_ideal,
    dual_ideal_basis_matrix,
    dual_ideal_element,
    get_Q_from_bounds,
    ideal_basis_matrix,
    ideal_coordinates,
    ideal_generator,
    is_tuple_zero,
    length_from_M,
    map_int_to_tuple,
    map_tuple_to_int,
    mongo_cache,
)

if TYPE_CHECKING:
    from .hilbert_maass_space import HilbertMaassFormSpace


# @cached_method
def get_pb_pts_set_params(
    space: "HilbertMaassFormSpace",
    spectral_parameter: tuple[ComplexNumber] = None,
    M: tuple[tuple[Integer_t]] = None,
    Y: tuple = None,
    Q_set: tuple[Integer_t, ...] = None,
    smax: float | RealNumber = None,
    prec: Integer_t = None,
    ideala: NumberFieldFractionalIdeal = None,
    use_symmetry: bool = False,
    use_shift: bool = False,
) -> tuple:
    r"""
    Compute pullback points and related parameters for coefficient computation.

    INPUT:

    - ``space`` -- Hilbert Maass form space
    - ``spectral_parameter`` -- tuple of complex numbers (optional)
    - ``M`` -- tuple of integer bounds (optional)
    - ``Y`` -- tuple of real numbers (optional)
    - ``Q_set`` -- tuple of integer Q values (optional)
    - ``smax`` -- maximum |s| (optional)
    - ``prec`` -- precision (optional)
    - ``ideala`` -- fractional ideal (optional)
    - ``use_symmetry`` -- use symmetry (default: False)
    - ``use_shift`` -- use shift (default: False)

    OUTPUT:

    - tuple (zpb, zm, Qs, M, Y)

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.compute_coefficients import get_pb_pts_set_params
        sage: from maass_forms_hilbert.modform.hilbert_maass_space import HilbertMaassFormSpace
        sage: space = HilbertMaassFormSpace(2)
        sage: s = (0.5+0.5j, 0.5+0.5j)
        sage: M = ((0, 1), (0,1))
        sage: Y = (0.32, 0.32)
        sage: Q = (2, 2)
        sage: t = get_pb_pts_set_params(space, spectral_parameter=s, M=M, Y=Y, Q_set=Q)
        sage: len(t) == 5
        True
        sage: len(t[0]) == len(t[1]) == 36
        True
        sage: t[2] == (3, 3)
        True
        sage: t[3] == M
        True
        sage: t[4] == Y
        True
    """
    if isinstance(spectral_parameter, tuple) and hasattr(spectral_parameter[0], "parent"):
        prec = spectral_parameter[0].parent().prec()
    elif not prec:
        prec = 53
    log.debug(
        f"in get_bb_pts_set_params: M={M}, Y={Y}, Q_set={Q_set}, smax={smax}, prec={prec}, smax={smax}"
    )
    if spectral_parameter is None and smax is None:
        raise ValueError("Need either spectral parameter or smax set.")
    if spectral_parameter and smax is None:
        smax = max(abs(s) for s in spectral_parameter)
    n = space.number_field().absolute_degree()
    if use_symmetry:
        # check if we need to swap the M-values
        ideala_dual = dual_ideal(ideala)
        delta = ideal_generator(ideala_dual)
        coords = ideal_coordinates(ideala_dual, delta)
        log.debug(f"coords={coords}")
        coord_signs = [-1 if x < 0 else 1 for x in coords]
        if coord_signs.count(-1) == 0 and coord_signs:
            coord_signs[0] = -1
    else:
        coord_signs = [-1] * n
    log.debug(f"coord_signs={coord_signs}")
    if not isinstance(M, tuple):
        M0 = M or ceil((smax + 12) / 6.28318530717959 + 1)
        M = tuple([(0 if coord_signs[i] == 1 else -M0, M0) for i in range(n)])
    if isinstance(M, tuple) and not isinstance(M[0], tuple):
        M = (M,) * space.number_field().absolute_degree()
    if use_symmetry and not all(m[0] == 0 for i, m in enumerate(M) if coord_signs[i] == 1):
        raise ValueError("If use_symmetry is True, M[0][0] must be 0")
    if not ideala:
        ideala = space.group().ideal_cusp_representatives()[0]
    Qs = get_Q_from_bounds(M)
    if Q_set and not isinstance(Q_set, (tuple, list)):
        Q_set = (Q_set,) * n
    if Q_set and min(Q_set) >= max(Qs):
        log.debug(f"Q_set {Q_set} is used instead of Qs={Qs}")
        Qs = Q_set
    # Try to find best value of Y
    Y = find_max_y(space, M, Qs=Qs, starting_Y=Y)
    Y = tuple([RealField(prec)(y) for y in Y])
    # We try with given Y and if it doesn't work we keep decreasing Y until it does.
    log.debug(
        f"2in get_bb_pts_set_params: M={M}, Y={Y}, Q_set={Q_set}, ideala={ideala}, "
        f"use_symmetry={use_symmetry}, use_shift={use_shift}"
    )

    try:
        zpb, zm = get_pb_pts(
            space,
            Qs,
            ideala,
            Y,
            use_symmetry=use_symmetry,
            use_shift=use_shift,
            prec=prec,
        )
    except ArithmeticError as e:
        msg = f"Could not find good pullback points. Error: {e}"
        log.debug(msg)
        raise ArithmeticError(msg) from e
    return zpb, zm, Qs, M, Y


@mongo_cache
def get_pb_pts(
    space: "HilbertMaassFormSpace",
    Q: tuple,
    ideala: NumberFieldFractionalIdeal,
    Y: tuple,
    prec: int = 53,
    check: bool = True,
    use_symmetry: bool = False,
    use_shift: bool = False,
) -> tuple:
    """
    Get the list of points in the scaled lattice together with the corresponding pullbacks.

    INPUT:

    - ``space`` -- Hilbert Maass form space
    - ``Q`` -- tuple of integers
    - ``ideala`` -- ideal
    - ``Y`` -- tuple of integers
    - ``prec`` -- precision
    - ``check`` -- if True, check that the points are in the lattice
    - ``use_symmetry`` -- if True, use symmetry
    - ``use_shift`` -- if True, use shift

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.compute_coefficients import get_pb_pts
        sage: from maass_forms_hilbert.modform.hilbert_maass_space import HilbertMaassFormSpace
        sage: space = HilbertMaassFormSpace(QuadraticField(5))
        sage: zpb, zm = get_pb_pts(space, (1,), space.number_field().ideal(1),
        ....:     (0.55,0.55))
        Traceback (most recent call last):
        ...
        ValueError: Q must be a tuple of length 2
        sage: zpb, zm = get_pb_pts(space, (1,1), space.number_field().ideal(1),
        ....:     (0.55,0.55))
        sage: RF = RealField(103)
        sage: zpb, zm = get_pb_pts(space, (1,1), space.number_field().ideal(1),
        ....:   (RF(0.55), RF(0.55)), prec=103)
        sage: len(zpb) == len(zm) == 4
        True
        sage: zm[0].prec()
        103
        sage: zpb[0][0].prec()
        103
    """
    P = space.pullback()
    n = P.number_field().degree()
    CF = RealField(prec)
    ideala_matrix = matrix(P.basis_matrix_ideal(ideala, prec=prec))
    if not isinstance(Q, (tuple, list)) or len(Q) != n:
        raise ValueError(f"Q must be a tuple of length {n}")
    basis_matrix_m = ideala_matrix * diagonal_matrix([CF(1) / CF(2 * q) for q in Q])
    if use_symmetry:
        Q_combination = [range(1, Q[0] + 1)] + [range(1 - q, q + 1) for q in Q[1:]]
    else:
        Q_combination = [range(1 - q, q + 1) for q in Q]
    zmpb = []
    zm = []
    log.info(f"Computing pullback for Q = {Q}, Y = {Y} idealamatrix={ideala_matrix}")
    if use_symmetry:
        use_shift = True
    if use_shift:
        half_vector = vector([CF(1) / CF(2)] * n)
    for m in cartesian_product(Q_combination):
        if use_shift:
            xm = basis_matrix_m * (vector(m) - half_vector)
        else:
            xm = basis_matrix_m * vector(m)
        zm_elt = UpperHalfPlaneProductElement([(xm[i], Y[i]) for i in range(n)], prec=prec)
        zm.append(zm_elt)
        pbpt = P.reduce(zm_elt)
        if check and any(y <= Y[i] for i, y in enumerate(pbpt.imag())):
            raise ArithmeticError(f"Point {pbpt} has imaginary part smaller than {Y}. zm={zm_elt}")
        zmpb.append(pbpt)
    return zmpb, zm


def compute_coefficients(
    space: "HilbertMaassFormSpace",
    spectral_parameter: tuple[complex | ComplexNumber],
    ideala: NumberFieldFractionalIdeal = None,
    idealb: NumberFieldFractionalIdeal = None,
    M: tuple[Integer_t] = None,
    Y: tuple[float | RealNumber] = None,
    Q: tuple[Integer_t] = None,
    returnV: bool = False,
    set_coefficients: dict = None,
    use_shift: bool = False,
    return_error_estimate=False,
    ncpus: int = 1,
) -> "HilbertMaassCoefficients":
    r"""

    INPUT:

    - ``space`` -- Hilbert Maass form space
    - ``ideala``  -- NumberField Fractional Ideal corresponding to cusp.
    - ``idealb``  -- ? not used at the moment
    - ``s``       -- tuple of complex numbers - spectral parameter
    - ``M``       -- tuple of tuples of integers (or integer) - truncation bound
    - ``Y``       -- tuple of real numbers - height of sampling points
    - ``Q``       -- tuple of integers - lattice size
    - ``returnV`` -- boolean - return the matrix V - only used for debugging
    - ``ncpus``   -- number of cpus to use (need to be less than $SAGE_NUM_THREADS)
    - ``set_coefficients`` -- dictionary of coefficients to be set (default: None)

    EXAMPLES::

        sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
        sage: from maass_forms_hilbert.modform.compute_coefficients import compute_coefficients
        sage: M = (2,2)
        sage: s = CC(0.5,5.12632439882674), CC(0.5,5.12632439882674)
        sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
        sage: X = compute_coefficients(H, s, Y=(0.32, 0.32),M = 2); X
        Coefficients of a Hilbert Maass form with M=((-2, 2), (-2, 2)) and 1 cusp
        sage: X[(0,0)] == 0
        True
        sage: X[(1,1)] == 1
        True
        sage: X[(0,1 )] # tol 1e-10
        -0.642563499325214 - 0.425608656688034*I
        sage: X = compute_coefficients(H, s, Y=(0.32, 0.32),M = 2, use_shift=True)
        sage: X[(0,0)] == 0
        True
        sage: X[(1,1)] == 1
        True
        sage: X[(1,0)] # abs tol 1e-10
        -0.685854898678222 + 8.65221920914551e-17*I
        sage: s = CC(0.5, 4.893781291438), CC(0.5, 4.893781291438)
        sage: H = HilbertMaassFormSpace(QuadraticField(5), cuspidal=True)
        sage: X = compute_coefficients(H, s, Y=(0.55, 0.55),M = 5); X
        Coefficients of a Hilbert Maass form with M=((-5, 5), (-5, 5)) and 1 cusp
        sage: X[(0,0)] == 0
        True
        sage: X[(1,-1)] == 1
        True
        sage: X[(0,-1)] # tol 1e-10
        1.00636761973816 + 0.00932996088880083*I
        sage: X[(0,1)] # tol 1e-10
        1.07092699890260 + 0.0127366413240728*I
        sage: X = compute_coefficients(H, s, Y=(0.55, 0.55),M = 5, use_shift=True)
        sage: X[(0,0)] == 0
        True
        sage: X[(1,-1)] == 1
        True
        sage: X[(0,-1)].real() # abs tol 1e-10
        1.00581013082305 + 4.00801061768804e-16*I
        sage: X[(0,1)].real() # abs tol 1e-10
        1.02379848615043 + 1.19775334209045e-15*I
    """
    if hasattr(spectral_parameter[0], "parent"):
        complex_field = spectral_parameter[0].parent()
    else:
        complex_field = ComplexField(prec=53)
    ideala = ideala or space.number_field().ideal(1)
    idealb = idealb or space.number_field().ideal(1)
    smax = max(ceil(abs(s0)) for s0 in spectral_parameter)
    prec = complex_field.precision()
    zpb, zm, Qs, M, Y = get_pb_pts_set_params(
        space,
        smax=smax,
        M=M,
        Y=Y,
        ideala=ideala,
        prec=prec,
        Q_set=Q,
        use_shift=use_shift,
    )
    use_iR = all(real(s - 0.5) == 0 for s in spectral_parameter)
    matrixV = {}
    if ncpus > 1:
        pass
    matrixV = setup_matrix(space, spectral_parameter, ideala, idealb, Y, M, Qs, zpb, zm)
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
            RHS[(V, W)] = matrix_element(spectral_parameter, Qs, v, w, zpb, zm, sgn="+")
            if V == W:
                tmp = bessel_prod(
                    v, tuple(Y), spectral_parameter, sgn="+", use_iR=use_iR, prec=prec
                )
                RHS[(V, W)] = RHS[(V, W)] - tmp
    else:
        t_0 = (0,) * space.number_field().absolute_degree()
        for V in cartesian_product_from_M(M):
            RHS[(V, t_0)] = 0
            for n, v in normalisation.items():
                t = map_int_to_tuple(n, M)
                RHS[(V, t_0)] += matrixV[(V, t)] * v
    n = length_from_M(M)
    Vmat = [
        [matrixV[map_int_to_tuple(r, M), map_int_to_tuple(k, M)] for k in range(n)]
        for r in range(n)
    ]
    Vmat = matrix(complex_field, n, n, Vmat)
    RHSmat = [RHS[(map_int_to_tuple(k, M), t_0)] for k in range(n)]
    RHSmat = matrix(complex_field, n, 1, RHSmat)
    if returnV and not space.is_cuspidal() or returnV:
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
                rows.append(X.rows()[n - skip_step])
        X = matrix(rows)
    else:
        X = Vmat.solve_right(-RHSmat)
    if return_error_estimate:
        err_est = (
            Vmat.norm(Infinity)
            * Vmat.inverse().norm(Infinity)
            * X.norm(Infinity)
            * RHSmat.norm(Infinity)
            * abs(M[0]) ** (n / 2)
        )
    else:
        err_est = None
    # Recreate the actual used set_coefficients dictionary
    set_coefficients_used = {map_int_to_tuple(k, M): v for k, v in normalisation.items()}
    return HilbertMaassCoefficients(
        X,
        M,
        spectral_parameter=spectral_parameter,
        space=space,
        coordinate_ideals=space.dual_ideals(),
        set_coefficients=set_coefficients_used,
        Y=Y,
        Q=Qs,
        error_est=err_est,
    )


def setup_matrix(
    space: "HilbertMaassFormSpace",
    spectral_parameter: tuple[complex | ComplexNumber],
    ideala: NumberFieldFractionalIdeal,
    idealb: NumberFieldFractionalIdeal,
    Y,
    M: tuple[tuple[Integer_t]],
    Qs: tuple,
    zpb: list,
    zm: list,
    Mcol: tuple[tuple[Integer_t]] = None,
    sub_diagonal: bool = True,
) -> dict[tuple[tuple[Integer_t]]]:
    """
    Setup the matrix V

    INPUT::

    - ``space`` -- Hilbert Maass form space
    - ``spectral_parameter`` -- tuple of complex or real numbers
    - ``ideala`` -- fractional ideal
    - ``idealb`` -- fractional ideal
    - ``Y`` -- tuple of real numbers
    - ``M`` -- tuple of integers giving the limits used to map the indices of the matrix to tuples
    - ``Qs`` -- tuple of positive integers. >= M+2
    - ``zpb`` -- list of complex numbers
    - ``zm`` -- list of complex numbers
    - ``Mcol`` -- same as ``M`` but for column indices if different from row indices

    OUTPUT:

    - dictionary representing the matrix V

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.compute_coefficients import setup_matrix
        sage: from maass_forms_hilbert.modform.hilbert_maass_space import HilbertMaassFormSpace
        sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
        sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
        sage: s = CC(0.5, 5.12632439882674), CC(0.5, 5.12632439882674)
        sage: ideala = H.dual_ideals()[0]
        sage: idealb = H.dual_ideals()[0]
        sage: Y = (0.32, 0.32)
        sage: M = ((-1, 1), (-1, 1))
        sage: Qs = (3, 3)
        sage: zpb = [[1.81818181818182*I, 1.81818181818182*I],
        ....: [-0.772673684476164 + 0.574707540859232*I, -0.158407197032179 + 1.38193905174673*I],
        ....: [0.0950226244343892 + 0.995475113122172*I, 0.0950226244343892 + 0.995475113122172*I],
        ....: [-0.841592802967824 + 1.38193905174674*I, -0.227326315523835 + 0.574707540859230*I]]
        sage: zpb = [UpperHalfPlaneProductElement(z) for z in zpb]
        sage: zm = [[0.550000000000000*I, 0.550000000000000*I],
        ....: [-0.809016994374945 + 0.550000000000000*I, 0.309016994374948 + 0.550000000000000*I],
        ....: [0.500000000000000 + 0.550000000000000*I, 0.500000000000000 + 0.550000000000000*I],
        ....: [-0.309016994374945 + 0.550000000000000*I, 0.809016994374947 + 0.550000000000000*I]]
        sage: zm = [UpperHalfPlaneProductElement(z) for z in zm]
        sage: A = setup_matrix(H, s, ideala, idealb, Y=Y, M=M, Qs=Qs, zpb=zpb, zm=zm)
        sage: len(A.keys())
        81
        sage: A[((1, 0), (0, 0))]
        0
        sage: A[(0, 0), (0, 0)]
        0
        sage: A[(-1, 0), (0, -1)] # tol 1e-10
        (3.7447884858466083e-05-0.0004907274543417816j)
    """
    matrixV = {}
    use_iR = all(real(s - 0.5) == 0 for s in spectral_parameter)
    bes_values = {}
    n = len(spectral_parameter)
    xms = [zmi.real() for zmi in zm]
    xpbs = [zpbi.real() for zpbi in zpb]
    ypbs = [zpbi.imag() for zpbi in zpb]
    if hasattr(spectral_parameter[0], "parent"):
        prec = spectral_parameter[0].parent().prec()
    else:
        prec = 53
    # Pre-compute dual ideal elements
    dual_ideal_elements = {0: {}, 1: {}}
    if not Mcol:
        Mcol = M
    for W in cartesian_product_from_M(Mcol):
        dual_ideal_elements[0][W] = tuple(dual_ideal_element(W, ideala))
        dual_ideal_elements[1][W] = tuple(dual_ideal_element(W, idealb))
    for W in cartesian_product_from_M(M):
        if W not in dual_ideal_elements[0]:
            dual_ideal_elements[0][W] = tuple(dual_ideal_element(W, ideala))
        if W not in dual_ideal_elements[1]:
            dual_ideal_elements[1][W] = tuple(dual_ideal_element(W, idealb))
    # Pre-compute the Bessel product values
    for m, ympb in enumerate(ypbs):  # m in cartesian_product(Q_combination):
        bes_values[m] = {}
        for W in cartesian_product_from_M(Mcol):
            # For cuspidal forms we don't need to compute the row corresponding to 0
            if space.is_cuspidal() and is_tuple_zero(W):
                bes_values[m][W] = 0
                continue
            w = dual_ideal_elements[1][W]
            if n == 2 and prec == 53:
                bes = bessel_prod_dp2(
                    w[0],
                    w[1],
                    ympb[0],
                    ympb[1],
                    spectral_parameter[0],
                    spectral_parameter[1],
                    sgn=0,
                )
            else:
                bes = bessel_prod(
                    tuple(w),
                    tuple(ympb),
                    spectral_parameter,
                    sgn="-",
                    use_iR=use_iR,
                    prec=prec,
                )
            exp_arg = (xpbs[m][0] * w[0], xpbs[m][1] * w[1])
            exp_val = exp_trace_prod_dp(exp_arg)
            bes_values[m][W] = bes * exp_val
    # Pre-compute the exponential values
    exp_values = {}
    for m, xm in enumerate(xms):
        exp_values[m] = {}
        for V in cartesian_product_from_M(M):
            v = dual_ideal_elements[0][V]
            exp_arg = (-xm[0] * v[0], -xm[1] * v[1])
            exp_values[m][V] = exp_trace_prod_dp(exp_arg)
    factor = len(xms)
    for V in cartesian_product_from_M(M):
        for W in cartesian_product_from_M(Mcol):
            # For cuspidal forms we don't need to compute the row corresponding to 0
            if space.is_cuspidal() and (is_tuple_zero(W) or is_tuple_zero(V)):
                matrixV[(V, W)] = 0
            else:
                summa = 0
                for m in range(len(xms) - 1, -1, -1):
                    bes = bes_values[m][W]
                    term = bes * exp_values[m][V]
                    summa += term
                matrixV[(V, W)] = summa / factor
        v = dual_ideal_elements[0][V]
        if (V, V) in matrixV and sub_diagonal:
            if space.is_cuspidal() and is_tuple_zero(V):
                continue
            if n == 2:
                bes = bessel_prod_dp2(
                    v[0],
                    v[1],
                    Y[0],
                    Y[1],
                    spectral_parameter[0],
                    spectral_parameter[1],
                    sgn=0,
                )
            else:
                bes = bessel_prod(v, tuple(Y), spectral_parameter, sgn="-", use_iR=use_iR)
            matrixV[(V, V)] = matrixV[(V, V)] - bes
    return matrixV


def compute_coefficients_symmetric(
    space: "HilbertMaassFormSpace",
    spectral_parameter: tuple[complex | ComplexNumber],
    ideala: NumberFieldFractionalIdeal = None,
    idealb: NumberFieldFractionalIdeal = None,
    M: tuple[Integer_t] = None,
    Y: tuple[float | RealNumber] = None,
    returnV: bool = False,
    set_coefficients: dict = None,
    eps: int = 1,
    ncpus: int = 1,
) -> "HilbertMaassCoefficients":
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

        sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
        sage: from maass_forms_hilbert.modform.compute_coefficients import compute_coefficients_symmetric
        sage: M = (2,2)
        sage: s = CC(0.5,5.12632439882674), CC(0.5,5.12632439882674)
        sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
        sage: X = compute_coefficients_symmetric(H, s, Y=(0.32, 0.32), eps=1, M = 2); X
        Coefficients of a Hilbert Maass form with M=((-2, 2), (0, 2)) and 1 cusp
        sage: X[(0,0)] == 0
        True
        sage: X[(1,1)] == 1
        True
        sage: X[(1,0)] # tol 1e-10
        -1.02751145868017
        sage: s = CC(0.5, 4.893781291438), CC(0.5, 4.893781291438)
        sage: H = HilbertMaassFormSpace(QuadraticField(5), cuspidal=True)
        sage: X = compute_coefficients_symmetric(H, s, Y=(0.55, 0.55),M = 5, eps=1); X # long time
        Coefficients of a Hilbert Maass form with M=((0, 5), (-5, 5)) and 1 cusp
        sage: X[(0,0)] == 0  # long time
        True
        sage: X[(1,-1)] == 1  # long time
        True
        sage: X[(0,1)] # tol 1e-10  long time
        0.996248597284956
        sage: X[(1,0)] # tol 1e-10 long time
        0.801225372213468
    """
    complex_field = spectral_parameter[0].parent()
    ideala = ideala or space.number_field().ideal(1)
    idealb = idealb or space.number_field().ideal(1)
    smax = max(ceil(abs(s0)) for s0 in spectral_parameter)
    zpb, zm, Qs, M, Y = get_pb_pts_set_params(
        space, smax=smax, M=M, Y=Y, ideala=ideala, use_symmetry=True
    )
    log.debug(f"M = {M}, Y = {Y}, Qs = {Qs}")
    use_iR = all(real(s - 0.5) == 0 for s in spectral_parameter)
    matrixV = {}
    if ncpus > 1:
        pass
    list_of_coordinates_orig = cartesian_product_from_M(M)
    list_of_coordinates = []
    m0_0 = 0 if M[0][0] == 0 else 1
    m0_1 = 1 if m0_0 == 0 else 0
    for W in list_of_coordinates_orig:
        if W[m0_0] == 0 and W[m0_1] < 0:
            continue
        list_of_coordinates.append(tuple(W))

    def map_int_to_tuple(index, limits):
        return list_of_coordinates[index]

    def map_tuple_to_int(tuple_index, limits):
        try:
            return list_of_coordinates.index(tuple_index)
        except IndexError:
            return None

    matrixV = setup_matrix_symmetric(
        space, spectral_parameter, ideala, idealb, Y, M, Qs, zpb, zm, sgn="-", eps=eps
    )
    RHS = {}
    normalisation = {}
    t_0 = (0,) * space.number_field().absolute_degree()
    n_0 = map_tuple_to_int(t_0, M)
    log.debug(f"list={list_of_coordinates}")
    log.debug(f"t_0,n_0={t_0, n_0}")
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
            RHS[(V, W)] = matrix_element(spectral_parameter, Qs, v, w, zpb, zm, sgn="+")
            if V == W:
                RHS[(V, W)] = RHS[(V, W)] - bessel_prod(
                    v, tuple(Y), spectral_parameter, sgn="+", use_iR=use_iR
                )
    else:
        t_0 = (0,) * space.number_field().absolute_degree()
        for V in list_of_coordinates:
            RHS[(V, t_0)] = 0
            for n, v in normalisation.items():
                t = map_int_to_tuple(n, M)
                RHS[(V, t_0)] += matrixV[(V, t)] * v
    n = len(list_of_coordinates)  # length_from_M(M)
    Vmat = [
        [matrixV[map_int_to_tuple(r, M), map_int_to_tuple(k, M)] for k in range(n)]
        for r in range(n)
    ]
    Vmat = matrix(complex_field, n, n, Vmat)
    RHSmat = [RHS[(map_int_to_tuple(k, M), t_0)] for k in range(n)]
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
                rows.append(X.rows()[n - skip_step])
        # rows = X.rows()[0:n_0] + [(c0,)] + X.rows()[n_0:n_1 - 1] \
        #                        + [(c1,)] + X.rows()[n_1 - 1:]
        X = matrix(rows)
    else:
        X = Vmat.solve_right(-RHSmat)
    # Recreate the actual used set_coefficients dictionary
    set_coefficients_used = {map_int_to_tuple(k, M): v for k, v in normalisation.items()}
    return HilbertMaassCoefficients(
        X,
        M,
        spectral_parameter=spectral_parameter,
        space=space,
        coordinate_ideals=space.dual_ideals(),
        set_coefficients=set_coefficients_used,
        index_tuples=[list_of_coordinates],
        check=False,
        Y=Y,
    )


def setup_matrix_symmetric(
    space: "HilbertMaassFormSpace",
    spectral_parameter: tuple[complex | ComplexNumber],
    ideala: NumberFieldFractionalIdeal,
    idealb: NumberFieldFractionalIdeal,
    Y,
    M: tuple[tuple[Integer_t]],
    Qs: tuple,
    zpb: list,
    zm: list,
    sgn: str = "-",
    eps: int = 1,
) -> dict[tuple[tuple[Integer_t]]]:
    """
    Setup the symmetric matrix V

    INPUT:

    - ``space`` -- Hilbert Maass form space
    - ``spectral_parameter`` -- tuple of complex or real numbers
    - ``ideala`` -- fractional ideal
    - ``idealb`` -- fractional ideal
    - ``Y`` -- tuple of real numbers
    - ``M`` -- tuple of integers giving the limits used to map the indices of the matrix to tuples
    - ``Qs`` -- tuple of positive integers. >= M+2
    - ``zpb`` -- list of complex numbers
    - ``zm`` -- list of complex numbers
    - ``sgn`` -- sign string, either "+" or "-"
    - ``eps`` -- sign of symmetry, either -1 or 1

    OUTPUT:

    - dictionary representing the symmetric matrix V

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.compute_coefficients import setup_matrix_symmetric
        sage: # This function is complex and used internally for symmetric matrix computation
        sage: # It takes many parameters and returns a dictionary structure
        sage: setup_matrix_symmetric.__name__
        'setup_matrix_symmetric'
    """
    if eps not in [-1, 1]:
        raise ValueError("eps must be -1 or 1")
    matrixV = {}
    use_iR = all(real(s - 0.5) == 0 for s in spectral_parameter)
    bes_values = {}
    n = len(spectral_parameter)
    xms = [zmi.real() for zmi in zm]
    xpbs = [zpbi.real() for zpbi in zpb]
    ypbs = [zpbi.imag() for zpbi in zpb]
    # Pre-compute dual ideal elements
    dual_ideal_elements = {0: {}, 1: {}}
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
                bes = bessel_prod_dp2(
                    w[0],
                    w[1],
                    ympb[0],
                    ympb[1],
                    spectral_parameter[0],
                    spectral_parameter[1],
                    sgn=0,
                )
            else:
                bes = bessel_prod(tuple(w), tuple(ympb), spectral_parameter, sgn="-", use_iR=use_iR)
            exp_arg = (xpbs[m][0] * w[0], xpbs[m][1] * w[1])
            exp_val = exp_trace_prod_dp(exp_arg, symmetry=eps)
            bes_values[m][W] = bes * exp_val
    # Pre-compute the exponential values
    exp_values = {}
    for m, xm in enumerate(xms):
        exp_values[m] = {}
        for V in list_of_coordinates:
            v = dual_ideal_elements[0][V]
            exp_arg = (-xm[0] * v[0], -xm[1] * v[1])
            exp_values[m][V] = exp_trace_prod_dp(exp_arg, symmetry=eps)
    factor = len(xms) / 2
    for V in list_of_coordinates:
        for W in list_of_coordinates:
            # For cuspidal forms we don't need to compute the row corresponding to 0
            if space.is_cuspidal() and (is_tuple_zero(W) or is_tuple_zero(V)):
                matrixV[(V, W)] = 0
            else:
                summa = 0
                for m in range(len(xms) - 1, -1, -1):
                    # for m in range(len(xms)):
                    bes = bes_values[m][W]
                    term = bes * exp_values[m][V]
                    summa += term
                matrixV[(V, W)] = summa / factor
        v = dual_ideal_elements[0][V]
        if n == 2:
            bes = bessel_prod_dp2(
                v[0],
                v[1],
                Y[0],
                Y[1],
                spectral_parameter[0],
                spectral_parameter[1],
                sgn=0,
            )
        else:
            bes = bessel_prod(v, tuple(Y), spectral_parameter, sgn="-", use_iR=use_iR)
        matrixV[(V, V)] = matrixV[(V, V)] - bes
    return matrixV


def matrix_element(
    s: tuple, Q: tuple, v: tuple, w: tuple, zpb_v: list, zm_v: list, sgn: str = "+"
) -> ComplexNumber:
    r"""
    Compute a matrix element for the coefficient system.

    INPUT:

    - ``s`` -- tuple of complex numbers (spectral parameter)
    - ``Q`` -- tuple of integer Q values
    - ``v`` -- tuple (row index)
    - ``w`` -- tuple (column index)
    - ``zpb_v`` -- list of pullback points
    - ``zm_v`` -- list of lattice points
    - ``sgn`` -- sign, '+' or '-'

    OUTPUT:

    - ComplexNumber, the matrix element

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.compute_coefficients import matrix_element
        sage: # Typical usage (requires setup):
        sage: # matrix_element(s, Q, v, w, zpb_v, zm_v)
    """
    factor = prod(2 * q for q in Q)
    summa = 0
    sgn_bool = bool(sgn == "+")
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
def find_max_y(
    space: "HilbertMaassFormSpace",
    M: tuple[tuple[Integer_t]],
    Qs: tuple[tuple[Integer_t]] = None,
    starting_Y: tuple[RealNumber] = None,
    max_iterations: int = 100,
) -> tuple[Real_t]:
    """
    Find max allowed Y.

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.compute_coefficients import find_max_y
        sage: from maass_forms_hilbert.modform.hilbert_maass_space import HilbertMaassFormSpace
        sage: H = HilbertMaassFormSpace(5)
        sage: find_max_y(H, ((-1, 1),(-1,1)))
        (0.55, 0.55)

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
        Qs = get_Q_from_bounds(M)
    log.debug(f"Trying: {Qs, Y}")
    for _i in range(max_iterations):
        ok_all_cusps = True
        for id in space.group().ideal_cusp_representatives():
            zp, zm = get_pb_pts(space, Qs, id, Y, check=False)
            miny = min(min(z.imag()) for z in zp)
            if miny <= Y[0] + 1e-10:
                log.debug(f"Y={Y} is too large")
                Y = (miny * 0.99,) * len(Y)
                ok_all_cusps = False
                break
        if ok_all_cusps:
            return Y
    raise ArithmeticError("Could not find max Y")


def actual_tail_sum(space, M: Integer_t, Y: Real_t | None = None, Q: Integer_t = 100):
    """
    Return the estimate of the truncated lattice sum from Equation 5.2 of the article.

    INPUT:

    - ``space`` - HilbertMaassFormSpace
    - ``M`` - Integer; truncation bound, the sum is over coordinates with
      ``||x||_inf > M``
    - ``Y`` - Real (optional); if ``None``, auto-tuned via :func:`find_max_y`
    - ``Q`` -  Integer (default: 100); lattice cutoff. Coordinates are drawn from
      ``[-Q, Q]^n``, so ``Q`` must be large enough that the contribution from
      ``||x||_inf > Q`` is negligible relative to the desired precision. In
      practice ``Q`` should be chosen so that the exponential decay dominates
      well before the cutoff; ``Q >> M`` is usually sufficient.

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.compute_coefficients import actual_tail_sum
        sage: from maass_forms_hilbert.modform.hilbert_maass_space import HilbertMaassFormSpace
        sage: space = HilbertMaassFormSpace(2)
        sage: actual_tail_sum(space, 20, 0.55) # tol 6e-15
        3.44069786315208e-21
        sage: space = HilbertMaassFormSpace(5)
        sage: actual_tail_sum(space, 30, 0.31) # tol 3e-15
        9.44541005162014e-11
    """
    n = space.number_field().absolute_degree()
    if Y is None:
        Y = find_max_y(space, ((-M, M),) * n)[0]
    ideala = space.number_field().ideal(1)
    two_pi_Y = 2 * RR.pi() * Y
    B_dual = dual_ideal_basis_matrix(ideala)
    coordinates = cartesian_product_from_M(((-Q, Q),) * n)
    return sum(
        (-(B_dual * vector(x)).norm(1) * two_pi_Y).exp()
        for x in coordinates
        if max(abs(xi) for xi in x) > M
    )


def tail_sum_formula_gen(space, M: Integer_t, Y: Real_t | None = None):
    """
    Return value of error estimate for general number field from Equation 5.3 of the article.

    INPUT:

    - ``space`` - HilbertMaassFormSpace
    - ``M`` - Integer
    - ``Y`` - Real

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.compute_coefficients import tail_sum_formula_gen
        sage: from maass_forms_hilbert.modform.hilbert_maass_space import HilbertMaassFormSpace
        sage: space = HilbertMaassFormSpace(2)
        sage: tail_sum_formula_gen(space, 20, 0.55)
        5.4188622597365765e-12
        sage: space = HilbertMaassFormSpace(5)
        sage: tail_sum_formula_gen(space, 30, 0.31)
        1.656644109567389e-08
    """
    n = space.number_field().absolute_degree()
    if Y is None:
        Y = find_max_y(space, ((-M, M),) * n)[0]
    delta = RR.pi() * 2 * Y
    ideala = space.number_field().ideal(1)
    B = ideal_basis_matrix(ideala).transpose().norm(1)
    B1 = 2 * B / delta
    val = max(B1, 1)
    return n * val * (M + 2) ** (n - 1) * (-delta * B**-1 * (M + 1)).exp()


def tail_sum_formula_quad(space, M: Integer_t, Y: Real_t | None = None):
    """
    Return value of error estimate for quadratic field from Equation 5.4 of the article.

    INPUT:

    - ``space`` - HilbertMaassFormSpace
    - ``M`` - Integer
    - ``Y`` - Real

    EXAMPLES::

        sage: from maass_forms_hilbert.modform.compute_coefficients import tail_sum_formula_quad
        sage: from maass_forms_hilbert.modform.hilbert_maass_space import HilbertMaassFormSpace
        sage: space = HilbertMaassFormSpace(2)
        sage: tail_sum_formula_quad(space, 20, 0.55)
        2.19166092241365e-19
        sage: space = HilbertMaassFormSpace(5)
        sage: tail_sum_formula_quad(space, 30, 0.31)
        3.80320923718002e-9
    """
    n = space.number_field().absolute_degree()
    if Y is None:
        Y = find_max_y(space, ((-M, M),) * n)[0]
    D = space.number_field().discriminant()
    d = D / 4 if D % 4 == 0 else D
    delta = RR.pi() * 2 * Y
    val = max(delta**-2, 1)
    return 4 * RR(d).sqrt() * (3 * M + 5) * val * (-delta * M / RR(d).sqrt()).exp()
