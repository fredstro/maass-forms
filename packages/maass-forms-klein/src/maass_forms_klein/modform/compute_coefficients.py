import logging

import numpy
from maass_form_core.utils.json_converters import dict_from_json
from sage.all import CC
from sage.functions.other import ceil, imag, real
from sage.matrix.constructor import matrix
from sage.matrix.matrix0 import Matrix
from sage.misc.cachefunc import cached_function as cached_function_default
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.real_mpfr import RealField
from sage.structure.sequence import Sequence

from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement
from maass_forms_klein.hyperbolic_space.utils import get_lattice_values
from maass_forms_klein.modform.coefficients import KleinianMaassFormCoefficients
from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
from maass_forms_klein.modform.utils import Complex_t, Integer_t, Real_t, bessel_function, get_prec

log = logging.getLogger(__name__)

try:
    from comp_manager.decorators import mongo_cache

    cached_function = mongo_cache()
except ImportError:
    cached_function = cached_function_default


@cached_function
def setup_matrix(
    space,
    spectral_parameter: Complex_t | Real_t,
    Y: Real_t,
    M: tuple[Integer_t] | Integer_t,
    zpb: tuple,
    zm: tuple,
) -> Matrix:
    """
    Compute the matrix to solve.

    INPUT:

    - ``space`` --  KleinianGroup
    - ``spectral_parameter`` --  Complex_t | Real_t
    - ``Y`` --  Real_t
    - ``M`` --  tuple[Integer_t] | Integer_t
    - ``Qs`` --  tuple
    - ``zpb`` --  list
    - ``zm`` --  list
    - ``dual_lattice_basis`` --  tuple[Real_t]

    EXAMPLES:
        sage: from maass_forms_klein.modform.compute_coefficients import setup_matrix
        sage: from maass_forms_klein.all import KleinianMaassFormSpace
        sage: from maass_forms_klein.modform.compute_coefficients import get_pb_pts
        sage: space = KleinianMaassFormSpace(-4)
        sage: Y = 0.5; M = 1; Q = 3
        sage: zpb, zm = get_pb_pts(space.group(), Q=Q, Y=Y)
        sage: s = CC(0.5,6.62211934)
        sage: V = setup_matrix(space, spectral_parameter=s, M=M, Y=Y, zpb=tuple(zpb), zm=tuple(zm))
        sage: V.ncols() == 9 and V.nrows() == 9
        True


    """

    bes_values = {}
    xms = [zmi.z() for zmi in zm]
    assert {zmi.y() for zmi in zm} == {Y}
    xpbs = [zpbi.z() for zpbi in zpb]
    ypbs = [zpbi.y() for zpbi in zpb]
    prec = get_prec(spectral_parameter)
    CF = ComplexField(prec)
    RF = RealField(prec)
    if real(spectral_parameter) == 0.5 and isinstance(spectral_parameter, Complex_t):
        spectral_parameter = RF(imag(spectral_parameter))
    dual_lattice_values = space.group().dual_translation_lattice_vectors(M)
    matrixV = matrix(CF, len(dual_lattice_values), len(dual_lattice_values))

    # Pre-compute the Bessel product values
    twopi = 2 * CF.pi()
    twopii = CF(0, twopi)
    for m, ympb in enumerate(ypbs):
        bes_values[m] = {}
        for nw, w in enumerate(dual_lattice_values):
            # For cuspidal forms we don't need to compute the row corresponding to 0
            if space.is_cuspidal() and abs(w) == 0:
                bes_values[m][nw] = 0
                continue
            bes = ympb * bessel_function(abs(w), twopi * ympb, spectral_parameter)
            exp_arg = xpbs[m][0] * w[0] - xpbs[m][1] * w[1]
            exp_val = (twopii * exp_arg).exp()
            bes_values[m][nw] = bes * exp_val
    # Pre-compute the exponential values
    exp_values = {}
    for m, xm in enumerate(xms):
        exp_values[m] = {}
        for nv, v in enumerate(dual_lattice_values):
            exp_arg = xm[0] * v[0] - xm[1] * v[1]
            exp_values[m][nv] = (-twopii * exp_arg).exp()
    factor = len(xms)
    for nv, v in enumerate(dual_lattice_values):
        if abs(v) == 0 and space.is_cuspidal():
            continue
        for nw, w in enumerate(dual_lattice_values):
            # For cuspidal forms we don't need to compute the row corresponding to 0
            if space.is_cuspidal() and abs(w) == 0:
                matrixV[nv, nw] = 0
            else:
                summa = 0
                for m in range(len(xms)):
                    summa += bes_values[m][nw] * exp_values[m][nv]
                matrixV[nv, nw] = summa / factor
        bes = Y * bessel_function(abs(v), Y * twopi, spectral_parameter)
        matrixV[nv, nv] = matrixV[nv, nv] - bes
    return matrixV


@cached_method()
def find_max_y(
    space: KleinianMaassFormSpace,
    M: Integer_t,
    Q: Integer_t,
    starting_Y: Real_t = None,
    max_iterations: int = 100,
) -> Real_t:
    r"""
    Find the maximum allowed height Y for pullback computations.

    Starting from ``starting_Y``, repeatedly decreases Y until pullback
    points can be computed without arithmetic errors.

    INPUT:

    - ``space`` -- KleinianMaassFormSpace; the space
    - ``M`` -- integer; truncation parameter
    - ``Q`` -- integer; number of sample points parameter
    - ``starting_Y`` -- real (default: 0.73); initial Y value to try
    - ``max_iterations`` -- integer (default: 100); maximum number of attempts

    OUTPUT:

    - Real; the largest valid Y found

    EXAMPLES::

        sage: from maass_forms_klein.modform.compute_coefficients import find_max_y
        sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
        sage: space = KleinianMaassFormSpace(-4)
        sage: Y = find_max_y(space, M=1, Q=2)  # doctest: +SKIP
    """
    if not starting_Y:
        starting_Y = 0.73
    Y = starting_Y
    for _i in range(max_iterations):
        try:
            get_pb_pts(space.group(), Q, Y)
        except ArithmeticError:
            log.debug("Arithmetic error, trying smaller Y")
            Y = Y * 0.98
        else:
            log.debug(f"Y={Y} is ok")
            break
    return Y


def compute_coefficients(
    space: "KleinianMaassFormSpace",
    spectral_parameter: Complex_t | Real_t,
    Q: tuple[Integer_t] | Integer_t | None = None,
    Y: Real_t | None = None,
    M: tuple[Integer_t] | Integer_t | None = None,
    set_coefficients: dict | None = None,
    return_mat: bool = False,
    use_numpy: bool = True,
) -> KleinianMaassFormCoefficients:
    r"""
    Compute Fourier coefficients of a Kleinian Maass form.

    Sets up and solves a linear system arising from the Fourier expansion
    evaluated at pullback points at two different heights.

    INPUT:

    - ``space`` -- KleinianMaassFormSpace; the ambient space
    - ``spectral_parameter`` -- complex or real; the eigenvalue parameter
    - ``Q`` -- integer or None (default: None); number of sample points.
      If None, determined from M
    - ``Y`` -- real or None (default: None); height parameter.
      If None, determined automatically
    - ``M`` -- integer or None (default: None); truncation parameter.
      If None, determined from spectral parameter
    - ``set_coefficients`` -- dict or None (default: None); coefficients
      to fix for normalization
    - ``return_mat`` -- bool (default: False); if True, return the matrix
      and right-hand side instead of solving
    - ``use_numpy`` -- bool (default: True); if True, use numpy for solving

    OUTPUT:

    - KleinianMaassFormCoefficients; the computed coefficients

    EXAMPLES::

        sage: from maass_forms_klein.modform.compute_coefficients import compute_coefficients
        sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
        sage: space = KleinianMaassFormSpace(-4)
        sage: s = CC(0.5, 6.62211934)
        sage: C = compute_coefficients(space, s, M=1, Q=3, Y=0.5)  # doctest: +SKIP
    """
    if not M:
        M = ceil((abs(spectral_parameter) + 12) / (6.28318530717959))
    zpb, zm, Q, M, Y = get_pb_pts_set_params(space, spectral_parameter, M=M, Y=Y, Q_set=Q)
    log.debug(f"Q={Q}, M={M}, Y={Y}")
    V = setup_matrix(space, spectral_parameter, Y, M, tuple(zpb), tuple(zm))
    # normalise:
    dual_lattice_values = space.group().dual_translation_lattice_vectors(M, return_indices=True)
    dual_lattice_values, dual_lattice_indices = dual_lattice_values
    delete_indices = []
    normalisation = {}
    if set_coefficients is None and space.is_cuspidal():
        if vector((1, 0)) in dual_lattice_values:
            set_coefficients = {(0, 0): 0, (1, 0): 1}
        elif vector((0, 1)) in dual_lattice_values:
            set_coefficients = {(0, 0): 0, (0, 1): 1}
        else:
            set_coefficients = {(0, 0): 0, tuple(dual_lattice_values[1]): 1}
    if isinstance(set_coefficients, str):
        set_coefficients = dict_from_json(set_coefficients)
    for c, value in set_coefficients.items():
        log.debug(f"c={c} type={type(c)} value={value}")
        if isinstance(c, Integer_t):
            ci = dual_lattice_values[c]
        else:
            ci = vector(c)
        if ci not in dual_lattice_values:
            raise ValueError(f"set_coefficients[{c}]={value} not in dual_lattice_values")
        ni = dual_lattice_values.index(ci)
        delete_indices.append(ni)
        normalisation[ni] = value
    V = V.delete_rows(delete_indices)
    RHS = 0
    for ni, value in normalisation.items():
        RHS += V.column(ni) * value
    V = V.delete_columns(delete_indices)
    if return_mat:
        return V, -RHS
    # Numpy is quicker solve
    if use_numpy:
        A = numpy.complex128(V.numpy())
        B = numpy.complex128(-RHS)
        X = numpy.linalg.solve(A, B)
    else:
        X = V.solve_right(-RHS)
    C = []
    skip_step = 0
    for n, _v in enumerate(dual_lattice_values):
        if n in normalisation:
            C.append(normalisation[n])
            skip_step += 1
        else:
            C.append(CC(X[n - skip_step]))
    coefficients = KleinianMaassFormCoefficients(
        C,
        M,
        Q,
        spectral_parameter,
        space,
        Y,
        coordinate_indices=dual_lattice_indices,
        coordinate_values=dual_lattice_values,
        set_coefficients=set_coefficients,
    )
    return coefficients


@cached_function
def get_pb_pts(
    group: KleinianGroup, Q: Integer_t = None, Y: Real_t = None, prec: int = 53, check: bool = True
) -> tuple:
    """
    Get the list of points in the scaled lattice together with the corresponding pullbacks.

    EXAMPLES:
        sage: from maass_forms_klein.modform.compute_coefficients import get_pb_pts
        sage: from maass_forms_klein.all import KleinianMaassFormSpace
        sage: space = KleinianMaassFormSpace(-4)
        sage: get_pb_pts(space.group(), 1, 0.5)
        ([0.000000000000000 + 0.000000000000000i + 2.00000000000000j,
          0.000000000000000 + 0.000000000000000i + 1.00000000000000j,
          0.000000000000000 + 0.000000000000000i + 1.00000000000000j,
          -0.500000000000000 - 0.500000000000000i + 1.00000000000000j],
         [0.000000000000000 + 0.000000000000000i + 0.500000000000000j,
          0.000000000000000 + 0.500000000000000i + 0.500000000000000j,
          0.500000000000000 + 0.000000000000000i + 0.500000000000000j,
          0.500000000000000 + 0.500000000000000i + 0.500000000000000j])
        sage: get_pb_pts(space.group(), Q=1, Y=1)
        Traceback (most recent call last):
        ...
        ArithmeticError: Point 0.000000000000000 + 0.000000000000000i + 1.00000000000000j has...
    """

    lattice_basis = group.translation_lattice().basis_matrix() / (2 * Q)
    lattice_basis = tuple(lattice_basis)
    zmpb = []
    zm = []
    # Ensure that elements have same parent
    s = Sequence(
        [lattice_basis[0][0], lattice_basis[0][1], lattice_basis[1][0], lattice_basis[1][1], Y]
    )
    parent = s.universe()
    lattice_basis = (
        (parent(lattice_basis[0][0]), parent(lattice_basis[0][1])),
        (parent(lattice_basis[1][0]), parent(lattice_basis[1][1])),
    )
    lattice_values = get_lattice_values(lattice_basis, 1 - Q, Q + 1, shift=True)
    for _n, v in enumerate(lattice_values):
        zm_elt = UpperHalfSpaceElement((v[0], v[1], Y))
        zm.append(zm_elt)
        try:
            pbpt, g = group.pullback(zm_elt)
        except ValueError as e:
            log.debug("zm=%s", zm_elt)
            raise e
        if check and pbpt.y() <= Y:
            raise ArithmeticError(
                f"Point {pbpt} has imaginary part not greater than {Y}. zm={zm_elt} g={g}"
            )
        zmpb.append(pbpt)
    return zmpb, zm


def get_pb_pts_set_params(
    space: "KleinianMaassFormSpace",
    spectral_parameter: ComplexNumber = None,
    M: Integer_t = None,
    Y: Real_t = None,
    Q_set: Integer_t = None,
    smax: float | Real_t = None,
    prec: Integer_t = None,
) -> tuple:
    r"""
    Determine pullback points and computation parameters.

    Computes suitable values of Q, M, and Y if not provided, then
    obtains pullback points for the coefficient computation.

    INPUT:

    - ``space`` -- KleinianMaassFormSpace; the ambient space
    - ``spectral_parameter`` -- complex or None; the eigenvalue parameter
    - ``M`` -- integer or None; truncation parameter
    - ``Y`` -- real or None; height parameter
    - ``Q_set`` -- integer or None; override for Q
    - ``smax`` -- real or None; maximum spectral parameter magnitude
    - ``prec`` -- integer or None; working precision in bits

    OUTPUT:

    - tuple ``(zpb, zm, Q, M, Y)`` where ``zpb`` are pullback points,
      ``zm`` are lattice points, and Q, M, Y are the determined parameters

    EXAMPLES::

        sage: from maass_forms_klein.modform.compute_coefficients import get_pb_pts_set_params
        sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
        sage: space = KleinianMaassFormSpace(-4)
        sage: s = CC(0.5, 6.62211934)
        sage: zpb, zm, Q, M, Y = get_pb_pts_set_params(space, s, M=1, Y=0.5)  # doctest: +SKIP
    """
    if hasattr(spectral_parameter, "parent"):
        prec = spectral_parameter.parent().prec()
    else:
        prec = 53
    if not spectral_parameter and not smax:
        raise ValueError("Need either spectral parameter or smax set.")
    if spectral_parameter and not smax:
        smax = abs(spectral_parameter)
    if not M:
        M = ceil((smax + 12) / (6.28318530717959) + 1)
    Q = M + 1
    if Q_set and Q_set >= Q:
        log.debug(f"Q_set {Q_set} is used instead of Qs={Q}")
        Q = Q_set
    # Try to find best value of Y
    Y = find_max_y(space, M, Q=Q, starting_Y=Y)
    Y = RealField(prec)(Y)
    # We try with given Y and if it doesn't work we keep decreasing Y until it does.
    try:
        zpb, zm = get_pb_pts(space.group(), Q, Y, prec=prec)
    except ArithmeticError as e:
        msg = f"Could not find good pullback points. Error: {e}"
        log.debug(msg)
        raise ArithmeticError(msg) from e
    return zpb, zm, Q, M, Y
