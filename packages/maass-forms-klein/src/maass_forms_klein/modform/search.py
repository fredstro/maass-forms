import logging
from typing import Optional, ParamSpec

import numpy
from sage.all import RR
from sage.functions.other import ceil, imag, real
from sage.misc.cachefunc import cached_function
from werkzeug.datastructures import ImmutableDict

from maass_forms_klein.modform.compute_coefficients import compute_coefficients
from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
from maass_forms_klein.modform.utils import Integer_t, Real_t
from maass_form_core.utils.json_converters import dict_from_json

P = ParamSpec("P")

log = logging.getLogger(__name__)


@cached_function
def h(
    r,
    space: KleinianMaassFormSpace,
    M: Integer_t,
    Q: Integer_t,
    Y1: Real_t,
    Y2: Real_t,
    set_coefficients: Optional[dict | str] = None,
    coefficient_indices: Optional[tuple[Integer_t, ...]] = None,
    real_imag: str = "real",
) -> tuple[Real_t, ...] | tuple[tuple[Real_t, Real_t], ...]:
    r"""
    Compute the difference of coefficients at two heights for eigenvalue search.

    This function computes coefficients at heights Y1 and Y2 and returns
    the difference. A true eigenvalue will give matching coefficients at
    both heights, so zero of this function indicates an eigenvalue.

    INPUT:

    - ``r`` -- real; the spectral parameter to test
    - ``space`` -- KleinianMaassFormSpace; the ambient space
    - ``M`` -- integer; truncation parameter
    - ``Q`` -- integer; number of sample points
    - ``Y1`` -- real; first height
    - ``Y2`` -- real; second height
    - ``set_coefficients`` -- dict or str or None (default: None);
      coefficients to fix for normalization
    - ``coefficient_indices`` -- tuple of integers or None (default: None);
      which coefficient indices to compare
    - ``real_imag`` -- str (default: ``'real'``); which part to return:
      ``'real'``, ``'imag'``, or ``'both'``

    OUTPUT:

    - tuple of reals or tuple of pairs of reals; the coefficient differences

    EXAMPLES::

        sage: from maass_forms_klein.modform.search import h
        sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
        sage: space = KleinianMaassFormSpace(-4)
        sage: h(RR(6.5), space, 3, 4, 0.5, 0.475)  # doctest: +SKIP
    """
    if isinstance(set_coefficients, str):
        set_coefficients = dict_from_json(set_coefficients)
    C0 = compute_coefficients(space, r, Q, Y1, M, set_coefficients)
    C1 = compute_coefficients(space, r, Q, Y2, M, set_coefficients)
    if not coefficient_indices:
        coefficient_indices = [1]
    for i in coefficient_indices:
        log.debug(
            f"Computing h({r})({i}) = {C0[i]} - {C1[i]}="
            f"{(C0[i] - C1[i]).real()}, {(C0[i] - C1[i]).imag()}"
        )
    if real_imag == "imag":
        return tuple([imag(C0[i] - C1[i]) for i in coefficient_indices])
    if real_imag == "real":
        return tuple([real(C0[i] - C1[i]) for i in coefficient_indices])
    return tuple([((C0[i] - C1[i]).real(), (C0[i] - C1[i]).imag()) for i in coefficient_indices])


def search_eigenvalues(
    space: KleinianMaassFormSpace, R1: Real_t, R2: Real_t, **kwargs: P.kwargs
) -> list:
    r"""
    Search for eigenvalues of Maass forms for the space.

    EXAMPLES:

        sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
        sage: from maass_forms_klein.modform.search import search_eigenvalues
        sage: space = KleinianMaassFormSpace(-4)
        sage: R1 = 6
        sage: R2 = 7
        sage: search_eigenvalues(space, R1, R2, M=3, Y=0.6) # long time
        sage: from snappy import Manifold
        sage: M = Manifold('4_1')
        sage: space = KleinianMaassFormSpace(M)
        sage: R1 = 0
        sage: R2 = 1
        sage: search_eigenvalues(space, R1, R2, M=2, Y=0.15)
        Traceback (most recent call last):
        ...
        ValueError:...
        sage: search_eigenvalues(space, R1, R2, M=2, Y=0.15,num_steps=5) # long time
        Traceback (most recent call last):
        ...
        ArithmeticError:...

    """
    if "step_size" in kwargs:
        step_size = kwargs.get("step_size", 0.05)
        num_steps = ceil((R2 - R1) / step_size + 1)
    elif "num_steps" not in kwargs:
        raise ValueError("Need either step size or number of steps.")
    else:
        num_steps = int(kwargs.get("num_steps", 11))
    range_of_r = numpy.linspace(R1, R2, num_steps)
    M = kwargs.get("M", 4)
    Q = kwargs.get("Q", M + 5)
    Y1 = kwargs.get("Y", 0.5)
    Y2 = kwargs.get("Y2", Y1 * 0.95)
    set_coefficients = kwargs.get("set_coefficients", None)
    coefficient_indices = kwargs.get("coefficient_indices", None)
    if set_coefficients:
        set_coefficients = ImmutableDict(set_coefficients)

    r0 = R1
    to_check = []
    for r in range_of_r[1:]:
        # Check for sign change in interval
        h0 = h(
            RR(r0),
            space,
            M,
            Q,
            Y1,
            Y2,
            set_coefficients=set_coefficients,
            coefficient_indices=coefficient_indices,
        )
        h1 = h(
            RR(r),
            space,
            M,
            Q,
            Y1,
            Y2,
            set_coefficients=set_coefficients,
            coefficient_indices=coefficient_indices,
        )
        tests = [h0[i] * h1[i] <= 0 or h0[i] * h1[i] <= 0 for i in range(len(h0))]
        if any(tests):
            to_check.append((r0, r, tests.count(True)))
        r0 = r
    return to_check


def newton_method_search(
    r0: Real_t,
    r1: Real_t,
    space: KleinianMaassFormSpace,
    M: Integer_t,
    Q: Integer_t,
    Y1: Real_t,
    Y2: Real_t,
    set_coefficients: dict,
    coefficient_indices: Optional[tuple[Integer_t, ...]] = None,
    tolerance: Real_t = 1e-10,
    verbose: bool = False,
    real_imag: str = "real",
) -> tuple[Real_t, tuple]:
    r"""
    Search for eigenvalues of Maass forms for the space.

    """
    log.debug(f"step: {r0} {r1}:")
    r = newton_method_step(
        r0,
        r1,
        space,
        M,
        Q,
        Y1,
        Y2,
        set_coefficients,
        coefficient_indices=coefficient_indices,
        tolerance=tolerance,
        real_imag=real_imag,
    )

    log.debug(f"first rnew: {r0} {r1} = {r}")
    if r < r0 or r > r1:
        msg = f"r_new={r} outside of the interval [{r0}, {r1}]"
        log.debug(msg)
        raise ValueError(msg)
    f = h(
        r,
        space,
        M,
        Q,
        Y1,
        Y2,
        set_coefficients,
        coefficient_indices=coefficient_indices,
        real_imag=real_imag,
    )
    log.debug(f"h({r}) = {f}")

    err = sum([abs(x) for x in f])
    if err < tolerance:
        return r, f
    while err > tolerance:
        try:
            r_new = newton_method_step(
                r0,
                r,
                space,
                M,
                Q,
                Y1,
                Y2,
                set_coefficients,
                coefficient_indices=coefficient_indices,
                tolerance=tolerance,
                real_imag=real_imag,
            )
            log.debug(f"1: r_new={r_new} r0={r0} r={r} r1={r1}")
            if r_new < r0 or r_new > r:
                msg = f"r_new={r_new} outside of the interval [{r0}, {r}]"
                log.debug(msg)
                raise ValueError(msg)
            r0 = r0
            r1 = r
        except ValueError:
            r_new = newton_method_step(
                r,
                r1,
                space,
                M,
                Q,
                Y1,
                Y2,
                set_coefficients,
                coefficient_indices=coefficient_indices,
                tolerance=tolerance,
                real_imag=real_imag,
            )
            log.debug(f"2: r_new={r_new} r0={r0} r={r} r1={r1}")
            if r_new < r or r_new > r1:
                msg = f"r_new={r_new} outside of the interval [{r}, {r1}]"
                log.debug(msg)
                raise ValueError(msg) from err
            r0 = r
            r1 = r1
        log.debug("r-diff: " + str(abs(r - r_new)))
        if abs(r - r_new) < tolerance:
            break
        f = h(
            r_new,
            space,
            M,
            Q,
            Y1,
            Y2,
            set_coefficients,
            coefficient_indices=coefficient_indices,
            real_imag=real_imag,
        )
        r = r_new
        err = sum([abs(x) for x in f])
        log.debug(f"f(r_new)={f} err={err}")
    return r, f


def newton_method_step(
    r0,
    r1,
    space,
    M,
    Q,
    Y1,
    Y2,
    set_coefficients: Optional[dict] = None,
    coefficient_indices: Optional[tuple[Integer_t, ...]] = None,
    tolerance: Real_t = 1e-10,
    real_imag: str = "real",
) -> Real_t:
    r"""
    Perform a single Newton method step for eigenvalue search.

    Computes the next approximation to a zero of :func:`h` using
    the secant method on the interval ``[r0, r1]``.

    INPUT:
    - ``r0`` -- initial guess
    - ``r1`` -- final guess
    - ``space`` -- space
    - ``M`` -- M
    - ``Q`` -- Q
    - ``Y1`` -- Y1
    - ``Y2`` -- Y2
    - ``set_coefficients`` -- coefficients to be set
    - ``coefficient_indices`` -- indices of the coefficients to be set
    - ``real_imag`` -- version of the function to be used (default: ``real``)
         - ``real`` use real part of C1 - C2
         - ``imag`` use imaginary part of C1 - C2
         - ``both`` use both real and imaginary parts
    """
    h0 = h(r0, space, M, Q, Y1, Y2, set_coefficients, coefficient_indices=coefficient_indices)
    log.debug(f"h(r0)({r0}) = {h0}")
    h1 = h(r1, space, M, Q, Y1, Y2, set_coefficients, coefficient_indices=coefficient_indices)
    log.debug(f"h(r1)({r1}) = {h1}")
    if any(h0[i] * h1[i] > 0 and min(h0[i], h1[i]) > tolerance for i in range(len(h0))):
        raise ValueError("No zero in this interval")
    if len(h0) == 1:  # real_imag in ['real', 'imag']:
        return r0 - h0[0] * (r0 - r1) / (h0[0] - h1[0])
    new_r = sum([r0 - h0[i] * (r0 - r1) / (h0[i] - h1[i]) for i in range(len(h0))])
    return new_r / len(h0)


def brute_force_search(
    r0,
    r1,
    numpts,
    space,
    M,
    Q,
    Y1,
    Y2,
    set_coefficients: Optional[dict] = None,
    coefficient_indices: Optional[tuple[Integer_t, ...]] = None,
    tolerance: Real_t = 1e-10,
    verbose: bool = False,
):
    pts = numpy.linspace(r0, r1, numpts)
    res = []
    for n, r00 in enumerate(pts):
        r01 = pts[n + 1] if n < len(pts) - 1 else r1
        try:
            search = newton_method_search(
                r00,
                r01,
                space,
                M,
                Q,
                Y1,
                Y2,
                set_coefficients,
                coefficient_indices=coefficient_indices,
                tolerance=tolerance,
                verbose=verbose,
            )
            log.debug(f"search={r00}, {r01}, {search}")
            res.append(search)
        except ValueError as e:
            log.debug(e)
            continue
    return res
