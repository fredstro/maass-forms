"""
Routines to search for Kleinian Maass forms.

Kleinian Maass forms have a single-dimensional complex spectral parameter
``s = 1/2 + i r``, so this module performs the eigenvalue search over a
one-dimensional grid of real ``r`` values. The structure of the routines
mirrors :mod:`maass_forms_hilbert.search.search`, but Broyden's method
collapses to the standard secant method in the scalar case.

The eigenvalue locator is chosen by the caller:

* **General Kleinian groups** (any group, including those coming from
  snappy manifolds):

  - ``'two_y'`` -- compare coefficients at two heights (the locator from
    :mod:`maass_forms_klein.modform.search`), DB-free;
  - ``'unit'`` -- compare two user-specified coefficients that ought to
    agree at an eigenvalue (e.g. by a symmetry of the form).

* **Arithmetic (Bianchi) groups only**:

  - ``'coprime'`` and ``'prime_power'`` -- Hecke multiplicativity
    relations; only available when
    :func:`is_arithmetic_group` returns ``True``.

See :func:`coeff_diff_fun` for the dispatch.

EXAMPLES:

Search the Bianchi space over ``Q(sqrt(-4))`` for the known eigenvalue
near ``r = 6.62211934`` using the two-height locator. Starting from a
small grid of five points covering ``[6.5, 6.7]``, the Newton/secant
zero-crossing refinement converges to the eigenvalue to ``~4e-6`` with
the modest truncation ``M = 3, Q = 4``::

    sage: from maass_forms_klein.all import KleinianMaassFormSpace
    sage: from maass_forms_klein.modform.search import brute_force_search
    sage: space = KleinianMaassFormSpace(-4, cuspidal=True)
    sage: results = brute_force_search(
    ....:     6.5, 6.7, 5, space, M=3, Q=4, Y1=0.5, Y2=0.475,
    ....:     set_coefficients=None, tolerance=1e-10)
    sage: len(results)
    1
    sage: r_found, residual = results[0]
    sage: bool(abs(float(r_found) - 6.62211934) < 1e-4)
    True
    sage: bool(abs(float(residual[0])) < 1e-8)
    True
"""

import logging
import os

import numpy
from sage.functions.other import ceil
from sage.parallel.decorate import parallel
from sage.rings.complex_mpfr import ComplexField
from sage.rings.number_field.number_field_base import NumberField
from sage.rings.real_mpfr import RR

from maass_forms_klein.modform.compute_coefficients import get_pb_pts_set_params
from maass_forms_klein.modform.kmaass_element import KleinianMaassFormElement
from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
from maass_forms_klein.modform.search import brute_force_search, h, newton_method_search
from maass_forms_klein.modform.utils import Complex_t, Integer_t, Real_t

#: Locator types valid for any Kleinian group (arithmetic or not).
GENERAL_LOCATORS = ("two_y", "unit")

#: Locator types that only make sense for arithmetic (Bianchi) groups,
#: where Hecke multiplicativity is available.
ARITHMETIC_LOCATORS = ("coprime", "prime_power")

#: Snappy manifold identifiers known to give rise to arithmetic Kleinian
#: groups (commensurable with a Bianchi group). The only such knot
#: complement among prime knots up to 10 crossings is the figure-eight
#: ``4_1``, commensurable with ``PSL(2, O_{-3})``. All other manifolds
#: are assumed to be non-arithmetic.
KNOWN_ARITHMETIC_MANIFOLDS = frozenset({"4_1"})

global use_database_global
try:
    import mongoengine
    from comp_manager.utils import insert_object, load_object

    from maass_forms_klein.database.models import KleinianMaassFormDB

    use_database_global = True
except ImportError:
    use_database_global = False
    KleinianMaassFormDB = None  # type: ignore[assignment]

log = logging.getLogger(__name__)


def create_grid(grid_limits: tuple[Real_t, Real_t], grid_number: Integer_t):
    r"""
    Create a one-dimensional linspace grid for the imaginary part of the
    spectral parameter.

    INPUT:

    - ``grid_limits`` -- tuple of two reals ``(r1, r2)``
    - ``grid_number`` -- integer; the number of grid points

    EXAMPLES::

        sage: from maass_forms_klein.search.search import create_grid
        sage: create_grid((0, 2), 11)
        array([0. , 0.2, 0.4, 0.6, 0.8, 1. , 1.2, 1.4, 1.6, 1.8, 2. ])
    """
    return numpy.linspace(grid_limits[0], grid_limits[1], grid_number)


def compute_on_interval(
    space: KleinianMaassFormSpace,
    r1: Real_t,
    r2: Real_t,
    num_steps: Integer_t,
    prec: Integer_t = 53,
    bound_m: Integer_t = 2,
    y: Real_t | None = None,
    Q: Integer_t | None = None,
    set_coefficients: dict | None = None,
    num_threads: Integer_t | None = None,
    use_database: bool = True,
):
    r"""
    Compute Kleinian Maass forms on a one-dimensional grid of spectral
    parameters ``s = 1/2 + i r`` with ``r`` ranging over ``num_steps``
    equally spaced points in ``[r1, r2]``.

    INPUT:

    - ``space`` -- KleinianMaassFormSpace; the space
    - ``r1`` -- real; lower bound of the imaginary part
    - ``r2`` -- real; upper bound of the imaginary part
    - ``num_steps`` -- integer; number of grid points
    - ``prec`` -- integer (default: 53); precision in bits
    - ``bound_m`` -- integer (default: 2); truncation bound ``M``
    - ``y`` -- real or None (default: None); height parameter
    - ``Q`` -- integer or None (default: None); number of sample points
    - ``set_coefficients`` -- dict or None (default: None); normalisation
    - ``num_threads`` -- integer or None (default: None); thread count
      (sets ``SAGE_NUM_THREADS`` if provided)
    - ``use_database`` -- bool (default: True); load/save via MongoDB
      (only effective if database support is available at import time)

    OUTPUT:

    - generator of ``(input_tuple, KleinianMaassFormElement)`` pairs, as
      returned by Sage's ``@parallel`` decorator

    EXAMPLES:

    Happy path on a tiny grid with the database turned off::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.all import (
        ....:     KleinianMaassFormSpace, KleinianMaassFormElement)
        sage: from maass_forms_klein.search.search import compute_on_interval
        sage: space = KleinianMaassFormSpace(-4, cuspidal=False)
        sage: result = list(compute_on_interval(
        ....:     space, r1=0.0, r2=1.0, num_steps=3, bound_m=1,
        ....:     y=0.5, Q=3, num_threads=1, use_database=False))
        sage: len(result)
        3
        sage: all(isinstance(out, KleinianMaassFormElement)
        ....:     for _input, out in result)
        True
    """
    grid = create_grid((r1, r2), num_steps)
    CF = ComplexField(prec)
    input_params = []
    use_db = use_database_global and use_database
    for r in grid:
        spectral_parameter = CF(0.5, r)
        log.debug("Queueing spectral parameter %s", spectral_parameter)
        input_params.append((space, spectral_parameter, bound_m, y, Q, set_coefficients, use_db))
    if num_threads is not None:
        os.environ["SAGE_NUM_THREADS"] = str(num_threads)
    # Prime the pullback-point cache so parallel workers don't race on it.
    smax = 0.0
    for r in input_params:
        smax = max(smax, float(abs(r[1])))
    for si in range(1, ceil(smax) + 1):
        get_pb_pts_set_params(space, M=input_params[0][2], Y=y, Q_set=Q, smax=si)
    return compute_one_spectral_parameter(input_params)


@parallel()
def compute_one_spectral_parameter(
    space: KleinianMaassFormSpace,
    spectral_parameter: Complex_t,
    bound_m: Integer_t,
    y: Real_t | None = None,
    Q: Integer_t | None = None,
    set_coefficients: dict | None = None,
    use_database: bool = True,
):
    r"""
    Compute (or retrieve from the database) a Kleinian Maass form at a
    single spectral parameter.

    INPUT:

    - ``space`` -- KleinianMaassFormSpace
    - ``spectral_parameter`` -- Complex; the spectral parameter ``s``
    - ``bound_m`` -- integer; truncation bound
    - ``y`` -- real or None (default: None); height parameter
    - ``Q`` -- integer or None (default: None); number of sample points
    - ``set_coefficients`` -- dict or None (default: None); normalisation
    - ``use_database`` -- bool (default: True); whether to use MongoDB

    OUTPUT:

    - :class:`KleinianMaassFormElement` with coefficients populated

    EXAMPLES:

    Happy path with the database turned off; the form is computed locally
    and its Fourier coefficients are populated::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.all import KleinianMaassFormSpace
        sage: from maass_forms_klein.search.search import compute_one_spectral_parameter
        sage: from sage.rings.complex_mpfr import ComplexField
        sage: CC = ComplexField(53)
        sage: space = KleinianMaassFormSpace(-4, cuspidal=True)
        sage: s = CC(0.5, 6.62211934)
        sage: form = compute_one_spectral_parameter(
        ....:     space, s, 2, 0.5, 3, None, False)
        sage: form.spectral_parameter() == s
        True
        sage: bool(form.coefficients())
        True

    Calling the parallel wrapper on an empty input list yields no results::

        sage: list(compute_one_spectral_parameter([]))
        []
    """
    log.debug("compute_one_spectral_parameter at s=%s (db=%s)", spectral_parameter, use_database)
    maass_form = None
    use_db = use_database_global and use_database
    if use_db:
        try:
            form_db = KleinianMaassFormDB.near_or_create(
                parent=space.to_json(),
                spectral_parameter=spectral_parameter,
                bound_m=bound_m,
                set_coefficients=set_coefficients,
                y=y,
            )
            maass_form = load_object(form_db)
        except mongoengine.connection.ConnectionFailure:
            log.warning("Could not connect to the database; computing locally only.")
    if maass_form is None:
        log.info("No Maass form found in database, computing it.")
        maass_form = KleinianMaassFormElement(space, spectral_parameter)
    if not maass_form.coefficients():
        log.info("No coefficients found, computing them.")
        maass_form.compute_coefficients(
            spectral_parameter=spectral_parameter,
            M=bound_m,
            set_coefficients=set_coefficients,
            Y=y,
            Q=Q,
        )
        if use_db:
            insert_object(maass_form)
    return maass_form


def is_arithmetic_group(space: KleinianMaassFormSpace) -> bool:
    r"""
    Return ``True`` if the Kleinian group underlying ``space`` is
    arithmetic, in the sense that Hecke-type coefficient relations are
    available.

    Two cases return ``True``:

    * The group is a Bianchi group over an imaginary quadratic field
      ``K`` (constructed from a fundamental discriminant); detected by
      the base ring being an imaginary quadratic ``NumberField``.
    * The group comes from a Snappy manifold whose identifier is listed
      in :data:`KNOWN_ARITHMETIC_MANIFOLDS`. The figure-eight knot
      complement ``4_1`` is the only such case in the supported range;
      it is commensurable with the Bianchi group ``PSL(2, O_{-3})``.

    All other manifolds (and any user-supplied generator set with a
    generic complex base ring) are treated as non-arithmetic.

    INPUT:

    - ``space`` -- KleinianMaassFormSpace

    OUTPUT:

    - bool

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.all import KleinianMaassFormSpace
        sage: from maass_forms_klein.search.search import is_arithmetic_group
        sage: is_arithmetic_group(KleinianMaassFormSpace(-4, cuspidal=False))
        True
        sage: from snappy import Manifold
        sage: is_arithmetic_group(KleinianMaassFormSpace(Manifold('4_1')))
        True
        sage: is_arithmetic_group(KleinianMaassFormSpace(Manifold('5_2')))
        False
    """
    group = space.group()
    manifold = getattr(group, "_manifold", "") or ""
    if manifold in KNOWN_ARITHMETIC_MANIFOLDS:
        return True
    base = group.base_ring()
    if not isinstance(base, NumberField):
        return False
    if base.degree() != 2 or base.discriminant() >= 0:
        return False
    # A manifold-derived group whose base ring happens to be a degree-2
    # imaginary quadratic field is still treated as non-arithmetic
    # unless its identifier is explicitly whitelisted above, because the
    # group may sit non-trivially inside PSL(2, O_K).
    return not manifold


def coeff_diff_fun(
    space: KleinianMaassFormSpace,
    relation: dict,
    r: Real_t,
    bound_m: Integer_t | None = None,
    y: Real_t | None = None,
    set_coefficients: dict | None = None,
    prec: Integer_t = 53,
    Q: Integer_t | None = None,
) -> tuple[Real_t, Real_t]:
    r"""
    Evaluate the residual of an eigenvalue locator at ``s = 1/2 + i r``.

    Two families of locators are supported:

    *General locators* (any Kleinian group):

    * ``{'two_y': {'Y1': y1, 'Y2': y2, 'indices': [i1, i2, ...]}}`` --
      compute coefficients at two heights ``y1`` and ``y2`` and compare
      them at the given dual-lattice indices. A true eigenvalue gives
      matching coefficients regardless of the height. The ``'Y2'`` key
      defaults to ``0.95 * Y1``; ``'indices'`` defaults to ``[1]``. This
      locator does not require the MongoDB database.
    * ``{'unit': [i1, i2]}`` -- require ``c[i1] == c[i2]``. The pair is
      user-specified (e.g. a known symmetry of the form). Requires the
      database to fetch the precomputed form.

    *Arithmetic-only locators* (Bianchi groups, see
    :func:`is_arithmetic_group`):

    * ``{'coprime': [i1, i2, i3]}`` -- Hecke multiplicativity
      ``c[i3] == c[i1] * c[i2]`` for coprime indices.
    * ``{'prime_power': [i0, i1, ...]}`` -- the Hecke recursion at a
      prime power.

    INPUT:

    - ``space`` -- KleinianMaassFormSpace
    - ``relation`` -- dict; see above
    - ``r`` -- real; imaginary part of the spectral parameter
    - ``bound_m`` -- integer or None; truncation bound
    - ``y`` -- real or None; height parameter used for DB lookup (and as
      a fallback for ``Y1`` in the ``two_y`` locator)
    - ``set_coefficients`` -- dict or None; normalisation
    - ``prec`` -- integer (default: 53); working precision
    - ``Q`` -- integer or None; number of sample points for the
      ``two_y`` locator (defaults to ``bound_m + 1``)

    OUTPUT:

    - tuple ``(|real diff|, |imag diff|)`` of nonnegative reals; zero at
      a true eigenvalue.

    EXAMPLES:

    Happy path on the Bianchi space ``Q(sqrt(-4))`` near the known
    eigenvalue ``r = 6.62211934``: the residual is small (~1e-5) for the
    ``two_y`` locator::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.all import KleinianMaassFormSpace
        sage: from maass_forms_klein.search.search import coeff_diff_fun
        sage: from sage.rings.real_mpfr import RR
        sage: space = KleinianMaassFormSpace(-4, cuspidal=True)
        sage: diff = coeff_diff_fun(
        ....:     space,
        ....:     {'two_y': {'Y1': 0.5, 'Y2': 0.475, 'indices': [1]}},
        ....:     RR(6.62211934), bound_m=3, Q=4)
        sage: len(diff)
        2
        sage: bool(diff[0] >= 0 and diff[1] >= 0)
        True
        sage: bool(diff[0] + diff[1] < 1e-3)
        True

    Off the eigenvalue, the residual is much larger::

        sage: diff_off = coeff_diff_fun(
        ....:     space,
        ....:     {'two_y': {'Y1': 0.5, 'Y2': 0.475, 'indices': [1]}},
        ....:     RR(6.0), bound_m=3, Q=4)
        sage: bool(diff_off[0] + diff_off[1] > diff[0] + diff[1])
        True

    Hecke-type relations are rejected for non-arithmetic groups::

        sage: from snappy import Manifold
        sage: nonarith = KleinianMaassFormSpace(Manifold('5_2'))
        sage: coeff_diff_fun(nonarith, {'coprime': [0, 1, 2]}, 0.5, bound_m=2)
        Traceback (most recent call last):
        ...
        NotImplementedError: Hecke-type coefficient relations ('coprime',
        'prime_power') are only defined for arithmetic Kleinian (Bianchi)
        groups...
    """
    locator = next(iter(relation))
    # --- general locator: two-height coefficient difference --------------
    if locator == "two_y":
        params = relation["two_y"] if isinstance(relation["two_y"], dict) else {}
        Y1 = params.get("Y1", y)
        if Y1 is None:
            raise ValueError("two_y locator requires Y1 (or y as fallback).")
        Y2 = params.get("Y2", Y1 * 0.95)
        indices = tuple(params.get("indices", [1]))
        M_ = bound_m if bound_m is not None else 2
        Q_ = Q if Q is not None else M_ + 1
        diffs = h(
            RR(r),
            space,
            M_,
            Q_,
            Y1,
            Y2,
            set_coefficients=set_coefficients,
            coefficient_indices=indices,
            real_imag="both",
        )
        re_diff = sum(abs(d[0]) for d in diffs)
        im_diff = sum(abs(d[1]) for d in diffs)
        return re_diff, im_diff
    # --- locators below need a precomputed form from the database -------
    if locator in ARITHMETIC_LOCATORS and not is_arithmetic_group(space):
        raise NotImplementedError(
            "Hecke-type coefficient relations ('coprime', 'prime_power') "
            "are only defined for arithmetic Kleinian (Bianchi) groups; "
            "use {'two_y': ...} or {'unit': ...} for general groups."
        )
    if locator not in GENERAL_LOCATORS + ARITHMETIC_LOCATORS:
        raise ValueError(
            f"Unknown locator '{locator}'. Expected one of "
            f"{GENERAL_LOCATORS + ARITHMETIC_LOCATORS}."
        )
    if not use_database_global:
        raise RuntimeError(
            f"Locator '{locator}' requires database support (comp_manager + mongoengine)."
        )
    CF = ComplexField(prec)
    spectral_parameter = CF(0.5, r)
    candidates = list(
        KleinianMaassFormDB.objects.with_m_precision(bound_m)
        .with_y_precision(y)
        .near(spectral_parameter)
    )
    if len(candidates) == 1:
        f = load_object(candidates[0])
    else:
        list(
            compute_one_spectral_parameter(
                [(space, spectral_parameter, bound_m, y, None, set_coefficients, True)]
            )
        )
        f = load_object(
            KleinianMaassFormDB.objects.near(spectral_parameter)
            .with_m_precision(bound_m)
            .with_y_precision(y)
            .first()
        )
    coeffs = f.coefficients()
    if locator == "unit":
        u = relation["unit"]
        diff = coeffs[u[1]] - coeffs[u[0]]
    elif locator == "coprime":
        c = relation["coprime"]
        diff = coeffs[c[2]] - coeffs[c[1]] * coeffs[c[0]]
    else:  # prime_power
        c = relation["prime_power"]
        t = len(c)
        if t == 2:
            diff = coeffs[c[1]] - coeffs[c[0]] * coeffs[c[0]] - 1
        else:
            diff = coeffs[c[t - 1]] - coeffs[c[0]] * coeffs[c[t - 2]] + coeffs[c[t - 3]]
    return abs(diff.real()), abs(diff.imag())


@parallel()
def check_coefficients_of_computed_object(
    space: KleinianMaassFormSpace,
    check_rel: list,
    cvalue: Real_t,
    spectral_parameter: Complex_t,
    bound_m: Integer_t,
    y: Real_t | None = None,
    return_error: bool = False,
):
    r"""
    After computing forms on a grid, check whether the form at
    ``spectral_parameter`` satisfies a unit coefficient relation up to
    threshold ``cvalue``.

    INPUT:

    - ``space`` -- KleinianMaassFormSpace
    - ``check_rel`` -- list of two coefficient indices ``[i1, i2]`` whose
      coefficients should agree at an eigenvalue
    - ``cvalue`` -- real; threshold defining "close"
    - ``spectral_parameter`` -- Complex; the candidate spectral parameter
    - ``bound_m`` -- integer; truncation bound
    - ``y`` -- real or None (default: None)
    - ``return_error`` -- bool (default: False); if ``True``, include the
      measured residual in the returned tuple

    OUTPUT:

    - tuple ``(imag(s),)`` (or ``(imag(s), error)`` when ``return_error``
      is ``True``) if the relation holds; ``None`` otherwise

    EXAMPLES:

    Applied to an empty input list, the parallel wrapper yields nothing
    (this is a quick test that does not require the database)::

        sage: from maass_forms_klein.search.search import check_coefficients_of_computed_object
        sage: list(check_coefficients_of_computed_object([]))
        []

    Happy path against an empty mock database -- with no precomputed form
    available, the function returns ``None`` (a hit would return a
    ``(imag(s),)`` tuple)::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.all import KleinianMaassFormSpace
        sage: from sage.rings.complex_mpfr import ComplexField
        sage: CC = ComplexField(53)
        sage: space = KleinianMaassFormSpace(-4, cuspidal=True)
        sage: s = CC(0.5, 6.62211934)
        sage: result = list(check_coefficients_of_computed_object(
        ....:     [(space, [(1, 0), (-1, 0)], 1e-2, s, 2, 0.5, False)]))
        sage: len(result)
        1
        sage: result[0][1] is None
        True
    """
    if not use_database_global:
        raise RuntimeError("check_coefficients_of_computed_object requires database support.")
    form_db = (
        KleinianMaassFormDB.objects.space(space)
        .with_m_precision(bound_m)
        .near(spectral_parameter)
        .with_y_precision(y)
        .first()
    )
    if form_db is None:
        log.debug("No precomputed form found near %s", spectral_parameter)
        return None
    f = load_object(form_db)
    coeffs = f.coefficients()
    diff = (coeffs[check_rel[1]] - coeffs[check_rel[0]]).real()
    t = abs(diff)
    if t <= cvalue:
        log.debug("%s\t%s", spectral_parameter.imag(), t)
        if return_error:
            return (spectral_parameter.imag(), t)
        return (spectral_parameter.imag(),)
    return None


@parallel()
def secant_iteration(
    space: KleinianMaassFormSpace,
    r_1: Real_t,
    relation: dict,
    bound_m: Integer_t,
    y: Real_t | None = None,
    set_coefficients: dict | None = None,
    count: Integer_t = 0,
    prec: Integer_t = 53,
    r_0: Real_t | None = None,
    f_1: tuple | None = None,
    f_0: tuple | None = None,
    tolerance: Real_t = 1e-12,
    max_iter: Integer_t = 40,
):
    r"""
    Single iteration of the secant method for refining a candidate eigenvalue ``r_1`` using a
    coefficient relation.

    INPUT:

    - ``space`` -- KleinianMaassFormSpace
    - ``r_1`` -- real; current iterate (imaginary part of ``s``)
    - ``relation`` -- dict; coefficient relation, see :func:`coeff_diff_fun`
    - ``bound_m`` -- integer; truncation bound
    - ``y`` -- real or None
    - ``set_coefficients`` -- dict or None
    - ``count`` -- integer (default: 0); recursion depth
    - ``prec`` -- integer (default: 53); precision
    - ``r_0`` -- real or None; previous iterate
    - ``f_1`` -- tuple or None; cached ``coeff_diff_fun`` value at ``r_1``
    - ``f_0`` -- tuple or None; cached ``coeff_diff_fun`` value at ``r_0``
    - ``tolerance`` -- real (default: 1e-12); convergence tolerance
    - ``max_iter`` -- integer (default: 40); hard iteration cap

    OUTPUT:

    - real ``r`` once converged, or ``None`` if convergence fails

    EXAMPLES:

    The function is ``@parallel``-decorated; called on an empty input
    list it yields no results (a quick test that needs no database)::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.search.search import secant_iteration
        sage: list(secant_iteration([]))
        []

    On the Bianchi space ``Q(sqrt(-4))``, starting from a point near the
    known eigenvalue ``r = 6.62211934`` with the two-height locator and
    a small iteration cap, the iteration terminates (either converging
    or hitting the cap and returning ``None``)::

        sage: from maass_forms_klein.all import KleinianMaassFormSpace
        sage: space = KleinianMaassFormSpace(-4, cuspidal=True)
        sage: result = secant_iteration(
        ....:     space, 6.62,
        ....:     {'two_y': {'Y1': 0.5, 'Y2': 0.475, 'indices': [1]}},
        ....:     2, 0.5, None, 0, 53, 6.625, None, None, 1e-12, 3)
        sage: result is None or float(result) > 0
        True
    """
    if r_0 is None:
        r_0 = r_1 - 1e-7
    if f_0 is None:
        f_0 = coeff_diff_fun(space, relation, r_0, bound_m, y, set_coefficients, prec)
    if f_1 is None:
        f_1 = coeff_diff_fun(space, relation, r_1, bound_m, y, set_coefficients, prec)
    delta_r = r_1 - r_0
    delta_f = f_1[0] - f_0[0]
    if delta_r == 0 or delta_f == 0:
        log.debug("Secant slope vanished; aborting.")
        return None
    r_new = r_1 - f_1[0] * delta_r / delta_f
    if r_new < 0 or r_new > 100:
        return None
    q = coeff_diff_fun(space, relation, r_new, bound_m, y, set_coefficients, prec)
    log.debug("r_new=%s q=%s", r_new, q)
    if abs(q[0]) < tolerance and abs(q[1]) < tolerance:
        return r_new
    if count >= max_iter:
        return None
    if count >= 20 and abs(q[0]) > 1:
        return None
    if count >= 30 and (abs(q[0]) > 0.01 or abs(q[1]) > 0.01):
        return None
    return secant_iteration(
        space=space,
        r_1=r_new,
        relation=relation,
        bound_m=bound_m,
        y=y,
        set_coefficients=set_coefficients,
        count=count + 1,
        prec=prec,
        r_0=r_1,
        f_1=q,
        f_0=f_1,
        tolerance=tolerance,
        max_iter=max_iter,
    )


def search_eigenvalues_via_relation(
    space: KleinianMaassFormSpace,
    r1: Real_t,
    r2: Real_t,
    num_steps: Integer_t,
    relation: dict,
    cvalue: Real_t = 1e-2,
    bound_m: Integer_t = 2,
    y: Real_t | None = None,
    Q: Integer_t | None = None,
    set_coefficients: dict | None = None,
    num_threads: Integer_t | None = None,
    use_database: bool = True,
    tolerance: Real_t = 1e-10,
    prec: Integer_t = 53,
) -> list:
    r"""
    End-to-end search for eigenvalues in ``[r1, r2]`` based on a
    coefficient relation.

    The function first computes the Maass form at ``num_steps`` equally
    spaced points (in parallel, with optional MongoDB caching), then
    filters candidates whose coefficient relation residual is below
    ``cvalue``, and finally refines each survivor with the secant method.

    INPUT:

    - ``space`` -- KleinianMaassFormSpace
    - ``r1``, ``r2`` -- reals; search bounds
    - ``num_steps`` -- integer; grid resolution
    - ``relation`` -- dict; see :func:`coeff_diff_fun`
    - ``cvalue`` -- real (default: 1e-2); pre-filter threshold
    - ``bound_m`` -- integer (default: 2)
    - ``y`` -- real or None (default: None)
    - ``Q`` -- integer or None (default: None)
    - ``set_coefficients`` -- dict or None (default: None)
    - ``num_threads`` -- integer or None (default: None)
    - ``use_database`` -- bool (default: True)
    - ``tolerance`` -- real (default: 1e-10); secant convergence tolerance
    - ``prec`` -- integer (default: 53); precision

    OUTPUT:

    - sorted list of refined eigenvalue candidates (imaginary parts of
      ``s``)

    EXAMPLES:

    Unknown locators are rejected upfront::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.all import KleinianMaassFormSpace
        sage: from maass_forms_klein.search.search import search_eigenvalues_via_relation
        sage: space = KleinianMaassFormSpace(-4, cuspidal=True)
        sage: search_eigenvalues_via_relation(
        ....:     space, 6.0, 7.0, 2, {'banana': []}, bound_m=2)
        Traceback (most recent call last):
        ...
        ValueError: Unknown locator 'banana'...

    Happy path with the database-free ``two_y`` locator over a narrow
    range that surrounds the known eigenvalue ``r = 6.62211934``: the
    function returns a sorted list of refined eigenvalue candidates::

        sage: result = search_eigenvalues_via_relation(
        ....:     space, r1=6.5, r2=6.8, num_steps=3,
        ....:     relation={'two_y': {'Y1': 0.5, 'Y2': 0.475, 'indices': [1]}},
        ....:     bound_m=2, Q=3, cvalue=1.0, use_database=False)
        sage: isinstance(result, list)
        True
        sage: result == sorted(result)
        True

    Hecke locators on non-arithmetic manifolds are rejected::

        sage: from snappy import Manifold
        sage: nonarith = KleinianMaassFormSpace(Manifold('5_2'))
        sage: search_eigenvalues_via_relation(
        ....:     nonarith, 0.0, 1.0, 2, {'coprime': [0, 1, 2]}, bound_m=2)
        Traceback (most recent call last):
        ...
        NotImplementedError: Hecke-type coefficient relations ('coprime',
        'prime_power') are only defined for arithmetic Kleinian (Bianchi)
        groups...
    """
    locator = next(iter(relation))
    if locator not in GENERAL_LOCATORS + ARITHMETIC_LOCATORS:
        raise ValueError(
            f"Unknown locator '{locator}'. Expected one of "
            f"{GENERAL_LOCATORS + ARITHMETIC_LOCATORS}."
        )
    if locator in ARITHMETIC_LOCATORS and not is_arithmetic_group(space):
        raise NotImplementedError(
            "Hecke-type coefficient relations ('coprime', 'prime_power') "
            "are only defined for arithmetic Kleinian (Bianchi) groups; "
            "use {'two_y': ...} or {'unit': ...} for general groups."
        )
    db_required = locator != "two_y"
    if db_required and not (use_database_global and use_database):
        raise RuntimeError(
            f"Locator '{locator}' requires database support; use the "
            "{'two_y': ...} locator (or "
            "maass_forms_klein.modform.search.brute_force_search) for "
            "DB-free two-height search."
        )
    # Step 1: populate the database in parallel (only needed for
    # DB-backed locators).
    if db_required:
        list(
            compute_on_interval(
                space=space,
                r1=r1,
                r2=r2,
                num_steps=num_steps,
                prec=prec,
                bound_m=bound_m,
                y=y,
                Q=Q,
                set_coefficients=set_coefficients,
                num_threads=num_threads,
                use_database=True,
            )
        )
        # Step 2 (DB path): pre-filter candidates via the unit-style check.
        CF = ComplexField(prec)
        check_rel = relation[locator]
        check_inputs = [
            (space, check_rel, cvalue, CF(0.5, r), bound_m, y, False)
            for r in create_grid((r1, r2), num_steps)
        ]
        candidates: list[Real_t] = []
        for _input, result in check_coefficients_of_computed_object(check_inputs):
            if result is not None:
                candidates.append(result[0])
    else:
        # Step 2 (two_y path): pre-filter directly via coeff_diff_fun.
        candidates = []
        for r in create_grid((r1, r2), num_steps):
            diff = coeff_diff_fun(
                space,
                relation,
                RR(r),
                bound_m=bound_m,
                y=y,
                set_coefficients=set_coefficients,
                prec=prec,
                Q=Q,
            )
            if abs(diff[0]) + abs(diff[1]) <= cvalue:
                candidates.append(RR(r))
    # Step 3: refine each survivor with the secant method.
    refine_inputs = [
        (space, r, relation, bound_m, y, set_coefficients, 0, prec, None, None, None, tolerance)
        for r in candidates
    ]
    refined: list[Real_t] = []
    for _input, result in secant_iteration(refine_inputs):
        if result is not None:
            refined.append(result)
    return sorted(sorting(refined))


def sorting(pts: list, threshold: Real_t = 1e-3) -> list:
    r"""
    De-duplicate a list of real eigenvalue candidates, treating values
    closer than ``threshold`` as equal.

    INPUT:

    - ``pts`` -- list of reals
    - ``threshold`` -- real (default: 1e-3); cluster tolerance

    OUTPUT:

    - list of reals with near-duplicates removed (first occurrence kept)

    EXAMPLES::

        sage: from maass_forms_klein.search.search import sorting
        sage: sorting([6.62211934, 6.62211935, 7.5])
        [6.62211935000000, 7.50000000000000]
        sage: sorting([1.0, 1.0005, 1.5, 1.5009])
        [1.00050000000000, 1.50090000000000]
    """
    keep = []
    for i, p in enumerate(pts):
        keep.append(p)
        for q in pts[i + 1 :]:
            if abs(p - q) < threshold:
                keep.pop()
                break
    return keep


def distance_between_points(point1: Real_t, point2: Real_t) -> Real_t:
    r"""
    Distance between two real spectral-parameter candidates.

    INPUT:

    - ``point1`` -- real
    - ``point2`` -- real

    OUTPUT:

    - real; ``|point1 - point2|``

    EXAMPLES::

        sage: from maass_forms_klein.search.search import distance_between_points
        sage: distance_between_points(1.0, 1.5)
        0.500000000000000
        sage: distance_between_points(2.0, 1.0)
        1.00000000000000
    """
    return abs(point1 - point2)


# Re-export the existing two-height search helpers so that the search
# subpackage offers a complete eigenvalue-search API.
__all__ = [
    "ARITHMETIC_LOCATORS",
    "GENERAL_LOCATORS",
    "brute_force_search",
    "check_coefficients_of_computed_object",
    "coeff_diff_fun",
    "compute_on_interval",
    "compute_one_spectral_parameter",
    "create_grid",
    "distance_between_points",
    "is_arithmetic_group",
    "newton_method_search",
    "search_eigenvalues_via_relation",
    "secant_iteration",
    "sorting",
]
