"""
Routines to search for Hilbert Maass forms
"""
import logging
import os
from cmath import asin

import mongoengine
import numpy

from hilbert_maass.modform.hilbert_maass_element import HilbertMaassForm
from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
from hilbert_maass.modform.compute_coefficients import get_pb_pts_set_params
from hilbert_maass.modform.utils import Integer_t, Real_t, Complex_t
from comp_manager.utils import insert_object, load_object
from hilbert_maass.database.models import HilbertMaassFormDB
from sage.all import sin, cos, pi
from sage.functions.other import ceil
from sage.categories.sets_cat import cartesian_product
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.parallel.decorate import parallel
from sage.rings.complex_mpfr import ComplexField
from sage.rings.integer import Integer

log = logging.getLogger(__name__)


def create_parallel_grid(grid_limits: tuple[Real_t], grid_number: Integer_t):
    """
     INPUT:
        - ``grid_limits`` -- tuple of Real Number
        - ``grid_numbers`` -- Integer

    EXAMPLES:

    sage: from hilbert_maass.search.search import create_parallel_grid
    sage: create_parallel_grid((0, 2), 11)
    array([0. , 0.2, 0.4, 0.6, 0.8, 1. , 1.2, 1.4, 1.6, 1.8, 2. ])
    """
    grids = numpy.linspace(grid_limits[0], grid_limits[1], grid_number)
    return grids


# Method 1
# For creating the grid points. It will work for real field of any degree.
def create_grid_non_circular(grid_limits: tuple[tuple[Real_t]],
                             grid_numbers: tuple[Integer_t]):
    """
    Create non-circular grid points.

    INPUT:

    - ``space`` -- space of Hilbert Maass forms
    - ``grid_limits`` -- tuple of tuples to represent the boundary of the grid
    - ``grid_numbers`` -- tuple of integers to represent the number of grid points in each dimension

    EXAMPLES::

        sage: from hilbert_maass.all import HilbertMaassFormSpace
        sage: from hilbert_maass.search.search import create_grid_non_circular
        sage: from hilbert_maass.search.search import compute_on_non_circular_grid
        sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
        sage: create_grid_non_circular(((0, 1),(0,1)), (2,2))
        ([array([[0., 1.],
                 [0., 1.]]),
          array([[0., 0.],
                 [1., 1.]])],
         The Cartesian product of ({0, 1}, {0, 1}))
    """
    # Check that sizes match
    if len(grid_limits) != len(grid_numbers):
        raise ValueError("Number of grid limits does not match number of grid numbers")

    # Set up the grid
    st = []
    for i in range(0, len(grid_limits)):
        st.append(create_parallel_grid(grid_limits[i], grid_numbers[i]))
    grids = numpy.meshgrid(*(st[i] for i in range(len(st))))
    return grids, cartesian_product([range(grids[0].shape[i]) for i in range(len(grids[0].shape))])


def compute_on_non_circular_grid(space: HilbertMaassFormSpace, grid_limits: tuple[tuple[Real_t]],
                                 grid_numbers: tuple[Integer_t], prec: Integer_t = 53,
                                 bound_m: tuple[tuple[Integer_t]] | Integer_t = 2,
                                 y: tuple[Real_t] | None = None,
                                 num_threads: Integer_t = None,
                                 spectral_symmetry: bool = True,
                                 set_coefficients: dict = None,
                                 use_database: bool = True,
                                 spectral_epsilon: Real_t = 1e-2):
    """
    Compute a Hilbert Maass form on the non-circular grid.

    INPUT:

    - ``space`` -- space of Hilbert Maass forms
    - ``grid_limits`` -- tuple of tuples to represent the boundary of the grid
    - ``grid_numbers`` -- tuple of integers to represent the number of grid points in each dimension
    - ``spectral_symmetry`` -- (default=True) only calculate for spectral parameters (r1,r2) with r1 <= r2 + eps
    - ``spectral_epsilon`` -- eps (see above)
    - ``num_threads`` -- integer (default=1) number of threads to use,
                         need to be less than or equal to SAGE_NUM_THREADS
    - ``use_database`` -- boolean(default: True) use the database

    EXAMPLES::

        sage: from hilbert_maass.all import HilbertMaassFormSpace
        sage: from hilbert_maass.search.search import compute_on_non_circular_grid
        sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal = False)
        sage: result = list(compute_on_non_circular_grid(H, ((0, 1),(0, 1)), (2,2), bound_m=1, num_threads=1, spectral_symmetry=True, use_database=False))
        sage: len(result)
        3
        sage: result = list(compute_on_non_circular_grid(H, ((0,1),(0,1)), (2,2), bound_m=1, num_threads=1,
        ....: spectral_symmetry=False, use_database=False))
        sage: len(result)
        4
        sage: spectral_parameters = [result[x][0][0][1] for x in range(4)]
        sage: sorted(spectral_parameters)
        [(0.500000000000000, 0.500000000000000),
         (0.500000000000000, 0.500000000000000 + 1.00000000000000*I),
         (0.500000000000000 + 1.00000000000000*I, 0.500000000000000),
         (0.500000000000000 + 1.00000000000000*I,
          0.500000000000000 + 1.00000000000000*I)]
        sage: result[0][1]
        Hilbert Maass form for HilbertMaassFormSpace(Hilbert Modular Group PSL(2) over Maximal...
    """
    if len(grid_limits) != len(grid_numbers):
        raise ValueError("Number of grid limits does not match number of grid numbers")
    if len(grid_numbers) != space.number_field().absolute_degree():
        raise ValueError("Number of number of grid points does not match number field degree")
    if isinstance(bound_m, (Integer, int)):
        bound_m = tuple([(-bound_m, bound_m)] * space.number_field().absolute_degree())
    grids, grid_indices = create_grid_non_circular(grid_limits, grid_numbers)
    CF = ComplexField(prec)
    input_params = []
    for m in grid_indices:
        spectral_parameter = tuple(CF(0.5, grid[tuple(m)]) for grid in grids)
        log.debug(f"Computing spectral parameter {spectral_parameter}")
        log.debug(
            f"Ineqs:{[spectral_parameter[i + 1].imag() > spectral_parameter[i].imag() + spectral_epsilon for i in range(len(spectral_parameter) - 1)]}")
        if spectral_symmetry and any(spectral_parameter[i + 1].imag() > spectral_parameter[i].imag() +
                                     spectral_epsilon for i in range(len(spectral_parameter) - 1)):
            log.debug(f"Skipping spectral parameter {spectral_parameter}")
            continue

        input_params.append((space, spectral_parameter, bound_m, y,
                             set_coefficients, use_database))
    if num_threads is not None:
        os.environ['SAGE_NUM_THREADS'] = str(num_threads)
    # Prepare the cache to avoid race errors
    smax = 0
    for r in input_params:
        smax = max(smax, max([abs(x) for x in r[1]]))
    for si in range(0, ceil(smax) + 1):
        get_pb_pts_set_params(space, M=input_params[0][2], Y=y, smax=si)
    return compute_one_spectral_parameter(input_params)


# Method 2 for creating grid points.
# We will work with this one while working with real field of degree 2.
def create_grid_on_circumference_from_r1_to_r2(r1: Real_t,
                                               r2: Real_t,
                                               min_distance: Real_t = 0.01,
                                               prec: Integer_t = 53):
    r"""
       Create the grid (tuple) point on the arcs of angle (0, pi/4) and different radiuses r lying between
       $r1 \le r \le r2$

        INPUT:
         - ``r1`` --  Real Number
         - ``r2`` --  Real Number
         - ``min_distance`` -- Real Number
         - ``prec``  --- Integer_t

       EXAMPLES::
       sage: from hilbert_maass.search.search import create_grid_on_circumference_from_r1_to_r2
       sage: grid_points=create_grid_on_circumference_from_r1_to_r2(r1=0.3, r2=0.4, min_distance=0.1)
       sage: grid_points
       [(0.300000000000000, 0.000000000000000),
        (0.294235584120969, 0.0585270966048385),
        (0.277163859753386, 0.114805029709527),
        (0.249440883690764, 0.166671069905881),
        (0.212132034355964, 0.212132034355964)]

     """

    CF = ComplexField(prec)
    grid_limits_on_angle = (0, CF(pi / 4))
    grid_numbers_on_radiuses = int(round(float((r2 - r1) / min_distance)))
    grids_on_r1_to_r2 = create_parallel_grid((r1, r2), grid_numbers_on_radiuses)
    CF = ComplexField(prec)
    center = (0, 0)
    grid_points = []
    for r in grids_on_r1_to_r2:
        if r == 0:
            grid = center
            grid_points.append(grid)
        else:
            grid_numbers_on_angle = int(round(float(pi * r / (2 * min_distance))))
            grids_on_angle = create_parallel_grid(grid_limits_on_angle, grid_numbers_on_angle)
            for m in grids_on_angle:
                grid = (CF(center[0] + r * cos(m)), CF(center[1] + r * sin(m)))
                grid_points.append(grid)
    return grid_points


def compute_on_circumference_grid(space: HilbertMaassFormSpace,
                                  r1: Real_t, r2: Real_t,
                                  min_distance: Real_t = 0.01,
                                  prec: Integer_t = 53,
                                  bound_m: tuple[tuple[Integer_t]] | Integer_t = 2,
                                  y: tuple[Real_t] | None = None,
                                  set_coefficients: dict = None,
                                  num_threads: Integer_t = None,
                                  spectral_symmetry: bool = True,
                                  use_database: bool = True,
                                  spectral_epsilon: Real_t = 1e-2):
    """
    Compute Hilbert Maass form type objects on grid points.

    INPUT:

    - ``space`` -- space of Hilbert Maass forms
    - ``r1`` -- initial radius
    - ``r2`` -- final radius
    - ``min_distance`` -- IF g= (a, b) is a grid point then there exist another grid g'=(a', b')
            whose distance from g will be less than or equal to min_distance + epsilon for any epsilon >0.
    - ``num_threads`` -- integer (default=1) number of threads to use,
                         need to be less than or equal to SAGE_NUM_THREADS
    - ``use_database`` -- boolean(default: True) use the database


    EXAMPLES::

        sage: from hilbert_maass.all import HilbertMaassFormSpace
        sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=True)
        sage: from hilbert_maass.search.search import compute_on_circumference_grid
        sage: from comp_manager import create_app
        sage: import logging
        sage: logging.captureWarnings(True)
        sage: result = list(compute_on_circumference_grid(H, r1=0.3, r2=0.4, min_distance= 0.1, bound_m=1,
        ....: num_threads=1, use_database= False))
        sage: len(result)
        5
        sage: spectral_parameters = [result[x][0][0][1] for x in range(len(result))]
        sage: sorted(spectral_parameters)
        [(0.500000000000000 + 0.212132034355964*I,
         0.500000000000000 + 0.212132034355964*I),
         (0.500000000000000 + 0.249440883690764*I,
         0.500000000000000 + 0.166671069905881*I),
         (0.500000000000000 + 0.277163859753386*I,
         0.500000000000000 + 0.114805029709527*I),
         (0.500000000000000 + 0.294235584120969*I,
         0.500000000000000 + 0.0585270966048385*I),
         (0.500000000000000 + 0.300000000000000*I, 0.500000000000000)]
        sage: result[0][0][0][2]
        ((-1, 1), (-1, 1))
        sage: result[0][0][0][3]

        sage: result[0][1]
        Hilbert Maass form for HilbertMaassFormSpace(Hilbert Modular Group PSL(2) over Maximal...
    """
    if isinstance(bound_m, (Integer, int)):
        bound_m = tuple([(-bound_m, bound_m)] * space.number_field().absolute_degree())
    grid_points = create_grid_on_circumference_from_r1_to_r2(r1, r2, min_distance, prec)
    input_params = []
    CF = ComplexField(prec)
    for r in grid_points:
        spectral_parameter = (CF(0.5, r[0]), CF(0.5, r[1]))
        if spectral_symmetry and any(spectral_parameter[i + 1].imag() > spectral_parameter[i].imag() +
                                     spectral_epsilon for i in range(len(spectral_parameter) - 1)):
            log.debug(f"Skipping spectral parameter {spectral_parameter}")
            continue
        log.debug(f"Computing spectral parameter {spectral_parameter}")
        input_params.append((space, spectral_parameter, bound_m, y,
                             set_coefficients, use_database))
    if num_threads is not None:
        os.environ['SAGE_NUM_THREADS'] = str(num_threads)
    smax = 0
    for r in input_params:
        smax = max(smax, max([abs(x) for x in r[1]]))
    for si in range(0, ceil(smax) + 1):
        get_pb_pts_set_params(space, M=input_params[0][2], Y=y, smax=si)
    return compute_one_spectral_parameter(input_params)


@parallel()
def compute_one_spectral_parameter(space: HilbertMaassFormSpace,
                                   spectral_parameter: tuple[Complex_t],
                                   bound_m: tuple[tuple[Integer_t]] | Integer_t,
                                   y: tuple[Real_t] | None = None,
                                   set_coefficients: dict = None,
                                   use_database: bool = True):
    """

    INPUT:
        - ``space`` -- HilbertMaassFormSpace
        - ``spectral_parameter`` -- tuple of complex numbers
        - ``bound_m`` -- tuple of tuples of integers or integer
        - ``y`` -- tuple of real numbers
        - ``set_coefficients`` -- dict (coefficients to set and values)
        - ``use_database`` -- bool (default: True) set to False to not use database.

    EXAMPLES::

    sage: from hilbert_maass.search.search import compute_one_spectral_parameter
    sage: from hilbert_maass.all import HilbertMaassFormSpace
    sage: space = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
    sage: # Define the input parameters
    sage: spectral_parameter = (0.5 + 0.5j, 0.5 + 1j)
    sage: bound_m = ((-1, 1), (-1, 1))
    sage: y = (0.1, 0.2)
    sage: set_coefficients = {(0,1): 1, (0,0): 2}
    sage: use_database = False
    sage: compute_one_spectral_parameter(space, spectral_parameter, bound_m, y,
    ....: set_coefficients, use_database)
    Hilbert Maass form for HilbertMaassFormSpace(Hilbert Modular Group PSL(2) over Maximal Order...
    """
    maass_form = None
    if use_database:
        try:
            maass_form_db = HilbertMaassFormDB.near_or_create(
                parent=space.to_json(),
                spectral_parameter=spectral_parameter,
                bound_m=bound_m,
                set_coefficients=set_coefficients,
                y=y)
            maass_form = load_object(maass_form_db)
        except mongoengine.connection.ConnectionFailure:
            log.warning(f"Could not connect to database. Compute locally only")
    if not maass_form:
        maass_form = HilbertMaassForm(space, spectral_parameter)
    if not maass_form.coefficients():
        maass_form.compute_coefficients(M=bound_m, set_coefficients=set_coefficients,
                                        Y=y)
        if use_database:
            insert_object(maass_form)
        log.debug(f"Computed Hilbert Maass form for s={spectral_parameter}")
    return maass_form


@parallel
def check_coefficients_of_computed_object(space: HilbertMaassFormSpace,
                                          check_rel: list,
                                          cvalue: Real_t,
                                          spectral_parameter: tuple[Complex_t],
                                          bound_m: tuple[tuple[Integer_t]],
                                          y: tuple[Real_t] | None = None):
    """
    After computing object in the database, we use this function to check the object which are near to
    an eigenvalue. We need to do this as Broyden method work fine if the initial point is close to
    an eigenvalue.

    INPUT:
        - ``space`` -- HilbertMaassFormSpace
        - ``check_rel`` -- unit relation. Check the unit relation.
        - ``cvalue`` -- Real number
        - ``spectral_parameter`` -- tuple of complex numbers
        - ``bound_m`` -- tuple of tuples of integers
        - ``y`` -- tuple of real numbers

    Examples::
    Note: We have removed sage below as it requires the use of database.
    from hilbert_maass.all import HilbertMaassFormSpace
    H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
    from hilbert_maass.search.search import compute_on_circumference_grid
    result = list(compute_on_circumference_grid(H, r1=0.3, r2=0.4, min_distance= 0.1, bound_m=1,
    ....: num_threads=1))
    ax = []
    for i in range (0, len(result)):
         ax.append((H, [(1, -1), (1, 1)], 0.1, result[i][0][0][1], result[i][0][0][2], result[i][0][0][3]))
    from hilbert_maass.search.search import check_coefficients_of_computed_object
    result1 = list(check_coefficients_of_computed_object(ax))
    sp = [result1[x][1] for x in range(len(result1))]
    [None, None, None, None, (0.212132034355964, 0.212132034355964)]
    """
    f = load_object(list(HilbertMaassFormDB.objects.space(space).with_m_precision(bound_m)
                         .near(spectral_parameter).with_y_precision(y))[0])
    t = abs((f.coefficients()[check_rel[1]] - f.coefficients()[check_rel[0]]).real())
    if t <= cvalue:
        print(tuple((spectral_parameter[0].imag(), spectral_parameter[1].imag())), t)
        return tuple((spectral_parameter[0].imag(), spectral_parameter[1].imag()))


@parallel()
def broyden_iteration(space: HilbertMaassFormSpace, r_1: tuple,
                      relation: dict,
                      bound_m: tuple[tuple[Integer_t]],
                      y: tuple[Real_t] | None = None,
                      set_coefficients: dict | None = None,
                      count: Integer_t = 0,
                      prec=53,
                      r_0: tuple = None,
                      f_1: tuple = None,
                      f_0: tuple = None,
                      J_0: matrix = None):
    """

    NOTE: See https://en.wikipedia.org/wiki/Broyden%27s_method for more information.

    INPUT:
        - ``space`` -- HilbertMaassFormSpace
        - ``r_1`` -- tuple[Real_t]. Starting point.
        - ``relation`` -- dictionary. Ex= {'unit': [(1, -1),(1, 1)]}
        - ``bound_m`` -- tuple of tuples of integers or integer
    Examples::
    Note: We have removed sage below as it requires the use of database.
    from hilbert_maass.all import HilbertMaassFormSpace
    from hilbert_maass.search.search import broyden_iteration
    from sage.rings.complex_mpfr import ComplexField
    H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
    bx = [(H, ( 2.23, 2.21), {'unit' : [(1, -1), (1, 1)] }, ((-1, 1) ,(-1, 1)))]
    result2 = list(broyden_iteration(bx))
    result2[0][1]
    (2.21206924454105, 2.19603006584595)
    """
    CF = ComplexField(prec)
    if (r_0 is None):
        r_0 = vector(r_1) - vector((CF(0.0000001), CF(0.0000001)))
        f_0 = coeff_diff_fun(space, relation, r_0, bound_m, y, set_coefficients)
    if (f_1 is None):
        f_1 = coeff_diff_fun(space, relation, r_1, bound_m, y, set_coefficients)
    delta__r_1 = vector([CF(r_1[0]) - CF(r_0[0]), CF(r_1[1]) - CF(r_0[1])])
    delta__f_1 = vector([CF(f_1[0]) - CF(f_0[0]), CF(f_1[1]) - CF(f_0[1])])
    if (J_0 is None):
        if delta__r_1[0] == 0 or delta__r_1[1] == 0:
            s = 'Zero divisor'
            raise ValueError("J_0 became zero. Try with some other entry.")
        J_0 = matrix([[delta__f_1[0] / delta__r_1[0], 0], [0, delta__f_1[1] / delta__r_1[1]]])
    delta__J = ((delta__f_1 - J_0 * delta__r_1) / delta__r_1.norm(2) ** 2).column() * (delta__r_1).row()
    J_1 = J_0 + delta__J
    if (J_1.det() == 0):
        raise ValueError("J_1 became zero. Try with some other entry.")
    r_1 = vector(r_1)
    f_1 = vector(f_1)
    r_new = r_1 - J_1 ** (-1) * f_1  # It should be  r_new = r1 - J1**(-1) * f1.
    if (r_new[0] < 0 or r_new[1] < 0):
        return None
    elif ((r_new[0] ** 2 + r_new[1] ** 2) > 100):
        return None
    q = coeff_diff_fun(space, relation, r_new, bound_m, y, set_coefficients)
    print('r_new=', r_new, 'value=', q)
    if (abs(q[0]) < 1e-12 and abs(q[1]) < 1e-12):
        return r_new
    elif (count >= 20 and (abs(q[0]) > 1)):
        return None
    elif (count >= 25 and (abs(q[0]) > 0.1)):
        return None
    elif (count >= 30 and (abs(q[0]) > 0.01 or abs(q[1]) > 0.01)):
        return None
    elif (count > 40):
        return None
    else:
        return broyden_iteration(space=space, r_1=r_new, relation=relation, bound_m=bound_m, y=y,
                                 set_coefficients=set_coefficients, count=count + 1, r_0=r_1, f_1=q,
                                 f_0=f_1, J_0=J_1)


def coeff_diff_fun(space, relation: dict,
                   r: tuple,
                   bound_m: tuple[tuple[Integer_t]] | Integer_t = None,
                   y: tuple[Real_t] | None = None,
                   set_coefficients: dict | None = None,
                   prec=53):
    """
    INPUT:
        - ``space`` -- HilbertMaassFormSpace
        - ``relation`` -- dictionary. Ex= {'unit': [(1, -1),(1, 1)]}
        - ``r` -- tuple[Real_t].
        - ``bound_m`` -- tuple of tuples of integers or integer
    Example::
    Note: We have removed sage below as it requires the use of database.
    from hilbert_maass.search.search import coeff_diff_fun
    from comp_manager import create_app
    ap = create_app()
    from hilbert_maass.all import HilbertMaassFormSpace
    H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
    from hilbert_maass.modform.utils import best_hecke_relation
    best_hecke_relation(H, 1, 'unit')
    [[(-1, -1), (-1, 1)]]
    coeff_diff_fun(H, {'unit': [(1, -1), (1, 1)]}, (2.21206924454105, 2.19603006584595), ((-1, 1),(-1, 1)))
    (9.99200722162641e-16, 0.000000000000000)
    coeff_diff_fun(H, {'unit': [(-1, -1), (-1, 1)]}, (2.21206924454105, 2.19603006584595), ((-1, 1), (-1, 1)))
    (0.000000000000000, 0.000000000000000)
    best_hecke_relation(H, 1, 'coprime')
    []
    """

    CF = ComplexField(prec)
    spectral_parameter = (CF(0.5, r[0]), CF(0.5, r[1]))
    t = len(list(HilbertMaassFormDB.objects.with_m_precision(bound_m)
                 .with_y_precision(y).near(spectral_parameter)))
    if (t == 1):
        f = load_object(list(HilbertMaassFormDB.objects.near(spectral_parameter)
                             .with_m_precision(bound_m).with_y_precision(y))[0])
    else:
        result1 = list(compute_one_spectral_parameter([(space, spectral_parameter,
                                                        bound_m, y, set_coefficients)]))
        f = load_object(list(HilbertMaassFormDB.objects.near(spectral_parameter)
                             .with_m_precision(bound_m).with_y_precision(y))[0])
    if 'unit' in relation:
        u = relation['unit']
        coef1 = abs((f.coefficients()[u[1]] - f.coefficients()[u[0]]).real())
        coef2 = abs((f.coefficients()[u[1]] - f.coefficients()[u[0]]).imag())
    elif 'coprime' in relation:
        c = relation['coprime']
        coef1 = abs((f.coefficients()[c[2]] - f.coefficients()[c[1]] * f.coefficients()[c[0]]).real())
        coef2 = abs((f.coefficients()[c[2]] - f.coefficients()[c[1]] * f.coefficients()[c[0]]).imag())
    else:
        c = relation['prime_power']
        t = len(c)
        if t == 2:
            coef1 = abs((f.coefficients()[c[1]] - f.coefficients()[c[0]] * f.coefficients()[c[0]] - 1).real())
            coef2 = abs((f.coefficients()[c[1]] - f.coefficients()[c[0]] * f.coefficients()[c[0]] - 1).imag())
        else:
            coef1 = abs((f.coefficients()[c[t - 1]] - f.coefficients()[c[0]] * f.coefficients()[c[t - 2]]
                         + f.coefficients()[c[t - 3]]).real())
            coef2 = abs((f.coefficients()[c[t - 1]] - f.coefficients()[c[0]] * f.coefficients()[c[t - 2]]
                         + f.coefficients()[c[t - 3]]).imag())

    return tuple((coef1, coef2))


def sorting(pts: list, prec=53):
    """
    INPUT:
        - ``pts`` -- list of tuples of Real number

    Example::
    sage: from hilbert_maass.search.search import sorting
    sage: cx = [(0.1, 0.5), (0.1001, 0.50001)]
    sage: sorting(cx)
    [(0.100100000000000, 0.500010000000000)]
    """

    ats = []
    for i in range(len(pts)):
        ats.append(pts[i])
        for j in range(i + 1, len(pts)):
            s = distance_between_points_two(pts[i], pts[j])
            if abs(s[0]) < 0.001:
                ats.pop()  # Remove the last added element if distance condition is met
                break

    return ats


def distance_between_points_two(point1: tuple, point2: tuple, prec=53):
    """
    INPUT:
        - ``point1`` -- tuple of Real Number
        - ``point2`` -- tuple of Real Number
    sage: from hilbert_maass.search.search import distance_between_points_two
    sage: distance_between_points_two((5.29138359828915, 3.66832128742010),(5.29138359828915, 3.66832128742007))
    (2.97539770599542e-14, -1.57079632679490)
    """
    x1, y1 = point1
    x2, y2 = point2
    CF = ComplexField(prec)
    d = CF(numpy.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2))
    theta = CF(asin((y2 - y1) / d))

    return d, theta