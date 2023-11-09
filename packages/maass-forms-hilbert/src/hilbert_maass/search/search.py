"""
Routines to search for Hilbert Maass forms
"""
import logging
import os

import mongoengine
import numpy
from comp_manager.utils import insert_object, load_object
from hilbert_maass.database.models import HilbertMaassFormDB
from hilbert_maass.modform.coefficients import get_pb_pts_set_params
from hilbert_maass.modform.hilbert_maass_element import HilbertMaassForm
from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
from hilbert_maass.modform.utils import Integer_t, complex_tuple_to_json, Real_t, Complex_t
from sage.categories.sets_cat import cartesian_product
from sage.functions.other import ceil
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.parallel.decorate import parallel
from sage.rings.complex_mpfr import ComplexField
from sage.rings.integer import Integer
from sage.rings.real_mpfr import RealNumber as RealNumber_class


log = logging.getLogger(__name__)


def create_grid(grid_limits: tuple[tuple[Real_t]],
                    grid_numbers: tuple[Integer_t]):
    """
    Compute a Hilbert Maass form on a grid.

    INPUT:

    - ``space`` -- space of Hilbert Maass forms
    - ``grid_limits`` -- tuple of tuples to represent the boundary of the grid
    - ``grid_numbers`` -- tuple of integers to represent the number of grid points in each dimension

    EXAMPLES::

        sage: from hilbert_maass.all import HilbertMaassFormSpace
        sage: from hilbert_maass.search.search import create_grid, compute_on_grid
        sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
        sage: create_grid(((0,1),(0,1)), (2,2))
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
    grids = numpy.meshgrid(
        *(numpy.linspace(x0, x1, grid_numbers[i]) for i, (x0, x1) in enumerate(grid_limits)))
    return grids, cartesian_product([range(grids[0].shape[i]) for i in range(len(grids[0].shape))])


def compute_on_grid(space: HilbertMaassFormSpace, grid_limits: tuple[tuple[Real_t]],
                    grid_numbers: tuple[Integer_t], prec: Integer_t = 53,
                    bound_m: tuple[tuple[Integer_t]] | Integer_t = 2,
                    num_threads: Integer_t = None):
    """
    Compute a Hilbert Maass form on a grid.

    INPUT:

    - ``space`` -- space of Hilbert Maass forms
    - ``grid_limits`` -- tuple of tuples to represent the boundary of the grid
    - ``grid_numbers`` -- tuple of integers to represent the number of grid points in each dimension

    EXAMPLES::

        sage: from hilbert_maass.all import HilbertMaassFormSpace
        sage: from hilbert_maass.search.search import compute_on_grid
        sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
        sage: result = list(compute_on_grid(H, ((0,1),(0,1)), (2,2), bound_m=1, num_threads=1))
        sage: len(result)
        4
        sage: spectral_parameters = [result[x][0][0][1] for x in range(4)]
        sage: spectral_parameters.sort()
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
    grids, grid_indices = create_grid(grid_limits, grid_numbers)
    CF = ComplexField(prec)
    input_params = []
    for m in grid_indices:
        spectral_parameter = tuple(CF(0.5, grid[tuple(m)]) for grid in grids)
        input_params.append((space, spectral_parameter, bound_m))
    if num_threads is not None:
        os.environ['SAGE_NUM_THREADS'] = str(num_threads)
    # Prepare the cache to avoid race errors
    smax = 0
    for r in input_params:
        smax = max(smax, max([abs(x) for x in r[1]]))
    for si in range(0, ceil(smax)+1):
        get_pb_pts_set_params(space, M=input_params[0][2], smax = si)
    return compute_one_spectral_parameter(input_params)


@parallel()
def compute_one_spectral_parameter(space: HilbertMaassFormSpace,
                                   spectral_parameter: tuple[Complex_t],
                                   bound_m: tuple[tuple[Integer_t]] | Integer_t):
    try:
        maass_form_db = HilbertMaassFormDB.near_or_create(
            parent=space.to_json(),
            spectral_parameter=spectral_parameter,
            bound_m=bound_m)
        maass_form = load_object(maass_form_db)
    except mongoengine.connection.ConnectionFailure:
        log.warning(f"Could not connect to database. Compute locally only")
        maass_form = HilbertMaassForm(space, spectral_parameter)
        maass_form.compute_coefficients(M=bound_m)
    if not maass_form.coefficients():
        maass_form.compute_coefficients(M=bound_m)
        insert_object(maass_form)
        log.debug(f"Computed Hilbert Maass form for s={spectral_parameter}")
    return maass_form

def broyden_iteration(previous_iterations: list):
    """
    Basic implementation of Broyden's method to find solution of 2x2 system of equations:
    a x + b y = f
    c x + d y = g


    INPUT:

    - previous_iterations -- tuple of length 3 consisting of tuples of length 4 (x,y,f,g).

    We use 3 points represented at tuples of length 4:
    v0 = (x_{n-1}, y_{n-1}, f_{n-1}, f_{n-1})
    v0 = (x_{n-1}, y_{n-1}, f_{n-1}, f_{n-1})
    v0 = (x_{n-1}, y_{n-1}, f_{n-1}, f_{n-1})
    v0 = (x_{n-1}, y_{n-1}, f_{n-1}, f_{n-1})

    Step n-1 = (x0,y0) step n = (x1,y1)

    NOTE: See https://en.wikipedia.org/wiki/Broyden%27s_method
    """
    if not isinstance(previous_iterations, tuple) or len(previous_iterations) != 3:
        raise ValueError("Input should be a tuple of lenght 3")
    v0, v1, v2 = previous_iterations
    delta_x0 = vector([v1[0] - v0[0], v1[1] - v0[1]])
    delta_y0 = vector([v1[2] - v0[2], v1[3] - v0[3]])
    delta_x1 = vector([v2[0] - v1[0], v2[1] - v1[1]])
    delta_y1 = vector([v2[2] - v1[2], v2[3] - v1[3]])
    # Initial approximation of the Jacobian
    def jacobian_approximation(x, f):
        return matrix([[f[0] / x[0], 0], [0, f[1] / x[1]]])
    J0 = jacobian_approximation(delta_x0, delta_y0)
    # Finite difference
    delta_J = (delta_y1 - J0 * delta_x1) / delta_x1.norm(2) ** 2 * delta_x1
    print("diff=", delta_y1 - J0 * delta_x1, delta_J)
    J1 = J0 + delta_J
    print("J1=", J1)
    x1 = vector([v2[0], v2[1]])
    f1 = vector([v2[2], v2[3]])
    x_new = x1 - J1 * f1
    print("x new=", x_new)
    return x_new
