"""
Routines to search for Hilbert Maass forms
"""
import logging
import os

import numpy
from comp_manager.utils import insert_object, load_object
from hilbert_maass.database.models import HilbertMaassFormDB
from hilbert_maass.modform.hilbert_maass_element import HilbertMaassForm
from hilbert_maass.modform.utils import Integer_t, complex_tuple_to_json, Real_t
from sage.categories.sets_cat import cartesian_product
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
        sage: from hilbert_maass.search.search import compute_on_grid
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


def compute_on_grid(space, grid_limits: tuple[tuple[Real_t]],
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
        sage: compute_on_grid(H, ((0,1),(0,1)), (2,2))

    """
    if len(grid_limits) != len(grid_numbers):
        raise ValueError("Number of grid limits does not match number of grid numbers")
    if len(grid_numbers) != space.number_field().absolute_degree():
        raise ValueError("Number of number of grid points does not match number field degree")
    if isinstance(bound_m, (Integer, int)):
        bound_m = [(-bound_m, bound_m)] * space.number_field().absolute_degree()
    grids, grid_indices = create_grid(grid_limits, grid_numbers)
    CF = ComplexField(prec)
    input_params = []
    for m in grid_indices:
        spectral_parameter = tuple(CF(0.5, grid[tuple(m)]) for grid in grids)
        input_params.append((space, spectral_parameter, bound_m))
    if num_threads is not None:
        os.environ['SAGE_NUM_THREADS'] = str(num_threads)
    compute_one_spectral_parameter(input_params)


@parallel()
def compute_one_spectral_parameter(space, spectral_parameter, bound_m):
    maass_form_db = HilbertMaassFormDB.near_or_create(
        parent=space.to_json(),
        spectral_parameter=spectral_parameter,
        bound_m=bound_m)
    if not maass_form_db.coefficients:
        maass_form = load_object(maass_form_db)
        maass_form.compute_coefficients(M=bound_m)
        insert_object(maass_form)
        log.debug(f"Computed Hilbert Maass form for s={spectral_parameter}")
