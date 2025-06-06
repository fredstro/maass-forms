try:
    # Need to import some symbols explicitly for command line imports when using Passagemath
    from sage.all__sagemath_symbolics import Category
except ImportError:
    pass
try:
    from .database.models import HilbertMaassFormDB
except ImportError as e:
    import logging

    logging.error(f"Cannot import HilbertMaassFormDB: {e}")
from .modform.utils import (
    get_Q_from_bounds,
    map_int_to_tuple,
    map_tuple_to_int,
    ideal_coordinates,
    ideal_basis_matrix,
)
from .functions.functions import bessel_prod, exp_trace_prod
from .functions.bessel.besselk_dp import besselk_dp
from .modform.eisenstein_series import HilbertEisensteinSeries
from .modform.hilbert_maass_space import HilbertMaassFormSpace
from .modform.hilbert_maass_element import HilbertMaassForm_Element, HilbertMaassForm
from .modform.coefficients import HilbertMaassCoefficients
from .modform.compute_coefficients import get_pb_pts, compute_coefficients
from .search.search import (
    compute_on_non_circular_grid,
    create_grid_non_circular,
    compute_one_spectral_parameter,
)
