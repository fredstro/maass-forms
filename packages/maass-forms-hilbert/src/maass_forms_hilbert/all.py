import contextlib

# Need to import some symbols explicitly for command line imports when using Passagemath
with contextlib.suppress(ImportError):
    from sage.all__sagemath_symbolics import Category
try:
    from maass_forms_hilbert.database.models import HilbertMaassFormDB
except ImportError as e:
    import logging

    logging.error(f"Cannot import HilbertMaassFormDB: {e}")
try:
    from maass_form_core.functions.bessel.besselk_dp import besselk_dp
except ImportError:
    besselk_dp = None
from maass_forms_hilbert.functions.functions import bessel_prod, exp_trace_prod
from maass_forms_hilbert.modform.coefficients import HilbertMaassCoefficients
from maass_forms_hilbert.modform.compute_coefficients import compute_coefficients, get_pb_pts
from maass_forms_hilbert.modform.eisenstein_series import HilbertEisensteinSeries
from maass_forms_hilbert.modform.hilbert_maass_element import (
    HilbertMaassForm,
    HilbertMaassForm_Element,
)
from maass_forms_hilbert.modform.hilbert_maass_space import HilbertMaassFormSpace
from maass_forms_hilbert.modform.utils import (
    get_Q_from_bounds,
    ideal_basis_matrix,
    ideal_coordinates,
    map_int_to_tuple,
    map_tuple_to_int,
)
from maass_forms_hilbert.search.search import (
    compute_on_non_circular_grid,
    compute_one_spectral_parameter,
    create_grid_non_circular,
)
