import logging

from comp_manager.utils import load_object
from mongoengine import get_connection
from sage.rings.cc import CC
from sage.rings.number_field.number_field import QuadraticField

from . import dbconnect

dbconnect()
conn = get_connection("default")

log = logging.getLogger(__name__)


def test_check_coefficients_of_computed_object(caplog):
    from maass_forms_hilbert.all import HilbertMaassFormSpace
    from maass_forms_hilbert.database.models import HilbertMaassFormDB
    from maass_forms_hilbert.search.search import check_coefficients_of_computed_object

    H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
    spectral_parameter = CC(0.5, 0.3), CC(0.5, 0)
    bound_m = (-1, 1), (-1, 1)
    maass_form_db = HilbertMaassFormDB.near_or_create(
        parent=H.to_json(),
        spectral_parameter=spectral_parameter,
        bound_m=bound_m,
        set_coefficients=None,
        y=None,
        q=None,
    )
    f = (H, [(1, -1), (1, 1)], 0.1, spectral_parameter, bound_m, None)
    result1 = list(check_coefficients_of_computed_object([f]))
    log.critical("result1: %s", result1)
    sp = [result1[x][1] for x in range(len(result1))]
    assert sp == [None]
    spectral_parameter = (
        CC(0.500000000000000, 0.212132034355964),
        CC(0.500000000000000, 0.212132034355964),
    )
    maass_form_db = HilbertMaassFormDB.near_or_create(
        parent=H.to_json(),
        spectral_parameter=spectral_parameter,
        bound_m=bound_m,
        set_coefficients=None,
        y=None,
        q=None,
    )
    f = load_object(maass_form_db)
    f = (H, [(1, -1), (1, 1)], 0.1, spectral_parameter, bound_m, None, True)
    result1 = list(check_coefficients_of_computed_object([f]))
    sp = result1[0][1]
    assert abs(sp[0] - 0.212132034355964) < 1e-10
    assert abs(sp[1] - 0.212132034355964) < 1e-10
    assert abs(sp[2]) < 1e-14
