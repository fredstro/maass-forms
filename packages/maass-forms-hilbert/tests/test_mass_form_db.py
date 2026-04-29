"""
Base test suite using pytest
"""


def test_insert_mass_form_db():
    from maass_forms_hilbert.database.models import HilbertMaassFormDB

    assert HilbertMaassFormDB.objects.count() == 0
