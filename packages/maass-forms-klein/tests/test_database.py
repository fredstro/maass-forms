"""
Base test suite using pytest for the maass_forms_klein database models.
"""

from mongoengine import get_connection

from . import dbconnect

dbconnect()
conn = get_connection("default")


def test_database_connection():
    """A mock database connection is established via mongomock."""
    assert conn is not None


def test_insert_kleinian_maass_form_db():
    """KleinianMaassFormDB collection is empty on a fresh mock database."""
    from maass_forms_klein.database.models import KleinianMaassFormDB

    assert KleinianMaassFormDB.objects.count() == 0


def test_insert_kleinian_group_db():
    """KleinianGroupDB collection is empty on a fresh mock database."""
    from maass_forms_klein.database.models import KleinianGroupDB

    assert KleinianGroupDB.objects.count() == 0


def test_kleinian_group_db_save_and_retrieve():
    """A KleinianGroupDB document can be saved and retrieved."""
    from maass_forms_klein.database.models import KleinianGroupDB

    group = KleinianGroupDB(_name_string="4_1")
    group.save()
    try:
        assert KleinianGroupDB.objects(name="4_1").count() == 1
        retrieved = KleinianGroupDB.objects(name="4_1").first()
        assert retrieved._name_string == "4_1"
        assert repr(retrieved) == "4_1"
    finally:
        group.delete()
