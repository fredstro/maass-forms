"""
Shared test fixtures and assertions for Maass forms packages.
"""

import logging

import mongomock
from mongoengine import get_db, register_connection


def patched_create_coll(self, name: str, **kwargs) -> object:
    """
    Patched create_collection method for mongomock compatibility.

    This function patches the create_collection method to work with mongomock
    since it doesn't fully implement capped collections.

    INPUT:

    - ``name`` -- string, the name of the collection to create
    - ``kwargs`` -- additional keyword arguments (ignored)

    OUTPUT:

    The created collection object.

    EXAMPLES::

        sage: from maass_form_core.testing import patched_create_coll
        sage: import mongomock
        sage: client = mongomock.MongoClient()
        sage: db = client.test_db
        sage: coll = patched_create_coll(db, "test_collection")
        sage: coll.name
        'test_collection'
    """
    self._store.create_collection(name)
    return self[name]


def connect_mockdb() -> None:
    """
    Connect to the mock database and patch the create_collection method.

    This function sets up a mongomock connection for testing purposes and
    applies a patch to handle capped collections which are not fully
    implemented in mongomock.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()  # Sets up mock database connection
        sage: from mongoengine import get_db
        sage: db = get_db("default")
        sage: hasattr(db, 'create_collection')
        True
    """
    register_connection(
        alias="default",
        host="localhost",
        port=27017,
        uuidRepresentation="standard",
        mongo_client_class=mongomock.MongoClient,
    )
    db = get_db("default")
    logging.debug(f"{db}: {id(db)} id(db['mongo_cache_d_b'])={id(db['mongo_cache_d_b'])}")
    db.create_collection = lambda name, **kwargs: patched_create_coll(db, name, **kwargs)
