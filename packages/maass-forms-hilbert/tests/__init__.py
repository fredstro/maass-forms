import logging
from unittest.mock import patch

import mongomock
from comp_manager import create_app
from mongoengine import ConnectionFailure, get_connection, get_db, register_connection


def patched_create_coll(self, name: str, **kwargs) -> object:
    """
    Patched create_collection method.
    """
    self._store.create_collection(name)
    return self[name]


def dbconnect() -> None:
    """
    Connect to the (mock) database and patch the create_collection method
    since it is not implemented for capped collections.
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
