import logging

import mongomock
from comp_manager.utils import insert_object
from mongoengine import get_db, register_connection
from sage.matrix.constructor import matrix
from sage.misc.mrange import cartesian_product_iterator
from sage.rings.cc import CC

from maass_forms_hilbert.modform.coefficients import HilbertMaassCoefficients

fixtures_inserted = False


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

        sage: from maass_forms_hilbert.database.tests import patched_create_coll
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

        sage: from maass_forms_hilbert.database.tests import connect_mockdb
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


def insert_fixtures(reset_db: bool = False):
    """
    Insert test fixtures into the database for testing purposes.

    This function creates sample Hilbert Maass forms with various parameters
    and coefficients to be used in tests. It only inserts fixtures once
    unless reset_db is True.

    INPUT:

    - ``reset_db`` -- boolean (default: False), whether to reset the database
      before inserting fixtures

    EXAMPLES::

        sage: from maass_forms_hilbert.database.tests import connect_mockdb, insert_fixtures
        sage: connect_mockdb()
        sage: insert_fixtures(reset_db=True)
        sage: from maass_forms_hilbert.all import HilbertMaassFormDB
        sage: len(HilbertMaassFormDB.objects) > 0
        True
    """
    from maass_forms_hilbert.all import HilbertMaassForm, HilbertMaassFormDB, HilbertMaassFormSpace

    global fixtures_inserted
    if fixtures_inserted:
        return
    setp1 = {(1, 1): 1, (-1, 1): 2}
    setp2 = {(1, 1): -1, (-1, 1): 2}

    parameter_list = [
        {"space": 5, "s": (CC(0.5, 1), CC(0.5, 1)), "y": (0.3, 0.3), "set_parameters": setp1},
        {"space": 5, "s": (CC(0.5, 1), CC(0.6, 1)), "y": (0.2, 0.2), "set_parameters": setp1},
        {"space": 2, "s": (CC(0.5, 1), CC(0.5, 1)), "y": (0.2, 0.2), "set_parameters": setp2},
    ]
    m_list = [1, 4, 5]
    Cmat = {}
    for m in m_list:
        Cmat[m] = matrix(
            CC,
            [
                [CC(a, b)]
                for a, b in cartesian_product_iterator([range(-m, m + 1), range(-m, m + 1)])
            ],
        )
    if reset_db:
        HilbertMaassFormDB.objects.delete()
    for parameters in parameter_list:
        space = HilbertMaassFormSpace(parameters["space"])
        maass_form = HilbertMaassForm(space, parameters["s"])
        for m in m_list:
            C = HilbertMaassCoefficients(
                Cmat[m],
                ((-m, m), (-m, m)),
                parameters["s"],
                space,
                Y=parameters["y"],
                Q=(m + 1, m + 1),
                set_coefficients=parameters["set_parameters"],
            )
            maass_form._coefficients = C
            insert_object(maass_form)
        insert_object(maass_form)
    fixtures_inserted = True


def return_db_parameters():
    """
    Return sample database parameters for testing.

    This function returns a dictionary containing sample parameters that
    match the structure expected by the database models, useful for
    testing serialization and deserialization.

    OUTPUT:

    A dictionary containing sample coefficients, spectral parameters,
    and space information.

    EXAMPLES::

        sage: from maass_forms_hilbert.database.tests import return_db_parameters
        sage: params = return_db_parameters()
        sage: 'coefficients' in params
        True
        sage: 'parent' in params
        True
        sage: params['coefficients']['prec']
        53
    """
    return {
        "coefficients": {
            "M": [[-1, 1], [-1, 1]],
            "coefficients": [
                ["-1.00000000000000 - 1.00000000000000*I"],
                ["-1.00000000000000"],
                ["-1.00000000000000 + 1.00000000000000*I"],
                ["-1.00000000000000*I"],
                ["0.000000000000000"],
                ["1.00000000000000*I"],
                ["1.00000000000000 - 1.00000000000000*I"],
                ["1.00000000000000"],
                ["1.00000000000000 + 1.00000000000000*I"],
            ],
            "prec": 53,
            "spectral_parameter": [
                {"prec": 53, "val": "0.500000000000000 + 1.00000000000000*I"},
                {"prec": 53, "val": "0.500000000000000 + 1.00000000000000*I"},
            ],
            "set_coefficients": {},
            "space": {"number_field": {"polynomial": "x^2 - 5", "names": ["a"]}, "cuspidal": False},
            "index_tuples": None,
            "Y": [0.3, 0.3],
            "Q": [2, 2],
        },
        "parent": {"number_field": {"polynomial": "x^2 - 5", "names": ["a"]}, "cuspidal": False},
    }
