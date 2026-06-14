"""
Common database patterns and abstractions for Maass forms.

Provides shared MongoDB/MongoEngine infrastructure:
- PointDB embedded document for complex number storage
- MaassFormQuerySet with SageMath Integer coercion

Requires the optional ``[db]`` extras::

    pip install 'maass_form_core[db]'
"""

import warnings

try:
    with warnings.catch_warnings():
        # connexion/flask_mongoengine still call jsonschema APIs that
        # were deprecated in jsonschema 4.18. The warnings are not
        # actionable from this package and would otherwise break
        # doctests that import database models.
        warnings.simplefilter("ignore", DeprecationWarning)
        import comp_manager  # noqa: F401
        import mongoengine  # noqa: F401
except ImportError as e:
    raise ModuleNotFoundError(
        "maass_form_core.database requires the optional [db] extras. "
        "Install with: pip install 'maass_form_core[db]'"
    ) from e

from maass_form_core.database.models import PointDB

__all__ = [
    "PointDB",
    "MaassFormQuerySet",
]


def __getattr__(name):
    if name == "MaassFormQuerySet":
        from maass_form_core.database.queryset import MaassFormQuerySet

        return MaassFormQuerySet
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
