"""
Database layer for Kleinian Maass forms.

Requires the optional ``[db]`` extras::

    pip install 'maass_forms_klein[db]'
"""

try:
    import comp_manager  # noqa: F401
    import mongoengine  # noqa: F401
except ImportError as _e:
    raise ModuleNotFoundError(
        "maass_forms_klein.database requires the optional [db] extras. "
        "Install with: pip install 'maass_forms_klein[db]'"
    ) from _e
