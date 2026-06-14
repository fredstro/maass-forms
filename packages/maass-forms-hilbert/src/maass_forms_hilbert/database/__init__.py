"""
Database layer for Hilbert Maass forms.

Requires the optional ``[db]`` extras::

    pip install 'maass_forms_hilbert[db]'
"""

try:
    import comp_manager  # noqa: F401
    import mongoengine  # noqa: F401
except ImportError as _e:
    raise ModuleNotFoundError(
        "maass_forms_hilbert.database requires the optional [db] extras. "
        "Install with: pip install 'maass_forms_hilbert[db]'"
    ) from _e
