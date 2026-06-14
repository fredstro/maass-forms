"""
Base QuerySet class with SageMath type coercion.

Provides a QuerySet base that handles Sage Integer indexing, which is
common across all Maass form domain packages.

EXAMPLES::

    sage: from sage.rings.integer import Integer
    sage: isinstance(Integer(3), int)
    False
"""

from comp_manager.core.queryset import QuerySetCompat
from sage.rings.integer import Integer


class MaassFormQuerySet(QuerySetCompat):
    """Base QuerySet with SageMath Integer coercion for indexing.

    Both HilbertMaassformQuerySet and KleinianMaassFormQuerySet share
    the same pattern of converting Sage Integers to Python ints for
    MongoDB compatibility.

    EXAMPLES::

        sage: import warnings
        sage: warnings.filterwarnings('ignore', category=DeprecationWarning)
        sage: warnings.filterwarnings('ignore', category=PendingDeprecationWarning)
        sage: from maass_form_core.database.queryset import MaassFormQuerySet
        sage: MaassFormQuerySet  # class exists
        <class 'maass_form_core.database.queryset.MaassFormQuerySet'>
    """

    def __getitem__(self, item):
        r"""
        Retrieve an item or slice, converting Sage Integers to Python ints.

        INPUT:

        - ``item`` -- an integer index (Python int or Sage Integer) or a slice

        OUTPUT:

        - The item at the specified index or a sliced QuerySet

        EXAMPLES::

            sage: from sage.rings.integer import Integer
            sage: # Integer coercion is handled transparently
            sage: Integer(3) == int(3)
            True
        """
        if isinstance(item, Integer):
            item = int(item)
        if isinstance(item, slice) and isinstance(item.stop, Integer):
            item = slice(int(item.start), int(item.stop))
        return super().__getitem__(item)
