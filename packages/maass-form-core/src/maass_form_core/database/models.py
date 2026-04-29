"""
Base database document classes for Maass form packages.

Provides shared embedded documents, QuerySet patterns, and abstract
base models used by both maass-forms-hilbert and maass-forms-klein.

EXAMPLES::

    sage: import mongoengine as me
    sage: from maass_form_core.database.models import PointDB
    sage: pt = PointDB(x=0.5, y=14.1)
    sage: str(pt)
    '(0.5, 14.1)'
"""

import logging

import mongoengine as me

log = logging.getLogger(__name__)


class PointDB(me.EmbeddedDocument):
    """Embedded document representing a point in the complex plane.

    Stores real and imaginary parts for spectral parameters and other
    complex-valued mathematical objects in the database.

    EXAMPLES::

        sage: from maass_form_core.database.models import PointDB
        sage: pt = PointDB(x=0.5, y=14.1)
        sage: str(pt)
        '(0.5, 14.1)'
        sage: pt = PointDB(x=-1.0, y=0.0)
        sage: str(pt)
        '(-1.0, 0.0)'
    """

    x = me.FloatField(required=True, help_text="Real part of the complex number")
    y = me.FloatField(required=True, help_text="Imaginary part of the complex number")

    def __str__(self) -> str:
        """String representation of the point.

        OUTPUT:

        - String in the format "(x, y)"

        EXAMPLES::

            sage: from maass_form_core.database.models import PointDB
            sage: str(PointDB(x=1.5, y=-2.3))
            '(1.5, -2.3)'
        """
        return f"({self.x}, {self.y})"
