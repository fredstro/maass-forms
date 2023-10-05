"""
Space of Maass waveforms for Knot complements.

"""
from typing import ParamSpec

from sage.structure.parent import Parent

from .kmaass_element import KnotMaassFormElement

P = ParamSpec('P')


class KnotMaassFormSpace(Parent):
    """
    Class for spaces of Maass waveform for Knot complements.
    """

    knot = None
    group = None
    Element = KnotMaassFormElement

    def __init__(self,  *args: P.args, **kwargs: P.kwargs) -> None:
        r"""

        INPUT:

        - `args`   -- arguments
        - `kwargs` -- keyword arguments

        EXAMPLES::

            sage: from knot_maass.all import KnotMaassFormSpace
            sage: S = KnotMaassFormSpace()
            sage: sage: TestSuite(S).run()

        """
        super(KnotMaassFormSpace, self).__init__(*args, **kwargs)

    def _element_constructor_(self, *args: P.args, **kwargs: P.kwargs) -> KnotMaassFormElement:
        r"""
        Construct an element of this space.

        EXAMPLES::

            sage: from knot_maass.all import KnotMaassFormSpace
        """
        return self.element_class(self, *args, **kwargs)


