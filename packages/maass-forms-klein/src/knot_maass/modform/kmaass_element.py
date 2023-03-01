"""
Elements of spaces of Maass waveforms for Knot complements.

"""
from typing import ParamSpec

from sage.structure.parent import Parent
from sage.structure.element import Element, Matrix, Vector

P = ParamSpec('P')


class KnotMaassFormElement(Element):
    r"""
    Element of a KnotMaassFormSpace.

    TODO: Decide whether to keep this as a Python class or make a cdef Cython class.
    """

    coefficients = []

    def __init__(self,  parent: Parent, *args: P.args, **kwargs: P.kwargs) -> None:
        r"""

        INPUT:

        - `parent` -- parent space of type KnotMaassFormSpace
        - `args`   -- arguments
        - `kwargs` -- keyword arguments

        EXAMPLES::

        sage: from knot_maass.all import KnotMaassFormSpace, KnotMaassFormElement
        sage: S = KnotMaassFormSpace()
        sage: F = S()
        sage: sage: TestSuite(F).run()

        """
        super(KnotMaassFormElement, self).__init__(parent, *args, **kwargs)

    def compute_coefficients(self, *args: P.args, **kwargs: P.kwargs) -> None:
        """
        Compute the Fourier coefficients of self.

        EXAMPLES::

            sage: from knot_maass.all import KnotMaassFormSpace, KnotMaassFormElement
            sage: S = KnotMaassFormSpace()
            sage: F = S.an_element()

        """
        matrixV = self.setup_matrix(*args, **kwargs)
        matrixV = self.normalise_matrix(matrixV, *args, **kwargs)
        coefficients = self.solve_system(matrixV, *args, **kwargs)
        self.coefficients = coefficients

    def setup_matrix(self, *args: P.args, **kwargs: P.kwargs) -> Matrix:
        r"""
        Set up the matrix for the linear system.

        EXAMPLES::

            sage: from knot_maass.all import KnotMaassFormSpace, KnotMaassFormElement
            sage: S = KnotMaassFormSpace()
            sage: F = S.an_element()
            sage: V = F.setup_matrix()
        """
        raise NotImplementedError

    def normalise_matrix(self, *args: P.args, **kwargs: P.kwargs) -> tuple:
        r"""
        Normalise the matrix, e.g. set first coefficient to 1 etc.

        EXAMPLES::

            sage: from knot_maass.all import KnotMaassFormSpace, KnotMaassFormElement
            sage: S = KnotMaassFormSpace()
            sage: F = S.an_element()
            sage: V = F.setup_matrix()
            sage: V = F.normalise_matrix(V)
        """
        raise NotImplementedError

    def solve_system(self, *args: P.args, **kwargs: P.kwargs) -> Vector:
        r"""
        Solve the system.

        EXAMPLES::

            sage: from knot_maass.all import KnotMaassFormSpace, KnotMaassFormElement
            sage: S = KnotMaassFormSpace()
            sage: F = S.an_element()
            sage: V = F.setup_matrix()
            sage: V, B = F.normalise_matrix(V)
            sage: X = F.solve_system(V, B)
        """
        raise NotImplementedError