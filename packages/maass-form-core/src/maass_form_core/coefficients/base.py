r"""
Base class for Maass form coefficient containers.

Provides shared storage, accessors, hashing, and comparison for coefficient
containers used by both Hilbert and Kleinian Maass form packages.

This class does not inherit from any SageMath base class, allowing domain
packages to mix it with their preferred SageMath parent (``ModuleElement``
for Hilbert, ``Parent`` for Klein) without MRO conflicts.

EXAMPLES::

    sage: from maass_form_core.coefficients.base import CoefficientManagerBase
    sage: hasattr(CoefficientManagerBase, 'M')
    True
    sage: hasattr(CoefficientManagerBase, 'spectral_parameter')
    True
"""

import logging

from sage.matrix.constructor import matrix as sage_matrix
from sage.structure.element import Matrix

log = logging.getLogger(__name__)


class CoefficientManagerBase:
    r"""Abstract base for coefficient storage across all Maass form domains.

    Provides shared infrastructure for coefficient management:

    - Storage of coefficient matrix, bounds (M), spectral parameter, space, Y, Q
    - Accessors for all stored data
    - ``__hash__()`` and ``__eq__()`` based on stored data
    - ``__repr__()`` skeleton (subclasses customize the label)

    Subclasses must:

    - Call ``_init_coefficients()`` from their ``__init__`` to set up storage
    - Implement ``to_json()`` and ``from_json()`` for domain-specific serialization
    - Implement ``__getitem__()`` for domain-specific indexing
    - Implement ``keys()`` for domain-specific key listing
    - Override ``_coefficient_label`` for repr customization

    INPUT:

    - ``coefficients`` -- matrix; the coefficient matrix
    - ``M`` -- bounds for coefficient indices (tuple of tuples or integer)
    - ``spectral_parameter`` -- spectral parameter value(s)
    - ``space`` -- the parent Maass form space
    - ``Y`` -- (optional) Y parameter(s) used in computation
    - ``Q`` -- (optional) Q parameter(s) used in computation
    - ``set_coefficients`` -- (optional) dict of pre-set coefficient values

    EXAMPLES::

        sage: from maass_form_core.coefficients.base import CoefficientManagerBase
        sage: from sage.matrix.constructor import matrix

        sage: class TestCoeffs(CoefficientManagerBase):
        ....:     _coefficient_label = "Test"
        ....:     def __init__(self, coefficients, M, s, space, **kw):
        ....:         self._init_coefficients(coefficients, M, s, space, **kw)
        ....:     def to_json(self):
        ....:         return {}
        ....:     def keys(self):
        ....:         return list(range(self._coefficients.nrows()))
        ....:     def __getitem__(self, key):
        ....:         return self._coefficients[key, 0]
        sage: C = TestCoeffs(matrix([[1],[2],[3]]), ((-1,1),), (0.5,), 'space')
        sage: C.M()
        ((-1, 1),)
        sage: C.spectral_parameter()
        (0.500000000000000,)
        sage: C.coefficient_matrix()
        [1]
        [2]
        [3]
    """

    # Subclasses should set this for __repr__
    _coefficient_label = "Maass form"

    def _init_coefficients(
        self,
        coefficients,
        M,
        spectral_parameter,
        space,
        Y=None,
        Q=None,
        set_coefficients=None,
        **kwargs,
    ):
        r"""
        Initialize shared coefficient storage.

        This method should be called from the subclass ``__init__`` after
        calling the SageMath parent's ``__init__``.

        INPUT:

        - ``coefficients`` -- matrix; the coefficient data
        - ``M`` -- bounds for coefficient indices
        - ``spectral_parameter`` -- spectral parameter
        - ``space`` -- parent Maass form space
        - ``Y`` -- (optional) Y parameters
        - ``Q`` -- (optional) Q parameters
        - ``set_coefficients`` -- (optional) dict of pre-set coefficients

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     _coefficient_label = "Test"
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self):
            ....:         return {}
            ....:     def keys(self):
            ....:         return []
            ....:     def __getitem__(self, key):
            ....:         return self._coefficients[key, 0]
            sage: C = TestCoeffs(matrix([[1]]), ((0,0),), CC(0.5,1), 'sp', Y=(0.5,), Q=(3,))
            sage: C.Y()
            (0.500000000000000,)
            sage: C.Q()
            (3,)
            sage: C.set_coefficients() is None
            True
        """
        if isinstance(coefficients, Matrix):
            self._coefficients = coefficients
        else:
            self._coefficients = sage_matrix(coefficients)
        self._M = M
        self._spectral_parameter = spectral_parameter
        self._space = space
        self._Y = Y or ()
        self._Q = Q
        self._set_coefficients = set_coefficients

    def M(self):
        r"""
        Return the M parameter (bounds for coefficient indices).

        OUTPUT:

        - The bounds parameter (tuple of tuples or integer, depending on domain)

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self): return {}
            ....:     def keys(self): return []
            ....:     def __getitem__(self, k): pass
            sage: C = TestCoeffs(matrix([[1],[2],[3]]), ((-1,1),), (0.5,), 'sp')
            sage: C.M()
            ((-1, 1),)
        """
        return self._M

    def Y(self):
        r"""
        Return the Y parameter used in computation.

        OUTPUT:

        - tuple of real values, or empty tuple if not set

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self): return {}
            ....:     def keys(self): return []
            ....:     def __getitem__(self, k): pass
            sage: TestCoeffs(matrix([[1]]), 1, 0, 'sp').Y()
            ()
            sage: TestCoeffs(matrix([[1]]), 1, 0, 'sp', Y=(0.5, 0.7)).Y()
            (0.500000000000000, 0.700000000000000)
        """
        return self._Y

    def Q(self):
        r"""
        Return the Q parameter.

        OUTPUT:

        - tuple of integers, or None if not set

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self): return {}
            ....:     def keys(self): return []
            ....:     def __getitem__(self, k): pass
            sage: TestCoeffs(matrix([[1]]), 1, 0, 'sp').Q() is None
            True
            sage: TestCoeffs(matrix([[1]]), 1, 0, 'sp', Q=(3, 3)).Q()
            (3, 3)
        """
        return self._Q

    def space(self):
        r"""
        Return the Maass form space associated with these coefficients.

        OUTPUT:

        - The parent space object

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self): return {}
            ....:     def keys(self): return []
            ....:     def __getitem__(self, k): pass
            sage: TestCoeffs(matrix([[1]]), 1, 0, 'my_space').space()
            'my_space'
        """
        return self._space

    def spectral_parameter(self):
        r"""
        Return the spectral parameter.

        OUTPUT:

        - The spectral parameter (tuple or single value, depending on domain)

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self): return {}
            ....:     def keys(self): return []
            ....:     def __getitem__(self, k): pass
            sage: TestCoeffs(matrix([[1]]), 1, (CC(0.5,1),), 'sp').spectral_parameter()
            (0.500000000000000 + 1.00000000000000*I,)
        """
        return self._spectral_parameter

    def coefficient_matrix(self):
        r"""
        Return the underlying coefficient matrix.

        OUTPUT:

        - A SageMath matrix

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self): return {}
            ....:     def keys(self): return []
            ....:     def __getitem__(self, k): pass
            sage: m = matrix(RR, [[1],[2],[3]])
            sage: TestCoeffs(m, 1, 0, 'sp').coefficient_matrix()
            [1.00000000000000]
            [2.00000000000000]
            [3.00000000000000]
        """
        return self._coefficients

    def set_coefficients(self):
        r"""
        Return the dictionary of pre-set coefficient values.

        These are the coefficients that were fixed during the computation
        (e.g., normalization conditions).

        OUTPUT:

        - dict or None

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self): return {}
            ....:     def keys(self): return []
            ....:     def __getitem__(self, k): pass
            sage: TestCoeffs(matrix([[1]]), 1, 0, 'sp').set_coefficients() is None
            True
            sage: d = {(0,0): 1}
            sage: TestCoeffs(matrix([[1]]), 1, 0, 'sp', set_coefficients=d).set_coefficients()
            {(0, 0): 1}
        """
        return self._set_coefficients

    def __hash__(self):
        r"""
        Return a hash value based on the coefficient data.

        OUTPUT:

        - integer hash value

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self): return {}
            ....:     def keys(self): return []
            ....:     def __getitem__(self, k): pass
            sage: m = matrix(RR, [[1],[2],[3]])
            sage: C1 = TestCoeffs(m, ((-1,1),), (0.5,), 'sp')
            sage: C2 = TestCoeffs(m, ((-1,1),), (0.5,), 'sp')
            sage: hash(C1) == hash(C2)
            True
        """
        self._coefficients.set_immutable()
        return hash(
            (
                self._coefficients,
                str(self._Y),
                str(self._Q),
                str(self._M),
                str(self._spectral_parameter),
                str(self._space),
                str(self._set_coefficients),
            )
        )

    def __eq__(self, other):
        r"""
        Check equality with another coefficient object.

        Two coefficient objects are equal if they have the same space,
        coefficient matrix, spectral parameter, M, Y, Q, and set_coefficients.

        INPUT:

        - ``other`` -- object to compare with

        OUTPUT:

        - bool

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self): return {}
            ....:     def keys(self): return []
            ....:     def __getitem__(self, k): pass
            sage: m = matrix(RR, [[1],[2],[3]])
            sage: C1 = TestCoeffs(m, ((-1,1),), (0.5,), 'sp')
            sage: C2 = TestCoeffs(m, ((-1,1),), (0.5,), 'sp')
            sage: C1 == C2
            True
            sage: C3 = TestCoeffs(m, ((-1,1),), (1.5,), 'sp')
            sage: C1 == C3
            False
            sage: C1 == 'not coefficients'
            False
        """
        if not isinstance(other, self.__class__):
            return False
        return (
            self._space == other._space
            and self._coefficients == other._coefficients
            and self._spectral_parameter == other._spectral_parameter
            and self._M == other._M
            and self._Y == other._Y
            and self._Q == other._Q
            and self._set_coefficients == other._set_coefficients
        )

    def __repr__(self):
        r"""
        Return string representation of the coefficients.

        Uses ``_coefficient_label`` for domain-specific naming.

        OUTPUT:

        - string

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     _coefficient_label = "Test"
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self): return {}
            ....:     def keys(self): return []
            ....:     def __getitem__(self, k): pass
            sage: TestCoeffs(matrix([[1],[2],[3]]), ((-1,1),), (0.5,), 'sp')
            Coefficients of a Test with M=((-1, 1),)
        """
        return f"Coefficients of a {self._coefficient_label} with M={self._M}"

    def __iter__(self):
        r"""
        Iterate over coefficient values (first column of the matrix).

        OUTPUT:

        - iterator over coefficient values

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self): return {}
            ....:     def keys(self): return []
            ....:     def __getitem__(self, k): pass
            sage: C = TestCoeffs(matrix(RR, [[1],[2],[3]]), 1, 0, 'sp')
            sage: list(C)
            [1.00000000000000, 2.00000000000000, 3.00000000000000]
        """
        return iter(self._coefficients.column(0))

    def to_json(self) -> dict:
        r"""
        Return a JSON-serializable representation of the coefficients.

        Subclasses must override for domain-specific serialization.

        OUTPUT:

        - dict

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: CoefficientManagerBase.to_json(None)
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclasses must implement to_json()
        """
        raise NotImplementedError("Subclasses must implement to_json()")

    def keys(self):
        r"""
        Return the list of coefficient keys/indices.

        Subclasses must override for domain-specific index enumeration.

        OUTPUT:

        - list

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: CoefficientManagerBase.keys(None)
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclasses must implement keys()
        """
        raise NotImplementedError("Subclasses must implement keys()")

    def __getitem__(self, key):
        r"""
        Get the coefficient at the given index.

        Subclasses must override for domain-specific indexing.

        INPUT:

        - ``key`` -- index (type depends on domain)

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: CoefficientManagerBase.__getitem__(None, 0)
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclasses must implement __getitem__()
        """
        raise NotImplementedError("Subclasses must implement __getitem__()")

    def __copy__(self):
        r"""
        Return a shallow copy of this coefficient object.

        Creates a new instance with the same data. The coefficient matrix
        is copied; other fields are shared references.

        OUTPUT:

        - A new instance of the same class

        EXAMPLES::

            sage: from maass_form_core.coefficients.base import CoefficientManagerBase
            sage: from sage.matrix.constructor import matrix
            sage: from copy import copy
            sage: class TestCoeffs(CoefficientManagerBase):
            ....:     def __init__(self, c, M, s, space, **kw):
            ....:         self._init_coefficients(c, M, s, space, **kw)
            ....:     def to_json(self): return {}
            ....:     def keys(self): return []
            ....:     def __getitem__(self, k): pass
            sage: m = matrix(RR, [[1],[2],[3]])
            sage: C = TestCoeffs(m, ((-1,1),), (0.5,), 'sp')
            sage: C2 = copy(C)
            sage: C == C2
            True
            sage: C is C2
            False
        """
        from copy import copy as _copy

        result = object.__new__(self.__class__)
        result._coefficients = _copy(self._coefficients)
        result._M = self._M
        result._spectral_parameter = self._spectral_parameter
        result._space = self._space
        result._Y = self._Y
        result._Q = self._Q
        result._set_coefficients = self._set_coefficients
        return result
