r"""
Base class for individual Maass form elements.

Provides the abstract SageMath Element structure that domain-specific
elements (Hilbert, Kleinian) inherit from. Handles spectral parameter
storage, coefficient access, serialization, comparison, and arithmetic.

EXAMPLES::

    sage: from maass_form_core.spaces.base_element import MaassFormElement
    sage: from sage.structure.element import Element
    sage: issubclass(MaassFormElement, Element)
    True
"""

import logging
from copy import copy

from sage.structure.element import Element

log = logging.getLogger(__name__)


# Minimal test helpers for doctests — avoids importing the full space module
def _make_test_space_and_element():
    r"""
    Create minimal test space and element classes for doctests.

    OUTPUT:

    - tuple (TestSpace, TestElement) classes

    EXAMPLES::

        sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
        sage: TestSpace, TestElement = _make_test_space_and_element()
        sage: S = TestSpace('G')
        sage: e = TestElement(S, int(0))
        sage: e.spectral_parameter()
        0
    """
    from maass_form_core.spaces.base_space import MaassFormSpace

    class TestElement(MaassFormElement):
        def to_json(self):
            return {
                "parent": self.parent().to_json(),
                "spectral_parameter": str(self.spectral_parameter()),
                "coefficients": str(self.coefficients()),
            }

        @classmethod
        def from_json(cls, data, **kwargs):
            raise NotImplementedError

    class TestSpace(MaassFormSpace):
        Element = TestElement

        def to_json(self):
            return {"group": str(self.group()), "cuspidal": self.is_cuspidal()}

    return TestSpace, TestElement


class MaassFormElement(Element):
    r"""Abstract base for individual Maass form instances.

    This class provides the shared infrastructure for Maass form elements
    across different mathematical domains. It implements the SageMath ``Element``
    interface for the Parent/Element coercion framework.

    Subclasses must:

    - Override ``to_json()`` and ``from_json()`` for domain-specific serialization
    - Override ``__repr__()`` for domain-specific string representation

    INPUT:

    - ``parent`` -- a MaassFormSpace; the parent space
    - ``spectral_parameter`` -- the spectral parameter (complex number or tuple)
    - ``coefficients`` -- (optional) coefficient object

    EXAMPLES::

        sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
        sage: TestSpace, TestElement = _make_test_space_and_element()
        sage: S = TestSpace('G', cuspidal=False)
        sage: e = TestElement(S, 0.5)
        sage: e.spectral_parameter()
        0.500000000000000
        sage: e.is_cuspidal()
        False
        sage: e.coefficients() is None
        True
        sage: e.parent() is S
        True
    """

    def __init__(self, parent, spectral_parameter, coefficients=None, **kwargs):
        r"""
        Initialize a Maass form element.

        INPUT:

        - ``parent`` -- MaassFormSpace; the parent space
        - ``spectral_parameter`` -- spectral parameter value(s)
        - ``coefficients`` -- (optional) coefficient object
        - ``**kwargs`` -- additional keyword arguments passed to Element

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e = TestElement(S, CC(0.5, 14.1))
            sage: e.spectral_parameter()
            0.500000000000000 + 14.1000000000000*I
            sage: e.parent() is S
            True
            sage: e.is_cuspidal()
            True
            sage: e2 = TestElement(S, CC(0, 1), coefficients='some_coeffs')
            sage: e2.coefficients()
            'some_coeffs'
        """
        super().__init__(parent, **kwargs)
        self.cuspidal = parent.is_cuspidal()
        self._spectral_parameter = spectral_parameter
        self._coefficients = coefficients

    def spectral_parameter(self):
        r"""
        Return the spectral parameter of this Maass form.

        OUTPUT:

        - The spectral parameter (complex number or tuple of complex numbers)

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e = TestElement(S, CC(0.5, 14.1))
            sage: e.spectral_parameter()
            0.500000000000000 + 14.1000000000000*I
            sage: e2 = TestElement(S, (CC(0.5, 1), CC(0.5, 2)))
            sage: e2.spectral_parameter()
            (0.500000000000000 + 1.00000000000000*I, 0.500000000000000 + 2.00000000000000*I)
        """
        return self._spectral_parameter

    def coefficients(self):
        r"""
        Return the coefficients of this Maass form.

        OUTPUT:

        - The coefficient object, or None if not yet computed

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e = TestElement(S, 0.5)
            sage: e.coefficients() is None
            True
            sage: e2 = TestElement(S, 0.5, coefficients='my_coeffs')
            sage: e2.coefficients()
            'my_coeffs'
        """
        return self._coefficients

    def is_cuspidal(self):
        r"""
        Return whether this form is cuspidal.

        Delegates to the parent space's ``is_cuspidal()`` method.

        OUTPUT:

        - bool

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: e1 = TestElement(TestSpace('G', cuspidal=True), 0.5)
            sage: e1.is_cuspidal()
            True
            sage: e2 = TestElement(TestSpace('G', cuspidal=False), 0.5)
            sage: e2.is_cuspidal()
            False
        """
        return self.cuspidal

    def to_json(self) -> dict:
        r"""
        Return a JSON-serializable representation of this element.

        Subclasses must override this method to include domain-specific data
        (e.g., spectral parameter encoding, coefficient serialization).

        OUTPUT:

        - dict; JSON-compatible dictionary

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e = TestElement(S, 0.5)
            sage: d = e.to_json()
            sage: 'parent' in d and 'spectral_parameter' in d
            True
        """
        raise NotImplementedError("Subclasses must implement to_json()")

    @classmethod
    def from_json(cls, data, **kwargs):
        r"""
        Create an element from a JSON-compatible dictionary.

        Subclasses must override this method for domain-specific deserialization.

        INPUT:

        - ``data`` -- dict or str; JSON data
        - ``**kwargs`` -- additional construction options

        OUTPUT:

        - An instance of the element class

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import MaassFormElement
            sage: MaassFormElement.from_json({})
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclasses must implement from_json()
        """
        raise NotImplementedError("Subclasses must implement from_json()")

    def __reduce__(self):
        r"""
        Return pickling data for this element.

        OUTPUT:

        - tuple (class, constructor_args)

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e = TestElement(S, 0.5)
            sage: cls, args = e.__reduce__()
            sage: cls == TestElement
            True
            sage: args[1]
            0.500000000000000
        """
        return self.__class__, (
            self.parent(),
            self.spectral_parameter(),
            self._coefficients,
        )

    def __hash__(self):
        r"""
        Return a hash value for this element.

        Hash is based on parent, spectral parameter, and coefficients.

        OUTPUT:

        - integer hash value

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e1 = TestElement(S, 0.5)
            sage: e2 = TestElement(S, 0.5)
            sage: hash(e1) == hash(e2)
            True
            sage: e3 = TestElement(S, 1.5)
            sage: hash(e1) == hash(e3)
            False
        """
        coeffs_str = str(self._coefficients) if self._coefficients else ""
        return hash((str(self.parent().to_json()), str(self.spectral_parameter()), coeffs_str))

    def __eq__(self, other):
        r"""
        Check equality with another element.

        Two elements are equal if they have the same parent, spectral parameter,
        and coefficients.

        INPUT:

        - ``other`` -- object to compare with

        OUTPUT:

        - bool

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e1 = TestElement(S, 0.5)
            sage: e2 = TestElement(S, 0.5)
            sage: e1 == e2
            True
            sage: e3 = TestElement(S, 1.5)
            sage: e1 == e3
            False
            sage: e1 == 'not an element'
            False
        """
        if not isinstance(other, self.__class__):
            return False
        return (
            self.parent() == other.parent()
            and self.spectral_parameter() == other.spectral_parameter()
            and self.coefficients() == other.coefficients()
        )

    def __ne__(self, other):
        r"""
        Check inequality with another element.

        INPUT:

        - ``other`` -- object to compare with

        OUTPUT:

        - bool

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e1 = TestElement(S, 0.5)
            sage: e2 = TestElement(S, 0.5)
            sage: e1 != e2
            False
            sage: e3 = TestElement(S, 1.5)
            sage: e1 != e3
            True
        """
        return not self.__eq__(other)

    def __copy__(self):
        r"""
        Return a copy of this element.

        OUTPUT:

        - A new element with copied attributes

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: from copy import copy
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e = TestElement(S, 0.5, coefficients='c')
            sage: e2 = copy(e)
            sage: e == e2
            True
            sage: e is e2
            False
        """
        coefficients = copy(self._coefficients)
        return self.__class__(self.parent(), self._spectral_parameter, coefficients)

    def _add_(self, other):
        r"""
        Add two Maass form elements.

        Both elements must have the same parent. If either has no coefficients,
        the other's coefficients are used.

        INPUT:

        - ``other`` -- MaassFormElement with the same parent

        OUTPUT:

        - MaassFormElement

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e1 = TestElement(S, 0.5)
            sage: e2 = TestElement(S, 0.5)
            sage: result = e1._add_(e2)
            sage: result.parent() is S
            True
            sage: result.spectral_parameter()
            0.500000000000000
        """
        if not isinstance(other, self.__class__):
            raise TypeError(f"Cannot add {type(self).__name__} and {type(other).__name__}")
        if other.parent() != self.parent():
            raise ValueError("Cannot add elements from different spaces")
        result = copy(self)
        c_self = self._coefficients
        c_other = other._coefficients
        if not c_self:
            result._coefficients = c_other
        elif not c_other:
            result._coefficients = c_self
        else:
            result._coefficients = c_self + c_other
        return result

    def __add__(self, other):
        r"""
        Add two Maass form elements.

        INPUT:

        - ``other`` -- MaassFormElement

        OUTPUT:

        - MaassFormElement

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e1 = TestElement(S, 0.5)
            sage: e2 = TestElement(S, 0.5)
            sage: result = e1 + e2
            sage: result.parent() is S
            True
        """
        return self._add_(other)

    def _sub_(self, other):
        r"""
        Subtract two Maass form elements.

        INPUT:

        - ``other`` -- MaassFormElement with the same parent

        OUTPUT:

        - MaassFormElement

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e1 = TestElement(S, 0.5)
            sage: e2 = TestElement(S, 0.5)
            sage: result = e1._sub_(e2)
            sage: result.parent() is S
            True
        """
        return self + other * -1

    def __sub__(self, other):
        r"""
        Subtract two Maass form elements.

        INPUT:

        - ``other`` -- MaassFormElement

        OUTPUT:

        - MaassFormElement

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e1 = TestElement(S, 0.5)
            sage: e2 = TestElement(S, 0.5)
            sage: result = e1 - e2
            sage: result.parent() is S
            True
        """
        return self._sub_(other)

    def __mul__(self, other):
        r"""
        Multiply this element by a scalar.

        Scales the coefficients by the given scalar value.

        INPUT:

        - ``other`` -- scalar (int, float, complex, or SageMath number)

        OUTPUT:

        - MaassFormElement with scaled coefficients

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e = TestElement(S, 0.5, coefficients=[1, 2, 3])
            sage: result = e * 2
            sage: result.coefficients()
            [1, 2, 3, 1, 2, 3]

        When coefficients are None, multiplication still returns an element::

            sage: e2 = TestElement(S, 0.5)
            sage: result2 = e2 * 2
            sage: result2.coefficients() is None
            True
        """
        result = copy(self)
        if result._coefficients is not None:
            result._coefficients = other * result._coefficients
        return result

    def _rmul_(self, other):
        r"""
        Right-multiply this element by a scalar (SageMath coercion framework).

        INPUT:

        - ``other`` -- scalar

        OUTPUT:

        - MaassFormElement with scaled coefficients

        EXAMPLES::

            sage: from maass_form_core.spaces.base_element import _make_test_space_and_element
            sage: TestSpace, TestElement = _make_test_space_and_element()
            sage: S = TestSpace('G')
            sage: e = TestElement(S, 0.5)
            sage: result = e._rmul_(2)
            sage: result.parent() is S
            True
        """
        return self.__mul__(other)
