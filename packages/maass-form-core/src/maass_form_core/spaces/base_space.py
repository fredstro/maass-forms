r"""
Base class for mathematical spaces of Maass forms.

Provides the abstract SageMath Parent structure that domain-specific
spaces (Hilbert, Kleinian) inherit from. Handles group storage,
cuspidality, serialization, and the Parent/Element coercion framework.

EXAMPLES::

    sage: from maass_form_core.spaces.base_space import MaassFormSpace
    sage: from maass_form_core.spaces.base_element import MaassFormElement
    sage: from sage.structure.parent import Parent
    sage: issubclass(MaassFormSpace, Parent)
    True
"""

import logging

from sage.structure.parent import Parent

log = logging.getLogger(__name__)


class MaassFormSpace(Parent):
    r"""Abstract base space for Maass forms.

    This class provides the shared infrastructure for spaces of Maass forms
    across different mathematical domains. It implements the SageMath ``Parent``
    interface, enabling the coercion framework.

    Subclasses must:

    - Set the ``Element`` class attribute to their element class
    - Override ``to_json()`` and ``from_json()`` for domain-specific serialization
    - Override ``__repr__()`` for domain-specific string representation

    INPUT:

    - ``group`` -- a group object (e.g., HilbertModularGroup, KleinianGroup)
    - ``cuspidal`` -- bool (default: True); whether forms are required to be cuspidal

    EXAMPLES::

        sage: from maass_form_core.spaces.base_space import MaassFormSpace
        sage: from sage.structure.element import Element
        sage: class TestElement(Element):
        ....:     def __init__(self, parent, s):
        ....:         super().__init__(parent)
        ....:         self._s = s
        sage: class TestSpace(MaassFormSpace):
        ....:     Element = TestElement
        ....:     def to_json(self):
        ....:         return {'group': str(self.group()), 'cuspidal': self.is_cuspidal()}
        ....:     def __repr__(self):
        ....:         return f"TestSpace({self.group()})"
        sage: S = TestSpace('my_group', cuspidal=False)
        sage: S.group()
        'my_group'
        sage: S.is_cuspidal()
        False
        sage: hash(S) == hash(S)
        True
    """

    Element = None  # Must be set by subclass

    def __init__(self, group, **kwargs):
        r"""
        Initialize a Maass form space.

        INPUT:

        - ``group`` -- a group object
        - ``cuspidal`` -- bool (default: True); restrict to cuspidal forms

        EXAMPLES::

            sage: from maass_form_core.spaces.base_space import MaassFormSpace
            sage: from sage.structure.element import Element
            sage: class TestElement(Element):
            ....:     def __init__(self, parent, s):
            ....:         super().__init__(parent)
            sage: class TestSpace(MaassFormSpace):
            ....:     Element = TestElement
            ....:     def to_json(self):
            ....:         return {'group': 'G', 'cuspidal': self.is_cuspidal()}
            sage: S = TestSpace('G')
            sage: S.group()
            'G'
            sage: S.is_cuspidal()
            True
            sage: S2 = TestSpace('G', cuspidal=False)
            sage: S2.is_cuspidal()
            False
        """
        self._group = group
        self._cuspidal = kwargs.pop("cuspidal", True)
        super().__init__(**kwargs)

    def group(self):
        r"""
        Return the group of this space.

        OUTPUT:

        - The group object associated with this space

        EXAMPLES::

            sage: from maass_form_core.spaces.base_space import MaassFormSpace
            sage: from sage.structure.element import Element
            sage: class TestElement(Element):
            ....:     def __init__(self, parent, s):
            ....:         super().__init__(parent)
            sage: class TestSpace(MaassFormSpace):
            ....:     Element = TestElement
            ....:     def to_json(self):
            ....:         return {'group': 'G'}
            sage: S = TestSpace('my_group')
            sage: S.group()
            'my_group'
        """
        return self._group

    def is_cuspidal(self):
        r"""
        Return whether this space consists of cuspidal forms.

        OUTPUT:

        - bool

        EXAMPLES::

            sage: from maass_form_core.spaces.base_space import MaassFormSpace
            sage: from sage.structure.element import Element
            sage: class TestElement(Element):
            ....:     def __init__(self, parent, s):
            ....:         super().__init__(parent)
            sage: class TestSpace(MaassFormSpace):
            ....:     Element = TestElement
            ....:     def to_json(self):
            ....:         return {}
            sage: TestSpace('G', cuspidal=True).is_cuspidal()
            True
            sage: TestSpace('G', cuspidal=False).is_cuspidal()
            False
        """
        return self._cuspidal

    def to_json(self) -> dict:
        r"""
        Return a JSON-serializable representation of this space.

        Subclasses must override this method to include domain-specific data.

        OUTPUT:

        - dict; JSON-compatible dictionary

        EXAMPLES::

            sage: from maass_form_core.spaces.base_space import MaassFormSpace
            sage: from sage.structure.element import Element
            sage: class TestElement(Element):
            ....:     def __init__(self, parent, s):
            ....:         super().__init__(parent)
            sage: class TestSpace(MaassFormSpace):
            ....:     Element = TestElement
            ....:     def to_json(self):
            ....:         return {'group': str(self.group()), 'cuspidal': self.is_cuspidal()}
            sage: S = TestSpace('G', cuspidal=False)
            sage: d = S.to_json()
            sage: d['group'], d['cuspidal']
            ('G', False)
        """
        raise NotImplementedError("Subclasses must implement to_json()")

    @classmethod
    def from_json(cls, data, **kwargs):
        r"""
        Create a space from a JSON-compatible dictionary.

        Subclasses must override this method for domain-specific deserialization.

        INPUT:

        - ``data`` -- dict; JSON data
        - ``**kwargs`` -- additional construction options

        OUTPUT:

        - An instance of the space class

        EXAMPLES::

            sage: from maass_form_core.spaces.base_space import MaassFormSpace
            sage: from sage.structure.element import Element
            sage: class TestElement(Element):
            ....:     def __init__(self, parent, s):
            ....:         super().__init__(parent)
            sage: class TestSpace(MaassFormSpace):
            ....:     Element = TestElement
            ....:     def to_json(self):
            ....:         return {'group': str(self.group()), 'cuspidal': self.is_cuspidal()}
            ....:     @classmethod
            ....:     def from_json(cls, data, **kwargs):
            ....:         return cls(data['group'], cuspidal=data.get('cuspidal', True))
            sage: S = TestSpace('G', cuspidal=False)
            sage: S2 = TestSpace.from_json(S.to_json())
            sage: S2.group()
            'G'
            sage: S2.is_cuspidal()
            False
        """
        raise NotImplementedError("Subclasses must implement from_json()")

    def __hash__(self):
        r"""
        Return a hash value for this space.

        Hash is based on the JSON representation of the space, ensuring
        that equal spaces have equal hashes.

        OUTPUT:

        - integer hash value

        EXAMPLES::

            sage: from maass_form_core.spaces.base_space import MaassFormSpace
            sage: from sage.structure.element import Element
            sage: class TestElement(Element):
            ....:     def __init__(self, parent, s):
            ....:         super().__init__(parent)
            sage: class TestSpace(MaassFormSpace):
            ....:     Element = TestElement
            ....:     def to_json(self):
            ....:         return {'group': str(self.group()), 'cuspidal': self.is_cuspidal()}
            sage: S1 = TestSpace('G', cuspidal=False)
            sage: S2 = TestSpace('G', cuspidal=False)
            sage: hash(S1) == hash(S2)
            True
            sage: S3 = TestSpace('G', cuspidal=True)
            sage: hash(S1) == hash(S3)
            False
        """
        return hash(str(self.to_json()))

    def __eq__(self, other):
        r"""
        Check equality with another space.

        Two spaces are equal if they are of the same type and have the same
        JSON representation.

        INPUT:

        - ``other`` -- object to compare with

        OUTPUT:

        - bool

        EXAMPLES::

            sage: from maass_form_core.spaces.base_space import MaassFormSpace
            sage: from sage.structure.element import Element
            sage: class TestElement(Element):
            ....:     def __init__(self, parent, s):
            ....:         super().__init__(parent)
            sage: class TestSpace(MaassFormSpace):
            ....:     Element = TestElement
            ....:     def to_json(self):
            ....:         return {'group': str(self.group()), 'cuspidal': self.is_cuspidal()}
            sage: S1 = TestSpace('G', cuspidal=False)
            sage: S2 = TestSpace('G', cuspidal=False)
            sage: S1 == S2
            True
            sage: S3 = TestSpace('H', cuspidal=False)
            sage: S1 == S3
            False
            sage: S1 == 'not a space'
            False
        """
        if not isinstance(other, self.__class__):
            return False
        return self.to_json() == other.to_json()

    def __ne__(self, other):
        r"""
        Check inequality with another space.

        INPUT:

        - ``other`` -- object to compare with

        OUTPUT:

        - bool

        EXAMPLES::

            sage: from maass_form_core.spaces.base_space import MaassFormSpace
            sage: from sage.structure.element import Element
            sage: class TestElement(Element):
            ....:     def __init__(self, parent, s):
            ....:         super().__init__(parent)
            sage: class TestSpace(MaassFormSpace):
            ....:     Element = TestElement
            ....:     def to_json(self):
            ....:         return {'group': str(self.group()), 'cuspidal': self.is_cuspidal()}
            sage: S1 = TestSpace('G')
            sage: S2 = TestSpace('G')
            sage: S1 != S2
            False
            sage: S3 = TestSpace('H')
            sage: S1 != S3
            True
        """
        return not self.__eq__(other)

    def _element_constructor_(self, *args, **kwargs):
        r"""
        Construct an element of this space.

        This is the SageMath coercion framework entry point. Subclasses
        may override for domain-specific construction logic.

        EXAMPLES::

            sage: from maass_form_core.spaces.base_space import MaassFormSpace
            sage: from sage.structure.element import Element
            sage: class TestElement(Element):
            ....:     def __init__(self, parent, s):
            ....:         super().__init__(parent)
            ....:         self._s = s
            sage: class TestSpace(MaassFormSpace):
            ....:     Element = TestElement
            ....:     def to_json(self):
            ....:         return {}
            sage: S = TestSpace('G')
            sage: e = S._element_constructor_(42)
            sage: e.parent() is S
            True
        """
        return self.element_class(self, *args, **kwargs)
