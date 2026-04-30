r"""
Space of Maass waveforms for Kleinian complements.

This module implements the parent structure for spaces of Maass waveforms
associated with Kleinian groups (discrete subgroups of PSL(2,C)) acting
on hyperbolic 3-space.

AUTHORS:

- Fredrik Strömberg (2024): initial version
"""

from typing import Any, Dict, Optional, ParamSpec

from sage.structure.parent import Parent

from maass_forms_klein.exceptions import InvalidGroupError, InvalidSpaceError
from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup, KleinianGroup_class

from .kmaass_element import KleinianMaassFormElement

P = ParamSpec("P")


class KleinianMaassFormSpace(Parent):
    r"""Space of Maass waveforms for Kleinian groups.

    This class represents the parent structure for spaces of Maass waveforms
    associated with discrete subgroups of PSL(2,C). Maass waveforms are
    eigenfunctions of the hyperbolic Laplacian on quotient spaces.

    A Kleinian Maass form is a smooth function f on the upper half-space
    satisfying:

    1. Automorphy: f(γz) = f(z) for all γ in the Kleinian group # ruff: noqa: RUF002
    2. Eigenvalue condition: Δf = λf for some spectral parameter λ
    3. Growth conditions: appropriate decay at cusps (for cuspidal forms)

    ATTRIBUTES:

    - ``group`` -- KleinianGroup; the underlying discrete group
    - ``Element`` -- class; the element class for this space
    - ``_is_cuspidal`` -- bool; whether forms are required to be cuspidal

    EXAMPLES::

        sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
        sage: # Create space for the figure-eight knot complement
        sage: S = KleinianMaassFormSpace('4_1')
        sage: S
        Kleinian cuspidal Maass Form Space (KleinianGroup from Knot: 4_1)
        sage: S.is_cuspidal()
        True
        sage: # Create non-cuspidal space
        sage: S_nc = KleinianMaassFormSpace('4_1', cuspidal=False)
        sage: S_nc
        Kleinian Maass Form Space (KleinianGroup from Knot: 4_1)
        sage: S_nc.is_cuspidal()
        False

    TESTS::

        sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
        sage: S = KleinianMaassFormSpace('4_1')
        sage: S.group() is not None
        True
        sage: hasattr(S, '_is_cuspidal')
        True
    """  # noqa: RUF002

    _group = None
    Element = KleinianMaassFormElement

    def __init__(self, *args: P.args, **kwargs: P.kwargs) -> None:
        r"""Initialize a Kleinian Maass form space.

        INPUT:

        - ``group_or_manifold`` -- KleinianGroup or manifold identifier;
          the underlying discrete group or a string/integer identifying a
          knot/manifold (e.g., '4_1' for figure-eight knot)
        - ``cuspidal`` -- bool (default: True); whether to restrict to cuspidal forms
        - ``base_ring`` -- ring (optional); coefficient ring for the space

        EXAMPLES::

            sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
            sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
            sage: # Create space from manifold identifier
            sage: S1 = KleinianMaassFormSpace('4_1')
            sage: S1
            Kleinian cuspidal Maass Form Space (KleinianGroup from Knot: 4_1)
            sage: # Create space from existing group
            sage: G = KleinianGroup('4_1')
            sage: S2 = KleinianMaassFormSpace(G)
            sage: S2.group() == G
            True

            sage: # Create non-cuspidal space
            sage: S3 = KleinianMaassFormSpace('4_1', cuspidal=False)
            sage: S3.is_cuspidal()
            False

        TESTS::

            sage: # Test invalid input
            sage: KleinianMaassFormSpace()
            Traceback (most recent call last):
            ...
            InvalidSpaceError: Must provide group or manifold identifier

            sage: # Test with invalid group
            sage: KleinianMaassFormSpace(None)
            Traceback (most recent call last):
            ...
            maass_forms_klein.exceptions.InvalidGroupError: Failed to create group from None...
        """
        if not args:
            raise InvalidSpaceError("Must provide group or manifold identifier")

        try:
            if isinstance(args[0], KleinianGroup_class):
                self._group = args[0]
            elif args[0] is None:
                raise InvalidGroupError("Invalid group specification: None")
            else:
                self._group = KleinianGroup(args[0])
        except Exception as e:
            raise InvalidGroupError(
                f"Failed to create group from {args[0]}: {e}", group_data=args[0]
            ) from e

        self._is_cuspidal = kwargs.pop("cuspidal", True)
        super(KleinianMaassFormSpace, self).__init__(*args, **kwargs)

    def __repr__(self) -> str:
        """String representation of the space.

        OUTPUT:
        - String describing the space and its underlying group

        EXAMPLES::

            sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
            sage: S = KleinianMaassFormSpace('4_1')
            sage: repr(S)
            'Kleinian cuspidal Maass Form Space (KleinianGroup from Knot: 4_1)'
        """
        cuspidal_str = "cuspidal " if self._is_cuspidal else ""
        return f"Kleinian {cuspidal_str}Maass Form Space ({self.group()})"

    def group(self) -> KleinianGroup_class:
        """Return the underlying Kleinian group.

        OUTPUT:
        - KleinianGroup; the discrete group for this space

        EXAMPLES::

            sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
            sage: S = KleinianMaassFormSpace('4_1')
            sage: S.group()
            KleinianGroup from Knot: 4_1
        """
        return self._group

    def is_cuspidal(self) -> bool:
        """Return whether this space consists of cuspidal forms.

        Cuspidal forms are those that vanish at all cusps of the quotient space.

        OUTPUT:
        - bool; True if space is restricted to cuspidal forms

        EXAMPLES::

            sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
            sage: S1 = KleinianMaassFormSpace('4_1')  # default cuspidal=True
            sage: S1.is_cuspidal()
            True
            sage: S2 = KleinianMaassFormSpace('4_1', cuspidal=False)
            sage: S2.is_cuspidal()
            False
        """
        return self._is_cuspidal

    def dimension_bound(self, spectral_parameter: Any) -> Optional[int]:
        """Return an upper bound for the dimension at given spectral parameter.

        This method provides theoretical bounds on the dimension of the
        eigenspace for a given spectral parameter.

        INPUT:
        - ``spectral_parameter`` -- complex number; the eigenvalue parameter

        OUTPUT:
        - int or None; upper bound on dimension, or None if unbounded

        EXAMPLES::

            sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
            sage: S = KleinianMaassFormSpace('4_1')
            sage: # For most parameters, expect dimension ≤ 1
            sage: bound = S.dimension_bound(0.5 + 14.1*I)

        .. NOTE::

            The actual computation of dimension bounds is complex and depends
            on the spectral theory of the underlying hyperbolic space.
        """
        # This is a placeholder - actual implementation would require
        # sophisticated spectral theory computations
        return 1 if self._is_cuspidal else None

    def to_json(self, **kwargs: P.kwargs) -> Dict[str, Any]:
        r"""Return JSON-compatible representation of the space.

        This method serializes the space to a dictionary format suitable
        for database storage and network transmission.

        INPUT:
        - ``**kwargs`` -- additional serialization options

        OUTPUT:
        - dict; JSON-compatible representation containing group and space data

        EXAMPLES::

            sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
            sage: S = KleinianMaassFormSpace('4_1')
            sage: data = S.to_json()
            sage: 'group' in data
            True
            sage: 'cuspidal' in data
            True
            sage: data['cuspidal']
            True

        TESTS::

            sage: S = KleinianMaassFormSpace('4_1', cuspidal=False)
            sage: data = S.to_json()
            sage: data['cuspidal']
            False
        """
        return {
            "group": self._group.to_json(**kwargs) if self._group else None,
            "cuspidal": self._is_cuspidal,
            "space_type": "KleinianMaassFormSpace",
        }

    @classmethod
    def from_json(cls, data: Dict[str, Any], **kwargs: P.kwargs) -> "KleinianMaassFormSpace":
        r"""Create space from JSON representation.

        INPUT:
        - ``data`` -- dict; JSON data from to_json()
        - ``**kwargs`` -- additional construction options

        OUTPUT:
        - KleinianMaassFormSpace; reconstructed space

        EXAMPLES::

            sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
            sage: S1 = KleinianMaassFormSpace('4_1')
            sage: data = S1.to_json()
            sage: S2 = KleinianMaassFormSpace.from_json(data)
            sage: S1.is_cuspidal() == S2.is_cuspidal()
            True

        TESTS::

            sage: KleinianMaassFormSpace.from_json({'invalid': 'data'})
            Traceback (most recent call last):
            ...
            InvalidSpaceError: JSON data must contain 'group' key
        """
        if not isinstance(data, dict):
            raise InvalidSpaceError("JSON data must be a dictionary", space_config=data)

        if "group" not in data:
            raise InvalidSpaceError("JSON data must contain 'group' key", space_config=data)

        # Reconstruct group from JSON
        group = KleinianGroup_class.from_json(data["group"]) if data["group"] else None
        if group is None:
            raise InvalidSpaceError("Failed to reconstruct group from JSON", space_config=data)

        # Extract space parameters
        cuspidal = data.get("cuspidal", True)

        return cls(group, cuspidal=cuspidal, **kwargs)

    def _test_space_properties(self) -> bool:
        """Internal method to test basic space properties.

        This method performs consistency checks on the space structure.
        Used primarily for debugging and testing.

        OUTPUT:
        - bool; True if all tests pass

        TESTS::

            sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
            sage: S = KleinianMaassFormSpace('4_1')
            sage: S._test_space_properties()
            True
        """
        # Check that group is properly set
        if self._group is None:
            return False

        # Check that cuspidal flag is boolean
        if not isinstance(self._is_cuspidal, bool):
            return False

        # Check that Element class is properly set
        if self.Element is None:
            return False

        return True

    def _an_element_(self):
        r"""
        Return a typical element of this space.

        This is used by the SageMath coercion framework to obtain a sample
        element for testing purposes.

        OUTPUT:

        - KleinianMaassFormElement; an element with spectral parameter 0

        EXAMPLES::

            sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
            sage: S = KleinianMaassFormSpace('4_1')
            sage: f = S._an_element_()
            sage: f.spectral_parameter()
            0
        """
        return KleinianMaassFormElement(self, 0)

    def _element_constructor_(self, *args: P.args, **kwargs: P.kwargs) -> KleinianMaassFormElement:
        r"""
        Construct an element of this space.

        EXAMPLES::

            sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
            sage: S = KleinianMaassFormSpace('4_1')
            sage: S._element_constructor_(0) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
            Kleinian cuspidal MaassFormElement... from Knot: 4_1), 0)
        """
        return self.element_class(self, *args, **kwargs)
