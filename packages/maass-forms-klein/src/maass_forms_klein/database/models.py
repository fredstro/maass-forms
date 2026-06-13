import logging

# Lazy imports to avoid circular dependencies
from typing import TYPE_CHECKING, ClassVar, ParamSpec, Union

import mongoengine as me
from comp_manager.core.models import DBObjectBase, DBObjectBaseAbstract
from comp_manager.core.queryset import QuerySetCompat
from comp_manager.utils import insert_object
from mongoengine import QuerySet
from sage.rings.complex_mpfr import ComplexField
from sage.rings.integer import Integer

from maass_forms_klein.exceptions import InvalidSpaceError, ValidationError
from maass_forms_klein.modform.utils import Complex_t, Integer_t, Real_t, map_tuple_to_int

if TYPE_CHECKING:
    from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace

P = ParamSpec("P")


log = logging.getLogger(__name__)


class Point(me.EmbeddedDocument):
    """Embedded document representing a point in the complex plane.

    This class stores coordinates for spectral parameters and other
    complex-valued mathematical objects in the database.

    ATTRIBUTES:
    - ``x`` -- Real; the real part (x-coordinate)
    - ``y`` -- Real; the imaginary part (y-coordinate)

    EXAMPLES::

        sage: from maass_forms_klein.all import Point
        sage: pt = Point(x=0.5, y=14.1)
        sage: str(pt)
        '(0.5, 14.1)'
    """

    x = me.FloatField(required=True, help_text="Real part of the complex number")
    y = me.FloatField(required=True, help_text="Imaginary part of the complex number")

    def __str__(self) -> str:
        """String representation of the point.

        OUTPUT:
        - String in the format "(x, y)"

        EXAMPLES::

            sage: from maass_forms_klein.all import Point
            sage: pt = Point(x=1.5, y=-2.3)
            sage: str(pt)
            '(1.5, -2.3)'
        """
        return f"({self.x}, {self.y})"


class ParallelogramDB(me.EmbeddedDocument):
    """Embedded document representing a parallelogram in the hyperbolic plane.

    This class stores the geometric data for fundamental domains and
    other parallelogram regions used in hyperbolic geometry computations.

    ATTRIBUTES:
    - ``base`` -- List[Real]; base point coordinates [x, y]
    - ``v1`` -- List[Real]; first spanning vector [x, y]
    - ``v2`` -- List[Real]; second spanning vector [x, y]

    EXAMPLES::

        sage: from maass_forms_klein.all import ParallelogramDB
        sage: pg = ParallelogramDB(base=[0,0], v1=[1,0], v2=[0,1])
        sage: pg.base
        [0.0, 0.0]
    """

    base = me.ListField(
        me.FloatField(), max_length=2, min_length=2, help_text="Base point coordinates [x, y]"
    )
    v1 = me.ListField(
        me.FloatField(), max_length=2, min_length=2, help_text="First spanning vector [x, y]"
    )
    v2 = me.ListField(
        me.FloatField(), max_length=2, min_length=2, help_text="Second spanning vector [x, y]"
    )


class Word(DBObjectBase):
    """
    Reduced words in the fundamental group of a manifold given by a know complement
    """

    label = me.StringField()
    words = me.ListField(me.StringField())
    reduced_to_fd = me.BooleanField()
    reduced_cover = me.BooleanField()
    fd = me.EmbeddedDocumentField(ParallelogramDB)
    max_length = me.IntField()
    group = me.DictField()
    meta: ClassVar[dict] = {
        "indexes": [{"fields": ("label", "reduced_to_fd", "max_length"), "unique": True}],
    }

    def save(self, **kwargs: P.kwargs) -> None:
        r"""
        Save this Word document, computing ``max_length`` if not set.

        INPUT:

        - ``**kwargs`` -- keyword arguments passed to the parent save method

        EXAMPLES::

            sage: from maass_forms_klein.database.models import Word  # doctest: +SKIP
            sage: w = Word(label='4_1', words=['aB', 'Ab'])  # doctest: +SKIP
            sage: w.save()  # doctest: +SKIP
        """
        if not self.max_length:
            self.max_length = max([len(w) for w in self.words])
        super(Word, self).save(**kwargs)

    def __repr__(self):
        r"""
        Return a string representation of this Word document.

        OUTPUT:

        - string; a representation showing label, max_length, and reduction status

        EXAMPLES::

            sage: from maass_forms_klein.database.models import Word  # doctest: +SKIP
            sage: w = Word(label='4_1', words=['aB'], max_length=2)  # doctest: +SKIP
            sage: repr(w)  # doctest: +SKIP
            "Words(4_1, 2, None, None)"
        """
        return f"Words({self.label}, {self.max_length}, {self.reduced_to_fd}, {self.fd})"


class KleinianMaassFormQuerySet(QuerySetCompat):
    """Enhanced QuerySet for KleinianMaassFormsDB with comprehensive validation.

    This QuerySet provides specialized query methods for Kleinian Maass forms
    with proper input validation, error handling, and type safety.

    EXAMPLES::

        sage: from maass_forms_klein.database.models import KleinianMaassFormDB
        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: qs = KleinianMaassFormDB.objects
        sage: # Query for forms in spectral range
        sage: forms = qs.spectral_range((0.5, 0.6), (14.0, 15.0))
    """

    def __getitem__(self, item):
        """
        Get an item from the QuerySet.

        INPUT:

        - ``item`` -- integer or slice

        EXAMPLES::

            sage: from maass_forms_klein.database.models import (  # doctest: +SKIP
            ....:     KleinianMaassFormDB)
            sage: from maass_form_core.testing import connect_mockdb  # doctest: +SKIP
            sage: connect_mockdb()  # doctest: +SKIP
            sage: qs = KleinianMaassFormDB.objects  # doctest: +SKIP
            sage: qs[0]  # doctest: +SKIP
        """
        if isinstance(item, Integer):
            item = int(item)
        if isinstance(item, slice) and isinstance(item.stop, Integer):
            item = slice(int(item.start), int(item.stop))
        return super().__getitem__(item)

    def space(self, space: Union["KleinianMaassFormSpace", dict]) -> QuerySet:
        """Filter KleinianMaassFormsDB objects by space.

        INPUT:
        - ``space`` -- KleinianMaassFormSpace or dict; space specification
          If dict, must contain 'group' key for filtering

        OUTPUT:
        - QuerySet filtered by the specified space parameters

        EXAMPLES::

            sage: from maass_forms_klein.all import KleinianMaassFormSpace, KleinianMaassFormDB
            sage: space = KleinianMaassFormSpace('4_1')
            sage: qs = KleinianMaassFormDB.objects.space(space)

        TESTS::
            sage: from maass_form_core.testing import connect_mockdb
            sage: connect_mockdb()
            sage: qs = KleinianMaassFormDB.objects
            sage: qs.space({'invalid': 'dict'})
            Traceback (most recent call last):
            ...
            InvalidSpaceError: Dict must contain 'group' key
        """
        # Import here to avoid circular imports
        from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace

        if isinstance(space, KleinianMaassFormSpace):
            space = space.to_json()
        elif isinstance(space, dict):
            if "group" not in space:
                raise InvalidSpaceError("Dict must contain 'group' key", space_config=space)
        else:
            raise ValidationError(
                f"space must be KleinianMaassFormSpace or dict with 'group', got {type(space)}",
                field_name="space",
                value=space,
            )
        return self(parent__group=space["group"], parent__cuspidal=space["cuspidal"])

    def spectral_range(
        self,
        range_real: tuple[Real_t, Real_t],
        range_imag: tuple[Real_t, Real_t] | None = None,
        eps: Real_t = 1e-15,
    ) -> QuerySet:
        """Filter for spectral parameters within given ranges.

        INPUT:
        - ``range_real`` -- Tuple[Real_t, Real_t]; bounds for real part (min, max)
        - ``range_imag`` -- Tuple[Real_t, Real_t] or None (default: None); bounds for imaginary part
        - ``eps`` -- Real_t (default: 1e-15); tolerance for range matching

        OUTPUT:
        - QuerySet filtered by spectral parameter ranges

        EXAMPLES::

            sage: from maass_forms_klein.all import KleinianMaassFormDB
            sage: from maass_form_core.testing import connect_mockdb
            sage: connect_mockdb()
            sage: qs = KleinianMaassFormDB.objects
            sage: # Find forms with real part in [0.4, 0.6]
            sage: forms = qs.spectral_range((0.4, 0.6))
            sage: # Find forms in complex rectangle
            sage: forms = qs.spectral_range((0.4, 0.6), (14.0, 15.0))

        TESTS::
            sage: from maass_form_core.testing import connect_mockdb
            sage: connect_mockdb()
            sage: qs = KleinianMaassFormDB.objects
            sage: qs.spectral_range((0.9, 0.1))  # Invalid range
            Traceback (most recent call last):
            ...
            ValidationError: Lower bound must be less than or equal to upper bound
        """
        # Validate input ranges
        if not isinstance(range_real, (tuple, list)) or len(range_real) != 2:
            raise ValidationError(
                "range_real must be a tuple/list of length 2",
                field_name="range_real",
                value=range_real,
            )

        if range_real[0] > range_real[1]:
            raise ValidationError(
                "Lower bound must be less than or equal to upper bound",
                field_name="range_real",
                value=range_real,
            )

        if range_imag is not None:
            if not isinstance(range_imag, (tuple, list)) or len(range_imag) != 2:
                raise ValidationError(
                    "range_imag must be a tuple/list of length 2",
                    field_name="range_imag",
                    value=range_imag,
                )
            if range_imag[0] > range_imag[1]:
                raise ValidationError(
                    "Lower bound must be less than or equal to upper bound",
                    field_name="range_imag",
                    value=range_imag,
                )

        # Build query conditions
        lower_bds_x = float(range_real[0] - eps)
        upper_bds_x = float(range_real[1] + eps)
        conditions = {
            "spectral_parameter_point.x": {"$gte": lower_bds_x, "$lte": upper_bds_x},
        }

        if range_imag:
            lower_bds_y = float(range_imag[0] - eps)
            upper_bds_y = float(range_imag[1] + eps)
            conditions["spectral_parameter_point.y"] = {"$gte": lower_bds_y, "$lte": upper_bds_y}

        return self(__raw__=conditions)

    def near(self, spectral_parameter: Complex_t, max_distance: Real_t = 1e-15) -> QuerySet:
        """Find KleinianMaassFormsDB objects near the given spectral parameter.

        INPUT:
        - ``spectral_parameter`` -- Complex_t; target spectral parameter
        - ``max_distance`` -- Real_t (default: 1e-15); maximum distance tolerance

        OUTPUT:
        - QuerySet containing forms within max_distance of the parameter

        EXAMPLES::

            sage: from maass_forms_klein.all import KleinianMaassFormDB
            sage: from sage.rings.complex_mpfr import ComplexField
            sage: from maass_form_core.testing import connect_mockdb
            sage: connect_mockdb()
            sage: CC = ComplexField(53)
            sage: s = CC(0.5, 14.1)
            sage: qs = KleinianMaassFormDB.objects
            sage: nearby_forms = qs.near(s, max_distance=1e-10)

        TESTS::
            sage: from maass_form_core.testing import connect_mockdb
            sage: connect_mockdb()
            sage: qs = KleinianMaassFormDB.objects
            sage: qs.near("invalid")
            Traceback (most recent call last):
            ...
            ValidationError: spectral_parameter must be a complex number
        """
        if not hasattr(spectral_parameter, "real") or not hasattr(spectral_parameter, "imag"):
            raise ValidationError(
                "spectral_parameter must be a complex number",
                field_name="spectral_parameter",
                value=spectral_parameter,
            )

        if max_distance < 0:
            raise ValidationError(
                "max_distance must be non-negative", field_name="max_distance", value=max_distance
            )

        return self.spectral_range(
            (spectral_parameter.real(), spectral_parameter.real()),
            (spectral_parameter.imag(), spectral_parameter.imag()),
            eps=max_distance,
        )

    def with_m_precision(self, m_bound: Integer_t | tuple[Integer_t, Integer_t] | None) -> QuerySet:
        """
        Find KleinianMaassFormsDB objects with coefficient precision at least
        ``m_bound``.

        INPUT:

        - ``m_bound`` -- integer (the truncation ``M`` for Kleinian Maass
          forms), a ``(low, high)`` tuple, or ``None`` to skip filtering.
          When an integer is given, the filter matches forms whose stored
          ``max_m`` is at least ``m_bound``.

        OUTPUT:

        - QuerySet filtered by coefficient bound, ordered by descending max_m

        EXAMPLES::

            sage: from maass_forms_klein.database.models import (  # doctest: +SKIP
            ....:     KleinianMaassFormDB)
            sage: from maass_form_core.testing import connect_mockdb  # doctest: +SKIP
            sage: connect_mockdb()  # doctest: +SKIP
            sage: qs = KleinianMaassFormDB.objects  # doctest: +SKIP
            sage: qs.with_m_precision(5)  # doctest: +SKIP
        """
        conditions = {}
        if m_bound is None:
            pass
        elif isinstance(m_bound, (tuple, list)):
            if len(m_bound) == 2:
                lo, hi = int(m_bound[0]), int(m_bound[1])
                if lo > hi:
                    lo, hi = hi, lo
                conditions = {"max_m": {"$gte": lo, "$lte": hi}}
        else:
            conditions = {"max_m": {"$gte": int(m_bound)}}
        return self(__raw__=conditions).order_by("-max_m")

    def with_set_coefficients(self, set_coefficients: dict | None) -> QuerySet:
        """
        Filter KleinianMaassFormsDB objects by the stored ``set_coefficients``
        normalisation.

        INPUT:

        - ``set_coefficients`` -- dict mapping coefficient indices to their
          fixed values, or ``None`` to skip filtering.

        OUTPUT:

        - QuerySet filtered by ``set_coefficients`` equality.

        EXAMPLES::

            sage: from maass_forms_klein.database.models import (  # doctest: +SKIP
            ....:     KleinianMaassFormDB)
            sage: from maass_form_core.testing import connect_mockdb  # doctest: +SKIP
            sage: connect_mockdb()  # doctest: +SKIP
            sage: qs = KleinianMaassFormDB.objects  # doctest: +SKIP
            sage: qs.with_set_coefficients({(0, 0): 0, (1, 0): 1})  # doctest: +SKIP
        """
        if not set_coefficients:
            return self
        return self(set_coefficients=set_coefficients)

    def with_y_precision(self, y: tuple[Real_t] | None = None, eps: Real_t = 1e-15) -> QuerySet:
        """
        Find KleinianMaassFormsDB objects with Y-value precision in the given range.

        INPUT:

        - ``y`` -- tuple of reals or None (default: None); bounds for Y value
        - ``eps`` -- real (default: 1e-15); tolerance

        OUTPUT:

        - QuerySet filtered by Y-value range

        EXAMPLES::

            sage: from maass_forms_klein.database.models import (  # doctest: +SKIP
            ....:     KleinianMaassFormDB)
            sage: from maass_form_core.testing import connect_mockdb  # doctest: +SKIP
            sage: connect_mockdb()  # doctest: +SKIP
            sage: qs = KleinianMaassFormDB.objects  # doctest: +SKIP
            sage: qs.with_y_precision((0.5, 1.0))  # doctest: +SKIP
        """
        if not y:
            return self
        if not isinstance(y, (tuple, list)):
            y = [y] * 2
        conditions = {
            "coefficients.Y": {"$lte": float(y[1]) + float(eps), "$gte": float(y[0]) - float(eps)}
        }
        return self(__raw__={"$and": conditions})


class KleinianMaassFormDB(DBObjectBaseAbstract):
    """
    Kleinian Maass form database object.
    """

    meta: ClassVar[dict] = {
        "collection": "kleinian_maass_forms",
        "object_class_name_base": "KleinianMaassForm",
        "queryset_class": KleinianMaassFormQuerySet,
        "indexes": [
            {"fields": ("hash",), "unique": True},
            {"fields": ("parent",), "unique": False},
        ],
    }
    # Properties matching those of KleinianMaassForm_Element
    # and in particular the output of the 'to_json' method
    spectral_parameter = me.ListField(me.DictField())
    # Storing the spectral parameters as list of floats
    # coressponding to r-values, only used for cusp forms
    r_value = me.FloatField()
    y_value = me.FloatField()
    # Describe which coefficients has been set in the normalisation
    set_coefficients = me.DictField()
    coefficients = me.DictField()
    parent = me.DictField()
    # Set manually (or automatically) to 'tentative' if the form is
    # close to a true eigenvalue, otherwise 'checked'.
    # If it is a known lift we mark it as 'lift'
    status = me.StringField(
        choices=["tentative", "unchecked", "checked", "lift"], default="unchecked"
    )
    comments = me.StringField()
    max_m = me.IntField()
    # Skip 'coefficients' since we only want to compare against the
    # input values, not the computed values.
    _skip_keys: ClassVar[list] = [
        "_id",
        "created_at",
        "updated_at",
        "hash",
        "comments",
        "coefficients.coefficients",
        "spectral_parameter_points",
    ]

    def save(self, **kwargs: P.kwargs):
        """
        Save this Kleinian Maass form to the database.

        Automatically computes spectral parameter points, r-values,
        y-values, and max_m if not already set.

        INPUT:

        - ``**kwargs`` -- keyword arguments passed to the parent save method

        EXAMPLES::

            sage: from maass_forms_klein.database.models import (  # doctest: +SKIP
            ....:     KleinianMaassFormDB)
            sage: form = KleinianMaassFormDB()  # doctest: +SKIP
            sage: form.save()  # doctest: +SKIP
        """
        if self.spectral_parameter and not self.r_values:
            complex_pts = [
                complex(s["val"].replace("*I", "j").replace(" ", ""))
                for s in self.spectral_parameter
            ]
            coords = [Point(**{"x": s.real, "y": s.imag}) for s in complex_pts]
            self.spectral_parameter_points = coords
            self.r_values = [float(s.imag) for s in complex_pts]
        if self.coefficients and not self.y_values:
            self.y_values = [float(y) for y in self.coefficients["Y"]]
        if not self.max_m and self.coefficients:
            self.max_m = max(max(m) for m in self.coefficients["M"])
        super(KleinianMaassFormDB, self).save(**kwargs)

    def __str__(self, *args: P.args, **kwargs: P.kwargs) -> str:
        """
        Return a string representation of this Kleinian Maass form.

        OUTPUT:

        - string; a description including the number field and spectral parameter

        EXAMPLES::

            sage: from maass_forms_klein.database.models import (  # doctest: +SKIP
            ....:     KleinianMaassFormDB)
            sage: form = KleinianMaassFormDB()  # doctest: +SKIP
            sage: str(form)  # doctest: +SKIP
        """
        poly = self.parent.get("number_field", {}).get("polynomial", "")
        spectral_parameter = [s.get("val", "") for s in self.spectral_parameter]
        return (
            f"Kleinian Maass form for NumberField({poly}) with spectral parameter"
            f" {spectral_parameter}"
        )

    def coefficient(self, t: tuple[Integer_t]) -> Complex_t:
        """
        Return the coefficient corresponding to the given index tuple.

        INPUT:

        - ``t`` -- tuple of integers; the coefficient index

        OUTPUT:

        - complex number; the Fourier coefficient at index ``t``

        EXAMPLES::

            sage: from maass_forms_klein.database.models import (  # doctest: +SKIP
            ....:     KleinianMaassFormDB)
            sage: form = KleinianMaassFormDB.objects.first()  # doctest: +SKIP
            sage: c = form.coefficient((0, 0))  # doctest: +SKIP
        """
        bound_tuple = ((-self.max_m, self.max_m),) * len(self.spectral_parameter)
        n = map_tuple_to_int(t, bound_tuple)
        prec = self.spectral_parameter[0]["prec"]
        return ComplexField(prec)(self.coefficients["coefficients"][n][0])

    @classmethod
    def near_or_create(
        cls,
        parent: "KleinianMaassFormSpace",
        spectral_parameter: tuple[Complex_t],
        max_distance: Real_t = 1e-15,
        bound_m: tuple[Integer_t] | None = None,
        y: tuple[Real_t] | None = None,
        set_coefficients: dict | None = None,
    ) -> "KleinianMaassFormDB":
        """
        Find or create KleinianMaassFormsDB objects near the given spectral parameter.

        INPUT:

        - ``parent`` -- KleinianMaassFormSpace or dict; the parent space
        - ``spectral_parameter`` -- tuple of complex numbers; the target spectral parameter
        - ``max_distance`` -- real (default: 1e-15); maximum distance tolerance
        - ``bound_m`` -- tuple of integers or None; coefficient bound
        - ``y`` -- tuple of reals or None; Y-value range
        - ``set_coefficients`` -- dict or None; normalisation coefficients

        OUTPUT:

        - KleinianMaassFormDB; existing or newly computed form

        EXAMPLES::

            sage: from maass_forms_klein.database.models import (  # doctest: +SKIP
            ....:     KleinianMaassFormDB)
            sage: from maass_forms_klein.all import KleinianMaassFormSpace  # doctest: +SKIP
            sage: space = KleinianMaassFormSpace('4_1')  # doctest: +SKIP
            sage: s = (0.5 + 14.1*I,)  # doctest: +SKIP
            sage: form = KleinianMaassFormDB.near_or_create(space, s)  # doctest: +SKIP
        """
        if not isinstance(parent, dict):
            parent = parent.to_json()
        maass_form_db = (
            cls.objects(parent=parent)
            .near(spectral_parameter, max_distance=max_distance)
            .with_m_precision(bound_m)
            .with_y_precision(y)
            .with_set_coefficients(set_coefficients)
            .first()
        )
        if not maass_form_db:
            from maass_forms_klein.modform.kmaass_element import (
                KleinianMaassFormElement,
            )
            from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace

            log.debug(f"Compute for s,m,y={spectral_parameter, bound_m, y}")
            space = KleinianMaassFormSpace.from_json(parent)
            maass_form = KleinianMaassFormElement(space, spectral_parameter)
            maass_form.compute_coefficients(M=bound_m, Y=y, set_coefficients=set_coefficients)
            maass_form_db = insert_object(maass_form)
        return maass_form_db


class MatrixDB(me.DynamicEmbeddedDocument):
    """
    Database representation of Matrix.
    """

    base_ring = me.DictField()
    entries = me.ListField(me.ListField())


class KleinianGroupDB(DBObjectBaseAbstract):
    """
    Database representation of KleinianGroup.
    """

    name = me.StringField()
    _name_string = me.StringField()
    _gens = me.ListField(me.DictField())  # me.EmbeddedDocumentField(MatrixDB))
    _covering_generators_words = me.ListField(me.StringField())
    _named_gens = me.DictField()
    _covering_generators = me.DictField()
    _manifold = me.StringField()
    _translation_lattice = me.DictField()
    _latex_string = me.StringField()
    type = me.StringField()
    meta: ClassVar[dict] = {
        "collection": "kleinian_group",
        "object_class_name_base": "KleinianGroup",
    }

    def __repr__(self):
        r"""
        Return a string representation of this KleinianGroupDB document.

        OUTPUT:

        - string; the name string of the group

        EXAMPLES::

            sage: from maass_forms_klein.database.models import KleinianGroupDB  # doctest: +SKIP
            sage: g = KleinianGroupDB(_name_string='4_1')  # doctest: +SKIP
            sage: repr(g)  # doctest: +SKIP
            '4_1'
        """
        return self._name_string

    def save(self, **kwargs: P.kwargs) -> None:
        """
        Save this KleinianGroupDB document, setting name from _name_string if needed.

        INPUT:

        - ``**kwargs`` -- keyword arguments passed to the parent save method

        EXAMPLES::

            sage: from maass_forms_klein.database.models import KleinianGroupDB  # doctest: +SKIP
            sage: g = KleinianGroupDB(_name_string='4_1')  # doctest: +SKIP
            sage: g.save()  # doctest: +SKIP
        """
        if not self.name:
            self.name = self._name_string
        super(KleinianGroupDB, self).save(**kwargs)
