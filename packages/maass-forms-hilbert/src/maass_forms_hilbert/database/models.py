"""
Database representation of Hilbert Maass forms.
"""

import logging
from typing import Any, ClassVar, ParamSpec

import mongoengine as me
from comp_manager.core.models import DBObjectBaseAbstract
from comp_manager.core.queryset import QuerySetCompat
from comp_manager.utils import insert_object
from mongoengine import QuerySet
from sage.all import Integer
from sage.functions.other import imag, real
from sage.rings.complex_mpfr import ComplexField
from sage.rings.number_field.number_field_base import NumberField as NumberField_class

from maass_forms_hilbert.modform.hilbert_maass_element import HilbertMaassForm
from maass_forms_hilbert.modform.hilbert_maass_space import HilbertMaassFormSpace
from maass_forms_hilbert.modform.utils import (
    Complex_t,
    Integer_t,
    Real_t,
    coefficient_dict_to_json,
    integer_to_bounds_tuple,
    map_tuple_to_int,
)

log = logging.getLogger(__name__)

P = ParamSpec("P")


class Point(me.EmbeddedDocument):
    x = me.FloatField()
    y = me.FloatField()

    def __str__(self):
        """
        Return the string representation of the Point object.

        OUTPUT:
        A string in the format ``(x, y)``, where ``x`` and ``y`` are the coordinates of the point.

        EXAMPLES::

            sage: import mongoengine as me
            sage: class Point(me.EmbeddedDocument):
            ....:     x = me.FloatField()
            ....:     y = me.FloatField()
            ....:     def __str__(self):
            ....:         return f"({self.x}, {self.y})"
            sage: point = Point(x=1.5, y=2.5)
            sage: str(point)
            '(1.5, 2.5)'
        """
        return f"({self.x}, {self.y})"


class HilbertMaassformQuerySet(QuerySetCompat):
    """
    Customised QuerySet for HilbertMaassFormsDB.
    """

    def __getitem__(self, item):
        r"""
        Retrieve an item or a slice of items from the QuerySet.

        INPUT:
        - ``item`` -- an integer index or a slice object.

        OUTPUT:
        The item at the specified index or a sliced QuerySet.

        EXAMPLES::

            sage: from maass_forms_hilbert.database.models import HilbertMaassFormDB
            sage: from maass_forms_hilbert.database.tests import connect_mockdb, insert_fixtures
            sage: connect_mockdb()
            sage: insert_fixtures()
            sage: qs = HilbertMaassFormDB.objects
            sage: qs[0] # doctest: +NORMALIZE_WHITESPACE
            <HilbertMaassFormDB: Hilbert Maass form for NumberField(x^2 - 5) with spectral
                parameter ['0.500000000000000 + 1.00000000000000*I',
                '0.500000000000000 + 1.00000000000000*I']>
            sage: qs[0].parent["number_field"]
             {'polynomial': 'x^2 - 5', 'names': ['a']}
            sage: qs[0:1] # doctest: +NORMALIZE_WHITESPACE
             [<HilbertMaassFormDB: Hilbert Maass form for NumberField(x^2 - 5) with spectral
                parameter ['0.500000000000000 + 1.00000000000000*I',
                '0.500000000000000 + 1.00000000000000*I']>]

        """
        if isinstance(item, Integer):
            item = int(item)
        if isinstance(item, slice) and isinstance(item.stop, Integer):
            item = slice(int(item.start), int(item.stop))
        return super().__getitem__(item)

    def space(self, space: HilbertMaassFormSpace | NumberField_class | dict) -> QuerySet:
        r"""
        Filter HilbertMaassFormsDB objects by space.

        INPUT:

        - ``space`` -- HilbertMaassFormSpace or dict with 'number_field'

        OUTPUT:

        - QuerySet filtered by the given space

        EXAMPLES::

            sage: from maass_forms_hilbert.database.models import HilbertMaassFormDB
            sage: from maass_forms_hilbert.database.tests import connect_mockdb, insert_fixtures
            sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
            sage: connect_mockdb()
            sage: insert_fixtures()
            sage: qs = HilbertMaassFormDB.objects
            sage: space = HilbertMaassFormSpace(7)
            sage: qs.space(space)
            []
            sage: space = HilbertMaassFormSpace(5)
            sage: len(qs.space(space))
            6
            sage: space = HilbertMaassFormSpace(2)
            sage: len(qs.space(space))
            3

        TESTS::

            sage: qs.space({'invalid': 'dict'})
            Traceback (most recent call last):
            ...
            TypeError: space must be HilbertMaassFormSpace or dict with 'number_field'

        """
        if isinstance(space, HilbertMaassFormSpace):
            space = space.to_json()
        elif not (isinstance(space, dict) and "number_field" in space):
            raise TypeError("space must be HilbertMaassFormSpace or dict with 'number_field'")
        return self(
            parent__number_field=space["number_field"],
            parent__cuspidal=space["cuspidal"],
        )

    def spectral_range(
        self,
        range_real: list[tuple[Real_t]],
        range_imag: list[tuple[Real_t]] = None,
        eps: Real_t = 1e-15,
    ) -> QuerySet:
        r"""
        Filter for spectral parameter in a given range.

        INPUT:

        - ``range_real`` -- tuple of tuples of real numbers
        - ``range_imag`` -- tuple of tuples of imaginary numbers
        - ``eps`` -- (default: 1e-15) tolerance for range

        OUTPUT:

        - QuerySet filtered by the given spectral range

        EXAMPLES::

            sage: from maass_forms_hilbert.database.models import HilbertMaassFormDB
            sage: from maass_forms_hilbert.database.tests import connect_mockdb, insert_fixtures
            sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
            sage: connect_mockdb()
            sage: insert_fixtures()
            sage: qs = HilbertMaassFormDB.objects
            sage: tol = 1e-15
            sage: len(qs.spectral_range(range_real=[(0.4, 0.55), (0.4, 0.65)],
            ....:                     range_imag=[(0.9, 1.1), (0.9, 1.1)]))
            9
            sage: len(qs.spectral_range(range_real=[(0.4, 0.5), (0.4, 0.55)],
            ....:                     range_imag=[(0.9, 1.1), (0.9, 1.1)],
            ....:                     eps=1e-16))
            6
            sage: len(qs.spectral_range(range_real=[(0.5, 0.5), (0.5, 0.5)],
            ....:                     range_imag=[(0.9, 1.1), (0.9, 1.1)],
            ....:                     eps=1e-16))
            6
            sage: len(qs.spectral_range(range_real=[(0.5-1e-10, 0.5-1e-10), (0.5, 0.5)],
            ....:                     range_imag=[(0.9, 1.1), (0.9, 1.1)],
            ....:                     eps=1e-16))
            0
            sage: len(qs.spectral_range(range_real=[(0.5-1e-10, 0.5-1e-10), (0.5, 0.5)],
            ....:                     range_imag=[(0.9, 1.1), (0.9, 1.1)],
            ....:                     eps=1e-8))
            6
        """
        lower_bds_x = [float(x[0] - eps) for x in range_real]
        upper_bds_x = [float(x[1] + eps) for x in range_real]
        if not range_imag:
            range_imag = [0.5] * len(range_real)
        lower_bds_y = [float(y[0] - eps) for y in range_imag]
        upper_bds_y = [float(y[1] + eps) for y in range_imag]

        conditions = [
            {
                f"spectral_parameter_points.{i}.x": {"$gt": lower_bds_x[i]},
                f"spectral_parameter_points.{i}.y": {"$gt": lower_bds_y[i]},
            }
            for i in range(len(range_real))
        ]
        conditions += [
            {
                f"spectral_parameter_points.{i}.x": {"$lt": upper_bds_x[i]},
                f"spectral_parameter_points.{i}.y": {"$lt": upper_bds_y[i]},
            }
            for i in range(len(range_real))
        ]
        return self(__raw__={"$and": conditions})

    def near(self, spectral_parameter: tuple[Complex_t], max_distance: Real_t = 1e-15) -> QuerySet:
        r"""
        Find HilbertMaassFormsDB objects near the given spectral parameter.

        INPUT:

        - ``spectral_parameter`` -- tuple of complex numbers
        - ``max_distance`` -- (default: 1e-15) maximum allowed distance

        OUTPUT:

        - QuerySet of objects near the given spectral parameter

        EXAMPLES::

            sage: from maass_forms_hilbert.database.models import HilbertMaassFormDB
            sage: from maass_forms_hilbert.database.tests import connect_mockdb, insert_fixtures
            sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
            sage: connect_mockdb()
            sage: insert_fixtures()
            sage: qs = HilbertMaassFormDB.objects
            sage: tol = 1e-15
            sage: len(qs.near((CC(0.5,1.0), CC(0.5,1.0)), max_distance=1e-10))
            6
            sage: len(qs.near((CC(0.5,1.00001), CC(0.5,1.0)), max_distance=1e-10))
            0
            sage: len(qs.near((CC(0.5,1.00001), CC(0.5,1.0)), max_distance=0.001))
            6

        """
        return self.spectral_range(
            [(real(x), real(x)) for x in spectral_parameter],
            [(imag(x), imag(x)) for x in spectral_parameter],
            eps=max_distance,
        )

    def with_m_precision(self, m_bound: tuple[Integer_t]) -> QuerySet:
        r"""
        Find HilbertMaassFormsDB objects with coefficient precision bounded by m_bound.

        INPUT:

        - ``m_bound`` -- tuple of integer bounds

        OUTPUT:

        - QuerySet filtered by m_bound

        EXAMPLES::

            sage: from maass_forms_hilbert.database.models import HilbertMaassFormDB
            sage: from maass_forms_hilbert.database.tests import connect_mockdb, insert_fixtures
            sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
            sage: connect_mockdb()
            sage: insert_fixtures()
            sage: qs = HilbertMaassFormDB.objects
            sage: qs.with_m_precision(((-10, 10), (-10, 10)))
            []
            sage: l = qs.with_m_precision(((-1, 1), (-1, 1)))
            sage: [x.max_m for x in l]
             [5, 5, 5, 4, 4, 4, 1, 1, 1]
            sage: l =qs.with_m_precision(((-2, 2), (-2, 2)))
            sage: [x.max_m for x in l]
            [5, 5, 5, 4, 4, 4]
        """
        conditions = []
        if m_bound:
            conditions += [
                {
                    f"coefficients.M.{i}.0": {"$lte": int(m_bound[i][0])},
                    f"coefficients.M.{i}.1": {"$gte": int(m_bound[i][1])},
                }
                for i in range(len(m_bound))
            ]
        return self(__raw__={"$and": conditions}).order_by("-max_m")

    def with_q_precision(self, q: tuple[Integer_t] = None) -> QuerySet:
        r"""
        Find HilbertMaassFormsDB objects with coefficient associated to Q.

        INPUT:

        - ``q`` -- tuple of integer Q values

        OUTPUT:

        - QuerySet filtered by Q

        EXAMPLES::

            sage: from maass_forms_hilbert.database.models import HilbertMaassFormDB
            sage: from maass_forms_hilbert.database.tests import connect_mockdb, insert_fixtures
            sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
            sage: connect_mockdb()
            sage: insert_fixtures()
            sage: qs = HilbertMaassFormDB.objects
            sage: qs.with_q_precision((10,10))
            []
            sage: l = qs.with_q_precision((2, 2))
            sage: sorted([x.q_values[0] for x in l])
             [2, 2, 2, 5, 5, 5, 6, 6, 6]
            sage: l = qs.with_q_precision((6, 6))
            sage: sorted([x.q_values[0] for x in l])
            [6, 6, 6]

        """
        if not q:
            return self
        conditions = [{f"coefficients.Q.{i}": {"$gte": int(q[i])}} for i in range(len(q))]
        return self(__raw__={"$and": conditions}).order_by("-Q")

    def with_y_precision(self, y: tuple[Real_t] = None, eps: Real_t = 1e-15) -> QuerySet:
        r"""
        Find HilbertMaassFormsDB objects with coefficient precision bounded by y.

        INPUT:

        - ``y`` -- tuple of real values
        - ``eps`` -- (default: 1e-15) tolerance

        OUTPUT:

        - QuerySet filtered by y

        EXAMPLES::

            sage: from maass_forms_hilbert.database.models import HilbertMaassFormDB
            sage: from maass_forms_hilbert.database.tests import connect_mockdb, insert_fixtures
            sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
            sage: connect_mockdb()
            sage: insert_fixtures()
            sage: qs = HilbertMaassFormDB.objects
            sage: l = qs.with_y_precision((0.3,0.3))
            sage: sorted(set(   [x.y_values[0] for x in l]))
            [0.3]
            sage: l = qs.with_y_precision((0.2, 0.2))
            sage: sorted(set([x.y_values[0] for x in l]))
            [0.2]
            sage: qs.with_y_precision((0.1, 0.1))
            []
        """
        if not y:
            return self
        if not isinstance(y, (tuple, list)):
            y = [y] * len(y)
        conditions = [
            {
                f"coefficients.Y.{i}": {
                    "$lte": float(y[i]) + float(eps),
                    "$gte": float(y[i]) - float(eps),
                }
            }
            for i in range(len(y))
        ]
        return self(__raw__={"$and": conditions})

    def with_set_coefficients(self, set_coefficients: dict) -> QuerySet:
        r"""
        Filter HilbertMaassFormsDB objects by a set of coefficients.

        INPUT:

        - ``set_coefficients`` -- dictionary of coefficients

        OUTPUT:

        - QuerySet filtered by set_coefficients

        EXAMPLES::

            sage: from maass_forms_hilbert.database.models import HilbertMaassFormDB
            sage: from maass_forms_hilbert.database.tests import connect_mockdb, insert_fixtures
            sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
            sage: connect_mockdb()
            sage: insert_fixtures()
            sage: qs = HilbertMaassFormDB.objects
            sage: len(qs.with_set_coefficients({(1,1): 1, (-1,1): 2}))
            6
            sage: len(qs.with_set_coefficients({(1,1): -1, (-1,1): 2}))
            3
            sage: len(qs.with_set_coefficients({(1,1): -1, (-1,1): -2}))
            0

        """
        set_coefficients_db = coefficient_dict_to_json(set_coefficients)
        return self(__raw__={"coefficients.set_coefficients": set_coefficients_db})


class HilbertMaassFormDB(DBObjectBaseAbstract):
    """
    Hilbert Maass form database object.
    """

    meta: ClassVar[dict[str, Any]] = {
        "collection": "hilbert_maass_forms",
        "object_class_name_base": "HilbertMaassForm",
        "queryset_class": HilbertMaassformQuerySet,
        "indexes": [
            {"fields": ("hash",), "unique": True},
            {"fields": ("parent",), "unique": False},
        ],
    }
    # Properties matching those of HilbertMaassForm_Element
    # and in particular the output of the 'to_json' method
    spectral_parameter = me.ListField(me.DictField())
    spectral_parameter_points = me.EmbeddedDocumentListField(Point, default=[])
    # Storing the spectral parameters as list of floats
    # corresponding to r-values, only used for cusp forms
    r_values = me.ListField(me.FloatField())
    y_values = me.ListField(me.FloatField())
    q_values = me.ListField(me.IntField())
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
    _skip_keys: ClassVar[list[str]] = [
        "_id",
        "created_at",
        "updated_at",
        "hash",
        "comments",
        "coefficients.coefficients",
        "spectral_parameter_points",
    ]

    def save(self, **kwargs: P.kwargs):
        r"""
        Save the Hilbert Maass form to the database with additional processing.

        This method computes derived fields like r_values, y_values, max_m, and q_values
        before saving the object to the database.

        INPUT:

        - ``**kwargs`` -- Keyword arguments passed to the parent save method

        EXAMPLES::

            sage: from maass_forms_hilbert.database.models import HilbertMaassFormDB
            sage: from maass_forms_hilbert.database.tests import connect_mockdb, return_db_parameters
            sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
            sage: connect_mockdb()
            sage: data = return_db_parameters()
            sage: F = HilbertMaassFormDB(**data)
            sage: F.save()
            sage: F.r_values
            [1.0, 1.0]
            sage: F.y_values
            [0.3, 0.3]
            sage: F.max_m
            1
            sage: F.q_values
            [2, 2]
        """
        if self.coefficients:
            if not self.spectral_parameter:
                self.spectral_parameter = self.coefficients["spectral_parameter"]
            if not self.y_values:
                self.y_values = [float(y) for y in self.coefficients["Y"]]
            if not self.max_m:
                self.max_m = max(max(m) for m in self.coefficients["M"])
            if not self.q_values:
                self.q_values = [int(q) for q in self.coefficients["Q"]]

        if self.spectral_parameter and not self.r_values:
            complex_pts = [
                complex(s["val"].replace("*I", "j").replace(" ", ""))
                for s in self.spectral_parameter
            ]
            coords = [Point(**{"x": s.real, "y": s.imag}) for s in complex_pts]
            self.spectral_parameter_points = coords
            self.r_values = [float(s.imag) for s in complex_pts]
        super().save(**kwargs)

    def __str__(self, *args: P.args, **kwargs: P.kwargs) -> str:
        r"""
        Return a string representation of the Hilbert Maass form.

        INPUT:

        - ``*args`` -- Variable length argument list
        - ``**kwargs`` -- Arbitrary keyword arguments

        OUTPUT:

        - String representation of the Hilbert Maass form

        EXAMPLES::

            sage: from maass_forms_hilbert.database.models import HilbertMaassFormDB
            sage: from maass_forms_hilbert.database.tests import connect_mockdb, insert_fixtures
            sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
            sage: connect_mockdb()
            sage: insert_fixtures()
            sage: F = HilbertMaassFormDB.objects[0]
            sage: str(F)
             "Hilbert Maass form for NumberField(x^2 - 5) with spectral parameter
             ['0.500000000000000 + 1.00000000000000*I', '0.500000000000000 + 1.00000000000000*I']"
        """
        poly = self.parent.get("number_field", {}).get("polynomial", "")
        spectral_parameter = [s.get("val", "") for s in self.spectral_parameter]
        return (
            f"Hilbert Maass form for NumberField({poly}) with spectral parameter"
            f" {spectral_parameter}"
        )

    def coefficient(self, t: tuple[Integer_t]) -> Complex_t:
        r"""
        Return the coefficient of the Hilbert Maass form at the specified index.

        This method maps the input tuple to an integer index in the coefficient array
        and returns the corresponding complex coefficient.

        INPUT:

        - ``t`` -- tuple of integers representing the coefficient index

        OUTPUT:

        - Complex number representing the coefficient at index t

        EXAMPLES::

            sage: from maass_forms_hilbert.database.models import HilbertMaassFormDB
            sage: from maass_forms_hilbert.database.tests import connect_mockdb, insert_fixtures
            sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
            sage: connect_mockdb()
            sage: insert_fixtures()
            sage: F = HilbertMaassFormDB.objects[0]
            sage: F.save()
            sage: F.coefficient((1, 1))
            1.00000000000000 + 1.00000000000000*I
        """
        bound_tuple = integer_to_bounds_tuple(self.max_m, len(self.spectral_parameter))
        n = map_tuple_to_int(t, bound_tuple)
        prec = self.spectral_parameter[0]["prec"]
        return ComplexField(prec)(self.coefficients["coefficients"][n][0])

    @classmethod
    def near_or_create(
        cls,
        parent: HilbertMaassFormSpace,
        spectral_parameter: tuple[Complex_t],
        max_distance: Real_t = 1e-15,
        bound_m: tuple[Integer_t] = None,
        y: tuple[Real_t] = None,
        q: tuple[Integer_t] = None,
        set_coefficients: dict = None,
    ) -> "HilbertMaassFormDB":
        r"""
        Find an existing Hilbert Maass form or create a new one if none exists.

        This method first tries to find a Hilbert Maass form in the database with the
        specified parameters. If none is found, it computes a new one and stores it.

        INPUT:

        - ``parent`` -- HilbertMaassFormSpace or dict; the parent space
        - ``spectral_parameter`` -- tuple of complex numbers; spectral parameter
        - ``max_distance`` -- real number (default: 1e-15); maximum distance for finding nearby forms
        - ``bound_m`` -- tuple of integers (default: None); bounds for the Fourier coefficients
        - ``y`` -- tuple of real numbers (default: None); y-values for evaluation
        - ``q`` -- tuple of integers (default: None); Q-values
        - ``set_coefficients`` -- dict (default: None); coefficients to set

        OUTPUT:

        - HilbertMaassFormDB; the found or created database object

        EXAMPLES::

            sage: from maass_forms_hilbert.database.models import HilbertMaassFormDB
            sage: from maass_forms_hilbert.database.tests import connect_mockdb, insert_fixtures
            sage: from maass_forms_hilbert.all import HilbertMaassFormSpace
            sage: connect_mockdb()
            sage: insert_fixtures()
            sage: HilbertMaassFormDB.objects.count()
            10
            sage: F = HilbertMaassFormDB.near_or_create(
            ....: parent=HilbertMaassFormSpace.from_json(
            ....: {"number_field": {"polynomial": "x^2 - 5", "names": ["a"]},
            ....:  "cuspidal": False}),
            ....: spectral_parameter=(CC(0.5, 1.0001), CC(0.5, 1.0)),
            ....: max_distance=1e-3,
            ....: bound_m=((-1, 1), (-1, 1)),
            ....: y=(0.3, 0.3),
            ....: q=(2, 2),
            ....: set_coefficients={(1,1): 1, (-1,1): 2}
            ....: )
            sage: HilbertMaassFormDB.objects.count()
            10
            sage: F = HilbertMaassFormDB.near_or_create(
            ....: parent=HilbertMaassFormSpace.from_json(
            ....: {"number_field": {"polynomial": "x^2 - 5", "names": ["a"]},
            ....:  "cuspidal": False}),
            ....: spectral_parameter=(CC(0.5, 2.0), CC(0.5, 1.0)),
            ....: max_distance=1e-15,
            ....: bound_m=((-1, 1), (-1, 1)),
            ....: y=(0.3, 0.3),
            ....: q=(2, 2),
            ....: set_coefficients={(1,1): 1, (-1,1): 2}
            ....: )
            sage: HilbertMaassFormDB.objects.count()
            11
        """
        if not isinstance(parent, dict):
            parent = parent.to_json()
        maass_form_db = (
            cls.objects(parent=parent)
            .near(spectral_parameter, max_distance=max_distance)
            .with_m_precision(bound_m)
            .with_y_precision(y)
            .with_q_precision(q)
            .with_set_coefficients(set_coefficients)
            .first()
        )
        if not maass_form_db:
            log.debug(f"Compute for s,m,y={spectral_parameter, bound_m, y}")
            space = HilbertMaassFormSpace.from_json(parent)
            maass_form = HilbertMaassForm(space, spectral_parameter)
            maass_form.compute_coefficients(M=bound_m, Y=y, Q=q, set_coefficients=set_coefficients)
            maass_form_db = insert_object(maass_form)
        return maass_form_db
