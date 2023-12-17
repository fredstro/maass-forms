"""
Database representation of Hilbert Maass forms.
"""
import hashlib
import logging
from json import dumps
from typing import ParamSpec

import mongoengine as me
from comp_manager.document.models import DBObjectBase
from comp_manager.document.queryset import QuerySetCompat
from comp_manager.utils import insert_object
from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
from hilbert_maass.modform.utils import Real_t, Integer_t, Complex_t
from mongoengine import QuerySet
from sage.all import Integer
from hilbert_maass.modform.hilbert_maass_element import HilbertMaassForm
from sage.rings.number_field.number_field_base import NumberField as NumberField_class

from hilbert_maass.modform.utils import coefficient_dict_to_json

log = logging.getLogger(__name__)

P = ParamSpec('P')


class Point(me.EmbeddedDocument):
    x = me.FloatField()
    y = me.FloatField()

    def __str__(self):
        """
        String representation of self.

        """
        return f"({self.x}, {self.y})"


class HilbertMaassformQuerySet(QuerySetCompat):
    """
    Customised QuerySet for HilbertMaassFormsDB.
    """

    def __getitem__(self, item):
        """
        Get an item from the QuerySet

        INPUT:

        - ``item`` -- integer or slice

        """
        if isinstance(item, Integer):
            item = int(item)
        if isinstance(item, slice) and isinstance(item.stop, Integer):
            item = slice(int(item.start), int(item.stop))
        return super().__getitem__(item)

    def space(self, space: HilbertMaassFormSpace | NumberField_class | dict) -> QuerySet:
        """
        Filter HilbertMaassFormsDB objects by space.

        INPUT:

        - ``space`` -- HilbertMaassFormSpace or dict with 'number_field'


        """
        if isinstance(space, HilbertMaassFormSpace):
            space = space.to_json()
        elif not (isinstance(space, dict) and 'number_field' in space):
            raise TypeError("space must be HilbertMaassFormSpace or dict with 'number_field'")
        return self(parent__number_field=space['number_field'], parent__cuspidal=space['cuspidal'])


    def spectral_range(self, range_real: tuple[tuple[Real_t]],
                       range_imag: tuple[tuple[Real_t]] = None,
                       eps: Real_t = 1e-10) -> QuerySet:
        """
        Filter for spectral parameter in a given range

        INPUT:

        - ``range_real`` -- tuple of tuples of real numbers
        - ``range_imag`` -- tuple of tuples of imaginary numbers


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
                f"spectral_parameter_points.{i}.y": {"$gt": lower_bds_y[i]}
            }
            for i in range(len(range_real))
        ]
        conditions += [
            {
                f"spectral_parameter_points.{i}.x": {"$lt": upper_bds_x[i]},
                f"spectral_parameter_points.{i}.y": {"$lt": upper_bds_y[i]}
            }
            for i in range(len(range_real))
        ]
        return self(__raw__={"$and": conditions})

    def near(self, spectral_parameter: tuple[Complex_t],
             max_distance: Real_t = 1e-10) -> QuerySet:
        """
        Find HilbertMaassFormsDB objects near the given spectral parameter.
        """
        return self.spectral_range([(x.real(),x.real()) for x in spectral_parameter],
                                   [(x.imag(),x.imag()) for x in spectral_parameter],
                                   eps=max_distance)

    # @queryset_manager
    def with_m_precision(self, m_bound: tuple[Integer_t]) -> QuerySet:
        """
        Find HilbertMaassFormsDB objects with coefficient precision bounded by m_bound.

        INPUT:

            queryset:
            min_m:


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
        return self(__raw__={"$and": conditions}).order_by('-max_m')

    def with_y_precision(self, y: tuple[Real_t] = None, eps: Real_t = 1e-10) -> QuerySet:
        """
        Find HilbertMaassFormsDB objects with coefficient precision bounded by m_bound.

        INPUT:

            queryset:
            min_m:


        """
        if not y:
            return self
        if not isinstance(y, (tuple, list)):
            y = [y] * len(y)
        conditions = [
                {
                    f"coefficients.Y.{i}": {
                        "$lte": float(y[i]) + float(eps),
                        "$gte": float(y[i]) - float(eps)
                    }
                }
            for i in range(len(y))
            ]
        return self(__raw__={"$and": conditions}).order_by('-max_m')


    def with_set_coefficients(self, set_coefficients: dict) -> QuerySet:
        set_coefficients_db = coefficient_dict_to_json(set_coefficients)
        return self(__raw__={"coefficients__set_coefficients": set_coefficients_db})


class HilbertMaassFormDB(DBObjectBase):
    """
    Hilbert Maass form database object.
    """
    meta = {
        'collection': 'hilbert_maass_forms',
        'object_class_name_base': 'HilbertMaassForm',
        'queryset_class': HilbertMaassformQuerySet,
        'indexes': [
            {'fields': ('hash',), 'unique': True},
            {'fields': ('parent',), 'unique': False},
        ],
    }
    # Properties matching those of HilbertMaassForm_Element
    # and in particular the output of the 'to_json' method
    spectral_parameter = me.ListField(me.DictField())
    # Storing the spectral parameters as list of floats
    # coressponding to r-values, only used for cusp forms
    r_values = me.ListField(me.FloatField())
    y_values = me.ListField(me.FloatField())
    # Describe which coefficients has been set in the normalisation
    set_coefficients = me.DictField()
    coefficients = me.DictField()
    parent = me.DictField()
    # Set manually (or automatically) to 'tentative' if the form is
    # close to a true eigenvalue, otherwise 'checked'.
    # If it is a known lift we mark it as 'lift'
    status = me.StringField(choices=['tentative', 'unchecked', 'checked',
                                     'lift'],
                            default='unchecked')
    comments = me.StringField()
    max_m = me.IntField()
    # Skip 'coefficients' since we only want to compare against the
    # input values, not the computed values.
    _skip_keys = ['_id', 'created_at', 'updated_at', 'hash', 'comments',
                  'coefficients']

    def save(self, **kwargs: P.kwargs):
        """
        Save self.

        INPUT:

        - ``**kwargs``  -- Keyword arguments

        """
        if self.spectral_parameter and not self.r_values:
            complex_pts = [complex(s['val'].replace('*I', 'j').replace(' ', ''))
                           for s in self.spectral_parameter]
            #coords = [Point(**{'x': s.real, 'y': s.imag}) for s in complex_pts]
            #self.spectral_parameter_points = coords
            self.r_values = [float(s.imag) for s in complex_pts]
        if self.coefficients and not self.y_values:
            self.y_values = [float(y) for y in self.coefficients['Y']]
        if not self.max_m and self.coefficients:
            self.max_m = max(max(m) for m in self.coefficients['M'])
        super(HilbertMaassFormDB, self).save(**kwargs)

    def __str__(self, *args: P.args, **kwargs: P.kwargs) -> str:
        """
        String representation of self.
        """
        poly = self.parent.get('number_field', {}).get('polynomial', '')
        spectral_parameter = [s.get('val', "") for s in self.spectral_parameter]
        return f"Hilbert Maass form for NumberField({poly}) with spectral parameter" \
           f" {spectral_parameter}"


    @classmethod
    def near_or_create(cls, parent: HilbertMaassFormSpace, spectral_parameter: tuple[Complex_t],
                       max_distance: Real_t=1e-10,
                       bound_m: tuple[Integer_t] = None,
                       y: tuple[Real_t] = None,
                       set_coefficients: dict = None) -> 'HilbertMaassFormDB':
        """
        Find or create HilbertMaassFormsDB objects near the given spectral parameter.
        """
        if not isinstance(parent, dict):
            parent = parent.to_json()
        maass_form_db = cls.objects(parent=parent).near(spectral_parameter,
                                                        max_distance=max_distance)\
            .with_m_precision(bound_m).with_y_precision(y).with_set_coefficients(set_coefficients).first()
        if not maass_form_db:
            log.debug(f"Compute for s,m,y={spectral_parameter, bound_m, y}")
            space = HilbertMaassFormSpace.from_json(parent)
            maass_form = HilbertMaassForm(space, spectral_parameter)
            maass_form.compute_coefficients(M=bound_m, Y=y, set_coefficients=set_coefficients)
            maass_form_db = insert_object(maass_form)
        return maass_form_db