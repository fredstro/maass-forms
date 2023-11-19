"""
Database representation of Hilbert Maass forms.
"""
import logging
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
        if isinstance(item, Integer):
            item = int(item)
        if isinstance(item, slice) and isinstance(item.stop, Integer):
            item = slice(int(item.start), int(item.stop))
        return super().__getitem__(item)

    def near(self, spectral_parameter: tuple[Complex_t],
             max_distance: Real_t = 1e-10) -> QuerySet:
        """
        Find HilbertMaassFormsDB objects near the given spectral parameter.
        """
        lower_bds_x = [float(x.real() - max_distance) for x in spectral_parameter]
        upper_bds_x = [float(x.real() + max_distance) for x in spectral_parameter]
        lower_bds_y = [float(x.imag() - max_distance) for x in spectral_parameter]
        upper_bds_y = [float(x.imag() + max_distance) for x in spectral_parameter]

        conditions = [
            {
                f"spectral_parameter_points.{i}.x": {"$gt": lower_bds_x[i]},
                f"spectral_parameter_points.{i}.y": {"$gt": lower_bds_y[i]}
            }
            for i in range(len(spectral_parameter))
        ]
        conditions += [
            {
                f"spectral_parameter_points.{i}.x": {"$lt": upper_bds_x[i]},
                f"spectral_parameter_points.{i}.y": {"$lt": upper_bds_y[i]}
            }
            for i in range(len(spectral_parameter))
        ]
        return self(__raw__={"$and": conditions})

    # @queryset_manager
    def with_precision(self, m_bound: tuple[Integer_t], y: tuple[Integer_t] = None) -> QuerySet:
        """
        Find HilbertMaassFormsDB objects with coefficient precision bounded by m_bound.

        INPUT:

            queryset:
            min_m:


        """
        conditions = [
            {
                f"coefficients.M.{i}.0": {"$lte": int(m_bound[i][0])},
                f"coefficients.M.{i}.1": {"$gte": int(m_bound[i][1])},
            }
            for i in range(len(m_bound))
        ]
        if y:
            conditions += [
                {
                    f"coefficients.Y.{i}.0": float(y[i]),
                    f"coefficients.Y.{i}.1": float(y[i]),
                }
                for i in range(len(y))
            ]
        return self(__raw__={"$and": conditions}).order_by('-max_m')

    def with_set_coefficients(self, set_coefficients: dict) -> QuerySet:
        return self(__raw__={"set_coefficient": True})


class HilbertMaassFormDB(DBObjectBase):
    """
    Hilbert Maass form database object.
    """
    meta = {
        'collection': 'hilbert_maass_forms',
        'object_class_name_base': 'HilbertMaassForm',
        'queryset_class': HilbertMaassformQuerySet,
    }
    # Properties matching those of HilbertMaassForm_Element
    # and in particular the output of the 'to_json' method
    spectral_parameter = me.ListField(me.DictField())
    # Storing the spectral parameters as list of points on a line enables geo searching
    spectral_parameter_points = me.EmbeddedDocumentListField(Point, default=[])
    y_values = me.ListField(me.FloatField())
    # Describe which coefficients has been set in the normalisation
    set_coefficient = me.DictField()
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
    _skip_keys = ['_id', 'created_at', 'updated_at', 'hash',
                  'coefficients']

    def save(self, **kwargs: P.kwargs):
        """
        Save self.

        INPUT:

        - ``**kwargs``  -- Keyword arguments

        """
        if self.spectral_parameter and not self.spectral_parameter_points:
            complex_pts = [complex(s['val'].replace('*I', 'j').replace(' ', ''))
                           for s in self.spectral_parameter]
            coords = [Point(**{'x': s.real, 'y': s.imag}) for s in complex_pts]
            self.spectral_parameter_points = coords
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
            .with_precision(bound_m, y).with_set_coefficients(set_coefficients).first()
        if not maass_form_db:
            log.debug(f"Compute for s,m,y={spectral_parameter, bound_m, y}")
            space = HilbertMaassFormSpace.from_json(parent)
            maass_form = HilbertMaassForm(space, spectral_parameter)
            maass_form.compute_coefficients(M=bound_m, Y=y, set_coefficients=set_coefficients)
            maass_form_db = insert_object(maass_form)
        return maass_form_db