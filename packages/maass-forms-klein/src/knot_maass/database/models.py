import logging
import mongoengine as me
from knot_maass.modform.kmaass_element import KleinianMaassFormElement
from mongoengine import QuerySet
from sage.rings.complex_mpfr import ComplexField
from sage.rings.integer import Integer

from comp_manager.core.models import DBObjectBase
from comp_manager.core.queryset import QuerySetCompat
from comp_manager.utils import insert_object
from knot_maass.hyperbolic_space.utils import P
from knot_maass.modform.kmaass_space import KleinianMaassFormSpace
from knot_maass.modform.utils import Real_t, Complex_t, Integer_t, map_tuple_to_int


log = logging.getLogger(__name__)


class Point(me.EmbeddedDocument):
    x = me.FloatField()
    y = me.FloatField()

    def __str__(self):
        """
        String representation of self.

        """
        return f"({self.x}, {self.y})"

class ParallelogramDB(me.EmbeddedDocument):
    base = me.ListField(me.FloatField())
    v1 = me.ListField(me.FloatField())
    v2 = me.ListField(me.FloatField())


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
    meta = {
        'indexes': [{'fields': ('label', 'reduced_to_fd', 'max_length'),
                    'unique': True}]
    }

    def save(self, **kwargs: P.kwargs) -> None:
        if not self.max_length:
            self.max_length = max([len(w) for w in self.words])
        super(Word, self).save(**kwargs)

    def __repr__(self):
        return f"Words({self.label}, {self.max_length}, {self.reduced_to_fd}, {self.fd})"

class KleinianMaassFormQuerySet(QuerySetCompat):
    """
    Customised QuerySet for KleinianMaassFormsDB.
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

    def space(self, space: KleinianMaassFormSpace | dict) -> QuerySet:
        """
        Filter KleinianMaassFormsDB objects by space.

        INPUT:

        - ``space`` -- KleinianMaassFormSpace or dict with 'number_field'


        """
        if isinstance(space, KleinianMaassFormSpace):
            space = space.to_json()
        elif not (isinstance(space, dict) and 'group' in space):
            raise TypeError("space must be KleinianMaassFormSpace or dict with 'group'")
        return self(parent__group=space['group'], parent__cuspidal=space['cuspidal'])

    def spectral_range(self, range_real: tuple[Real_t, Real_t],
                       range_imag: tuple[Real_t, Real_t] = None,
                       eps: Real_t = 1e-15) -> QuerySet:
        """
        Filter for spectral parameter in a given range

        INPUT:

        - ``range_real`` -- tuple of tuples of real numbers
        - ``range_imag`` -- tuple of tuples of imaginary numbers


        """
        lower_bds_x = float(range_real[0] - eps)
        upper_bds_x = float(range_real[1] + eps)
        conditions = {
            f"spectral_parameter_point.x": {"$gte": lower_bds_x, "$lte": upper_bds_x},
        }
        if range_imag:
            lower_bds_y = float(range_imag[0] - eps)
            upper_bds_y = float(range_imag[1] + eps)
            conditions[f"spectral_parameter_point.y"] = {"$gte": lower_bds_y, "$lte": upper_bds_y}
        return self(__raw__=conditions)

    def near(self, spectral_parameter: Complex_t,
             max_distance: Real_t = 1e-15) -> QuerySet:
        """
        Find KleinianMaassFormsDB objects near the given spectral parameter.
        """
        return self.spectral_range((spectral_parameter.real(), spectral_parameter.real()),
                                   (spectral_parameter.imag(), spectral_parameter.imag()),
                                   eps=max_distance)

    def with_m_precision(self, m_bound: tuple[Integer_t]) -> QuerySet:
        """
        Find KleinianMaassFormsDB objects with coefficient precision bounded by m_bound.

        INPUT:

            queryset:
            min_m:


        """
        conditions = {}
        if m_bound:
            conditions = {
                    f"coefficients.M": {"$lte": int(m_bound[0]), "$gte": int(m_bound[1])}
            }
        return self(__raw__=conditions).order_by('-max_m')

    def with_y_precision(self, y: tuple[Real_t] = None, eps: Real_t = 1e-15) -> QuerySet:
        """
        Find KleinianMaassFormsDB objects with coefficient precision bounded by m_bound.

        INPUT:

            queryset:
            min_m:


        """
        if not y:
            return self
        if not isinstance(y, (tuple, list)):
            y = [y] * 2
        conditions = {f"coefficients.Y": {
                        "$lte": float(y[1]) + float(eps),
                        "$gte": float(y[0]) - float(eps)
                    }}
        return self(__raw__={"$and": conditions})


    # def with_set_coefficients(self, set_coefficients: dict) -> QuerySet:
    #     set_coefficients_db = coefficient_dict_to_json(set_coefficients)
    #     return self(__raw__={"coefficients.set_coefficients": set_coefficients_db})

class KleinianMaassFormDB(DBObjectBase):
    """
    Kleinian Maass form database object.
    """
    meta = {
        'collection': 'kleinian_maass_forms',
        'object_class_name_base': 'KleinianMaassForm',
        'queryset_class': KleinianMaassFormQuerySet,
        'indexes': [
            {'fields': ('hash',), 'unique': True},
            {'fields': ('parent',), 'unique': False},
        ],
    }
    # Properties matching those of KleinianMaassForm_Element
    # and in particular the output of the 'to_json' method
    spectral_parameter = me.ListField(me.DictField())
    # spectral_parameter_points = me.EmbeddedDocumentListField(Point, default=[])
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
    status = me.StringField(choices=['tentative', 'unchecked', 'checked',
                                     'lift'],
                            default='unchecked')
    comments = me.StringField()
    max_m = me.IntField()
    # Skip 'coefficients' since we only want to compare against the
    # input values, not the computed values.
    _skip_keys = ['_id', 'created_at', 'updated_at', 'hash', 'comments',
                  'coefficients.coefficients', 'spectral_parameter_points']

    def save(self, **kwargs: P.kwargs):
        """
        Save self.

        INPUT:

        - ``**kwargs``  -- Keyword arguments

        """
        if self.spectral_parameter and not self.r_values:
            complex_pts = [complex(s['val'].replace('*I', 'j').replace(' ', ''))
                           for s in self.spectral_parameter]
            coords = [Point(**{'x': s.real, 'y': s.imag}) for s in complex_pts]
            self.spectral_parameter_points = coords
            self.r_values = [float(s.imag) for s in complex_pts]
        if self.coefficients and not self.y_values:
            self.y_values = [float(y) for y in self.coefficients['Y']]
        if not self.max_m and self.coefficients:
            self.max_m = max(max(m) for m in self.coefficients['M'])
        super(KleinianMaassFormDB, self).save(**kwargs)

    def __str__(self, *args: P.args, **kwargs: P.kwargs) -> str:
        """
        String representation of self.
        """
        poly = self.parent.get('number_field', {}).get('polynomial', '')
        spectral_parameter = [s.get('val', "") for s in self.spectral_parameter]
        return f"Kleinian Maass form for NumberField({poly}) with spectral parameter" \
           f" {spectral_parameter}"

    def coefficient(self, t: tuple[Integer_t]) -> Complex_t:
        """
        Return the coefficient corresponding to the given tuple.

        """
        bound_tuple = integer_to_bounds_tuple(self.max_m, len(self.spectral_parameter))
        n = map_tuple_to_int(t, bound_tuple)
        prec = self.spectral_parameter[0]['prec']
        return ComplexField(prec)(self.coefficients['coefficients'][n][0])

    @classmethod
    def near_or_create(cls, parent: KleinianMaassFormSpace, spectral_parameter: tuple[Complex_t],
                       max_distance: Real_t=1e-15,
                       bound_m: tuple[Integer_t] = None,
                       y: tuple[Real_t] = None,
                       set_coefficients: dict = None) -> 'KleinianMaassFormDB':
        """
        Find or create KleinianMaassFormsDB objects near the given spectral parameter.
        """
        if not isinstance(parent, dict):
            parent = parent.to_json()
        maass_form_db = cls.objects(parent=parent).near(spectral_parameter,
                                                        max_distance=max_distance)\
            .with_m_precision(bound_m).with_y_precision(y).with_set_coefficients(set_coefficients).first()
        if not maass_form_db:
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


class KleinianGroupDB(DBObjectBase):
    """
    Database representation of KleinianGroup.
    """
    name = me.StringField()
    _name_string = me.StringField()
    _gens = me.ListField(me.DictField())  #me.EmbeddedDocumentField(MatrixDB))
    _covering_generators_words = me.ListField(me.StringField())
    _named_gens = me.DictField()
    _covering_generators = me.DictField()
    _manifold = me.StringField()
    _translation_lattice = me.DictField()
    _latex_string = me.StringField()
    type = me.StringField()
    meta = {
        'collection': 'kleinian_group',
        'object_class_name_base': 'KleinianGroup',
    }

    def __repr__(self):
        return self._name_string

    def save(self, **kwargs: P.kwargs) -> 'BaseDocument':
        """
        Save self.

        INPUT:

        - ``**kwargs``  -- Keyword arguments

        """
        if not self.name:
            self.name = self._name_string
        super(KleinianGroupDB, self).save(**kwargs)
