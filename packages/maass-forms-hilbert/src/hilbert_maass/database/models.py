"""
Database representation of Hilbert Maass forms.
"""
import mongoengine as me
from comp_manager.document.models import DBObjectBase

class HilbertMaassFormDB(DBObjectBase):
    """
    Hilbert Maass form database object
    """
    meta = {
        'collection': 'hilbert_maass_forms',
        'object_class_name_base' : 'HilbertMaassForm'
    }
    # Properties matching those of HilbertMaassForm_Element
    # and in particular the output of the 'to_json' method
    spectral_parameters = me.ListField(me.DictField())
    coefficients = me.DictField()
    parent = me.DictField()

    def __str__(self, *args, **kwargs):
        """
        String representation of self.
        """
        return self.name

    def save(self, **kwargs):
        pass

