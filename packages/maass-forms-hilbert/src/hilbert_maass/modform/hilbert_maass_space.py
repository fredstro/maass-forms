"""
Class for spaces of Hilbert Maass forms.
"""
from hilbert_modgroup.all import HilbertModularGroup
from sage.modules.module import Module
from sage.rings.complex_mpfr import ComplexField
from typing import ParamSpec
from .hilbert_maass_element import HilbertMaassForm_Element

P = ParamSpec('P')


class HilbertMaassFormSpace(Module):

    Element = HilbertMaassForm_Element

    def __init__(self, group: HilbertModularGroup, **kwargs: P.kwargs) -> None:
        self._group = group
        self._numerical_precision = kwargs.pop('numerical_precision', 53)
        self._complex_field = ComplexField(self._numerical_precision)
        self._cuspidal = kwargs.pop('cuspidal', False)
        self._number_field = group.base_ring().number_field()
        super(HilbertMaassFormSpace, self).__init__(group.base_ring(), **kwargs)

    def group(self):
        return self._group

    def number_field(self):
        return self._number_field

    def is_cuspidal(self):
        return self._cuspidal