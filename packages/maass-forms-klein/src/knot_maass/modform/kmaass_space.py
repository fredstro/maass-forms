"""
Space of Maass waveforms for Kleinian complements.

"""
from typing import ParamSpec

from knot_maass.hyperbolic_space.kleinian_group import KleinianGroup_class, KleinianGroup
from sage.structure.parent import Parent

from .kmaass_element import KleinianMaassFormElement

P = ParamSpec('P')


class KleinianMaassFormSpace(Parent):
    """
    Class for spaces of Maass waveform for Kleinian groups.
    """

    group = None
    Element = KleinianMaassFormElement

    def __init__(self,  *args: P.args, **kwargs: P.kwargs) -> None:
        r"""

        INPUT:

        - `args`   -- arguments
        - `kwargs` -- keyword arguments

        EXAMPLES::

            sage: from knot_maass.all import KleinianMaassFormSpace
            sage: S = KleinianMaassFormSpace(-4)
            sage: # TestSuite(S).run()

        TODO: finish this and write examples. The TestSuite(S).run() must pass.
        """
        if isinstance(args[0], KleinianGroup_class):
            self._group = args[0]
        else:
            self._group = KleinianGroup(args[0])
        self._is_cuspidal = kwargs.get('cuspidal', True)
        super(KleinianMaassFormSpace, self).__init__(*args, **kwargs)

    def __repr__(self):
        return f'Kleinian Maass Form Space ({self.group()})'

    def to_json(self, **kwargs: P.kwargs) -> dict:
        r"""
        JSON compatible representation of self.

        EXAMPLES::

            sage: from knot_maass.all import KleinianMaassFormSpace
            sage: S = KleinianMaassFormSpace(-4)
            sage: S.to_json()
            {'cuspidal': True,
             'group': {'_covering_generators': {},
              '_covering_generators_words': [],
              '_gens': [{'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['1.00000000000000', '1.00000000000000'],
                 ['0.000000000000000', '1.00000000000000']]},
               {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['1.00000000000000', '-1.00000000000000*I'],
                 ['0.000000000000000', '1.00000000000000']]},
               {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['0.000000000000000', '-1.00000000000000'],
                 ['1.00000000000000', '0.000000000000000']]},
               {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['1.00000000000000', '-1.00000000000000'],
                 ['0.000000000000000', '1.00000000000000']]},
               {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['1.00000000000000', '1.00000000000000*I'],
                 ['0.000000000000000', '1.00000000000000']]},
               {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['0.000000000000000', '-1.00000000000000'],
                 ['1.00000000000000', '0.000000000000000']]}],
              '_latex_string': 'Bianchi Group: $\\mathbb{Q}(\\sqrt{-4})$',
              '_manifold': None,
              '_name': 'Bianchi Group: Q(sqrt(-4))',
              '_named_gens': {'A': {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['0.000000000000000', '-1.00000000000000'],
                 ['1.00000000000000', '0.000000000000000']]},
               'B': {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['0.000000000000000', '-1.00000000000000'],
                 ['1.00000000000000', '0.000000000000000']]},
               'L': {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['1.00000000000000', '1.00000000000000'],
                 ['0.000000000000000', '1.00000000000000']]},
               'M': {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['1.00000000000000', '-1.00000000000000*I'],
                 ['0.000000000000000', '1.00000000000000']]},
               'a': {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['0.000000000000000', '1.00000000000000'],
                 ['-1.00000000000000', '0.000000000000000']]},
               'b': {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['0.000000000000000', '1.00000000000000'],
                 ['-1.00000000000000', '0.000000000000000']]},
               'l': {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['1.00000000000000', '-1.00000000000000'],
                 ['0.000000000000000', '1.00000000000000']]},
               'm': {'__type__': 'matrix',
                'base_ring': {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53},
                'entries': [['1.00000000000000', '1.00000000000000*I'],
                 ['0.000000000000000', '1.00000000000000']]}},
              '_translation_lattice': {'__type__': 'matrix',
               'base_ring': {'__type__': 'ring', 'name': 'RealField', 'prec': 53},
               'entries': [['1.00000000000000', '0.000000000000000'],
                ['0.000000000000000', '1.00000000000000']]},
              'type': 'KleinianGroup'}}
        """
        return {
            'group': self.group().to_json(),
            'cuspidal': self.is_cuspidal(),
        }

    def _an_element_(self):
        return KleinianMaassFormElement(self, 0)

    def _element_constructor_(self, *args: P.args, **kwargs: P.kwargs) -> KleinianMaassFormElement:
        r"""
        Construct an element of this space.

        EXAMPLES::

            sage: from knot_maass.all import KleinianMaassFormSpace
        """
        return self.element_class(self, *args, **kwargs)

    def is_cuspidal(self) -> bool:
        """

        EXAMPLES::

            sage: from knot_maass.all import KleinianMaassFormSpace
            sage: S = KleinianMaassFormSpace(-4)
            sage: S.is_cuspidal()
            True

        """
        return self._is_cuspidal
    def group(self):
        return self._group