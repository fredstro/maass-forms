from knot_maass.modform.utils import Integer_t, Real_t, Complex_t
from sage.categories.sets_cat import cartesian_product
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.structure.element import Vector, Matrix
from sage.structure.parent import Parent

from knot_maass.hyperbolic_space.utils import P


class KleinianMaassFormCoefficients(Parent):
    """
    Class for coefficients of Kleinian Maass forms.
    """

    def __init__(self, coefficients: list | Matrix, M: Integer_t,
                 spectral_parameter: Complex_t | Real_t,
                 space: 'KleinianMaassFormSpace',
                 Y: Real_t = None,
                 coordinate_indices: list = None,
                 coordinate_values: list = None,
                 set_coefficients: dict = None,
                 check: bool = True, **kwargs: P.kwargs) -> None:
        r"""
        INPUT:

        - `args`   -- arguments
        - `kwargs` -- keyword arguments

        EXAMPLES::

            sage: from knot_maass.modform.kmaass_space import KleinianMaassFormSpace
            sage: from knot_maass.modform.coefficients import KleinianMaassFormCoefficients
            sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
            sage: H = KleinianMaassFormSpace(-4)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = Matrix(RR, [1,2,3,4,5,6,7,8,9])
            sage: C = KleinianMaassFormCoefficients(Cmat, 1, spectral_parameter, H)
            sage: #TestSuite(C).run()

        """
        super(KleinianMaassFormCoefficients, self).__init__(**kwargs)
        if isinstance(coefficients, (Vector, list)):
            coefficients = matrix(coefficients).transpose()
        self._coefficients = coefficients
        self._M = M
        self._spectral_parameter = spectral_parameter
        self._space = space
        self._Y = Y
        if not coordinate_indices:
            coordinate_indices = cartesian_product([range(-M, M+1), range(-M, M+1)])
        self._coordinate_indices = [vector(x, immutable=True) for x in coordinate_indices]
        if not coordinate_values:
            coordinate_values = space.group().dual_translation_lattice_vectors(M)
        self._coordinate_values = [vector(x, immutable=True) for x in coordinate_values]
        self._set_coefficients = set_coefficients or {}
        self._check = check

    def to_json(self, **kwargs: P.kwargs) -> dict:
        r"""
        JSON compatible representation of self.

        EXAMPLES:

            sage: from knot_maass.modform.kmaass_space import KleinianMaassFormSpace
            sage: from knot_maass.modform.coefficients import KleinianMaassFormCoefficients
            sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
            sage: H = KleinianMaassFormSpace(-4)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = Matrix(RR, [1,2,3,4,5,6,7,8,9])
            sage: C = KleinianMaassFormCoefficients(Cmat, 1, spectral_parameter, H)
            sage: C.to_json()
            {'M': 1,
             'Y': 'None',
             'coefficients': ['(1.00000000000000, 2.00000000000000, 3.00000000000000, ...
             'coordinate_indices': ['(-1, -1)',
              '(-1, 0)',
              '(-1, 1)',
              '(0, -1)',
              '(0, 0)',
              '(0, 1)',
              '(1, -1)',
              '(1, 0)',
              '(1, 1)'],
             'coordinate_values': ['(0.000000000000000, 0.000000000000000)',
              '(-1.00000000000000, 0.000000000000000)',
              '(0.000000000000000, -1.00000000000000)',
              '(0.000000000000000, 1.00000000000000)',
              '(1.00000000000000, 0.000000000000000)',
              '(-1.00000000000000, -1.00000000000000)',
              '(-1.00000000000000, 1.00000000000000)',
              '(1.00000000000000, -1.00000000000000)',
              '(1.00000000000000, 1.00000000000000)'],
             'prec': 53,
             'set_coefficients': {},
             'space': {'cuspidal': True, 'group': 'Bianchi Group: Q(sqrt(-4))'},
             'spectral_parameter': '(0.500000000000000 + 1.00000000000000*I, 0.50000000000000...
        """
        json_dict = {
                'M': self._M,
                'Y': str(self._Y),
                'coefficients': [str(x) for x in self._coefficients],
                'prec': int(self._coefficients.base_ring().prec()),
                'spectral_parameter': str(self._spectral_parameter),
                'set_coefficients': {str(k): str(v) for k, v in self._set_coefficients.items()},
                'coordinate_indices': [str(x) for x in self._coordinate_indices],
                'coordinate_values': [str(x) for x in self._coordinate_values],
            }
        if kwargs.get('include_space', False):
            json_dict['space'] = self._space.to_json()
        else:
            json_dict['space'] = {
                'cuspidal': self._space.is_cuspidal(),
                'group': self._space.group().name()
            }
        return json_dict

    def coordinate_indices(self):
        r"""
        Return coordinate indices.

        EXAMPLES::

            sage: from knot_maass.all import KleinianMaassFormSpace, KleinianMaassFormCoefficients
            sage: H = KleinianMaassFormSpace(-4)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = Matrix(RR, [1,2,3,4,5,6,7,8,9])
            sage: C = KleinianMaassFormCoefficients(Cmat, 1, spectral_parameter, H)
            sage: C.coordinate_indices()
            [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]
        """
        return self._coordinate_indices

    def coordinate_values(self):
        r"""
        Return coordinate values.

        EXAMPLES::

            sage: from knot_maass.all import KleinianMaassFormSpace, KleinianMaassFormCoefficients
            sage: H = KleinianMaassFormSpace(-4)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = Matrix(RR, [1,2,3,4,5,6,7,8,9])
            sage: C = KleinianMaassFormCoefficients(Cmat, 1, spectral_parameter, H)
            sage: C.coordinate_values()
             [(0.000000000000000, 0.000000000000000),
             (-1.00000000000000, 0.000000000000000),
             (0.000000000000000, -1.00000000000000),
             (0.000000000000000, 1.00000000000000),
             (1.00000000000000, 0.000000000000000),
             (-1.00000000000000, -1.00000000000000),
             (-1.00000000000000, 1.00000000000000),
             (1.00000000000000, -1.00000000000000),
             (1.00000000000000, 1.00000000000000)]
        """
        return self._coordinate_values

    def __repr__(self):
        """
        String representation of self.

        EXAMPLES::

            sage: from knot_maass.all import KleinianMaassFormSpace, KleinianMaassFormCoefficients
            sage: H = KleinianMaassFormSpace(-4)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = Matrix(RR, [1,2,3,4,5,6,7,8,9])
            sage: C = KleinianMaassFormCoefficients(Cmat, 1, spectral_parameter, H)
            sage: C.__repr__()
            'Coefficients of a Kleinian Maass form with M=1'


        """
        return f"Coefficients of a Kleinian Maass form with M={self._M}"

    def __hash__(self):
        self._coefficients.set_immutable()
        return hash((self._coefficients, self._Y, self._M, str(self._index_tuples),
                     self._spectral_parameter, self._space,
                     str(self._set_coefficients)))

    def __getitem__(self, item) -> Complex_t:
        """
        Return the coefficient of the given index or slice.

        INPUT:

        - ``item`` -- index or slice

        EXAMPLES:

            sage: from knot_maass.all import KleinianMaassFormSpace, KleinianMaassFormCoefficients
            sage: H = KleinianMaassFormSpace(-4)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = vector(RR, [1,2,3,4,5,6,7,8,9])
            sage: C = KleinianMaassFormCoefficients(Cmat, 1, spectral_parameter, H)
            sage: C[0]
            1.00000000000000
            sage: C[0:3]
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: C[(1.0,1.0)]
            9.00000000000000
        """
        if isinstance(item, (Integer_t, slice)):
            return self._coefficients.column(0)[item]
        index = None
        if isinstance(item, (tuple, list)):
            item = vector(item)
        if item in self._coordinate_indices:
            index = self._coordinate_indices.index(item)
        elif item in self._coordinate_values:
            index = self._coordinate_values.index(item)
        else:
            raise KeyError(f"No index {item} in {self}")
        return self._coefficients.column(0)[index]

    def __iter__(self):
        yield from self._coefficients

    def keys(self, as_elements=False) -> list:
        """
        Return the list of keys.
        :param as_elements:
        :return:

        EXAMPLES:

            sage: from knot_maass.all import KleinianMaassFormSpace, KleinianMaassFormCoefficients
            sage: H = KleinianMaassFormSpace(-4)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = vector(RR, [1,2,3,4,5,6,7,8,9])
            sage: C = KleinianMaassFormCoefficients(Cmat, 1, spectral_parameter, H)
            sage: len(C.keys())
            9
            sage: C.keys()
            [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]
            sage: C.keys(as_elements=True)
             [(0.000000000000000, 0.000000000000000),
             (-1.00000000000000, 0.000000000000000),
             (0.000000000000000, -1.00000000000000),
             (0.000000000000000, 1.00000000000000),
             (1.00000000000000, 0.000000000000000),
             (-1.00000000000000, -1.00000000000000),
             (-1.00000000000000, 1.00000000000000),
             (1.00000000000000, -1.00000000000000),
             (1.00000000000000, 1.00000000000000)]
            sage: dict(C)
            {(-1, -1): 1.00000000000000,
             (-1, 0): 2.00000000000000,
             (-1, 1): 3.00000000000000,
             (0, -1): 4.00000000000000,
             (0, 0): 5.00000000000000,
             (0, 1): 6.00000000000000,
             (1, -1): 7.00000000000000,
             (1, 0): 8.00000000000000,
             (1, 1): 9.00000000000000}
        """
        if as_elements:
            return self._coordinate_values
        return self._coordinate_indices

    def values(self) -> list:
        r"""
        Return the values.

        EXAMPLE:

            sage: from knot_maass.all import KleinianMaassFormSpace, KleinianMaassFormCoefficients
            sage: H = KleinianMaassFormSpace(-4)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = vector(RR, [1,2,3,4,5,6,7,8,9])
            sage: C = KleinianMaassFormCoefficients(Cmat, 1, spectral_parameter, H)
            sage: list(C.values())
            [1.00000000000000,
             2.00000000000000,
             3.00000000000000,
             4.00000000000000,
             5.00000000000000,
             6.00000000000000,
             7.00000000000000,
             8.00000000000000,
             9.00000000000000]
        """
        return list(self._coefficients.column(0))
