import json
from sage.matrix.constructor import matrix
from sage.rings.complex_mpfr import ComplexNumber, ComplexField
from sage.rings.real_mpfr import RealField
from typing import ParamSpec

import logging
from hilbert_maass.modform.utils import ideal_coordinates, map_tuple_to_int, map_int_to_tuple
from sage.rings.integer import Integer
from sage.rings.number_field.number_field_element import NumberFieldElement
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.structure.element import Matrix
from sage.structure.sage_object import SageObject
from .utils import Integer_t, length_from_M, \
    complex_tuple_to_json, complex_tuple_from_json, Real_t, dual_ideal, ideal_generator, \
    Complex_t, coefficient_dict_to_json, coefficient_dict_from_json

P = ParamSpec('P')
log = logging.getLogger(__name__)


class HilbertMaassCoefficients(SageObject):
    """
    Coefficients of Hilbert Maass Forms

    """
    def __init__(self, coefficients: Matrix,
                 M: tuple[Integer_t | tuple[Integer_t]],
                 spectral_parameter: tuple[Complex_t | Real_t],
                 space: 'HilbertMaassFormSpace', Y: tuple[Real_t] = None,
                 Q: tuple[Integer_t] = None,
                 coordinate_ideals: tuple[NumberFieldFractionalIdeal] = None,
                 set_coefficients: dict = None,
                 index_tuples: list[list] = None,
                 check: bool = True,
                 **kwargs: P.kwargs) -> None:
        r"""
        Init self.

        INPUT:

        - ``coefficients`` -- matrix of coefficients
        - ``M`` -- tuple of integers giving the limits used to map the indices of the matrix to tuples
        - ``spectral_parameter`` -- tuple of complex or real numbers
        - `` Q`` - tuple of positive intgers. >= M+2
        - ``space`` -- Hilbert Maass form space
        - ``coordinate_ideals`` -- list of fractional ideals indexing the different components
        - ``check`` -- boolean, if True, check that the coefficients are of the correct dimension.
        - ``set_coefficients`` -- dictionary with coefficients to set.
        - ``index_tuples`` -- list of lists of all tuples corresponding to indices in coefficients (in the same order)

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = matrix(RR, [[1],[2],[3],[4],[5],[6],[7],[8],[9]])
            sage: C = HilbertMaassCoefficients(Cmat, ((-1,1),(-1,1)), spectral_parameter, H); C
            Coefficients of a Hilbert Maass form with M=((-1, 1), (-1, 1)) and 1 cusp

        """
        super(HilbertMaassCoefficients, self).__init__(**kwargs)
        if not coordinate_ideals:
            representatives = space.group().ideal_cusp_representatives()
            coordinate_ideals = [dual_ideal(ideal) for ideal in representatives]
        if len(coordinate_ideals) != 1:
            raise NotImplementedError("Only one cusp supported for now")
        if check:
            if any(M0 == M1 == 0 for M0, M1 in M):
                raise ValueError('M must be non-zero')
            if coefficients.nrows() != length_from_M(M):
                raise ValueError(f'coefficients has wrong number of rows: {coefficients.nrows()}'
                                 f' != {length_from_M(M)}')
            if coefficients.ncols() != len(coordinate_ideals):
                raise ValueError('coefficients has number of cols: {coefficients.ncols()}'
                                 f' != {len(coordinate_ideals)}')
        from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
        if not coordinate_ideals and not space or not isinstance(space, HilbertMaassFormSpace):
            raise ValueError("Need to either specify coordinate ideals or a Maass form")

        self._coordinate_ideals = coordinate_ideals
        self._coefficients = coefficients
        self._coefficients.set_immutable()
        self._set_coefficients = set_coefficients
        self._M = M
        self._Q = Q
        self._spectral_parameter = spectral_parameter
        self._space = space
        self._Y = Y or tuple()
        if index_tuples:
            if not isinstance(index_tuples, list):
                raise ValueError('index_tuples must be a dict')
            if len(index_tuples) != len(self._coordinate_ideals):
                raise ValueError('index_tuples and coordinate_ideals have different lengths')
            if not isinstance(index_tuples[0], list):
                raise ValueError('index_tuples must be a list of lists')
        self._index_tuples = index_tuples

    def coordinate_ideals(self):
        return self._coordinate_ideals

    def coordinate_ideal_generators(self):
        return [ideal_generator(ida) for ida in self._coordinate_ideals]

    def to_json(self) -> dict:
        """
            JSON representation of self.

            EXAMPLES::

                sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
                sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
                sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
                sage: Y = (0.5, 0.5)
                sage: Q = (3, 3)
                sage: Cmat = matrix(RR, [[1],[2],[3],[4],[5],[6],[7],[8],[9]])
                sage: set_coefficients = {(0, 0): 0, (0, 1): 1}
                sage: C = HilbertMaassCoefficients(Cmat, ((-1,1),(-1,1)), spectral_parameter,
                ....: H,  Y, Q, set_coefficients=set_coefficients)
                sage: C.to_json()
                {'M': ((-1, 1), (-1, 1)),
                 'Q': (3, 3),
                 'Y': (0.5, 0.5),
                 'coefficients': [['1.00000000000000'],
                  ['2.00000000000000'],
                  ['3.00000000000000'],
                  ['4.00000000000000'],
                  ['5.00000000000000'],
                  ['6.00000000000000'],
                  ['7.00000000000000'],
                  ['8.00000000000000'],
                  ['9.00000000000000']],
                 'index_tuples': None,
                 'prec': 53,
                 'set_coefficients': {'[0, 0]': {'prec': 53, 'val': '0.000000000000000'},
                  '[0, 1]': {'prec': 53, 'val': '1.00000000000000'}},
                 'space': {'cuspidal': False,
                  'number_field': {'names': ['a'], 'polynomial': 'x^2 - 2'}},
                 'spectral_parameter': [{'prec': 53,
                   'val': '0.500000000000000 + 1.00000000000000*I'},
                  {'prec': 53, 'val': '0.500000000000000 + 1.00000000000000*I'}]}
        """
        return {
            'M': tuple((int(M0[0]), int(M0[1])) for M0 in self._M),
            'coefficients': [[str(x) for x in r] for r in self._coefficients],
            'prec': int(self._coefficients.base_ring().prec()),
            'spectral_parameter': complex_tuple_to_json(self._spectral_parameter),
            'set_coefficients': coefficient_dict_to_json(self._set_coefficients),
            'space': self._space.to_json(),
            'index_tuples': self._index_tuples,
            'Y': tuple(float(y) for y in self._Y),
            'Q': tuple(int(q) for q in self._Q)
        }

    @classmethod
    def from_json(cls, data: str | dict) -> 'HilbertMaassCoefficients':
        """
            Create an instance of HilbertMaassCoefficients from json formatted dict or string.

            EXAMPLES::

                sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
                sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
                sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
                sage: Y = (0.5, 0.5)
                sage: Q = (3, 3)
                sage: Cmat = matrix(RR, [[1],[2],[3],[4],[5],[6],[7],[8],[9]])
                sage: set_coefficients = {(0, 0): 0, (0, 1): 1}
                sage: C = HilbertMaassCoefficients(Cmat, ((-1,1),(-1,1)), spectral_parameter,
                ....: H,  Y, Q, set_coefficients=set_coefficients)
                sage: HilbertMaassCoefficients.from_json(C.to_json()) == C
                True
                sage: import json
                sage: json_string = json.dumps(C.to_json())
                sage: HilbertMaassCoefficients.from_json(json_string) == C
                True

        """
        if isinstance(data, str):
            data = json.loads(data)
        existing_keys = {'prec', 'coefficients', 'M', 'Y', 'Q', 'spectral_parameter', 'index_tuples',
            'space', 'set_coefficients'}
        if existing_keys.difference(data.keys()) not in [set({}), {'index_tuples'}]:
            raise ValueError("Not a valid JSON representation of HilbertMaassCoefficients")
        CF = ComplexField(data['prec'])
        coefficients = matrix(CF, data['coefficients'])
        from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
        M = tuple(tuple(x) for x in data['M'])
        Y = tuple(RealField(data['prec'])(x) for x in data['Y'])
        Q = tuple(x for x in data['Q'])
        return cls(
            coefficients, M=M,
            spectral_parameter=complex_tuple_from_json(data['spectral_parameter']),
            space=HilbertMaassFormSpace.from_json(data['space']),
            Y=Y,
            Q=Q,
            set_coefficients=coefficient_dict_from_json(data['set_coefficients']),
            index_tuples=data.get('index_tuples', [])
        )

    def __hash__(self):
        self._coefficients.set_immutable()
        return hash((self._coefficients, self._Y, self._Q, self._M, str(self._index_tuples),
                     self._spectral_parameter, self._space,
                     str(self._set_coefficients)))

    def M(self):
        return self._M

    def Y(self):
        return self._Y
    def Q(self):
        return self._Q
    def space(self):
        return self._space

    def coefficient_matrix(self):
        return self._coefficients

    def set_coefficients(self):
        return self._set_coefficients

    def spectral_parameter(self):
        return self._spectral_parameter

    def __eq__(self, other):
        if not isinstance(other, HilbertMaassCoefficients):
            return False
        return self.space() == other.space() and \
            self.coefficient_matrix() == other.coefficient_matrix() and \
            self.spectral_parameter() == other.spectral_parameter() and \
            self._Y == other._Y and \
            self._Q == other._Q and \
            self._M == other._M and \
            self._set_coefficients == other._set_coefficients

    def __repr__(self):
        """
        String representation of self.

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = matrix(RR, [[1],[2],[3],[4],[5],[6],[7],[8],[9]])
            sage: C = HilbertMaassCoefficients(Cmat, ((-1,1),(-1,1)), spectral_parameter, H)
            sage: C.__repr__()
            'Coefficients of a Hilbert Maass form with M=((-1, 1), (-1, 1)) and 1 cusp'


        """
        num_cusps = len(self._coordinate_ideals)
        return f"Coefficients of a Hilbert Maass form with M={self._M} and" \
               f" {num_cusps} cusp{'s' if num_cusps > 1 else ''}"

    def __getitem__(self, key: tuple | Integer_t) -> ComplexNumber:
        """
        Get one coefficient of self.

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassCoefficients
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: spectral_parameter = (CC(0.5,1),CC(0.5,1))
            sage: Cmat = matrix(RR, [[1],[2],[3],[4],[5],[6],[7],[8],[9]])
            sage: C = HilbertMaassCoefficients(Cmat, ((-1,1),(-1,1)), spectral_parameter, H)
            sage: C[0]
            1.00000000000000
            sage: C[(0, 1)]
            8.00000000000000
            sage: C[0, (0, -1)]
            2.00000000000000
            sage: C[0, (0, 1)]
            8.00000000000000
            sage: C[0, (1, 1)]
            9.00000000000000

            TESTS::

            sage: C[0,(0, 2)]
            Traceback (most recent call last):
            ...
            IndexError: Tuple element (0, 2) is out of bounds!
            sage: C[1,(0, 0)]
            Traceback (most recent call last):
            ...
            ValueError: Cusp index 1 out of bounds.
        """
        cusp = 0
        v = None
        if isinstance(key, tuple) and len(key) == 2 and isinstance(key[1], tuple):
            cusp, v = key
        elif isinstance(key, (int, Integer)) or isinstance(key, tuple):
            v = key
        elif isinstance(key, NumberFieldElement):
            cusp = 0
            v = key
        if cusp < 0 or cusp >= len(self._coordinate_ideals):
            raise ValueError(f"Cusp index {cusp} out of bounds.")
        if isinstance(v, NumberFieldElement):
            v = ideal_coordinates(self._coordinate_ideals[cusp], v)
        if isinstance(v, tuple):
            if v in self.keys():
                v = self.keys().index(v)
            else:
                raise IndexError(f"Tuple element {v} is out of bounds!")
        try:
            return self._coefficients.column(cusp)[v]
        except IndexError:
            raise ValueError(f'Can not get coefficient for {cusp}:{v}')

    def get_coefficient_index(self, v: tuple, cusp: int = 0) -> int:
        if self._index_tuples:
            return self._index_tuples[cusp].index(v)
        return map_tuple_to_int(v, self._M)

    # @cached_method
    def keys(self, as_elements=False) -> list:
        if self._index_tuples:
            keys_tuple = self._index_tuples[0]
        else:
            keys_tuple = [map_int_to_tuple(v, self._M)
                      for v in range(self._coefficients.nrows())]
            self._index_tuples = [keys_tuple]
        if not as_elements:
            return keys_tuple
        else:
            basis = self._coordinate_ideals[0].integral_basis()
            return [sum(t[i] * basis[i] for i in range(len(basis))) for t in keys_tuple]

    def __iter__(self) -> iter:
        return iter(self._coefficients.column(0))

    def norms(self, use_abs=False):
        r"""
        Add up coefficients of the same norm.
        This is mainly a way to see if a form is "close to 0".
        """
        elements_by_norm = {}
        for k in self.keys(as_elements=True):
            if k.norm() not in elements_by_norm:
                elements_by_norm[k.norm()] = 0
            elements_by_norm[k.norm()] += abs(self[k]) if use_abs else self[k]
        return elements_by_norm