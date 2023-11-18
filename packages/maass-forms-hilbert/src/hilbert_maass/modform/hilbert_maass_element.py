"""
Classes For Hilbert-Maass forms

"""
import json
import logging
from typing import ParamSpec

# from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
from hilbert_modgroup.pullback import HilbertPullback
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.rings.real_mpfr import RealNumber as RealNumber_class
from sage.structure.element import ModuleElement, Matrix

from .coefficients import HilbertMaassCoefficients, compute_coefficients
from hilbert_maass.modform.utils import Integer_t, complex_tuple_to_json

P = ParamSpec('P')
log = logging.getLogger(__name__)


class HilbertMaassForm_Element(ModuleElement):

    def __init__(self, parent: 'HilbertMaassFormSpace',
                 spectral_parameter: tuple[ComplexNumber | RealNumber_class],
                 coefficients: Matrix | HilbertMaassCoefficients = None,
                 **kwargs: P.kwargs) -> None:
        super(HilbertMaassForm_Element, self).__init__(parent, **kwargs)
        self.cuspidal = parent.is_cuspidal()
        self.spectral_parameter = spectral_parameter
        self._number_field = parent.number_field()
        self._complex_field = self.spectral_parameter[0].parent().prec()
        self.has_coefficients = False
        if coefficients is None:
            self._coefficients = None
            # {c: {} for c in range(self.parent().group().ncusps())}
        elif isinstance(coefficients, HilbertMaassCoefficients):
            self._coefficients = coefficients
        else:
            self._coefficients = HilbertMaassCoefficients(coefficients, self.parent().group())
        self._pullback = HilbertPullback(self.parent().group())
        different = self.parent().number_field().different()
        representatives = self.parent().group().ideal_cusp_representatives()
        self._dual_ideals = [ideal**-1*different**-1 for ideal in representatives]
        self._dual_ideal_matrix = {
            ideal:  matrix(self._pullback._get_lattice_and_ideal_basis(ideal**-1*different**-1)[0])
            for ideal in representatives
        }

    def __reduce__(self):
        return self.__class__, (self.parent(), self.spectral_parameter, self._coefficients)

    def to_json(self):
        return {
            'parent': self.parent().to_json(),
            'spectral_parameter': complex_tuple_to_json(self.spectral_parameter),
            'coefficients':
                self.coefficients().to_json()
        }

    @classmethod
    def from_json(cls, data):
        from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
        if isinstance(data, str):
            data = json.loads(data)
        parent = HilbertMaassFormSpace.from_json(data=data['parent'])
        spectral_parameter = tuple(ComplexField(x['prec'])(x['val'])
                              for x in data['spectral_parameter'])
        coefficients = HilbertMaassCoefficients.from_json(data['coefficients'])
        return cls(parent, spectral_parameter, coefficients)

    def is_cuspidal(self):
        return self.cuspidal

    def __eq__(self, other):
        if not isinstance(other, HilbertMaassForm_Element):
            return False
        return self.parent() == other.parent() and \
            self.spectral_parameter() == other.spectral_parameter() and \
            self.coefficients() == other.coefficients()

    def __repr__(self):
        return f"Hilbert Maass form for {self.parent()} with spectral parameter" \
               f" {self.spectral_parameter}"

    def dual_ideal_element(self, coordinates: tuple[Integer_t] or vector,
                           ideal: NumberFieldFractionalIdeal):
        return self._dual_ideal_matrix[ideal]*vector(coordinates)

    def spectral_parameter(self):
        return self._spectral_parameter

    def pullback(self):
        return self._pullback

    def coefficients(self):
        return self._coefficients

    def compute_coefficients(self, s: tuple = None,
                             ideala: NumberFieldFractionalIdeal = None,
                             idealb: NumberFieldFractionalIdeal = None,
                             M: tuple[tuple[Integer_t]] = None,
                             Y: tuple = None,
                             prec: int = 53,
                             sgn: str = '-',
                             returnV: bool = False) -> 'HilbertMaassCoefficients' or tuple:
        r"""

        INPUT:

        - ``ideala``  -- NumberField Fractional Ideal corresponding to cusp.
        - ``idealb``  --
        - ``s``       --
        - ``M``  --
        - ``Y``  --
        - ``prec``  --
        - ``cuspidal``  --
        - ``sgn``  --
        - ``returnV``  --

        EXAMPLES::

            sage: from hilbert_maass.all import HilbertMaassForm
            sage: M = (2,2)
            sage: spectral_parameter = (CC(1.5,1.5),)*2
            sage: F = HilbertMaassForm(QuadraticField(2), cuspidal=False)
            Traceback (most recent call last):
            ...
            TypeError: HilbertMaassForm() missing 1 required positional argument: 'spectral...
            sage: F = HilbertMaassForm(QuadraticField(2), spectral_parameter, cuspidal=False)
            sage: C = F.compute_coefficients(spectral_parameter, M = (-1,1))
            Cannot connect to cache collection.
            sage: C[(0,0)] # abs tol 1e-10
            0.167572243136260 - 0.590945405474222*I
            sage: F = HilbertMaassForm(QuadraticField(2), spectral_parameter, cuspidal=False)
            sage: C = F.compute_coefficients(spectral_parameter, M = (-3,3)) # long time (100s)
            sage: C[(0,0)] # abs tol 1e-10 # long time (100s)
            0.245942691776149 - 0.595428499664962*I

        """
        s = s or self.spectral_parameter
        if not s:
            raise ValueError("Spectral parameter must be set in the HilbertMaassForm or "
                             "passed as parameter")
        C = compute_coefficients(space=self.parent(),
                                 spectral_parameter=s,
                                 ideala=ideala,
                                 idealb=idealb,
                                 M=M,
                                 Y=Y,
                                 )
        self._coefficients = C
        return C


def HilbertMaassForm(group: 'HilbertModularGroup' or 'HilbertMaassFormSpace',
                     spectral_parameter: tuple[ComplexNumber | RealNumber_class],
                     **kwargs: P.kwargs) -> HilbertMaassForm_Element:
    from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
    if isinstance(group, HilbertMaassFormSpace):
        space = group
    else:
        space = HilbertMaassFormSpace(group, **kwargs)
    return HilbertMaassForm_Element(space, spectral_parameter=spectral_parameter)
