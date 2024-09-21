"""
Classes For Hilbert-Maass forms

"""
import json
import logging
from copy import copy, deepcopy
from typing import ParamSpec

from hilbert_maass.functions.functions import bessel_prod
from hilbert_maass.functions.functions_cy import exp_trace_prod_dp
from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement__class
from matplotlib import pyplot as plt
from sage.all import CC
from sage.arith.srange import xsrange
# from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
from sage.functions.other import real, imag
from sage.matrix.constructor import matrix
from sage.plot.animate import animate, Animation
from sage.plot.misc import setup_for_eval_on_grid
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.rings.number_field.number_field_base import NumberField as NumberField_class
from sage.rings.real_mpfr import RealNumber as RealNumber_class
from sage.structure.element import ModuleElement, Matrix

from .coefficients import HilbertMaassCoefficients
from .compute_coefficients import compute_coefficients
from hilbert_maass.modform.utils import Integer_t, complex_tuple_to_json, cartesian_product_from_M, \
    Complex_t, Real_t, ideal_coordinates

P = ParamSpec('P')
log = logging.getLogger(__name__)


class HilbertMaassForm_Element(ModuleElement):

    def __init__(self, parent: 'HilbertMaassFormSpace',
                 spectral_parameter: tuple[ComplexNumber | RealNumber_class],
                 coefficients: Matrix | HilbertMaassCoefficients = None,
                 **kwargs: P.kwargs) -> None:
        super(HilbertMaassForm_Element, self).__init__(parent, **kwargs)
        self.cuspidal = parent.is_cuspidal()
        self._spectral_parameter = spectral_parameter
        self._number_field = parent.number_field()
        if hasattr(self._spectral_parameter[0], 'parent') and \
                hasattr(self._spectral_parameter[0].parent(), 'prec'):
            self._complex_field = self._spectral_parameter[0].parent()
        else:
            self._complex_field = ComplexField(prec=53)
        self.has_coefficients = False
        if coefficients is None:
            self._coefficients = None
            # {c: {} for c in range(self.parent().group().ncusps())}
        elif isinstance(coefficients, HilbertMaassCoefficients):
            self._coefficients = coefficients
        else:
            self._coefficients = HilbertMaassCoefficients(coefficients, self.parent().group())

    def __reduce__(self):
        return self.__class__, (self.parent(), self.spectral_parameter(), self._coefficients)

    def __hash__(self):
        return hash((self.parent(),
                    self.spectral_parameter(),
                    self.coefficients()))

    def __repr__(self):
        return f"Hilbert Maass form for {self.parent()} with spectral parameter" \
               f" {self.spectral_parameter()}"

    def to_json(self):
        """
        Json representation of self.

        EXAMPLES:

            sage: from hilbert_maass.all import HilbertMaassFormSpace, HilbertMaassForm_Element
            sage: H = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
            sage: element = HilbertMaassForm_Element(H, (CC(0,1),CC(0,1)))
            sage: json_data = element.to_json()
            sage: json_data
            {'coefficients': {},
             'parent': {'cuspidal': False,
             'number_field': {'names': ['a'], 'polynomial': 'x^2 - 2'}},
             'spectral_parameter': [{'prec': 53, 'val': '1.00000000000000*I'},
             {'prec': 53, 'val': '1.00000000000000*I'}]}

            Convert the JSON representation back to an instance of `HilbertMaassElement`:

            sage: new_element = HilbertMaassForm_Element.from_json(json_data)
            sage: new_element == element
            True
        """
        return {
            'parent': self.parent().to_json(),
            'spectral_parameter': complex_tuple_to_json(self.spectral_parameter()),
            'coefficients': self.coefficients().to_json() if self.coefficients() else {}
        }

    @classmethod
    def from_json(cls, data):
        from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
        if isinstance(data, str):
            data = json.loads(data)
        parent = HilbertMaassFormSpace.from_json(data=data['parent'])
        spectral_parameter = tuple(ComplexField(x['prec'])(x['val'])
                                   for x in data['spectral_parameter'])
        if not data['coefficients']:
            return cls(parent, spectral_parameter)
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

    def __call__(self, z: list | tuple, **kwargs: P.kwargs) -> ComplexNumber:
        if self._coefficients is None:
            raise ValueError("Coefficients must be computed first")
        if self.parent().group().ncusps() > 1:
            # TODO: support multiple cusps: find closest cusp and use correct Fourier expansion
            raise NotImplementedError("Only one cusp supported for now")
        C = self._coefficients
        if isinstance(z, UpperHalfPlaneProductElement__class):
            x = z.real()
            y = z.imag()
        else:
            x = [real(zi) for zi in z]
            y = [imag(zi) for zi in z]
            if not all([yi > 0 for yi in y]):
                raise ValueError("y must be positive")
        summa = 0
        ideala = 0
        for V in cartesian_product_from_M(self.coefficients().M()):
            V = tuple(V)
            v = self.parent().dual_ideal_element(V, ideala)
            bes = bessel_prod(v, tuple(y), self.spectral_parameter())
            exp_arg = tuple([x[i] * v for i, v in enumerate(v)])
            term = bes * exp_trace_prod_dp(exp_arg)
            summa += self.coefficients()[V] * term
        return summa

    def spectral_parameter(self):
        return self._spectral_parameter

    def pullback(self):
        return self._pullback

    def coefficients(self):
        return self._coefficients

    def __mul__(self, other):
        """
        """
        if not isinstance(other, (Real_t, Complex_t, Integer_t)):
            raise ValueError("Multiplication is only defined for real or complex numbers")
        result = copy(self)
        result._coefficients._coefficients = other * result._coefficients._coefficients
        return result

    def __copy__(self):
        coefficients = copy(self._coefficients)
        return self.__class__(self.parent(), self._spectral_parameter, coefficients)

    def _lmul_(self, other):
        """
        """
        if not isinstance(other, (Real_t, Complex_t, Integer_t)):
            raise ValueError("Multiplication is only defined for real or complex numbers")
        result = copy(self)
        result._coefficients._coefficients = other * result._coefficients._coefficients
        return result

    def _add_(self, other):
        """
        """
        if not isinstance(other, self.__class__):
            raise ValueError("Addition is only defined for HilbertMaassForms_Elements objects")
        if other.parent() != self.parent():
            raise ValueError("Addition is only defined for HilbertMaassForms_Elements objects "
                             "with the same parent")
        # Note that the sum will be supported on the intersection
        # of the indices of the individual forms.
        result = copy(self)
        # Note that indices may differ...
        used_indices = []
        coefficients = []
        for k, v in dict(self._coefficients).items():
            if k in dict(other._coefficients):
                coefficients.append((v + other._coefficients[k],))
                used_indices.append(k)

        result._coefficients._coefficients = matrix(coefficients)
        result._coefficients._index_tuples = [used_indices]
        return result

    def _sub_(self, other):
        """
        """
        return self + other * -1

    def galois_conjugate(self, i):
        r"""
        Return Galois conjugate no. i of self.
        """
        if not isinstance(i, Integer_t):
            raise ValueError("Conjugate no. must be an integer")
        indices_used = []
        galois_group = self._number_field.galois_group()
        galois_map = galois_group[i]
        coordinate_ideal = self.coefficients().coordinate_ideals()[0]
        coefficients = []
        for k in self._coefficients.keys(as_elements=True):
            mapped_index = galois_map(k)
            coordinates_mapped = ideal_coordinates(coordinate_ideal, mapped_index)
            indices_used.append(coordinates_mapped)
            coefficients.append((self._coefficients[k],))

        coeffs = HilbertMaassCoefficients(matrix(coefficients),
                                          M=self.coefficients().M(),
                                          Y=self.coefficients().Y(),
                                          spectral_parameter=self.spectral_parameter(),
                                          space=self.parent(),
                                          coordinate_ideals=self.coefficients().coordinate_ideals(),
                                          set_coefficients=self.coefficients().set_coefficients(),
                                          index_tuples=[indices_used],
                                          check=False)
        return HilbertMaassForm(self.parent(),
                                self.spectral_parameter(),
                                coefficients=coeffs)


    def action_by_unit(self, u):
        r"""
        Act on self by unit u through action on the coefficients.
        """
        if u not in self.parent().number_field():
            raise ValueError("Unit must be in the number field")
        if u not in self.parent().number_field().unit_group() and not u.is_unit():
            raise ValueError("Unit must be a unit")
        indices_used = []
        coordinate_ideal = self.coefficients().coordinate_ideals()[0]
        coefficients = []
        for k in self._coefficients.keys(as_elements=True):
            mapped_index = k * u
            coordinates_mapped = ideal_coordinates(coordinate_ideal, mapped_index)
            indices_used.append(coordinates_mapped)
            coefficients.append((self._coefficients[k],))

        coeffs = HilbertMaassCoefficients(matrix(coefficients),
                                          M=self.coefficients().M(),
                                          Y=self.coefficients().Y(),
                                          spectral_parameter=self.spectral_parameter(),
                                          space=self.parent(),
                                          coordinate_ideals=self.coefficients().coordinate_ideals(),
                                          set_coefficients=self.coefficients().set_coefficients(),
                                          index_tuples=[indices_used],
                                          check=False)
        return HilbertMaassForm(self.parent(),
                                self.spectral_parameter(),
                                coefficients=coeffs)


    def compute_coefficients(self, s: tuple = None,
                             ideala: NumberFieldFractionalIdeal = None,
                             idealb: NumberFieldFractionalIdeal = None,
                             M: tuple[tuple[Integer_t]] = None,
                             Y: tuple = None,
                             Q: tuple = None,
                             set_coefficients: dict = None,
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
            sage: C = F.compute_coefficients(spectral_parameter, M = (-1,1), Q=(10,10))
            sage: C[(0,0)] # abs tol 1e-10
            0.245523867043680 - 0.593475166148685*I
            sage: F = HilbertMaassForm(QuadraticField(2), spectral_parameter, cuspidal=False)
            sage: C = F.compute_coefficients(spectral_parameter, M = (-3,3)) # long time (100s)
            sage: C[(0,0)] # abs tol 1e-10 # long time (100s)
            0.245942691776149 - 0.595428499664962*I

        """
        s = s or self.spectral_parameter()
        if not s:
            raise ValueError("Spectral parameter must be set in the HilbertMaassForm or "
                             "passed as parameter")
        C = compute_coefficients(space=self.parent(),
                                 spectral_parameter=s,
                                 ideala=ideala,
                                 idealb=idealb,
                                 M=M,
                                 Y=Y,
                                 Q=Q,
                                 set_coefficients=set_coefficients
                                 )
        self._coefficients = C
        return C

    def animation(self, num_steps: Integer_t = 100,
                  y_start: Real_t = 0, y_stop: Real_t = 1,
                  x_start: Real_t = 0, x_stop: Real_t = 0,
                  **kwargs: P.kwargs) -> Animation:
        """
        Create an animation of the Hilbert Maass form.

        INPUT:


        - kwargs:
        """
        num_steps = kwargs.get('num_steps', 100)
        h = (y_stop - y_start) / num_steps
        glist = []
        for i in range(num_steps):
            g = self.plot(yset=[y_start + h * (i+1)], **kwargs)
            g.save_image = g.savefig
            glist.append(g)
        return animate(glist)

    def plot(self, **kwargs):
        """
        Density plot of self along one copy of the hyperbolic upper half-plane with
        other parameters set to fixed values (by default set to i).

        Note: You need to compute Fourier coefficients before plotting.

        PLOT OPTIONS:

        - ``plot_points`` -- (default: `200`); the minimal number of plot points.

        - ``xmin`` -- starting x value.
        - ``xmax`` -- ending x value.
        - ``ymin`` -- starting y value.
        - ``ymax`` -- ending y value.
        - ``xset`` -- list of fixed x values (default: 0).
        - ``yset`` -- list of fixed y values (default: 1).
        - ``cmap`` -- color map (default: `jet`).

        EXAMPLES::

            sage: from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
            sage: F = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False).an_element()
            sage: F.plot()
            Traceback (most recent call last):
            ...
            ValueError: Coefficients must be computed first
            sage: F.compute_coefficients(M=1)
            Coefficients of a Hilbert Maass form with M=((-1, 1), (-1, 1)) and 1 cusp
            sage: F.plot()
            <Figure size 800x399.99 with 1 Axes>
        """
        n = self.parent().number_field().degree()
        xset = kwargs.get("xset", [0] * (n - 1))
        yset = kwargs.get("yset", [1] * (n - 1))
        xmin = kwargs.get("xmin", -4)
        xmax = kwargs.get("xmax", 4)
        ymin = kwargs.get("ymin", 0.0001)
        ymax = kwargs.get("ymax", 4)
        show_axis = kwargs.get("show_axis", False)
        plot_points_x = kwargs.get("plot_points_x", 50)
        plot_points_y = kwargs.get("plot_points_y", 50)
        cmap = kwargs.get('cmap', ['jet'])
        # Create grid points
        fixed_zs = [CC(x, y) for x, y in zip(xset, yset)]

        def function_to_eval(x, y):
            return abs(self(fixed_zs + [CC(x, y)]))

        g, ranges = setup_for_eval_on_grid([function_to_eval],
                                           [[xmin, xmax], [ymin, ymax]],
                                           [plot_points_x, plot_points_y])
        g = g[0]
        xy_data_array = [[g(x, y) for x in xsrange(*ranges[0], include_endpoint=True)] for y in
                         xsrange(*ranges[1], include_endpoint=True)]
        res = []
        for cmapi in cmap:
            g = plt.figure(figsize=(xmax - xmin, ymax - ymin))
            ax = g.add_subplot(111)
            t = ax.imshow(xy_data_array, origin='lower',
                      cmap=cmapi,
                      extent=(xmin, xmax, ymin, ymax),
                      interpolation='catrom')
            if not show_axis:
                ax.set_frame_on(False)
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
            res.append(g)
        if len(res) == 1:
            return res[0]
        return res


def HilbertMaassForm(group: 'HilbertModularGroup' or 'HilbertMaassFormSpace' or NumberField_class,
                     spectral_parameter: tuple[ComplexNumber | RealNumber_class],
                     **kwargs: P.kwargs) -> HilbertMaassForm_Element:
    """
    Create a Hilbert Maass form

    INPUT:

    - ``group``  -- Hilbert modular group or space of Hilbert maass forms
    - ``spectral_parameter`` -- tuple of complex numbers


    EXAMPLES::

        sage: from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
        sage: space = HilbertMaassFormSpace(QuadraticField(2), cuspidal=False)
        sage: from hilbert_maass.modform.hilbert_maass_element import HilbertMaassForm
        sage: spectral_parameter = (0.5 + 0.5j, 0.5 + 1j)
        sage: HilbertMaassForm(space, (0.5, 1.5))
        Hilbert Maass form for HilbertMaassFormSpace(Hilbert Modular Group PSL(2) over Maximal Order...
        sage: HilbertMaassForm(QuadraticField(2), spectral_parameter)
        Hilbert Maass form for HilbertMaassFormSpace(Hilbert Modular Group PSL(2) over Maximal Order...

    """
    from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
    coefficients = kwargs.pop('coefficients', None)
    if isinstance(group, HilbertMaassFormSpace):
        space = group
    else:
        space = HilbertMaassFormSpace(group, **kwargs)
    return HilbertMaassForm_Element(space, spectral_parameter=spectral_parameter,
                                    coefficients=coefficients)

