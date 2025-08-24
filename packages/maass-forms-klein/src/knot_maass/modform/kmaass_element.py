"""
Elements of spaces of Maass waveforms for Kleinian complements.

"""
import json
from typing import ParamSpec
import mongoengine as me
from matplotlib import pyplot as plt
from sage.all import CC
from sage.arith.srange import xsrange
from sage.functions.other import imag, real
from sage.misc.cachefunc import cached_method
from sage.plot.animate import Animation, animate
from sage.plot.misc import setup_for_eval_on_grid
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.real_mpfr import RealField
from sage.structure.parent import Parent
from sage.structure.element import Element, Matrix, Vector
from hilbert_maass.modform.utils import Real_t, Complex_t
from knot_maass.hyperbolic_space.upper_half_space import UpperHalfSpaceElement__class
from knot_maass.modform.coefficients import KleinianMaassFormCoefficients
from knot_maass.modform.utils import bessel_function
from knot_maass.hyperbolic_space.upper_half_space import UpperHalfSpaceElement

from knot_maass.modform.utils import Integer_t

P = ParamSpec('P')


class KleinianMaassFormElement(Element):
    r"""
    Element of a KleinianMaassFormSpace.

    TODO: Decide whether to keep this as a Python class or make a cdef Cython class.
    """

    coefficients = []

    def __init__(self,  parent: Parent, spectral_parameter, coefficients=None,
                 **kwargs: P.kwargs) -> None:
        r"""

        INPUT:

        - `parent` -- parent space of type KleinianMaassFormSpace
        - `spectral_parameter` -- spectral parameter
        - `coefficients` -- coefficients
        - `kwargs` -- keyword arguments

        EXAMPLES::

        sage: from knot_maass.all import KleinianMaassFormSpace, KleinianMaassFormElement
        sage: S = KleinianMaassFormSpace(-4)
        sage: F = KleinianMaassFormElement(S, CC(0.5, 1.0))
        sage: #TestSuite(F).run()

        """
        super(KleinianMaassFormElement, self).__init__(parent, **kwargs)
        self.cuspidal = parent.is_cuspidal()
        self._spectral_parameter = spectral_parameter
        if coefficients and not isinstance(coefficients, KleinianMaassFormCoefficients):
            raise ValueError("Coefficients must be of type KleinianMaassFormCoefficients.")
        elif not coefficients:
            coefficients = KleinianMaassFormCoefficients([], 0,
                                                         self._spectral_parameter,
                                                         self.parent())
        self._coefficients = coefficients

    def __repr__(self):
        return f"KleinianMaassFormElement({self.parent()}, {self.spectral_parameter()})"

    def spectral_parameter(self) -> Complex_t:
        return self._spectral_parameter

    def __call__(self, z: list | tuple, **kwargs: P.kwargs) -> ComplexNumber:
        if self._coefficients is None:
            raise ValueError("Coefficients must be computed first")
        C = self._coefficients
        if not isinstance(z, UpperHalfSpaceElement__class):
            z = UpperHalfSpaceElement(z)
        x = z.z()
        y = z.y()
        complex_values = self._complex_value_vector(x)
        bessel_values = self._bessel_value_vector(y)
        summa = 0
        for n, v in enumerate(C.coordinate_values()):
            summa += C[v] * complex_values[n] * bessel_values[n]
        return summa

    @cached_method
    def _bessel_value_vector(self, y):
        twopi = RealField(53).pi() * 2
        return [y * bessel_function(abs(v), twopi * y, self.spectral_parameter()) for
                v in self.coefficients().coordinate_values()]
    @cached_method
    def _complex_value_vector(self, x):
        CF = ComplexField(53)
        twopii = CF(0, RealField(53).pi() * 2)
        return [(twopii * (v[0] * x[0] - v[1] * x[1])).exp() for
                v in self.coefficients().coordinate_values()]


    def coefficients(self) -> dict:
        """
        Return the coefficients of self.

        EXAMPLES::

            sage: from knot_maass.all import KleinianMaassFormSpace, KleinianMaassFormElement
            sage: S = KleinianMaassFormSpace(-4)
            sage: F = S.an_element()
            sage: F.coefficients()
            Coefficients of a Kleinian Maass form with M=0
        """
        return self._coefficients

    def compute_coefficients(self, spectral_parameter: Real_t | ComplexNumber = None,
                             **kwargs: P.kwargs) -> None:
        """
        Compute the Fourier coefficients of self.

        EXAMPLES::

            sage: from knot_maass.all import KleinianMaassFormSpace, KleinianMaassFormElement
            sage: S = KleinianMaassFormSpace(-4)
            sage: F = S.an_element()

        """
        from knot_maass.modform.compute_coefficients import compute_coefficients
        if not spectral_parameter:
            spectral_parameter = self._spectral_parameter
        self._coefficients = compute_coefficients(self.parent(), spectral_parameter, **kwargs)

    def setup_matrix(self, *args: P.args, **kwargs: P.kwargs) -> Matrix:
        r"""
        Set up the matrix for the linear system.

        EXAMPLES::

            sage: from knot_maass.all import KleinianMaassFormSpace, KleinianMaassFormElement
            sage: S = KleinianMaassFormSpace(-4)
            sage: F = S.an_element()
            sage: V = F.setup_matrix()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def normalise_matrix(self, *args: P.args, **kwargs: P.kwargs) -> tuple:
        r"""
        Normalise the matrix, e.g. set first coefficient to 1 etc.

        EXAMPLES::

            sage: from knot_maass.all import KleinianMaassFormSpace, KleinianMaassFormElement
            sage: S = KleinianMaassFormSpace(-4)
            sage: F = S.an_element()
            sage: V = F.setup_matrix()
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: V = F.normalise_matrix(None)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def solve_system(self, *args: P.args, **kwargs: P.kwargs) -> Vector:
        r"""
        Solve the system.

        EXAMPLES::

            sage: from knot_maass.all import KleinianMaassFormSpace, KleinianMaassFormElement
            sage: S = KleinianMaassFormSpace(-4)
            sage: F = S.an_element()
            sage: V = F.setup_matrix()
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: V, B = F.normalise_matrix(None)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: X = F.solve_system(None, None)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

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
            'spectral_parameter': {'prec': self.spectral_parameter().prec(),
                                   'val': str(self.spectral_parameter())},
            'coefficient_keys': self.coefficient_keys(),
            'coefficient': {
                'keys': self.coefficients().keys(),
                'values': self.coefficients().values()
            },
            'r_value': imag(self._spectral_parameter),
        }

    @classmethod
    def from_json(cls, data):
        if isinstance(data, str):
            data = json.loads(data)
        from knot_maass.modform.kmaass_space import KleinianMaassFormSpace
        parent = KleinianMaassFormSpace.from_json(data=data['parent'])
        s = data['spectral_parameter']
        spectral_parameter = ComplexField(s['prec'])(s['val'])
        if not data['coefficients']:
            return cls(parent, spectral_parameter)
        coefficients = dict(zip(data['coefficient_keys'], data['coefficient']['values']))
        return cls(parent, spectral_parameter, coefficients)

    def to_json(self):
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

    def animation(self, num_steps: Integer_t = 100,
                  y_start: Real_t = 0, y_stop: Real_t = 1,
                  x_start: Real_t = 0, x_stop: Real_t = 0,
                  **kwargs: P.kwargs) -> Animation:
        """
        Create an animation of the Hilbert Maass form.

        INPUT:


        - kwargs:
        """
        y_is_set = x_is_set = False
        if y_start != y_stop:
            h = (y_stop - y_start) / num_steps
            x_is_set = True
        elif x_start != x_stop:
            h = (x_stop - x_start) / num_steps
            y_is_set = True
        else:
            raise ValueError("Either y_start or x_start must be different from y_stop or x_stop")
        glist = []
        for i in range(num_steps):
            if y_is_set:
                g = self.plot(y0set=[y_start + h * i], **kwargs)
            else:
                g = self.plot(x0set=[x_start + h * i], **kwargs)
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

        # We fix either x0 or x1 in the complex plane
        x0set = kwargs.get("x0set", 0)
        x1set = kwargs.get("x1set", 0)
        xmin = kwargs.get("xmin", -4)
        xmax = kwargs.get("xmax", 4)
        ymin = kwargs.get("ymin", 0.0001)
        ymax = kwargs.get("ymax", 4)
        show_axis = kwargs.get("show_axis", False)
        plot_points_x = kwargs.get("plot_points_x", 50)
        plot_points_y = kwargs.get("plot_points_y", 50)
        cmap = kwargs.get('cmap', ['jet'])
        # Create grid points

        def function_to_eval(x, y):
            if 'x1set' in kwargs:
                return abs(self([x, x1set, y]))
            else:
                return abs(self([x0set, x, y]))

        g, ranges = setup_for_eval_on_grid([function_to_eval],
                                           [[xmin, xmax], [ymin, ymax]],
                                           [plot_points_x, plot_points_y])
        g = g[0]
        xy_data_array = [[g(x, y) for x in xsrange(*ranges[0], include_endpoint=True)] for y in
                         xsrange(*ranges[1], include_endpoint=True)]
        res = []
        for cmapi in cmap:
            g = plt.figure(figsize=(5*(xmax - xmin), 5*(ymax - ymin)))
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


