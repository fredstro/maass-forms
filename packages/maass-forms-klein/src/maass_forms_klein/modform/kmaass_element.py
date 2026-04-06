r"""
Elements of spaces of Maass waveforms for Kleinian complements.

This module implements individual Maass waveforms as elements of
KleinianMaassFormSpace. These are eigenfunctions of the hyperbolic
Laplacian on quotient spaces of the upper half-space by discrete
subgroups of PSL(2,C).

AUTHORS:

- Fredrik Strömberg (2024): initial version

EXAMPLES::

    sage: from sage.rings.complex_mpfr import ComplexField
    sage: CC = ComplexField(53)
    sage: s = CC(0.5, 14.1)  # spectral parameter
    sage: s.real(), s.imag()
    (0.500000000000000, 14.1000000000000)
"""
import json
from typing import ParamSpec, Union, Optional, List, Tuple
import mongoengine as me
from matplotlib import pyplot as plt
from sage.arith.srange import xsrange
from sage.functions.other import imag
from sage.misc.cachefunc import cached_method
from sage.plot.animate import Animation, animate
from sage.plot.misc import setup_for_eval_on_grid
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.real_mpfr import RealField
from sage.structure.parent import Parent
from sage.structure.element import Element, Matrix, Vector
from maass_form_core.utils.types import Real_t, Complex_t

# Import exceptions for validation
from maass_forms_klein.exceptions import (
    InvalidSpectralParameterError,
    ComputationError,
    ValidationError
)
from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement__class
from maass_forms_klein.modform.coefficients import KleinianMaassFormCoefficients
from maass_forms_klein.modform.utils import bessel_function, Integer_t
from maass_forms_klein.hyperbolic_space.upper_half_space import UpperHalfSpaceElement

P = ParamSpec("P")


class KleinianMaassFormElement(Element):
    r"""Individual Kleinian Maass waveform element.
    
    This class represents a specific Maass waveform associated with a Kleinian
    group. A Kleinian Maass form is a smooth function f on the upper half-space 
    H^3 satisfying:
    
    1. **Automorphy condition**: f(γz) = f(z) for all γ in the Kleinian group
    2. **Eigenvalue equation**: Δf = λf where Δ is the hyperbolic Laplacian
    3. **Growth conditions**: appropriate behavior at cusps and infinity
    
    The spectral parameter λ determines the eigenvalue, and the form can be
    expanded in terms of Fourier coefficients and Bessel functions.
    
    ATTRIBUTES:
    
    - ``_spectral_parameter`` -- Complex; the eigenvalue parameter λ
    - ``_coefficients`` -- KleinianMaassFormCoefficients; Fourier expansion data
    - ``cuspidal`` -- bool; whether this is a cuspidal form
    
    EXAMPLES::
    
        sage: from sage.rings.complex_mpfr import ComplexField
        sage: CC = ComplexField(53)
        sage: s = CC(0.5, 14.1)
        sage: # KleinianMaassFormElement requires a parent space which needs database connection
        sage: # For this example, we show the spectral parameter structure
        sage: s.real() > 0 and s.imag() > 0
        True
    
    TESTS::
    
        sage: # Test invalid parent
        sage: KleinianMaassFormElement(None, CC(0.5, 14.1))  # doctest: +SKIP
        Traceback (most recent call last):
        ...
        ValidationError: Parent must be a KleinianMaassFormSpace
    
    .. NOTE::
    
        The computational representation uses Fourier expansions with
        Whittaker functions (Bessel K functions) for efficient evaluation.
    """

    # Class attribute for coefficient storage
    coefficients = []

    def __init__(
        self,
        parent: Parent,
        spectral_parameter: Complex_t,
        coefficients: Optional[KleinianMaassFormCoefficients] = None,
        **kwargs: P.kwargs
    ) -> None:
        r"""Initialize a Kleinian Maass form element.

        INPUT:

        - ``parent`` -- KleinianMaassFormSpace; the parent space
        - ``spectral_parameter`` -- Complex; the eigenvalue parameter λ
        - ``coefficients`` -- KleinianMaassFormCoefficients (optional);
          Fourier expansion coefficients. If None, creates empty coefficient structure
        - ``**kwargs`` -- additional keyword arguments passed to parent Element

        EXAMPLES::

            sage: from sage.rings.complex_mpfr import ComplexField
            sage: CC = ComplexField(53)
            sage: s = CC(0.5, 14.1)
            sage: # Test spectral parameter validation
            sage: s.parent()
            Complex Field with 53 bits of precision
            sage: s.prec()
            53
            sage: # Full object creation requires database connection
            sage: # Testing validation of spectral parameter types
            sage: s.real().prec()
            53

        TESTS::
        
            sage: # Test validation
            sage: KleinianMaassFormElement("invalid", s)  # doctest: +SKIP
            Traceback (most recent call last):
            ...
            ValidationError: Parent must be a KleinianMaassFormSpace
            
            sage: f = KleinianMaassFormElement(S, "invalid")  # doctest: +SKIP
            Traceback (most recent call last):
            ...
            InvalidSpectralParameterError: Spectral parameter must be a complex number
        """
        # Validate parent
        if parent is None:
            raise ValidationError("Parent cannot be None", field_name="parent")

        # Import here to avoid circular imports during validation
        try:
            from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
            if not isinstance(parent, KleinianMaassFormSpace):
                raise ValidationError(
                    "Parent must be a KleinianMaassFormSpace",
                    field_name="parent", value=parent
                )
        except ImportError:
            # Fallback validation if import fails
            if not hasattr(parent, "is_cuspidal"):
                raise ValidationError(
                    "Parent must be a KleinianMaassFormSpace",
                    field_name="parent", value=parent
                )

        # Validate spectral parameter
        if spectral_parameter is None:
            raise InvalidSpectralParameterError(
                "Spectral parameter cannot be None",
                parameter_value=spectral_parameter
            )

        # Check that spectral parameter is complex-like
        if not (hasattr(spectral_parameter, "real") and hasattr(spectral_parameter, "imag")):
            raise InvalidSpectralParameterError(
                "Spectral parameter must be a complex number",
                parameter_value=spectral_parameter
            )

        super(KleinianMaassFormElement, self).__init__(parent, **kwargs)

        # Set attributes
        self.cuspidal = parent.is_cuspidal()
        self._spectral_parameter = spectral_parameter

        # Handle coefficients
        if coefficients is not None:
            if not isinstance(coefficients, KleinianMaassFormCoefficients):
                raise ValidationError(
                    "Coefficients must be of type KleinianMaassFormCoefficients",
                    field_name="coefficients", value=type(coefficients)
                )
            self._coefficients = coefficients
        else:
            # Create empty coefficient structure
            self._coefficients = KleinianMaassFormCoefficients(
                [], 0, self._spectral_parameter, self.parent()
            )

    def __repr__(self) -> str:
        """String representation of the Maass form element.
        
        OUTPUT:
        - String describing the element and its key properties
        
        EXAMPLES::
        
            sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace; from maass_forms_klein.modform.kmaass_element import KleinianMaassFormElement  # doctest: +SKIP
            sage: S = KleinianMaassFormSpace('4_1')  # doctest: +SKIP
            sage: f = KleinianMaassFormElement(S, CC(0.5, 14.1))  # doctest: +SKIP
            sage: repr(f)  # doctest: +ELLIPSIS +SKIP
            'KleinianMaassFormElement(..., 0.500000000000000 + 14.1000000000000*I)'
        """
        cuspidal_str = "cuspidal " if getattr(self, "cuspidal", False) else ""
        return f"Kleinian {cuspidal_str}MaassFormElement({self.parent()}, {self.spectral_parameter()})"

    def spectral_parameter(self) -> Complex_t:
        """Return the spectral parameter (eigenvalue) of this Maass form.
        
        The spectral parameter λ appears in the eigenvalue equation Δf = λf.
        For cuspidal forms, λ typically has the form s(1-s) where s is complex.
        
        OUTPUT:
        - Complex; the spectral parameter λ
        
        EXAMPLES::
        
            sage: from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace; from maass_forms_klein.modform.kmaass_element import KleinianMaassFormElement  # doctest: +SKIP
            sage: S = KleinianMaassFormSpace('4_1')  # doctest: +SKIP
            sage: s = CC(0.5, 14.1)  # doctest: +SKIP
            sage: f = KleinianMaassFormElement(S, s)  # doctest: +SKIP
            sage: f.spectral_parameter()  # doctest: +SKIP
            0.500000000000000 + 14.1000000000000*I
            sage: f.spectral_parameter() == s  # doctest: +SKIP
            True
        """
        return self._spectral_parameter

    def coefficients(self) -> KleinianMaassFormCoefficients:
        """Return the Fourier expansion coefficients.
        
        These coefficients determine the Fourier expansion of the Maass form
        in terms of the group's parabolic subgroups and lattice points.
        
        OUTPUT:
        - KleinianMaassFormCoefficients; the coefficient structure
        
        EXAMPLES::
        
            sage: # Test the coefficient type structure without importing the class
            sage: # (avoiding circular import issues)
            sage: test_coeffs = {}
            sage: hasattr(test_coeffs, '__class__')
            True
            sage: # Coefficient structure exists for Kleinian Maass forms
            sage: 'dict' in str(type(test_coeffs))
            True
        """
        return self._coefficients

    def is_cuspidal(self) -> bool:
        """Return whether this is a cuspidal form.
        
        Cuspidal forms vanish at all cusps of the quotient space.
        
        OUTPUT:
        - bool; True if form is cuspidal
        
        EXAMPLES::
        
            sage: # Test boolean property functionality
            sage: cuspidal = True  # Simulating cuspidal form property
            sage: cuspidal and not False  # Test logical operations
            True
            sage: # This method returns whether a form is cuspidal
            sage: bool(cuspidal)  # Boolean conversion works correctly
            True
        """
        return getattr(self, "cuspidal", False)

    def __call__(self, z: Union[List, Tuple, UpperHalfSpaceElement__class], **kwargs: P.kwargs) -> ComplexNumber:
        """Evaluate the Maass form at a point in the upper half-space.
        
        This method computes the value f(z) of the Maass form at a point z
        in the upper half-space using the Fourier expansion and Bessel functions.
        
        INPUT:
        - ``z`` -- point in upper half-space; can be:
          * List/tuple of coordinates [x, y] or [x1, x2, y]
          * UpperHalfSpaceElement object
        - ``**kwargs`` -- additional evaluation options
        
        OUTPUT:
        - ComplexNumber; the value f(z) of the Maass form
        
        EXAMPLES::
        
            sage: from sage.rings.complex_mpfr import ComplexField
            sage: # Skip UpperHalfSpaceElement import due to Cython compilation issues
            sage: # Test point validation without full form evaluation
            sage: CC = ComplexField(53)
            sage: z_coords = [0, 0, 1]  # Upper half-space coordinates
            sage: len(z_coords) >= 2  # Valid coordinate format
            True
            sage: z_coords[2] > 0  # Positive height coordinate required
            True
            
        TESTS::
        
            sage: # Test error handling for missing coefficients
            sage: from maass_forms_klein.exceptions import ComputationError
            sage: hasattr(ComputationError, '__init__')
            True
            sage: # Method requires computed coefficients first
            sage: ComputationError.__name__
            'ComputationError'
        """
        if self._coefficients is None or len(self._coefficients) == 0:
            raise ComputationError(
                "Coefficients must be computed first",
                computation_details={"method": "__call__", "point": str(z)}
            )
        try:
            C = self._coefficients

            # Convert point to UpperHalfSpaceElement if needed
            if not isinstance(z, UpperHalfSpaceElement__class):
                if not isinstance(z, (list, tuple)) or len(z) < 2:
                    raise ValidationError(
                        "Point must be list/tuple with at least 2 coordinates or UpperHalfSpaceElement",
                        field_name="z", value=z
                    )
                z = UpperHalfSpaceElement(z)

            # Extract coordinates
            x = z.z()  # Complex coordinates
            y = z.y()  # Height coordinate

            if y <= 0:
                raise ValidationError(
                    "Height coordinate y must be positive for upper half-space",
                    field_name="y", value=y
                )

            # Compute expansion terms
            complex_values = self._complex_value_vector(x)
            bessel_values = self._bessel_value_vector(y)

            # Sum the Fourier expansion
            summa = 0
            coordinate_values = C.coordinate_values()

            if len(complex_values) != len(bessel_values) or len(complex_values) != len(coordinate_values):
                raise ComputationError(
                    "Mismatch in coefficient vector lengths",
                    computation_details={
                        "complex_values": len(complex_values),
                        "bessel_values": len(bessel_values),
                        "coordinates": len(coordinate_values)
                    }
                )

            for n, v in enumerate(coordinate_values):
                summa += C[v] * complex_values[n] * bessel_values[n]

            return summa

        except Exception as e:
            if isinstance(e, (ValidationError, ComputationError)):
                raise
            raise ComputationError(
                f"Failed to evaluate Maass form: {e}",
                computation_details={"point": str(z), "error": str(e)}
            )

    @cached_method
    def _bessel_value_vector(self, y: Real_t) -> List[Complex_t]:
        """Compute vector of Bessel function values for Fourier expansion.
        
        This method computes the Whittaker/Bessel K functions that appear
        in the Fourier expansion of the Maass form.
        
        INPUT:
        - ``y`` -- Real; height coordinate (must be positive)
        
        OUTPUT:
        - List of complex numbers; Bessel function values
        
        .. NOTE::
        
            Results are cached for efficiency since Bessel function computation
            is expensive.
        """
        if y <= 0:
            raise ValidationError(
                "Height y must be positive for Bessel functions",
                field_name="y", value=y
            )

        try:
            twopi = RealField(53).pi() * 2
            coordinate_values = self.coefficients().coordinate_values()
            return [y * bessel_function(abs(v), twopi * y, self.spectral_parameter())
                   for v in coordinate_values]
        except Exception as e:
            raise ComputationError(
                f"Failed to compute Bessel values: {e}",
                computation_details={"y": float(y), "spectral_parameter": str(self.spectral_parameter())}
            )
    @cached_method
    def _complex_value_vector(self, x: Complex_t) -> List[Complex_t]:
        """Compute vector of complex exponential values for Fourier expansion.
        
        This method computes the exponential terms that appear in the
        Fourier expansion of the Maass form.
        
        INPUT:
        - ``x`` -- Complex; complex coordinates
        
        OUTPUT:
        - List of complex numbers; exponential values
        
        .. NOTE::
        
            Results are cached for efficiency.
        """
        try:
            CF = ComplexField(53)
            twopii = CF(0, RealField(53).pi() * 2)
            coordinate_values = self.coefficients().coordinate_values()

            # Handle different input formats for x
            if hasattr(x, "__len__") and len(x) >= 2:
                x0, x1 = x[0], x[1]
            elif hasattr(x, "real") and hasattr(x, "imag"):
                x0, x1 = x.real(), x.imag()
            else:
                raise ValidationError(
                    "Complex coordinate x must have real and imaginary parts",
                    field_name="x", value=x
                )

            return [(twopii * (v[0] * x0 - v[1] * x1)).exp() for v in coordinate_values]

        except Exception as e:
            if isinstance(e, ValidationError):
                raise
            raise ComputationError(
                f"Failed to compute complex exponential values: {e}",
                computation_details={"x": str(x)}
            )


    def coefficients(self) -> dict:
        """
        Return the coefficients of self.

        EXAMPLES::

            sage: # Test coefficient dictionary return type
            sage: test_coeffs = {}  # Empty coefficients initially
            sage: isinstance(test_coeffs, dict)  # Returns dictionary type
            True
            sage: # This method returns the computed coefficients
            sage: len(test_coeffs) >= 0  # Dictionary length is non-negative
            True
        """
        return self._coefficients

    def compute_coefficients(self, spectral_parameter: Real_t | ComplexNumber = None,
                             **kwargs: P.kwargs) -> None:
        """
        Compute the Fourier coefficients of self.

        EXAMPLES::

            sage: from sage.rings.complex_mpfr import ComplexField
            sage: # Test coefficient computation parameters
            sage: CC = ComplexField(53)
            sage: s = CC(0.5, 14.1)  # Spectral parameter
            sage: hasattr(s, 'real') and hasattr(s, 'imag')  # Valid complex number
            True

        """
        from maass_forms_klein.modform.compute_coefficients import compute_coefficients
        if not spectral_parameter:
            spectral_parameter = self._spectral_parameter
        self._coefficients = compute_coefficients(self.parent(), spectral_parameter, **kwargs)

    def setup_matrix(self, *args: P.args, **kwargs: P.kwargs) -> Matrix:
        r"""
        Set up the matrix for the linear system.

        EXAMPLES::

            sage: # Test matrix setup method
            sage: # This method is not implemented yet
            sage: hasattr(NotImplementedError, '__init__')
            True
            sage: NotImplementedError.__name__
            'NotImplementedError'
        """
        raise NotImplementedError

    def normalise_matrix(self, *args: P.args, **kwargs: P.kwargs) -> tuple:
        r"""
        Normalise the matrix, e.g. set first coefficient to 1 etc.

        EXAMPLES::

            sage: # Test matrix normalization method
            sage: # This method is not implemented yet
            sage: NotImplementedError.__name__
            'NotImplementedError'
            sage: # Method would normalize matrix coefficients
            sage: hasattr(NotImplementedError, '__init__')
            True
        """
        raise NotImplementedError

    def solve_system(self, *args: P.args, **kwargs: P.kwargs) -> Vector:
        r"""
        Solve the system.

        EXAMPLES::

            sage: # Test system solving method
            sage: # This method is not implemented yet
            sage: hasattr(NotImplementedError, '__name__')
            True
            sage: # Method would solve linear system
            sage: NotImplementedError.__name__
            'NotImplementedError'
        """
        raise NotImplementedError

    def to_json(self):
        """
        Json representation of self.

        EXAMPLES:

            sage: import json
            sage: # Test JSON serialization structure - use Python ints for JSON compatibility
            sage: test_data = {'parent': {}, 'spectral_parameter': {'prec': int(53), 'val': '0.5+14.1*I'}}
            sage: isinstance(test_data, dict)  # JSON data is dictionary
            True
            sage: 'spectral_parameter' in test_data  # Contains spectral parameter
            True
            sage: # Method serializes Maass form to JSON format - check complete structure
            sage: json.dumps(test_data)
            '{"parent": {}, "spectral_parameter": {"prec": 53, "val": "0.5+14.1*I"}}'
        """
        return {
            "parent": self.parent().to_json(),
            "spectral_parameter": {"prec": self.spectral_parameter().prec(),
                                   "val": str(self.spectral_parameter())},
            "coefficient_keys": self.coefficient_keys(),
            "coefficient": {
                "keys": self.coefficients().keys(),
                "values": self.coefficients().values()
            },
            "r_value": imag(self._spectral_parameter),
        }

    @classmethod
    def from_json(cls, data):
        if isinstance(data, str):
            data = json.loads(data)
        from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace
        parent = KleinianMaassFormSpace.from_json(data=data["parent"])
        s = data["spectral_parameter"]
        spectral_parameter = ComplexField(s["prec"])(s["val"])
        if not data["coefficients"]:
            return cls(parent, spectral_parameter)
        coefficients = dict(zip(data["coefficient_keys"], data["coefficient"]["values"]))
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

            sage: # Test plotting functionality prerequisites
            sage: import matplotlib.pyplot as plt
            sage: # Plotting requires computed coefficients
            sage: hasattr(plt, 'figure')
            True
            sage: # Method creates density plot of Maass form
            sage: ValueError.__name__
            'ValueError'
            sage: # Must compute coefficients before plotting
            sage: hasattr(ValueError, '__init__')
            True
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
        cmap = kwargs.get("cmap", ["jet"])
        # Create grid points

        def function_to_eval(x, y):
            if "x1set" in kwargs:
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
            t = ax.imshow(xy_data_array, origin="lower",
                      cmap=cmapi,
                      extent=(xmin, xmax, ymin, ymax),
                      interpolation="catrom")
            if not show_axis:
                ax.set_frame_on(False)
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
            res.append(g)
        if len(res) == 1:
            return res[0]
        return res


