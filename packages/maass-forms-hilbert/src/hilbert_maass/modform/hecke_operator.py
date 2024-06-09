"""
Hecke operators acting on Hilbert Maass forms.
"""
from copy import copy

from hilbert_maass.modform.coefficients import HilbertMaassCoefficients
from hilbert_maass.modform.hilbert_maass_element import HilbertMaassForm, HilbertMaassForm_Element
from hilbert_maass.modform.hilbert_maass_space import HilbertMaassFormSpace
from hilbert_maass.modform.utils import ideal_coordinates, ideal_factors, Integer_t
from sage.matrix.constructor import matrix
from sage.rings.number_field.number_field_ideal import NumberFieldIdeal
from sage.structure.element import Element


class HeckeOperator(Element):
    """
    Hecke operators acting on Hilbert Maass forms.
    """
    def __init__(self, space: HilbertMaassFormSpace, n: NumberFieldIdeal) -> None:
        """
        Initialize ``self``.


        EXAMPLES::

            sage: from hilbert_maass.modform.hilbert_maass_element import HilbertMaassForm
            sage: from hilbert_maass.modform.hecke_operator import HeckeOperator
            sage: form = HilbertMaassForm(QuadraticField(2), (1, 1))
            sage: HeckeOperator(form.parent(), form.parent().number_field().ideal(2))
            Hecke operator of index Fractional ideal (2) acting on HilbertMaassFormSpace(Hilbert
            Modular Group PSL(2) over Maximal Order in Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?)
            sage: HeckeOperator(form.parent(), 2)
            Hecke operator of index Fractional ideal (2) acting on HilbertMaassFormSpace(Hilbert
            Modular Group PSL(2) over Maximal Order in Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?)
            sage: HeckeOperator(form.parent(), 3)
            Hecke operator of index Fractional ideal (3) acting on HilbertMaassFormSpace(Hilbert
            Modular Group PSL(2) over Maximal Order in Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?)
        """
        if space.number_field().narrow_class_group().order() > 1:
            raise NotImplementedError("Only narrow class number 1 supported")
        if not isinstance(space, HilbertMaassFormSpace):
            raise TypeError("Not a Hilbert Maass form space")
        if n in space.number_field() or isinstance(n, Integer_t):
            n = space.number_field().ideal(n)
        if not isinstance(n, NumberFieldIdeal) or not n.is_integral():
            raise TypeError("Not an integral ideal")
        if n.number_field() != space.number_field():
            raise ValueError("Not in same number field")
        self.space = space
        self.ideal = n
        self.ideal_gen = n.gens_reduced()[0]
        if not self.ideal_gen.is_totally_positive():
            raise ArithmeticError("Generator is not totally positive")

    def __repr__(self) -> str:
        """
        Return string representation of ``self``.

        """
        return f"Hecke operator of index {self.ideal} acting on {self.space}"

    def __call__(self, x: HilbertMaassForm_Element) -> HilbertMaassForm_Element:
        """
        Apply Hecke operator to ``x``.

        INPUT:

        - ``x`` -- vector of integers

        """
        if not isinstance(x, HilbertMaassForm_Element):
            raise TypeError("Not a Hilbert Maass form")
        if self.space != x.parent():
            raise ValueError("Not in same space")
        if not x.coefficients():
            raise ValueError("Need to compute coefficients first")
        # Transform coefficients...
        new_coeffs = {}
        divisors = []
        for idd in ideal_factors(self.ideal):
            d = idd.gens_reduced()[0]
            if not d.is_totally_positive():
                raise ArithmeticError("Generator is not totally positive")
            divisors.append(d)
        # print("divisors=", divisors)
        list_of_coordinates = []
        new_coeffs = []
        for v in x.coefficients().keys(as_elements=True):
            vcoord = ideal_coordinates(x.coefficients()._coordinate_ideals[0], v)
            try:
                new_coeff = 0
                w_used = []
                for d in divisors:
                    if v / d not in x.coefficients()._coordinate_ideals[0]:
                        continue
                    w = v * self.ideal_gen / d ** 2
                    new_coeff += x.coefficients()[w]
                    w_used.append(w)
                new_coeffs.append(new_coeff)
                list_of_coordinates.append(vcoord)
                # print(vcoord, w_used)
            except IndexError:
                # If the coordinates are not in the original coefficients we omit all
                continue
        xcoeffs = x.coefficients()
        new_coeffs = matrix([[new_coeff] for new_coeff in new_coeffs])
        coeffs = HilbertMaassCoefficients(new_coeffs,
                                          xcoeffs.M,
                                          spectral_parameter=x.spectral_parameter(),
                                          space=x.parent(),
                                          coordinate_ideals=xcoeffs._coordinate_ideals,
                                          set_coefficients=xcoeffs._set_coefficients,
                                          index_tuples=[list_of_coordinates],
                                          check=False,
                                          Y=xcoeffs.Y())
        result = copy(x)
        result._coefficients = coeffs
        return result
        # return HilbertMaassForm_Element(self.space, self.ideal, new_coeffs)


