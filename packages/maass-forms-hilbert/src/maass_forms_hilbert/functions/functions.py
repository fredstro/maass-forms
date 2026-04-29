from sage.misc.cachefunc import cached_function
from sage.misc.misc_c import prod

from maass_forms_hilbert.functions.functions_cy import bessel_prod_dp2
from maass_forms_hilbert.modform.utils import Complex_t, Integer_t, Real_t

try:
    from maass_forms_hilbert.functions.bessel.besselk_dp import besselk_dp
except ImportError:
    besselk_dp = None
from sage.functions.bessel import bessel_K
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.real_mpfr import RealNumber as RealNumber_class


@cached_function
def bessel_prod(
    v: tuple[Real_t, ...],
    y: tuple[Real_t, ...],
    s: tuple[Complex_t, ...],
    sgn: str = "+",
    use_iR: bool = False,
    prec: Integer_t = 53,
) -> RealNumber_class:
    r"""
    A product of scaled K-Bessel functions: sqrt(y_i) e^{pi R_i/2}K_{iR_i}(2pi |v_i|y_i)
    where s_i = 1/2 + Ri

    INPUT:

    - ``v`` -- tuple of real numbers
    - ``y`` -- tuple of real numbers
    - ``s`` -- tuple of complex numbers
    - ``sgn`` -- '+' or '-'
    - ``use_iR`` -- boolean, if True compute K_{iR}(x) where s_i = 1/2 + Ri
    - ``prec`` -- precision

    OUTPUT:

    - The product of K-Bessel functions as a real number

    EXAMPLES::

        sage: from maass_forms_hilbert.functions.functions import bessel_prod
        sage: bessel_prod((1.0,), (1.0,), (1.0,))
        0.000933721365853991
    """
    n = len(v)
    if len(y) != n or len(s) != n:
        raise ValueError("Need tuples of same length")
    CF = ComplexField(prec)
    if n == 2 and use_iR and prec == 53:
        return bessel_prod_dp2(v[0], v[1], y[0], y[1], s[0], s[1], sgn == "+")
    if all(vi == 0 for vi in v):
        if sgn == "+":
            return prod(CF(y[i] ** s[i]) for i in range(n))
        else:
            return prod(CF(y[i] ** (CF(1) - s[i])) for i in range(n))
    s_minus_half = [si - CF(1) / CF(2) for si in s]
    twopi = CF(2) * CF.pi()
    if not all(si.real() == 0 for si in s_minus_half) or not besselk_dp or CF.prec() > 53:
        bessels = [
            CF(
                y[i].sqrt()
                * (CF.pi() / 2 * s_minus_half[i].imag()).exp()
                * bessel_K(s_minus_half[i], twopi * abs(v[i]) * y[i])
            )
            for i in range(n)
        ]
    else:
        R = [si.imag() for si in s_minus_half]
        bessels = [
            CF(y[i]).sqrt() * besselk_dp(R[i], twopi * abs(v[i]) * y[i], pref=1) for i in range(n)
        ]
    return prod(bessels)


@cached_function
def exp_trace_prod(x: tuple, prec: int = 53) -> ComplexNumber:
    """
    Return e( trace(x) )

    EXAMPLES::

        sage: from maass_forms_hilbert.functions.functions import exp_trace_prod
        sage: exp_trace_prod((1.0,1.5)) # abs tol 1e-10
        -1.00000000000000 - 3.31384565490311e-14*I
    """
    summa = 0
    CF = ComplexField(prec)
    twopi = 2 * CF.pi()
    for xi in x:
        summa += xi
    return CF(0, twopi * summa).exp()
