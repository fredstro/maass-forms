from sage.misc.cachefunc import cached_function
from sage.misc.misc_c import prod
try:
    from hilbert_maass.functions.bessel.besselk_dp import besselk_dp
except ImportError:
    besselk_dp = None
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.real_mpfr import RealNumber as RealNumber_class
from sage.functions.bessel import bessel_K


@cached_function
def bessel_prod(v: tuple, y: tuple, s: tuple, sgn: str = '+') -> RealNumber_class:
    """
    A product of scaled K-Bessel functions: sqrt(y_i) e^{pi R_i/2}K_{iR_i}(2pi |v_i|y_i)
    where s_i = 1/2 + Ri
    """
    n = len(v)
    if len(y) != n:
        raise ValueError("Need vectors of same length")
    i = 0
    CF = s[0].parent()
    if all(vi == 0 for vi in v):
        if sgn == '+':
            return prod(CF(y[i] ** s[i]) for i in range(n))
        else:
            return prod(CF(y[i] ** (CF(1) - s[i])) for i in range(n))
    s_minus_half = [si - CF(1) / CF(2) for si in s]
    twopi = CF(2) * CF.pi()
    if not all(si.real() == 0 for si in s_minus_half) or not besselk_dp:
        bessels = [
            CF(y[i].sqrt() * bessel_K(s_minus_half[i], twopi * abs(v[i]) * y[i])) for i in range(n)
        ]
    else:
        R = [si.imag() for si in s_minus_half]
        bessels = [
            CF(y[i]).sqrt()*besselk_dp(R[i], twopi * abs(v[i]) * y[i], pref=1) for i in range(n)
        ]
    return prod(bessels)

@cached_function
def exp_trace_prod(x: tuple, prec: int = 53) -> ComplexNumber:
    """
    Return e( trace(x) )
    """
    summa = 0
    CF = ComplexField(prec)
    twopi = 2*CF.pi()
    for xi in x:
        summa += xi
    return CF(0, twopi*summa).exp()
