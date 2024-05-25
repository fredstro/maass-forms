r"""
Functions.
"""
from .bessel.besselk_dp cimport besselk_dp_c
from .bessel.besselk_dp import besselk_dp
from cysignals.memory cimport sig_free, check_malloc
cdef extern from "complex.h":
    cdef double complex _Complex_I
    cdef double complex cexp(double complex)
    cdef double creal(double complex)
    cdef double cimag(double complex)
    cdef double complex cpow(double complex, double complex)

cdef extern from "math.h":
    cdef double sqrt(double)
    cdef double cos(double)
    cdef double sin(double)

cdef double twopi = 6.28318530717958647692528676656
cdef double complex twopii = _Complex_I * 6.28318530717958647692528676656
from sage.functions.bessel import bessel_K

cpdef exp_trace_prod_dp(x, symmetry=0):
    """
    Return e( trace(x) )

    EXAMPLES::

        sage: from hilbert_maass.functions.functions_cy import exp_trace_prod_dp
        sage: exp_trace_prod_dp((1.0,1.5)) # abs tol 1e-15
        (-1+6.123233995736766e-16j)
    """
    cdef double summa = 0.0
    cdef double xi
    for xi in x:
        summa += <double>xi
    if symmetry == 0:
        return cexp(twopii*summa)
    if symmetry == 1:
        return cos(twopi*summa)
    if symmetry == -1:
        return sin(twopi*summa)

cpdef bessel_prod_dp2(double v0, double v1, double y0, double y1,
                      double complex s0, double complex s1, int sgn):
    """
    Special case of degree 2. A product of scaled K-Bessel functions: sqrt(y_i) e^{pi R_i/2}K_{iR_i}(2pi |v_i|y_i)
    where s_i = 1/2 + Ri
    
    EXAMPLES::

        sage: from hilbert_maass.functions.functions_cy import bessel_prod_dp2
        sage: bessel_prod_dp2(1,1, 1,1, CC(0.5,1.0), CC(0.5,1.0), 1) # abs tol 1e-14
        1.6758678170356153e-05
        sage: bessel_prod_dp2(1,1, 1,1, CC(0.6,1.0), CC(0.6,1.0), 0)
        (7.249717552777861e-07+2.1566552189937963e-08j)

    """
    cdef int i = 0
    cdef int all_v_zero = 1
    cdef double complex product = 1
    if v0 == 0 and v1 == 0:
        if sgn:
            return cpow(y0, s0) * cpow(y1, s1)
        return  cpow(y0, 1.0 - s0) * cpow(y1, 1.0 - s1)
    #s_minus_half = [si - CF(1) / CF(2) for si in s]
    cdef int all_si_real = 1
    cdef double R0, R1
    cdef double bessel_prod = 1
    cdef double complex bessel_prodc = 1
    if creal(s0 - 0.5) != 0 or creal(s1 - 0.5) != 0:
        bessel_prodc = sqrt(y0) * bessel_K(s0 - 0.5, twopi * abs(v0) * y0)
        bessel_prodc = bessel_prodc * sqrt(y1) * bessel_K(s1 - 0.5, twopi * abs(v1) * y1)
        return bessel_prodc
    R0 = cimag(s0 - 0.5)
    R1 = cimag(s1 - 0.5)
    cdef double prec = 1e-14
    cdef double * bes = NULL
    cdef int res
    bes = <double *> check_malloc(sizeof(double) * 1)
    res = besselk_dp_c(bes, R0, twopi * abs(v0) * y0, prec, 1)
    bessel_prod = sqrt(y0) * bes[0]
    res = besselk_dp_c(bes, R1, twopi * abs(v1) * y1, prec, 1)
    bessel_prod = bessel_prod * sqrt(y1) * bes[0]
    sig_free(bes)
    return bessel_prod

cpdef bessel_prod_dp_gen(tuple v, tuple y, tuple s, str sgn='+'):
    """
    A product of scaled K-Bessel functions: sqrt(y_i) e^{pi R_i/2}K_{iR_i}(2pi |v_i|y_i)
    where s_i = 1/2 + Ri
    
        EXAMPLES::

        sage: from hilbert_maass.functions.functions_cy import bessel_prod_dp_gen
        sage: bessel_prod_dp_gen((1,1), (1,1), (CC(0.5,1.0), CC(0.5,1.0)), '+') # abs tol 1e-14
        1.6758678170356153e-05

    """
    cdef int n = len(v)
    if len(y) != n:
        raise ValueError("Need vectors of same length")
    cdef int i = 0
    cdef int all_v_zero = 1
    cdef double complex product = 1
    for i in range(n):
        if v[i]:
            all_v_zero = 0
            break
    if all_v_zero == 1:
        if sgn == '+':
            for i in range(n):
              product = product * cpow(y[i], s[i])
        else:
            for i in range(n):
                product = product * cpow(y[i], 1.0 - s[i])
        return product
    #s_minus_half = [si - CF(1) / CF(2) for si in s]
    cdef int all_si_real = 1
    for i in range(n):
        if creal(s[i]-0.5) != 0:
            all_si_real = 0
            break
    cdef double R
    cdef double bessel_prod = 1
    cdef double complex bessel_prodc = 1
    if not all_si_real or not besselk_dp:
        for i in range(n):
            bessel_prodc = bessel_prodc * sqrt(y[i]) * bessel_K(s[i] - 0.5, twopi * abs(v[i]) * y[i])
        return bessel_prodc
    else:
        for i in range(n):
            R = (s[i] - 0.5).imag()
            bessel_prod = bessel_prod * sqrt(y[i])*besselk_dp(R, twopi * abs(v[i]) * y[i], pref=1)
        return bessel_prod
