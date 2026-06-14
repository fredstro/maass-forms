# -*- coding=utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010  Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from libc.stdio cimport *
import cython

r"""
Low-precision (fast) algorithms for the K-Bessel function.
Also algorithms for incomplete gamma function.

    AUTHOR:

    - Fredrik Strömberg (April 2010)


    EXAMPLES::

        sage: from maass_form_core.functions.bessel.besselk_dp import besselk_dp, besselk_dp_pow
        sage: from maass_form_core.functions.bessel.besselk_dp import besselk_dp_rec
        sage: from sage.all import RealField, ComplexField
        sage: besselk_dp(10.0, 5.0) # tol 2e-14
        -0.7183327166568183
        sage: besselk_dp_rec(10.0, 3.0) # tol 2e-14
        -6.3759939798738967e-08
        sage: besselk_dp_rec(10.0, 3.0, pref=1) # tol 2e-14
        -0.42308698672505796
        sage: besselk_dp_pow(10.0, 3.0) # tol 2e-14
        -6.3759939798739404e-08
        sage: besselk_dp_pow(10.0, 3.0, pref=1) # tol 2e-14
        -0.42308698672506084
        sage: a = besselk_dp_pow(10.0, 3.0, pref=1); a  # tol 2e-14
        -0.42308698672506084
        sage: b = besselk_dp_rec(10.0, 3.0, pref=1); b  # tol 2e-14
        -0.42308698672505796
        sage: abs(a - b) # abs tol 2e-14
        0
        sage: besselk_dp_rec(100.0, 3.0, pref=1)  # tol 3e-13
        0.082457701468151401
        sage: RF = RealField(150)
        sage: CF = ComplexField(150)
        sage: bessel_K(CF(100 * I), RF(3)) * exp(RF.pi() * RF(50))  # tol 1e-14
        0.0824577014681303
        sage: _ - besselk_dp_rec(100.0, 3.0, pref=1)  # tol 5e-14
        0


"""

cdef extern from "complex.h":
    cdef double complex _Complex_I

cdef complex CMPLXF(float x, float y):
    return x + _Complex_I * y


cdef double complex CMPLX(double x, double y):
    return x + _Complex_I * y
cdef extern from "math.h" nogil:
    double log(double)
    double exp(double)
    double cos(double)
    double sin(double)
    double atan(double)
    double sinh(double)
    double fabs(double)
    double sqrt(double)
    int ceil(double)
    double cabs(double complex)
    double pow(double, double)
    double complex clog(double complex)
    double cimag(double complex)
    double creal(double complex)
    double carg(double complex)
    double complex cexp(double complex)
    

cdef double cppi = <double>3.14159265358979323846264338327950288419716939  # =pi
cdef double pihalf = <double>1.5707963267948966192313216916397514420985847  # =pi/2
cdef double ln2PI2 = <double>0.91893853320467274178032973640561763986139747  # =ln(2pi)/2
cdef double d_one = <double>1.0
cdef double d_half = <double>0.5
cdef double d_two = <double>2.0
cdef double d_ten = <double>10.0
cdef double d_20 = <double>20.0
cdef double d_100 = <double>100.0
cdef double complex c_one = CMPLX(1.0, 0.0)
cdef double complex c_zero = CMPLX(0.0, 0.0)


@cython.cdivision(True)
cpdef besselk_dp(double R, double x, double prec=1e-14, int pref=1, algorithm="default"):
    r"""
    Modified K-Bessel function in double precision. Chooses the most appropriate algorithm.

    INPUT:

        - `R` -- parameter (double)
        - `x` -- argument (double)
        - `prec` -- precision (default 1E-14, double)
        - `pref` -- use prefactor (integer, default 1)
               = 1 => computes K_iR(x)*exp(pi*R/2)
               = 0 => computes K_iR(x)
               = -1 => computes K_R(x)
        - `algorithm` -- use specific algorithm instead of "best"
               = 'default' => use automatic choice
               = 'pow' => use power series
               = 'rec' => use recursion
    OUTPUT:
    
     - exp(Pi*R/2)*K_{i*R}(x)  -- double

    EXAMPLES::

        sage: from maass_form_core.functions.bessel.besselk_dp import besselk_dp
        sage: besselk_dp(0.0, 1.0) # tol 2e-14
        0.42102443824070833
        sage: besselk_dp(10.0, 5.0) # tol 2e-14
        -0.7183327166568183
        sage: besselk_dp(10.0, 3.0, pref=0) # tol 2e-14
        -6.3759939798738967e-08
        sage: besselk_dp(10.0, 3.0, pref=1) # tol 2e-14
        -0.42308698672505796
        

    """
    RR = fabs(R)
    if x <= 0:
        raise ValueError(f" Need x>0! Got x={x}")
    if pref == -1:
        return besselk_real_dp(R, x, prec)
    cdef double xc, xcral, S
    xc = sqrt((x + R) * (x - R))
    S = R / x
    xcral = xc + R * 2.0 * atan(S / (1.0 + (xc / x)))
    cdef double kbes
    if R * pihalf - xcral < -125.0:
        return 0.0
    cdef int res = -1
    if algorithm != "default":
        if algorithm == "pow":
            res = besselk_dp_pow_c(RR, x, &kbes, prec, pref)
        elif algorithm == "rec":
            res = besselk_dp_rec_c(RR, x, &kbes, prec, pref)
        if res == 1:
            raise ValueError(
                f"The K-bessel routine failed (k large) for x,R={x},{RR}, value={kbes}"
            )
        elif res == 2:
            raise ValueError(
                f"The K-bessel routine failed (too many iterations) for x,R={x},{RR}, value={kbes}"
            )
        return kbes
    if x < R * 0.7:
        res = besselk_dp_pow_c(RR, x, &kbes, prec, pref)
        if res != 0:
            res = besselk_dp_rec_c(RR, x, &kbes, prec, pref)
    else:
        res = besselk_dp_rec_c(RR, x, &kbes, prec, pref)
    if res == 1:
        raise ValueError(f"The K-bessel routine failed (k large) for x,R={x},{RR}, value={kbes}")
    elif res == 2:
        raise ValueError(
            f"The K-bessel routine failed (too many iterations) for x,R={x},{RR}, value={kbes}"
        )
    elif res == 0:
        return kbes
    else:
        raise ValueError(f"The K-bessel routine failed (unknown error) for x,R={x},{RR}, "
                         f"value={kbes}")

cdef int besselk_dp_c(
    double *kbes, double R, double x, double prec, int pref, int verbose=0
) nogil:  # double prec=1e-14,int pref=0):
    r"""
    Modified K-Bessel function in double precision. Chooses the most appropriate algorithm.

    INPUT:

       - `R` -- parameter (double)
       - `x` -- argument (double)
       - `prec` -- precision (default 1E-14, double)
       - `pref` -- use prefactor (integer, default 0)
               = 1 => computes K_iR(x)*exp(pi*R/2)
               = 0 => computes K_iR(x)
       - ``verbose`` -- int (default: 0) set to 1 to print debug output.

    OUTPUT:
    
     - exp(Pi*R/2)*K_{i*R}(x)  -- double

    """
    RR = fabs(R)
    if x <= 0:
        if verbose:
            printf("Need x>0! Got x=%f", x)
        return -1
    cdef int res = 0
    if x > 0.5 * cppi * RR - 2.0 * log(prec):
        kbes[0] = 0.0
        return 0
    if x < R * 0.7:
        res = besselk_dp_pow_c(RR, x, kbes, prec, pref)
        if res != 0:
            res = besselk_dp_rec_c(RR, x, kbes, prec, pref)
    else:
        res = besselk_dp_rec_c(RR, x, kbes, prec, pref)

    if res == 1:
        printf(
            "The K-bessel routine failed (k large) for x,R=%g,%g, value=%g \n",
            x,
            RR,
            kbes[0],
        )
    elif res == 2:
        printf(
            "the K-bessel routine failed (too many iterations) for x,R=%g,%g, value=%g",
            x,
            RR,
            kbes[0],
        )
    return res

@cython.cdivision(True)
cpdef double besselk_dp_rec(
    double R, double x, double prec=1e-14, int pref=0
):  # |r| <<1000, x >> 0 !
    r"""
    Modified K-Bessel function in double precision using the backwards Miller-recursion algorithm. 

    INPUT:
     - ''R'' -- double
     - ''x'' -- double
     - ''prec'' -- (default 1E-12) double
     - ''pref'' -- (default 1) int
                =1 => computes K_iR(x)*exp(pi*R/2)
                =0 => computes K_iR(x)

    OUTPUT:
    
     - exp(Pi*R/2)*K_{i*R}(x)  -- double

    EXAMPLES::

        sage: from maass_form_core.functions.bessel.besselk_dp import besselk_dp_rec
        sage: besselk_dp_rec(10.0, 5.0, prec=1e-15, pref=1) # tol 2e-13
        -0.7183327166568183
        sage: besselk_dp_rec(10.0, 3.0, prec=1e-15) # tol 2e-14
        -6.375993979873876e-08
        sage: besselk_dp_rec(10.0, 3.0, prec=1e-15, pref=1) # tol 2e-14
        -0.42308698672505657
        sage: besselk_dp_rec(100.0, 3.0, prec=1e-15, pref=1) # tol 4e-13
        0.08245770146815011

    REFERENCES:
    - Temme, N. M. (1975). On the numerical evaluation of the modified Bessel functions of the third kind. Journal of Computational Physics, 19(3), 324 - 337.

    """
    cdef int res = 0
    cdef double value
    res = besselk_dp_rec_c(R, x, &value, prec, pref)
    if res == 0:
        return value
    if res == 1:
        msg = f"besselk_dp_rec_c failed (k large) for x,R={x}, {R}, value={value}"
    elif res == -1:
        msg = f"Must have x > 0"
    elif res == 2:
        msg = f"besselk_dp_rec_c failed (too many iterations) for x,R={x}, {R}, value={value}"
    elif res == -2:
        msg = f"Precision requested is smaller than machine epsilon."
    else:
        msg = f"Error in besselk_dp_rec. Code: {res}"
    raise ArithmeticError(msg)

@cython.cdivision(True)
cdef int besselk_dp_rec_c(
    double R, double x, double *val, double prec=1e-14, int pref=0
) nogil:  # |r| <<1000, x >> 0 !
    r"""
    Modified K-Bessel function in double precision using the backwards Miller-recursion algorithm. 

    INPUT:

     - ``R'' -- double
     - ''x'' -- double
     - ''prec'' -- (default 1E-12) double
     - ''pref'' -- (default 1) int
                =1 => computes K_iR(x)*exp(pi*R/2)
                =0 => computes K_iR(x)

    OUTPUT:
    
     - exp(Pi*R/2)*K_{i*R}(x)  -- double

    """
    if x <= 0:
        return -1
    if sizeof(double) == 8 and prec < 2**-52:
        return -2
    cdef double p = <double>0.25 + R * R
    cdef double q = <double>2.0 * (x - d_one)
    cdef double t = 0.0  # /*arbitrary*/
    cdef double k = d_one  # /*arbitrary*/;
    cdef double err = d_one
    cdef int NMAX = 10000
    cdef int n_start = 40  # 128
    cdef double ef1 = log(d_two * x / cppi)
    cdef double mr
    cdef int n = n_start
    cdef double tmp, tmp2, ef, nr, mr_p1, mr_m1, efarg
    cdef int nn
    cdef double den
    cdef double exp_Pih_R
    cdef double err_old = <double>0
    if pref == 1:
        exp_Pih_R = exp(pihalf * R)

    for nn in range(1, NMAX+1): #from 1 <= nn <= NMAX:
        err = fabs(t - k)
        # If the error gets worse, but it is not too far off the requested precision
        # we return the value
        if (err < prec or (err_old < err and err < 100 * prec)) and n > n_start + 40:
            break
        err_old = err
        n = n + 20
        t = k
        y = d_one
        k = d_one
        d = d_one
        tmp = d_two * x - R * cppi
        nr = <double> n
        if tmp > 1300.0:
            val[0] = 0.0
            return 1
        ef = exp((ef1 + tmp) / (d_two * nr))
        mr_p1 = <double>n + 1  # m + 1
        mr = <double>n  # m
        for m in range(n, 0, -1):  # from n >=m>=1:
            mr_m1 = <double>m - 1  # m-1
            den = (q + mr_p1 * (d_two - y))
            y = (mr_m1 + p / mr) / den
            k = ef * (d + y * k)
            d = d * ef
            mr_p1 = mr
            mr = mr_m1
        if k == 0.0:
            val[0] = k
            return 1
        k = d_one / k
    if k > 1E30:
        val[0] = k
        return 1
    if nn >= NMAX:
        return 2
    if pref == 1:
        val[0] = k
    else:
        val[0] = k * exp(-pihalf * R)
    return 0

cpdef double besselk_dp_pow(double R, double x, double prec=1e-12, int pref=0):
    r"""
    Computes the modified K-Bessel function: K_iR(x) using power series.

    INPUT:

        - `R` -- parameter (double)

        - `x` -- argument (double)

        - `prec` -- precision (double)

        - `pref` -- use prefactor (integer, default 0)
               = 1 => computes K_iR(x)*exp(pi*R/2)
               = 0 => computes K_iR(x)

    OUTPUT:
        - `val`  value of the function
        - `res` -- int = 0 on success otherwise non-zero
        
    EXAMPLES::

        sage: from maass_form_core.functions.bessel.besselk_dp import besselk_dp_pow
        sage: besselk_dp_pow(10.0, 3.0) # tol 2e-14
        -6.3759939798739404e-08
        sage: besselk_dp_pow(10.0, 3.0, pref=1) # tol 2e-14
        -0.423086986725055826
        sage: besselk_dp_pow(10.0, 3.0, pref=1, prec=1e-15) # abs tol 8e-15
        -0.423086986725055826
        
    
    NOTE: The actual error might be larger than the estimated (requested) error in ``prec``
           due to loss of precision.

    REFERENCES:
    - A. Gil, J. Segura, and N. M. Temme, Evaluation of the modified Bessel function of the third kind of imaginary orders, J. Comput. Phys. 175 (2002), no. 2, 398-411.
    """
    cdef double value
    cdef int res = 0
    res = besselk_dp_pow_c(R, x, &value, prec, pref)
    if res != 0:
        if res == 2:
            msg = f"Maximum number of iterations reached x,R={x}, {R}, value = {value}"
        elif res == -2:
            msg = f"Precision smaller than machine epsilon requested"
        else:
            msg = f"Error in Besselk_dp_pow. Code: {res}"
        raise ValueError(msg)
    return value


@cython.cdivision(True)
cdef int besselk_dp_pow_c(
    double R, double x, double *val, double prec=1e-12, int pref=0
) nogil:
    r"""
    Computes the modified K-Bessel function: K_iR(x) using power series.

    INPUT:

        - `R` -- parameter (double)
        - `x` -- argument (double)
        - `val` -- pointer to result (double*)
        - `prec` -- precision (double)
        - `pref` -- use prefactor (integer, default 0) 
               = 1 => computes K_iR(x)*exp(pi*R/2)
               = 0 => computes K_iR(x)

    OUTPUT:

        - int = 0 -- if the algorithm is successful.
              =-2 -- if too large precision is requested.
              = 2 -- if too many terms are required.
        
    NOTE: The actual error might be larger than the estimated (requested) error in ``prec``
           due to loss of precision.
           
    REFERENCES:
    - Temme, N. M. (1975). On the numerical evaluation of the modified Bessel functions of the third kind. Journal of Computational Physics, 19(3), 324 - 337.
    - Gil, A., Segura, J., & Temme, N. M. (2002). Evaluation of the modified Bessel function of the third kind of imaginary orders. Journal of Computational Physics, 175(2), 398-411.


    """
    if sizeof(double) == 8 and prec < 2**-52:
        return -2
    cdef double xh = d_half * x
    cdef double xh2 = xh * xh
    cdef double complex gamma0 = loggamma_dp_c(d_one, R)
    cdef double sigma0 = cimag(gamma0)
    cdef double rsigma0 = creal(gamma0)
    cdef double th = R * log(xh) - sigma0
    cdef double tmp_sin = sin(th)
    cdef double tmp_cos = cos(th)
    cdef double tmp_factor = sqrt(cppi / (R * sinh(cppi * R)))
    cdef double f0 = -tmp_factor * tmp_sin
    cdef double Rsq = R * R
    cdef double pi_halfR = pihalf * R
    cdef double r0 = R * tmp_factor * tmp_cos
    cdef double r1 = R * tmp_factor / (d_one + Rsq) * (tmp_cos + R * tmp_sin)
    cdef double c0 = d_one
    cdef double rk1 = r0  # r(k-1)
    cdef double fk1 = f0  # f(k-1)
    cdef double summa = f0
    cdef double exp_Pih_R = <double>1.0
    cdef double fk = (fk1 + rk1) / (d_one + Rsq)
    cdef double rk = r1
    cdef double ck = xh2
    if pref == 1:
        exp_Pih_R = exp(pi_halfR)
    summa = summa + ck * fk
    fk1 = fk
    rk1 = r1
    cdef double rk2 = r0
    cdef double ck1 = ck
    cdef double test
    cdef int N_max = 1000
    cdef int k
    for k in range(2, N_max + 1):
        kk = <double>k
        den = kk * kk + Rsq
        fk = (kk * fk1 + rk1) / den
        rk = ((d_two * kk - d_one) * rk1 - rk2) / den
        ck = xh2 * ck1 / kk
        summa = summa + ck * fk
        test = ck * fk * exp_Pih_R / summa
        test = fabs(test)
        if test < prec:
            break
        fk1 = fk
        rk2 = rk1
        rk1 = rk
        ck1 = ck
    if k >= N_max:
        stat = 1
        return 2
    if pref == 1:
        val[0] = summa * exp_Pih_R
    else:
        val[0] = summa
    return 0



### Scaled Bernoulli numbers
cdef int num_ber = 50
cdef double Ber[51]  # Ber[n]=B[2n]/(2n*(2n-1))

Ber[1] =0.083333333333333333333333333333333333333333333
Ber[2] =-0.0027777777777777777777777777777777777777777778
Ber[3] =0.00079365079365079365079365079365079365079365079
Ber[4] =-0.00059523809523809523809523809523809523809523810
Ber[5] =0.00084175084175084175084175084175084175084175084
Ber[6] =-0.0019175269175269175269175269175269175269175269
Ber[7] =0.0064102564102564102564102564102564102564102564
Ber[8] =-0.029550653594771241830065359477124183006535948
Ber[9] =0.17964437236883057316493849001588939669435025
Ber[10] =-1.3924322169059011164274322169059011164274322
Ber[11] =13.402864044168391994478951000690131124913734
Ber[12] =-156.84828462600201730636513245208897382810426
Ber[13] =2193.1033333333333333333333333333333333333333
Ber[14] =-36108.771253724989357173265219242230736483610
Ber[15] =691472.26885131306710839525077567346755333407
Ber[16] =-1.5238221539407416192283364958886780518659077e7
Ber[17] =3.8290075139141414141414141414141414141414141e8
Ber[18] =-1.0882266035784391089015149165525105374729435e10
Ber[19] =3.4732028376500225225225225225225225225225225e11
Ber[20] =-1.2369602142269274454251710349271324881080979e13
Ber[21] =4.8878806479307933507581516251802290210847054e14
Ber[22] =-2.1320333960919373896975058982136838557465453e16
Ber[23] =1.0217752965257000775652876280535855003940110e18
Ber[24] =-5.3575472173300203610827709191969204484849041e19
Ber[25] =3.0615782637048834150431510513296227581941868e21
Ber[26] =-1.8999917426399204050293714293069429029473425e23
Ber[27] =1.2763374033828834149234951377697825976541634e25
Ber[28] =-9.2528471761204163072302423483476227795193312e26
Ber[29] =7.2188225951856102978360501873016379224898404e28
Ber[30] =-6.0451834059958569677431482387545472860661444e30
Ber[31] =5.4206704715700945451934778148261000136612022e32
Ber[32] =-5.1929578153140819467001947643918576846997063e34
Ber[33] =5.3036588551197005966548392430697586436992926e36
Ber[34] =-5.7633253481649640138944358507809925551907376e38
Ber[35] =6.6511557148484539375165201458105559510397394e40
Ber[36] =-8.1373783581366805387161726320935756918406892e42
Ber[37] =1.0536966953357141803754804927641810189648373e45
Ber[38] =-1.4418180599962206261805377801511812809570332e47
Ber[39] =2.0817356522089565462424808241263562311317343e49
Ber[40] =-3.1670226634886661827413495567742561342918070e51
Ber[41] =5.0700064612111373431792648153174876567629628e53
Ber[42] =-8.5299728203005518816208400522162278887807045e55
Ber[43] =1.5064172809340598576695117360379879076101931e58
Ber[44] =-2.7893494703831636871288381686312781712347569e60
Ber[45] =5.4093504352860415005763561871884152582336290e62
Ber[46] =-1.0975337821508519855016788726170795167099903e65
Ber[47] =2.3274876202618479173478641032052184930193480e67
Ber[48] =-5.1539291620653213901946552121717171770391159e69
Ber[49] =1.1906210230890226457684816176871380462750227e72
Ber[50] =-2.8668938960296673696226420541900772462913819e74

cpdef loggamma_dp(double x, double R, double prec=1e-16):
    r"""
    Logarithm (principal branch) of the Gamma function for the argument x + i*R.
    
    ALGORITHM: Stirling's formula for gamma(z+N) and the recursion gamma(z+1)=zgamma(z).
        
    INPUT:

     - ``x`` -- double
     - ``R`` -- double
     - ``prec`` -- double (default: 1e-16)  
     
    OUTPUT:
    
    - double complex -- log(Gamma(x + i*R))

    EXAMPLES::

        sage: from maass_form_core.functions.bessel.besselk_dp import loggamma_dp
        sage: loggamma_dp(1.0, 3.0) # tol 1e-14
        (-3.2441442995897556 + 1.053350771068613j)
        sage: import mpmath
        sage: a = mpmath.mpc(loggamma_dp(1.0, 3.0))
        sage: b = mpmath.loggamma(mpmath.mpc(1, 3))
        sage: abs(a - b) # abs tol 2e-15
        mpf('0.0')
        sage: abs(a - b) < 2e-15
        True

    """
    return loggamma_dp_c(x, R, prec)

@cython.cdivision(True)
cdef double complex loggamma_dp_c(double x, double R, double prec=1e-16) nogil:
    """
    Logarithm (principal branch) of the Gamma function for the argument x + i*R.
    
    ALGORITHM: Stirling's formula for gamma(z+N) and the recursion gamma(z+1)=zgamma(z).
        
    """
    cdef double a, jr, lnw, argw, Ims, Res, R2, S, r
    cdef double complex cr, tmp, res, iR, d_c_i
    cdef double complex z, w, w2, ww, summa, stirling
    cdef int i, j, N, m, M, k
    R2 = R * R
    d_c_i = (<double complex>_Complex_I)
    z = d_c_i * R + x
    y = R
    iR = d_c_i * R
    N = 15
    Nr = <double>N  # d_ten
    u = Nr + x
    v = y
    w = iR + u
    w2 = w * w  # d_100-R*R + d_20*iR  #w*w
    summa = c_zero
    ww = c_one / w
    for i in range(1, num_ber):
        tmp = (<double complex>Ber[i]) * ww
        summa = summa + tmp
        if cabs(tmp) < prec:
            break
        ww = ww / w2
    a = cabs(w)  # d_100 + R2  # cabs(w) # sqrt(u*u + y*y)
    lnw = log(a)
    argw = carg(w)
    Res = (u - d_half) * lnw - R * argw - u + ln2PI2
    Ims = R * (lnw - d_one) + (u - d_half) * argw
    stirling = <double complex>Res + (<double complex>Ims) * d_c_i
    #  this was really for logGAMMA(z + N)
    res = stirling + summa
    for j in range(N-1, -1, -1):
        cr = <double complex>j
        tmp = iR + cr + x
        res = res - clog(tmp)
    return res


@cython.cdivision(True)
cpdef besselk_real_dp(double r, double x, double eps=0, int verbose=0):
    r"""
    Modified Bessel function K_r(x) for real parameters with |r| <= 0.5.

    This function computes the modified Bessel function K_r(x) through the series expansion 
    using a recurrence relation when the parameter r satisfies |r| <= 0.5.

    INPUT:

    - ``r`` -- double, the order parameter with |r| <= 0.5
    - ``x`` -- double, the argument, must be positive
    - ``eps`` -- double (default: 0), convergence tolerance, if 0 uses machine precision
    - ``verbose`` -- int (default: 0), verbosity level for debugging output

    OUTPUT:

    The value of K_r(x) as a double.

    EXAMPLES::

        sage: from maass_form_core.functions.bessel.besselk_dp import besselk_real_dp
        sage: from sage.functions.bessel import bessel_K
        sage: # Test with r=0, should give K_0(1)
        sage: result = besselk_real_dp(0.0, 1.0)  # doctest: +ELLIPSIS
        sage: result # tol 2e-14
        0.42102443824071234
        sage: abs(result - bessel_K(0.0, 1.0)) < 2e-14
        True
        sage: # Test with small r
        sage: result = besselk_real_dp(0.25, 2.0)
        sage: result # tol 3e-13 
        0.11537827684090246
        sage: abs(result - bessel_K(0.25, 2.0)) # tol 2e-14
        0
        sage: # Test with large r
        sage: result = besselk_real_dp(0.4, 2.0)
        sage: result # tol 3e-14
        0.11772913317039947
        sage: abs(result - bessel_K(0.4, 2.0)) # tol 3e-14
        0
        sage: besselk_real_dp(0.6, 1.0) # tol 1e-13
        0.479715694892866 
        
    TESTS::
    
        sage: # Test error conditions
        sage: besselk_real_dp(1.1, 1.1)
        Traceback (most recent call last):
        ...
        ValueError: Use only for r in ]-1,1[
        sage: besselk_real_dp(0.0, 0.0)
        Traceback (most recent call last):
        ...
        ValueError: Use only for x>0

    REFERENCES:
    - Gil, A., Segura, J., & Temme, N. M. (2002). Evaluation of the modified Bessel function of the third kind of imaginary orders. Journal of Computational Physics, 175(2), 398-411.

    NOTE: This is essentially the same as besselk_dp_pow but specialized for real r in [-0.5,0.5].

    """
    if r < 0:
        r = -r
    if abs(r) >= 1:
        raise ValueError("Use only for r in ]-1,1[")
    if x <= 0.0:
        raise ValueError("Use only for x>0")
    cdef double kbval
    if r == 0:
        besselk_dp_rec_c(0.0, x, &kbval, prec=1e-15, pref=1)
        return kbval
    if eps == 0:
        eps = 2.0 ** (2 - 53)
    if eps < 2.0 ** (1 - 53):
        raise ValueError("Can not request more than double precision!")
    cdef int N = 200
    cdef double one, two, four, rpi, sqrt_two
    cdef double f0, f1
    cdef double t1, t2, t3, xtwo, xtwo_by_four
    cdef double ck0, ck1
    cdef double kk, k2
    cdef double rmax = 0.5
    cdef double ef1, ef2, ef3, err_est
    if verbose > 0:
        import logging
        logging.debug(f"r={r}, x={x}")
    cdef double r_square
    one = d_one
    two = d_two
    r_square = r ** 2
    t1 = (x / two) ** (-r)
    t2 = cppi / two / sin(cppi * r)
    t3 = t1 * creal(cexp(loggamma_dp_c(1 + r,0.0)))
    # f0 = t2 * (t1 / creal(cexp(loggamma_dp_c(1-r,0.0))) -
    #                    t1 ** -1 / creal(cexp(loggamma_dp_c(1+r,0.0))))
    f0 = (t3 / 2 / r - t2 * t3 ** -1)
    # r0 = t2 * r * (t1 / creal(cexp(loggamma_dp_c(1-r,0.0))) +
    #                    t1 ** -1 / creal(cexp(loggamma_dp_c(1+r,0.0))))
    r0 = (t3 / 2 + t2 * r * t3 ** -1)
    # r1 = t2 * r * (t1 / creal(cexp(loggamma_dp_c(2-r,0.0))) +
    #                    t1 ** -1 / creal(cexp(loggamma_dp_c(2+r,0.0))))
    t3 = t1 * creal(cexp(loggamma_dp_c(2 + r,0.0)))
    r1 = (t3 / (1 - r*r) / 2 + t2 * r * t3 ** -1)
    f1 = (f0 + r0) / (one - r * r)
    ck0 = one
    xtwo = x ** 2
    xtwo_by_four = xtwo / 4
    ck1 = xtwo_by_four
    s = ck0 * f0 + ck1 * f1
    if abs(x) < 2:
        ef1 = 2.0 * (x / 2.0) ** rmax
    else:
        ef1 = 2.0 * (2.0 / x) ** rmax
    ef1 = ef1 / (1.0 - xtwo)

    if verbose > 1:
        logging.debug(f"s= {s}")
        logging.debug(f"f0 = {f0}")
        logging.debug(f"r0 = {r0}")
        logging.debug(f"f1 = {f1}")
        logging.debug(f"r1 = {r1}")

    cdef int kmin = int(max(r + 1, xtwo_by_four)) + 1
    # for k in range(2, N + 1):
    kk = 2
    while kk <= N:
        k2 = kk ** 2
        denom = k2 - r_square
        rnew = ((2 * kk - 1) * r1 - r0 )/ denom
        fnew = (kk * f1 + r1) / denom
        ck1 = ck1 * xtwo_by_four / kk
        term = ck1 * fnew
        s += term
        if verbose > 1:
            logging.debug(f"rnew[{kk}] = {rnew}")
            logging.debug(f"fnew[{kk}] = {fnew}")
            logging.debug(f"ck[{kk}] = {ck1}")
            logging.debug(f"term= {kk}, {term}, {abs(term) / abs(s)}")
            logging.debug(f"s= {s}")

        if kk > kmin:
            # Get a rigorous error term for truncation
            ef2 = (2.0 * kk + 2.0) ** rmax
            ef3 = (3.0 / kk) ** (kk + 1)
            err_est = ef1 * ef2 * ef3 * ck1
            if verbose > 1:
                logging.debug(f"ef1,ef2,ef3,ck1={ef1},{ef2},{ef3},{ck1}")
                logging.debug(f"error est= {err_est}")
            if abs(err_est) < eps:  # abs(term)/abs(s)<eps:
                break
        r0 = r1
        r1 = rnew
        f0 = f1
        f1 = fnew
        kk += 1
    return s
