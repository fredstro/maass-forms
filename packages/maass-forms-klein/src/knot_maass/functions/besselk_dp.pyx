# cython: profile=True
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
#*****************************************************************************
# from psage.rings.mp_cimports cimport *
from libc.stdio cimport *
import cython

r"""
Low-precision (fast) algorithms for the K-Bessel function.
Also algorithms for incomplete gamma function.

    AUTHOR:
    
    - Fredrik Strömberg (April 2010)


    EXAMPLES::

        sage: from knot_maass.functions.besselk_dp import besselk_dp, besselk_dp_pow
        sage: from knot_maass.functions.besselk_dp import besselk_dp_rec
        sage: from sage.all import RealField, ComplexField
        sage: besselk_dp(10.0,5.0) # tol 1e-14
        -0.7183327166568183
        sage: besselk_dp_rec(10.0,3.0) # tol 1e-14
        -6.3759939798738967e-08
        sage: besselk_dp_rec(10.0,3.0,pref=1) # tol 1e-14
        -0.42308698672505796
        sage: besselk_dp_pow(10.0,3.0) # tol 1e-14
        -6.3759939798739404e-08
        sage: besselk_dp_pow(10.0,3.0,pref=1) # tol 1e-14
        -0.42308698672506084
        sage: a=besselk_dp_pow(10.0,3.0,pref=1);a  # tol 1e-14
        -0.42308698672506084
        sage: b=besselk_dp_rec(10.0,3.0,pref=1);b  # tol 1e-14
        -0.42308698672505796
        sage: abs(a-b) # tol 6e-15
        0
        sage: besselk_dp_rec(100.0,3.0,pref=1)  # tol 3e-13
        0.082457701468151401
        sage: RF=RealField(150)
        sage: CF=ComplexField(150)
        sage: bessel_K(CF(100*I),RF(3))*exp(RF.pi()*RF(50))  # tol 1e-14
        0.0824577014681303
        sage: _-besselk_dp_rec(100.0,3.0,pref=1)  # tol 5e-14
        0


"""

cdef extern from "complex.h":
    cdef double complex _Complex_I

cdef complex CMPLXF(float x,float y):
    return x+_Complex_I*y

cdef double complex CMPLX(double x,double y):
    return x+_Complex_I*y
    
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
    double pow(double,double)
    double complex clog(double complex)
    double cimag(double complex)
    double creal(double complex)
    double carg(double complex)
    

cdef double cppi  =<double> 3.14159265358979323846264338327950288419716939 # =pi
cdef double pihalf=<double> 1.5707963267948966192313216916397514420985847  # =pi/2
cdef double ln2PI2= <double> 0.91893853320467274178032973640561763986139747 #=ln(2pi)/2
cdef double d_one =<double> 1.0
cdef double d_half =<double> 0.5
cdef double d_two =<double> 2.0
cdef double d_ten =<double> 10.0
cdef double d_20 =<double> 20.0
cdef double d_100 =<double> 100.0
cdef double complex c_one = CMPLX(1.0,0.0)
cdef double complex c_zero = CMPLX(0.0,0.0)


cpdef besselk_dp(double R,double x,double prec=1e-14,int pref=1,algorithm='default'):
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

        sage: from knot_maass.functions.besselk_dp import besselk_dp
        sage: besselk_dp(10.0,5.0) # tol 2e-14
        -0.7183327166568183
        sage: besselk_dp(10.0,3.0,pref=0) # tol 2e-14
        -6.3759939798738967e-08
        sage: besselk_dp(10.0,3.0,pref=1) # tol 2e-14
        -0.42308698672505796
        

    """
    RR = fabs(R)
    if x <= 0:
        raise ValueError(f" Need x>0! Got x={x}")
    if pref == -1:
        return besselk_real_rec_dp(R,x,prec)
    cdef double xc,xcral,S
    xc = sqrt((x + R) * (x - R))
    S=R / x
    xcral=xc + R * 2.0 * atan(S / (1.0 + (xc / x)))
    cdef double kbes
    if R * pihalf - xcral < -125.0:
        return 0.0
    cdef int res = -1
    if algorithm != 'default':
        if algorithm == 'pow':
            res=besselk_dp_pow_c(RR, x, &kbes, prec, pref)
        elif algorithm == 'rec':
            res=besselk_dp_rec_c(RR, x, &kbes, prec, pref)
        if res == 1:
            raise ValueError,'The K-bessel routine failed (k large) for x,R={0},{1}, value={2}'.format(x,RR,kbes)
        elif res == 2:
            raise ValueError,'The K-bessel routine failed (too many iterations)  for x,R={0},{1}, value={2}'.format(x,RR,kbes)
        return kbes
    if x < R * 0.7:
        res=besselk_dp_pow_c(RR,x,&kbes,prec,pref)
        if res != 0:
            res=besselk_dp_rec_c(RR,x,&kbes,prec,pref)
    else:
        res=besselk_dp_rec_c(RR,x,&kbes,prec,pref)
    if res == 1:
        raise ValueError(f"The K-bessel routine failed (k large) for x,R={x},{RR}, value={kbes}")
    elif res==2:
        raise ValueError("The K-bessel routine failed (too many iterations)  for x,R={x},{RR}, "
                         "value={kbes}")
    elif res==0:
        return kbes
    else:
        raise ValueError(f"The K-bessel routine failed (unknown error) for x,R={x},{RR}, "
                         f"value={kbes}")

cdef int besselk_dp_c(double *kbes,double R,double x,double prec,int pref,int verbose=0) nogil: #double prec=1e-14,int pref=0):
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
            printf('Need x>0! Got x=%f',x)
        return -1
    cdef int res=0
    if x > 0.5*cppi*RR-2.0*log(prec):
        kbes[0] = 0.0
        return 0
    if x<R*0.7:
        res=besselk_dp_pow_c(RR,x,kbes,prec,pref)
        if res != 0:
            res=besselk_dp_rec_c(RR,x,kbes,prec,pref)
    else:
        res = besselk_dp_rec_c(RR,x,kbes,prec,pref)

    if res == 1:
        printf('The K-bessel routine failed (k large) for x,R=%g,%g, value=%g \n',x,RR,kbes[0])
    elif res==2:
        printf('the K-bessel routine failed (too many iterations) for x,R=%g,%g, value=%g',x,RR,kbes[0])
    return res

@cython.cdivision(True)
cpdef double besselk_dp_rec(double R, double x, double prec=1e-14,
                        int pref=0):  # |r| <<1000, x >> 0 !
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

        sage: from knot_maass.functions.besselk_dp import besselk_dp_rec
        sage: besselk_dp_rec(10.0,5.0,prec=1e-16,pref=1) # tol 2e-13
        -0.7183327166568183
        sage: besselk_dp_rec(10.0,3.0,prec=1e-16) # tol 2e-15
        -6.3759939798738967e-08
        sage: besselk_dp_rec(10.0,3.0,prec=1e-16,pref=1) # tol 3e-15
        -0.42308698672505796
        sage: besselk_dp_rec(100.0,3.0,prec=1e-16,pref=1) # tol 4e-13
        0.08245770146816264
        sage: besselk_dp_rec(100.0,3.0,prec=1e-16,pref=1) # tol 4e-13
        0.0824577014681302


    """
    cdef res = 0
    cdef double value
    res = besselk_dp_rec_c(R, x, &value, prec, pref)
    if res != 0:
        if res == 1:
            msg = f"besselk_dp_rec_c failed (k large) for x,R={x}, {R}, value={value}"
        elif res == -1:
            msg = f"Must have x > 0"
        elif res == 2:
            msg = f"besselk_dp_rec_c failed (too many iterations) for x,R={x}, {R}, value={value}"
        else:
            msg = f"Error in besselk_dp_rec. Code: {res}"
        raise ArithmeticError(msg)
    return value

@cython.cdivision(True) 
cdef int besselk_dp_rec_c(double R,double x,double *val,double prec=1e-14,int pref=0) nogil:   # |r| <<1000, x >> 0 !
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

    """
    if x < 0:
        return -1
    if sizeof(double) == 8 and prec < 2**-52:
        return -2
    cdef double p = <double>0.25 + R * R
    cdef double q = <double>2.0 * (x - d_one)
    cdef double t = 0.0 #/*arbitrary*/
    cdef double k = d_one #/*arbitrary*/;
    cdef double err = d_one
    cdef int NMAX = 10000
    cdef int n_start = 40 #128
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
        #printf('n=%g, err=%g, k=%g\n', <double>nn, err, k)
        # If the error gets worse, but it is not too far off the requested precision
        # we return the value
        if (err < prec or (err_old < err and err < 100*prec)) and n > n_start + 40:
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
        ef = exp( (ef1 + tmp) / (d_two * nr))
        mr_p1 = <double> n+1   # m+1
        mr = <double> n        # m
        for m in range(n, 0, -1): #from n >=m>=1:
            mr_m1 = <double> m-1 # m-1
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
        val[0] =  k * exp(-pihalf * R)
    return 0

cpdef double besselk_dp_pow(double R,double x, double prec=1E-12,int pref=0):
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

        sage: from knot_maass.functions.besselk_dp import besselk_dp_pow
        sage: besselk_dp_pow(10.0,3.0) # tol 7e-15
        -6.3759939798739404e-08
        
        # Value that is expected to be returned
        
        sage: besselk_dp_pow(10.0,3.0,pref=1) # tol 1e-13
        -0.42308698672506334
        
        # Compared to correct value

        sage: besselk_dp_pow(10.0,3.0,pref=1, prec=1e-15) # abs tol 8e-15 
        -0.42308698672505557
        
    
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
cdef int besselk_dp_pow_c(double R, double x, double *val, double prec=1E-12, int pref=0) nogil:
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
    """
    if sizeof(double) == 8 and prec < 2 ** -52:
        return -2
    cdef double xh=d_half*x
    cdef double xh2=xh*xh
    cdef double complex gamma0 = loggamma_dp_c(d_one,R)
    cdef double sigma0 = cimag(gamma0)
    cdef double rsigma0 = creal(gamma0)
    cdef double th=R * log(xh) - sigma0
    cdef double tmp_sin = sin(th)
    cdef double tmp_cos = cos(th)      
    cdef double tmp_factor=sqrt(cppi / (R*sinh(cppi*R)))
    cdef double f0 = -tmp_factor * tmp_sin
    cdef double Rsq = R * R
    cdef double pi_halfR = pihalf * R
    cdef double r0=R * tmp_factor * tmp_cos
    cdef double r1=R * tmp_factor / (d_one+Rsq) * (tmp_cos + R * tmp_sin)
    cdef double c0=d_one
    cdef double rk1=r0  # r(k-1)
    cdef double fk1=f0  # f(k-1)
    cdef double summa=f0
    cdef double exp_Pih_R = <double>1.0
    cdef double fk=(fk1+rk1)/(d_one+Rsq)
    cdef double rk=r1
    cdef double ck=xh2
    if pref == 1:
        exp_Pih_R = exp(pi_halfR)
    summa = summa + ck * fk
    fk1 = fk
    rk1 = r1
    cdef double rk2=r0
    cdef double ck1=ck
    cdef double test
    cdef int N_max=1000
    cdef int k
    for k in range(2, N_max+1): #from  2<=k<=N_max:
        kk = <double> k
        den = kk * kk + Rsq
        fk = (kk * fk1 + rk1) / den
        rk = ((d_two * kk - d_one) * rk1 - rk2) / den
        ck = xh2 * ck1 / kk
        summa = summa + ck * fk
        test = ck * fk * exp_Pih_R / summa
        test = fabs(test)
        if test < prec:
            break
        fk1=fk
        rk2=rk1
        rk1=rk
        ck1=ck
    if k >= N_max:
        stat=1  #! We reached end of loop
        return 2
    if pref == 1:
        val[0] = summa * exp_Pih_R
    else:
        val[0] = summa
    return 0



### Scaled Bernoulli numbers
cdef int num_ber=50
cdef double Ber[51] # Ber[n]=B[2n]/(2n*(2n-1))

Ber[ 1 ] = 0.083333333333333333333333333333333333333333333
Ber[ 2 ] = -0.0027777777777777777777777777777777777777777778
Ber[ 3 ] = 0.00079365079365079365079365079365079365079365079
Ber[ 4 ] = -0.00059523809523809523809523809523809523809523810
Ber[ 5 ] = 0.00084175084175084175084175084175084175084175084
Ber[ 6 ] = -0.0019175269175269175269175269175269175269175269
Ber[ 7 ] = 0.0064102564102564102564102564102564102564102564
Ber[ 8 ] = -0.029550653594771241830065359477124183006535948
Ber[ 9 ] = 0.17964437236883057316493849001588939669435025
Ber[ 10 ] = -1.3924322169059011164274322169059011164274322
Ber[ 11 ] = 13.402864044168391994478951000690131124913734
Ber[ 12 ] = -156.84828462600201730636513245208897382810426
Ber[ 13 ] = 2193.1033333333333333333333333333333333333333
Ber[ 14 ] = -36108.771253724989357173265219242230736483610
Ber[ 15 ] = 691472.26885131306710839525077567346755333407
Ber[ 16 ] = -1.5238221539407416192283364958886780518659077e7
Ber[ 17 ] = 3.8290075139141414141414141414141414141414141e8
Ber[ 18 ] = -1.0882266035784391089015149165525105374729435e10
Ber[ 19 ] = 3.4732028376500225225225225225225225225225225e11
Ber[ 20 ] = -1.2369602142269274454251710349271324881080979e13
Ber[ 21 ] = 4.8878806479307933507581516251802290210847054e14
Ber[ 22 ] = -2.1320333960919373896975058982136838557465453e16
Ber[ 23 ] = 1.0217752965257000775652876280535855003940110e18
Ber[ 24 ] = -5.3575472173300203610827709191969204484849041e19
Ber[ 25 ] = 3.0615782637048834150431510513296227581941868e21
Ber[ 26 ] = -1.8999917426399204050293714293069429029473425e23
Ber[ 27 ] = 1.2763374033828834149234951377697825976541634e25
Ber[ 28 ] = -9.2528471761204163072302423483476227795193312e26
Ber[ 29 ] = 7.2188225951856102978360501873016379224898404e28
Ber[ 30 ] = -6.0451834059958569677431482387545472860661444e30
Ber[ 31 ] = 5.4206704715700945451934778148261000136612022e32
Ber[ 32 ] = -5.1929578153140819467001947643918576846997063e34
Ber[ 33 ] = 5.3036588551197005966548392430697586436992926e36
Ber[ 34 ] = -5.7633253481649640138944358507809925551907376e38
Ber[ 35 ] = 6.6511557148484539375165201458105559510397394e40
Ber[ 36 ] = -8.1373783581366805387161726320935756918406892e42
Ber[ 37 ] = 1.0536966953357141803754804927641810189648373e45
Ber[ 38 ] = -1.4418180599962206261805377801511812809570332e47
Ber[ 39 ] = 2.0817356522089565462424808241263562311317343e49
Ber[ 40 ] = -3.1670226634886661827413495567742561342918070e51
Ber[ 41 ] = 5.0700064612111373431792648153174876567629628e53
Ber[ 42 ] = -8.5299728203005518816208400522162278887807045e55
Ber[ 43 ] = 1.5064172809340598576695117360379879076101931e58
Ber[ 44 ] = -2.7893494703831636871288381686312781712347569e60
Ber[ 45 ] = 5.4093504352860415005763561871884152582336290e62
Ber[ 46 ] = -1.0975337821508519855016788726170795167099903e65
Ber[ 47 ] = 2.3274876202618479173478641032052184930193480e67
Ber[ 48 ] = -5.1539291620653213901946552121717171770391159e69
Ber[ 49 ] = 1.1906210230890226457684816176871380462750227e72
Ber[ 50 ] = -2.8668938960296673696226420541900772462913819e74

cpdef loggamma_dp(double x, double R, double prec=1E-16):
    r"""
    Logarithm of Gamma function for the argument x+i*R
    (using principal branch of the logarithm).
    
    INPUT:
    
     - ''x'' -- double 
     - ''R'' -- double 
     
    OUTPUT:
    
    - double complex -- log(Gamma(x+i*R))

    EXAMPLES::

        sage: from knot_maass.functions.besselk_dp import loggamma_dp
        sage: loggamma_dp(1.0, 3.0) # tol 5e-16
        (-3.2441442995897556+1.053350771068613j) 
        sage: import mpmath
        sage: a=mpmath.mpc(loggamma_dp(1.0, 3.0))
        sage: b=mpmath.loggamma(mpmath.mpc(1,3))
        sage: abs(a-b) # abs tol 3e-15
        mpf('0.0')
        sage: abs(a-b) < 3e-15
        True

    """
    return loggamma_dp_c(x,R,prec)

@cython.cdivision(True) 
cdef double complex loggamma_dp_c(double x, double R,double prec=1E-16) nogil:
    """
    Logarithm of Gamma function for the argument x+i*R
    (using principal branch of the logarithm).
    
    """
    cdef double a,jr,lnw,argw,Ims,Res,R2,S,r
    cdef double complex cr,tmp, res,iR,d_c_i
    cdef double complex z,w ,w2,ww, summa , stirling 
    cdef int i,j ,N ,m,M,k
    R2 = R * R
    d_c_i = (<double complex>_Complex_I)
    z = d_c_i * R + x
    y = R
    iR = d_c_i*R
    N = 11
    Nr = <double> N # d_ten
    u = Nr + x
    v = y
    w = iR + u
    w2 = w * w   #d_100-R*R+d_20*iR  #w*w
    summa = c_zero
    ww = c_one / w
    #ns=50
    for i in range(1,num_ber):
        tmp = (<double complex>Ber[i]) * ww
        summa = summa + tmp
        if cabs(tmp) < prec:
            break
        ww = ww / w2
    a = cabs(w) #d_100+R2  # cabs(w) # sqrt(u*u+y*y)
    lnw = log(a)
    argw = carg(w)
    Res = (u - d_half) * lnw - R * argw - u + ln2PI2
    Ims = R * (lnw - d_one) + (u - d_half) * argw
    stirling = <double complex>Res + (<double complex>Ims) * d_c_i
    #  this was really for GAMMA(z+N)
    res = stirling + summa
    for j in range(N): #from 0 <=j<=N-1:
        cr = <double complex> j
        tmp = iR + cr + x #d_one
        res = res - clog(tmp)
    return res



## Bessel function for real parameter |r|<=0.5
cpdef besselk_real_rec_dp(double r,double x,double eps=0,int verbose=0):
    r"""
    K_r(x) with x real, positive and r real with |r|<=0.5
    """
    if r<0:
        r = -r
    if abs(r)>0.5:
        raise ValueError,"Use only for r in [-0.5,0.5]"
    if abs(x)==0.0:
        raise ValueError,"Use only for x>0"
    if eps == 0:
        eps = 2.0**(2-53)
#    cdef dict rk,fk,ck
    cdef int N=100
    cdef double one,two,pi,four,rpi
    cdef double t1,f11,f111,xtwo,xtwo_by_four
    cdef double rk0,rk1,rk_new
    cdef double fk0,fk1,fk_new
    cdef double ck0,ck1,ck_new
    cdef double R2
    cdef double kk,k2
    #    rk = {}; fk={}; ck={}
    one = 1.0; two = 2.0
    pi = cppi
    four=4.0
    k2=1.0; kk=1.0
    t1 = exp( creal(loggamma_dp_c(r,0.0)))
    t1*=(x/two)**(-r)
    rpi=r*pi
    f1 = sin(rpi)/pi    
    f11= f1*t1
    f111=r*f11
    f2 = t1**-1
    rk0 = f111+f2
    fk0 = f11 - f2/r
    rk1 = f111/(one-r) + f2/(one+r)
    R2 = r*r
    fk1 = (fk0 + rk0)/(one-R2)
    s = 0.0
    xtwo = x**2
    xtwo_by_four = xtwo/four
    ck0 = one
    ck1 = xtwo_by_four
    s = ck0*fk0+ck1*fk1
    if eps < 2.0**(1-53):
        return 1 #raise ValueError,"Need higher precision input"
    ##
    cdef double rmax = 0.5
    if verbose>1:
        printf("s= %f",s)
        printf("r[0] = %f",rk0)
        printf("r[1] = %f",rk1)
        printf("f[0] = %f",fk0)
        printf("f[1] = %f",fk1)
        printf("c[0] = %f",ck0)
        printf("c[1] = %f",ck1)
    cdef int kmin,k

    cdef double ef1,ef2,ef3,err_est
    if abs(x)<2:
        ef1 = 2.0*(<double>x/2.0)**rmax
    else:
        ef1 = 2.0*(2.0/x)**rmax
    ef1=ef1/(1.0-xtwo)
    
    kmin = int( max(r+1,xtwo_by_four))+1
    # fk1 = fk[k-1]
    for k in range(2,N+1):
        kk = <double>k
        k2 = kk**2
        fk1 = (kk*fk1 + rk1)/(k2-R2)
        ck1 = ck1*xtwo_by_four/kk
        term = ck1*fk1
        if verbose>1:
            printf("f[%d] = %f",k,fk1)
            printf("c[%d] = %f",k,ck1)
        s+=term
        if k > kmin:
            # Get a rigorous error term for truncation
            ef2 = (2.0*k+2.0)**rmax
            ef3 = (3.0/kk)**(k+1)
            err_est = ef1*ef2*ef3*ck1
            # Also add numerical error
#            err_est+=k*meps
            if verbose>0:
                printf("term= %d, %f, %f",k,term,abs(term)/abs(s))
                printf("error est= %f",err_est)
                printf("s= %f",s)
            if abs(err_est)< eps: #abs(term)/abs(s)<eps:
                break
        rk_new = ((2*kk-1)*rk1 - rk0)/(k2 - R2)
        rk0 = rk1
        rk1 = rk_new
        if verbose>2:
            printf("r[%d] = %f",k,rk1)
    return s*cppi/two
