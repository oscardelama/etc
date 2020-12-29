from mpmath import mp, mpf
mp.dps = 30
from random import random

def poly_cheby(x):
    poly_cheby.p2 = 4*x*x-2
    poly_cheby.p_nm1 = x
    poly_cheby.p_nm2 = x
    yield x
    while True:
        p_np1 = poly_cheby.p2*poly_cheby.p_nm1 \
            - poly_cheby.p_nm2
        poly_cheby.p_nm2 = poly_cheby.p_nm1
        poly_cheby.p_nm1 = p_np1
        yield p_np1


SQRT2_M1 = mp.sqrt(2) - mpf(1)


def atan_coef():
    atan_coef.k = mpf(0)
    atan_coef.sqrt2_pow = mpf(-1)/SQRT2_M1

    while True:
        atan_coef.k += mpf(1)
        atan_coef.sqrt2_pow *= - SQRT2_M1*SQRT2_M1
        coef = 2/(2*atan_coef.k-1) * atan_coef.sqrt2_pow
        yield coef


def atan_approx(x):
    result = 0
    cheby = poly_cheby(x)
    coef = atan_coef()

    for i in range(37):
        result += next(coef) * next(cheby)
        
    return result

def main():
    pi_approx = atan_approx(1)*4
    print("*** Approximation to pi ****")
    print("Pi Approx: ", pi_approx)
    print("Pi Exact : ",mp.pi)
    print("Pi Approx Error: ", pi_approx-mp.pi)
    
    x = mpf(random()*2 - 1)
    approx = atan_approx(x)
    exact = mp.atan(x)
    
    print("*** atan of random number ****")
    print("x: ", x)
    print("Approx atan: ", approx)
    print("Exact atan: ", mp.atan(x))
    print("Atan Approx Error: ", approx-exact)

main()

# *** Approximation to pi ****
# Pi Approx:  3.14159265358979323846264338328
# Pi Exact :  3.14159265358979323846264338328
# Pi Approx Error:  3.54987407349455312435277854377e-30
# *** atan of random number ****
# x:  0.500529222910760296372245647945
# Approx atan:  0.464070897698497645502372515528
# Exact atan:  0.464070897698497645502372515528
# Atan Approx Error:  -3.45126646034192664867631247311e-31

