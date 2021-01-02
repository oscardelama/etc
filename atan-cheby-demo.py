from mpmath import mp, mpf
mp.dps = 30

def atan_poly_cheb(x):
    atan_poly_cheb.p2 = 4*x*x-2
    atan_poly_cheb.p_nm1 = x
    atan_poly_cheb.p_nm2 = x
    yield x
    while True:
        p_np1 = atan_poly_cheb.p2*atan_poly_cheb.p_nm1 \
            - atan_poly_cheb.p_nm2
        atan_poly_cheb.p_nm2 = atan_poly_cheb.p_nm1
        atan_poly_cheb.p_nm1 = p_np1
        yield p_np1


def atan_cheby_coef():
    SQRT2_M1 = mp.sqrt(2) - mpf(1)
    SQRT2_M1_2 = SQRT2_M1*SQRT2_M1
    atan_cheby_coef.k = mpf(0)
    atan_cheby_coef.rho_factor = mpf(-2)/SQRT2_M1

    while True:
        atan_cheby_coef.k += mpf(1)
        atan_cheby_coef.rho_factor *= - SQRT2_M1_2
        coef = atan_cheby_coef.rho_factor/(2*atan_cheby_coef.k-1)
        yield coef


def atan_approx(x):
    result = 0
    cheby = atan_poly_cheb(x)
    coef = atan_cheby_coef()

    for i in range(37):
        result += next(coef) * next(cheby)
        
    return result

def main():
    pi_approx = atan_approx(1)*4
    print("*** Approximation to pi ****")
    print("Pi Approx: ", pi_approx)
    print("Pi Exact : ",mp.pi)
    print("Pi Approx Error: ", pi_approx-mp.pi)
    
    from random import random
    x = mpf(random()*2 - 1)
    approx = atan_approx(x)
    exact = mp.atan(x)
    
    print("\n*** atan of random number ****")
    print("x: ", x)
    print("Approx atan: ", approx)
    print("Exact atan : ", mp.atan(x))
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

