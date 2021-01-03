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


def atan_coef():
    SQRT2_M1 = mp.sqrt(2) - mpf(1)
    SQRT2_M1_2 = SQRT2_M1*SQRT2_M1
    atan_coef.k = mpf(0)
    atan_coef.rho_factor = mpf(-2)/SQRT2_M1

    while True:
        atan_coef.k += mpf(1)
        atan_coef.rho_factor *= - SQRT2_M1_2
        coef = atan_coef.rho_factor/(2*atan_coef.k-1)
        yield coef


def main():
    result = 0
    cheby = atan_poly_cheb(1)
    coef = atan_coef()

    for i in range(37):
        result += next(coef) * next(cheby)
    
    print("*** Approximation to pi ****")
    print("Pi Approx: ", result)
    print("Pi Exact : ",mp.pi)
    print("Pi Approx Error: ", result-mp.pi)
    
main()

# *** Approximation to pi ****
# Pi Approx:  0.785398163397448309615660845821
# Pi Exact :  3.14159265358979323846264338328
# Pi Approx Error:  -2.35619449019234492884698253746



