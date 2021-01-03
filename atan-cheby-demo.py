from mpmath import mp, mpf
mp.dps = 30

def poly_cheb_odd(x):
    poly_cheb_odd.p2 = 4*x*x-2
    poly_cheb_odd.p_nm1 = x
    poly_cheb_odd.p_nm2 = x
    yield x
    while True:
        p_np1 = poly_cheb_odd.p2*poly_cheb_odd.p_nm1 \
            - poly_cheb_odd.p_nm2
        poly_cheb_odd.p_nm2 = poly_cheb_odd.p_nm1
        poly_cheb_odd.p_nm1 = p_np1
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


def main():
    result = 0
    cheby = poly_cheb_odd(1)
    coef = atan_cheby_coef()

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

