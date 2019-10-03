#!/usr/bin/env python

import polypy
import polypy_test

import random

polypy_test.init()

def check_factorization(p, debug=False, check_product = True):
    global polypy_test
    if debug:
        print("check_factorization:")
        print("p = {0}".format(p))
    # Do the factorization
    factorization = p.factor()
    C, factors = factorization[0], factorization[1:]
    if debug:
        print("C = {0}".format(C))
        print("factors = {0}".format(factors))
    # The constant should always be positive
    if C <= 0:
        print("Wrong factorization (constant)")
        print("p = {0}".format(p))
        print("C = {0}".format(C))
        print("factors = {0}".format(factors))
        polypy_test.check(False)
        return
    # Check if factorization multiplies to the input
    if check_product:
        product = C
        for (f, d) in factors:
            if (f.degree() == 0):
                print("Wrong factorization (constant factor)")
                print("p = {0}".format(p))
                print("factors = {0}".format(factors))
                polypy_test.check(False)
            product = product * f**d
        if (p != product):
            print("Wrong factorization (product mismatch)")
            print("p       =  {0}".format(p))
            print("product =  {0}".format(product))
            print("factors =  {0}".format(factors))
            polypy_test.check(False)
            return
    # Done, we're OK
    polypy_test.check(True)


"""
 Returns a list of cyclotomic polynomials F1...Fn
"""
def cyclotomic(n):
    if (n == 1):
        P1 = polypy.x - 1
        return [P1]
    else:
        # Get up to 1 polynomials
        L = cyclotomic(n-1)
        # Compute Pn
        Pn = polypy.x**n - 1
        for d, Pd in enumerate(L):
            if (n % (d+1) == 0):
                Pn = Pn / Pd
        L.append(Pn)
        return L

polypy_test.start("Factorization in Z (Regressions)")

# polypy.trace_enable("factorization")
# polypy.trace_enable("hensel")
# polypy.trace_enable("arithmetic")

x = polypy.x

p = (x - 1)*(x - 2)*(x - 3)*(x - 4)*(x - 5);
check_factorization(p)

p = 1*x**4 + 94*x**3 + 107*x**2 + 771*x + (-690)
check_factorization(p)

p= x**4 + 1
check_factorization(p)

p = x**8 - x**6 + x**4 - x**2 + 1
check_factorization(p)

p = (x**2 + 1)*(x**2 + 3)
check_factorization(p)

p = x**2 + 1
check_factorization(p)

p = (x**2 + 1)*(x**2 + 2)
check_factorization(p)

p = (x - 1)*(x**2 + 1)
check_factorization(p)

p = (x - 1)*(x - 2)
check_factorization(p)

polypy_test.start("Berlekamp in Z_13 (Knuth)")

K = polypy.CoefficientRing(13)
x = polypy.x.to_ring(K);

p = x**8 + x**6 + 10*x**4 + 10*x**3 + 8*x**2 + 2*x + 8

check_factorization(p)

polypy_test.start("Square-free modular (simple)")

K = polypy.CoefficientRing(13)
x = polypy.x.to_ring(K);

p = x**13
check_factorization(p);

p = x**13 - x
check_factorization(p)

p = x*(x-1)**2*(x-2)**3
check_factorization(p)

p = 1*x**6 + 3*x**5 + (-5*x**4) + 4*x**3 + 6*x**2 + 6*x + (-1)
check_factorization(p)

polypy_test.start("Square-free modular (random)")

for factors_count in range(1, 11):
    for test in range(1, 10):
        # Constant
        p = 1
        # Factors
        deg = random.randint(1, 2)
        for k in range(0, factors_count):
            degree =  random.randint(1, 3)
            f_k = polypy_test.random_upolynomial(K, degree = 3, M = 10, lc = 1)
            p *= f_k**deg
            deg += random.randint(1, 2)
        # Check it
        check_factorization(p)

polypy_test.start("Cyclotomic (Z)")

#for p in cyclotomic(100):
#    print(p)

# http://mathworld.wolfram.com/Swinnerton-DyerPolynomial.html
polypy_test.start("Swinnerton-Dyer (Z)")
