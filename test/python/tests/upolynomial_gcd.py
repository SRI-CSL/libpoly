#!/usr/bin/env python

import polypy
import polypy_test

polypy_test.init()

def check_gcd_single(p, q):
    global polypy_test
    gcd = p.gcd(q)
    ok = polypy_test.sympy_checker.check_gcd(p, q, gcd)
    polypy_test.check(ok)

def check_gcd(p, q):
    check_gcd_single(p, q)
    check_gcd_single(p, -q)
    check_gcd_single(-p, q)
    check_gcd_single(-p, -q)

def check_extended_gcd_single(p, q):
    global polypy_test
    (gcd, u, v) = p.extended_gcd(q)
    ok = polypy_test.sympy_checker.check_extended_gcd(p, q, gcd, u, v)
    polypy_test.check(ok)

def check_extended_gcd(p, q):
    check_extended_gcd_single(p, q)
    check_extended_gcd_single(p, -q)
    check_extended_gcd_single(-p, q)
    check_extended_gcd_single(-p, -q)
            
polypy_test.start("GCD in Z_13 (Knuth)")
        
# Example from Knuth (p 424)
K = polypy.CoefficientRing(13)
x = polypy.x.to_ring(K);

p = x**8 + x**6 + 10*x**4 + 10*x**3 + 8*x**2 + 2*x + 8
q = 3*x**6 + 5*x**4 + 9*x**2 + 4*x + 8
check_gcd(p, q)

polypy_test.start("Regressions (modular)")

K = polypy.CoefficientRing(7);
x = polypy.x.to_ring(K)

p = 2*x
q = 4*x
check_gcd(p, q)

p = 1*x**4 + 2*x**2 + (-2*x)
q = 1*x**4 + 1*x**2 + 1*x + 3
check_gcd(p, q)

p = (-2*x**2) + 1*x
q = (-2*x**2) + (-1*x)
check_gcd(p, q)

p = x
q = x + 1
check_extended_gcd(p, q)

p = 1*x**2 + 1*x + (-2)
q = 1*x**2 + (-3*x) + 2
check_extended_gcd(p, q)

p = 1*x**3 + (-2*x**2) + 1*x + (-1) 
q = 1*x**3 + (-1*x**2) + (-1*x) + (-2)
check_extended_gcd(p, q)

K = polypy.CoefficientRing(11);
x = polypy.x.to_ring(K)

p = 1*x**6 + (-3*x**5) + 3*x**4 + (-4*x**3) + 1*x**2 + (-1)
q = 1*x**6 + (-3*x**5) + 3*x**4 + (-5*x**3) + (-4*x**2) + (-5*x) + (-2)
check_gcd(p, q)

polypy_test.start("Random (modular)")

for prime in [7, 11, 13]:
    K = polypy.CoefficientRing(prime)
    for d in range(1, 10):
        gcd_gold = polypy_test.random_upolynomial(K, d, 20, 1)
        tmp = polypy_test.random_upolynomial(K, d, 20);
        p = tmp*gcd_gold
        q = (tmp-1)*gcd_gold
        check_gcd(p, q)
    for d in range(1, 10):
        p = polypy_test.random_upolynomial(K, d, 20, 20)
        q = polypy_test.random_upolynomial(K, d, 20, 20)
        check_gcd(p, q);

polypy_test.start("Random extended (modular)")
 
for prime in [7, 11, 13]:
    K = polypy.CoefficientRing(prime)
    for d in range(1, 10):
        p = polypy_test.random_upolynomial(K, d, 20, 20)        
        q = polypy_test.random_upolynomial(K, d, 20, 20)
        d = polypy_test.random_upolynomial(K, 1, 20, 20)
        check_extended_gcd(p*d, q*d);

polypy_test.start("Regressions (Z)")

x = polypy.x

p = 5*x**6 + (-89*x**5) + 297*x**4 + 91*x**3 + (-298*x**2) + 12*x + 10
q = 5*x**6 + (-89*x**5) + 297*x**4 + 90*x**3 + (-284*x**2) + 2*x + 8
check_gcd(p, q)

p = x**2 + 2*x + (-63)
q = x**2 + 1*x + (-72)
check_gcd(p, q)

check_gcd(6*p, 9*q)

p = 3*(x + 1)
q = 6 + x - x
check_gcd(p, q)

polypy_test.start("Random (Z)")

for d in range(1, 10):
    gcd_gold = polypy_test.random_upolynomial(polypy.Z, d, 20)
    tmp = polypy_test.random_upolynomial(polypy.Z, d, 20);
    p = tmp*gcd_gold
    q = (tmp-1)*gcd_gold
    check_gcd(p, q)
for d in range(1, 10):
    p = polypy_test.random_upolynomial(polypy.Z, d, 20, 20)
    q = polypy_test.random_upolynomial(polypy.Z, d, 20, 20)
    check_gcd(p, q);
    
