#!/usr/bin/env python

import polypy
import polypy_test

polypy_test.init()

def check_factorization(p, p_factors, expected_size):
    mul = 1
    for (factor, multiplicity) in p_factors:
        mul = mul * (factor**multiplicity)
    ok = (p == mul) and (expected_size == len(p_factors))
    if (not ok):
        print("p = {0}".format(p))
        print("p_factors = {0}".format(p_factors))
        print("expected_size = {0}".format(expected_size))
        print("mul = {0}".format(mul))
    polypy_test.check(ok)

polypy_test.start("Square-Free Factorization")

# polypy.trace_enable("polynomial")
# polypy.trace_enable("factorization")
# polypy.trace_enable("coefficient::gcd")
#
# x00 = polypy.Variable("x00")
# x01 = polypy.Variable("x01")
# x10 = polypy.Variable("x10")
# x11 = polypy.Variable("x11")
# x20 = polypy.Variable("x20")
# x21 = polypy.Variable("x21")
# x30 = polypy.Variable("x30")
# x31 = polypy.Variable("x31")
#
# polypy.variable_order.set([x31, x30, x21, x20, x11])
# p = (16*x30**4 + (32*x31**2)*x30**2 + (16*x31**4))*x20**8 + ((64*x30**4 + (128*x31**2)*x30**2 + (64*x31**4))*x21**2 + ((-64*x31)*x30**4 + (-128*x31**3)*x30**2 + (-64*x31**5))*x21 + (-32*x30**6 + (-32*x31**2)*x30**4 + (32*x31**4 - 128*x31**2)*x30**2 + (32*x31**6 - 128*x31**4)))*x20**6 + ((96*x30**4 + (192*x31**2)*x30**2 + (96*x31**4))*x21**4 + ((-192*x31)*x30**4 + (-384*x31**3)*x30**2 + (-192*x31**5))*x21**3 + (-32*x30**6 + (96*x31**2 - 128)*x30**4 + (288*x31**4 - 384*x31**2)*x30**2 + (160*x31**6 - 256*x31**4))*x21**2 + ((-64*x31)*x30**6 + (-192*x31**3 + 512*x31)*x30**4 + (-192*x31**5 + 768*x31**3)*x30**2 + (-64*x31**7 + 256*x31**5))*x21 + (16*x30**8 + (64*x31**2)*x30**6 + (96*x31**4 - 128*x31**2)*x30**4 + (64*x31**6 - 256*x31**4)*x30**2 + (16*x31**8 - 128*x31**6 + 256*x31**4)))*x20**4 + ((64*x30**4 + (128*x31**2)*x30**2 + (64*x31**4))*x21**6 + ((-192*x31)*x30**4 + (-384*x31**3)*x30**2 + (-192*x31**5))*x21**5 + (32*x30**6 + (288*x31**2 - 256)*x30**4 + (480*x31**4 - 384*x31**2)*x30**2 + (224*x31**6 - 128*x31**4))*x21**4 + ((-128*x31)*x30**6 + (-384*x31**3 + 768*x31)*x30**4 + (-384*x31**5 + 1024*x31**3)*x30**2 + (-128*x31**7 + 256*x31**5))*x21**3 + (32*x30**8 + (128*x31**2 - 128)*x30**6 + (192*x31**4 - 384*x31**2)*x30**4 + (128*x31**6 - 384*x31**4 - 512*x31**2)*x30**2 + (32*x31**8 - 128*x31**6))*x21**2)*x20**2 + ((16*x30**4 + (32*x31**2)*x30**2 + (16*x31**4))*x21**8 + ((-64*x31)*x30**4 + (-128*x31**3)*x30**2 + (-64*x31**5))*x21**7 + (32*x30**6 + (160*x31**2 - 128)*x30**4 + (224*x31**4 - 128*x31**2)*x30**2 + (96*x31**6))*x21**6 + ((-64*x31)*x30**6 + (-192*x31**3 + 256*x31)*x30**4 + (-192*x31**5 + 256*x31**3)*x30**2 + (-64*x31**7))*x21**5 + (16*x30**8 + (64*x31**2 - 128)*x30**6 + (96*x31**4 - 256*x31**2 + 256)*x30**4 + (64*x31**6 - 128*x31**4)*x30**2 + (16*x31**8))*x21**4)
# p_factors = p.factor_square_free()

x = polypy.Variable("x");
y = polypy.Variable("y");
z = polypy.Variable("z");

polypy.variable_order.set([z, y, x])

# polynomials and expected factors
test1 = [(x + 1 - x, 0),
         (x + 3 - x, 1),
         (x + 1 - 1, 1),
         (2*x      , 2),
         (x + 1    , 1),
         (2*x + 2  , 2),
         (2*x + 1  , 1)]

for (p, expected) in test1:
    p_factors = p.factor_square_free()
    check_factorization(p, p_factors, expected)

test2 = [(x*y        , 2),
         (x*y + y    , 2),
         (x*y + 1    , 1),
         (2*x*y + 1  , 1),
         (2*x*y + 2  , 2),
         (2*x*y + 2*y, 3)]

for (p, expected) in test2:
    p_factors = p.factor_square_free()
    check_factorization(p, p_factors, expected)

test3 = [( (x + 1)*(x + 2)**2*(x + 3)**3, 3),
         ( (x + 1)*(y + 2)**2*(z + 3)**3, 3),
         ( (z + 1)*(y + 2)**2*(x + 3)**3, 3)]

for (p, expected) in test3:
    p_factors = p.factor_square_free()
    check_factorization(p, p_factors, expected)
