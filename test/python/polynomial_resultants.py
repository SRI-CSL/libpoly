#!/usr/bin/python

import polypy
import polypy_test

x = polypy.Variable("x");
y = polypy.Variable("y");
z = polypy.Variable("z"); 
w = polypy.Variable("w");

polypy.variable_order.set([w, z, y, x]);

def check_psc(p, q, expected):
    psc = p.psc(q)
    ok = cmp(psc, expected) == 0
    if (not ok):
        print "p =", p
        print "q =", q
        print "psc      =", psc
        print "expected =", expected
    polypy_test.check(ok)

def check_resultant(p, q, expected):
    resultant = p.resultant(q)
    ok = resultant == expected
    if (not ok):
        print "p =", p
        print "q =", q
        print "resultant =", resultant
        print "expected  =", expected

polypy_test.init()

polypy_test.start("Speed")

A = w - 2*x**3*y**2*z - x**2*y*z**3 - x*y**3*z**2 - z - 4
B = 2304*x**12 - 5184*x**10 + 2304*x**9 - 3996*x**8 - 2592*x**7 + 5376*x**6 - 3456*x**5 + 3132*x**4 + 2832*x**3 - 2736*x**2 + 169
resultant = A.resultant(B)

polypy_test.start("Resultants")

p_sqrt2 = x**2 - 2
p_sqrt3 = y**2 - 3

# polypy.trace_enable("polynomial")
# polypy.trace_enable("coefficient::resultant")

add = z - (x + y)
res1 = add.resultant(p_sqrt2)
check_resultant(res1, p_sqrt3, 1*z**4 - 10*z**2 + 1)

sub = z - (x - y)
res1 = sub.resultant(p_sqrt2)
check_resultant(res1, p_sqrt3, 1*z**4 - 10*z**2 + 1)

mul = z - x*y
res1 = mul.resultant(p_sqrt2)
check_resultant(res1, p_sqrt3, 1*z**4 - 12*z**2 + 36)

polypy_test.start("Principal Sub-resultant Coefficients")

p = (y-1)*z*x**2 + y*(z-1)*x + y*z
q =                z*(y-1)*x + y*z
expected = [-y**2*z**2 + y**3*z**2 + y*z**3 - 2*y**2*z**3 + y**3*z**3, -z + y*z]
check_psc(p, q, expected)

p = (y-1)*z*x**2 + y*(z-1)*x + y*z
q = (z-1)*y*x**2 + z*(y-1)*x + y*z
expected = [-y**4*z + y**3*z**2 + 3*y**4*z**2 + y**2*z**3 - 6*y**3*z**3 - y*z**4 + 3*y**2*z**4, 
            -y**2 + 2*y**2*z + z**2 - 2*y*z**2, 
            1]
check_psc(p, q, expected)

p = (y-3)*x**3 + (y-2)*x**2 + (y-1)*x + y
q = (z-3)*x**3 + (z-2)*x**2 + (z-1)*x + z
expected = [16*y**3 - 48*y**2*z + 48*y*z**2 - 16*z**3, 
            0, 
            y - z, 
            1]
check_psc(p, q, expected)
expected = [-16*y**3 + 48*y**2*z - 48*y*z**2 + 16*z**3, 
            0, 
            z - y, 
            1]
check_psc(q, p, expected)

p = x**2 + y**2 - 1
q = p.derivative()
expected = [-4 + 4*y**2, 2]
check_psc(p, q, expected)
check_psc(q, p, expected);

p = (y-1)*x**2 
q = (y-2)*x**2
expected = [0, 0, 1]
check_psc(p, q, expected)

p = x**2
q = x**2
expected = [0, 0, 1]
check_psc(p, q, expected)
 
p = x**2
q = x
expected = [0, 1]
check_psc(p, q, expected)

p = x**2 + y**2 + z**2 - 1
q = p.derivative()
expected = [-4 + 4*y**2 + 4*z**2, 2]
check_psc(p, q, expected)
 
p = (y-3)*x**3 + (y-2)*x**2 + (y-1)*x + y
q =              (z-2)*x**2 + (z-1)*x + z
expected = [-3*y - 5*y**2 + 3*z - y*z + 7*y**2*z + 6*z**2 - y*z**2 - 3*y**2*z**2 + 3*z**3 - 3*y*z**3 + y**2*z**3, 
            -3 + 3*y - 2*z - y*z + z**2, 
            -2 + z]
check_psc(p, q, expected)

p = z*x**3 + (y - 1)*x**2 + z*y
q = y*x**3 + (z - 1)*x**2 + z*y
expected = [y**5*z**2 - y**6*z**2 - 3*y**4*z**3 + 2*y**5*z**3 - y**6*z**3 + 3*y**3*z**4 + 3*y**5*z**4 - y**2*z**5 - 2*y**3*z**5 - 3*y**4*z**5 + y**2*z**6 + y**3*z**6, 
            y**3*z - y**4*z - 2*y**2*z**2 + y**3*z**2 + y*z**3 + y**2*z**3 - y*z**4, 
            y - y**2 - z + z**2, 
            1]
check_psc(p, q, expected)

p = z*x**4 + (y - 1)*x**3 + z*y
q = y*x**4 + (z - 1)*x**3 + z*y
expected = [-y**7*z**3 + y**8*z**3 + 4*y**6*z**4 - 3*y**7*z**4 + y**8*z**4 - 6*y**5*z**5 + 2*y**6*z**5 - 4*y**7*z**5 + 4*y**4*z**6 + 2*y**5*z**6 + 6*y**6*z**6 - y**3*z**7 - 3*y**4*z**7 - 4*y**5*z**7 + y**3*z**8 + y**4*z**8, 
            y**5*z**2 - y**6*z**2 - 3*y**4*z**3 + 2*y**5*z**3 + 3*y**3*z**4 - y**2*z**5 - 2*y**3*z**5 + y**2*z**6, 
            0, 
            y - y**2 - z + z**2, 
            1]
check_psc(p, q, expected)
