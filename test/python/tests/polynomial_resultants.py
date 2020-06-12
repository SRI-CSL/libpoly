#!/usr/bin/env python

import polypy
import polypy_test

x = polypy.Variable("x");
y = polypy.Variable("y");
z = polypy.Variable("z");
w = polypy.Variable("w");

polypy.variable_order.set([w, z, y, x]);

def cmp3(lhs, rhs):
    return ((lhs > rhs) - (lhs < rhs))

def check_psc(p, q, expected):
    psc = p.psc(q)
    ok = cmp3(psc, expected) == 0
    if (not ok):
        print("p = {0}".format(p))
        print("q = {0}".format(q))
        print("psc      = {0}".format(psc))
        print("expected = {0}".format(expected))
    polypy_test.check(ok)

def check_resultant(p, q, expected):
    resultant = p.resultant(q)
    ok = resultant == expected
    if (not ok):
        print("p = {0}".format(p))
        print("q = {0}".format(q))
        print("resultant = {0}".format(resultant))
        print("expected  = {0}".format(expected))

polypy_test.init()

# polypy.trace_enable("coefficient::resultant")
# polypy.trace_enable("polynomial::expensive")

polypy_test.start("Speed")

x5 = polypy.Variable("x5")
x9 = polypy.Variable("x9")
x10 = polypy.Variable("x10")
x11 = polypy.Variable("x11")
x12 = polypy.Variable("x12")

polypy.variable_order.set([x9, x12, x11, x10, x5]);

A = (((1*x12**6 + (4*x9)*x12**5 + (10*x9**2)*x12**4 + (14*x9**3)*x12**3 + (13*x9**4)*x12**2 + (6*x9**5)*x12 + (1*x9**6))*x11**2)*x10**2 + (((2*x9**2)*x12**6 + (4*x9**3)*x12**5 + (6*x9**4)*x12**4 + (2*x9**5)*x12**3)*x11)*x10 + ((1*x9**4)*x12**6))*x5**2 + ((((-1*x9)*x12**3 + (-1*x9**2)*x12**2 + (-1*x9**3)*x12)*x11**4)*x10**4 + ((-2*x12**6 + (-9*x9)*x12**5 + (-17*x9**2)*x12**4 + (-21*x9**3)*x12**3 + (-12*x9**4)*x12**2 + (-4*x9**5)*x12)*x11**3)*x10**3 + (((-2*x9**2)*x12**6 + (-6*x9**3 - 2*x9)*x12**5 + (-4*x9**4 - 8*x9**2)*x12**4 + (-3*x9**5 - 16*x9**3)*x12**3 + (-1*x9**6 - 18*x9**4)*x12**2 + (-1*x9**7 - 10*x9**5)*x12 + (-2*x9**6))*x11**2)*x10**2 + (((-1*x9**5 - 2*x9**3)*x12**5 + (-1*x9**6 - 4*x9**4)*x12**4 + (-1*x9**7 - 2*x9**5)*x12**3)*x11)*x10)*x5 + ((((1*x9)*x12**3 + (1*x9**2)*x12**2)*x11**5)*x10**5 + ((1*x12**6 + (5*x9)*x12**5 + (7*x9**2)*x12**4 + (6*x9**3)*x12**3 + (3*x9**4)*x12**2 + (1*x9**3)*x12)*x11**4)*x10**4 + (((2*x9**3 + 2*x9)*x12**5 + (2*x9**4 + 7*x9**2)*x12**4 + (1*x9**5 + 13*x9**3)*x12**3 + (1*x9**6 + 8*x9**4)*x12**2 + (4*x9**5)*x12)*x11**3)*x10**3 + (((-1*x9**2)*x12**6 + (1*x9**5)*x12**5 + (1*x9**6 - 2*x9**4 + 1*x9**2)*x12**4 + (2*x9**5 + 4*x9**3)*x12**3 + (6*x9**4)*x12**2 + (1*x9**7 + 4*x9**5)*x12 + (1*x9**6))*x11**2)*x10**2 + (((-2*x9**4)*x12**6 + (-1*x9**6)*x12**4 + (1*x9**7)*x12**3)*x11)*x10 + ((-1*x9**6)*x12**6))
B = (((2*x12**6 + (8*x9)*x12**5 + (20*x9**2)*x12**4 + (28*x9**3)*x12**3 + (26*x9**4)*x12**2 + (12*x9**5)*x12 + (2*x9**6))*x11**2)*x10**2 + (((4*x9**2)*x12**6 + (8*x9**3)*x12**5 + (12*x9**4)*x12**4 + (4*x9**5)*x12**3)*x11)*x10 + ((2*x9**4)*x12**6))*x5 + ((((-1*x9)*x12**3 + (-1*x9**2)*x12**2 + (-1*x9**3)*x12)*x11**4)*x10**4 + ((-2*x12**6 + (-9*x9)*x12**5 + (-17*x9**2)*x12**4 + (-21*x9**3)*x12**3 + (-12*x9**4)*x12**2 + (-4*x9**5)*x12)*x11**3)*x10**3 + (((-2*x9**2)*x12**6 + (-6*x9**3 - 2*x9)*x12**5 + (-4*x9**4 - 8*x9**2)*x12**4 + (-3*x9**5 - 16*x9**3)*x12**3 + (-1*x9**6 - 18*x9**4)*x12**2 + (-1*x9**7 - 10*x9**5)*x12 + (-2*x9**6))*x11**2)*x10**2 + (((-1*x9**5 - 2*x9**3)*x12**5 + (-1*x9**6 - 4*x9**4)*x12**4 + (-1*x9**7 - 2*x9**5)*x12**3)*x11)*x10)
A.psc(B)

polypy.variable_order.set([w, z, y, x]);

A = w - 2*x**3*y**2*z - x**2*y*z**3 - x*y**3*z**2 - z - 4
B = 2304*x**12 - 5184*x**10 + 2304*x**9 - 3996*x**8 - 2592*x**7 + 5376*x**6 - 3456*x**5 + 3132*x**4 + 2832*x**3 - 2736*x**2 + 169
resultant = A.resultant(B)

polypy_test.start("Bugs")

polypy.variable_order.set([w, z, y, x]);

p = 7*x**6 + 14*x**5 + (7*y - 1)
q = 42*x**5 + 70*x**4
expected = [92254156521408*y**5-461361174686720*y**4+244807578081920*y**3-51113953954560*y**2+4803956911040*y-170197631744,
            174327582240*y**3-74711820960*y**2+10673117280*y-508243680,
            0,
            0,
            -6860,
            42]
check_psc(p, q, expected)

polypy_test.start("Resultants")

p_sqrt2 = x**2 - 2
p_sqrt3 = y**2 - 3

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
