#!/usr/bin/env python

import polypy
import polypy_test

polypy_test.init()

x = polypy.Variable("x");
y = polypy.Variable("y");
z = polypy.Variable("z");

polypy.variable_order.set([z, y, x])

def check_gcd(p, q, expected):
    gcd = p.gcd(q)
    ok = gcd == expected
    if (not ok):
        print("p = {0}".format(p))
        print("q = {0}".format(q))
        print("expected = {0}".format(expected))
        print("gcd = {0}".format(gcd))
    polypy_test.check(ok)

def check_lcm(p, q, expected):
    gcd = p.lcm(q)
    ok = gcd == expected
    if (not ok):
        print("p = {0}".format(p))
        print("q = {0}".format(q))
        print("expected = {0}".format(expected))
        print("gcd = {0}".format(gcd))
    polypy_test.check(ok)

polypy_test.start("GCD")

# polypy.trace_enable("coefficient::gcd")
# polypy.trace_enable("polynomial")
# polypy.trace_enable("coefficient")
# polypy.trace_enable("coefficient::order")
# polypy.trace_enable("coefficient::reduce")

polypy.variable_order.set([z, y, x])

# p = ((1*z**68)*y**71 + (-2*z**48)*y**52 + (26896*z**118 + 26896*z**116 + 6724*z**114)*y**48 + (1*z**28)*y**33 + (-172200*z**59 - 86100*z**57)*y**24 + 275625)
# q = ((-328*z**129 + 164*z**127 + 82*z**125)*y**95 + (656*z**109 - 328*z**107 - 164*z**105)*y**76 + (4410944*z**177 + 6616416*z**175 + 3308208*z**173 + 551368*z**171)*y**72 + (-525*z**68)*y**71 + (-328*z**89 + 164*z**87 + 82*z**85)*y**57 + (1050*z**48)*y**52 + (-42361200*z**118 - 42361200*z**116 - 10590300*z**114)*y**48 + (-525*z**28)*y**33 + (135607500*z**59 + 67803750*z**57)*y**24 - 144703125)
# p.gcd(q)

# Wrong gcd
p = 4*z**4 + 4*z**2 + 1
q = 2*z**2 + 1
expected = 2*z**2 + 1;
check_gcd(p, q, expected)

# Wrong gcd
p = 172200*z**39 + 86100*z**37
q = 26896*z**78 + 26896*z**76 + 6724*z**74
expected = 328*z**39 + 164*z**37
check_gcd(p, q, expected)

p = x*y + 1
q = x*z + 1
expected = 1
check_gcd(p, q, expected)

p = ((-1*z**23)*y**24 + (1*z**3)*y**5)
q = ((1*z**68)*y**71 + (-2*z**48)*y**52 + (26896*z**118 + 26896*z**116 + 6724*z**114)*y**48 + (1*z**28)*y**33 + (-172200*z**59 - 86100*z**57)*y**24 + 275625)
expected = 1
check_gcd(p, q, expected)

p = (1*z**68)*y**71 + (-2*z**48)*y**52 + (26896*z**118 + 26896*z**116 + 6724*z**114)*y**48 + (1*z**28)*y**33 + (-172200*z**59 - 86100*z**57)*y**24 + 275625
q = (328*z**106)*y**71 + (-328*z**86)*y**52
expected = 1
check_gcd(p, q, expected)

# An improved EZ-GCD algorithm for multivariate polynomials, Example 2
p = y**4*z**2*x**4 + (y**2*z**2 + y**3*z**2 + 2*y**2*z**4 + y**3*z**4 + y**4*z**4)*x**3 + (2*y**2*z + y**3*z + 2*z**4 + 3*y*z**4 + 2*y**2*z**4 + y**3*z**4)*x**2 + (2*z + 2*y*z + 2*y*z**3 + y**2*z**3 + y**3*z**3)*x + 2*y
q = y**2*z*x**5 + (2*z**3 + y*z**3 + y**2*z**3)*x**4 + (2 + 2*y**7*z**6)*x**3 + (4*y**5*z**5 + 4*y**5*z**8 + 2*y**6*z**8 + 2*y**7*z**8)*x**2 + (4*y**5*z**5 + 8*y**3*z**7 + 4*y**4*z**7 + 4*y**5*z**7)*x + 8*y**3*z**4
expected =  ((1*z)*y**2)*x**2 + ((1*z**3)*y**2 + (1*z**3)*y + (2*z**3))*x + 2
check_gcd(p, q, expected)

# polypy.trace_enable("polynomial")
# polypy.trace_enable("coefficient")
# polypy.trace_enable("coefficient::reduce")

p = ((1*z**4)*y**5 + (2*z**4))
q = ((1*z**7)*y**8 + (1*z**7)*y**7 + (2*z**7)*y**3 + (2*z**7)*y**2)
expected = y**5*z**4 + 2*z**4
check_gcd(p, q, expected)

p = (1*z)*x**2 + ((1*z**4)*y**3 + (1*z**4)*y**2)*x + ((1*z**6)*y**4)
q = ((1*z**4)*y**5 + (2*z**4))*x**2 + ((1*z**7)*y**8 + (1*z**7)*y**7 + (2*z**7)*y**3 + (2*z**7)*y**2)*x + ((1*z**9)*y**9 + (2*z**9)*y**4)
expected = y**4*z**6 + x*y**3*z**4 + x*y**2*z**4 + x**2*z
check_gcd(p, q, expected)

p = ((1*z**5)*y**3 - 1*y**2 + (1*z)*y)*x**3 + ((1*z**8)*y**6 + (1*z**8)*y**5 + (1*z**4 - 1*z**3)*y**4 + (1*z**4)*y**3 + (2*z**3))*x**2 + ((1*z**6)*y**8 + (1*z**10 + 1*z**6)*y**7 + (-1*z**5)*y**6 + (1*z**6)*y**5 + (2*z**6)*y**3 + (2*z**6)*y**2)*x + ((1*z**8)*y**9 + (2*z**8)*y**4)
q = ((1*z**7 + 1*z)*y**9 + (-2*z**2)*y**8 + (1*z**3)*y**7 + (-2*z**11)*y**6 + (6*z**6)*y**5 + (-2*z**7)*y**4 + (4*z**6))*x**2 + ((1*z**10 + 1*z**4)*y**12 + (1*z**10 - 2*z**5 + 1*z**4)*y**11 + (1*z**6 - 2*z**5)*y**10 + (-2*z**14 + 1*z**6)*y**9 + (-2*z**14 + 6*z**9)*y**8 + (-2*z**10 + 6*z**9)*y**7 + (-2*z**10)*y**6 + (4*z**9)*y**3 + (4*z**9)*y**2)*x + ((1*z**12 + 1*z**6)*y**13 + (-2*z**7)*y**12 + (1*z**8)*y**11 + (-2*z**16)*y**10 + (6*z**11)*y**9 + (-2*z**12)*y**8 + (4*z**11)*y**4)
expected = y**4*z**5 + x*y**3*z**3 + x*y**2*z**3 + x**2
check_gcd(p, q, expected)

# An improved EZ-GCD algorithm for multivariate polynomials, Example 1
p = (x**2 + (y**2*z**3 + y**3*z**3)*x + y**4*z**5)*(x**3 + y**4*z**2*x + 2*y**3*z**4)
q = (x**2 + (y**2*z**3 + y**3*z**3)*x + y**4*z**5)*(x**2 + y**3*z**3*x + y**5*z)
expected = y**4*z**5 + x*y**3*z**3 + x*y**2*z**3 + x**2
check_gcd(p, q, expected)

p = (1*z**6)*y**6 + (-1*z)*y**5 + (1*z**2)*y**4
q = (-1*z**9)*y**12 + (-2*z**9)*y**7
expected = y**4*z
check_gcd(p, q, expected)

p = (2*z**11 - 2*z**5)
q = (1*z**5)
expected = (2*z**11 - 2*z**5)
check_lcm(p, q, expected)
