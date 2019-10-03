#!/usr/bin/env python

import polypy
import polypy_test

import sys

polypy_test.init()

x = polypy.Variable("x");
y = polypy.Variable("y");

polypy.variable_order.set([y, x])

polypy_test.start("model-based GCD")

# polypy.trace_enable("coefficient::gcd")
# polypy.trace_enable("polynomial")
# polypy.trace_enable("coefficient")
# polypy.trace_enable("coefficient::order")
# polypy.trace_enable("coefficient::reduce")
polypy.trace_enable("coefficient::mgcd")

x0 = polypy.Variable("x0");
x1 = polypy.Variable("x1");
x2 = polypy.Variable("x2");
x3 = polypy.Variable("x3");
x4 = polypy.Variable("x4");

polypy.variable_order.set([x0, x1, x2, x3, x4])

m = polypy.Assignment()
m.set_value(x0, 2)
m.set_value(x1, 3)
m.set_value(x2, 3)
m.set_value(x3, 9)


A = (1*x3)*x4**4 + ((-2*x0)*x3 + (2*x1))*x4**3 + ((1*x0**2)*x3 + (2*x2 + ((-2*x0)*x1 - 40)))*x4**2 + (80*x0)*x4 + (-40*x0**2)
B = A.derivative()

PSC = A.psc(B)
print("PSC = {0}".format(PSC))
MGCD = A.mgcd(B, m)
print("MGCD = {0}".format(MGCD))

sys.exit()



m = polypy.Assignment()
m.set_value(y, -1)
A = 1*x**11 - 110*x**9 + 7920*x**7 - 332640*x**5 + 6652800*x**3 + (110*y**8 - 7920*y**6 + 332640*y**4 - 6652800*y**2)*x
B = 72*x**6 - 3024*x**4 + 60480*x**2 + (1*y**8 - 72*y**6 + 3024*y**4 - 60480*y**2)
print(A.mgcd(B, m))

sys.exit()

m = polypy.Assignment()
m.set_value(y, -1)

A = 1*x
B = 72*x**6 - 3024*x**4 + 60480*x**2 + (1*y**8 - 72*y**6 + 3024*y**4 - 60480*y**2)
C = 1*x**11 - 110*x**9 + 7920*x**7 - 332640*x**5 + 6652800*x**3 + (110*y**8 - 7920*y**6 + 332640*y**4 - 6652800*y**2)*x

m = polypy.Assignment()
m.set_value(y, 0)

A = 72*x**6 - 3024*x**4 + 60480*x**2 + (y**8 - 72*y**6 + 3024*y**4 - 60480*y**2)
B = A.derivative()

PSC = A.psc(B)
MGCD = A.mgcd(B, m)

print("A = {0}".format(A))
print("B = {0}".format(B))
print("PSC = {0}".format(PSC))
print("MGCD = {0}".format(MGCD))

sys.exit(0)

x5 = polypy.Variable("x5")
x9 = polypy.Variable("x9")
x10 = polypy.Variable("x10")
x11 = polypy.Variable("x11")
x12 = polypy.Variable("x12")

polypy.variable_order.set([x9, x12, x11, x10, x5]);

m = polypy.Assignment()
m.set_value(x9, 0)
m.set_value(x12, 0)
m.set_value(x11, 0)
m.set_value(x10, 0)

A = (((1*x12**6 + (4*x9)*x12**5 + (10*x9**2)*x12**4 + (14*x9**3)*x12**3 + (13*x9**4)*x12**2 + (6*x9**5)*x12 + (1*x9**6))*x11**2)*x10**2 + (((2*x9**2)*x12**6 + (4*x9**3)*x12**5 + (6*x9**4)*x12**4 + (2*x9**5)*x12**3)*x11)*x10 + ((1*x9**4)*x12**6))*x5**2 + ((((-1*x9)*x12**3 + (-1*x9**2)*x12**2 + (-1*x9**3)*x12)*x11**4)*x10**4 + ((-2*x12**6 + (-9*x9)*x12**5 + (-17*x9**2)*x12**4 + (-21*x9**3)*x12**3 + (-12*x9**4)*x12**2 + (-4*x9**5)*x12)*x11**3)*x10**3 + (((-2*x9**2)*x12**6 + (-6*x9**3 - 2*x9)*x12**5 + (-4*x9**4 - 8*x9**2)*x12**4 + (-3*x9**5 - 16*x9**3)*x12**3 + (-1*x9**6 - 18*x9**4)*x12**2 + (-1*x9**7 - 10*x9**5)*x12 + (-2*x9**6))*x11**2)*x10**2 + (((-1*x9**5 - 2*x9**3)*x12**5 + (-1*x9**6 - 4*x9**4)*x12**4 + (-1*x9**7 - 2*x9**5)*x12**3)*x11)*x10)*x5 + ((((1*x9)*x12**3 + (1*x9**2)*x12**2)*x11**5)*x10**5 + ((1*x12**6 + (5*x9)*x12**5 + (7*x9**2)*x12**4 + (6*x9**3)*x12**3 + (3*x9**4)*x12**2 + (1*x9**3)*x12)*x11**4)*x10**4 + (((2*x9**3 + 2*x9)*x12**5 + (2*x9**4 + 7*x9**2)*x12**4 + (1*x9**5 + 13*x9**3)*x12**3 + (1*x9**6 + 8*x9**4)*x12**2 + (4*x9**5)*x12)*x11**3)*x10**3 + (((-1*x9**2)*x12**6 + (1*x9**5)*x12**5 + (1*x9**6 - 2*x9**4 + 1*x9**2)*x12**4 + (2*x9**5 + 4*x9**3)*x12**3 + (6*x9**4)*x12**2 + (1*x9**7 + 4*x9**5)*x12 + (1*x9**6))*x11**2)*x10**2 + (((-2*x9**4)*x12**6 + (-1*x9**6)*x12**4 + (1*x9**7)*x12**3)*x11)*x10 + ((-1*x9**6)*x12**6))
B = (((2*x12**6 + (8*x9)*x12**5 + (20*x9**2)*x12**4 + (28*x9**3)*x12**3 + (26*x9**4)*x12**2 + (12*x9**5)*x12 + (2*x9**6))*x11**2)*x10**2 + (((4*x9**2)*x12**6 + (8*x9**3)*x12**5 + (12*x9**4)*x12**4 + (4*x9**5)*x12**3)*x11)*x10 + ((2*x9**4)*x12**6))*x5 + ((((-1*x9)*x12**3 + (-1*x9**2)*x12**2 + (-1*x9**3)*x12)*x11**4)*x10**4 + ((-2*x12**6 + (-9*x9)*x12**5 + (-17*x9**2)*x12**4 + (-21*x9**3)*x12**3 + (-12*x9**4)*x12**2 + (-4*x9**5)*x12)*x11**3)*x10**3 + (((-2*x9**2)*x12**6 + (-6*x9**3 - 2*x9)*x12**5 + (-4*x9**4 - 8*x9**2)*x12**4 + (-3*x9**5 - 16*x9**3)*x12**3 + (-1*x9**6 - 18*x9**4)*x12**2 + (-1*x9**7 - 10*x9**5)*x12 + (-2*x9**6))*x11**2)*x10**2 + (((-1*x9**5 - 2*x9**3)*x12**5 + (-1*x9**6 - 4*x9**4)*x12**4 + (-1*x9**7 - 2*x9**5)*x12**3)*x11)*x10)

print(A.mgcd(B, m))

m = polypy.Assignment()
m.set_value(y, 0)

polypy.variable_order.set([y, x])

p = (x - y)*(x + y) - 1
q = x - 1

print(p.gcd(q))
print(p.mgcd(q, m))
