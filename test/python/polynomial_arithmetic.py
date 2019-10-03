#!/usr/bin/env python

import polypy
import polypy_test

def check_binary(op, op_name, p, q, expected):
    result = op(p, q)
    ok = result == expected
    if (not ok):
        print("p = {0}".format(p))
        print("q = {0}".format(q))
        print("{0} = {1}".format(op_name, result))
        print("expected = {0}".format(expected))
    polypy_test.check(ok)

def check_unary(op, op_name, p, expected):
    result = op(p)
    ok = result == expected
    if (not ok):
        print("p = {0}".format(p))
        print("{0} = {1}".format(op_name, result))
        print("expected = {0}".format(expected))
    polypy_test.check(ok)

polypy_test.init()

polypy_test.start("Addition")

x = polypy.Variable("x");
y = polypy.Variable("y");
z = polypy.Variable("z");

polypy.variable_order.set([z, y, x])

# polypy.trace_enable("polynomial")
# polypy.trace_enable("coefficient")
# polypy.trace_enable("coefficient::reduce")

def poly_plus(p, q):
    return p + q

p = x**2 + 2*x + 1
q = y**2 + 2*y + 1
expected = x**2 + y**2 + 2*x + 2*y + 2
check_binary(poly_plus, "plus", p, q, expected)

p = x**2
q = 1
expected = x**2 + 1
check_binary(poly_plus, "plus", p, q, expected)
check_binary(poly_plus, "plus", q, p, expected)

p = x**2
q = x
expected = x**2 + x
check_binary(poly_plus, "plus", p, q, expected)
check_binary(poly_plus, "plus", q, p, expected)

polypy_test.start("Negation")

def poly_neg(p):
    return -p

p = x**2 + 2*x + 1
expected = -x**2 -2*x - 1
check_unary(poly_neg, "neg", p, expected)

q = y**2 + 2*y + 1
expected = -y**2 -2*y - 1
check_unary(poly_neg, "neg", q, expected)

polypy_test.start("Subtraction")

def poly_sub(p, q):
    return p - q

p = x**2 + 2*x + 1
q = y**2 + 2*y + 1
expected = x**2 + 2*x - y**2 - 2*y
check_binary(poly_sub, "sub", p, q, expected)

polypy_test.start("Multiplication")

def poly_mul(p, q):
    return p*q

p = 3*y + (-3*z)
q = 1*y - 3
expected = 3*y**2 + (-3*z - 9)*y + 9*z
check_binary(poly_mul, "mul", p, q, expected)

p = (1*y + (-1*z))*x**2 + (2*y + (-2*z))*x + (3*y + (-3*z))
q = (1*y - 3)*x
expected = (y**2 + (-z - 3)*y + 3*z)*x**3 + (2*y**2 + (-2*z - 6)*y + 6*z)*x**2 + (3*y**2 + (-3*z - 9)*y + 9*z)*x
check_binary(poly_mul, "mul", p, q, expected)

p = x + 1
q = y + 1
expected = x*y + x + y + 1
check_binary(poly_mul, "mul", p, q, expected)

expected = 3*x + 3
check_binary(poly_mul, "mul", 3, p, expected)
check_binary(poly_mul, "mul", p, 3, expected)

expected = x*z + z
check_binary(poly_mul, "mul", p, z, expected)

expected = x*y*z + x*z + y*z + z;
check_binary(poly_mul, "mul", (p*q), z, expected)

p = (1*y - 3)*x**3 + (1*y - 2)*x**2 + (1*y - 1)*x + (1*y)
q = 1*y + (-1*z)
expected = (y**2 + (-z - 3)*y + 3*z)*x**3 + (y**2 + (-z - 2)*y + 2*z)*x**2 + (y**2 + (-z - 1)*y + z)*x + y**2 - z*y
check_binary(poly_mul, "mul", p, q, expected)

polypy_test.start("Power")

def poly_pow(p, k):
    return p**k

p = (x + y)
expected = x**2 + 2*x*y + y**2
check_binary(poly_pow, "pow", p, 2, expected)

expected = [1,
            x + 1,
            x**2 + 2*x + 1,
            x**3 + 3*x**2 + 3*x + 1,
            x**4 + 4*x**3 + 6*x**2 + 4*x + 1,
            x**5 + 5*x**4 + 10*x**3 + 10*x**2 + 5*x + 1,
            x**6 + 6*x**5 + 15*x**4 + 20*x**3 + 15*x**2 + 6*x + 1,
            x**7 + 7*x**6 + 21*x**5 + 35*x**4 + 35*x**3 + 21*x**2 + 7*x + 1,
            x**8 + 8*x**7 + 28*x**6 + 56*x**5 + 70*x**4 + 56*x**3 + 28*x**2 + 8*x + 1,
            x**9 + 9*x**8 + 36*x**7 + 84*x**6 + 126*x**5 + 126*x**4 + 84*x**3 + 36*x**2 + 9*x + 1,
            x**10 + 10*x**9 + 45*x**8 + 120*x**7 + 210*x**6 + 252*x**5 + 210*x**4 + 120*x**3 + 45*x**2 + 10*x + 1,
            x**11 + 11*x**10 + 55*x**9 + 165*x**8 + 330*x**7 + 462*x**6 + 462*x**5 + 330*x**4 + 165*x**3 + 55*x**2 + 11*x + 1,
            x**12 + 12*x**11 + 66*x**10 + 220*x**9 + 495*x**8 + 792*x**7 + 924*x**6 + 792*x**5 + 495*x**4 + 220*x**3 + 66*x**2 + 12*x + 1,
            x**13 + 13*x**12 + 78*x**11 + 286*x**10 + 715*x**9 + 1287*x**8 + 1716*x**7 + 1716*x**6 + 1287*x**5 + 715*x**4 + 286*x**3 + 78*x**2 + 13*x + 1,
            x**14 + 14*x**13 + 91*x**12 + 364*x**11 + 1001*x**10 + 2002*x**9 + 3003*x**8 + 3432*x**7 + 3003*x**6 + 2002*x**5 + 1001*x**4 + 364*x**3 + 91*x**2 + 14*x + 1,
            x**15 + 15*x**14 + 105*x**13 + 455*x**12 + 1365*x**11 + 3003*x**10 + 5005*x**9 + 6435*x**8 + 6435*x**7 + 5005*x**6 + 3003*x**5 + 1365*x**4 + 455*x**3 + 105*x**2 + 15*x + 1]

for k in range(0, 16):
    check_binary(poly_pow, "pow", (x + 1), k, expected[k])

polypy_test.start("Derivative");

def poly_d(p):
    return p.derivative()

polypy.variable_order.set([z, y, x]);

p = (x + y + z + 1)**3
expected = 3*x**2 + 6*x*y + 3*y**2 + 6*x*z + 6*y*z + 3*z**2 + 6*x + 6*y + 6*z + 3
check_unary(poly_d, "derivative", p, expected)

polypy_test.start("Division")

polypy.variable_order.set([z, y, x]);

def poly_div(p, q):
    return p / q

def poly_rem(p, q):
    return p % q

p = 2*x + 2
q = 2
expected = x + 1
check_binary(poly_div, "div", p, q, expected)

p = x**2 - 1
q = x - 1
expected = x + 1
check_binary(poly_div, "div", p, q, expected)

p = x**2 - 1
q = x - 1
expected = 0
check_binary(poly_rem, "rem", p, q, expected)

p = (1*z**11 - 1*z**5)
q = (z**5)
expected = z**6 - 1
check_binary(poly_div, "div", p, q, expected)

p = x*y*z
q = (x*y)
expected = z
check_binary(poly_div, "div", p, q, expected)

p = 6*y*(x + 1)
q = (2*y)
expected = 3*(x + 1)
check_binary(poly_div, "div", p, q, expected)

p = 6*y*(x + 1)
q = 3*(x+1)
expected = 2*y
check_binary(poly_div, "div", p, q, expected)

q = x + 1
expected = 6*y
check_binary(poly_div, "div", p, q, expected)

p = 2*x**2 + 4*x + 6
q = 2
expected = x**2 + 2*x + 3
check_binary(poly_div, "div", p, q, expected)

p = 2*x**2 + 4*x + 6
q = (x + 2) - x
expected = x**2 + 2*x + 3
check_binary(poly_div, "div", p, q, expected)

p = (x + 6) - x
q = (x + 3) - x
expected = 2
check_binary(poly_div, "div", p, q, expected)
