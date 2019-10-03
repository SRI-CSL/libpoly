#!/usr/bin/env python

import polypy
import polypy_test

polypy_test.init()

def check_str(o, expected):
    ok = str(o) == expected
    if (not ok):
        print("o = {0}".format(o))
        print("expected = {0}".format(expected))
    polypy_test.check(ok)

x = polypy.Variable("x")
y = polypy.Variable("y")

polypy_test.start("Variable")

order = polypy.VariableOrder([x])
check_str(order, "[x]")

order = polypy.VariableOrder([x, y])
check_str(order, "[x, y]")

order = polypy.variable_order
check_str(order, "[]")

order.push(x)
check_str(order, "[x]")

order.push(y)
check_str(order, "[x, y]")

order.pop()
check_str(order, "[x]")

order.pop()
check_str(order, "[]")

polypy_test.start("Variable Addition")

p = x + 1
check_str(p, "1*x + 1")

p = y + 1
check_str(p, "1*y + 1")

p = 1 + x
check_str(p, "1*x + 1")

p = 1 + y
check_str(p, "1*y + 1")

p = x + y
check_str(p, "1*y + (1*x)")

p = y + x
check_str(p, "1*y + (1*x)")

order.set([x, y])
check_str(p, "1*y + (1*x)")

order.set([y, x])
check_str(p, "1*x + (1*y)")

p = x + x
check_str(p, "2*x")

polypy_test.start("Variable Negation")

p = -x
check_str(p, "-1*x")

p = -y
check_str(p, "-1*y")

polypy_test.start("Variable Subtraction")

p = x - 1
check_str(p, "1*x - 1")

p = 1 - x
check_str(p, "-1*x + 1")

p = y - 1
check_str(p, "1*y - 1")

p = 1 - y
check_str(p, "-1*y + 1")

p = x - y
check_str(p, "1*x + (-1*y)")

p = y - x
check_str(p, "-1*x + (1*y)")

p = x - x
check_str(p, "0")

p = y - y
check_str(p, "0")

polypy_test.start("Variable Multiplication")

p = x*2
check_str(p, "2*x")

p = y*2
check_str(p, "2*y")

p = 2*x
check_str(p, "2*x")

p = 2*y
check_str(p, "2*y")

p = x*y
check_str(p, "(1*y)*x")

p = y*x
check_str(p, "(1*y)*x")

order.set([x, y])
check_str(p, "(1*x)*y")

p = x*x
check_str(p, "1*x**2")

p = y*y
check_str(p, "1*y**2")

polypy_test.start("Variable Power")

p = x**0
check_str(p, "1")

p = x**1
check_str(p, "1*x")

p = x**2
check_str(p, "1*x**2")
