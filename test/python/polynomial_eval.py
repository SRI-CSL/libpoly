#!/usr/bin/env python

import polypy
import polypy_test

polypy_test.init()

[x, y, z] = [polypy.Variable(name) for name in ['x', 'y', 'z']]
polypy.variable_order.set([z, y, x])

def check_value(p, assignment, value, expected_value):
    value_double = value.to_double()
    ok = abs(value_double - expected_value) < 0.000001
    if (not ok):
        print("p = {0}".format(p))
        print("assignment = {0}".format(assignment))
        print("value = {0}".format(value))
        print("value_double = {0}".format(value_double))
        print("expected_value = {0}".format(expected_value))
    polypy_test.check(ok)

polypy_test.start("Polynomial Evaluation")

# polypy.trace_enable("polynomial")
# polypy.trace_enable("factorization")
# polypy.trace_enable("algebraic_number")
# polypy.trace_enable("coefficient")
# polypy.trace_enable("coefficient::sgn")
# polypy.trace_enable("coefficient::roots")
# polypy.trace_enable("coefficient::arith")

sqrt2 = polypy.AlgebraicNumber(x**2 - 2, 1)
sqrt2_4 = polypy.AlgebraicNumber(x**4 - 2, 1)
sqrt3 = polypy.AlgebraicNumber(x**2 - 3, 1)

assignment = polypy.Assignment()
assignment.set_value(x, sqrt2)
assignment.set_value(y, sqrt2_4)
assignment.set_value(z, sqrt3)

# print assignment

p = x + y + z
p_value = p.evaluate(assignment)
check_value(p, assignment, p_value, 4.335471485)

p = x**2
p_value = p.evaluate(assignment)
check_value(p, assignment, p_value, 2)

p = x**2 + y**4
p_value = p.evaluate(assignment)
check_value(p, assignment, p_value, 4)

p = x**2 + y**4 + z**2
p_value = p.evaluate(assignment)
check_value(p, assignment, p_value, 7)
