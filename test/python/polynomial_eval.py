#!/usr/bin/python

import polypy
import polypy_test
 
polypy_test.init()
 
[x, y, z] = [polypy.Variable(name) for name in ['x', 'y', 'z']]
polypy.variable_order.set([z, y, x])

def check_sgn(p, assignment, expected_sgn):
    sgn = p.sgn(assignment)
    ok = (sgn > 0 and expected_sgn > 0) or (sgn < 0 and expected_sgn < 0) or (sgn == 0 and expected_sgn == 0)
    if (not ok):
        print "p =", p
        print "assignment =", assignment
        print "sgn =", sgn
        print "expected_sgn =", expected_sgn
    polypy_test.check(ok)

polypy_test.start("Polynomial Evaluation")

polypy.trace_enable("polynomial")
# polypy.trace_enable("factorization")
# polypy.trace_enable("algebraic_number")
polypy.trace_enable("coefficient")
polypy.trace_enable("coefficient::sgn")
# polypy.trace_enable("coefficient::arith")

sqrt2 = polypy.AlgebraicNumber(x**2 - 2, 1)
sqrt2_4 = polypy.AlgebraicNumber(x**4 - 2, 1)
sqrt3 = polypy.AlgebraicNumber(x**2 - 3, 1)

assignment = polypy.Assignment()
assignment.set_value(x, sqrt2_4)
assignment.set_value(y, sqrt2)
assignment.set_value(z, sqrt3)

print assignment

p = x**2
p_value = p.evaluate(assignment)
print p_value

p = x**2 + y**4
p_value = p.evaluate(assignment)
print p_value

p = x**2 + y**4 + z**2
p_value = p.evaluate(assignment)
print p_value

