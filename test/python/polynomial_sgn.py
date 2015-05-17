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

polypy_test.start("Sign Determination")

polypy.trace_enable("polynomial")
# polypy.trace_enable("factorization")
# polypy.trace_enable("algebraic_number")
polypy.trace_enable("coefficient")
polypy.trace_enable("coefficient::sgn")
# polypy.trace_enable("coefficient::arith")

sqrt2 = polypy.AlgebraicNumber(x**2 - 2, 1)

assignment = polypy.Assignment()
assignment.set_value(x, sqrt2)
assignment.set_value(y, -sqrt2)
assignment.set_value(z, 0)

print assignment

p = x + y + z
check_sgn(p, assignment, 0)

p = x + y + z + 1
check_sgn(p, assignment, 1)

p = x + y + z - 1
check_sgn(p, assignment, -1)

assignment = polypy.Assignment()
assignment.set_value(x, sqrt2)
assignment.set_value(y, 1)
assignment.set_value(z, -1)

# print assignment

p = x + y + z
check_sgn(p, assignment, 1)

p = x + y + z + 1
check_sgn(p, assignment, 1)

p = x + y + z - 2
check_sgn(p, assignment, -1)

assignment = polypy.Assignment()
assignment.set_value(x, 0)
assignment.set_value(y, 1)
assignment.set_value(z, -1)

# print assignment

p = x + y + z
check_sgn(p, assignment, 0)

p = x + y + z + 1
check_sgn(p, assignment, 1)

p = x + y + z - 2
check_sgn(p, assignment, -1)

assignment.set_value(x, -1)
assignment.set_value(y, -1)
assignment.set_value(z, -1)

# print assignment

p = (y+z)*x**3
check_sgn(p, assignment, 2)

p = (y-z)*x**2
check_sgn(p, assignment, 0)

p = y*x
check_sgn(p, assignment, 1)

p = (y+z)*x**3 + (y-z)*x**2 + y*x + z
check_sgn(p, assignment, 2)
