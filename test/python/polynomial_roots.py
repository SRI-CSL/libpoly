#!/usr/bin/python

import polypy
import polypy_test
 
polypy_test.init()

def check_roots(p, assignment, roots, expected_roots):
    ok = len(roots) == len(expected_roots)
    if ok: 
        for root, expected_root in zip(roots, expected_roots):
            root_double = root.to_double()
            ok = abs(root_double - expected_root) < 0.000001
            if (not ok):
                print "p =", p
                print "assignment =", assignment
                print "root =", root
                print "expected_root =", expected_root
                break
    else:
        print "p =", p
        print "assignment =", assignment
        print "roots =", roots
        print "expected_roots =", expected_roots
    polypy_test.check(ok)
 
[x, y, z] = [polypy.Variable(name) for name in ['x', 'y', 'z']]
polypy.variable_order.set([z, y, x])

polypy_test.start("Polynomial Root Isolation")

# polypy.trace_enable("polynomial")
# polypy.trace_enable("factorization")
# polypy.trace_enable("algebraic_number")
# polypy.trace_enable("coefficient")
# polypy.trace_enable("coefficient::sgn")
# polypy.trace_enable("coefficient::roots")
# polypy.trace_enable("coefficient::arith")

sqrt2 = polypy.AlgebraicNumber(x**2 - 2, 1)
sqrt3 = polypy.AlgebraicNumber(x**2 - 3, 1)

assignment = polypy.Assignment()
assignment.set_value(y, sqrt2)
assignment.set_value(z, sqrt3)

p = ((1*z)*y**2 + (1*z**2)*y)*x**2 + ((1*z**2)*y**2)*x + 1
roots = p.roots_isolate(assignment)
check_roots(p, assignment, roots, [-0.536830573034912, -0.241708498946619])

# print assignment

p = (x**2 - y**2)*(x**2 - z**2)
roots = p.roots_isolate(assignment);
check_roots(p, assignment, roots, [-1.732050808, -1.414213562, 1.414213562, 1.732050808])

assignment = polypy.Assignment()
assignment.set_value(y, 0)
assignment.set_value(z, 0)

# print assignment

p = 10*(x-1)*(x**2 - y - z - 1)
roots = p.roots_isolate(assignment)
check_roots(p, assignment, roots, [-1, 1])

p = x**2 - y - z - 1
roots = p.roots_isolate(assignment)
check_roots(p, assignment, roots, [-1, 1])

p = x + y + z
roots = p.roots_isolate(assignment)
check_roots(p, assignment, roots, [0])

