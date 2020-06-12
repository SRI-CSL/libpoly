#!/usr/bin/env python

import sys
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
                print("p = {0}".format(p))
                print("assignment = {0}".format(assignment))
                print("root = {0}".format(root))
                print("expected_root = {0}".format(expected_root))
                break
    else:
        print("p = {0}".format(p))
        print("assignment = {0}".format(assignment))
        print("roots = {0}".format(roots))
        print("expected_roots = {0}".format(expected_roots))
    polypy_test.check(ok)

# polypy.trace_enable("polynomial")
# polypy.trace_enable("polynomial::expensive")
# polypy.trace_enable("factorization")
# polypy.trace_enable("algebraic_number")
# polypy.trace_enable("coefficient")
# polypy.trace_enable("coefficient::sgn")
# polypy.trace_enable("coefficient::roots")
# polypy.trace_enable("roots")


[x, y, z, w] = [polypy.Variable(name) for name in ['x', 'y', 'z', 'w']]
[x0, x1, x2, x3, x4, x5, x6, x7] = [polypy.Variable(name) for name in ['x0', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7']]

# polypy_test.start("Polynomial Root Isolation: Speed")
#
# polypy.variable_order.set([y, x])
#
# y_value_poly = 8192*x**6 + (-14336*x**5) + 75376*x**4 + (-32736*x**3) + 109496*x**2 + 133752*x + (-32441)
# y_value = polypy.AlgebraicNumber(y_value_poly, 1)
#
# assignment = polypy.Assignment()
# assignment.set_value(y, y_value)
#
# p = (1*y**2 + 3)*x**6 + (4*y**2 + 24*y + 20)*x**5 + (36*y**2 + 112*y + 44)*x**4 + (8*y**3 + 114*y**2 + 144*y + 30)*x**3 + (1*y**4 + 24*y**3 + 142*y**2 + 40*y - 15)*x**2 + (2*y**4 + 40*y**3 + 64*y**2 - 32*y - 26)*x + (5*y**4 + 16*y**3 + 7*y**2 - 16*y - 8)
# roots = p.roots_isolate(assignment)

polypy_test.start("Polynomial Root Isolation: Bugs");

polypy.variable_order.set([x0, x1, x2, x3, x4, x5, x6, x7])
p = ((16*x0**2 - 32*x0 + 16)*x5**2 + ((32*x0 - 32)*x1)*x5 + (16*x1**2))*x7**2 + ((8*x0**2 - 16*x0 + 8)*x5**2 + ((16*x0**2 - 32*x0 + 16)*x2**2 + (-8*x1**2)))*x7 + ((1*x0**2 - 2*x0 + 1)*x5**2 + ((-2*x0 + 2)*x1)*x5 + (1*x1**2))

p_lc = ((16*x0**2 - 32*x0 + 16)*x5**2 + ((32*x0 - 32)*x1)*x5 + (16*x1**2))

assignment = polypy.Assignment()
v = polypy.AlgebraicNumber(32*x**2 + (-64*x) + 15, 0)
assignment.set_value(x0, v)
assignment.set_value(x1, 28);
assignment.set_value(x2, -1, 4)
assignment.set_value(x3, 3)
assignment.set_value(x4, 277771, 4096)
assignment.set_value(x5, 786743, 20480)
assignment.set_value(x6, -6437)

roots = p.roots_isolate(assignment);
check_roots(p, assignment, roots, [-166318.899794716301953293454440, -8961.48539760581031126063921576])

polypy.variable_order.set([x, w, y, z])

sqrt180 = polypy.AlgebraicNumber(x**2 - 180, 0)
assignment = polypy.Assignment()
assignment.set_value(x, 9)
assignment.set_value(w, sqrt180)
assignment.set_value(y, 0)

p = ((1*x)*y)*z + (-1*w**2 + (20*x))

roots = p.roots_isolate(assignment)
check_roots(p, assignment, roots, []);


polypy_test.start("Polynomial Root Isolation: Basic")

polypy.variable_order.set([z, y, x])


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
