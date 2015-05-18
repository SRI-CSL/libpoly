#!/usr/bin/python

import polypy
import polypy_test
 
polypy_test.init()
 
[x, y, z] = [polypy.Variable(name) for name in ['x', 'y', 'z']]
polypy.variable_order.set([z, y, x])

polypy_test.start("Polynomial Root Isolation")

polypy.trace_enable("polynomial")
# polypy.trace_enable("factorization")
# polypy.trace_enable("algebraic_number")
polypy.trace_enable("coefficient")
polypy.trace_enable("coefficient::sgn")
polypy.trace_enable("coefficient::roots")
# polypy.trace_enable("coefficient::arith")

assignment = polypy.Assignment()
assignment.set_value(y, 0)
assignment.set_value(z, 0)

print assignment

p = 10*(x-1)*(x**2 - y - z - 1)
roots = p.roots_isolate(assignment)
print "p_roots =", roots

p = x**2 - y - z - 1
roots = p.roots_isolate(assignment)
print "p_roots =", roots


p = x + y + z
roots = p.roots_isolate(assignment)
print "p_roots =", roots

