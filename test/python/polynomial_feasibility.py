#!/usr/bin/python

import polypy
import polypy_test
 
polypy_test.init()

[x, y, z] = [polypy.Variable(name) for name in ['x', 'y', 'z']]
polypy.variable_order.set([z, y, x])

polypy_test.start("Polynomial Feasibility Intervals")

# polypy.trace_enable("polynomial")
# polypy.trace_enable("factorization")
# polypy.trace_enable("algebraic_number")
# polypy.trace_enable("coefficient")
# polypy.trace_enable("coefficient::sgn")
# polypy.trace_enable("coefficient::roots")
# polypy.trace_enable("value::pick");
# polypy.trace_enable("coefficient::arith")

def get_feasible(p, var, assignment, sgn):
    p_feasible = p.feasible_intervals(assignment, sgn)
    print p_feasible
    for I in p_feasible:
        print "Picking value in", I 
        v = I.pick_value()
        print v

sgns = [polypy.SGN_LT_0, polypy.SGN_LE_0, polypy.SGN_EQ_0, polypy.SGN_NE_0, polypy.SGN_GT_0, polypy.SGN_GE_0]

assignment = polypy.Assignment()
assignment.set_value(y, 1)
assignment.set_value(z, 1)

p = x**2 - y - z
for sgn in sgns:
    get_feasible(p, x, assignment, sgn)


assignment = polypy.Assignment()
assignment.set_value(y, 0)
assignment.set_value(z, 1)

p = (x - y)*(x - z)*(y - z)
for sgn in sgns:
    get_feasible(p, x, assignment, sgn)


assignment = polypy.Assignment()
assignment.set_value(y, 0)
assignment.set_value(z, 1)

p = x**2 + y**2 + z**2
for sgn in sgns:
    get_feasible(p, x, assignment, sgn)
