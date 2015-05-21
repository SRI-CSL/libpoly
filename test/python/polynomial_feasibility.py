#!/usr/bin/python

import polypy
import polypy_test
 
polypy_test.init()

[x, y, z] = [polypy.Variable(name) for name in ['x', 'y', 'z']]
polypy.variable_order.set([z, y, x])

polypy_test.start("Polynomial Feasibility Intervals")

polypy.trace_enable("polynomial")
# polypy.trace_enable("factorization")
# polypy.trace_enable("algebraic_number")
polypy.trace_enable("coefficient")
# polypy.trace_enable("coefficient::sgn")
polypy.trace_enable("coefficient::roots")
# polypy.trace_enable("coefficient::arith")

assignment = polypy.Assignment()
assignment.set_value(y, 0)
assignment.set_value(z, 1)

p = x**2 + y**2 + z**2
p_feasible = p.feasible_intervals(assignment, polypy.SGN_LT_0)
print p_feasible
p_feasible = p.feasible_intervals(assignment, polypy.SGN_LE_0)
print p_feasible
p_feasible = p.feasible_intervals(assignment, polypy.SGN_EQ_0)
print p_feasible
p_feasible = p.feasible_intervals(assignment, polypy.SGN_NE_0)
print p_feasible
p_feasible = p.feasible_intervals(assignment, polypy.SGN_GT_0)
print p_feasible
p_feasible = p.feasible_intervals(assignment, polypy.SGN_GE_0)
print p_feasible
