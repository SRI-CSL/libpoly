#!/usr/bin/python

import polypy
import polypy_test
 
polypy_test.init()
 
x = polypy.Variable("x");
y = polypy.Variable("y");

polypy.variable_order.set([y, x])

polypy_test.start("model-based GCD")

# polypy.trace_enable("coefficient::gcd")
# polypy.trace_enable("polynomial")
# polypy.trace_enable("coefficient")
# polypy.trace_enable("coefficient::order")
# polypy.trace_enable("coefficient::reduce")

# Wrong gcd
m = polypy.Assignment()
m.set_value(y, 0)

p = (x - y)*(x + y) - 1
q = x - 1

print p.gcd(q)
print p.mgcd(q, m)
