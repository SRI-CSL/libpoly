#!/usr/bin/python

import polypy
import polypy_test
 
polypy_test.init()
 
polypy_test.start("Addition")

x = polypy.Variable("x");
y = polypy.Variable("y");
z = polypy.Variable("z"); 

polypy.variable_order.set([z, y, x])

def check_factorization(p, p_factors):
    mul = 1
    for (factor, multiplicity) in p_factors:
        mul = mul * (factor**multiplicity)
    ok = p == mul;
    if (not ok):
        print "p =", p
        print "p_factors =", q
        print "mul =", mul
    polypy_test.check(ok)

polypy_test.start("Square-Free Factorization")

polypy.trace_enable("polynomial")
polypy.trace_enable("factorization")
polypy.trace_enable("coefficient")

test1 = [x + 1 - x,
         x + 3 - x,
         x + 1 - 1, 
         2*x,
         x + 1,
         2*x + 2,
         2*x + 1]

for p in test1:
    p_factors = p.factor_square_free()
    print p, p_factors
    check_factorization(p, p_factors)
    
test2 = [x*y,
         x*y + y,
         x*y + 1,
         2*x*y + 1,
         2*x*y + 2,
         2*x*y + 2*y]

for p in test2:
    p_factors = p.factor_square_free()
    print p, p_factors
    check_factorization(p, p_factors)

test3 = [(x + 1)*(x + 2)**2*(x + 3)**3,
         (x + 1)*(y + 2)**2*(z + 3)**3,
         (z + 1)*(y + 2)**2*(x + 3)**3]


