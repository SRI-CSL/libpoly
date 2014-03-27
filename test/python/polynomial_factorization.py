#!/usr/bin/python

import polypy
import polypy_test
 
polypy_test.init()
 
x = polypy.Variable("x");
y = polypy.Variable("y");
z = polypy.Variable("z"); 

polypy.variable_order.set([z, y, x])

def check_factorization(p, p_factors, expected_size):
    mul = 1
    for (factor, multiplicity) in p_factors:
        mul = mul * (factor**multiplicity)
    ok = (p == mul) and (expected_size == len(p_factors))
    if (not ok):
        print "p =", p
        print "p_factors =", q
        print "expected_size =", expected_size
        print "mul =", mul
    polypy_test.check(ok)

polypy_test.start("Square-Free Factorization")

# polypy.trace_enable("polynomial")
# polypy.trace_enable("factorization")
# polypy.trace_enable("coefficient")

# polynomials and expected factors
test1 = [(x + 1 - x, 0),
         (x + 3 - x, 1),
         (x + 1 - 1, 1), 
         (2*x      , 2),
         (x + 1    , 1),
         (2*x + 2  , 2),
         (2*x + 1  , 1)]

for (p, expected) in test1:
    p_factors = p.factor_square_free()
    check_factorization(p, p_factors, expected)
    
test2 = [(x*y        , 2),
         (x*y + y    , 2),
         (x*y + 1    , 1),
         (2*x*y + 1  , 1),
         (2*x*y + 2  , 2),
         (2*x*y + 2*y, 3)]

for (p, expected) in test2:
    p_factors = p.factor_square_free()
    check_factorization(p, p_factors, expected)

test3 = [( (x + 1)*(x + 2)**2*(x + 3)**3, 3),
         ( (x + 1)*(y + 2)**2*(z + 3)**3, 3),
         ( (z + 1)*(y + 2)**2*(x + 3)**3, 3)]

for (p, expected) in test3:
    p_factors = p.factor_square_free()
    check_factorization(p, p_factors, expected)

