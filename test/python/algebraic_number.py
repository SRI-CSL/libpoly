#!/usr/bin/python

import polypy
import polypy_test

from polypy import x

polypy_test.init()
polypy_test.start("Construction")

a1 = polypy.AlgebraicNumber((x**2 - 2)*(x - 10), 1)
a2 = polypy.AlgebraicNumber((x**2 - 2), 1)
print a1 < a2
print a1
print a2

p = x**2
zero = polypy.AlgebraicNumber(p, 0);
print zero

p = x**2 - 1
one_neg = polypy.AlgebraicNumber(p, 0);
one_pos = polypy.AlgebraicNumber(p, 1);
print one_neg, one_pos

p = x**2 - 2
sqrt2_neg = polypy.AlgebraicNumber(p, 0);
sqrt2_pos = polypy.AlgebraicNumber(p, 1);
print sqrt2_neg, sqrt2_pos

p = x**2 - 3
sqrt3_neg = polypy.AlgebraicNumber(p, 0);
sqrt3_pos = polypy.AlgebraicNumber(p, 1);
print sqrt3_neg, sqrt3_pos 

a = [zero, one_neg, sqrt2_pos, one_pos, sqrt2_neg]
print a;
a.sort()
print a;

a = []
for n in range(1, 10):
    p = x**2 - n
    a.append(polypy.AlgebraicNumber(p, 0))
    a.append(polypy.AlgebraicNumber(p, 1))
a.sort()
print a
a_double = [b.to_double() for b in a]
print a_double

p1 = (x**2 - 2)*x        # -sqrt(2), 0, sqrt(2)
p2 = (x**2 - 2)*(x-1)    # -sqrt(2), 1, sqrt(2)
p3 = (x**2 - 2)*(x+2)    # -sqrt(2), -1, sqrt(2)

a1 = polypy.AlgebraicNumber(p1, 2)
a2 = polypy.AlgebraicNumber(p2, 2)
a3 = polypy.AlgebraicNumber(p3, 2)

a = [a1, a2, a3]
print a
a.sort()
print a

a_double = [x.to_double() for x in a]
print a_double

polypy_test.start("Arithmetic")

# polypy.trace_enable("coefficient")
# polypy.trace_enable("coefficient::resultant")
# polypy.trace_enable("algebraic_number")

x = polypy.x

a = polypy.AlgebraicNumber(x**2 - 2, 1)
b = polypy.AlgebraicNumber(x**2 - 3, 1)

mul = a * a
print mul, mul.to_double()

mul = a * b
print mul

a = polypy.AlgebraicNumber(x**2 - 2, 0)
b = polypy.AlgebraicNumber(x**2 - 2, 1)
add = a + b
print add, add.to_double()

a = polypy.AlgebraicNumber(x**2 - 2, 1)
b = polypy.AlgebraicNumber(x**2 - 3, 1)
add = a + b
print add, add.to_double()



