#!/usr/bin/python

import polypy
import polypy_test

from polypy import x

def check_comparison(a1, a2, cmp, result, expected):
    if (result != expected):
        print "a1 =", a1
        print "a2 =", a2
        print "cmp =", cmp
        print "result =", result
        print "expected =", expected
        polypy_test.check(False)
    else:
        polypy_test.check(True)
        
polypy_test.init()
polypy_test.start("Construction and order")

a1 = polypy.AlgebraicNumber((x**2 - 2)*(x - 10), 1)
a2 = polypy.AlgebraicNumber((x**2 - 2), 1)

result = (a1 < a2)
check_comparison(a1, a2, "<", result, False)

result = (a1 == a2)
check_comparison(a1, a2, "==", result, True)

p = x**2
zero = polypy.AlgebraicNumber(p, 0);
result = zero.to_double() == 0.0
check_comparison(zero, 0, ".to_double ==", result, True)

p = x**2 - 1
one_neg = polypy.AlgebraicNumber(p, 0);
one_pos = polypy.AlgebraicNumber(p, 1);
result = (one_neg < one_pos)
check_comparison(one_neg, one_pos, "<", result, True)

p = x**2 - 2
sqrt2_neg = polypy.AlgebraicNumber(p, 0);
sqrt2_pos = polypy.AlgebraicNumber(p, 1);
result = (sqrt2_neg < sqrt2_pos)
check_comparison(sqrt2_neg, sqrt2_pos, "<", result, True)

p = x**2 - 3
sqrt3_neg = polypy.AlgebraicNumber(p, 0);
sqrt3_pos = polypy.AlgebraicNumber(p, 1);
result = (sqrt3_neg < sqrt3_pos)
check_comparison(sqrt3_neg, sqrt3_pos, "<", result, True)

a = [sqrt3_pos, zero, sqrt3_neg, one_neg, sqrt2_pos, one_pos, sqrt2_neg]
a_sorted = [sqrt3_neg, sqrt2_neg, one_neg, zero, one_pos, sqrt2_pos, sqrt3_pos]
a.sort()

result = (a == a_sorted)
check_comparison(a, a_sorted, "==", result, True)

a = []
a_sorted = []
for n in range(1, 10):
    p = x**2 - n
    neg = polypy.AlgebraicNumber(p, 0)
    pos = polypy.AlgebraicNumber(p, 1)
    a.append(pos)
    a.append(neg)
    a_sorted.append(pos)
    a_sorted.insert(0, neg)
a.sort()
check_comparison(a, a_sorted, "==", result, True)

p1 = (x**2 - 2)*x        # -sqrt(2), 0, sqrt(2)
p2 = (x**2 - 2)*(x-1)    # -sqrt(2), 1, sqrt(2)
p3 = (x**2 - 2)*(x+2)    # -sqrt(2), -1, sqrt(2)

a1 = polypy.AlgebraicNumber(p1, 2)
a2 = polypy.AlgebraicNumber(p2, 2)
a3 = polypy.AlgebraicNumber(p3, 2)

check_comparison(a1, a2, "==", (a1 == a2), True)
check_comparison(a1, a3, "==", (a1 == a3), True)
check_comparison(a2, a3, "==", (a2 == a3), True)

polypy_test.start("Arithmetic")

# polypy.trace_enable("coefficient")
# polypy.trace_enable("coefficient::resultant")
# polypy.trace_enable("algebraic_number")

x = polypy.x

a = polypy.AlgebraicNumber(x**2 - 2, 1)
a2 = polypy.AlgebraicNumber(x**2 - 4, 1)

mul = a * a
check_comparison(mul, a2, "==", (mul == a2), True)

b = polypy.AlgebraicNumber(x**2 - 3, 1)
ab = polypy.AlgebraicNumber(x**2 - 6, 1)
mul = a * b
check_comparison(mul, ab, "==", (mul == ab), True)

a = polypy.AlgebraicNumber(x**2 - 2, 0)
b = polypy.AlgebraicNumber(x**2 - 2, 1)
add = a + b
check_comparison(add, zero, "==", (add == zero), True)

a = polypy.AlgebraicNumber(x**2 - 2, 1)
b = polypy.AlgebraicNumber(x**2 - 3, 1)
ab = polypy.AlgebraicNumber(x**4 + -10*x**2 + 1, 3)
add = a + b
check_comparison(add, ab, "==", (add == ab), True)
