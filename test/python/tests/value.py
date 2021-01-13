import polypy
import polypy_test

from polypy import x
from itertools import combinations_with_replacement 

polypy_test.init()

def check_comparison(x, expected, precision = 0.0000001):
    diff = x.to_double() - expected
    if (diff < -precision or diff > precision):
        print("x = {0}".format(x))
        print("expected = {0}".format(expected))
        print("precision = {0}".format(precision))
        polypy_test.check(False)
    else:
        polypy_test.check(True)

polypy_test.start("Construction")

zero_int = polypy.Value(0)
one_int = polypy.Value(1)
two_int = polypy.Value(2)
three_int = polypy.Value(3)

check_comparison(zero_int, 0)
check_comparison(one_int, 1)
check_comparison(two_int, 2)
check_comparison(three_int, 3)

p = x**2 - 2
sqrt2_neg = polypy.Value(polypy.AlgebraicNumber(p, 0))
sqrt2_pos = polypy.Value(polypy.AlgebraicNumber(p, 1))

check_comparison(sqrt2_neg, -1.41421356237)
check_comparison(sqrt2_pos, 1.41421356237)

p = x**2 - 3
sqrt3_neg = polypy.Value(polypy.AlgebraicNumber(p, 0))
sqrt3_pos = polypy.Value(polypy.AlgebraicNumber(p, 1))

check_comparison(sqrt3_neg, -1.73205080757)
check_comparison(sqrt3_pos, 1.73205080757)

half_rat = one_int / two_int
third_rat = one_int / three_int

check_comparison(half_rat, 0.5)
check_comparison(third_rat, 0.333333333333)

values = [zero_int, one_int, two_int, three_int, half_rat, third_rat, sqrt2_pos, sqrt2_neg, sqrt3_pos, sqrt3_neg]

polypy_test.start("Addition")

for v1, v2 in combinations_with_replacement(values, 2):
    result_value = v1 + v2
    result_float = v1.to_double() + v2.to_double()
    check_comparison(result_value, result_float)

polypy_test.start("Subtraction")

for v1, v2 in combinations_with_replacement(values, 2):
    result_value = v1 - v2
    result_float = v1.to_double() - v2.to_double()
    check_comparison(result_value, result_float)

polypy_test.start("Multiplication")

for v1, v2 in combinations_with_replacement(values, 2):
    result_value = v1 * v2
    result_float = v1.to_double() * v2.to_double()
    check_comparison(result_value, result_float)

polypy_test.start("Division")

for v1, v2 in combinations_with_replacement(values, 2):
    if v2 == zero_int:
        continue
    result_value = v1 / v2
    result_float = v1.to_double() / v2.to_double()
    check_comparison(result_value, result_float)

polypy_test.start("Power")

for v in values:
    for p in range(0, 10):
        result_value = v**p
        result_float = v.to_double() ** p
        check_comparison(result_value, result_float)

