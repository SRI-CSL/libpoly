#!/usr/bin/env python

import polypy
import polypy_test

polypy_test.init()

def check_roots_isolate(p, expected):
    roots = p.roots_isolate()
    sorted = all(roots[i] < roots[i+1] for i in range(len(roots)-1))
    count = len(roots)
    if ((not sorted) or count != expected):
        print("p = {0}".format(p))
        print("sorted = {0}".format(sorted))
        print("count = {0}".format(count))
        print("expected = {0}".format(expected))
        polypy_test.check(False)
    else:
        polypy_test.check(True)

def check_roots_count(p, expected, lb = None, ub = None):
    if (lb is None):
        count = p.roots_count()
    else:
        count = p.roots_count(lb, ub)

    if (count != expected):
        polypy_test.check(False)
        print("p = {0}".format(p))
        print("lb = {0}".format(lb))
        print("ub = {0}".format(ub))
        print("count = {0}".format(count))
        print("expected = {0}".format(expected))
    else:
        polypy_test.check(True)

# polypy.trace_enable("roots")
# polypy.trace_enable("algebraic_number")
# polypy.trace_enable("division")
# polypy.trace_enable("sturm_sequence_check")

x = polypy.x

polypy_test.start("Root isolation")

p = x**2 - 1
check_roots_isolate(p, 2)

p = (x - 1)*(x + 1)*(x - 2)
check_roots_isolate(p, 3)

p = 1*x**5 + (-23*x**4) + (-57*x**3) + 85*x**2 + 64*x + (-63)
check_roots_isolate(p, 5)

p = 41*x**5 + (-79*x**4) + 44*x**3 + 56*x**2 + (-10*x) + 41
check_roots_isolate(p, 1)

p = 9*x**13 - 18*x**11 - 33*x**10 + 102*x**8 + 7*x**7 - 36*x**6 - 122*x**5 + 49*x**4 + 93*x**3 - 42*x**2 - 18*x** + 9
check_roots_isolate(p, 6)

polypy_test.start("Root counting")

p = x + 1 - x
check_roots_count(p, 0)
check_roots_count(p, 0, -0.5, 0.5)
check_roots_count(p, 0, -1, 1)
check_roots_count(p, 0, -1.5, 1.5)
check_roots_count(p, 0, -2, 2)
check_roots_count(p, 0, -2.5, 2.5)
check_roots_count(p, 0, -3, 3)

p = (x - 1)
check_roots_count(p, 1)
check_roots_count(p, 0, -0.5, 0.5)
check_roots_count(p, 0, -1, 1)
check_roots_count(p, 1, -1.5, 1.5)
check_roots_count(p, 1, -2, 2)
check_roots_count(p, 1, -2.5, 2.5)
check_roots_count(p, 1, -3, 3)

p = (x - 1)*(x + 1)
check_roots_count(p, 2)
check_roots_count(p, 0, -0.5, 0.5)
check_roots_count(p, 0, -1, 1)
check_roots_count(p, 2, -1.5, 1.5)
check_roots_count(p, 2, -2, 2)
check_roots_count(p, 2, -2.5, 2.5)
check_roots_count(p, 2, -3, 3)

p = (x - 1)*(x + 1)*(x - 2)
check_roots_count(p, 3)
check_roots_count(p, 0, -0.5, 0.5)
check_roots_count(p, 0, -1, 1)
check_roots_count(p, 2, -1.5, 1.5)
check_roots_count(p, 2, -2, 2)
check_roots_count(p, 3, -2.5, 2.5)
check_roots_count(p, 3, -3, 3)

p = (x - 1)*(x + 1)*(x - 2)*(x + 2)
check_roots_count(p, 4)
check_roots_count(p, 0, -0.5, 0.5)
check_roots_count(p, 0, -1, 1)
check_roots_count(p, 2, -1.5, 1.5)
check_roots_count(p, 2, -2, 2)
check_roots_count(p, 4, -2.5, 2.5)
check_roots_count(p, 4, -3, 3)

p = (x - 1)**3 * (x + 1)** 2 * (x - 2)**2 * (x + 2)
check_roots_count(p, 4)
check_roots_count(p, 0, -0.5, 0.5)
check_roots_count(p, 0, -1, 1)
check_roots_count(p, 2, -1.5, 1.5)
check_roots_count(p, 2, -2, 2)
check_roots_count(p, 4, -2.5, 2.5)
check_roots_count(p, 4, -3, 3)
