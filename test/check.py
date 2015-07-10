#!/usr/bin/env python

import argparse
import polypy
import polypy_test

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--sympy', action="store_true")
parser.add_argument('--stats', action="store_true")

args = parser.parse_args()

tests = ["python/upolynomial_gcd.py",
         "python/upolynomial_factor.py",
         "python/upolynomial_roots.py",
         "python/algebraic_number.py",
         "python/variable.py",
         "python/polynomial_arithmetic.py",
         "python/polynomial_sgn.py",
         "python/polynomial_gcd.py",
         "python/polynomial_factorization.py",
         "python/polynomial_eval.py",
         "python/polynomial_roots.py",
         "python/polynomial_resultants.py",
         "python/polynomial_feasibility.py"]

if (args.sympy):
    print "Sympy checking enabled"
    polypy_test.sympy_checker.enabled = True
else:
    print "Sympy checking disabled"
    polypy_test.sympy_checker.enabled = False
        
for test in tests:    
    print "Running", test, ":"
    context = dict()
    execfile(test, context, context)
    module = context["polypy_test"]
    print "PASS:", module.PASS
    print "FAIL:", module.FAIL

if (args.stats):    
    print "Statistics:"
    polypy.stats_print()
