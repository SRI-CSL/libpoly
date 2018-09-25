#!/usr/bin/env python
from __future__ import print_function
from io import open

import sys

import argparse
import polypy
import polypy_test

print('Using Python', sys.version)

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
    print("Sympy checking enabled")
    polypy_test.sympy_checker.enabled = True
else:
    print("Sympy checking disabled")
    polypy_test.sympy_checker.enabled = False

def compile_file(filename):
    """Compile a file for use with exec() for both Python 2 and 3"""
    # See: https://portingguide.readthedocs.io/en/latest/builtins.html#removed-execfile
    with open(filename, encoding='utf-8') as f:
        return compile(f.read(), filename, 'exec')

for test in tests:
    print("Running", test, ":")
    context = dict()

    exec(compile_file(test), context, context)

    module = context["polypy_test"]
    print("PASS:", module.PASS)
    print("FAIL:", module.FAIL)

if (args.stats):    
    print("Statistics:")
    polypy.stats_print()
