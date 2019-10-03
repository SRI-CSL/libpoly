#!/usr/bin/env python

import polypy

import itertools
from builtins import zip

# Get all the reductums (including the polynomial itself), but not the constants
def get_reductums(f, x):
  R = []
  while f.var() == x: R.append(f); f = f.reductum()
  return R

# Add polynomials to projection map
def add_polynomial(poly_map, f):
  # Factor the polynomial
  for (f_factor, f_factor_multiplicity) in f.factor_square_free():
    # Add non-constant polynomials
    if (f_factor.degree() > 0):
      # Add to set of the top variable
      x = f_factor.var()
      if (x not in poly_map): poly_map[x] = set()
      poly_map[x].add(f_factor)
          
# Add a collection of polynomials to projection map
def add_polynomials(poly_map, polys):
  for f in polys: add_polynomial(poly_map, f) 
                       
# Project. Go down the variable stack and project:
def project(poly_map, vars):
  for x in reversed(vars):
    # Project variable x
    x_polys = poly_map[x]
    # [1] coeff(f) for f in poly[x]
    for f in x_polys:
      add_polynomials(poly_map, f.coefficients())
    # [2] psc(g, g') for f in poly[x], g in R(f, x)
    for f in x_polys:
      for g in get_reductums(f, x):
        g_d = f.derivative()
        if (g_d.var() == x):
          add_polynomials(poly_map, g.psc(g_d))
    # [3] psc(g1, g2) for f1, f2 in poly[x], g1 in R(f1, x), g2 in R(f2, x)
    for (f1, f2) in itertools.combinations(x_polys, 2):
      f1_R = get_reductums(f1, x)
      f2_R = get_reductums(f2, x)
      for (g1, g2) in itertools.product(f1_R, f2_R):
        add_polynomials(poly_map, g1.psc(g2))        
            
# Lift the first variable, update the assignment and lift recursively
def lift_first_var(poly_map, vars, m):
  # We've tried all vars, print the asignment
  if len(vars) == 0:
    print(m) # Evaluation point
    return        
  # The variable we're assigning
  x = vars[0] 
  # Get the roots
  roots = set()  # Set of root values
  for f in poly_map[x]:
    f_roots = f.roots_isolate(m)
    roots.update(f_roots)
  # Sort the roots and add infinities
  roots = [polypy.INFINITY_NEG] + sorted(roots) + [polypy.INFINITY_POS]
  # Select values in the regions, and go recursive
  r_i, r_j = itertools.tee(roots)
  next(r_j)
  for r1, r2 in zip(r_i, r_j):
    # Get the sector (r1, r2)          
    v = r1.get_value_between(r2);
    m.set_value(x, v)
    # Go recursive if m doesn't invalidate any constraints
    lift_first_var(poly_map, vars[1:], m)            
    # Get the section [r2]
    if r2 != polypy.INFINITY_POS:
      m.set_value(x, r2)
      lift_first_var(poly_map, vars[1:], m)                        
    m.unset_value(x)
    
# Do the lifting
def lift(poly_map, vars):  
  m = polypy.Assignment()
  lift_first_var(poly_map, vars, m)
                 
# Run the CAD construction
def cad(polys, vars):
  polypy.variable_order.set(vars)
  poly_map = {}
  add_polynomials(poly_map, polys)
  project(poly_map, vars)
  lift(poly_map, vars) 
      
if __name__ == "__main__":
  x = polypy.Variable("x");
  y = polypy.Variable("y");
  vars = [x, y]
  polys = [x**2 + y**2 - 1]    
  cad(polys, vars)
