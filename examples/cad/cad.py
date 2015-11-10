#!/usr/bin/python

import polypy

import itertools

#
# Get all the reductums (including the polynomial itself), but not the constant
#
def get_reductums(f, x):
    R = []
    while f.var() == x:
        R.append(f)
        f = f.reductum()
    return R

#
# Do CAD
# 
class CAD:    
  
    # Map from variables to polynomials
    projection_map = None
    
    # Set of variables we're working with
    variables = None
    
    # Initialize
    def __init__(self, variable_list):
        # Set the variable order
        polypy.variable_order.set(variable_list)
        # Map from variables to sets of polynomials
        self.projection_map = {}
        for x in variable_list:
            self.projection_map[x] = set()
        # Variables 
        self.variables = variable_list
    
    #
    # Add a polynomial to CAD
    # 
    def add_polynomial(self, f):
        # Factor the polynomial
        for (f_factor, f_facto_multiplicity) in f.factor_square_free():
            # Add non-constant polynomials
            if (f_factor.degree() > 0):
                # Top variable
                x = f_factor.var()
                # Add to projection map
                self.projection_map[x].add(f_factor)
    
    #
    # Add a collection of polynomials to CAD
    #
    def add_polynomials(self, fs):
        for f in fs:
            self.add_polynomial(f) 
                       
    # 
    # Project. Go down the variable stack and project:
    # [1] coeff(f) for f in poly[x]
    # 
    # [3] psc(g1, g2) for f1, f2 in poly[x], g1 in R(f1, x), g2 in R(f2, x)  
    # 
    def project(self):
        for x in reversed(self.variables):
            # Project variable x
            x_poly_set = self.projection_map[x]
            # [1] coeff(f) for f in poly[x]
            for f in self.projection_map[x]:
                self.add_polynomials(f.coefficients())
            # [2] psc(g, g') for f in poly[x], g in R(f, x)
            for f in self.projection_map[x]:
                for g in get_reductums(f, x):
                    g_d = f.derivative()
                    if (g_d.var() == x):
                        self.add_polynomials(g.psc(g_d))
            # psc(g1, g2) for f1, f2 in poly[x], g1 in R(f1, x), g2 in R(f2, x)
            for (f1, f2) in itertools.combinations(self.projection_map[x], 2):
                f1_R = get_reductums(f1, x)
                f2_R = get_reductums(f2, x)
                for (g1, g2) in itertools.product(f1_R, f2_R):
                    self.add_polynomials(g1.psc(g2))

    #
    # Print internal state
    # 
    def print_state(self):
        print "Variables:", self.variables
        for x in reversed(self.variables):
            print x, ":", self.projection_map[x]
      
if __name__ == "__main__":
    # Some variables
    x = polypy.Variable("x");
    y = polypy.Variable("y");
    # Setup CAD  
    cad = CAD([x, y]) 
    cad.add_polynomial(x**2 + y**2 - 1)
    cad.add_polynomial((x-1)**2 + y**2 - 1)
    # Project
    cad.project()
    cad.print_state()