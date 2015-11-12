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
# An element of the CAD assignment tree
#



#
# Do CAD
# 
class CAD:    
  
    # Map from variables to polynomials
    projection_map = None
    
    # Signs of relevant polynomials, mapping top variables to polynomials (for lifting)
    sign_condition_map = None
    
    # Set of variables we're working with
    variables = None
    
    # Initialize
    def __init__(self, variable_list):
        # Set the variable order
        polypy.variable_order.set(variable_list)
        # Map from variables to sets of polynomials
        self.projection_map = {}
        self.sign_condition_map = {}
        for x in variable_list:
            self.projection_map[x] = set()
            self.sign_condition_map[x] = set()
            
        # Variables 
        self.variables = variable_list
    
    #
    # Add a polynomial to CAD
    # 
    def add_polynomial(self, f, sign_condition = None):
        # Factor the polynomial
        for (f_factor, f_facto_multiplicity) in f.factor_square_free():
            # Add non-constant polynomials
            if (f_factor.degree() > 0):
                # Top variable
                x = f_factor.var()
                # Add to projection map
                self.projection_map[x].add(f_factor)
        # If sign give, remember it
        if sign_condition is not None:
            self.sign_condition_map[x].add((f, sign_condition))            
    
    # 
    # Check if the assignmnent respects the polynomials with top variable x  
    #
    def check_assignment(self, x, assignment):
        for p, sgn_condition in self.sign_condition_map[x]:
            # Check
            if not p.sgn_check(assignment, sgn_condition):
                print "Discarding ", assignment, "due to", p
                return False
        return True
    
    #
    # Add a collection of polynomials to CAD
    #
    def add_polynomials(self, fs):
        for f in fs:
            self.add_polynomial(f) 
                       
    # 
    # Project. Go down the variable stack and project:
    # [1] coeff(f) for f in poly[x]
    # [2] psc(g, g') for f in poly[x], g in R(f, x)
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
            # [3] psc(g1, g2) for f1, f2 in poly[x], g1 in R(f1, x), g2 in R(f2, x)
            for (f1, f2) in itertools.combinations(self.projection_map[x], 2):
                f1_R = get_reductums(f1, x)
                f2_R = get_reductums(f2, x)
                for (g1, g2) in itertools.product(f1_R, f2_R):
                    self.add_polynomials(g1.psc(g2))

    
    def lift_first_var(self, variables, assignment):
        # We've tried all variables, this assignment checks out
        if len(variables) == 0:
            print "Full model:", assignment
            return        
        # Lift first variable
        x = variables[0] 
        # Get the roots
        roots = set()
        for f in self.projection_map[x]:
            roots.update(f.roots_isolate(assignment))
        roots_sorted = sorted(roots)            
        # Get the evaluation points
        if len(roots_sorted) == 0:
            # No roots, just pick 0
            cad_points = [0]
        else:
            # A point smaller than all roots
            first = roots_sorted[0].get_value_between(polypy.INFINITY_NEG)
            cad_points = [first, roots_sorted[0]]
            # Points between roots
            root_i, root_j = itertools.tee(roots_sorted)
            next(root_j)
            for r1, r2 in itertools.izip(root_i, root_j):
                v = r1.get_value_between(r2);
                cad_points.extend([v, r2])
            # A point larger than all roots
            last = roots_sorted[-1].get_value_between(polypy.INFINITY_POS)
            cad_points.append(last)      
        # Go recursive
        for v in cad_points:
            # Set the assignment
            assignment.set_value(x, v)            
            # Go recursive if assignment doesn't invalidate any constraints
            if self.check_assignment(x, assignment):
                self.lift_first_var(variables[1:], assignment)            
            # Unset assignment
            assignment.unset_value(x)                    
    
    # 
    # Do the lifting
    #
    def lift(self):
        assignment = polypy.Assignment()
        self.lift_first_var(self.variables, assignment)
                                    
    #
    # Print internal state
    # 
    def print_state(self):
        print "Variables:", self.variables
        print "Projection map:"
        for x in reversed(self.variables):
            print x, ":", self.projection_map[x]
        print "Sign map:"
        for x in self.variables:
            print x, ":", self.sign_condition_map[x]
      
if __name__ == "__main__":
    # Some variables
    x = polypy.Variable("x");
    y = polypy.Variable("y");
    # Setup CAD  
    cad = CAD([x, y]) 
    cad.add_polynomial(2*x - 1, polypy.SGN_EQ_0)
    cad.add_polynomial((x-1)**2 + y**2 - 1, polypy.SGN_EQ_0)
    # Project
    cad.project()
    # Lift
    cad.lift()