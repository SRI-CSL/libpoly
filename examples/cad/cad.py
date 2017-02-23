#!/usr/bin/python

import polypy
import itertools
        
# Notifications for the cylinders of CAD
class CylinderNotify(object):        
    # A cylinder_nofity and an assignment in the cylinder_nofity
    def cylinder_notify(self, cylinder, assignment):
        print "Assignmnent: ", assignment
        print "Cylinder:\n", cylinder
        pass


# Get all the reductums (including the polynomial itself), but not the constant
def get_reductums(f, x):
    R = []
    while f.var() == x:
        R.append(f)
        f = f.reductum()
    return R

# Polynomials with sign conditions are considered a conjunction
CONJUNCTIVE  = 0
# Polynomials with sign conditions are considered a disjunction
DISJUNCTIVE = 1
    
# Do cylindrical algebraic decomposition (CAD)
class CAD:    
          
    # Initialize
    def __init__(self, variable_list, type = CONJUNCTIVE):
        # Init variables
        self.projection_map = {}
        self.sign_condition_map = {}
        self.variables = variable_list
        self.cylinder_notify = CylinderNotify()
        self.type = type
        # Initialize the maps
        for x in variable_list:
            self.projection_map[x] = set()
            self.sign_condition_map[x] = set()
        # Set the variable order
        polypy.variable_order.set(variable_list)
        
    # Add a polynomial to CAD
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
    
    # Check if the assignment respects the polynomials with top variable x  
    def check_assignment(self, x, assignment):        
        if self.type == CONJUNCTIVE:
            for p, sgn_condition in self.sign_condition_map[x]:
                if not p.sgn_check(assignment, sgn_condition):
                    return False
            return True
        else:
            for p, sgn_condition in self.sign_condition_map[x]:
                if p.sgn_check(assignment, sgn_condition):
                    return True
            return False
    
    # Add a collection of polynomials to CAD
    def add_polynomials(self, polynomials):
        for f in polynomials:
            self.add_polynomial(f) 
                       
    # Project. Go down the variable stack and project:
    def project(self):
        for x in reversed(self.variables):
            # Project variable x
            x_poly_set = self.projection_map[x]
            # [1] coeff(f) for f in poly[x]
            for f in x_poly_set:
                self.add_polynomials(f.coefficients())
            # [2] psc(g, g') for f in poly[x], g in R(f, x)
            for f in x_poly_set:
                for g in get_reductums(f, x):
                    g_d = f.derivative()
                    if (g_d.var() == x):
                        self.add_polynomials(g.psc(g_d))
            # [3] psc(g1, g2) for f1, f2 in poly[x], g1 in R(f1, x), g2 in R(f2, x)
            for (f1, f2) in itertools.combinations(x_poly_set, 2):
                f1_R = get_reductums(f1, x)
                f2_R = get_reductums(f2, x)
                for (g1, g2) in itertools.product(f1_R, f2_R):
                    self.add_polynomial(g1.psc(g2))        
            
    # Lift the first variable, update the assignment and lift recursively
    def lift_first_var(self, variables, assignment, cylinder):
        # We've tried all variables, this assignment checks out
        if len(variables) == 0:
            self.cylinder_notify.cylinder_notify(cylinder, assignment)
            return        
        # Lift first variable
        x = variables[0] 
        # Get the roots
        roots_set = set()  # Set of root values
        roots_list = [] # List of (root_value, value_definition)
        for f in self.projection_map[x]:
            f_roots = f.roots_isolate(assignment)
            roots_list.extend([(v, (f, k)) for (k, v) in enumerate(f_roots) if v not in roots_set])
            roots_set.update(f_roots)
        # Sort the roots and eliminate duplicates
        roots_sorted = sorted(roots_list, key = lambda (v,d): v)
        # Go recursive on the sectors
        roots_sorted = [(polypy.INFINITY_NEG, polypy.INFINITY_NEG)] + \
                       roots_sorted + \
                       [(polypy.INFINITY_POS, polypy.INFINITY_POS)]
        root_i, root_j = itertools.tee(roots_sorted)
        next(root_j)
        for r1, r2 in itertools.izip(root_i, root_j):
            # Get the sector (r1, r2)          
            v = r1[0].get_value_between(r2[0]);
            cylinder.append({ "var" : x, "sector" : (r1[1], r2[1]), "value" : v})
            assignment.set_value(x, v)
            # Go recursive if assignment doesn't invalidate any constraints
            if self.check_assignment(x, assignment):
                self.lift_first_var(variables[1:], assignment, cylinder)            
            assignment.unset_value(x)
            cylinder.pop()
            # Get the section [r2]
            if r2[0] != polypy.INFINITY_POS:
                v = r2[0]
                cylinder.append({ "var" : x, "section" : r2[1], "value": v})
                assignment.set_value(x, v)
                # Go recursive if assignment doesn't invalidate any constraints
                if self.check_assignment(x, assignment):
                    self.lift_first_var(variables[1:], assignment, cylinder)                        
                assignment.unset_value(x)
                cylinder.pop()
    
    # Do the lifting
    def lift(self):
        assignment = polypy.Assignment()
        cylinder = []
        self.lift_first_var(self.variables, assignment, cylinder)
                 
    # Run the CAD construction
    def run(self):
        self.project()
        self.lift() 
                                   
    # Print internal state
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
    cad.add_polynomial(x**2 + y**2 - 1, polypy.SGN_GE_0)
    # Run CAD
    cad.run()
