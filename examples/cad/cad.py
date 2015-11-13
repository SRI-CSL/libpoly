#!/usr/bin/python

import polypy
import itertools

# Definition of a cylindrical region
class Cylinder:
    
    # Variables (ordered)    
    variables = []    
    # List of sections/sector per variable
    definitions = []
    # Variable values in the region
    values = []

    # Add a on top of the cylinder definition
    def push(self, x, definition, value):
        self.variables.append(x)
        self.definitions.append(definition)
        self.values.append(value)
        
    # Pop the top of cylinder definition
    def pop(self):
        self.values.pop()
        self.definitions.pop()
        self.variables.pop()

    # For printouts
    def __str__(self):
        out = ""
        for (x, definition, value) in zip(self.variables, self.definitions, self.values):
            out = out + "%s -> %s in %s\n" % (x, value, definition)
        return out 
        
# Notifications for the cylinders of CAD
class CylinderNotify(object):
        
    # A cylinder and an assignment in the cylinder
    def notify(self, cylinder, assignment):
        print "Assignmnent: ", assignment
        print "Cylinder:\n", cylinder
        pass

    
# Do cylindrical algebraic decomposition (CAD)
class CAD:    
  
    # Map from variables to polynomials
    projection_map = {}
    
    # Signs of relevant polynomials, mapping top variables to polynomials (for lifting)
    sign_condition_map = {}
    
    # Set of variables we're working with
    variables = []
    
    # Notification
    region_notify = CylinderNotify()
    
    # Polynomials with sign conditions are considered a conjunction
    CONJUCTIVE  = 0
    # Polynomials with sign conditions are considered a disjunction
    DISJUNCTIVE = 1
    
    # Type of CAD
    type = CONJUCTIVE
    
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
    
    # Get all the reductums (including the polynomial itself), but not the constant
    @staticmethod
    def get_reductums(f, x):
        R = []
        while f.var() == x:
            R.append(f)
            f = f.reductum()
        return R
    
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
        if self.type == CAD.CONJUCTIVE:
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
    # [1] coeff(f) for f in poly[x]
    # [2] psc(g, g') for f in poly[x], g in R(f, x)
    # [3] psc(g1, g2) for f1, f2 in poly[x], g1 in R(f1, x), g2 in R(f2, x)  
    def project(self):
        for x in reversed(self.variables):
            # Project variable x
            x_poly_set = self.projection_map[x]
            # [1] coeff(f) for f in poly[x]
            for f in self.projection_map[x]:
                self.add_polynomials(f.coefficients())
            # [2] psc(g, g') for f in poly[x], g in R(f, x)
            for f in self.projection_map[x]:
                for g in CAD.get_reductums(f, x):
                    g_d = f.derivative()
                    if (g_d.var() == x):
                        self.add_polynomials(g.psc(g_d))
            # [3] psc(g1, g2) for f1, f2 in poly[x], g1 in R(f1, x), g2 in R(f2, x)
            for (f1, f2) in itertools.combinations(self.projection_map[x], 2):
                f1_R = CAD.get_reductums(f1, x)
                f2_R = CAD.get_reductums(f2, x)
                for (g1, g2) in itertools.product(f1_R, f2_R):
                    self.add_polynomials(g1.psc(g2))
    
    # Lift the first variable, update the assignment and lift recursively
    def lift_first_var(self, variables, assignment, cylinder):
        # We've tried all variables, this assignment checks out
        if len(variables) == 0:
            self.region_notify.notify(cylinder, assignment)
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
        roots_sorted = [(polypy.INFINITY_NEG, polypy.INFINITY_NEG)] + roots_sorted + [(polypy.INFINITY_POS, polypy.INFINITY_POS)]
        root_i, root_j = itertools.tee(roots_sorted)
        next(root_j)
        for r1, r2 in itertools.izip(root_i, root_j):
            # Get the sector (r1, r2)          
            v = r1[0].get_value_between(r2[0]);
            cylinder.push(x, (r1[1], r2[1]), v)
            assignment.set_value(x, v)
            # Go recursive if assignment doesn't invalidate any constraints
            if self.check_assignment(x, assignment):
                self.lift_first_var(variables[1:], assignment, cylinder)            
            assignment.unset_value(x)
            cylinder.pop()
            # Get the section [r2]
            if r2[0] != polypy.INFINITY_POS:
                cylinder.push(x, (r2[1],), v)
                assignment.set_value(x, r2[0])
                # Go recursive if assignment doesn't invalidate any constraints
                if self.check_assignment(x, assignment):
                    self.lift_first_var(variables[1:], assignment, cylinder)                        
                assignment.unset_value(x)
                cylinder.pop()
    
    # Do the lifting
    def lift(self):
        assignment = polypy.Assignment()
        cylinder = Cylinder()
        self.lift_first_var(self.variables, assignment, cylinder)
                                    
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
    cad.add_polynomial(x**2 + y**2 - 1, polypy.SGN_LE_0)
    # Project
    cad.project()
    # Lift
    cad.lift()