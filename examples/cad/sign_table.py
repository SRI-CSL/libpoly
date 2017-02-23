#!/usr/bin/python
import sys
import polypy
import itertools

def sign_table(x, polys, assignment, width = 15):
    # Get the roots    
    roots_set = set()  # Set of root values
    roots_list = []    # List of (root_value, value_definition)
    sys.stdout.write("{}".format("poly/int").center(width))
    for f in polys:
        sys.stdout.write("{}".format(f).center(width))
        f_roots = f.roots_isolate(assignment)
        roots_list.extend([v for (k, v) in enumerate(f_roots) if v not in roots_set])
        roots_set.update(f_roots)
    sys.stdout.write("\n")
    # Sort the roots and eliminate duplicates
    roots_sorted = [polypy.INFINITY_NEG] + sorted(roots_list) + [polypy.INFINITY_POS]                
    root_i, root_j = itertools.tee(roots_sorted)
    next(root_j)
    for r1, r2 in itertools.izip(root_i, root_j):        
        sys.stdout.write("({:.2f} {:.2f})".format(r1.to_double(), r2.to_double()).center(width))
        # The interval (r1, r2)          
        v = r1.get_value_between(r2);
        assignment.set_value(x, v)
        for f in polys: sys.stdout.write("{}".format(f.sgn(m)).center(width))
        sys.stdout.write("\n") 
        assignment.unset_value(x)
        # The interval [r2]
        if r2 != polypy.INFINITY_POS:
            sys.stdout.write("{:.2f}".format(r2.to_double()).center(width))
            assignment.set_value(x, r2)
            for f in polys: sys.stdout.write("{}".format(f.sgn(m)).center(width))
            sys.stdout.write("\n") 
            assignment.unset_value(x)
                  
if __name__ == "__main__":
    # Some variables
    x = polypy.Variable("x");
    m = polypy.Assignment()
    # Print sign table
    polys = [x**2 - 2, x**2 - 3]
    sign_table(x, polys, m)
