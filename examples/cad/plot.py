#!/usr/bin/python

import polypy
import cad

import matplotlib.pyplot as plt

# 2D plotting of polynomials
class PolyPlot2D:
        
    # The variables
    x = None
    y = None
    
    # List of polynomials with sign conditions
    polynomials = None
    
    # The CAD 
    cad = None
    
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.cad = cad.CAD([x, y])
        self.polynomials = []

    # Add a polynomial
    def add_polynomial(self, f, sign_condition):
        self.polynomials.append((f, sign_condition))
        self.cad.add_polynomial(f, sign_condition)

    # Show the plot
    def show(self):
        # Do CAD
        self.cad.project()
        self.cad.lift()
        

if __name__ == "__main__":
    # Some variables
    x = polypy.Variable("x");
    y = polypy.Variable("y");
    # Setup   
    plot = PolyPlot2D(x, y) 
    plot.add_polynomial(x**2 + y**2 - 1, polypy.SGN_EQ_0)
    # Show
    plot.show()
        
        
