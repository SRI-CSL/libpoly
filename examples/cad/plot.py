#!/usr/bin/env python

import polypy
import cad

import matplotlib.pyplot as plt

# 2D plotting of polynomials
class PolyPlot2D(cad.CylinderNotify):
            
    # Initialize   
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.cad = cad.CAD([x, y])
        self.polynomials = []
        self.cylinders = []

    # Add a polynomial
    def add_polynomial(self, f, sign_condition):
        self.polynomials.append((f, sign_condition))

    # Show the plot
    def show(self):
        # Run the CAD and collect cyllinders
        mycad = cad.CAD([x, y])
        mycad.cylinder_notify = self
        for (f, sign_condition) in self.polynomials: 
            mycad.add_polynomial(f, sign_condition)
        mycad.run()
        
    # Notifications on sylinders
    def cylinder_notify(self, cylinder, assignment):
        self.cylinders.append(cylinder)
        print("Cylinder:\n", cylinder)

if __name__ == "__main__":
    # Some variables
    x = polypy.Variable("x");
    y = polypy.Variable("y");
    # Setup   
    plot = PolyPlot2D(x, y) 
    plot.add_polynomial(x**2 + y**2 - 1, polypy.SGN_EQ_0)
    # Show
    plot.show()
        
        
