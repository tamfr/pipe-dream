# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 14:19:08 2014

@author: Mott

Numerical root finder using Newton's method for systems of equations with form Ax+Bx^2+C.
"""
# Imports necessary math modules for linear algebra

from __future__ import division
from math import floor, log10, pi
from numpy import matrix, zeros, ones, multiply
from pysnatch import xlsnatch
import os
mypath = os.path.dirname(__file__) # Makes available the absolute path to the directory of this file such that if it is imported and gets executed outside its original directory, any subsequently imported files can be found.

#import matplotlib.pyplot as plt

def pipe():
    # Define coefficients
    #a = matrix('[1,0,0;0,1,0;0,0,1]')
    #b = matrix('[0,0,0;1,0,0;0,1,0]')
    #c = matrix('[-15;20;10]')
        
    A = matrix(zeros((68,68)))
    B = matrix(zeros((68,68)))
    C = matrix(zeros((len(A),1)))
    
    # Define constants
    Diameters = xlsnatch(os.path.join(mypath,'coefficients.xlsx'),0,12,2,1,44) # Snatch diameters for later calculations
    rho = 1000 # density of water [kg/m^3]
    #mu32 = 32*0.000798 # viscosity of water @ ~30 centigrade [mPa*s] multiplied by 32 
    #g = 9.81 # Gravitational acceleration [m/s^2]
    
    num_pressure_unknowns = 25 # Number of unknown pressures
     
    # Import matrices of coefficients
    A = xlsnatch(os.path.join(mypath,'coefficients.xlsx'),1,3,2,68,68) # Imports data from Excel file by using absolute path to file.
    B = xlsnatch(os.path.join(mypath,'coefficients.xlsx'),2,3,2,68,68)
    C = xlsnatch(os.path.join(mypath,'coefficients.xlsx'),3,3,2,68,1)
    
    # Set coefficients for equation Ax+Bx^2+C
    A = A
    B = B
    C = C
    
    # Deinfe function f to solve for roots
    f=lambda x: A*x+B*multiply(x,x)+C
    
    # Initiate Jacobian based on number of equations as defined by matrix A
    J = matrix(zeros((A.shape[0],A.shape[1])))
    
    # Define parameters to run loop
    n = 0 # Define iteration counter
    max_iterations = 600 # Define max number of iterations
    tol = 2e-11 # Define tolerance upon which to interate
    
    # Initial guesses
    vel_guess = .0000003 # Good guest found to be .0000003 to yield all positive roots.
    pressure_guess = 103421.25
    
    X = matrix(ones((A.shape[0],1))*pressure_guess)
    X[num_pressure_unknowns:X.shape[0],0] = vel_guess
    #X = matrix('[-5;7;-3]')
    
    # Newton's method for finding roots
    while any(abs(v) > tol for v in f(X)):
    
        # propagate jacobian (A+Bx)
        for i in range(0,J.shape[0]):
            for j in range(0,J.shape[1]):
                J[i,j] = A[i,j] + B[i,j]*X[j]
            
        X = X - J.I*f(X)
        
        #X = matrix([(abs(X[x,0]) if X[x,0] < 0 else X[x,0]) for x in range(0,len(X))]).T
        
        if n < max_iterations:
            n = n + 1
        else:
            tol = 999999999999999999
    pressure_precision = 10 # Define number of sig figs  
    velocity_precision = 5 # Define number of sig figs 
    ans = []    
    ans[0:num_pressure_unknowns] = [round(x,-int(floor(log10(abs(x)/10**pressure_precision)))) for x in X[0:num_pressure_unknowns,0]] # Round numbers to certain number of sig figs based on precision      
    ans[num_pressure_unknowns:X.shape[0]] = [round(x,-int(floor(log10(abs(x)/10**velocity_precision)))) for x in X[num_pressure_unknowns:X.shape[0],0]] # Round numbers to certain number of sig figs based on precision
    ans.append(round(C[len(C)-1,0]/Diameters[0,len(C)-num_pressure_unknowns]**2,-int(floor(log10(abs(C[len(C)-1,0]/Diameters[0,len(C)-num_pressure_unknowns]**2)/10**velocity_precision))))) # Append the given velocity
    
    m_dot = [round(x,-int(floor(log10(abs(x)/10**velocity_precision)))) for x in multiply(multiply(Diameters[0,:-1].T,Diameters[0,:-1].T),X[num_pressure_unknowns:,0])*pi*rho/4] # Calculate mass flow rates in each branch.
    m_dot.append(round(ans[len(ans)-1]*Diameters[0,len(C)-num_pressure_unknowns]**2*pi*rho/4,-int(floor(log10(abs(ans[len(ans)-1]*Diameters[0,len(C)-num_pressure_unknowns]**2*pi*rho/4)/10**velocity_precision))))) # Append the calculated mass flow rate from given velocity
    
    #print X
    #print f(X)
    #V = matrix(range(0,43))
    #plt.plot(V, X[25:,0].T, linestyle='--',marker='o',color='g')
    #plt.show()
    variables = {}
    variable_names = xlsnatch(os.path.join(mypath,'coefficients.xlsx'),1,2,2,1,69) #
    
    for i in range(0,len(ans)):
        variables[variable_names[0,i].encode('ascii','ignore')] = ans[i]
        if i >= num_pressure_unknowns:
            variables['m_'+variable_names[0,i].encode('ascii','ignore')[2:]] = m_dot[i-num_pressure_unknowns]

    return variables