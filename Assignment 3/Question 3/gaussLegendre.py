#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 03:02:58 2016

@author: Sharhad
"""
from numpy import *
import math

def legPoly (n,xi):
    if (n == 0):
        Pn = 1.0
    elif (n == 1):
        Pn = xi
    else:
        n = float(n)
        Pn = ((2 * n - 1) * xi * legPoly(n - 1,xi)) / n - ((n - 1) * legPoly(n - 2,xi)) / n
    return Pn

def dlegPoly (n,xi):
    dPn = (n / (xi ** 2 - 1)) * (xi * legPoly(n,xi) - legPoly(n - 1,xi))  
    return dPn  

def weights (n,xi):
    n = int(n)
    xi = float(xi)
    wi = 2.0/((1.0 - xi**2) * (dlegPoly(n,xi))**2)
    return wi

def roots (n):
    root = []
    weight = []
    tolerance = 1e-20
    for i in range (1,n + 1):
        n = float(n)
        i = float(i)
        xi = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        error = 10 * tolerance
        iters = 0
        while (error > tolerance) and (iters < 5):
            dx = legPoly(n,xi) / dlegPoly(n,xi)
            xi = xi - dx
            iters = iters + 1
            error = abs(dx)
        wi = weights(n,xi)
        root.append(xi)
        weight.append(wi)
    return (root,weight)

def integralCalc (func,n,a,b):
    a = float(a)
    b = float(b)
    (xi,wi) = roots(n)
    integral = 0
    for i in range (1,n + 1):
        integral += wi[i - 1] * func(((b - a) / 2.0) * xi[i - 1] + ((b - a) / 2.0)) 
    return (((b - a) / 2) * integral)
print (roots(20))
    
#print ('Integral of sin(x):')
#for n in range (1,21):
#    print ('N = ' + str(n) + ', Integral = ' + str(integralCalc(math.sin,n,0,1)))
#print ('')
#print ('Integral of ln(x):')    
#for n in range (1,21):
#    print ('N = ' + str(n * 10) + ', Integral = ' + str(integralCalc(math.log,n*,0,1)))


















