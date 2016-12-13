#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Sharhad Bashar
ECSE 543
Assignment 3
interpolate.py
Contains the functions to perform Lagrange Interpolation and Cubic Hermite
"""
from sympy import *
def diff(f,x):
    n = len(x)
    y = []
    for i in range(n):
        if (i == n - 1):
            y.append((f[i]-f[i-1])/(x[2]-x[1]))
        else:
            y.append((f[i+1]-f[i]/(x[2]-x[1])))
    return y

def interLag(pointsX, pointsY):
    x = symbols('x')    
    n = len(pointsX)
    Ljx = [None for j in range (n)]    
    for j in range (n):
        Fjx = 1
        Fjxj = 1
        yx = 0
        for k in range (n):
            if (k != j):
                Fjx *= (x - pointsX[k]) 
                Fjxj *= (pointsX[j] - pointsX[k])
                Ljx[j] = Fjx.expand()/Fjxj
    # Now we have al the Lj(x)
    for j in range (n):
        yx += pointsY[j] * Ljx[j]             
    return yx 

def cubicHer(pointsX, pointsY):
    x = symbols('x')    
    n = len(pointsX)
    Ljx = [None for j in range (n)]
    LjxPrime = [None for j in range (n)]
    Ujx = [None for j in range (n)]
    Vjx = [None for j in range (n)] 
    b = [None for j in range (n)] 
    yx = 0
    for j in range (n):
        Fjx = 1
        Fjxj = 1
        for k in range (n):
            if (k != j):
                Fjx *= (x - pointsX[k]) 
                Fjxj *= (pointsX[j] - pointsX[k])
                Ljx[j] = Fjx.expand()/Fjxj
                #now we have Lj(x)
    for j in range (n):
        LjxPrime[j] = Ljx[j].diff(x)
        #now we have Lj'(x)
    for j in range(n):
        Ujx[j] = ((1 - 2 * (LjxPrime[j] * (x - pointsX[j]))) * (Ljx[j] * Ljx[j])).expand() 
        #now we have Uj(x)
        Vjx[j] = ((x - pointsX[j]) * (Ljx[j] * Ljx[j])).expand()
        #now we have Uj(x)
        if (j < 5):    
            b[j] = (pointsY[j + 1] - pointsY[j]) / (pointsX[j + 1] - pointsX[j])
        elif (j == 5):
            b[j] = pointsY[j]/pointsX[j]
        #now we have the bj
    for j in range (n):
        yx += pointsY[j] * Ujx[j] + b[j] * Vjx[j] 
    return yx    
    
B = [0.0, 0.2,  0.4,  0.6,  0.8,   1.0,   1.1,   1.2,   1.3,   1.4,    1.5,    1.6,    1.7,    1.8,     1.9]
H = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4, 256.2, 348.7, 540.6, 1062.8, 2318.0, 4781.9, 8687.4, 13924.3, 22650.2]
Bb = [0.0, 1.3,   1.4,    1.7,    1.8,     1.9]
Hb = [0.0, 540.6, 1062.8, 8687.4, 13924.3, 22650.2]
    
a = interLag (B[:6],H[:6])
b = interLag (Bb,Hb)

print ("First Six points: " + str(a))
print ('\n')
print ("Part b: " + str(b))

cubic = (cubicHer(Bb, Hb))
print cubic



