#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Sharhad Bashar
260519664
ECSE 543
Assignment 3
gaussLegendreOnePoint.py
Computing the Integral using Gauss Legendre One Point method
Dec 5th, 2016
"""
import math

def integral(func1,func2,n,a,b):
    segments = n
    n,a,b = float(n),float(a),float(b)
    wi = 2.0
    xi = 0.0
    summation = 0.0
    segSize = (b - a) / float(segments)
    if (func2 == 1):
        while (a < b):
            lowLim = a
            a += segSize
            highLim = a  
            summation += (highLim - lowLim) * func1((highLim + lowLim)/2.0 * xi + (highLim + lowLim)/2.0)         
    else:
        while (a < b):
            lowLim = a
            a += segSize
            highLim = a  
            summation += (highLim - lowLim) * func2 (0.2 * func1((highLim + lowLim)/2.0 * xi + (highLim + lowLim)/2.0))         
    return summation
    
def integralUneven(func1,func2, a, b):
    relativeWidths = [2**i for i in range (0,10)]
    #[1,2,4,8,16,32,64,128,256,512]
    a,b = float(a),float(b)
    scale = (b - a) / sum(relativeWidths)
    widths = [width * scale for width in relativeWidths]
    runningWidth = 0
    runningSum = 0
    if (func2 == 1):
        for w in widthso:
            runningSum += func1(w/2 + runningWidth) * w
            runningWidth += w
    else:
        for w in widths:
            runningSum += func2(0.2*abs(func1(w/2 + runningWidth))) * w
            runningWidth += w
    return runningSum
    
print ("Integral of sin(x): 0.45970")
for i in range (1,21): 
    error = 0.45970
    integ =  integral (math.sin,1,i,0,1)
    print ("i = " + str(i) + " Integral = " + str(integ) + " Error: " + str(integ - error))
print ('')
print ("Integral of ln(x): -1") 
for i in range (1,21):  
    error = -1
    integ =  integral (math.log,1,i * 10,0,1)
    oprint ("i = " + str(i*10) + " Integral = " + str(integ) + " Error: " + str(integ - error))
print ('')
print ("Integral of ln(0.2(|sin(x)|): -2.66616") 
for i in range (1,21):    
    error = -2.66616
    integ =  integral(math.sin,math.log,i,0,1)
    print ("i = " + str(i*10) + " Integral = " + str(integ) + " Error: " + str(integ - error))    
print ('')
print ('Unever integral')
print ("Integral of ln(x): ")    
integ =  integralUneven (math.log,1,0,1)
error = -1
print ("Integral = " + str(abs(integ)) + " Error: " + str(integ - error))
print ('')
print ("Integral of ln(0.2(|sin(x)|): ") 
integ = integralUneven(math.sin,math.log,0,1)
error = -2.66616
print ("Integral = " + str(integ) + " Error: " + str(integ - error))
print ('')



