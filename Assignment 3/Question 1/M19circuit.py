#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Sharhad Bashar
ECSE 543
Assignment 3
M19circuit.py
Contains the functions to calculate the flux of the M19 circuit
"""

HvB = [[0.0, 0.0],[0.2, 14.7],[0.4,36.5],[0.6,71.7],[0.8,121.4],[1,197.4],[1.1,256.2],[1.2,348.7],[1.3,540.6],[1.4,1062.8],[1.5,2318.0],[1.6,4781.9],[1.7,8687.4],[1.8,13924.3],[1.9,22650.2]]

def Hval(flux):
    B = flux/(1.0/pow(100,2))
    # interpolate for values outside domain
    if B > HvB[-1][0]:
        slope = (HvB[-1][1] - HvB[-2][1]) / (HvB[-1][0] - HvB[-2][0])
        return (B - HvB[-1][0]) * slope + HvB[-1][1]
	
    for i in range(len(HvB)):
        if HvB[i][0] > B:
            slope = (HvB[i][1] - HvB[i-1][1]) / (HvB[i][0] - HvB[i-1][0])
            return (B - HvB[i-1][0]) * slope + HvB[i-1][1]
    else: # must be smaller
        slope = (HvB[1][1] - HvB[0][1]) / (HvB[1][0] - HvB[0][0])
        return (B - HvB[0][0]) * slope + HvB[0][1]
	

def Hder(flux):
	B = flux/(1.0/pow(100,2))
	if B > HvB[-1][0]:
		return (HvB[-1][1] - HvB[-2][1]) / (HvB[-1][0] - HvB[-2][0])
	
	for i in range(len(HvB)):
		if HvB[i][0] > B:
			slope = (HvB[i][1] - HvB[i-1][1]) / (HvB[i][0] - HvB[i-1][0])
			return slope
	else: # must be smaller
		slope = (HvB[1][1] - HvB[0][1]) / (HvB[1][0] - HvB[0][0])
		return slope
		
def fFlux(flux):
	return 3.978873577e7 * flux + 0.3 * Hval(flux) - 8000
	
def fFluxDer(flux):
	return 3.978873577e7 + 0.3 * Hderivative(flux)/(1.0/pow(100,2))

	
def newRap(x, tolerance):
	i = 0
	while abs(fFlux(x)/fFlux(0)) > tolerance:
		i += 1
		x -= fFlux(x)/fFluxDer(x)
	
	print('Iterations: ' + str(i) + ' Flux: ' + str(x))
	return x

def fSubstitution(flux):
    return 8000/(39.78873577e6 + 0.3 * Hval(flux)/flux)
 
def succesSub (x,tolerance):
    i = 0
    while abs(fFlux(x)/fFlux(0)) > tolerance:
        i += 1
		x = fSubstitution(x) 
    print('Iterations: ' + str(i) + ' Flux: ' + str(x))
    return x
	
NR = newtonRaph(0.0, 1e-6)

