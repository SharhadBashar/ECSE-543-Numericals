 #!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Sharhad Bashar
260519664
ECSE 543
Assignment 2
findCapacitance.py
Computes the Capacitance per legnth of the cable
"""
#Function Imports
import numpy as np #numpy, for reading the output file
from basicDefinitions import matTranspose, matrixMult #multiplication and transpose functions

#variable declaration
totalEnergy = 0.0
nodeHeight = 6
nodeWidth = 6
pot = 110.0
UCon = [0 for x in range (4)]
#S, which was calculated in Q1
S = [[1.0,-0.5,0.0,-0.5],
     [-0.5,1.0,-0.5,0.0],
     [0.0,-0.5,1.0,-0.5],
     [-0.5,0.0,-0.5,1.0]]
###############################################################################
#Create the Mesh of voltages from the output file of Simple2d_m
mesh = [[0 for x in range (nodeHeight)] for y in range(nodeWidth)]
values = np.loadtxt('outfile.txt', dtype = float, skiprows = 3, unpack= False)
for i in range (len(values)):
    xNode = int(values[i][1]/0.02)
    yNode = int(values[i][2]/0.02)
    mesh[yNode][xNode] = values[i][3]
#Voltage of the inner conductor
mesh[0][0] = pot
mesh[0][1] = pot
###############################################################################
#Capacitance per lenght calculator
#c = eplison0*|divU|^2/V^2
#Total energy: |divU|^2 = Ucon^T * S * Ucon:
for y in range (nodeHeight - 1):
    for x in range(nodeWidth - 1):
        UCon[0] = mesh[y][x]
        UCon[1] = mesh[y][x + 1]
        UCon[2] = mesh[y + 1][x + 1]
        UCon[3] = mesh[y + 1][x]
        UConT = matTranspose(UCon)       
        potential = matrixMult(matrixMult(UCon,S),UConT)
        print potential
        totalEnergy += potential[0][0] 
#Calculates the Capacitance per length
#multiplies by 4, to include all 4 quadrants
epsilon = 8.854187817620e-12
Vsquared = pot*pot
capPerLen = (totalEnergy * epsilon * 4)/ Vsquared
print (str(capPerLen) + ' F/m')





    
    