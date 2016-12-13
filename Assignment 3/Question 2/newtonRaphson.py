
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Sharhad Bashar
260519664
ECSE 543
Assignment 3
newtonRaphson.py
Computing the voltages using Newton Raphson
Dec 5th, 2016
"""
import math
from basicDefinitions import matrixAddorSub, matrixMult, matInverse,is1Dor2D,matTranspose, scalarMult     
############################################################################################              
def newRap (E,R,V1,V2,Isa,Isb,k):    
    f1 = V1 - E + R * Isa * (math.exp((V1 - V2) / 0.025) - 1.0)
    f2 = Isa * ((math.exp((V1 - V2) / 0.025) - 1.0)) - Isb * (math.exp(V2 / 0.025) - 1.0)   
    print ('f1 : ' + str(f1))
    print ('f2 : ' + str(f2))
    print ('Number of iterations : ' + str(k))
    print('')
    f1V1Prime = 1.0 + (R * Isa / 0.025) * (math.exp((V1 - V2) / 0.025))
    f1V2Prime = -1.0 * (R * Isa / 0.025) * (math.exp((V1 - V2) / 0.025))
    f2V1Prime = (Isa / 0.025) * (math.exp((V1 - V2) / 0.025))
    f2V2Prime = -1.0 * (Isa / 0.025) * (math.exp((V1 - V2) / 0.025)) - (Isb / 0.025) * (math.exp(V2 / 0.025))    
    V = [[V1],[V2]]
    f = [[f1],[f2]]
    J = [[None for x in range (2)] for y in range(2)]    
    J[0][0] = f1V1Prime
    J[0][1] = f1V2Prime
    J[1][0] = f2V1Prime
    J[1][1] = f2V2Prime    
    invJ = matInverse(J)
    V = matrixAddorSub(matTranspose(scalarMult(-1.0,matrixMult(invJ,f))),V,'a') 
    print ('V1 : ' + str(V[0][0])+ ' V')
    print ('V2 : ' + str(V[1][0])+ ' V')    
    return V
############################################################################################    
E = 0.2    
R = 512.0
Isa = 0.0000006
Isb = 0.0000012
V1 = 0.0
V2 = 0.0
k = 0


print  ('V1 : ' + str(V1)+ ' V')
print ('V2 : ' + str(V2)+ ' V')

V = newRap(E,R,V1,V2,Isa,Isb,k)

f1 = V[0][0] - E + R * Isa *(math.exp((V[0][0] - V[1][0])/0.025) - 1.0)

while (f1 > 0):
    k += 1
    V = newRap(E,R,V[0][0],V[1][0],Isa,Isb,k)
    f1 = V[0][0] - E + R * Isa *(math.exp((V[0][0] - V[1][0])/0.025) - 1.0)
f1 = V[0][0] - E + R * Isa * (math.exp((V[0][0] - V[1][0]) / 0.025) - 1.0)
f2 = Isa * ((math.exp((V[0][0] - V[1][0]) / 0.025) - 1.0)) - Isb * (math.exp(V[1][0] / 0.025) - 1.0)  
print ('f1 : ' + str(f1))
print ('f2 : ' + str(f2))
print ('Number of iterations : ' + str(k + 1))




