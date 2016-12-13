# -*- coding: utf-8 -*-
"""
Sharhad Bashar
ECSE 543
Assignment 1
OCt 17th, 2016
question_2.py
"""

import time
from basicDefinitions import voltageSolver
from generateMesh import genMesh, EVector, JVector, RVector

#N = int(input("Please enter the size of mesh: ")) #Lets the user enter the desired dimention of the mesh
for N in range(2,21,1):
    incidentMatrix = genMesh(N) #Generates the incident matrix for the input N
    #Generates the E, J and R vectors
    E = EVector(N)
    J = JVector(N)
    R = RVector(N)
    #Solves for the voltage across the mesh
    startTime = time.clock() #Starts the clock
    V = voltageSolver(incidentMatrix,E,J,R)
    endTime = time.clock() #Stops the clock
    #Vm = Vs*(Req/1+Req)
    Rreq = V[0]/(1-V[0]) #Solves for the equivalent resistance using voltage dividor
    
    timeTaken = endTime - startTime #Gives the execution time
    print("Equivaent Resistance for size "+ str(N) +" is: " + str(Rreq) + " ohms")
    print("Execution time: " + str(timeTaken) +"s")
    print('')

