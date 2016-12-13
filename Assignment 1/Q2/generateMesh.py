# -*- coding: utf-8 -*-
"""
Sharhad Bashar
ECSE 543
Assignment 1
OCt 17th, 2016
generateMesh.py
"""
##########################################################################################
#Generates the incident matrix
def genMesh(meshDim):    
    N = (meshDim + 1) #N is the number of nodes in one line of the mesh
    totalNodes = N ** 2 #totla nodes in the mesh circuit
    totalBranches = 2 * N * (N - 1) #total branches in the mesh circuit
    #Createing incident matrix, and filling it with 0's
    incMat = [[0 for rows in range (totalBranches + 1)] for cols in range (totalNodes)]
    #i is the horizontal rows, j is the vertical rows    
    for i in range (1,(N + 1),1):
        for j in range (1,(N + 1),1):
            node = N * (j - 1) + i  #Numbering the nodes          
            bUp = node + (N - 1) * (j - 1) #Branch above the node  
            bDown = bUp - 1 #Branch below the node
            bLeft = bUp - N #Branch to the left of the node 
            bRight = bUp + N - 1 #Branch to the right of the node            
            #Populating the Incident Matrix
            #Leaving node = +1
            #Entering node = -1
            #Taking into account of the voltage source connected to the bottom left and 
            #top right of the mesh 
            incMat[0][totalBranches] = -1
            incMat[totalNodes - 1][totalBranches] = 1            
            #Rest of the Mesh
            if (j == 1): #left most vertical branch            
                incMat[node - 1][bRight - 1] = 1
                if (i == 1):
                    incMat[node - 1][bUp - 1] = 1
                elif (i == N):
                    incMat[node - 1][bDown - 1] = -1
                else:
                    incMat[node - 1][bUp - 1] = 1
                    incMat[node - 1][bDown - 1] = -1                    
            elif (j == N): #right most vertical branch  
                incMat[node - 1][bLeft - 1] = -1
                if(i == 1):
                    incMat[node - 1][bUp - 1] = 1
                elif (i == N):
                    incMat[node - 1][bDown - 1] = -1
                else: 
                    incMat[node - 1][bUp - 1] = 1
                    incMat[node - 1][bDown - 1] = -1                    
            else:
                incMat[node - 1][bLeft - 1] = -1
                incMat[node - 1][bRight - 1] = 1
                if (i == 1):
                    incMat[node - 1][bUp - 1] = 1
                elif (i == N):
                    incMat[node - 1][bDown - 1] = -1
                else:
                    incMat[node - 1][bUp - 1] = 1
                    incMat[node - 1][bDown - 1] = -1
    incMatR = [[0 for rows in range (totalBranches + 1)] for cols in range (totalNodes - 1)]
    for i in range (totalBranches + 1):
        for j in range (totalNodes - 1):
            incMatR[j][i] = incMat[j][i] 
    return incMatR                

##########################################################################################           
#create E vector
def EVector(meshDim):    
    N = (meshDim + 1) #N is the number of nodes in one line of the mesh
    totalBranches = 2 * N * (N - 1)    
    E = [[0.00 for rows in range(1)]for rows in range (totalBranches + 1)]
    E[totalBranches][0] = 1.00 #voltage source of 1V in the last branch   
    return E
##########################################################################################
#create J vector
def JVector(meshDim):    
    N = (meshDim + 1) #N is the number of nodes in one line of the mesh
    totalBranches = 2 * N * (N - 1)    
    J = [[0.00 for rows in range(1)]for rows in range (totalBranches + 1)]
    return J
##########################################################################################  
#create R vector
def RVector(meshDim):    
    N = (meshDim + 1) #N is the number of nodes in one line of the mesh
    totalBranches = 2 * N * (N - 1)    
    R = [[1000 for rows in range(1)]for rows in range (totalBranches + 1)]
    R[totalBranches][0] = 1 #i ohm resistance connecting the source to the mesh   
    return R
##########################################################################################
    
    
    
    
    
    