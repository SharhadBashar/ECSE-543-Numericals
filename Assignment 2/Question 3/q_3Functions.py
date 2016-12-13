# -*- coding: utf-8 -*-
"""
Sharhad Bashar
ECSE 543
Assignment 2
q_3Functions.py
Functons for generating the matrix from Assign 1
Nov 14th, 2016
"""

#using the top right quater of the cable, due to symmetry
import math
#######################################################################################################
#Generates the initial mesh, taking into considering the boundary conditions
def genMesh (h):
    cableHeight = 0.1
    cableWidth = 0.1
    coreHeight = 0.02
    coreWidth = 0.04
    corePot = 110.0
    nodeHeight = (int)(cableHeight/h + 1)
    nodeWidth = (int)(cableWidth/h + 1)    
    #Create the mesh, with Dirchlet conditions 
    mesh = [[(None,corePot) if x <= coreWidth/h and y <= coreHeight/h \
             else ((None,0.0) if x == nodeHeight - 1 or y == nodeWidth - 1 \
             else (0,0.0)) for x in range(nodeWidth)] for y in range(nodeHeight)]    
    #update the mesh to take into account the Neuman conditions
    nodeForA = 0
    for y in range (nodeHeight):
        for x in range(nodeWidth):
            if (mesh[y][x][0] == 0):
                mesh[y][x] = (nodeForA,0.0)
                nodeForA += 1            
    return mesh 
#######################################################################################################
#The Equation that calculates SOR    
def SOR(mesh,h,w):
    cableHeight = 0.1
    cableWidth = 0.1
    coreHeight = 0.02
    coreWidth = 0.04
    nodeHeight = (int)(cableHeight/h + 1)
    nodeWidth = (int)(cableWidth/h + 1)    
    for y in range (1,nodeHeight - 1):
        for x in range (1,nodeWidth - 1):
            if (x > (int)(coreWidth/h) or y > (int)(coreHeight/h)):
                mesh[y][x] = (1 - w) * mesh[y][x] + (w/4) * (mesh[y][x-1] + mesh[y][x+1] + mesh[y-1][x] + mesh[y+1][x])    
    return mesh            
#######################################################################################################
#The Equation that calculates Jacobian
def jacobian(mesh,h):
    cableHeight = 0.1
    cableWidth = 0.1
    coreHeight = 0.02
    coreWidth = 0.04
    nodeHeight = (int)(cableHeight/h + 1)
    nodeWidth = (int)(cableWidth/h + 1) 
    for y in range (1,nodeHeight - 1):    
        for x in range (1,nodeWidth - 1):
            if (x > (int)(coreWidth/h) or y > (int)(coreHeight/h)):
                mesh[y][x] = (1/4) * (mesh[y][x-1] + mesh[y][x+1] + mesh[y-1][x] + mesh[y+1][x])
    return mesh

####################################################################################################### 
#Equation that computes the residue 
def computeMaxRes(mesh,h):
    cableHeight = 0.1
    cableWidth = 0.1
    coreHeight = 0.02
    coreWidth = 0.04
    nodeHeight = (int)(cableHeight/h + 1)
    nodeWidth = (int)(cableWidth/h + 1)
    maxRes = 0
    for y in range(1, nodeHeight - 1):
        for x in range(1, nodeWidth - 1):
            if (x > coreWidth/h or y > coreHeight/h):
                #calculate the residue of each free point
                res = mesh[y][x-1] + mesh[y][x+1] + mesh[y-1][x] + mesh[y+1][x] - 4 * mesh[y][x]
                res = math.fabs(res)
                if (res > maxRes):
                    #Updates variable with the biggest residue amongst the free point
                    maxRes = res
    return maxRes
#######################################################################################################
#Function that computes the number of iterations    
def numIteration (initialMesh,h,w,method):       
    minRes = 0.0001
    iteration = 1
    if (method == 's'):    
        mesh = SOR(initialMesh,h,w)
        while (computeMaxRes(mesh,h) >= minRes):
            mesh = SOR(mesh,h,w)
            iteration += 1
    elif (method == 'j'):
        mesh = jacobian(initialMesh,h) 
        while (computeMaxRes(mesh,h) >= minRes):
            mesh = jacobian(mesh,h)
            iteration += 1
    print ('Number of iterations: '+ str(iteration))   
    return(mesh)    
#######################################################################################################
#Function that returns the potential at a free node    
def getPot(mesh, x, y, h):
    cableHeight = 0.1
    cableWidth = 0.1
    nodeHeight = (int)(cableHeight/h + 1)
    nodeWidth = (int)(cableWidth/h + 1)
    xNode = int(nodeWidth - x/h - 1)
    yNode = int(nodeHeight - y/h - 1)
    return mesh[yNode][xNode]    
#######################################################################################################    


