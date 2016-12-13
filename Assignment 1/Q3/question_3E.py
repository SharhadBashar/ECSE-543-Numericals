# -*- coding: utf-8 -*-
"""
Sharhad Bashar
ECSE 543
Assignment 1
Oct 17th, 2016
question_3E.py
"""
import math
#######################################################################################################
#Generates the initial mesh, taking into considering the boundary conditions
def genMesh (verLine,horLine):
    cableHeight = 0.1
    cableWidth = 0.1
    coreHeight = 0.02
    coreWidth = 0.04
    corePot = 110.0    
    #Create the mesh, with Dirchlet conditions 
    mesh = [[corePot if x <= coreWidth and y <= coreHeight else 0.0 for x in verLine] for y in horLine]    
    #update the mesh to take into account the Neuman conditions
    rateofChangeX = 110/(cableWidth - coreWidth)
    rateofChangeY = 110/(cableHeight - coreHeight)    
    for x in range (len(verLine)):
        if (verLine[x] > coreWidth):        
            mesh[0][x] = 110 - rateofChangeX * (verLine[x] - coreWidth)            
    for y in range (len(horLine)):
        if (horLine[y] > coreHeight):        
            mesh[y][0] = 110 - rateofChangeY * (horLine[y] - coreHeight)
    return mesh 
####################################################################################################### 
#The Equation that calculates SOR    
def SOR(mesh,verLine,horLine):
    coreHeight = 0.02
    coreWidth = 0.04
    for y in range (1,len(horLine) - 1):
        for x in range(1,len(verLine) - 1):
            if (verLine[x] > coreWidth or horLine[y] > coreHeight):
                a1 = verLine[x] - verLine[x-1]
                a2 = verLine[x+1] - verLine[x]
                b1 = horLine[y+1] - horLine[y]
                b2 = horLine[y] - horLine[y-1]
                mesh[y][x] = (mesh[y][x-1]/(a1 * (a1 + a2)) + mesh[y][x+1]/(a2 * (a1 + a2)) + \
                             mesh[y-1][x]/(b1 * (b1 + b2)) + mesh[y+1][x]/(b2 * (b1 + b2))) / \
                             (1/(a1 * a2) + 1/(b1 * b2)) 
    return mesh
####################################################################################################### 
#Equation that computes the residue
def computeMaxRes(mesh,horLine,verLine):
    coreHeight = 0.02
    coreWidth = 0.04
    maxRes = 0 
    for y in range (1,len(horLine) - 1):
        for x in range (1,len(verLine) - 1):
            if (verLine[x] > coreWidth or horLine[y] > coreHeight):
                a1 = verLine[x] - verLine[x-1]
                a2 = verLine[x+1] - verLine[x]
                b1 = horLine[y+1] - horLine[y]
                b2 = horLine[y] - horLine[y-1]
                res = (mesh[y][x-1]/(a1 * (a1 + a2)) + mesh[y][x+1]/(a2 * (a1 + a2)) + mesh[y-1][x]/(b1 * (b1 + b2)) + mesh[y+1][x]/(b2 * (b1 + b2))) - (1/(a1 * a2) + 1/(b1 * b2))*mesh[y][x]
                res = math.fabs(res)
                if (res > maxRes):
                    #Updates variable with the biggest residue amongst the free point
                    maxRes = res
    return maxRes          
####################################################################################################### 
def numIteration (initialMesh,horLine,verLine):       
    minRes = 0.0001
    mesh = SOR(initialMesh,horLine,verLine)
    iteration = 1
    while (computeMaxRes(mesh,horLine,verLine) >= minRes):
        mesh = SOR(mesh,horLine,verLine)
        iteration += 1
    print ('Number of iterations: '+ str(iteration))   
    return(mesh)  
#######################################################################################################     
def getPot(mesh, x, y,verLine,horLine):
    xNode = verLine.index(x)
    yNode = horLine.index(y)
    return mesh[yNode][xNode]    
       
horLine = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
verLine = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]

x = 0.06
y = 0.04
print('SOR')
initialMesh = genMesh(horLine,verLine)
finalMesh = numIteration(initialMesh,horLine,verLine)
print ('Potential: ' + str(getPot(finalMesh,x,y,verLine,horLine)) + ' V')





    