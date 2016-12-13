#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Sharhad Bashar
260519664
ECSE 543
Assignment 2
conjGradMethod.py
Calculates the conjugate gradient for the cable
Nov 14th, 2016
""" 

import numpy as np
from basicDefinitions import matrixMult, matTranspose, cholesky, forwElim, backSub, matrixAddorSub, scalarMult
from q_3Functions import genMesh
import math

def genAb(h,freeNode):
    
    A = [[-4, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, -4, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
    [0, 1, -4, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
    [1, 0, 0, -4, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
    [0, 1, 0, 1, -4, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
    [0, 0, 1, 0, 1, -4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
    [0, 0, 0, 1, 0, 0, -4, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,],
    [0, 0, 0, 0, 1, 0, 1, -4, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,],
    [0, 0, 0, 0, 0, 1, 0, 1, -4, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 2, 0, 0, 0, 1, 0, 0, 0, 0,],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -4, 1, 0, 0, 0, 1, 0, 0, 0,],
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -4, 1, 0, 0, 0, 1, 0, 0,],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -4, 1, 0, 0, 0, 1, 0,],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -4, 0, 0, 0, 0, 1,],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -4, 2, 0, 0, 0,],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -4, 1, 0, 0,],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -4, 1, 0,],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -4, 1,],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -4,]]
    

    b = [-110,0,0,-110,0,0,-110,0,0,-110,-110,0,0,0,0,0,0,0,0]    
    cableHeight = 0.1
    cableWidth = 0.1
    coreHeight = 0.02
    coreWidth = 0.04
    corePot = 110.0
    mesh = genMesh(0.02)
    
    nodeHeight = (int)(cableHeight/h + 1)
    nodeWidth = (int)(cableWidth/h + 1)
    coreHeightNode = (int)(coreHeight/h + 1)
    codeWidthNode = (int)(coreWidth/h + 1)
    freeNode = (nodeHeight * nodeWidth) - (coreHeightNode * codeWidthNode) - (nodeHeight + nodeWidth) + 1
    
    A = [[0 for x in range (freeNode)] for y in range (freeNode)]
    b = [0 for x in range(freeNode)]    
         
    for y in range (nodeHeight):
        for x in range (nodeWidth):
            if (mesh[y][x][0] != None):
                cord = mesh[y][x][0]
                cordRight = mesh[y][x + 1][0]
                cordLeft = mesh[y][x - 1][0]
                cordUp = mesh[y + 1][x][0]
                cordDown = mesh[y - 1][x][0]
                A[cord][cord] = -4.0
                #print (cord, cordLeft, cordRight, cordUp, cordDown)
                #right
                if (cordRight != None):
                    A[cord][cordRight] = 1
                    if (cord == 0):
                        A[cord][cordRight] = 2
                        b[cord] = -corePot           
                #left
                if (cordLeft != None):
                    A[cord][cordLeft] = 1
                    if (cord == 1):
                        A[cord][cordLeft] = 2
                        b[cord + 1] = -corePot
                    if (cordLeft == 4 or cordLeft == 9 or cordLeft == 14):
                        A[cord][cordLeft] = 2
                #up
                if (cordUp != None):
                    A[cord][cordUp] = 1
                    if (cord == 4 or cord == 9):
                        A[cord][cordUp] = 2
                        if (cord == 4):
                            b[cord] = -corePot
                            b[cord + 1] = -corePot
                            b[cord + 2] = -corePot
                #down
                if (cordDown != None):
                    A[cord][cordDown] = 1
                    if (cord == 14 or cord == 9):
                        A[cord][cordDown] = 2 
                    if (cordDown == 0 or cordDown == 1):
                        A[cord][cordDown] = 2
    
    return (A,b)
##################################################
#Choleski
def Cholesky (A, b):
    AT = matTranspose(A)
    b = matTranspose(b)
    newA = matrixMult(AT,A)
    newb = matrixMult(AT,b)
    length = 19#len(newA)
    L = cholesky(newA,length)
    y = forwElim (L, newb, length)
    X = backSub (L, y, length) 
    print ('Choleski: ')
    print X
##################################################    
#Conjugate Gradient
def conjGrad(A,b,freeNode):
    AT = matTranspose(A)
    b = matTranspose(b)
    newA = matrixMult(AT,A)
    newb = matrixMult(AT,b)
    X = [0 for x in range (freeNode)] 
    X = matTranspose(X)     
    r = matrixAddorSub(newb, matrixMult(newA,X),'s')
    p = r
    
    for i in range(freeNode):
        alpha = float(matrixMult(matTranspose(p),r)[0][0])/float((matrixMult(matrixMult(matTranspose(p),newA),p))[0][0])
        X = matrixAddorSub(X, matTranspose(scalarMult(alpha,p)),'a')
        r = matrixAddorSub(newb, matrixMult(newA,X),'s')
        beta = (-1)*float(matrixMult((matrixMult(matTranspose(p),newA)),r)[0][0])/float((matrixMult(matrixMult(matTranspose(p),newA),p))[0][0])
        p = matrixAddorSub(r, matTranspose(scalarMult(beta,p)),'a')
        
        #finding the norms
        infNorm = 0
        twoNorm = 0
        for j in range (freeNode):
            value = abs(r[j][0])
            if (value > infNorm):
                infNorm = value
            twoNorm += r[j][0]**2
        twoNorm = math.sqrt(twoNorm)
        print (i,twoNorm,infNorm)
    print X
##################################################      
cableHeight = 0.1
cableWidth = 0.1
coreHeight = 0.02
coreWidth = 0.04
corePot = 110.0

h = 0.02    
nodeHeight = (int)(cableHeight/h + 1)
nodeWidth = (int)(cableWidth/h + 1)
coreHeightNode = (int)(coreHeight/h + 1)
codeWidthNode = (int)(coreWidth/h + 1)
freeNode = (nodeHeight * nodeWidth) - (coreHeightNode * codeWidthNode) - (nodeHeight + nodeWidth) + 1
 #Creates A and b   
(A,b) = genAb(0.02,freeNode)
#Cholesky(A,b)
conjGrad (A,b,freeNode)





 

