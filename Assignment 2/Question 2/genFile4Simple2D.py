#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Sharhad Bashar
260519664
ECSE 543
Assignment 2
genFile4Simple2D.py
Creates the input file for the Simple2d_M program
"""

output = open("fileForSimple2D.txt", 'w') #opens the file
#All code below will go in here
#Also has print statements, so we can see in the IDE if proper numbers were printed
horNodes = 6
verNodes = 5
pot = 110.0
botBound = [(x+1) for x in range (horNodes - 2)]
mesh = [[((x + 5) + y * horNodes) for x in range(horNodes)] for y in range(verNodes)]
#########################################################################################
# First part of input
# Node -- xVal -- yVal
for i in range (horNodes - 2):
    print("{:<3} {:>6} {:>6}\n".format(botBound[i], (i + 2) * 0.02, i * 0.00))
    output.write("{:<3} {:>6} {:>6}\n".format(botBound[i], (i + 2) * 0.02, i * 0.00))
         
for y in range (verNodes):
    for x in range (horNodes):
        print("{:<3} {:>6} {:>6}\n".format(mesh[y][x], x * 0.02, (y + 1) * 0.02))
        output.write("{:<3} {:>6} {:>6}\n".format(mesh[y][x], x * 0.02, (y + 1) * 0.02))
#########################################################################################        
print("\n")
output.write("\n")
#########################################################################################      
#Second part of the input 
# Node1 -- Node2 -- Node3 -- 0.000
for i in range (horNodes - 3):
    print("{:<3} {:<3} {:<3} {}\n".format(botBound[i],botBound[i + 1],mesh[0][i + 2],'0.000'))
    print("{:<3} {:<3} {:<3} {}\n".format(botBound[i + 1],mesh[0][i + 3],mesh[0][i + 2],'0.000'))
    output.write("{:<3} {:<3} {:<3} {}\n".format(botBound[i],botBound[i + 1],mesh[0][i + 2],'0.000'))
    output.write("{:<3} {:<3} {:<3} {}\n".format(botBound[i + 1],mesh[0][i + 3],mesh[0][i + 2],'0.000'))

for y in range (verNodes - 1):
    for x in range (horNodes - 1): 
        print("{:<3} {:<3} {:<3} {:>6}\n".format(mesh[y][x],mesh[y][x + 1],mesh[y + 1][x],'0.000')) 
        print("{:<3} {:<3} {:<3} {:>6}\n".format(mesh[y][x + 1],mesh[y + 1][x + 1],mesh[y + 1][x],'0.000')) 
        output.write("{:<3} {:<3} {:<3} {:>6}\n".format(mesh[y][x],mesh[y][x + 1],mesh[y + 1][x],'0.000')) 
        output.write("{:<3} {:<3} {:<3} {:>6}\n".format(mesh[y][x + 1],mesh[y + 1][x + 1],mesh[y + 1][x],'0.000')) 
########################################################################################  
print("\n")
output.write("\n")
######################################################################################### 
#Third part of the input
#Boundary Conditions
#Boundary Node -- Boundary Conditions
#Dirchlet Condition:
print ("{:<3} {}\n".format(1,pot)) 
output.write ("{:<3} {}\n".format(1,pot))
for i in range (horNodes):
    if (mesh[0][i] < 8):
        print ("{:<3} {}\n".format(mesh[0][i],pot))
        output.write ("{:<3} {}\n".format(mesh[0][i],pot))
for i in range (horNodes):
    print ("{:<3} {}\n".format(mesh[verNodes-1][i],'0.000')) 
    output.write ("{:<3} {}\n".format(mesh[verNodes-1][i],'0.000')) 
for i in range (verNodes - 1):
    print ("{:<3} {}\n".format(mesh[i][horNodes - 1],'0.000'))
    output.write ("{:<3} {}\n".format(mesh[i][horNodes - 1],'0.000')) 
######################################################################################### 
output.close() #closes the file
         
