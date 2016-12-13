# -*- coding: utf-8 -*-
"""
Sharhad Bashar
ECSE 543
Assignment 1
Oct 17th, 2016
question_3.py
"""

from q_3Functions import genMesh, numIteration, getPot
######################################################################
#Fixed value of h, w changing
x = 0.06
y = 0.04
h = 0.02

print('SOR')
for w in range (10,20,1):
    w = float(w/10)
    initialMesh = (genMesh(h))
    finalMesh = numIteration(initialMesh,h,w,'s')
    print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
    print('')
print('Jacobian')
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'j')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')
######################################################################
print('#######################################################')
#w = 1.25 is the value that gives least number of iterations
w = 1.25
h = 0.025
print('SOR')
print('h: ' + str(h))
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'s')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')
print('Jacobian')
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'j')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')

h = 0.02
print('SOR')
print('h: ' + str(h))
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'s')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')
print('Jacobian')
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'j')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')

h = 0.015
print('SOR')
print('h: ' + str(h))
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'s')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')
print('Jacobian')
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'j')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')

h = 0.01
print('SOR')
print('h: ' + str(h))
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'s')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')
print('Jacobian')
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'j')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')

h = 0.005
print('SOR')
print('h: ' + str(h))
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'s')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')
print('Jacobian')
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'j')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')

h = 0.002
print('SOR')
print('h: ' + str(h))
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'s')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')
print('Jacobian')
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'j')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')

h = 0.001
print('SOR')
print('h: ' + str(h))
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'s')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')
print('Jacobian')
initialMesh = (genMesh(h))
finalMesh = numIteration(initialMesh,h,w,'j')
print ('Potential: ' + str(getPot(finalMesh,x,y,h)) + ' V')
print('')








