# -*- coding: utf-8 -*-
"""
Sharhad Bashar
ECSE 543
Assignment 1
OCt 17th, 2016
question_1.py
"""
from basicDefinitions import cholesky, forwElim, backSub, readCell, randomSPD, voltageSolver
#######################################################################################################################
#Sets up the A,E,J,R Matricies for the five example circuits provided to us
incidentMatrix_1 = [readCell(1, 1),readCell(2, 1)]
E_1 = [[readCell(7, 1)],[readCell(7, 2)]]
J_1 = [[readCell(8, 1)],[readCell(8, 2)]]
R_1 = [[readCell(9, 1)],[readCell(9, 2)]]

incidentMatrix_2 = [readCell(1, 4),readCell(2, 4)]
E_2 = [[readCell(7, 4)],[readCell(7, 5)]]
J_2 = [[readCell(8, 4)],[readCell(8, 5)]]
R_2 = [[readCell(9, 4)],[readCell(9, 5)]]

incidentMatrix_3 = [readCell(1, 7),readCell(2, 7)]
E_3 = [[readCell(7, 7)],[readCell(7, 8)]]
J_3 = [[readCell(8, 7)],[readCell(8, 8)]]
R_3 = [[readCell(9, 7)],[readCell(9, 8)]]

incidentMatrix_4 = [[readCell(1,10),readCell(2,10),readCell(3,10),readCell(4,10)],
                    [readCell(1,11),readCell(2,11),readCell(3,11),readCell(4,11)]]
E_4 = [[readCell(7,10)],[readCell(7,11)],[readCell(7,12)],[readCell(7,13)]]
J_4 = [[readCell(8,10)],[readCell(8,11)],[readCell(8,12)],[readCell(8,13)]]
R_4 = [[readCell(9,10)],[readCell(9,11)],[readCell(9,12)],[readCell(9,13)]]

incidentMatrix_5 = [[readCell(1,15),readCell(2,15),readCell(3,15),readCell(4,15),readCell(5,15),readCell(6,15)],
                    [readCell(1,16),readCell(2,16),readCell(3,16),readCell(4,16),readCell(5,16),readCell(6,16)],
                    [readCell(1,17),readCell(2,17),readCell(3,17),readCell(4,17),readCell(5,17),readCell(6,17)],]
E_5 = [[readCell(7,15)],[readCell(7,16)],[readCell(7,17)],[readCell(7,18)],[readCell(7,19)],[readCell(7,20)]]
J_5 = [[readCell(8,15)],[readCell(8,16)],[readCell(8,17)],[readCell(8,18)],[readCell(8,19)],[readCell(8,20)]]
R_5 = [[readCell(9,15)],[readCell(9,16)],[readCell(9,17)],[readCell(9,18)],[readCell(9,19)],[readCell(9,20)]]
####################################################################################################################### 
#for a test run
#A = [[4,12,-16],[12,37,-43],[-16,-43,98]]
#b = [100,200,150]
#x = [2541.667, -683.334, 116.667]    
###################################################################################### 
#Main     
#A,b,L,y,x = 0,0,0,0,0 #initialize the values to zero
#lengthInput = int(input("Please enter the length of A: "))
#A = randomSPD(lengthInput) # creates a random SPD matrix of any requested length
#b = []
#for i in range (lengthInput):
#    count = i + 1
#    if (count == 1):
#        abbv = 'st'
#    elif (count == 2):
#        abbv = 'nd'
#    elif (count == 3):
#        abbv = 'rd'
#    else:
#        abbv = 'th'
#    bVal = input("Please enter the " + str(count) + abbv + " value of b: ")
#    b.append(float(bVal))
    
######################################################################################  
#function calls
chol = cholesky(A, lengthInput)
y = forwElim(chol, b, lengthInput)
x = backSub(chol,y,lengthInput)
######################################################################################
#prints the input A
print ("")
print('A: ')
for i in range (lengthInput):
    print(A[i])
print("")
#prints the input b
print("b:")
print(b)
print("")
#prints the output x
print ('x:')
print (x)
print ('')
#######################################################################################
#1_D 
print('Voltage for circuit 1 is: ' + str(voltageSolver(incidentMatrix_1,E_1,J_1,R_1)) + 'V')
print('Voltage for circuit 2 is: ' + str(voltageSolver(incidentMatrix_2,E_2,J_2,R_2)) + 'V')
print('Voltage for circuit 3 is: ' + str(voltageSolver(incidentMatrix_3,E_3,J_3,R_3)) + 'V')
print('Voltages for circuit 4 are: ' + str(voltageSolver(incidentMatrix_4,E_4,J_4,R_4)) + 'V')
print('Voltages for circuit 5 are: ' + str(voltageSolver(incidentMatrix_5,E_5,J_5,R_5)) + 'V')
#######################################################################################

