"""
This module contains the Pauli and Gell-Mann matrices. 
We put define the constants i and a for reading/writing convenience.

These objects are complex128
"""

import numpy as np

i = 1j
a = (1 / np.sqrt([3]))
p = [0,0,0,0]
l = [0,0,0,0,0,0,0,0,0]

#Pauli-Matrices
p[0] = np.matrix([[1, 0], [0, 1]], complex)
p[1] = np.matrix([[0, 1], [1, 0]], complex)
p[2] = np.matrix([[0,-i], [i, 0]], complex)
p[3] = np.matrix([[1, 0], [0,-1]], complex)

#Gell-Mann Matrices
l[0] = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]],   complex)
l[1] = np.matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]],   complex)
l[2] = np.matrix([[0,-i, 0], [i, 0, 0], [0, 0, 0]],   complex)
l[3] = np.matrix([[1, 0, 0], [0,-1, 0], [0, 0, 0]],   complex)
l[4] = np.matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]],   complex)
l[5] = np.matrix([[0, 0,-i], [0, 0, 0], [i, 0, 0]],   complex)
l[6] = np.matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]],   complex)
l[7] = np.matrix([[0, 0, 0], [0, 0,-i], [0, i, 0]],   complex)
l[8] = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0,-2]]*a, complex)  #a = 1/sqrt(3)

#Expand a matrix in a given basis
def expand(matrix, basis):
    a = []
    for k in range(len(basis)):
        a.append(np.ravel(basis[k]))
    A = (np.matrix(a)).getT()
    B = np.ravel(matrix)
    return np.linalg.solve(A, B) 

#Returns the structure constant f_abc for the algebra m
def strConst(m, a, b, c):
    com = m[a]*m[b] - m[b]*m[a]   #Commutator
    s = expand(com, m)
    return -i*s[c]/2              #Divide by 2 to get the GENERATORS s.c