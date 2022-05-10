# -*- coding: utf-8 -*-
"""
Created on Tue May 10 22:03:13 2022

@author: suliki
"""

import sys
import numpy as np


def gauss_elimination(matrix):
    matrix = np.array(matrix, dtype=np.float64)
    matrix.dtype = np.float64
    for x in range(len(matrix)):
        if matrix[x][x] == 0.0:
            sys.exit("dzielenie przez 0")
        for y in range (x+1, len(matrix)):
            ratio = matrix[y][x]/matrix[x][x]
            for z in range(len(matrix)+1):
                matrix[y][z] = matrix[y][z] - ratio * matrix[x][z]
    return matrix
matrix = [[2,-2,5,0,34],[3,1,3,5,34],[3,2,3,2,17],[5,2,3,2,1]]
print(gauss_elimination(matrix))

print(matrix)
print(matrix[0][2])

