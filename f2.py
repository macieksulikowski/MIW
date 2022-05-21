# -*- coding: utf-8 -*-
"""
Created on Sat May 21 08:24:42 2022

@author: suliki
"""

import numpy as np
import math

x = [1,1,1,1,1,1,1,1]
y = [1,1,1,1,-1,-1,-1,-1]
z = [1,1,-1,-1,0,0,0,0]
a = [0,0,0,0,1,1,-1,-1]
b = [1,-1,0,0,0,0,0,0]
c = [0,0,1,-1,0,0,0,0]
d = [0,0,0,0,1,-1,0,0]
e = [0,0,0,0,0,0,1,-1] 

macierz = np.array([x,y,z,a,b,c,d,e])

q = np.array([8,6,2,3,4,6,6,5])

# znormalizować wektory następnie dzielimy wszystkie elementy i wyiwetlamy

def len_v(v):
    return math.sqrt(np.dot(v,v))


baza = np.dot(macierz, macierz.T) # AA^T

print("baza:",baza)


ortogonalna = []
zm = 0
for row in macierz:
    zm = 0
    for i in row:
        zm += abs(i) # wartoć bezwzględna
    ortogonalna.append(zm)

print(ortogonalna) #macierz ortogonalna


ortonormalna = []

for i in range(len(macierz[1])):
    zm = macierz[i] * (1 / math.sqrt(ortogonalna[i]))
    ortonormalna.append(zm)
    
ortonormalna = np.array(ortonormalna) #macierz ortonormalna
print(ortonormalna)

print("test:",np.dot(ortonormalna,ortonormalna))

print(ortonormalna.T) # transponowana ortonormalna
print(np.linalg.inv(ortonormalna)) #funkcja odwracająca 
print(np.dot(ortonormalna,np.linalg.inv(ortonormalna))) # AA do -1

#zamiana baz

zamiana = np.dot(ortonormalna,q) #bo d już jest transponowane
print(zamiana)
