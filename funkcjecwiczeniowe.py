# -*- coding: utf-8 -*-
"""
Created on Sat Mar 5 20:30:19 2022

@author: suliki
"""
import sys
import math as m
import numpy as np
import random


def pobierz_plik():
    kontener = []
    with open("australian.dat", "r") as file:
        for line in file:
            kontener.append(list(map(lambda x: float(x), line.replace("\n", "").split())))
    return kontener

def metryka_euklidesowa(lista1, lista2):
    suma = 0
    for i in range(max(len(lista1), len(lista2))-1):
        suma = suma + (lista1[i]- lista2[i]) ** 2
    return np.sqrt(suma)

def mieszanie(zbior):
    zm = [i.copy() for i in zbior]
    for i in range(len(zbior)):
        zm[i][-1] = random.choice([0,1])
    return zm

def sort_by_key(zbior):
    kontener = {}
    for i in zbior:
        if i[-1] not in kontener.keys():
            kontener[i[-1]] = [i]
        else:
            kontener[i[-1]].append(i)
    return kontener

def srodek_masy(posortowany_zbior, n=False):
    kontener = {}
    for key in posortowany_zbior.keys():
        zm = {}
        for x in range(len(posortowany_zbior[key])):
            if not n:
                d = 0
                for y in range(len(posortowany_zbior[key])):
                    if x != y:
                        d += metryka_euklidesowa(posortowany_zbior[key][x], posortowany_zbior[key][y])
                zm[tuple(posortowany_zbior[key][x])] = d
            else:
                d = []
                for y in range(len(posortowany_zbior[key])):
                    if x != y:
                        d.append(metryka_euklidesowa(posortowany_zbior[key][x], posortowany_zbior[key][y]))
                zm[tuple(posortowany_zbior[key][x])] = sum(d[:n]) / n
        kontener[key] = list(list(zm.keys())[list(zm.values()).index(min(zm.values()))])
    return kontener

def decyzja(i, srodek_masy):
    zm = srodek_masy.copy()
    for key in srodek_masy.keys():
        zm[key] = metryka_euklidesowa(i, srodek_masy[key])
    kontener = [x for x,y in zm.items() if y == min(zm.values())]
    if len(kontener) > 1:
        return None
    return kontener[0]

def segregacja(zbior, zbior2):
    swaps = 1
    while swaps > 0:
        swaps = 0
        posortowany_zbior = sort_by_key(zbior)
        srodek_masy2 = srodek_masy(posortowany_zbior, 10)
        for x in range(len(zbior)):
            if decyzja(zbior[x], srodek_masy2) != zbior[x][-1]:
                zbior[x][-1] = decyzja(zbior[x], srodek_masy2)
                swaps += 1
                print(swaps)
        return zbior
    
def monte_carlo(n, xp, xk, yp, yk):
    p = (abs(yk - yp)*abs(xk - xp))
    u = lambda x: x
    pod = 0
    for r in range(n):
        x = random.uniform(0, 1)
        y = random.uniform(0, 1)
        if y < u(x):
            pod += 1
    integral = p * (pod/n)
    print("monte carlo =>",integral)
    
def metoda_kwadratowa(n, xp, xk):
    dx = ((xk - xp)/n)
    u = lambda x: x**2
    res = 0
    for i in range(n):
        res += u(xp + i * dx) * dx
    print("square =>",res)

def skalar(v1,v2):
    v3 = 0
    for i in range(len(v1)):
        v3 += v1[i] * v2[i]
    return v3
    
def srednia(v1):
    l = len(v1)
    return skalar(v1, np.ones(l))/l

def wariancja(v1):
    l = len(v1)
    sr = srednia(v1)
    x = v1 -(np.ones(l) * sr)
    return skalar(x, x)/l

def odchylenie(v1):
    return m.sqrt(wariancja(v1))

def dlugosc_wektora(v):
    return pow(np.dot(v,v), 1/2)

def wskaznik(macierz):
    return np.linalg.det(macierz)

def odwrotnosc(macierz):
    return np.linalg.inv(macierz)

def projekcja(u, v):
    return (np.dot(v,u) / np.dot(u,u)) * u

def qr(macierz):
    macierz = np.transpose(macierz)
    q, u_list = [], [macierz[0]]
    e = np.array(u_list[0]) / dlugosc_wektora(u_list[0])
    q.append(e)
    for x in range(1, len(macierz)):
        u = macierz[x]
        for y in range(len(u_list)):
            u = u - projekcja(u_list[y], macierz[x])
        u_list.append(u)
        e = u / dlugosc_wektora(u)
        q.append(e)
    r = np.dot(q, np.transpose(macierz))
    q = np.transpose(q)
    return (q, r)

def ak(macierz, k):
    for x in range(len(macierz)):
        if (macierz == np.transpose(macierz)).all():
            for x in range(k):
                q, r = qr(macierz)
                macierz = np.dot(r,q)
            return macierz
        else:
            print("Macierz nie jest symetryczna")
            return

def czy_trojkatna(macierz, s = 0.001):
    for x in range(len(macierz)):
        for y in range(len(macierz[x])):
            if y < x:
                if macierz[x][y] > s:
                    return False
    return True

def wartosci_wlasne(macierz):
    for x in range(len(macierz)):
        if (macierz == np.transpose(macierz)).all():
            while not czy_trojkatna(macierz):
                q, r = qr(macierz)
                macierz = np.dot(r,q)
            return macierz
        else:
            print("Macierz niesymetryczna")
            return

def eliminacja_gaussa(macierz):
    macierz = np.array(macierz, dtype=np.float64)
    macierz.dtype = np.float64
    for x in range(len(macierz)):
        if macierz[x][x] == 0.0:
            sys.exit("dzielenie przez 0")
        for y in range (x+1, len(macierz)):
            zm = macierz[y][x]/macierz[x][x]
            for z in range(len(macierz)+1):
                macierz[y][z] = macierz[y][z] - zm * macierz[x][z]
    return macierz



# zbior = pobierz_plik()
# zbior_zmieszany = mieszanie(zbior)
# posortowany_zbior = sort_by_key(zbior)
# srodek_masy2 = srodek_masy(posortowany_zbior)
# zbior_koncowy = segregacja(zbior_zmieszany, zbior)
# monte_carlo(1000,0,1,0,1)
# metoda_kwadratowa(1000, 0, 1)
# print(zbior_koncowy)
# ar = [[1,2],[0,2]]
# ar2 = [[2,-2,5,0,34],[3,1,3,5,34],[3,2,3,2,17],[5,2,3,2,1]]
# print(qr(ar))
# print(czy_trojkatna(ar))
# print(wartosci_wlasne(ar))
# print(eliminacja_gaussa(ar2))
# print(wskaznik(ar))
# print(odwrotnosc(ar))







