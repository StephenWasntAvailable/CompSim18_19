# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 12:25:11 2018

@author: Stephen
"""


import matplotlib.pyplot as plt
import math
import numpy as np

#Function to calculate the value x*exp(x)
def xexpx(x):
    return (x * math.exp(x))

def ord2diffexact(x):
    return 2 * math.exp(x) + xexpx(x)

#Function to calculate the 2nd order central difference approximation of x*exp(x)
def ord2centdiffapprox1(xval,hval):
    return (  ( xexpx(xval + hval) - (2 * xexpx(xval)) + xexpx(xval - hval) ) / hval ** 2 )

def createarrays(x):
    xvalues = np.linspace(0, 2, x)
    numpoints = len(xvalues)
    return xvalues, numpoints

def richextrap(i, x, h):
    if i>1:
        dip1 = (((4 ** (i-1)) * richextrap(i-1, x, h) ) - richextrap(i-1, x, 2 * h)) / (4 ** (i-1) - 1)
        return dip1
    elif i==1:
        return ( 4 * ord2centdiffapprox1(x, h) - ord2centdiffapprox1(x, 2 * h)) / 3
    else:
        raise ValueError('Invalid value of i used')



def hval(x):
    return 0.8 / (2 * x)

def initarrays(kmax, lmax):
    array1 = np.zeros((kmax, lmax-1))
    array2 = np.zeros((kmax, lmax-1))
    return array1, array2
    
#function for finding approximation Di(h) when varying h evaluated at x
def varyh(k, lmax, x, hvalues):
    approxval = np.zeros(lmax - 1)
    relerror = np.zeros(lmax - 1)
    for l in range(1, lmax): #varying h
        approxval[l-1] = richextrap(k, x, hvalues[l-1])
        relerror[l-1] = abs((ord2diffexact(x) - approxval[l-1]) / ord2diffexact(x))
    return approxval, relerror

def varyi(x, kmax, lmax):
    approxarray, errorarray = initarrays(kmax, lmax)
    hvalues = np.zeros(lmax - 1)
    for l in range(1, lmax):
        hvalues[l-1] = hval(l)
    for k in range(1, kmax):
        approxarray[k-1], errorarray[k-1] = varyh(k, lmax, x, hvalues)
    print(approxarray)
    print(errorarray)

varyi(2, 5, 20)
#some initial testing of the functions
#y = richextrap(2, 2, 0.0001)
#appr1 = ord2centdiffapprox1(2,0.000001)
#appr2 = ord2centdiffapprox1(2,0.000002)
#ex = ord2diffexact(2)
#print(y)
#print(appr1)
#print(appr2)
#print(ex)
#bigarray = initarrays(5, 10)
#hvalues = np.zeros(10)
#for i in range(1,10):
#    hvalues[i-1] = hval(i)
#    
#bigarray[1] = hvalues
#print(bigarray)
    
