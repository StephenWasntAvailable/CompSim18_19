# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 11:38:17 2018

@author: Stephen
"""

import matplotlib.pyplot as plt
import math
import numpy as np

def fx(x):
    return math.sin(x)

def intervals(a, b, N):
    interval = np.linspace(a, b, N)
    weights = np.zeros(len(interval))
    weights[0] = calch(a, b, N) / 2
    weights[len(interval)-1] = calch(a, b, N) / 2
    for i in range(1, len(interval)-2):
        weights[i] = calch(a, b, N)
    return interval, weights    

def calch(a, b, N):
    return (b - a) / N

def trapezoid(interval, weights):
    approx = 0
    for i in range(0, len(interval)-2):
        approx += fx(interval[i]) * weights[i]
    return approx
    
def traparray(a, b):
    N = 4;
    interval, weights = intervals(a, b, N)
    intapprox = []
    intapprox.append(trapezoid(interval, weights))
    for i in range(3, 10):
        N = 2 ** i
        interval, weights = intervals(a, b, N)
        intapprox.append(trapezoid(interval, weights))
    return intapprox

def errorcalc(approximations):
    return approximations
    
traparray(0, math.pi)