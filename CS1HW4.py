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

def intanalytical(a, b):
    return (math.cos(a) - math.cos(b))

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
    
def traparray(a, b, nmax):
    N = 4;
    interval, weights = intervals(a, b, N)
    intapprox = []
    intapprox.append(trapezoid(interval, weights))
    for i in range(3, nmax):
        N = 2 ** i
        interval, weights = intervals(a, b, N)
        intapprox.append(trapezoid(interval, weights))
    return intapprox

def relerrorcalc(approximations, a, b):
    analyticvalue = intanalytical(a, b)
    relerror = np.zeros(len(approximations))
    abserror = np.zeros(len(approximations))
    if analyticvalue == 0:
        print('Analytic Value for this integral is 0, relative error is undefined, absolute error will be used instead.')
        for i in range(len(relerror)-1):
            abserror[i] = abs(analyticvalue - approximations[i])
        return abserror
    else:  
        for i in range(len(relerror)-1):
            relerror[i] = abs((analyticvalue - approximations[i]) / analyticvalue)
        return relerror

    
#traparray(0, math.pi)

def main(a, b, nmax):
    approximations = traparray(a, b, nmax)
    relerrors = relerrorcalc(approximations, a, b)
    nvalues = np.zeros(len(relerrors))
    for i in range(2, nmax):
        nvalues[i-2] = 2 ** i
    plt.loglog(nvalues, relerrors, 'bo')
    plt.show()

    
    
main(0, 2 * math.pi, 15)
    