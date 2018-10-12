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
    intapprox = []
    for i in range(2, nmax):
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
    
    baseapprox, mac1approx, mac2approx = traparrayimp(a, b, nmax)
    relerrorbase = relerrorcalc(baseapprox, a, b)
    relerrormac1 = relerrorcalc(mac1approx, a, b)
    relerrormac2 = relerrorcalc(mac2approx, a, b)
    
    plt.figure()
    plt.loglog(nvalues, relerrorbase, 'bo')
    plt.loglog(nvalues, relerrormac1, 'ro')
    plt.loglog(nvalues, relerrormac2, 'go')
#    plt.ylim(ymin = 1e-1, ymax = 1e0)
#    plt.xlim(xmax = 1e1)
    plt.show()
    
    
def traparrayimp(a, b, nmax):
    intapprox = []
    mac1approx = []
    mac2approx = []
    for i in range(2, nmax):
        N = 2 ** i
        interval, weights = intervals(a, b, N)
        initialapprox = trapezoid(interval, weights)
        #add the first maclaurin expansion term to the approximation = (h^2 / 12) * (f'x1 - f'xN)
        #where f' is the derivative of f, will use the analytic values here i.e (sin x)' = cos x
        #and cos 0 = 1, cos pi = -1
        macapprox1 = initialapprox + ( ( (calch(a, b, N) ** 2) / 12 ) * 2 ) 
        #same again but for the second maclaurin expansion term as well
        #(sin x)''' - cos x
        macapprox2 = macapprox1 - ( ( (calch(a, b, N) ** 4) / 720 ) * (-2))
        intapprox.append(initialapprox)
        mac1approx.append(macapprox1)
        mac2approx.append(macapprox2)
    return intapprox, mac1approx, mac2approx

    
    
main(0, math.pi, 25)


print("Stephen O'Shea - SN = 13321762")
    