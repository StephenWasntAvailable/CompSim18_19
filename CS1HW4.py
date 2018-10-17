# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 11:38:17 2018

@author: Stephen
"""

import matplotlib.pyplot as plt
import math
import numpy as np
import timeit


#Function to be integrated
def fx(x):
    return np.sin(x)

#Analytic result of the integral on the region a,b
def intanalytical(a, b):
    return (np.cos(a) - np.cos(b))

#Function to find the evaluation points and calculate the weights
def intervals(a, b, N):
    interval = np.linspace(a, b, N)
    weights = [calch(a, b, N)] * N
    weights[0] = calch(a, b, N) / 2
    weights[len(interval)-1] = calch(a, b, N) / 2
    return interval, weights    

#Returns the value of h, used in the weights
def calch(a, b, N):
    return (b - a) / N

#Function to evaluate the numerical integration using trapezoid rule
#Commented out portion is the earlier implementation, which I changed to improve performance
#as it was causing both runtime and memory issues for very large values of N
def trapezoid(interval, weights):
    approx = 0
    interval = np.array(interval)
    weights = np.array(weights)
    fxinterval = np.sin(interval)
    approx = np.sum(np.multiply(fxinterval, weights))
#    for i in range(0, len(interval)-1):
#        approx += fx(interval[i]) * weights[i]
    return approx

#Function to iterate through several value of N and return the approximation  
def traparray(a, b, nmax):
    intapprox = []
    for i in range(2, nmax):
        N = 2 ** i
        interval, weights = intervals(a, b, N)
        intapprox.append(trapezoid(interval, weights))
    print(intapprox)
    return intapprox

#Function to find the relative error in an approximation, or the absolute error
# if the analytic value of the integral is 0
def relerrorcalc(approximations, a, b):
    analyticvalue = intanalytical(a, b)
    relerror = np.zeros(len(approximations))
    abserror = np.zeros(len(approximations))
    if analyticvalue == 0:
        print('Analytic Value for this integral is 0, relative error is undefined, absolute error will be used instead.')
        for i in range(len(relerror)):
            abserror[i] = abs(analyticvalue - approximations[i])
        return abserror
    else:  
        for i in range(len(relerror)):
            relerror[i] = abs((analyticvalue - approximations[i]) / analyticvalue)
        return relerror

    
#traparray(0, math.pi)


#Produces plots and fits to the data produced
def main(a, b, nmax):
    approximations = traparray(a, b, nmax)
    relerrors = relerrorcalc(approximations, a, b)
    nvalues = np.zeros(len(relerrors))
    for i in range(2, nmax):
        nvalues[i-2] = 2 ** i
    lognvalues = np.log(nvalues)
    logerrors = np.log(relerrors)
    fit = np.polyfit(lognvalues, logerrors, 1)
    print(fit)
    s = 'Slope = %f' %fit[0]
    plt.figure()
    plt.title('LogLog Plot of Relative Errors against Number of Intervals')
    plt.loglog(nvalues, relerrors, 'bo', label = s)
    plt.legend(loc = 1)
    plt.xlabel('N: Number of Intervals')
    plt.ylabel('Relative Error in Approximation')
    plt.show()
    
    baseapprox, mac1approx, mac2approx = traparrayimp(a, b, nmax)
    relerrorbase = relerrorcalc(baseapprox, a, b)
    relerrormac1 = relerrorcalc(mac1approx, a, b)
    relerrormac2 = relerrorcalc(mac2approx, a, b)
    
    plt.figure()
    plt.title('LogLog Plot of Relative Errors against Number of Intervals')
    plt.loglog(nvalues, relerrorbase, 'b.', label = 'Base Approximation')
    plt.loglog(nvalues, relerrormac1, 'r.', label = 'With First Maclaurin Correction Term')
    plt.loglog(nvalues, relerrormac2, 'g.', label = 'With First Two Maclaurin Correction Terms')
    plt.legend(loc = 1)
    plt.xlabel('N: Number of Intervals')
    plt.ylabel('Relative Error in Approximation')
    #use below for zooming on first approximation with N = 4
#    plt.ylim(ymin = 2.5e-1, ymax = 3.5e-1)
#    plt.yticks(np.arange(2.5e-1, 3.5e-1, 1e-2))
#    plt.xlim(xmax = 1e1)
    #use below for zooming in on last approximation with N = 1024
    plt.ylim(ymin = 0.976e-3, ymax = 0.979e-3)
    plt.yticks(np.arange(0.976e-3, 0.979e-3, 2e-7))
    plt.xlim(xmin = 1e2, xmax = 1e4)
    plt.show()
    logbaseerror = np.log(relerrorbase)
    logmac1error = np.log(relerrormac1)
    logmac2error = np.log(relerrormac2)
    lognvalues = np.log(nvalues)
#    print(logbaseerror, logmac1error, logmac2error)
    basefit = np.polyfit(lognvalues, logbaseerror, 1)
    mac1fit = np.polyfit(lognvalues, logmac1error, 1)
    mac2fit = np.polyfit(lognvalues, logmac2error, 1)
    print(basefit, mac1fit, mac2fit)
    
#Function for finding the approximations given the inclusion of the Maclaurin series correction terms
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
        #(sin x)''' =  - cos x
        macapprox2 = macapprox1 - ( ( (calch(a, b, N) ** 4) / 720 ) * (-2))
        intapprox.append(initialapprox)
        mac1approx.append(macapprox1)
        mac2approx.append(macapprox2)
    return intapprox, mac1approx, mac2approx

    
    
main(0, math.pi, 11)


#Equivalent to main() for executing the checking for rounding errors
def roundingcheck(a, b, nmax):
    approximations = traparraylargen(a, b, nmax)
    relerrors = relerrorcalc(approximations, a, b)
    nvalues = np.zeros(len(relerrors))
    for i in range(2, nmax):
        nvalues[i-2] = 100000 + i
    lognvalues = np.log(nvalues)
    logerrors = np.log(relerrors)
    fit = np.polyfit(lognvalues, logerrors, 1)
    print(fit)
    s = 'Slope = %f' %fit[0]
    plt.figure()
    plt.title('LogLog Plot for finding Rounding Errors')
    plt.loglog(nvalues, relerrors, 'b', label = s)
    plt.legend(loc = 1)
    plt.xlabel('N: Number of Intervals')
    plt.ylabel('Relative Error in Approximation')
    plt.show()


#equivalent to traparray(), just created as a separate function for testing for rounding issues  

def traparraylargen(a, b, nmax):
    intapprox = []
    for i in range(2, nmax):
        N = 100000 +  i
        interval, weights = intervals(a, b, N)
        intapprox.append(trapezoid(interval, weights))
    return intapprox

#Be careful with running this, as for larger N this will use significant amounts of RAM
#roundingcheck(0, math.pi, 20)

print("Stephen O'Shea - SN = 13321762")
    