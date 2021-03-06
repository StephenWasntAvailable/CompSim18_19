# -*- coding: utf-8 -*-
"""
@author: Stephen
"""


import matplotlib.pyplot as plt
import math
import numpy as np

#Function to calculate the value x*exp(x)
def xexpx(x):
    return (x * math.exp(x))

#Analytic value of the second derivative of x*exp(x)
def ord2diffexact(x):
    return 2 * math.exp(x) + xexpx(x)

#Function to calculate the 2nd order central difference approximation of x*exp(x)
def ord2centdiffapprox1(xval,hval):
    return (  ( xexpx(xval + hval) - (2 * xexpx(xval)) + xexpx(xval - hval) ) / hval ** 2 )

#Function to create an array of x-values for evaluating functions
def createarrays(x):
    xvalues = np.linspace(0, 2, x)
    numpoints = len(xvalues)
    return xvalues, numpoints


#Recursive function to execute the Richardson Extrapolation procedure
#Arguments are the order of the extrapolation i.e D_i, x the value of x for the derivative to be evaluated at 
#and h the step size for performing the central difference approximation at
#For i == 1, just returns the standard central difference approximation
def richextrap(i, x, h):
    if i>1:
        dip1 = (((4 ** (i-1)) * richextrap(i-1, x, h) ) - richextrap(i-1, x, 2 * h)) / (4 ** (i-1) - 1)
        return dip1
    elif i==1:
        return ord2centdiffapprox1(x, h) 
    else:
        raise ValueError('Invalid value of i used')

#Function for stepping down the value of h slowly
def hval(x):
    return 0.4 / (1.0 + 0.01 * (x-1))

def hvalfast(x):
    return 0.4 / (x ** 2)

#Initialise arrays for storing the approximations and relative errors
def initarrays(kmax, lmax):
    array1 = np.zeros((kmax-2, lmax-1))
    array2 = np.zeros((kmax-2, lmax-1))
    return array1, array2
    
#function for finding approximation Di(h) when varying h evaluated at x
def varyh(k, lmax, x, hvalues):
    approxval = np.zeros(lmax - 1)
    relerror = np.zeros(lmax - 1)
    for l in range(1, lmax): #varying h
        approxval[l-1] = richextrap(k, x, hvalues[l-1])
        relerror[l-1] = abs((ord2diffexact(x) - approxval[l-1]) / ord2diffexact(x))
    return approxval, relerror

#Function that adds the approximations and errors into the arrays
def varyi(x, kmax, lmax):
    approxarray, errorarray = initarrays(kmax, lmax)
    hvalues = np.zeros(lmax - 1)
    for l in range(1, lmax):
        hvalues[l-1] = hval(l)
    for k in range(2, kmax):
        approxarray[k-2], errorarray[k-2] = varyh(k, lmax, x, hvalues)
    return hvalues, errorarray, approxarray

#Function that executes the evaluation of the approximations where we vary i and h
#x the value of x to evaluate at, kmax-2 gives the number of approximations orders to check, starting at D_2
#so kmax = 6 will check D_2 through D_5
#lmax the number of hvalues to check ranging from 0.4 to 0.4 / (1 + (0.01 * (lmax - 1) ) ) 
#Produces plots after calculation
def varyivaryh(x, kmax, lmax):
    hvalues, errorarray, approxarray = varyi(x, kmax, lmax)
    exactval = np.zeros(lmax - 1)
    for j in range(lmax-1):
        exactval[j] = ord2diffexact(x)
    plt.figure(1)
    plt.title('Approximations and Analytic Value at x=2 with varying h')
    plt.xlabel('Step size: h')
    plt.ylabel('f''(x)')
    plt.plot(hvalues, exactval, label = 'Analytical Value')
    for i in range(kmax-2):
        s = 'Approximation D$_%a $' % (i+2)
        plt.plot(hvalues, approxarray[i], label = s)
    plt.legend(loc = 3)
    plt.show()
    plt.figure(2)
    plt.title('Relative Error of Approximations at x=2 with varying h')
    plt.xlabel('Step size: h')
    plt.ylabel('Relative Error')
    for i in range(kmax-2):
        s = 'Approximation D$_%a $' % (i+2)
        plt.plot(hvalues, errorarray[i], label = s)
    plt.legend(loc = 2)
    plt.show()
    plt.figure(3)
    plt.title('LogLog Plot of Relative Errors')
    plt.xlabel('Step size: h')
    plt.ylabel('Relative Error')
#    plt.ylim(ymin = 10**(-13),ymax = 10**(-9))
    for i in range(kmax-2):
        s = 'Approximation D$_%a $' % (i+2)
        plt.loglog(hvalues, errorarray[i], label = s)
    plt.legend(loc = 4)
    plt.show()



#Function for checking for the effects of rounding errors as h gets smaller
#x is the x value to evaluate the approximation at
#k is the order of the Richardson Extrapolation to evaluate the approximation at
#lmax is the number of h values to check, from 0.4 to 0.4 / lmax ** 2 for illustrating how it changes at low h
#Also used this for finding best approximation, using hval() instead of hvalfast()
def rounding(x, k, lmax):
    hvalues = np.zeros(lmax - 1)
    for l in range(1, lmax):
        hvalues[l-1] = hval(l)
#    for l in range(1, lmax):
#        hvalues[l-1] = hvalfast(l)
    approxarray, errorarray = varyh(k, lmax, x, hvalues)
#    plt.figure(4)
#    plt.plot(hvalues, approxarray, 'x')
#    plt.ylim(ymin = 29.4, ymax = 29.6)
#    plt.show()
    plt.figure(5)
    plt.plot(hvalues, errorarray, 'x')
    plt.show()
    plt.figure(6)
    plt.title('LogLog plot for illustrating loss of accuracy at small h')
    plt.ylabel('Relative Error')
    plt.xlabel('Step size: h')
    s = 'D$_%d$' % k
    plt.loglog(hvalues, errorarray, label=s)
    plt.legend(loc = 2)
    plt.show()    
    bestapprox = np.argmin(errorarray)
    bestapproxvalue = approxarray[bestapprox]
    s1 = 'Best Approximation for Second Derivative at x=2 : ' + str(bestapproxvalue)
    s2 = 'Exact value of Second Derivative at x=2 : ' + str(ord2diffexact(x))
    s3 = 'Relative Error in Approximation : ' + str(errorarray[bestapprox])
    print(s1)
    print(s2)
    print(s3)
    


varyivaryh(2, 10, 400)
rounding(2, 5, 322)

print("Stephen O'Shea")
print('SN = 13321762')
    

    
