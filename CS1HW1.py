#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Stephen
"""

import matplotlib.pyplot as plt
import math
import numpy as np

#Function, 2nd deriv. of which to be approximated, in this case cos(x)
def cosx(x):
    return (math.cos(x))

#Function to gradually decrease the step value h 
#for use in the central difference approximation
def hsizeslow(x):
    h = math.pi / (10.0 * x)
    return h

#Function to decrease the step value h faster
def hsizefast(x):
    h = math.pi / (x ** x)
    return h

#Function to calculate the 2nd order central difference approximation of cos(x)
def ord2centdiffapprox1(xval,hval):
    return (  ( cosx((xval + hval)) - (2 * cosx(xval)) + cosx((xval - hval)) ) / hval ** 2 )

#Function, algebraically the same as the other function for approximating 
#the second derivative, that might result in slightly different results due 
#to order of operations
def ord2centdiffapprox2(xval, hval):
    a1 = cosx((xval + hval)) - cosx(xval)
    a2 = cosx((xval - hval)) - cosx(xval)
    a3 = a1 + a2
    return a3 / (hval ** 2)

#Analytical value of the 2nd derivative of cos(X) == -cos(x)
def ord2diffcosx(x):
    return -1 * math.cos(x)

#Creating an array of the x values the approximation will be evaluated at,
#in the range -4pi to 4pi, i.e 4 periods of cos(x)
def createarrays(x):
    xvalues = np.linspace(0, 4 * math.pi, x)
    numpoints = len(xvalues)
    return xvalues, numpoints
    
#Evaluating both approximation methods, the analytical value, and the relative
#errors of both methods, given input of x values
def calcderivanderror(arrayx, arraysize, h):
    approxvalues1 = np.zeros(arraysize)
    approxvalues2 = np.zeros(arraysize)
    actualvalues = np.zeros(arraysize)
    relativeerror1 = np.zeros(arraysize)
    relativeerror2 = np.zeros(arraysize)
    errordiff = np.zeros(arraysize)
    for i in range(arraysize):
        approxvalues1[i] = ord2centdiffapprox1(arrayx[i], h)
        approxvalues2[i] = ord2centdiffapprox2(arrayx[i], h)
        actualvalues[i] = ord2diffcosx(arrayx[i])
        relativeerror1[i] = abs( ( actualvalues[i] - approxvalues1[i] ) / actualvalues[i] )
        relativeerror2[i] = abs( ( actualvalues[i] - approxvalues2[i] ) / actualvalues[i] )
        errordiff[i] = abs(relativeerror1[i]-relativeerror2[i])
    return approxvalues1, approxvalues2, actualvalues, relativeerror1, relativeerror2, errordiff

#Function that creates some reused plots
def StandardPlots(x, y, z):
    plt.figure(1)
    plt.title('Approximation of 2nd Derivative')
    plt.plot(x, y, 'r')
    plt.xlabel('X')
    plt.ylabel('[cos(x)]'' Central Difference Approximation')
    plt.show()
    plt.figure(2)
    plt.title('Relative Error compared to Analytical Result')
    plt.plot(x, z, 'b+')
    plt.ylim(ymin = 0, ymax = (1.5 * np.amax(z)))
    plt.xlabel('X')
    plt.ylabel('Relative Error of Approximation')
    plt.show()

#Main function for producing plots, which calls the previously defined functions
#Loops through j = 1 to maxj which steadily decreases the step size h for the 
#approximation
def main(x, maxj):
    averror1 = np.zeros(maxj - 1)
    averror2 = np.zeros(maxj - 1)
    hvalues = np.zeros(maxj - 1)
    errordiff = np.zeros(x)
    xaxis, numpoints = createarrays(x)
    for j in range(1, maxj):
        hvalue = hsizeslow(j) 
        yapprox1, yapprox2, yactual, relerror1, relerror2, errordiff = calcderivanderror(xaxis, numpoints, hvalue)
        averror1[j-1] = np.average(relerror1)
        averror2[j-1] = np.average(relerror2)
        hvalues[j-1] = hvalue
        #Plots for given values of j
        if (j == 1):
            s = str(hvalue)
            StandardPlots(xaxis, yapprox1, relerror1)
            print('Value of h: ' + s)
        elif (j % 10000 == 0):
            s = str(hvalue)
            StandardPlots(xaxis, yapprox1, relerror1)
            print('Value of h: ' + s)
    minerrortemp = averror1[1]
    smoothuptemp = abs((averror1[2] - averror1[1]) / averror1[1])
    smoothdowntemp = abs((averror1[1] - averror1[0]) / averror1[1])
    indextemp = 1
    for k in range(2, maxj-2):
        temp1 = abs((averror1[k+1] - averror1[k]) / averror1[k])
        temp2 = abs(averror1[k] - averror1[k-1] / averror1[k])
        if((temp1 <= smoothuptemp and temp2 <= smoothdowntemp) and averror1[k] < minerrortemp):
            minerrortemp = averror1[k]
            smoothuptemp = temp1
            smoothdowntemp = temp2
            indextemp = k
        
    hindexminerror = indextemp
    minerror = str(minerrortemp)
    hminerror = str(hsizeslow(hindexminerror + 1))
    s1 = 'H which produced the lowest* average error: ' + hminerror
    s2 = 'Minimum* relative average error produced: ' + minerror
    print(s1)
    print(s2)
            

    plt.figure(5)
    plt.title('Standard Plot: Av. Error against h')
    plt.plot(hvalues, averror1, 'r+')
    plt.xlabel('Step value: h')
    plt.ylabel('Average error in Approximation')
    plt.grid()
    plt.show()
    plt.figure(6)
    plt.title('LogLog plot: Av. Error against h')
    plt.loglog(hvalues, averror1)
    plt.xlabel('Step value: h')
    plt.ylabel('Average error in Approximation')
    plt.grid()
    plt.show()
    plt.figure(7)
    plt.title('LogLog plot: Both Approximation methods')
    plt.loglog(hvalues, averror1, label ='Method 1')
    plt.loglog(hvalues, averror2, label ='Method 2')
    plt.xlabel('Step Size: h')
    plt.ylabel('Approximation of 2nd Derivative')
    plt.legend(loc = 4)
    plt.grid()
    plt.show()
    plt.figure(8)
    plt.title('Zoomed in graph to see difference between approximation methods')
    plt.loglog(hvalues, averror1, 'x', label ='Method 1')
    plt.loglog(hvalues, averror2, '+', label ='Method 2')
    plt.legend(loc = 4)
    plt.xlim(xmin = 3e-3, xmax = 4e-3)
    plt.figure(9)
    plt.plot(xaxis, errordiff, 'ro')
    plt.show()
    #plt.plot(arrayx, approxvalues, 'ro')
    #plt.show()
    #plt.plot(arrayx, relativeerror, 'bo')
    #plt.ylim(ymin = 0, ymax = (1.5 * np.amax(relativeerror)))
    #plt.show()
    
    
#Checking for subtraction cancellation effect where x -> x' : f(x) - > 0
def subcancelchecklower():
    relerrorl1 = np.zeros(20)
    relerrorl2 = np.zeros(20)
    hvalues = np.zeros(20)
    for k in range(1, 21):
        h = 1 / (k * 1000)
        hvalues[k-1] = h
        relerrorl1[k-1] = abs( (ord2diffcosx(1.57079616972) - ord2centdiffapprox1(1.57079616972, h) )/ ord2diffcosx(1.57079616972))
        relerrorl2[k-1] = abs( (ord2diffcosx(1.57079616972) - ord2centdiffapprox2(1.57079616972, h) )/ ord2diffcosx(1.57079616972))
    plt.figure(9)
    plt.title('Relative Error of different approximation constructions')
    plt.plot(hvalues, relerrorl1, 'r+', label ='Method 1')
    plt.plot(hvalues, relerrorl2, 'bx', label ='Method 2')
    plt.xlabel('Value of h: step size')
    plt.ylabel('Relative Error')
    plt.legend(loc=1)
    plt.show()

    
#Checking for subtraction cancellation effect where x -> x' : f(x) - > 0
def subcancelcheckupper():
    relerroru1 = np.zeros(20)
    relerroru2 = np.zeros(20)
    hvalues = np.zeros(20)
    for k in range(1, 21):
        h = 1 / (k * 1000)
        hvalues[k-1] = h
        relerroru1[k-1] = abs( (ord2diffcosx(1.57079648387) - ord2centdiffapprox1(1.57079648387, h) )/ ord2diffcosx(1.57079648387))
        relerroru2[k-1] = abs( (ord2diffcosx(1.57079648387) - ord2centdiffapprox2(1.57079648387, h) )/ ord2diffcosx(1.57079648387))
    plt.figure(10)
    plt.title('Relative Error of different approximation constructions')
    plt.plot(hvalues, relerroru1, 'r+', label ='Method 1')
    plt.plot(hvalues, relerroru2, 'bx', label ='Method 2')
    plt.xlabel('Value of h: step size')
    plt.ylabel('Relative Error')
    plt.legend(loc=1)
    plt.show()

#Using accelerated h decrease to show at what point the approximation breaks
#down due to numerical precision
def machprectest(x, maxj):
    #averror1 = np.zeros(maxj - 1)
    #averror2 = np.zeros(maxj - 1)
    #hvalues = np.zeros(maxj - 1)
    errordiff = np.zeros(x)
    xaxis, numpoints = createarrays(x)
    for j in range(1, maxj):
        hvalue = hsizefast(j) 
        yapprox1, yapprox2, yactual, relerror1, relerror2, errordiff = calcderivanderror(xaxis, numpoints, hvalue)
        if (j % 3 == 0):
            StandardPlots(xaxis, yapprox1, relerror1)

#Calling the functions            
#main(200, 10000)      #
subcancelchecklower()
subcancelcheckupper()
#machprectest(100, 10)


    
    

    
    