#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 13:45:51 2018

@author: Stephen
"""
import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.optimize as opt

#Function to plot the provided data, as well as the eye-estimated fit equation
def initialplots():
    arrayx = np.array((5,15,25,35,45,55,65,75,85,95,105,115))
    arrayy = np.array((32,17,21,7,8,6,5,3,4,1,5,1))
    logarrayy = np.zeros((len(arrayy)))
    for i in range(len(arrayy)):
        logarrayy[i] = math.log(arrayy[i])
    plt.plot(arrayx, arrayy, 'bo', label = 'Provided Data')
    plt.title('Standard plot of Decays against Time')
    plt.xlabel('Time')
    plt.ylabel('Decays')
    plt.legend(loc = 1)
    plt.show()
    plt.plot(arrayx, logarrayy, 'ro', label = 'Provided Data')
    plt.title('Plot of ln(Decays) against Time')
    plt.xlabel('Time')
    plt.ylabel('Natural Log of Decays')
    plt.legend(loc = 1)
    plt.show()
#estimate of equation of line of best fit
#  ~ -3.5 / 110 = -7/220
# y = mx + c : y = -7/220*x + c
# want y = 0 when x ~ 110
# c = 770/220 = 3.5
#so y = -7/220 * x + 3.5
    approxlinex = np.linspace(0, 120, 121)
    approxliney = np.zeros(len(approxlinex))
    for i in range(len(approxlinex)):
        approxliney[i] = ((-7/220) * approxlinex[i]) + 3.5
    plt.figure()
    plt.plot(approxlinex, approxliney, 'b', label = 'Approximate Fit')
    plt.plot(arrayx, logarrayy, 'ro', label = 'Data')
    plt.xlabel('Time')
    plt.ylabel('Natural Log of Decays')
    plt.title('Data with Eye-Approximated Fit Line')
    plt.legend(loc = 3)
    plt.show()
    
#initialplots()


#Function which holds the data to be fitted, calls the evaluation function for finding
#the optimal a and b values in a given range amin -> amax and bmin -> bmax, returns
#the best values of a, b and the associated least square sum, for comparison between different runs
#of the evaluation() function
#also plots the current fit against the data provided
def lssqfit(amin, amax, bmin, bmax, size):
    datax = np.array((5,15,25,35,45,55,65,75,85,95,105,115))
    datay = np.array((32,17,21,7,8,6,5,3,4,1,5,1))
    datalogy = np.zeros((len(datay)))
    svalues = []
    for i in range(len(datay)):
        datalogy[i] = math.log(datay[i])
    abminval = []
    abminval = evaluate(datax, datalogy, amin, amax, bmin, bmax, size, svalues)
#    print(abminval)
    a = abminval[0][0]
    b = abminval[0][1]
    smin = abminval[1]
#    print(a, b, smin)
    approxlinex = np.linspace(0, 120, 121)
    approxliney = np.zeros(len(approxlinex))
    for i in range(len(approxlinex)):
        approxliney[i] = a + b * approxlinex[i]
#    plt.figure()
#    plt.plot(approxlinex, approxliney, 'b', label = 'Approximate Fit')
#    plt.plot(datax, datalogy, 'ro', label = 'Data')
#    plt.xlabel('Time')
#    plt.ylabel('Natural Log of Decays')
#    plt.title('Data with Calculated Fit Line')
#    plt.legend(loc = 3)
#    plt.show()
    s1 = 'The value of the offset a is %f' % a
    s2 = 'The value of the slope b is %f' % b
    s3 = 'The value of the least squares sum is %f' % smin
    print(s1, s2, s3)
    return a,b,smin
    

#function for executing the least square fitting of the data, iterates over test values
# for a and b, and checks the sum of the squares of the differences from the dataset
#a and b values are kept until the completion of the function, and the best pair is returned
#after their index is found by finding the minimum squarre root of the sum of the squares
def evaluate(xaxis, yaxis, mina, maxa, minb, maxb, size, svalues):
    avalues = np.linspace(mina, maxa, size)
    bvalues = np.linspace(minb, maxb, size)
    fx = np.zeros((len(xaxis)))
    squarevalues = np.zeros((len(xaxis)))
    diff = np.zeros((len(xaxis)))
    abval = []
    for i in range(len(avalues)-1):
        for j in range(len(bvalues)-1):
            for k in range(len(xaxis)-1):
                fx[k] = avalues[i] + bvalues[j] * xaxis[k]
                diff[k] = abs(yaxis[k]- fx[k])
            squarevalues = np.square(diff)
            svalues.append(math.sqrt(np.sum(squarevalues)))
            abval.append([avalues[i],bvalues[j]])
    npsvalues = np.array(svalues)
    lssqfitarg = np.argmin(npsvalues)
    
    return abval[lssqfitarg], npsvalues[lssqfitarg]
  
#init function, idea is to gradually increase the accuracy over iterations of the fitting
#by decreasing the search space of a and b values to be centred around the most recently 
#found good value. New values for a and b are accepted if the sum of the squares returned
#by the most recent call is less than the previous best, otherwise, the function continues
#without updating a and b
#increasing maxrec and initsize parameters can greatly increase the running time          
def recursive(amininit, amaxinit, bmininit, bmaxinit, maxrec, initsize):
    a, b, smin = lssqfit(amininit, amaxinit, bmininit, bmaxinit, initsize)
    i = 1
    size = initsize
    for i in range(1, maxrec):
        j = float(i)
        anewmin = a - (1 / j ** 2) 
        anewmax = a + (1 / j ** 2) 
        bnewmin = b - (1 / j ** 2) 
        bnewmax = b + (1 / j ** 2) 
        newa, newb, newsmin = lssqfit(anewmin, anewmax, bnewmin, bnewmax, size)
        if newsmin < smin:
            a = newa; b = newb; smin = newsmin;
    return a, b
 
#Function which returns S for given values of a and b           
def lssqeval(a, b):
    datax = np.array((5,15,25,35,45,55,65,75,85,95,105,115))
    datay = np.array((32,17,21,7,8,6,5,3,4,1,5,1))
    datalogy = np.zeros((len(datay)))
    diff = np.zeros((len(datax)))
    fx = np.zeros((len(datax)))
    for i in range(len(datay)):
        datalogy[i] = math.log(datay[i])
    for k in range(len(datax)-1):
        fx[k] = a + b * datax[k]
        diff[k] = abs(datalogy[k]- fx[k])
    squarevalues = np.square(diff)
    svalues = (math.sqrt(np.sum(squarevalues)))  
    print(svalues)
        
        
#lssqfit(0, 5, -1, 0, 300)
#lssqeval(3.5, (-7/220))      
#recursive(0, 5, -1, 0, 30, 300)

def func(x, a, b):
    return np.exp(a * x + b)

#Function which uses existing Python library functions to perform a least squares fit
def existingfunctions():
    #inbuilt linear regression
    datax = np.array((5,15,25,35,45,55,65,75,85,95,105,115))
    datay = np.array((32,17,21,7,8,6,5,3,4,1,5,1))
    datalogy = np.zeros((len(datay)))
    for i in range(len(datay)):
        datalogy[i] = math.log(datay[i])
    numpypolyfit = np.polyfit(datax, datalogy, 1)
    print(numpypolyfit)
    lssqeval(numpypolyfit[1], numpypolyfit[0])
    
    #inbuilt non-linear best fit
    x0 = np.array([-7.0/220, 3.5])
    nonlinfit = opt.curve_fit(func, datax, datay, x0)
    fitparams = nonlinfit[0]
    print(fitparams)
    lssqeval(fitparams[1], fitparams[0])
    return numpypolyfit, fitparams

def fitplotting(amininit, amaxinit, bmininit, bmaxinit, maxrec, initsize):
    datax = np.array((5,15,25,35,45,55,65,75,85,95,105,115))
    datay = np.array((32,17,21,7,8,6,5,3,4,1,5,1))
    myleastsqfit = recursive(amininit, amaxinit, bmininit, bmaxinit, maxrec, initsize)
    numpypolyfit, scipypolyfit = existingfunctions()
    approxlinex = np.linspace(0, 120, 121)
    myapproxliney = np.zeros(len(approxlinex))
    numpyapproxliney = np.zeros(len(approxlinex))
    scipyapproxliney = np.zeros(len(approxlinex))
    for i in range(len(approxlinex)):
        myapproxliney[i] = np.exp(myleastsqfit[0] + myleastsqfit[1] * approxlinex[i])
    for i in range(len(approxlinex)):
        numpyapproxliney[i] = np.exp(numpypolyfit[1] + numpypolyfit[0] * approxlinex[i])
    for i in range(len(approxlinex)):
        scipyapproxliney[i] = np.exp(scipypolyfit[1] + scipypolyfit[0] * approxlinex[i])
    plt.figure()
    plt.plot(datax, datay, 'b.', label = 'Data')
    plt.plot(approxlinex, myapproxliney, 'r', label = 'My Functions Fit')
    plt.plot(approxlinex, numpyapproxliney, 'g', label = 'Numpy Fit')
    plt.plot(approxlinex, scipyapproxliney, 'k', label = 'Scipy Fit')
    plt.xlabel('Time')
    plt.ylabel('Decays')
    plt.title('Data with Calculated Fit Lines')
    plt.legend(loc = 'upper right')
    plt.show()

#Plots exponential curve based on eye fit parameters
def plotexpfit():
    datax = np.array((5,15,25,35,45,55,65,75,85,95,105,115))
    datay = np.array((32,17,21,7,8,6,5,3,4,1,5,1))
    approxlinex = np.linspace(0, 120, 121)
    approxliney = np.zeros(len(approxlinex))
    for i in range(len(approxlinex)):
        approxliney[i] = np.exp(3.5 + (-7.0/220) * approxlinex[i])
    plt.figure(0)
    plt.plot(datax, datay, 'b.', label = 'Data')
    plt.plot(approxlinex, approxliney, 'r', label = 'Estimated Fit')
    plt.xlabel('Time')
    plt.ylabel('Decays')
    plt.title('Data with Eye Estimated Fit')
    plt.legend(loc = 3)
    plt.show()

fitplotting(0, 5, -1, 0, 30, 300)
#plotexpfit()

    


print("Stephen O'Shea - JSPY - SN: 13321762")
        


    