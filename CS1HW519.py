# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 15:13:06 2019

@author: soshe
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy
import random
import time

def factorial(n):
    if (n == 0):
        return 1
    else:
        return n * factorial(n-1)
    
def poisson(avnval, n):
    return (avnval ** n) * np.exp(-avnval) / factorial(n)

def plot_poisson(avnval):
    nvalues = np.linspace(0, 50, 51)
    poissonarray = np.zeros(len(nvalues))
    for i in range(len(nvalues)-1):
        poissonarray[i] = poisson(avnval, nvalues[i])
    plt.figure()
    plt.plot(nvalues, poissonarray, 'r')
    plt.show()
    
def compute_sum_poisson(avnval):
    
        
plot_poisson(50)