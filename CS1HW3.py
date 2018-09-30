#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 13:45:51 2018

@author: Stephen
"""
import matplotlib.pyplot as plt
import math
import numpy as np

def initialplots():
    arrayx = np.array((5,15,25,35,45,55,65,75,85,95,105,115))
    arrayy = np.array((32,17,21,7,8,6,5,3,4,1,5,1))
    logarrayy = np.zeros((len(arrayy)))
    for i in range(len(arrayy)):
        logarrayy[i] = math.log(arrayy[i])
    plt.plot(arrayx, arrayy, 'bo')
    plt.show()
    plt.plot(arrayx, logarrayy, 'ro')
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
    plt.Figure()
    plt.plot(approxlinex, approxliney, 'b')
    plt.plot(arrayx, logarrayy, 'ro')
    plt.show()
    
initialplots()





    