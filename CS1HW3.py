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

    
initialplots()

#estimate of equation of line of best fit
#  ~ -3.5 / 115 = -7/230
# y = mx + c : y = -7/230*x + c
# want y = 0 when x ~ 130
# c = 910/230 ~ 4
#so y = -7/230 * x + 4
