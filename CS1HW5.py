#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 16:05:03 2018

@author: Stephen
"""

import numpy as np
import matplotlib.pyplot as plt

def function(x, t):
    return (1 + t) * x + 1 - (3 * t) + t ** 2

def slope_field():
    Tt, Xx = (25, 25)
    trange = np.linspace(-3, 3, Tt)
    xrange = np.linspace(0, 5, Xx)
    grid = np.meshgrid(trange, xrange)
    fxt = function(grid[0] , grid[1])
#    print(fxt)
    plt.figure()
    plt.quiver(grid[0], grid[1], fxt, fxt)
    plt.show()
    
slope_field()


def simple_euler(xi, ti, deltat):
    xiplus1 = xi + function(xi, ti) * deltat
    return xiplus1

def simple_euler_routine(x0, t0, maxt, step):
    currentt = t0
    currentx = x0
    tvals = []
    xvals = []
    tvals.append(t0)
    xvals.append(x0)
    if t0 < maxt:
        while currentt < maxt:
            currentx = simple_euler(currentx, currentt, step)
            xvals.append(currentx)
            currentt += step
            tvals.append(currentt)
        return tvals, xvals
    elif t0 > maxt:
        while currentt > maxt:
            currentx = simple_euler(currentx, currentt, step)
            xvals.append(currentx)
            currentt += step
            tvals.append(currentt)
        return tvals, xvals
    else:
        raise ValueError('Invalid choice of t0, tmax')

def simple_euler_plotting():
    ttest1, xtest1 = simple_euler_routine(0.0655, 0, 2.5, 0.01)
    ttest2, xtest2 = simple_euler_routine(0.0663, 2.5, 0, -0.01)
    tsforcritx = np.linspace(0, 3, 20)
    critx = np.ones((len(tsforcritx)))
    critx = np.multiply(critx, 0.065923)
    plt.figure()
    plt.plot(ttest1, xtest1, 'r', label ='Approaching critical x from below')
    plt.plot(ttest2, xtest2, 'b', label ='Approaching critical x from above')
    plt.plot(tsforcritx, critx, 'g-', label ='Critical x')
    plt.legend(loc=2)
    plt.plot()
    
simple_euler_plotting()
    
    
    
    