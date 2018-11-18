#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 16:05:03 2018

@author: Stephen
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.integrate import odeint

def function(x, t):
    return (1 + t) * x + 1 - (3 * t) + t ** 2

def slope_field():
    Tt, Xx = (50, 50)
    trange = np.linspace(-3, 3, Tt)
    xrange = np.linspace(0, 5, Xx)
    X, T = np.meshgrid(xrange, trange)
    fxt = function(X , T)
#    print(fxt)
    plt.figure()
    plt.quiver(T, X, fxt, fxt, units = 'width')
    plt.show()
    plt.figure()
    plt.pcolormesh(X, T, fxt, cmap=cm.PuOr, vmax=0.25*fxt.max())
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
    ttest1, xtest1 = simple_euler_routine(0.0655, 0, 3.0, 0.02)
    ttest2, xtest2 = simple_euler_routine(0.0663, 3.0, 0, -0.02)
    tsforcritx = np.linspace(0, 3, 20)
    critx = np.ones((len(tsforcritx)))
    critx = np.multiply(critx, 0.065923)
    plt.figure()
    plt.plot(ttest1, xtest1, 'r', label ='Approaching critical x from below')
    plt.plot(ttest2, xtest2, 'b', label ='Approaching critical x from above')
    plt.xlabel('t')
    plt.ylabel('f(x, t)')
    plt.legend(loc=2)
    plt.plot()
    
simple_euler_plotting()

def improved_euler(xi, ti, deltat):
    xiplus1 = xi + 0.5*deltat*( function(xi,ti) + function(xi + deltat * function(xi, ti), ti + deltat) )
    return xiplus1

def improved_euler_routine(x0, t0, maxt, step):
    currentt = t0
    currentx = x0
    tvals = []
    xvals = []
    tvals.append(t0)
    xvals.append(x0)
    if t0 < maxt:
        while currentt < maxt:
            currentx = improved_euler(currentx, currentt, step)
            xvals.append(currentx)
            currentt += step
            tvals.append(currentt)
        return tvals, xvals
    elif t0 > maxt:
        while currentt > maxt:
            currentx = improved_euler(currentx, currentt, step)
            xvals.append(currentx)
            currentt += step
            tvals.append(currentt)
        return tvals, xvals
    else:
        raise ValueError('Invalid choice of t0, tmax')
        
def improved_euler_plotting():
    ttest1, xtest1 = improved_euler_routine(0.0655, 0, 3.0, 0.02)
    ttest2, xtest2 = improved_euler_routine(0.0663, 3.0, 0, -0.02)
    tsforcritx = np.linspace(0, 3, 20)
    critx = np.ones((len(tsforcritx)))
    critx = np.multiply(critx, 0.065923)
    plt.figure()
    plt.plot(ttest1, xtest1, 'r', label ='Approaching critical x from below')
    plt.plot(ttest2, xtest2, 'b', label ='Approaching critical x from above')
    plt.xlabel('t')
    plt.ylabel('f(x, t)')
    plt.legend(loc=3)
    plt.plot()
    
improved_euler_plotting()


def runge_kutta(xi, ti, deltat):
    k1 = function(xi,ti)
    k2 = function(xi + 0.5 * deltat * k1, ti + 0.5 * deltat)
    k3 = function(xi + 0.5 * deltat * k2, ti + 0.5 * deltat)
    k4 = function(xi + deltat * k3, ti + deltat)
    xiplus1 = xi + deltat / 6.0 * ( k1 + 2.0 * k2 + 2.0 * k3 + k4 ) 
    return xiplus1

def runge_kutta_routine(x0, t0, maxt, step):
    currentt = t0
    currentx = x0
    tvals = []
    xvals = []
    tvals.append(t0)
    xvals.append(x0)
    if t0 < maxt:
        while currentt < maxt:
            currentx = runge_kutta(currentx, currentt, step)
            xvals.append(currentx)
            currentt += step
            tvals.append(currentt)
        return tvals, xvals
    elif t0 > maxt:
        while currentt > maxt:
            currentx = runge_kutta(currentx, currentt, step)
            xvals.append(currentx)
            currentt += step
            tvals.append(currentt)
        return tvals, xvals
    else:
        raise ValueError('Invalid choice of t0, tmax')

def runge_kutta_plotting():
    ttest1, xtest1 = runge_kutta_routine(0.0655, 0.0, 3.0, 0.02)
    ttest2, xtest2 = runge_kutta_routine(0.0663, 3.0, 0.0, -0.02)
    tsforcritx = np.linspace(0, 3, 20)
    critx = np.ones((len(tsforcritx)))
    critx = np.multiply(critx, 0.065923)
    plt.figure()
    plt.plot(ttest1, xtest1, 'r', label ='Approaching critical x from below')
    plt.plot(ttest2, xtest2, 'b', label ='Approaching critical x from above')
    plt.xlabel('t')
    plt.ylabel('f(x, t)')
    plt.legend(loc=3)
    plt.plot()
    
runge_kutta_plotting()

    

print('''Stephen O'Shea - SN 13321762''')  
    
    