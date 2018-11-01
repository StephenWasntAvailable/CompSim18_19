#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 12:49:37 2018

@author: Stephen
"""
from scipy import constants
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import scipy
h = constants.Planck
c = constants.speed_of_light
kB = constants.Boltzmann

def planck(wl, temp):
    u1 = ((2 * h * (c ** 2)) / (wl ** 5))
    u2 = (1 / (np.exp((h * c ) / ((wl * kB * temp))) - 1))
    u3 = u1 * u2
    return u3

def planck2(wl, temp):
    return -planck(wl, temp)

def find_max_wl(temp):
    fmin = minimize_scalar(planck2, bracket=(0.2E-7, 2E-6), args=(temp))
    return fmin.x, -fmin.fun


def main():
    wls = np.linspace(1e-8, 2e-6, 200)
    Ts = np.linspace(3000, 8000, 6)

    fig, ax = plt.subplots()
    for temp in Ts:
        ax.plot(wls, planck(wls, temp), label=r'$T='+str(int(temp))+'\,\mathrm{K}$')
        ax.ticklabel_format(axis='both', style='sci', scilimits=(-3,3))
        ax.set_ylabel('Spectral radiance (kJ)')
        ax.set_xlabel('Wavelength (nm)')
    #
    #       
        
        
    ax.legend()
    plt.show()



    x = find_max_wl(8000)
    print(x)
main()
#wls = np.linspace(1e-8, 2e-6, 200)
#Ts = np.linspace(3000, 8000, 6)
#wmaxima = np.zeros(len(Ts))