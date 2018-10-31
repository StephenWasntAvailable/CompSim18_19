#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 11:50:05 2017

@author: Stephen
"""

import numpy as np
import math 
import random
import matplotlib.pyplot as plt


# Some constants: lattice dimensions, exchange energy, temperature (K), 
# size scaling factor, tolerance, minimum number of iterations before 
# considering equilibrium, external magnetic field
xy = 20
JInt= 1.0
T = 1.5
kB = 1.38065e-23
sizescalefac = float(xy * xy)
tol = 0.05
min_iter = 10000.0
mag_iter = 500.0
h = 0

# Initialise some matrices with all 0 entries
lattice = [[0 for i in range(xy)] for i in range(xy)]
newlattice = [[0 for i in range(xy)] for i in range(xy)]


# Randomly assigning spins to each lattice position
for i in range(xy):
    for j in range(xy):
        lattice[i][j] = random.choice((-1.0, 1.0))
        


def bc(i):
    if i+1 > xy-1:
        return 0
    if i-1 < 0:
        return xy-1
    else:
        return i


def flipcheck(system, i, j):
    EB = -1 * system[i][j] * (system[bc(i-1)][j] + system[bc(i+1)][j] + system[i][bc(j-1)] + system[i][bc(j+1)])
    EA = -1 * EB
    deltaE = -2 * EB
    if deltaE<=0.0:
        return (system[i][j] * -1), EA
    else:
        randnum = random.random() # Generate a random number between 0 and 1
        if np.exp(((-1.0 * deltaE) / T)) > randnum: # Prob. test for flip
            return system[i][j] * -1.0, EA
        else:
            return system[i][j], EB
        
def timestep(system):
    newsystem = [[0 for i in range(xy)] for i in range(xy)]
    newenergy = [[0 for i in range(xy)] for i in range(xy)]
               
    for i in range(xy):
        for j in range(xy):
            newsystem[i][j] = flipcheck(system, i, j)[0]
            newenergy[i][j] = flipcheck(system, i, j)[1]
    '''
    count = 0
    while count < 10 ** 6:
        rand1 = random.randint(0, xy-1)
        rand2 = random.randint(0, xy-1)
        newsystem[rand1][rand2] = flipcheck(system, rand1, rand2)[0]
        newenergy[rand1][rand2] = flipcheck(system, rand1, rand2)[1]
        count += 1
     '''      
    return newsystem, newenergy
   
# Function to calculate the magnetisation (scaled to range of [-1, 1] by using 
def magnet(system):     # a scaling factor based on lattice size)
    return (np.sum(system)) / sizescalefac    

    
# Function to evolve a system through multiple 'time' steps in order to bring
# the system to equilibrium, with equilibrium being defined as when the 
# change in magnetisation between steps is below a threshold, after a number of
# iterations have been completed  
def toequil(system): 
    iter = 0
    M_old = magnet(system)  
    newsystem = timestep(system)[0] # Evolve system through one step
    iter += 1
    system = newsystem
    M = magnet(system)
    while abs(M-M_old) > tol: #Criterion for equil.
        while iter < min_iter: #Go through a number of iteration before we 
            M = M_old          #check for equilibrium
            newsystem = timestep(system)[0]
            system = newsystem
            M = magnet(system)
            iter +=1
            print(iter)
        M = M_old
        newsystem = timestep(system)[0]
        system = newsystem
        M = magnet(system)
        iter += 1
        print(iter)
    return system
    
equillattice = toequil(lattice)
p = magnet(equillattice)
print(p)

plt.imshow(equillattice, shape = 'circle',interpolation = 'nearest') 
# Produces a figure depicting the lattice as red and blue squares corresponding
# to the spins  
      
