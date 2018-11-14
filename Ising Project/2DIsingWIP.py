#!/usr/bin/env python3  
# -*- coding: utf-8 -*-  
"""  
@author: Stephen 
"""  
  
import numpy as np  
import math   
import random  
import matplotlib.pyplot as plt  
  
  
  
  
# Some constants: lattice dimensions, exchange energy, temperature (K),   
# size scaling factor, tolerance, minimum number of iterations before   
# considering equilibrium, external magnetic field  
xy = 41  
T = 1.5  
sizescalefac = float(xy * xy)  
h = 0  
mag_iter = 300  
  
  
  
  
# Function to calculate the magnetisation (scaled to range of [-1, 1] by using   
def magnet(system):     # a scaling factor based on lattice size)  
    return (np.sum(system)) / sizescalefac      
  
  
# Initialise some matrices with all 0 entries  
lattice = [[0 for i in range(xy)] for i in range(xy)]  
  
  
  
# Randomly assigning spins to each lattice position  
for i in range(xy):  
    for j in range(xy):  
        lattice[i][j] = random.choice((-1.0, 1.0))  
          
q = magnet(lattice)         
  
#  
def bc(i):  
    if i+1 > xy-1:  
        return 0  
    if i-1 < 0:  
        return xy-1  
    else:  
        return i  
  
  
def flipcheck(system, i, j):  
      
    EB = -1. * system[i][j] * (system[bc(i-1)][j] + system[bc(i+1)][j] + system[i][bc(j-1)] + system[i][bc(j+1)]) - h * system[i][j]  
    EA = -1. * EB  
    deltaE = EA - EB  
    if deltaE<=0.0:  
        return system[i][j] * -1.0, EA  
    else:  
        randnum = random.random() # Generate a random number between 0 and 1  
        if np.exp((-1.0 *deltaE) / T) > randnum: # Prob. test for flip  
            return system[i][j] * -1.0, EA  
        else:  
            return system[i][j], EB  
  
def enerlattice(system):  
    newenergy = [[0 for i in range(xy)] for i in range(xy)]  
    for i in range(xy):  
        for j in range(xy):  
            newenergy[i][j] =  -1. * system[i][j] * (system[bc(i-1)][j] + system[bc(i+1)][j] + system[i][bc(j-1)] + system[i][bc(j+1)])  
    return newenergy         
  
def timestep(system):  
    newsystem = system  
    newenergy = [[0 for i in range(xy)] for i in range(xy)]                 
    for i in range(xy):  
        for j in range(xy):  
            newsystem[i][j] = flipcheck(system, i, j)[0]  
            system = newsystem  
    newenergy = enerlattice(newsystem)  
    return newsystem, newenergy  
  
def toequil1(system):  
    newsystem = system  
    newenergy = enerlattice(system)  
    count = 0  
    while count < 10 ** 6:  
        rand1 = random.randint(0, xy-1)  
        rand2 = random.randint(0, xy-1)  
        newsystem[rand1][rand2] = flipcheck(system, rand1, rand2)[0]  
        newenergy[rand1][rand2] = flipcheck(system, rand1, rand2)[1]  
        system = newsystem   
        count += 1  
        #print(count)       
    return newsystem, newenergy  
     
def toequil2(system):  
    newsystem = system  
    newenergy = enerlattice(system)  
    count = 0  
    while count < 30000:  
        for i in range(xy):  
            for j in range(xy):  
                newsystem[i][j] = flipcheck(system, i, j)[0]  
                newenergy[i][j] = flipcheck(system, i, j)[1]  
                system = newsystem  
        count += 1  
    return newsystem, newenergy  
'''''     
# Function to evolve a system through multiple 'time' steps in order to bring 
# the system to equilibrium, with equilibrium being defined as when the  
# change in magnetisation between steps is below a threshold, after a number of 
# iterations have been completed   
def toequil(system):  
    iteration = 0 
    M_old = magnet(system)   
    newsystem = timestep(system)[0] # Evolve system through one step 
    iteration += 1 
    system = newsystem 
    M = magnet(system) 
    while abs(M-M_old) > tol: #Criterion for equil. 
        while iteration < min_iter: #Go through a number of iteration before we  
            M = M_old          #check for equilibrium 
            newsystem = timestep(system)[0] 
            system = newsystem 
            M = magnet(system) 
            print(iteration) 
            iteration +=1 
        M = M_old 
        newsystem = timestep(system)[0] 
        system = newsystem 
        M = magnet(system) 
        iteration += 1 
        print(iteration) 
    return system 
'''     
  
'''''  
equillattice = toequil1(lattice)[0] 
sum1 = np.sum(enerlattice(equillattice)) / sizescalefac 
print(sum1) 
p = magnet(equillattice) 
print(q) 
print(p) 
plt.imshow(equillattice, shape = 'circle',interpolation = 'nearest')  
plt.show() 
 
# Produces a figure depicting the lattice as red and blue squares corresponding 
# to the spins   
'''  
  
def netmag(system):  
    iteration = 0  
    netmagarray = np.zeros(mag_iter)  
    netmag_sq = np.zeros(mag_iter)  
    while iteration < mag_iter:  
        M = abs(magnet(system))  
        M_sq = M ** 2  
        netmagarray[iteration] = M  
        netmag_sq[iteration] = M_sq  
        newsystem = timestep(system)[0]  
        system = newsystem  
        iteration += 1  
    return netmagarray, netmag_sq  
      
def matrixsum(a,b):  
    res = []  
    for i in range(len(a)):  
        row = []  
        for j in range(len(a[0])):  
            row.append(a[i][j]+b[i][j])  
        res.append(row)  
    return res  
  
def matrixmult(a,b):  
    res = []  
    for i in range(len(a)):  
        row = []  
        for j in range(len(a[0])):  
            row.append(a[i][j] * b[i][j])  
        res.append(row)  
    return res  
      
# A function that returns the average of the average energies of a lattice as   
# evolves through timestep() after already being brought to equilibrium  
def avenergy(system):  
    iteration = 0.0  
    # Some dummy matrices that will be filled with values  
    energymatrix = [[0 for i in range(xy)] for i in range(xy)]  
    cumulativeenergy = [[0 for i in range(xy)] for i in range(xy)]  
    cum_en_sq = [[0 for i in range(xy)]for i in range(xy)]  
    av_cumenergy = [[0 for i in range(xy)] for i in range(xy)]  
    av_cum_en_sq = [[0 for i in range(xy)] for i in range(xy)]  
    while iteration < mag_iter: # No. of times to repeat  
        newsystem = timestep(system)[0]  
        energymatrix = timestep(system)[1]  
        e_m_sq = matrixmult(energymatrix, energymatrix) #Matrix of energies sq  
        cum_en_sq = matrixsum(cum_en_sq, e_m_sq) #Matrix holding the   
        system = newsystem      #accumulating values of energy before averaging  
        cumulativeenergy = matrixsum(cumulativeenergy, energymatrix)  
        iteration += 1.0  
    for i in range(xy):  
        for j in range(xy):# Averaging the energy and energy squared  
            av_cumenergy[i][j] = cumulativeenergy[i][j] / float(mag_iter)  
            av_cum_en_sq[i][j] = cum_en_sq[i][j] / float(mag_iter)  
    av_energy = (np.average(av_cumenergy)) / sizescalefac  # Scaling the   
    av_en_sq = (np.average(av_cum_en_sq)) / sizescalefac  # average by lattice  
    return av_energy, av_en_sq  
  
    
# Routine to gather data about energy at varying temperatures, and then plots  
# heat capacity of the lattice over that temperature range, can also be altered 
# slightly to produce a plot of the average energy over the temperature range    
n = 60
av_en_arr = np.zeros(n) 
av_en_arr_sq = np.zeros(n) 
av_en_arr_sq2 = np.zeros(n) 
CV = np.zeros(n) 
for k in range(n): 
    T = 1.0 + 0.05 * k 
    print(T)  # Print statement included mostly just to keep track of progress  
    # while  it is running, as there were some issues with the code getting stuck 
    # unable to reach equilibrium 
    equillattice = toequil1(lattice)[0]  #Bring a lattice to equilibrium 
    av_en_arr[k] = avenergy(equillattice)[0] #Calc using equil. lattice 
    av_en_arr_sq[k] = (av_en_arr[k]) ** 2 
    av_en_arr_sq2[k] = avenergy(equillattice)[1]  
    CV[k] = (av_en_arr_sq2[k] - av_en_arr_sq[k]) /(T**2) #Eq for heatcap 
linspace = np.arange(1, 4, 0.05) #Generating an array to plot with 
plt.figure() 
plt.plot(linspace, av_en_arr, 'ro') 
plt.xlabel('Temperature') 
plt.ylabel('Average Energy per Site') 
plt.show()
 
'''
# Routine to gather data about net magnetisation at varying temperatures, and   
# then plots the magnetic susceptivity over that temp. range, can also be   
# altered to produce a plot of net magnetisation  
n = 60 
temp = np.zeros(n)  
av_mag_arr = np.zeros(n)  
av_mag_arr_sq = np.zeros(n)  
X = np.zeros(n)  
for k in range(n):  
    T = 1 + 0.05 * k  
    print(T)  
    equillattice = toequil1(lattice)[0]  
    magarray = netmag(equillattice)[0]  
    mag_sq_arr = netmag(equillattice)[1]  
    #numpy.histogram(magarray)  
    avmag = np.average(magarray)  
    avmag_sq = np.average(mag_sq_arr)  
    avmag_sq2 = avmag ** 2  
    av_mag_arr[k] = avmag  
    X[k] = (avmag_sq - avmag_sq2) / (T) # Eq for susceptibility  
linspace = (np.arange(1, 4, 0.05))  #Generating an array to plot with  
plt.figure()  
plt.plot(linspace, X, 'ro')  
plt.xlabel('Temperature (K)')  
plt.ylabel('Magnetic Susceptibility') 
'''

