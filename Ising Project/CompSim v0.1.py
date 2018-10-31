# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 12:23:44 2016

@author: Stephen
"""

import numpy 
import math 
import random
import matplotlib.pyplot as plt


# Some constants: lattice dimensions, exchange energy, temperature (K), 
# size scaling factor, tolerance, minimum number of iterations before 
# considering equilibrium
x=20
y=20
J= 1.0
T = 1.5
kB = 1.38065e-23
sizescalefac = float(x * y)
tol = 0.05
min_iter = 1000.0
mag_iter = 500.0

# Initialise some matrices with all 0 entries
lattice = [[0 for i in range(x)] for i in range(y)]
newlattice = [[0 for i in range(x)] for i in range(y)]


# Randomly assigning spins to each lattice position
for i in range(x):
    for j in range(y):
        lattice[i][j] = random.choice((-1.0, 1.0))




# Function that will check for spin flip based on energy difference, returns
# the single element either unchanged or flipped, as well as the energy  
def flipcheck(system, up , down , left , right):
    # Energy before flip (EB), Energy after a flip (EA)
    EB = -1.0 * J * (system[i][j] * up + system[i][j] * down + system[i][j] * 
    left + system[i][j] * right)
    
    EA = J * (system[i][j] * up + system[i][j] * down + system[i][j] * left + 
    system[i][j] * right)
    
    deltaE = EA - EB
    if deltaE<=0: # If flip reduces energy, accept flip immediately
        return system[i][j] * -1.0, EA
    else:
        randnum = random.random() # Generate a random number between 0 and 1
        if numpy.exp(-1.0 * (deltaE /  T)) > randnum: # Prob. test for flip
            return system[i][j] * -1.0, EA
        else:
            return system[i][j], EB
            

# Function to evolve the lattice through one time step where 
# each lattice site is checked for a spin flip, returns the new lattice as well
# as a matrix of energies      
def timestep(system):
    # Initialise some arrays
    newsystem = [[0 for i in range(x)] for i in range(y)]
    newenergy = [[0 for i in range(x)] for i in range(y)]
    for i in range (x):
        for j in range(y):
            # Accounting for all boundary conditions
            if i==0 and j==0:
                up = system[0][y-1]
                down = system[i][j+1]
                left = system[x-1][0]
                right = system[i+1][j]
                newsystem[i][j] = flipcheck(system, up , down , left , right)[0]
                newenergy[i][j] = flipcheck(system, up , down , left , right)[1]
            elif i==x-1 and j==0:
                up = system[x-1][y-1]
                down = system[i][j+1]
                left = system[i-1][j]
                right = system[0][0]
                newsystem[i][j] = flipcheck(system, up , down , left , right)[0]
                newenergy[i][j] = flipcheck(system, up , down , left , right)[1]
            elif i==0 and j==y-1:
                up = system[i][j-1]
                down = system[0][0]
                left = system[x-1][y-1]
                right = system[i+1][j]
                newsystem[i][j] = flipcheck(system, up , down , left , right)[0]
                newenergy[i][j] = flipcheck(system, up , down , left , right)[1]
            elif i==x-1 and j==y-1:
                up = system[i][j-1]
                down = system[x-1][0]
                left = system[i-1][j]
                right = system[0][y-1]
                newsystem[i][j] = flipcheck(system, up , down , left , right)[0]
                newenergy[i][j] = flipcheck(system, up , down , left , right)[1]
            elif j==0:
                up = system[i][y-1]
                down = system[i][j+1]
                left = system[i-1][j]
                right = system[i+1][j]
                newsystem[i][j] = flipcheck(system, up , down , left , right)[0]
                newenergy[i][j] = flipcheck(system, up , down , left , right)[1]
            elif j==x-1:
                up = system[i][j-1]
                down = system[i][0]
                left = system[i-1][j]
                right = system[i+1][j]
                newsystem[i][j] = flipcheck(system, up , down , left , right)[0]
                newenergy[i][j] = flipcheck(system, up , down , left , right)[1]
            elif i==0:
                up = system[i][j-1]
                down = system[i][j+1]
                left = system[x-1][j]
                right = system[i+1][j]
                newsystem[i][j] = flipcheck(system, up , down , left , right)[0]
                newenergy[i][j] = flipcheck(system, up , down , left , right)[1]
            elif i==x-1:
                up = system[i][j-1]
                down = system[i][j+1]
                left = system[i-1][j]
                right = system[0][j]
                newsystem[i][j] = flipcheck(system, up , down , left , right)[0]
                newenergy[i][j] = flipcheck(system, up , down , left , right)[1]
            else:
                up = system[i][j-1]
                down = system[i][j+1]
                left = system[i-1][j]
                right = system[i+1][j]
                newsystem[i][j] = flipcheck(system, up , down , left , right)[0]
                newenergy[i][j] = flipcheck(system, up , down , left , right)[1]
    return newsystem, newenergy
'''    
iter = 0

print "t = ", iter

for i in range(y):
    print lattice[i]


while iter < 100:                
    newlattice = timestep(lattice)[0]
    iter += 1
    
    print " "
    
    print "t = ", iter
    
    for i in range(y):
        print newlattice[i]
    print " "
    
    lattice = newlattice 
'''

# Function to calculate the magnetisation (scaled to range of [-1, 1] by using 
def magnet(system):     # a scaling factor based on lattice size)
    return (numpy.sum(system)) / sizescalefac
    

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
        M = M_old
        newsystem = timestep(system)[0]
        system = newsystem
        M = magnet(system)
        iter += 1
    return system
  
# Function to calculate an array of net magnetisms as the lattice evolves after
# reaching equilibrium, these need to be averaged and scaled to calculate the
# magnetic susceptibilty of the lattice. 
def netmag(system):
    iter = 0
    # initialise some dummy arrays to fill values into as the function runs 
    netmagarray = numpy.zeros(mag_iter)
    net_mag_sq = numpy.zeros(mag_iter)
    while iter < mag_iter:
        M = abs(magnet(system))
        M_sq = M ** 2
        netmagarray[iter] = M    # Place values into the the arrays with 
        net_mag_sq[iter] = M_sq  # positions corresponding to times it has 
        newsystem = timestep(system)[0]   #iterated
        iter += 1
        system = newsystem
    return netmagarray, net_mag_sq
   
# Function to sum two matrices entry by entry, found online at (1)    
def matrixsum(a, b):
    res = []
    for i in range(len(a)): # Where len(a) will return the no. of rows in a
        row = []
        for j in range(len(a[0])): #len(A[0]) returns the no. of columns in a
            row.append(a[i][j]+b[i][j])
        res.append(row)
    return res  
    
# The same as above, just modified to multiply entry by entry rather than add    
def matrixmult(a, b):
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
    iter = 0.0
    # Some dummy matrices that will be filled with values
    energymatrix = [[0 for i in range(x)] for i in range(y)]
    cumulativeenergy = [[0 for i in range(x)] for i in range(y)]
    cum_en_sq = [[0 for i in range(x)]for i in range(y)]
    av_cumenergy = [[0 for i in range(x)] for i in range(y)]
    av_cum_en_sq = [[0 for i in range(x)] for i in range(y)]
    while iter < mag_iter: # No. of times to repeat
        newsystem = timestep(system)[0]
        energymatrix = timestep(system)[1]
        e_m_sq = matrixmult(energymatrix, energymatrix) #Matrix of energies sq
        cum_en_sq = matrixsum(cum_en_sq, e_m_sq) #Matrix holding the 
        system = newsystem      #accumulating values of energy before averaging
        cumulativeenergy = matrixsum(cumulativeenergy, energymatrix)
        iter += 1.0
    for i in range(x):
        for j in range(y):# Averaging the energy and energy squared
            av_cumenergy[i][j] = cumulativeenergy[i][j] / float(mag_iter)
            av_cum_en_sq[i][j] = cum_en_sq[i][j] / float(mag_iter)
    av_energy = (numpy.average(av_cumenergy)) / sizescalefac  # Scaling the 
    av_en_sq = (numpy.average(av_cum_en_sq)) / sizescalefac  # average by lattice
    return av_energy, av_en_sq                               # size


# Routine to gather data about energy at varying temperatures, and then plots 
# heat capacity of the lattice over that temperature range, can also be altered
# slightly to produce a plot of the average energy over the temperature range   
n = 20
av_en_arr = numpy.zeros(n)
av_en_arr_sq = numpy.zeros(n)
av_en_arr_sq2 = numpy.zeros(n)
CV = numpy.zeros(n)
for k in range(n):
    T = 1.5 + 0.5 * k
    print T  # Print statement included mostly just to keep track of progress 
    # while  it is running, as there were some issues with the code getting stuck
    # unable to reach equilibrium
    equillattice = toequil(lattice)  #Bring a lattice to equilibrium
    av_en_arr[k] = avenergy(equillattice)[0] #Calc using equil. lattice
    av_en_arr_sq[k] = (av_en_arr[k]) ** 2
    av_en_arr_sq2[k] = avenergy(equillattice)[1] 
    CV[k] = (av_en_arr_sq2[k] - av_en_arr_sq[k]) / (kB * (T**2)) #Eq for heatcap
linspace = (numpy.arange(0, n, 1) /2.0 ) + 1.5 #Generating an array to plot with
plt.figure()
plt.plot(linspace, av_en_arr, 'ro')
plt.xlabel('Temperature')
plt.ylabel('Average energy scaled')
plt.show()
       
        
'''
# Routine to gather data about net magnetisation at varying temperatures, and 
# then plots the magnetic susceptivity over that temp. range, can also be 
# altered to produce a plot of net magnetisation
n = 20
temp = numpy.zeros(n)
av_mag_arr = numpy.zeros(n)
av_mag_arr_sq = numpy.zeros(n)
X = numpy.zeros(n)
for k in range(n):
    T = 1.5 + 0.5 * k
    print T
    equillattice = toequil(lattice)
    magarray = netmag(equillattice)[0]
    mag_sq_arr = netmag(equillattice)[1]
    #numpy.histogram(magarray)
    avmag = numpy.average(magarray)
    avmag_sq = numpy.average(mag_sq_arr)
    avmag_sq2 = avmag ** 2
    av_mag_arr[k] = avmag
    X[k] = (avmag_sq - avmag_sq2) / (kB * T) # Eq for susceptibility
linspace = (numpy.arange(0, n, 1)/ 2.0) + 1.5  #Generating an array to plot with
plt.figure()
plt.plot(linspace, X, 'ro')
plt.xlabel('Temperature (K)')
plt.ylabel('Magnetic Susceptibility')
plt.show()
'''
'''  
#average = [avmag] * mag_iter
#linspace = numpy.arange(0, mag_iter, 1)
'''

equillattice = toequil(lattice)
plt.imshow(equillattice, shape = 'circle',interpolation = 'nearest') 
# Produces a figure depicting the lattice as red and blue squares corresponding
# to the spins 
