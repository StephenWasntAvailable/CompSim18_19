# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 00:00:14 2016

@author: Stephen
"""

import numpy as np
import matplotlib.pyplot as plt
import random

x = 5
y = 5
z = 5
J = 1.0
kB = 1.38065e-23
sizescalefac = float(x * y * z)
T = 1.0
min_iter = 1000.0
mag_iter = 500.0
tol = 0.05

lattice = [[[0 for i in range(x)]for i in range (y)]for i in range(z)]
# Very similar to the function of the same name in the 2D code, despite 
# similarities, this function caused IndexErrors           
def flipcheck(system, up , down , left, right, forward, back):
    EB = -1 * J * (system[i][j][k] * (up + down + left + right + back + forward))
    EA = J * (system[i][j][k] * (up + down + left + right + back + forward))            
    deltaE = EA - EB
    if deltaE <= 0:
        return system[i][j][k] * -1.0, EA
    else:
        randnum = random.random()
        if np.exp(-1.0 * (deltaE)/T) > randnum:
            return system[i][j][k] * -1.0, EA
        else:
            return system[i][j][k], EB
# Essentially exactly the same as the function of the same name in the 2D code,
# just far longer to account for all the boundary conditions            
def timestep(system):
    newsystem = [[[0 for i in range(x)]for i in range (y)]for i in range(z)]
    newenergy = [[[0 for i in range(x)]for i in range (y)]for i in range(z)]
    for i in range(x):
        for j in range(y):
            for k in range(z):
                if i==0 and j==0 and k==0:
                    up=system[i+1][j][k]
                    down=system[x-1][j][k]
                    left=system[i][y-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][z-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==0 and j==0 and k==z-1:
                    up=system[i+1][j][k]
                    down=system[x-1][j][k]
                    left=system[i][y-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][0]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==0 and j==y-1 and k==0:
                    up=system[i+1][j][k]
                    down=system[x-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][0][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][z-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==x-1 and j==0 and k==0:
                    up=system[0][j][k]
                    down=system[i-1][j][k]
                    left=system[i][y-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][z-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==0 and j==y-1 and k==z-1:
                    up=system[i+1][j][k]
                    down=system[x-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][0][k]
                    forward=system[i][j][0]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==x-1 and j==0 and k==z-1:
                    up=system[0][j][k]
                    down=system[i-1][j][k]
                    left=system[i][y-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][0]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==x-1 and j==y-1 and k==0:
                    up=system[0][j][k]
                    down=system[i-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][0][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][z-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==x-1 and j==y-1 and k==z-1:
                    up=system[0][j][k]
                    down=system[i-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][0][k]
                    forward=system[i][j][0]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==0 and j==0:
                    up=system[i+1][j][k]
                    down=system[x-1][j][k]
                    left=system[i][y-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==0 and j==y-1:
                    up=system[i+1][j][k]
                    down=system[x-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][0][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==0 and k==0:
                    up=system[i+1][j][k]
                    down=system[x-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][z-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==0 and k==z-1:
                    up=system[i+1][j][k]
                    down=system[x-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][0]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif j==0 and k==0:
                    up=system[i+1][j][k]
                    down=system[i-1][j][k]
                    left=system[i][y-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][z-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif j==0 and k==z-1:
                    up=system[i+1][j][k]
                    down=system[i-1][j][k]
                    left=system[i][y-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][0]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif j==0 and i==x-1:
                    up=system[0][j][k]
                    down=system[i-1][j][k]
                    left=system[i][y-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif k==0 and i==x-1:
                    up=system[0][j][k]
                    down=system[i-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][z-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif k==0 and j==y-1:
                    up=system[i+1][j][k]
                    down=system[i-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][0][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][z-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif k==z-1 and i==x-1:
                    up=system[0][j][k]
                    down=system[i-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][0]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif k==z-1 and j==y-1:
                    up=system[i+1][j][k]
                    down=system[i-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][0][k]
                    forward=system[i][j][0]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==x-1 and j==y-1:
                    up=system[0][j][k]
                    down=system[i-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][0][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==0:
                    up=system[i+1][j][k]
                    down=system[x-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif j==0:
                    up=system[i+1][j][k]
                    down=system[i-1][j][k]
                    left=system[i][y-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif k==0:
                    up=system[i+1][j][k]
                    down=system[i-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][z-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif i==x-1:
                    up=system[0][j][k]
                    down=system[i-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif j==y-1:
                    up=system[i+1][j][k]
                    down=system[i-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][0][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                elif k==z-1:
                    up=system[i+1][j][k]
                    down=system[i-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][0]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
                    
                else:
                    up=system[i+1][j][k]
                    down=system[i-1][j][k]
                    left=system[i][j-1][k]
                    right=system[i][j+1][k]
                    forward=system[i][j][k+1]
                    back=system[i][j][k-1]
                    newsystem[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[0]
                    newenergy[i][j][k] = flipcheck(system, up , down , left , right, forward, back)[1]
               
    return newsystem, newenergy
#No reference to dimension, so no change from 2D
def magnet(system):
    return (np.sum(system)) / sizescalefac
    
#No reference to dimension, so no change from 2D   
def toequil(system): 
    iter = 0
    M_old = magnet(system)  
    newsystem = timestep(system)[0]
    iter += 1
    system = newsystem
    M = magnet(system)
    while abs(M-M_old) > tol:
        while iter < min_iter:
            M = M_old
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
#No reference to dimension, so no change from 2D    
def netmag(system):
    iter = 0
    netmagarray = np.zeros(mag_iter)
    net_mag_sq = np.zeros(mag_iter)
    while iter < mag_iter:
        M = abs(magnet(system))
        M_sq = M ** 2
        netmagarray[iter] = M
        net_mag_sq[iter] = M_sq
        newsystem = timestep(system)[0]
        iter += 1
        system = newsystem
    return netmagarray, net_mag_sq

# Some changes to bring into 3D, such as adding a dimension to the 2D arrays,
# and replacing the matrix sum and muliplicaiton functions with for loops, 
# maintains the same functionality    
def avenergy(system):
    iter = 0.0
    energymatrix = [[[0 for i in range(x)] for i in range(y)]for i in range(z)]
    cumulativeenergy = [[[0 for i in range(x)] for i in range(y)]for i in range(z)]
    av_cumenergy = [[[0 for i in range(x)] for i in range(y)]for i in range(z)]
    av_cum_en_sq = [[[0 for i in range(x)] for i in range(y)]for i in range(z)]
    while iter < mag_iter:
        newsystem = timestep(system)[0]
        energymatrix = timestep(system)[1]
        e_m_sq = [[[0 for i in range(x)] for i in range(y)]for i in range(z)]
        for i in range(x):
            for j in range(y):
                for k in range(z):
                    e_m_sq[i][j][k] = energymatrix[i][j][k] ** 2
                    cumulativeenergy[i][j][k] = cumulativeenergy[i][j][k] + energymatrix[i][j][k]
        system = newsystem
        iter += 1.0
    for i in range(len(cumulativeenergy[0])):
        for j in range(len(cumulativeenergy[1])):
            for k in range(len(cumulativeenergy[2])):
                av_cumenergy[i][j][k] = cumulativeenergy[i][j][k] / float(mag_iter)
                av_cum_en_sq[i][j][k] = e_m_sq[i][j][k] / float(mag_iter)
    av_energy = (np.average(av_cumenergy)) / sizescalefac
    av_en_sq = (np.average(av_cum_en_sq)) / sizescalefac
    return av_energy, av_en_sq

#No reference to dimension, so no change from 2D                
n = 20
av_en_arr = np.zeros(n)
av_en_arr_sq = np.zeros(n)
av_en_arr_sq2 = np.zeros(n)
CV = np.zeros(n)
for k in range(n):
    T = 1.6 + 0.2 * k
    print T
    equillattice = toequil(lattice)
    av_en_arr[k] = avenergy(equillattice)[0]
    av_en_arr_sq[k] = (av_en_arr[k]) ** 2
    av_en_arr_sq2[k] = avenergy(equillattice)[1] 
    CV[k] = (av_en_arr_sq[k] - av_en_arr_sq2[k]) / (kB * (T**2))
linspace = (np.arange(0, n, 1) /5.0 ) + 1.6
plt.figure()
plt.plot(linspace, CV, 'ro')
plt.show()

#No reference to dimension, so no change from 2D
n = 20
temp = np.zeros(n)
av_mag_arr = np.zeros(n)
av_mag_arr_sq = np.zeros(n)
X = np.zeros(n)
for k in range(n):
    temp[k] = 1.5 + 0.5 * k
    T = temp[k]
    print T
    equillattice = toequil(lattice)
    magarray = netmag(equillattice)[0]
    mag_sq_arr = netmag(equillattice)[1]
    #numpy.histogram(magarray)
    avmag = np.average(magarray)
    avmag_sq = np.average(mag_sq_arr)
    avmag_sq2 = avmag ** 2
    av_mag_arr[k] = avmag
    X[k] = (avmag_sq - avmag_sq2) / (kB * T)
linspace = (np.arange(0, n, 1)/ 2.0) + 1.5
plt.figure()
plt.plot(linspace, X, 'r')
plt.show()