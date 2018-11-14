# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 14:33:51 2018

@author: Stephen
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy
import random



class IsingSimple:
    """class representing the ising lattice for the case of nearest neighbour interactions
    and colinear spins, in 2 dimensions"""
    
    def __init__(self, spins, constants, modellattice, didflip, observables):
        '''initial properties used for generating the state and for later processing
        spins: possible spin values, for collinear case, 1 and -1'''
        self.spins = spins
        '''constants: 0 = lattice size, 1 = temperature, 2 = external magnetic field, 3 interaction energy J'''
        self.constants = constants
        
        '''modellattice: lattice to use as a starting point for the new lattice i.e
        if using a lattice from a previous run as a starting point
        modellattice should be 0 if no model being used'''
        if modellattice == 0:
            self.modellattice = np.array([[random.choice(self.spins) for i in range(self.constants[0])] for i in range(self.constants[0])])
        else:
            self.modellattice = modellattice
            
        '''didflip: an array of the same size and dim as the modellattice, which will store whether the lattice 
        site flipped it's spin on the previous sweep, as well as te cause of the flip
        i.e. whether the flip was required due to being a lower energy state, or if it was random
        represented by 
        0,1  flip due to energy and random flip respectively resulting in spin up
        2    for no flip
        3,4  flip due to energy and random flip respectively resulting in spin down'''
        self.didflip = np.ones((self.constants[0], self.constants[0]), dtype = int)
        self.didflip = np.multiply(self.didflip, 2)
        
        '''observables: physical properties of interest of the lattice, currently including
        average energy, magnetisation and magnetic susceptability. Also including 'time' for tracking
        number of sweeps of the algorithm a lattice has gone through. Will have time as the 0th index
        time 0
        net magnetisation per site 1
        average energy per site 2'''
        if observables != [0, 0, 0]:
            self.observables = observables
        else:
            self.observables = observables
            IsingSimple.update_observables(self)
            
        '''poss_energies: for storing the possible energy values for the given constants defining the lattice'''
        self.poss_energies = 0
        '''calculate the possible energy values
        for h!=0 there are up to 10 unique values when considering only nearest neighbour interaction
        in general, first index will be for up spin, second for down
        index 0, 1 - for 0 same spins nearby
        index 2, 3 - for 1 same spins nearby
        index 4, 5 - for 2 same spins nearby
        index 6, 7 - for 3 same spins nearby
        index 8, 9 - for 4 same spins nearby''' 
        if self.constants[2] == 0:
            self.poss_energies = np.zeros(10, dtype = float)
            self.poss_energies[0] = 4 * self.constants[3] - self.constants[2]
            self.poss_energies[1] = 4 * self.constants[3] + self.constants[2]
            self.poss_energies[2] = 2 * self.constants[3] - self.constants[2]
            self.poss_energies[3] = 2 * self.constants[3] + self.constants[2]
            self.poss_energies[4] = -1 * self.constants[2]
            self.poss_energies[5] = self.constants[2]
            self.poss_energies[6] = -2 * self.constants[3] - self.constants[2]
            self.poss_energies[7] = -2 * self.constants[3] + self.constants[2]
            self.poss_energies[8] = -4 * self.constants[3] - self.constants[2]
            self.poss_energies[9] = -4 * self.constants[3] + self.constants[2]
            '''similar to above, except there are only 5 unique values when h == 0'''
        else:
            self.poss_energies = np.zeros(5, dtype = float)
            self.poss_energies[0] = 4 * self.constants[3]
            self.poss_energies[1] = 2 * self.constants[3]
            self.poss_energies[2] = 0 
            self.poss_energies[3] = -2 * self.constants[3]
            self.poss_energies[4] = -4 * self.constants[3]
        
        '''poss_flips: object for holding which flips are possible given the constants provided
        i.e some flips are or are not possible depending on whether there is an external field'''
        self.poss_flips_up = 0
        self.poss_flips_down = 0
        if self.constants[2] > 0 and self.constants[2] < 2 * self.constants[3]:
            self.poss_flips_up = [1,1,0,0,0]
            self.poss_flips_down = [1,1,1,0,0]
        elif self.constants[2] < 0 and abs(self.constants[2]) < 2 * self.constants[3]:
            self.poss_flips_up = [1,1,1,0,0]
            self.poss_flips_down = [1,1,0,0,0]
        elif self.constants[2] == 0:
            self.poss_flips_up = [1,1,0,0,0]
            self.poss_flips_down = [1,1,0,0,0]
        else:
            raise ValueError('Currently not allowing for abs(h) greater than or equal to 2 * J')
            
            
            
    def lattice_grid(self):
        '''produces a simple plot showing the lattice, different colours indicating different spins'''
        plt.figure()
        plt.imshow(self.modellattice, shape = 'circle',interpolation = 'nearest')
        plt.show()
        
    def l_b_cnd(self, i):
        '''implements periodic boundary conditions on the lattice'''
        '''keeping this here as legacy, should I run into issues with the modulo method'''
#        if i+1 > self.constants[0]-1:  
#            return 0  
#        if i-1 < 0:  
#            return self.constants[0]-1  
#        else:  
#            return i 
        
        '''using modulo instead, here constants[0] is the size of the lattice'''
        return i % (self.constants[0])
    
    def update_observables(self):
        '''given the current lattice, updates the observables'''
       
        '''average energy per site'''
        net_energy = 0
        for i in range((self.constants[0]-1)):
            for j in range((self.constants[0]-1)):
                neighbour_values = self.neighbouring_sites(i, j)[0]
                net_energy +=  -1 * self.constants[3] * self.modellattice[i][j] * (np.sum(neighbour_values)) - self.constants[2] * self.modellattice[i][j]
        self.observables[2] = net_energy / (self.constants[0] ** 2)
        
        '''to be finished - add routine for updating based on energy from previous sweep
        if no pre-existing value for av. energy'''
#        if self.observables[2] == 0 and self.observables[0] == 0:
#            net_energy = 0
#            for i in range(self.constants[0]-1):
#                for j in range(self.constants[0]-1):
#                    net_energy +=  -1 * self.modellattice[i][j] * (self.modellattice[i][IsingSimple.l_b_cnd(j+1)] + self.modellattice[i][IsingSimple.l_b_cnd(j-1)] + self.modellattice[IsingSimple.l_b_cnd(i+1)][j] + self.modellattice[IsingSimple.l_b_cnd(i-1)][j]) - self.constants[2] * self.modelattice[i][j]
#            self.observables[2] = net_energy / (self.constants[0] ** 2)
        '''if there is an existing energy value and a record of the sites that flipped            '''
#        else:
        '''do not need to fully recalc for neighbouring sites, should be sufficient to adjust by +-2J 
        depending on whether the flipped spin is now the same or different to the neighbour'''
#            net_energy_change = 0
#            flipped_sites = np.argwhere(self.didflip = 1)
#            for k in range(len(flipped_sites)-1):
#                nearby_spins, nearby_sites = self.neighbouring_spins(flipped_sites[k][0],flipped_sites[k][1])
                
#                
#            net_energy_change = net_energy_change
        
        '''magnetisation
        if no pre-existing value for net mag.'''
        if self.observables[1] == 0:
            self.observables[1] = np.sum(self.modellattice) / (self.constants[0] ** 2)
            '''if existing value for net mag. exists + record of sites that flipped'''
        else:
            numdown = np.greater_equal(self.didflip, 3)
            numup = np.less_equal(self.didflip, 1)
            self.observables[1] = ( self.observables[1] + ( ( np.sum(numup) - np.sum(numdown) ) / ( self.constants[0] ** 2 ) ) )
        
        '''updating time'''
        self.observables[0] += 1
            
    def flip_check(self, i, j):
        '''function to check whether a site should be flipped'''
        '''This should be changed for better performance, precalculating possible energy values for the given lattice or below'''
        '''possible implementation for energy check, build np array of neighbouring spins, use np.equal to check how many are the same, if 
        number that are the same is greater than or equal to 3, do not flip, then proceed to random check'''
        h = self.constants[2]
        T = self.constants[1]
        nearby_spins, nearby_sites = self.neighbouring_sites([i],[j])
        number_same = np.sum(np.equal(nearby_spins, self.modellattice[i][j]))
        EB = -1.0 * self.constants[3] * self.modellattice[i][j] * (np.sum(nearby_spins)) - h * self.modellattice[i][j]  
        EA = -1.0 * EB  
        deltaE = EA - EB  
#        if deltaE<=0.0:  
#            self.modellattice[i][j] *= -1.0
#            if self.modellattice[i][j] == 1.0:
#                self.didflip[i][j] = 0
#            else:
#                self.didflip[i][j] = 3
#        else:  
#            randnum = random.random() 
#            if np.exp((-1.0 * deltaE) / T) > randnum: 
#                self.modellattice[i][j] *= -1.0 
#                if self.modellattice[i][j] == 1.0:
#                    self.didflip[i][j] = 1
#                else:
#                    self.didflip[i][j] = 4
#            else:  
#                self.modellattice[i][j] *= 1.0
#                self.didflip[i][j] = 2
        '''alternate implementation'''
        if self.modellattice[i][j] == 1 and self.poss_flips_up[number_same] == 1:
            self.didflip[i][j] = 0
        elif self.modellattice[i][j] == -1 and self.poss_flips_down[number_same] == 1:
            self.didflip[i][j] = 3
        else:  
            randnum = random.random() 
            if np.exp((-1.0 * deltaE) / T) > randnum: 
                self.modellattice[i][j] *= -1.0 
                if self.modellattice[i][j] == 1.0:
                    self.didflip[i][j] = 1
                else:
                    self.didflip[i][j] = 4
            else:  
                self.modellattice[i][j] *= 1.0
                self.didflip[i][j] = 2
        
            
            
    def time_step(self):
        '''function to implement a single metropolis sweep through the lattice
        evolving it though a single 'time' step'''
        '''might need some further modification depending on changes made to other functions
        but seems finished for the moment'''
        xy = self.constants[0]
        for i in range(xy-1):  
            for j in range(xy-1):  
                IsingSimple.flip_check(self, i, j)  
        IsingSimple.update_observables(self)
        
               
    def neighbouring_sites(self, i, j):
        '''function to return neighbouring site indexes and spins'''
        '''indexes currently unused'''
#        neighbour_sites = np.zeros((4,2))
#        s = np.array([[(self.l_b_cnd(i+1)),j], [(self.l_b_cnd(i-1)),j], [i,(self.l_b_cnd(j+1))], [i,(self.l_b_cnd(j-1))]])
#        neighbour_sites = s
#        neighbour_values = np.zeros(4)
#        v = np.array([self.modellattice[self.l_b_cnd(i+1)][j], self.modellattice[self.l_b_cnd(i-1)][j], self.modellattice[i][self.l_b_cnd(j+1)], self.modellattice[i][self.l_b_cnd(j-1)]])
#        neighbour_values = v
        neighbour_sites=np.array([[(self.l_b_cnd(i+1)),j], [(self.l_b_cnd(i-1)),j], [i,(self.l_b_cnd(j+1))], [i,(self.l_b_cnd(j-1))]])
        neighbour_values=np.array([self.modellattice[self.l_b_cnd(i+1)][j], self.modellattice[self.l_b_cnd(i-1)][j], self.modellattice[i][self.l_b_cnd(j+1)], self.modellattice[i][self.l_b_cnd(j-1)]])
        return neighbour_values, neighbour_sites



testlattice = IsingSimple([1.0,-1.0], [100, 2.0, 0.0, 1.0], 0, 0, [0, 0, 0])
#testlattice.lattice_grid()
#testlattice.time_step()
#testlattice.lattice_grid()
z = 0
#mags = []
#ts = []
while z <= 10:
#    ts.append(testlattice.observables[0])
#    mags.append(testlattice.observables[1])
    testlattice.time_step()
    testlattice.lattice_grid()
    
    z += 1
#plt.figure()
#plt.plot(ts, mags, 'ro')
#plt.show()

#nbsites, nbvals = neighbouring_sites(testlattice.modellattice, 3, 5)
#print(nbsites)

#testlattice.lattice_grid()
#testlattice = IsingSimple(100, [1.0, -1.0], 2.0, [[0,0], [0,0]])
#print(testlattice.modellattice)
#testlattice.update_observables()
#print(testlattice.observables[1])
        
        
"""
TODO
1. Time evolution to equilibrium
    1.1 Evolve lattice through single time step
        1.1.1 Define function for checking whether a given lattice point should flip
        1.1.2 Evaluate each lattice point using 1.1.1 DONE
    1.2 Evolve lattice to equilibrium
        1.2.1 Determine appropriate property which converges at equilibrium
        1.2.2 Function to evaluate lattice for property convergence as time proceeds
    1.3 Update observables after full sweep
        1.3.1 Average energy
        1.3.2 Magnetisation DONE
        1.3.3 Magnetic Suceptability
        
2. Clean up plot of lattice
    2.1 Remove ticks
"""     