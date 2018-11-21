# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 14:33:51 2018

@author: Stephen
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy
import random
import time



class IsingSimple:
    """class representing the ising lattice for the case of nearest neighbour interactions
    and colinear spins, in 2 dimensions"""
    
    def __init__(self, spins, constants, modellattice, didflip, observables):
        '''initial properties used for generating the state and for later processing
        spins: possible spin values, for collinear case, 1 and -1'''
        self.spins = spins
        '''constants: 0 = lattice size, 1 = temperature, 2 = external magnetic field, 3 interaction energy J, 4 boltzmann const'''
        self.constants = constants
        '''this is not used anywhere, just to confirm that I'm treating Boltzmann const as equal to unity'''
        self.constants[4] = 1
        
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
        if observables != [0, 0, 0, 0, 0]:
            self.observables = observables
        else:
            self.observables = observables
            IsingSimple.update_observables(self)
            
        '''poss_energies: for storing the possible energy values for the given constants defining the lattice'''
        ''' not currently using these, as they are causing problems with how the lattice evolves'''
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
        '''parameters determining how the to_equilibrium method functions'''
        '''index 0: the tolerance when calculating the new mag - old mag, values below this tolerance
        are taken to indicate a system near equilibrium
        index 1: the minimum number of metropolis sweeps before the above mentioned tolerance checking 
        starts
        index 2: number of consecutive sub threshold values in order to be considered at equil'''
        '''to change these parameters, see the change_equil_params method'''
        self.equilibrium_parameters = [3.0 / (self.constants[0] ** 2), 1000, 3]
        '''arrays for storing data as the lattice evolves'''
        self.energy_per_site=[]
        self.netmag_per_site=[]
        self.energy_sq_per_site=[]
        self.netmag_sq_per_site=[]    
            
            
    def lattice_grid(self):
        '''produces a simple plot showing the lattice, different colours indicating different spins'''
        plt.figure()
        plt.imshow(self.modellattice, shape = 'circle',interpolation = 'nearest')
        plt.show()
        
    
    def update_observables(self):
        '''given the current lattice, updates the observables'''
        '''full recalc of energy after every sweep'''
        '''average energy per site'''
        net_energy = 0
        for i in range((self.constants[0]-1)):
            for j in range((self.constants[0]-1)):
                neighbour_values, neighbour_sites = neighbouring_sites(self, i, j)
                net_energy +=  -1 * self.constants[3] * self.modellattice[i][j] * (np.sum(neighbour_values)) - self.constants[2] * self.modellattice[i][j]
        self.observables[2] = net_energy / (self.constants[0] ** 2)
        
        '''if there is an existing energy value and a record of the sites that flipped'''
        '''do not need to fully recalc for neighbouring sites, should be sufficient to adjust by +-2J 
        depending on whether the flipped spin is now the same or different to the neighbour'''
#        if self.observables[2] == 0 and self.observables[0] == 0:
#            net_energy = 0
#            for i in range(self.constants[0]-1):
#                for j in range(self.constants[0]-1):
#                    neighbour_values, neighbour_sites = neighbouring_sites(self, i, j)
#                    net_energy +=  -1 * self.constants[3] * self.modellattice[i][j] * (np.sum(neighbour_values)) - self.constants[2] * self.modellattice[i][j]
#            self.observables[2] = net_energy / (self.constants[0] ** 2)
#        else:
#            net_energy_change = 0
#            flipped_sites = np.argwhere(self.didflip!=2)
#            for k in range(len(flipped_sites)-1):
#                neighbour_values, neighbour_sites = neighbouring_sites(self, flipped_sites[k][0],flipped_sites[k][1])
#                number_same = np.sum(np.equal(neighbour_values, self.modellattice[flipped_sites[k][0]][flipped_sites[k][1]])) 
#                deltaE1 = 2.0 * self.constants[3] * self.modellattice[flipped_sites[k][0]][flipped_sites[k][1]] * (np.sum(neighbour_values)) 
#                deltaE2 = 2.0 * self.constants[2] * self.modellattice[flipped_sites[k][0]][flipped_sites[k][1]]
#                deltaE3 = -2.0 * self.constants[3] * number_same
#                deltaE4 = 2.0 * self.constants[3] * (4 - number_same)
#                net_energy_change = net_energy_change + deltaE1 + deltaE2 + deltaE3 + deltaE4
#
#            self.observables[2] = self.observables[2] + (net_energy_change / (self.constants[0] ** 2))
#            print(net_energy_change)
        
        '''magnetisation
        if no pre-existing value for net mag.'''
        if self.observables[1] == 0:
            self.observables[1] = np.sum(self.modellattice) / (self.constants[0] ** 2)
            '''if existing value for net mag. exists + record of sites that flipped'''
        else:
            numdown = np.greater_equal(self.didflip, 3)
            numup = np.less_equal(self.didflip, 1)
            self.observables[1] = self.observables[1] + ( ( np.sum(numup) - np.sum(numdown) ) / ( self.constants[0] ** 2 ) ) 
        
        '''magnetisation squared'''
        if self.observables[3] == 0:
            self.observables[3] = (np.sum(self.modellattice)/ (self.constants[0] ** 2) ) ** 2            
        else:
            numdown = np.greater_equal(self.didflip, 3)
            numup = np.less_equal(self.didflip, 1)
            self.observables[3] =  self.observables[3] + ( ( (np.sum(numup) - np.sum(numdown)) ** 2 ) / ( self.constants[0] ** 2 ) ) 
            
        '''energy squared'''
        net_energy_sq = 0
        for i in range((self.constants[0]-1)):
            for j in range((self.constants[0]-1)):
                neighbour_values, neighbour_sites = neighbouring_sites(self, i, j)
                net_energy_sq +=  -1 * self.constants[3] * self.modellattice[i][j] * (np.sum(neighbour_values)) - self.constants[2] * self.modellattice[i][j]
        self.observables[4] = ( (net_energy_sq) / (self.constants[0] ** 2) ) ** 2
        
        '''updating time'''
        self.observables[0] += 1
            
    def flip_check(self, i, j):
        '''function to check whether a site should be flipped'''
        '''This should be changed for better performance, precalculating possible energy values for the given lattice or below'''
        '''possible implementation for energy check, build np array of neighbouring spins, use np.equal to check how many are the same, if 
        number that are the same is greater than or equal to 3, do not flip, then proceed to random check'''
        h = self.constants[2]
        T = self.constants[1]
        nearby_spins, nearby_sites = neighbouring_sites(self,i,j)
        number_same = np.sum(np.equal(nearby_spins, self.modellattice[i][j]))
        deltaE = 2.0 * self.constants[3] * self.modellattice[i][j] * (np.sum(nearby_spins)) + 2 * h * self.modellattice[i][j]  
#        if h == 0:
#            EB = self.poss_energies[0 + number_same]
#        elif self.modellattice == 1.0:
#            EB = self.poss_energies[0 + 2 * number_same]
#        else:
#            EB = self.poss_energies[1 + 2 * number_same] 
        if deltaE<0.0:  
            self.modellattice[i][j] *= -1.0
            if self.modellattice[i][j] == 1.0:
                self.didflip[i][j] = 0
            else:
                self.didflip[i][j] = 3
        elif deltaE > 0:  
            randnum = random.random() 
            if np.exp((-1.0 * deltaE) / T) > randnum: 
                self.modellattice[i][j] *= -1.0 
                if self.modellattice[i][j] == 1.0:
                    self.didflip[i][j] = 1
                else:
                    self.didflip[i][j] = 4
            else:  
                self.didflip[i][j] = 2
        else:
            self.didflip[i][j] = 2
        '''alternate implementation'''
        '''possibly slower than above method'''
#        if self.modellattice[i][j] == 1 and self.poss_flips_up[number_same] == 1:
#            self.modellattice[i][j] *= -1.0 
#            self.didflip[i][j] = 0
#        elif self.modellattice[i][j] == -1 and self.poss_flips_down[number_same] == 1:
#            self.modellattice[i][j] *= -1.0 
#            self.didflip[i][j] = 3
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
        
    
    def to_equilibrium(self, displayinfo = 0):
        '''function to attempt to bring the lattice to equilibrium, based on the net magnetisation
        after some warm up sweeps, the method checks for a series of consecutive sweeps wherein
        the net magnetisation of the lattice only changes by a small amount
        the number of warm up sweeps, definition of 'small amount', and number of consecutives
        can all be changed using the update_equil_params method'''
        tol = self.equilibrium_parameters[0]
        min_consec_below_tol = self.equilibrium_parameters[2]
        min_iter = self.equilibrium_parameters[1]
        a = 'Threshold between steps: %f' %tol
        b = 'Minimum number of iterations before checking for equilibrium: %i' %min_iter
        c = 'Minimum number of iterations with change below threshold to consider at equilibrium: %i' %min_consec_below_tol
        d = 'To alter any of these parameters before running again, use the update_equil_params(newmin, newtol, newminconsec) method'
        if displayinfo == 1:
            print(b)
            print(a)
            print(c)
            print(d)
        i = 0
        consec_below_tol = 0
        while i < min_iter:
            self.time_step()
            i += 1
        M_new = self.observables[1]
        while consec_below_tol <= min_consec_below_tol:
            M_old = M_new
            self.time_step()
            M_new = self.observables[1]
            if abs(M_new - M_old) <= tol:
                consec_below_tol += 1
#                print('+1')
            else:
                consec_below_tol = 0
#                print('0')
    
    def update_equil_params(self, newmin, newtol, minconsec):
        '''method for updating the equilibrium parameters for the current lattice'''
        '''first calculated in terms of the size of the lattice, others just a number'''
        self.equilibrium_parameters[0] = newtol / self.constants[0] ** 2
        self.equilibrium_parameters[1] = newmin 
        self.equilibrium_parameters[2] = minconsec
        a = 'Threshold between steps: %f' %self.equilibrium_parameters[0]
        b = 'Minimum number of iterations before checking for equilibrium: %i' %self.equilibrium_parameters[1]
        c = 'Minimum number of iterations with change below threshold to consider at equilibrium: %i' %self.equilibrium_parameters[2]
        print(b, a, c)
        
    
    def equilibrium_evolution(self, num_sweeps):
        '''evolves the lattice while at equilibrium, and stores values of observables 
        inn arrays defined in the __init__, for use with graphing etc.'''
        i = 0
        self.energy_per_site.append(self.observables[2])
        self.energy_sq_per_site.append(self.observables[4])
        self.netmag_per_site.append(self.observables[1])
        self.netmag_sq_per_site.append(self.observables[3])
        while i < num_sweeps:
            self.time_step()
            self.energy_per_site.append(self.observables[2])
            self.energy_sq_per_site.append(self.observables[4])
            self.netmag_per_site.append(self.observables[1])
            self.netmag_sq_per_site.append(self.observables[3])
            
            
            
            i += 1
            
    def update_constants(self, newt, newh):
        '''method to update the 'constants' of the lattice, in particular the 
        temperature and external magnetic field, can be used for seeing how
        observables change with the change of these parameters'''
        if newt != self.constants[1]:
            self.constants[1] = newt
            a = 'New Temperature: %f' %newt
            print(a)
        if newh != self.constants[2]:
            self.constants[2] = newh
            b = 'New External Field: %f' %newh
            print(b)
        
            

def neighbouring_sites(ising, i, j):
    '''function to return neighbouring site indexes and spins'''
    '''indexes currently unused'''
    '''syntax a bit of a mess here, but was running into issues with simpler versions
    and this method works reasonably well'''
    
    neighbour_sites = np.zeros((4,2))
    s = [[(i+1)%(ising.constants[0]),j], [(i-1)%(ising.constants[0]),j], [i,(j+1)%(ising.constants[0])], [i,(j-1)%(ising.constants[0])]]
    neighbour_sites = np.array(s)
    neighbour_values = np.zeros(4)
    v = [ising.modellattice[(i+1)%(ising.constants[0])][j], ising.modellattice[(i-1)%(ising.constants[0])][j], ising.modellattice[i][(j+1)%(ising.constants[0])], ising.modellattice[i][(j-1)%(ising.constants[0])]]
    neighbour_values = np.array(v)
#    neighbour_sites=np.array([[(self.l_b_cnd(i+1)),j], [(self.l_b_cnd(i-1)),j], [i,(self.l_b_cnd(j+1))], [i,(self.l_b_cnd(j-1))]])
#    neighbour_values=np.array([self.modellattice[self.l_b_cnd(i+1)][j], self.modellattice[self.l_b_cnd(i-1)][j], self.modellattice[i][self.l_b_cnd(j+1)], self.modellattice[i][self.l_b_cnd(j-1)]])
    return neighbour_values, neighbour_sites

def l_b_cnd(size, i):
    '''implements periodic boundary conditions on the lattice'''
    '''using modulo instead, here size is the size of the lattice'''
    '''no longer really using this, the boundary conditions are built into the neighbouring_sites function'''
    return i % size




"""
will come back to this if I have time / the quicker to code method is too slow and an alternative is required
class ToEquilibrium(IsingSimple):
    '''class to execute simulated annealing in order to bring the lattice to equilibrium'''
    def __init__(self, ising_ob, T=-1, alpha=-1, stopping_T=-1, stopping_iter=-1):
        
        super().__init__(ising_ob.spins, ising_ob.constants, ising_ob.modellattice, ising_ob.didflip, ising_ob.observables)
        
        self.currentgrid = self.modellattice
        
        self.bestgrid = np.copy(self.modellattice)
        
        #'fitness' or energy of solutions
        self.cur_fitness = self.fitness()
        self.initial_fitness = self.cur_fitness
        self.best_fitness = self.cur_fitness
        self.fitness_list = [self.cur_fitness]

        #annealing parameters (use defaults if not set...)
        self.T = 1.0E6 if T == -1 else T
        self.alpha = 0.999 if alpha == -1 else alpha
        self.stopping_temperature = 0.0001 if stopping_T == -1 else stopping_T
        self.stopping_iter = 100000 if stopping_iter == -1 else stopping_iter
        self.iteration = 1
        
    def fitness(self):
        net_energy = 0
        for i in range((self.constants[0]-1)):
            for j in range((self.constants[0]-1)):
                neighbour_values = neighbouring_sites(self, i, j)[0]
                net_energy +=  -1 * self.constants[3] * self.modellattice[i][j] * (np.sum(neighbour_values)) - self.constants[2] * self.modellattice[i][j]
        fitness_value = net_energy / (self.constants[0] ** 2)
        return fitness_value
"""

'''testing area'''
#start=time.time()
#
#
#testlattice = IsingSimple([1.0,-1.0], [20, 2.0, 0.0, 1.0, 0], 0, 0, [0, 0, 0])
#testlattice.to_equilibrium()
#print(testlattice.observables[0])
#testlattice.lattice_grid()
#testlattice.equilibrium_evolution(100)
#print(testlattice.observables[0])
#
#
#
#end = time.time()    
#print(end - start)
'''end testing area'''       
"""
TODO
1. Time evolution to equilibrium
    1.1 Evolve lattice through single time step
        1.1.1 Define function for checking whether a given lattice point should flip DONE
        1.1.2 Evaluate each lattice point using 1.1.1 DONE
    1.2 Evolve lattice to equilibrium
        1.2.1 Determine appropriate property which converges at equilibrium DONE
        1.2.2 Function to evaluate lattice for property convergence as time proceeds DONE
    1.3 Update observables after full sweep
        1.3.1 Average energy DONE
        1.3.2 Magnetisation DONE
        1.3.3 Magnetic Suceptability
        
2. Clean up plot of lattice
    2.1 Remove ticks
"""     