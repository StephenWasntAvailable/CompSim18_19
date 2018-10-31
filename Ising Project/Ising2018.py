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
    
    def __init__(self, size, spins, temperature, modellattice, observables):
        """initial properties used for generating the state and for later processing
        size: integer n, gives a lattice of size nxn
        spins: possible spin values, for collinear case, 1 and -1
        temperature: temperature of the lattice at which to evaluate the model
        modellattice: lattice to use as a starting point for the new lattice i.e if using a lattice from a previous run as a starting point
        modellattice should be 0 if no model being used
        observables: physical properties of interest of the lattice, currently including
        average energy, magnetisation and magnetic susceptability"""
        self.size = size
        self.spins = spins
        self.temp = temperature
        if modellattice == 0:
            self.modellattice = np.array([[random.choice(self.spins) for i in range(self.size)] for i in range(self.size)])
        else:
            self.modellattice = modellattice
        self.observables = observables
            
    def lattice_grid(self):
        """produces a simple plot showing the lattice, different colours indicating different spins"""
        plt.imshow(self.modellattice, shape = 'circle',interpolation = 'nearest')
        
    def update_observables(self):
        """given the current lattice, updates the observables"""
        #average energy
        
        #magnetisation
        
    
               



testlattice = IsingSimple(100, [1.0,-1.0], 2.0, 0, 0)
testlattice.lattice_grid()
#testlattice = IsingSimple(100, [1.0, -1.0], 2.0, [[0,0], [0,0]])
#print(testlattice.modellattice)

        
        
"""
TODO
1. Time evolution to equilibrium
    1.1 Evolve lattice through single time step
        1.1.1 Define function for checking whether a given lattice point should flip
        1.1.2 Evaluate each lattice point using 1.1.1
    1.2 Evolve lattice to equilibrium
        1.2.1 Determine appropriate property which converges at equilibrium
        1.2.2 Function to evaluate lattice for property convergence as time proceeds
    1.3 Update observables after full sweep
        1.3.1 Average energy
        1.3.2 Magnetisation
        1.3.3 Magnetic Suceptability
        
2. Clean up plot of lattice
    2.1 Remove ticks
"""     