# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 21:22:36 2018

@author: Stephen
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
import random



class IsingPointSimple:
    """Class representing each individual lattice point for the case of nearest neighbour interactions
    and colinear spins, in 2 dimensions"""
    
    def __init__(self, x, y, spin, didflip, flipcause):
        self.x = x
        self.y = y
        self.spin = spin
        self.didflip = didflip
        self.flipcause = flipcause
        
        
class IsingLatticeSimple:
    """Class representing the overall lattice, built up of a set of IsingPointSimple objects"""
    
    def __init__(self, size, temp, spins, lattice):
        if lattice == 0:
            self.size = size
            self.temp = temp
            self.spins = spins
            self.lattice = np.zeros((self.size, self.size))
            for i in range(size):
                for j in range(size):
                    self.lattice[i][j] = IsingPointSimple(i, j, random.choice(self.spins), 0, 0)
        else:
            self.size = lattice.size
            self.temp = lattice.temp
            self.spins = lattice.spins
            self.lattice = lattice.lattice
#        selfnp.array([[random.choice(self.spins) for i in range(self.size)] for i in range(self.size)])
            
    def show_lattice(self):
        plt.imshow(self.lattice, shape = 'circle',interpolation = 'nearest')
        
        
testlattice = IsingLatticeSimple(10, 2, [1, -1], 0)
testlattice.show_lattice()