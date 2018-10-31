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
        self.temp = temp
        self.didflip = didflip
        self.flipcause = flipcause
        
        
class IsingLatticeSimple:
    """Class representing the overall lattice, built up of a set of IsingPointSimple objects"""
    
    def __init__(self, lattice, temp):
        self.lattice = lattice
        self.temp = temp