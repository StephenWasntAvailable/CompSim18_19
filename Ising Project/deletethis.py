#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 16:56:36 2018

@author: Stephen
"""

import numpy as np
i = 7
j = 22
xy = 22
neighbour_sites = np.zeros((4,2))

def bc(i):  
    return i % (xy+1)
    

s = np.array(([bc(i+1),j], [bc(i-1),j], [i,bc(j+1)], [i,bc(j-1)]))
neighbour_sites = s

print(neighbour_sites)