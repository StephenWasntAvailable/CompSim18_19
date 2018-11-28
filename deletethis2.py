#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 10:20:19 2018

@author: Stephen
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as random


import csv
import smopy

#smopy tile server and basic options
smopy.TILE_SERVER = "http://tile.basemaps.cartocdn.com/light_all/{z}/{x}/{y}@2x.png"
smopy.TILE_SIZE = 512
worldmap = smopy.Map((42., -1., 55., 3.), z=4)
plt.imshow(worldmap.img)