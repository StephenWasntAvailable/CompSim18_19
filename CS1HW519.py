# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 15:13:06 2019

@author: soshe
"""
import numpy as np
import matplotlib.pyplot as plt
import random
import collections as col


def factorial(n):
    if (n == 0):
        return 1
    else:
        return n * factorial(n-1)
    
def poisson(avnval, n):
    return (avnval ** n) * np.exp(-avnval) / factorial(n)

def plot_poisson(avnval):
    nvalues = np.linspace(0, 50, 51)
    poissonarray = np.zeros(len(nvalues))
    for i in range(len(nvalues)-1):
        poissonarray[i] = poisson(avnval, nvalues[i])
    plt.figure()
    plt.plot(nvalues, poissonarray, 'r')
    plt.show()
    
def compute_sum_poisson(avnval):
    nvalues = np.linspace(0, 50, 51)
    poissonarray = np.zeros(len(nvalues))
    for i in range(len(nvalues)-1):
        poissonarray[i] = poisson(avnval, nvalues[i])
    return np.sum(poissonarray)
    
def compute_sum_n_poisson(avnval):
    nvalues = np.linspace(0, 50, 51)
    poissonarray = np.zeros(len(nvalues))
    for i in range(len(nvalues)-1):
        poissonarray[i] = poisson(avnval, nvalues[i]) * nvalues[i]
    return np.sum(poissonarray)

def compute_sum_nsq_poisson(avnval):
    nvalues = np.linspace(0, 50, 51)
    poissonarray = np.zeros(len(nvalues))
    for i in range(len(nvalues)-1):
        poissonarray[i] = poisson(avnval, nvalues[i]) * nvalues[i] * nvalues[i]
    return np.sum(poissonarray)
        
def print_sums():
    avnvalues = np.array([1, 5, 10])
    for i in range(0, len(avnvalues)):
        x1 = compute_sum_poisson(avnvalues[i])
        s1 = "The value of the sum P(n) is %f" %x1
        print(s1)
        x2 = compute_sum_n_poisson(avnvalues[i])
        s2 = "The value of the sum nP(n) is %f" %x2
        print(s2)
        x3 = compute_sum_nsq_poisson(avnvalues[i])
        s3 = "The value of the sum n^2P(n) is %f" %x3
        print(s3)

#print_sums()

def dartboard_sizel_dartsn(l, n):
    hits_array = np.zeros(l)
    i = 1
    while (i <= n):
        x = random.randint(0, l-1)
        hits_array[x] += 1
        i += 1
    return hits_array

def compute_hn(l, n):
    hits_array = dartboard_sizel_dartsn(l, n)
    return np.unique(hits_array, return_counts = True)
    
        
def multiple_trials_dartboard(ntrials, l, n):
    uniques_np, hns_np = compute_hn(l, n)
    uniques = [uniques_np.tolist()]
    hns = [hns_np.tolist()]
    i = 2
    while (i <= ntrials):
        temp_np1, temp_np2 = compute_hn(l, n)
        temp1 = temp_np1.tolist()
        temp2 = temp_np2.tolist()
        uniques.append(temp1)
        hns.append(temp2)
        i += 1
    max_length = 0
    for i in range(0, len(hns)-1):
        if (len(hns[i]) > max_length):
            max_length = len(hns[i])
    for i in range(0, max_length - 1):
        while (len(hns[i]) < max_length):
            hns[i].append(0)
    hn = np.zeros(max_length)
    for i in range(0, max_length):
        for j in range(0, len(hns)):
            hn[i] += hns[j][i]
    print(hn)
    total_throws = np.sum(np.sum(hns))
    print(total_throws)
    pn = np.zeros(max_length)
    for i in range(0, max_length):
        pn[i] = hn[i] / total_throws
    print(pn)
        
multiple_trials_dartboard(5, 100, 50) 