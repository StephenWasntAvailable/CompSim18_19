# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 15:13:06 2019

@author: soshe
"""
print("Stephen O'Shea - JSPY")
print("SN : 13321762")

import numpy as np
import matplotlib.pyplot as plt
import random



def factorial(n):
    if (n == 0):
        return 1
    else:
        return n * factorial(n-1)
    
def poisson(avnval, n):
    return (avnval ** n) * np.exp(-avnval) / factorial(n)

def calc_poisson_array(avnval):
    nvalues = np.linspace(0, 50, 51)
    poissonarray = np.zeros(len(nvalues))
    for i in range(len(nvalues)-1):
        poissonarray[i] = poisson(avnval, nvalues[i])
    return poissonarray

def plot_poisson(avnvalues):
    nvalues = np.linspace(0, 50, 51)
    colours = ['blue' , 'green', 'red', 'cyan', 'magenta']
    j = 0
    plt.figure()
    for i in avnvalues:
        temp_array = calc_poisson_array(i)
        plt.plot(nvalues, temp_array, label = "<n> = %d" %i, color = colours[j], marker = '.', lw = 0)
        j += 1
    plt.xlabel('n')
    plt.ylabel('P(n)')
    plt.title('Poisson Distributions for differing average n <n>')
    plt.legend(loc = 'upper right')
    plt.show()
        
plot_poisson([1, 5, 10, 15, 20])      
        
 
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
        poissonarray[i] = poisson(avnval, nvalues[i]) * (nvalues[i] ** 2)
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

print_sums()

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
            hn[i] += boundary(hns, i, j)
    print(hn)
    total_throws = np.sum(np.sum(hns))
    print(total_throws)
    pn = np.zeros(max_length)
    for i in range(0, max_length):
        pn[i] = hn[i] / total_throws
    print("Probabilities: ", pn)
    average_hits_per_site = 0
    for i in range(0, max_length - 1):
        average_hits_per_site += hn[i] * i
    average_hits_per_site = average_hits_per_site / total_throws
    print("Average hits per site: %f" %average_hits_per_site)
    pn_theoretical = np.zeros(max_length)
    for i in range(0, max_length):
        pn_theoretical[i] = poisson(average_hits_per_site, i)
    print("Probabilities from Poisson Distribution: ", pn_theoretical)
    
def boundary(hns, i, j):
    if (j >= len(hns)):
        return 0
    elif (i >= len(hns[j])):
        return 0
    else:
        return hns[j][i]
    
     
multiple_trials_dartboard(10, 100, 50) 