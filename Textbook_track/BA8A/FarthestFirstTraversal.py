# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 19:32:35 2018

@author: Richard
"""

import math
import random


# function to compute the euclidian distance
def EuclidianDist(L1, L2):
    D = 0
    for i in range(len(L1)):
        D += (L1[i] - L2[i])**2
    return math.sqrt(D)

# function to find the max distance between centers and data points 
def GetMaxDist(centers, data):
    '''
    (list, list) -> tuple
    Take a list of centers and a list of data points and return a tuple with
    the index of the center with largest distance and the largest distance within a cluster
    '''
    
    # compute the minimum distance between centers and all data points
    # group centers and data based on minimum distance
    N ={}
    for i in range(len(data)):
        D = []
        for j in range(len(centers)):
            D.append([EuclidianDist(data[i], centers[j]), i, j])
        D.sort()
        assert D[0][1] not in N
        N[D[0][1]] = D[0][2]
    # reverse dict to get all data points most nearest per center, ie identify clusters
    R = {}
    for i in N:
        if N[i] in R:
            R[N[i]].append(i)
        else:
            R[N[i]] = [i]
    k, dmax = '', 0      
    for i in R:
        for j in R[i]:
            d = EuclidianDist(centers[i], data[j])
            if d > dmax:
                dmax = d
                k = j
    return k, dmax                         
    
    
def FarthestFirstTraversal(k, data):
    centers = []
    while len(centers) < k:
        if len(centers) == 0:
            centers.append(data.pop(0))
        else:
            i, dmax = GetMaxDist(centers, data)
            centers.append(data.pop(i))
    return centers               
                
def PrintResults(file):
    infile = open(file)
    k = int(infile.readline().rstrip()[0])
    data = []
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            data.append(list(map(lambda x: float(x), line)))
    infile.close()

    d = FarthestFirstTraversal(k, data)
    for i in d:
        print(' '.join(list(map(lambda x: str(x), i))))


