# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 16:38:46 2018

@author: rjovelin
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
def SqrErrDis(centers, data):
    '''
    (list, list) -> tuple
    Take a list of centers and a list of data points and return a tuple with
    the index pf the center with largest distance and the largest distance within a cluster
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
    # compute squared error distortion  
    d = 0
    for i in R:
        for j in R[i]:
            d += (EuclidianDist(centers[i], data[j]))**2
    return d / len(data)                         
    
    
def PrintResults(file):
    infile = open(file)
    k = int(infile.readline().rstrip()[0])
    data, centers = [], []
    for i in range(k):
        L = infile.readline().rstrip().split()
        centers.append(list(map(lambda x: float(x), L)))
    sep = infile.readline()
    assert '-' in sep
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            data.append(list(map(lambda x: float(x), line)))
    infile.close()

    print(SqrErrDis(centers, data))
    