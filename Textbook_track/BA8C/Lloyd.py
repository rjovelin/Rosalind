# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 19:35:19 2018

@author: Richard
"""


import math
import random
import copy

# function to compute the euclidian distance
def EuclidianDist(L1, L2):
    D = 0
    for i in range(len(L1)):
        D += (L1[i] - L2[i])**2
    return math.sqrt(D)

# function to compute the center of gravity
def CenterGravity(data):
    G = []
    for i in range(len(data[0])):
        x = 0
        for j in range(len(data)):
            x += data[j][i]
        G.append(x / len(data))
    return G


# function to assign data points to clusters 
def CentersToClusters(centers, data):
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
    return R

# function to assign centers of gravity to new centers
def ClustersToCenters(R, data):
    centers = []
    for i in R:
        D = [data[j] for j in R[i]]
        centers.append(CenterGravity(D))
    centers.sort()
    return centers
        
# function to implement heuristic search
def ImplementLloyd(k, data):
    
    Data = copy.deepcopy(data)
    # initialize with the first k points in data
    centers = []
    for i in range(k):
        centers.append(Data.pop(0))
    R= CentersToClusters(centers, Data)
    newcenters = ClustersToCenters(R, Data)
    for item in centers:
        Data.append(item)
    centers = newcenters
            
    previous = []
    # do a heuristic search, stop when centers are the same
    while previous != centers:
        print(i)
        previous = centers
        R= CentersToClusters(centers, Data)
        centers = ClustersToCenters(R, Data)

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

    for i in ImplementLloyd(k, data):
    #for i in ImplementLloyd(k, data):
        print(' '.join(list(map(lambda x: str(round(x, 3)), i))))
         