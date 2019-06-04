# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 22:23:59 2019

@author: rjovelin
"""

import math
import random
import copy



# convert list into dict
def ReadMatrix(L):
    '''
    (list) -> dict
    Take a distance matrix in the form of list of lists and return a dictionary
    of leaf1: {leafi: distance} value pairs
    '''
    
    M = {}
    for i in range(len(L)):
        for j in range(len(L[i])):
            if i not in M:
                M[i] = {}
            M[i][j] = L[i][j]
    return M
 

def GetMatrixFromFile(PbFile):
    '''
    (file) -> list
    Read file and return a list with distances among leafs 
    '''
    infile = open(PbFile)
    infile.readline()
    L = []
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            L.append(list(map(lambda x: float(x), line)))
    infile.close()
    return L


def FindLowest(M):
    '''
    (dict) --> int
    Return the lowest distance in the distance matrix M
    '''

    # pick a distance in the matrix to initialize Lowest
    leafs = list(M.keys())
    for i in M[leafs[0]]:
        if i != leafs[0]:
            Lowest = M[leafs[0]][i]
    # find the lowest distance
    for i in M:
        for j in M[i]:
            if i != j and M[i][j] < Lowest:
                Lowest = M[i][j] 
    return Lowest
    
def FindClusters(M, Lowest):
    '''
    (dict, int) --> list
    Return a list of clusters with the Lowest distance between them in the distance matrix M
    '''

    # make a list of clusters with lowest distances
    L = []
    for i in M:
        for j in M[i]:
            if M[i][j] == Lowest:
                L.extend([i, j])
    L = list(set(L))
    #L.sort()
    
    return L


def UpdateMatrix(M, Closest, Cnew, A):
    '''
    (dict, list, str) -> dict
    Return a modified distance matrix M adding node Cnew and average distance
    between clusters in Closest and removing these clusters
    '''
    
    # add node Cnew in columns
    for i in M:
        # get the distance of i from each point in the cluster
        D = (M[i][Closest[0]] * A[Closest[0]] + M[i][Closest[1]] * A[Closest[1]]) / (A[Closest[0]] + A[Closest[1]])
        M[i][Cnew] = D
    
    # add node Cnew in rows
    for i in list(M.keys()):
        if Cnew not in M:
            M[Cnew] = {}
        M[Cnew][i] = M[i][Cnew] 
    # add self-distance
    M[Cnew][Cnew] = 0
    
    # trim matrix and removes elements in Closest
    to_remove = [i for i in Closest]
    for i in to_remove:
        del M[i]
    ToRemove = {}
    for i in M:
        to_remove = [j for j in M[i] if j in Closest]
        ToRemove[i] = to_remove
    for i in ToRemove:
        for j in ToRemove[i]:
            del M[i][j]
    return M



# implement breadth-first-search for finding the path between 2 nodes
def BFS(tree, f, goal):
    '''
    (dict, str, str) -> list
    Take a dictionary representing with interactions among node, and return
    a list of nodes representing the path betwen nodes f and goal
    '''
    # create a queue to add the discovered nodes that are not yet visited
    # keep track of the path visited to reach discovered node
    queue = [(f, [f])]
    visited = []
    while len(queue) != 0:
        # get a node to visit. keep track of the path up to node
        (node, path) = queue.pop(0)
        if node not in visited:
            visited.append(node)
            # add node children to queue
            for child in tree[node]:
                if child == goal:
                    path.append(child)
                    return path
                else:
                    queue.append((child, path + [child]))



def HierarchicalClustering(M):
    '''
    Input: A space separated n x n distance matrix.
    Output: An adjacency list for the ultrametric tree returned by UPGMA.
    Edge weights should be accurate to two decimal places
    '''
    
    # initialize tree
    T = {}
    # get the last leaf in the matrix
    Last = list(M.keys())
    Last.sort()
    Last = Last[-1]
    
    # construct a list of single-element clusters
    Clusters = [i for i in M]
    # make a list of leafs
    leafs = [i for i in M]
    leafs.sort()
    
    # keep track of the weight of each cluster, ie the number of leafs 
    A = {}
    for i in leafs:
        A[i] = 1
    
    # construct a graph T with n isolated nodes labeled by single-elements
    T = {}
    for i in Clusters:
        T[i] = {}
    # create a graph to record the age of every node
    Age = {}
    for i in Clusters:
        Age[i] = 0
    
    while len(Clusters) > 1:
        # find the two closest clusters Ci and Cj
        Lowest = FindLowest(M)
        Closest = FindClusters(M, Lowest)
              
        # merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
        Cnew = str(Closest[0]) + ',' + str(Closest[1])
        
             
        
        # use 1-based numbering to print output clusters
        print(' '.join(list(map(lambda x: str(x), list(map(lambda x: int(x) + 1, Cnew.split(',')))))))
        
        
        
        
               
        # keep track of leaf size at each merged node
        A[Cnew] = 0
        for i in Closest:
            A[Cnew] += A[i]
        #add a new node labeled by cluster Cnew to T
        T[Cnew] = {}
        #connect node Cnew to Ci and Cj by directed edges
        for i in Closest:
            T[Cnew][i] = Lowest / len(Closest)
            T[i][Cnew] = Lowest / len(Closest)
        #Age(Cnew) ← DCi, Cj / 2
        Age[Cnew] = Lowest / 2
        # remove Ci and Cj from Clusters
        for i in Closest:
            Clusters.remove(i)
        # remove the rows and columns of D corresponding to Ci and Cj
        # add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters
        M = UpdateMatrix(M, Closest, Cnew, A)
        # add Cnew to Clusters
        Clusters.append(Cnew)
       
    # root ← the node in T corresponding to the remaining cluster
    Root = Clusters[0]
    
    # for each edge (v, w) in T
    # length of (v, w) ← Age(v) - Age(w)
    # use breadth-first-search to find the path between root and each leaf
    for i in leafs:
        if i in T[Root]:
            path = [Root, i]
        else:
            path = BFS(T, Root, i)
            # compute edge lenth along the path
        for j in range(len(path)-1):
            edgeLength = Age[path[j]] - Age[path[j+1]]
            # update Tree length
            T[path[j]][path[j+1]] = edgeLength
            T[path[j+1]][path[j]] = edgeLength
    return T  
    
    
def SolvePb(PbFile):
    '''
    
    The adjacency list must have consecutive integer node labels starting from 0.
    The n leaves must be labeled 0, 1, ..., n - 1 in order of their appearance
    in the distance matrix. Labels for internal nodes may be labeled in any order
    but must start from n and increase consecutively
    
    '''
    L = GetMatrixFromFile(PbFile)
    M = ReadMatrix(L)
    Tree = HierarchicalClustering(M)
    return Tree

#    nodes = list(Tree.keys())
#    nodes.sort()
#    for i in nodes:
#        for j in Tree[i]:
#            print(str(i) + '->' + str(j) + ':' + str(round(Tree[i][j], 3)))


