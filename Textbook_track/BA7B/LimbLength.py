# -*- coding: utf-8 -*-
"""
Created on Wed May  1 17:13:43 2019

@author: rjovelin
"""


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
 


# compute the distance between leaf i and its parent node
def DistanceToParent(i, j, k, M):
    '''
    (int, int, int, dict) --> num
    Return the distance between leaf i and its parent node
    using the leaf distances between leafs i, j and k
    '''
    
    D = (M[j][i] + M[i][k] - M[j][k]) / 2
    return D


# compute the limbdistance of leaf to its parent node
def ComputeLimbLength(leaf, M):
    '''
    (int, dict) --> num
    Return the limb distance of leaf, ie, the distance between leaf and its parent
    node using the distance matrix M    
    '''
    
    D = []
    # make a list of leafs other than i
    Leafs = [i for i in M[leaf] if i != leaf]
    for i in range(len(Leafs)-1):
        for j in range(i+1, len(Leafs)):
            D.append(DistanceToParent(leaf, Leafs[i], Leafs[j], M))
    return min(D)    




K = [[0, 13, 21, 22],
[13, 0, 12, 13],
[21, 12, 0, 13],
[22, 13, 13, 0]]


M = ReadMatrix(K)
print(ComputeLimbLength(1, M))




def SolvePb(PbFile):
    infile = open(PbFile)
    infile.readline()
    Leaf = int(infile.readline().rstrip())
    L = []
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = list(map(lambda x: int(x), line.split()))
            L.append(line)
    infile.close()
    M = ReadMatrix(L)
    print(int(ComputeLimbLength(Leaf, M)))