# -*- coding: utf-8 -*-
"""
Created on Fri May 17 18:27:31 2019

@author: rjovelin
"""



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
            L.append(list(map(lambda x: int(x), line)))
    infile.close()
    return L



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



def ComputeTotalDistance(M):
    '''
    (dict) -> dict
    Return a dictionary with total distance (distance from i to all other leafs) 
    for each leaf in distance matrix M
    '''
    
    TotalDist = {}
    for i in M:
        TotalDist[i] = 0
        for j in M[i]:
            TotalDist[i] += M[i][j]
    return TotalDist


def ComputeDistanceStar(M, TotalDist):
    '''
    (dict, dict) -> dict
    Take a distance matix and a dict with total distance for each leaf in M    
    and return the distance matrix star 
    '''
    
    n = len(M)
    
    Mstar = {}
    for i in M:
        Mstar[i] = {}
        for j in M[i]:
            if i == j:
                Mstar[i][j] = 0
            else:
                d = (n - 2) * M[i][j] - TotalDist[i] - TotalDist[j]
                Mstar[i][j] = d
    return Mstar


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
    
def FindElementsLowestDist(M, Lowest):
    '''
    (dict, int) --> list
    Return a list with leafs having the lowest distance in the distance matrix
    Take a single pair if multiple pairs exist
    '''

    # make a list of clusters with lowest distances
    L = []
    for i in M:
        for j in M[i]:
            if M[i][j] == Lowest:
                return [i, j]
    return L


def UpdateMatrix(M, i, j, m):
    '''
    (dict, list, str) -> dict
    Return a modified distance matrix M adding a new row/column m so that
    Dk,m = Dm,k = (1/2)(Dk,i + Dk,j - Di,j) for any k in M
    '''
    
    # add node m in rows
    M[m] = {}
    for k in M:
        if k != m:
            d = (M[k][i] + M[k][j] - M[i][j]) / 2
            M[k][m] = d
            
        
    # add m  in columns
    for k in list(M.keys()):
        if k != m:
            M[m][k] = M[k][m] 
        # add self-distance
        else:
            M[m][k] = 0
    
    return M


def TrimMatrix(M, i, j):
    '''
    (dict, int, int) -> 
    Return a modified distance matrix M without elements i and j
    '''
    
    # trim matrix and removes elements i and j
    del M[i]
    del M[j]
    ToRemove = {}
    for k in M:
        ToRemove[k] = [i,j]
    for k in ToRemove:
        for m in ToRemove[k]:
            del M[k][m]
    return M




def NeighborJoining(M):
    '''
    (dict) --> dict
    
    '''
    
    # n  ← number of rows in D
    n = len(list(M.keys()))
    if n == 2:
       clusters =list( M.keys())
       # T ← tree consisting of a single edge of length D1,2
       T = {clusters[0]:{clusters[1]:M[clusters[0]][clusters[1]]}, clusters[1]: {clusters[0]:M[clusters[1]][clusters[0]]}}
       return T
    
    # D* ← neighbor-joining matrix constructed from the distance matrix D
    TotalDist = ComputeTotalDistance(M)
    Mstar = ComputeDistanceStar(M, TotalDist)
    
    # find elements i and j such that D*i,j is a minimum non-diagonal element of D*
    Lowest = FindLowest(Mstar)
    i, j = FindElementsLowestDist(Mstar, Lowest)

    # Δ ← (TotalDistanceD(i) - TotalDistanceD(j)) /(n - 2)
    delta = (TotalDist[i] - TotalDist[j]) / (n -2)

    # limbLengthi ← (1/2)(Di,j + Δ)
    limblength_i = (M[i][j] + delta) / 2
    
    # limbLengthj ← (1/2)(Di,j - Δ)
    limblength_j = (M[i][j] - delta) / 2
    
    # add a new row/column m to D so that Dk,m = Dm,k = (1/2)(Dk,i + Dk,j - Di,j) for any k
    m = str(len(M))
    M = UpdateMatrix(M, i, j, m)
    
    # D ← D with rows i and j removed
    # D ← D with columns i and j removed
    M = TrimMatrix(M, i, j) 

    # T ← NeighborJoining(D)
    T = NeighborJoining(M)

    # add two new limbs (connecting node m with leaves i and j) to the tree T
    # assign length limbLengthi to Limb(i)
    # assign length limbLengthj to Limb(j)

    for node in [m, i, j]:
        if node not in T:
            T[node] = {}
    T[m][i] = limblength_i
    T[m][j] = limblength_j
    T[i][m] = limblength_i
    T[j][m] = limblength_j
    
    return T


def Format(Tree):
    
    leafs = [i for i in Tree if type(i) == int]
    Last = max(leafs) + 1
    nodes = [int(i) for i in Tree if type(i) == str]
    nodes.sort()
    nodes = list(map(lambda x: str(x), nodes))
    Map = {}
    for i in nodes:
        Map[i] = Last
        Last += 1
    for i in Map:
        Tree[Map[i]] = Tree[i]
    to_remove = [i for i in Tree if i in Map]
    for i in to_remove:
        del Tree[i]
    L = []
    for i in Tree:
        for j in Tree[i]:
            if j in Map:
                L.append([i, Map[j], j])
    for i in L:
        Tree[i[0]][i[1]] = Tree[i[0]][i[2]]             
    for i in Tree:
        to_remove = [j for j in Tree[i] if type(j) == str]
        for j in to_remove:
            del Tree[i][j]
    L = list(Tree.keys())
    L.sort()
    for i in L:
        for j in Tree[i]:
            print(str(i) + '->' + str(j) + ':' + str(float(Tree[i][j]))) 


def SolvePb(PbFile):
    '''
    (file) -> None
    Print the adjacency list of neibhor-joining tree T 
    '''
       
    L = GetMatrixFromFile(PbFile)
    M = ReadMatrix(L)
    T = NeighborJoining(M)
    Format(T)

