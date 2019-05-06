# -*- coding: utf-8 -*-
"""
Created on Fri May  3 17:12:32 2019

@author: rjovelin
"""


# The following recursive algorithm, called AdditivePhylogeny, finds the simple
# tree fitting an n x n additive distance matrix D. It relies on a function Limb(D, j)
# that computes LimbLength(j) for a leaf j based on the distance matrix D.
# Rather than selecting an arbitrary leaf j from Tree(D) for trimming,
# AdditivePhylogeny selects leaf n (corresponding to the last row and column of D).


import copy
import random

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


def SelectLeaf(M):
    '''
    (dict) -> int
    Return the leaf corresponding to the last row and column of M
    Precondition: leafs are labeled with integers
    '''
    
    a = M.keys()
    a.sort()
    return a[-1]




def ComputeDBald(M, n, d):
    for i in M[n]:
        if i != n:
            M[n][i] = M[n][i] - d
    for i in M:
        if i != n:
            M[i][n] = M[i][n] -d
    return M





def AdditivePhylogeny(M, n):
    '''
    
    Given: n and a tab-delimited n x n additive matrix.
    Return: A weighted adjacency list for the simple tree fitting this matrix.

    
    
    
    '''
    
    print('n', n)
    
    
    if n == 1:
        # return the tree consisting of a single edge of length D1,2
        return {0:{1:M[0][1]}, 1: {0:M[1][0]}}
    
    # compute distance between n and its parent node
    d = ComputeLimbLength(n, M)
    # compute Dbald in which n has bald limb length of 0
    # by substracting d to each element in row n and column n 
    M = ComputeDBald(M, n, d)
        
    print(n, d)
    #print(M)
    
    
    # find node i and k such that Di,k = Di,n + Dn,k    
    x = -1
    trio = []
    
    Found = False
    for i in M:
        if Found == False:
            for k in M:
                if i != n and k != n and i != k:
                    if M[i][k] == M[i][n] + M[n][k]:
                        trio.append([i, n, k])
                        Found = True
                        break
    print('trio', trio)
    trio = trio[0]
    
    #print(trio)
    
    # get the distance between i and n
    x, y, i, k = M[i][n], M[n][k], trio[0], trio[-1]
    
    print('node', 'distance', i, x)
    
    # compute Dtrim by removing row n and column n from M
    Dtrim = M
    del Dtrim[n]
    for i in Dtrim:
        del Dtrim[i][n]
        
    #$print(Dtrim)    
    
    # compute the tree corresonding to Dtrim
    T = AdditivePhylogeny(Dtrim, n - 1)
    
    print('tree', T)
    
    
    # check if a node exists at distance x from i
    # if not, create a new node and attach n
    # if exists, attach n to that node
    
#    # make a list of internal nodes
#    internal = [node for node in T if len(T[node]) > 1]
#    internal.sort()
#    if len(internal) == 0:
#        # new node label
#        v = 0
#    else:
#        v = internal[-1] + 1
    
    
#    # find the node where to attach n
#    if x == 0:
#        # node is i
#        T[i][n] = d
#        if n not in T:
#            T[n] = {}
#        T[n][i] = d
##    else:
##        Found = False
##        for j in Dtrim[i]:
##            if Dtrim[i][j] == x:
##                # add n to j
##                T[j][n] = d
##                if n not in T:
##                    T[n] = {}
##                T[n][j] = d
##                Found = True
##                break
##        if Found == False:
##            v = n + 1
##            # create new node
##            T[v] = {n:d}
##            if n not in T:
##                T[n] = {}.add({v:d})
##            T[n][i] = d
#
#
#    else:
#        v = n + 1
#        # create new node
#        T[v] = {n:d}
#        if n not in T:
#            T[n] = {}
#        T[n][i] = d

    v = str(n + 1)
    print('new node', n)
    # create new node
    T[v] = {n:d}
    if n not in T:
        T[n] = {}
    T[n][i] = d
    T[i][v] = x
    T[v][i] = x
    T[v][k] = y
    T[k][v] = y
    
    
    
   

    
#    T[v] = {n:d}
    
    print('tree after attach', T)        
    return T


#Additive Phylogeny Problem
#Construct the simple tree fitting an additive matrix.
#
#
#Note on formatting: The adjacency list must have consecutive integer node labels starting from 0. The n leaves must be labeled 0, 1, ..., n-1 in order of their appearance in the distance matrix. Labels for internal nodes may be labeled in any order but must start from n and increase consecutively.


#Sample Dataset
#4
#0   13  21  22
#13  0   12  13
#21  12  0   13
#22  13  13  0
#Sample Output
#0->4:11
#1->4:2
#2->5:6
#3->5:7
#4->0:11
#4->1:2
#4->5:4
#5->4:4
#5->3:7
#5->2:6
    

L = [[0, 13, 21, 22],
[13, 0, 12, 13],
[21, 12, 0, 13],
[22, 13, 13, 0]]


M = ReadMatrix(L)
Tree = AdditivePhylogeny(M, 3)
print(Tree)