# -*- coding: utf-8 -*-
"""
Created on Tue May 28 21:39:00 2019

@author: rjovelin
"""

# Solve the Nearest Neighbors of a Tree Problem.
# Input: Two internal nodes a and b specifying an edge e,
# followed by an adjacency list of an unrooted binary tree.
# Output: Two adjacency lists representing the nearest neighbors of the tree
# with respect to e. Separate the adjacency lists with a blank line.



import copy



def ReadTreeFromFile(PbFile):
    '''
    (file) -> tuple    
    Return a tuple with a dictionary representing the Tree written as adjacency
    list in file, and nodes a and b in Tree  
    '''
    
    
    infile = open(PbFile)
    # get the nodes specifying the edge where to do nearest neighbors interchange
    a, b = infile.readline().rstrip().split()
    a, b = int(a), int(b)
    # convert adjency list into dict
    T = {}
    for line in infile:
        if '-' in line:
            line = line.rstrip().split('->')
            node1, node2 = int(line[0]), int(line[1])
            if node1 not in T:
                T[node1] = set()
            if node2 not in T:
                T[node2] = set()
            T[node1].add(node2)
            T[node2].add(node1)
    infile.close()
    return (T, a, b)


def NearestNeighborsInterchange(T, a, b):
    '''
    (dict, int, int) -> tuple

    '''

    # make a list of children for a and b, exclusing a and b
    x, y = [i for i in T[a] if i != b]
    v, w = [i for i in T[b] if i != a]
   
    # make a copy of T
    D = copy.deepcopy(T)
    K = copy.deepcopy(T)
   
    # remove link between a and x
    D[a].remove(x)
    D[x].remove(a)
    # remove link between b and v
    D[b].remove(v)
    D[v].remove(b)
    # add link between a and v
    D[a].add(v)
    D[v].add(a)
    # add link between b and x
    D[b].add(x)
    D[x].add(b)

    # remove link between a and x
    K[a].remove(x)
    K[x].remove(a)
    # remove link between b and w
    K[b].remove(w)
    K[w].remove(b)
    # add link between a and w
    K[a].add(w)
    K[w].add(a)
    # add link between b and x
    K[b].add(x)
    K[x].add(b)

    return D, K




def FormatTree(T):
    '''
    
    
    '''
    
    
    for i in T:
        for j in T[i]:
            print(str(i) + '->' + str(j))
            
    


def SolvePb(PbFile):
    '''
    
    
    
    '''
    
    # get Tree from file
    T, a, b = ReadTreeFromFile(PbFile)
    # perform nearest neighbors interchange
    D, K = NearestNeighborsInterchange(T, a, b)
    # print Trees
    FormatTree(D)
    print('\n')
    FormatTree(K)



