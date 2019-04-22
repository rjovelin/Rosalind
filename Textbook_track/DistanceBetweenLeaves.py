# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:51:32 2019

@author: rjovelin
"""




def ComputeDistance(PbFile):
    '''
    (file) -> dict
    Take a file with An integer n followed by the adjacency list
    of a weighted tree with n leaves and return a dictionary with interactions
    among nodes
    '''
    
    infile = open(PbFile)
    NLeaves = int(infile.readline().rstrip())
    Distance = {}
    for line in infile:
        if ':' in line:
            line = line.rstrip()
            i = int(line[:line.index('-')])
            j = int(line[line.index('>')+1: line.index(':')])
            distance = int(line[line.index(':')+1:])
            if i not in Distance:
                Distance[i] = {}
            if j not in Distance:
                Distance[j] = {}
            Distance[i][j] = distance
            Distance[j][i] = distance
    infile.close()
    return Distance

    
    


def ConstructTree(PbFile):
    '''
    (file) -> dict
    file) -> dict
    Take a file with An integer n followed by the adjacency list
    of a weighted tree with n leaves and return a dictionary with interactions
    among nodes
    '''

    infile = open(PbFile)
    NLeaves = int(infile.readline().rstrip())
    D = {}
    for line in infile:
        if ':' in line:
            line = line.rstrip()
            i = int(line[:line.index('-')])
            j = int(line[line.index('>')+1: line.index(':')])
            if i not in D:
                D[i] = set()
            if j not in D:
                D[j] = set()
            D[i].add(j)
            D[j].add(i)
    infile.close()
    return D


def FindLeaves(tree):
    '''
    (dict) -> list
    Take a dict with interactions amond nodes and return a list
    of leaf nodes (interaction = 1)
    '''
    
    Leaves = [i for i in tree if len(tree[i]) == 1]
    return Leaves



    
# implement breadth-first-search for traversing the tree from f to goal
def BS(tree, f, goal, leaves):
    '''
    
    
    '''
    # create a queue to add the discovered nodes that are not yet visited
    queue = []
    queue.append(f)
    visited = []
    path = []
    while len(queue) != 0:
        # get a node to visit
        node = queue.pop(0)
        if node not in visited:
            path.append(node)
            visited.append(node)
            # add node children to queue
            if node in tree:
               for children in tree[node]:
                   if children == goal:
                       queue.append(children)
                       break
                   else:
                       if children not in leaves:
                           queue.append(children)
    return path        
    
def ComputeDistanceBetweenLeaves(tree, f, goal, leaves, distance):

    path = BS(tree, f, goal, leaves)
    d = 0
    for i in range(len(path)-1):
        d = d + distance[path[i]][path[i+1]]
    print(d)
    
    
