# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 17:12:41 2019

@author: rjovelin
"""

# The task is to use breadth-first search to compute single-source shortest distances in an unweighted directed graph.

import math


#### 1st approach, identify paths

def ConstructTree(PbFile):
    '''
    (file) -> dict
    A simple directed graph with n≤10**3 vertices in the edge list format
    and return a dictionary with directed connections
    '''

    infile = open(PbFile)
    infile.readline()
    D = {}
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            i = int(line[0].rstrip())
            j = int(line[1].rstrip())
            if i not in D:
                D[i] = set()
            D[i].add(j)
    infile.close()
    
    return D


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
            if node in tree:
                for child in tree[node]:
                    if child == goal:
                        path.append(child)
                        #print('path', path)
                        return path
                    else:
                        queue.append((child, path + [child]))
                        #print(f, goal, node, path + [child])

 
def PrintArray(PBFile):
    '''
    (file) -> None
    A simple directed graph with n≤10**3 vertices in the edge list format
    Print an array D[1..n] where D[i] is the length of a shortest path from
    the vertex 1 to the vertex i (D[1]=0). If i is not reachable from 1 set D[i] to −1.
    '''    
    
    # make a directed tree
    tree = ConstructTree(PBFile)
    #
    
    # collect all the nodes
    nodes = []
    for i in tree:
        nodes.append(i)
        for j in tree[i]:
            nodes.append(j)
    # remove duplicate nodes
    nodes = list(set(nodes))
    nodes.sort()
        
    # loop over nodes in tree
    Distances = []
    for i in nodes:
        if i == 1:
            # set distance to self to 0
            d = 0
        else:
            path = BFS(tree, 1, i)
            if path == None:
                d = -1
            else:
                d = len(path) -1
        # collect distances
        Distances.append(str(d))
    print(' '.join(Distances))        
    
    
### 2nd approach. use BFS to compute distance on the fly
    
    
def Breadth_First_Search_Dist(PbFile):
    
    s = 1
    
    infile = open(PbFile)
    infile.readline()
    E = infile.read().strip().split('\n')
    for i in range(len(E)):
        E[i] = list(map(lambda x: int(x), E[i].split()))
    
    V = []
    for i in E:
        V.extend(i)
    V = list(set(V))
    V.sort()
    
    infile.close()
    
    dist = {}
    for i in V:
        dist[str(s)+',' + str(i)] = math.inf
    
    dist[str(s)+','+str(s)] = 0
    Q = [s]
    while len(Q) != 0:
        u = Q.pop(0)
        for e in E:
            if e[0] == u:
                v = e[1]
                if dist[str(s)+','+str(v)] == math.inf:
                    Q.append(v)
                    dist[str(s)+','+str(v)] = dist[str(s)+','+str(u)] + 1
    D = {}
    for i in dist:
        j = i.split(',')
        a = int(j[1])
        D[a] = dist[i]
    L = [i for i in D]
    L.sort()
    K = []
    for i in L:
        if D[i] == math.inf:
            K.append(-1)
        else:
            K.append(D[i])
    K = ' '.join(list(map(lambda x: str(x), K)))
    print(K)