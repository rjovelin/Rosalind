# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 17:12:41 2019

@author: rjovelin
"""

# The task is to use breadth-first search to compute single-source shortest distances in an unweighted directed graph.

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
            i = int(line[0])
            j = int(line[1])
            if i not in D:
                D[i] = set()
            D[i].add(j)
    infile.close()
    return D


# implement breadth-first-search for finding the path between 2 nodes
#def BFS(tree, f, goal):
#    '''
#    (dict, str, str) -> list
#    Take a dictionary representing with interactions among node, and return
#    a list of nodes representing the path betwen nodes f and goal
#    '''
#    # create a queue to add the discovered nodes that are not yet visited
#    # keep track of the path visited to reach discovered node
#    queue = [(f, [f])]
#    visited = []
#    while len(queue) != 0:
#        # get a node to visit. keep track of the path up to node
#        (node, path) = queue.pop(0)
#        
#        
#        if node not in visited:
#            
#            visited.append(node)
#            # add node children to queue
#            for child in tree[node]:
#                
#                if child == goal:
#                    path.append(child)
#                    return path
#                else:
#                    queue.append((child, path + [child]))
        

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
        if node in tree:
            if node not in visited:
                visited.append(node)
                # add node children to queue
                for child in tree[node]:
                    if child == goal:
                        path.append(child)
                        return len(path) -1
                    else:
                        queue.append((child, path + [child]))
        else:
            if node == goal:
                return -1










   
  
#def PrintArray(PBFile):
#    '''
#    (file) -> None
#    A simple directed graph with n≤10**3 vertices in the edge list format
#    Print an array D[1..n] where D[i] is the length of a shortest path from
#    the vertex 1 to the vertex i (D[1]=0). If i is not reachable from 1 set D[i] to −1.
#    '''    
#    
#    # make a directed tree
#    tree = ConstructTree(PBFile)
#    #
#    
#    # collect all the nodes
#    nodes = []
#    for i in tree:
#        nodes.append(i)
#        for j in tree[i]:
#            nodes.append(j)
#    # remove duplicate nodes
#    nodes = list(set(nodes))
#    nodes.sort()
#        
#    # loop over nodes in tree
#    Distances = []
#    for i in nodes:
#        if i == 1:
#            # set distance to self to 0
#            d = 0
#        else:
#            try:
#                # find short path between nodes 1 and i
#                # compute distance = len(path) -1 
#                # all edges have same weight = 1
#                path = BFS(tree, 1, i)
#                d = len(path) -1
#            except:
#                # set d = -1 if node is not reachable
#                d = -1
#        # collect distances
#        Distances.append(str(d))
#    print(' '.join(Distances))    
    
    
    
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
            d = BFS(tree, 1, i)
            if d == None:
                d = -1
        # collect distances
        Distances.append(str(d))
    print(' '.join(Distances))        
    
    
