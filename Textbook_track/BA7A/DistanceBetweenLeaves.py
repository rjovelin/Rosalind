# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:51:32 2019

@author: rjovelin
"""


def ExtractDistance(PbFile):
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
        
   
def ComputeDistanceBetweenLeaves(tree, f, goal, distance):
    '''
    (dict, str, str, list, dict) -> int
    Return the distance between leaves f and goal by traversing tree
    '''
    
    path = BFS(tree, f, goal)
    d = 0
    if f == goal:
        d = 0
    else:
        for i in range(len(path)-1):
            d = d + distance[path[i]][path[i+1]]
    return d
    
    
def PrintMatrix(PBFile):
    '''
    (file) -> None
    Take a file with an integer n followed by the adjacency list of a weighted tree with
    n leaves, and print a n x n matrix (di,j), where di,j is the length of the path
    between leaves i and j.
    '''    
    
    distance = ExtractDistance(PBFile)
    tree = ConstructTree(PBFile)
    leaves = FindLeaves(tree)
    leaves.sort()
    matrix = []
    for i in range(len(leaves)):
        line = []
        for j in range(len(leaves)):
            line.append(str(ComputeDistanceBetweenLeaves(tree, leaves[i], leaves[j], distance)))
        matrix.append(' '.join(line))
    for i in matrix:
        print(i)
    