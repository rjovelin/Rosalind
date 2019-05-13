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


def ComputeDBald(M, n, d):
    '''
    (dict, int, int) -> dict
    Return a modified distance matrix M in which the limb distance d of n is
    substracted from the distance of all leafs
    '''
    
    for i in M[n]:
        if i != n:
            M[n][i] = M[n][i] - d
    for i in M:
        if i != n:
            M[i][n] = M[i][n] -d
    return M


def FindLeaves(M, n):
    '''
    (list) -> list
    Find leafs i and k such that Di,k = Di,n + Dn,k
    '''
    trio = []
    for i in M:
        for k in M:
            if M[i][k] == M[i][n] + M[n][k]:
                trio.append([i, k])
    return trio[0]
    

def TrimMatrix(M, n):
    '''
    (dict, int) -> dict
    Return a modified distance matrix M trimmed of n in columns and rows
    '''
    
    del M[n]
    for i in M:
        del M[i][n]
    return M

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


def CumulativeDistanceDistanceBetweenLeaves(tree, f, goal):
    '''
    (dict, int, int) -> list
    Return a list of distances from leaf f to all nodes on the path between f
    and goal in Tree. 
    Precondition: there is at least 1 internal node between f and goal
    '''
    
    path = BFS(tree, f, goal)
    D = []
    d = 0
    if f == goal:
        D = [d]
    else:
        for i in range(len(path)-1):
            d = d + tree[path[i]][path[i+1]]
            D.append(d)
    return D


def FindAttachmentNode(Distances, path, d):
    '''
    (list, list, int) -> list
    Take the list of Distances between i and all nodes on the path between i and k,
    a path with all nodes between leafs i and k included, and the distance between i and n
    Return a list with node representing the path node-k where to attach a new node,
    the distance between i and node, and the distance between node and k
    '''
    
    nodes = [[path[0], path[i]] for i in range(1, len(path))]
    
    i = 0
    
    if d < Distances[i]:
        return [nodes[0][0], d, Distances[0]-d]
    else:
        while Distances[i] < d:
            i += 1
        
        # get the node where to attach n between node and k
        # get the distance between i and node and distance between node and k
        return [nodes[i-1][-1], d - Distances[i-1], Distances[i] - d]

def AdditivePhylogeny(M, n):
    '''
    Given: n and a tab-delimited n x n additive matrix.
    Return: A weighted adjacency list for the simple tree fitting this matrix.
    '''
    
    if n == 1:
        # return the tree consisting of a single edge of length D1,2
        return {0:{1:M[0][1]}, 1: {0:M[1][0]}}
    
    # compute distance between n and its parent node
    d = ComputeLimbLength(n, M)
    # compute Dbald in which n has bald limb length of 0
    # by substracting d to each element in row n and column n 
    M = ComputeDBald(M, n, d)
        
    # find node i and k such that Di,k = Di,n + Dn,k    
    i, k = FindLeaves(M, n)
       
    # get the distance between i and n, and between n and k
    x, y = M[i][n], M[n][k] 
    
    # compute Dtrim by removing row n and column n from M
    Dtrim = TrimMatrix(M, n)
            
    # compute the tree corresonding to Dtrim
    T = AdditivePhylogeny(Dtrim, n - 1)
    
    # check if need to create new node
    if x != 0:
        v = str(n)
        # insert v between i and k at distance x from i
        # insert v between i and k if no node in between
        if k in T[i]:
            del T[i][k]
            del T[k][i]
            T[v] = {n:d}
            if n not in T:
                T[n] = {}
            T[n][v] = d
            T[v][i] = x
            T[i][v] = x
            T[v][k] = y
            T[k][v] = y
        else:
            # find the path between i and k to determine where to insert n
            path = BFS(T, i, k)
            # find where to attach v between i and k
            # comptute distance between i and each node in the path between i and k
            D = CumulativeDistanceDistanceBetweenLeaves(T, i, k)
            node, X, Y = FindAttachmentNode(D, path, x)
            # check if there are other nodes between node and k
            if node == path[-2]:
                # insert between node and k
                del T[node][k]
                del T[k][node]
                T[v] = {n:d}
                if n not in T:
                    T[n] = {}
                T[n][v] = d
                T[v][k] = Y
                T[k][v] = Y
                T[v][node] = X
                T[node][v] = X
            else:
                # insert between node and node immediately after
                # remove link between node and node immediately after in the path
                pos = path.index(node)
                node2 = path[pos+1]
                del T[node][node2]
                del T[node2][node]
                T[v] = {n:d}
                if not n in T:
                    T[n] = {}
                T[n][v] = d
                T[v][node] = X
                T[node][v] = X
                T[v][node2] = Y
                T[node2][v] = Y
    else:
        # attach n to i 
        T[i][n] = d
        T[n][i] = d
    
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
            print(str(i) + '->' + str(j) + ':' + str(int(Tree[i][j]))) 


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


def SolvePb(PbFile):
    # get matrix from file
    L = GetMatrixFromFile(PbFile)
    # convert list to dictionary
    M = ReadMatrix(L)
    # find the last node in matrix to start with
    leafs = list(M.keys())
    leafs.sort()
    n = leafs[-1]
    # use additive phylogeny to reconstruct the tree
    Tree = AdditivePhylogeny(M, n)
    #Format(Tree)
    Format(Tree)


