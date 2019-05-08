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
    a path with all nodes between leafs i and k included, and the limb distance
    of leaf n to its parent node
    Return a list with node representing the path node-k where to attach a new node,
    the distance between i and node, and the distance between node and k
    '''
    
    
    for i in range(len(Distances)):
        nodes = [path[0], path[i+1]]
        if d > Distances[i]:
            # get the node where to attach n between node and k
            # get the distance between i and node and distance between node and k
            return [nodes[-1], d - Distances[i], Distances[-1] - d]
    

def AdditivePhylogeny(M, n):
    '''
    Given: n and a tab-delimited n x n additive matrix.
    Return: A weighted adjacency list for the simple tree fitting this matrix.
    '''
    
    print('n:', n)
    
    
    if n == 1:
        # return the tree consisting of a single edge of length D1,2
        return {0:{1:M[0][1]}, 1: {0:M[1][0]}}
    
    # compute distance between n and its parent node
    d = ComputeLimbLength(n, M)
    # compute Dbald in which n has bald limb length of 0
    # by substracting d to each element in row n and column n 
    M = ComputeDBald(M, n, d)
        
    #print('n,d:', n, d)
    #print(M)
    
    # find node i and k such that Di,k = Di,n + Dn,k    
    i, k = FindLeaves(M, n)
       
    #print('i, k:', [i, k])
    
    # get the distance between i and n, and between n and k
    x, y = M[i][n], M[n][k] 
    
    #print('distance {0}_{1}'.format(i, n), x)
    #print('distance {0}_{1}'.format(n,k), y)
    
    # compute Dtrim by removing row n and column n from M
    Dtrim = TrimMatrix(M, n)
            
    #$print(Dtrim)    
    
    # compute the tree corresonding to Dtrim
    T = AdditivePhylogeny(Dtrim, n - 1)
    
    #print('tree', T)
    
    # check if need to create new node
    if x != 0:
        
        #print('i, k for insertion', i, k)
        #print('distance to attach {0} from i'.format(n), x)
        
        # create new node
        v = str(n+2)
        #print('new node', v)
        #print('tree before insertion', T)
        # insert v between i and k at distance x from i
        # find the path between i and k to determine where to insert n
        path = BFS(T, i, k)
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
            
            path = BFS(T, i, k)
            # find where to attach v between i and k
            # comptute distance between i and each node in the path between i and k
            D = CumulativeDistanceDistanceBetweenLeaves(T, i, k)
            #print('=====')
            #print('trying to find attachement point')
            #print('n=', n)
            #print('path', path)
            #print('distances', D)
            #print('d', d)
            #print('=====')
            
            
            node, X, Y = FindAttachmentNode(D, path, x)
            
            print('=====')
            print('trying to find attachement point')
            print('n=', n)
            print('path', path)
            print('distances', D)
            print('d', d)
            print('node, distances', node, X, Y)
            print('=====')
            
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
        # attach n to i 
        T[i][n] = d
        T[n][i] = d
    
    
       
    print('tree after attach', T)        
    return T






def Format(Tree):
    L = list(Tree.keys())
    L.sort()
    for i in L:
        for j in Tree[i]:
            print(str(i) + '->' + str(j) + ':' + str(int(Tree[i][j]))) 









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
    

#L = [[0, 13, 21, 22],
#[13, 0, 12, 13],
#[21, 12, 0, 13],
#[22, 13, 13, 0]]


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
    Format(Tree)



