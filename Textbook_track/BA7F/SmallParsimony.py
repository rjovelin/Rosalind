# -*- coding: utf-8 -*-
"""
Created on Fri May 24 17:10:58 2019

@author: rjovelin
"""

import math


#Small Parsimony Problem. Find the most parsimonious labeling of the internal nodes of a rooted tree. 
#
#Input: A rooted binary tree with each leaf labeled by a string of length m.
#Output: A labeling of all other nodes of the tree by strings of length m that minimizes the parsimony score of the tree.
#


def IsRoot(T, node):
    '''
    (dict, int) -> bool    
    Return True if node is root
    '''
    
    if len(T[node]) == 2:
        return True
    else:
        return False


def IsLeaf(T, node):
    '''
    (dict, int) -> bool
    Return True if node is a leaf
    '''
    if len(T[node]) == 1:
        return True
    else:
        return False


def CountRipeNodes(T, Tag):
    '''
    (dict, dict) -> int
    Count the ripe nodes in T. An internal node of T ripe if its tag is 0 but
    its children’s tags are both 1
    '''
    
    C = 0
    
    for i in T:
        if IsLeaf(T, i) == False:
            if Tag[i] == 0:
                childrenTags = [Tag[j] for j in T[i] if j < i]
                if len(list(set(childrenTags))) == 1 and list(set(childrenTags))[0] == 1:
                    # found a ripe node
                    C += 1
    return C            

            
def GetRipeNode(T, Tag):
    '''
    (dict, dict) -> int
    Return a ripe node v of T. An internal node of T ripe if its tag is 0 but
    its children’s tags are both 1
    '''
    
    for i in T:
        # check if internal node
        if IsLeaf(T,i) == False and Tag[i] == 0:
            childrenTags = [Tag[j] for j in T[i] if j < i]
            if len(list(set(childrenTags))) == 1 and list(set(childrenTags))[0] == 1:
                return i            

def HammingDistance(seq1, seq2):
    '''
    (str, str) -> int
    Return the hamming distances between sequences 1 and 2
    Pre-condition: seq1 and se2 have same length
    '''
    
    d = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            d += 1
    return d



def ComputeMinimumScore(Score, v, k):
    '''
    (dict, int, str) -> tuple
    Return a tuple with the minimum score of v among all characters i 
    when parent node of v is labaled with character k, and the corresponding character i
    '''
    
    S = {}
    for i in ['A', 'C', 'G', 'T']:
        d = Score[v][i] + HammingDistance(i, k)
        if d in S:
            S[d].append(i)
        else:
            S[d] = [i]
    # find minimum score
    m = min(list(S.keys()))
    return m, S[m][0] 


def FindCharacterMinScore(Score, node, Minimum):
    '''
    (dict, int, int) -> str
    Return the nucleotide at node v with Minimum Score
    '''
    
    for i in Score[node]:
        if Score[node][i] == Minimum:
            return i
    

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


def SmallParsimony(T, Character):
    '''
    (dict, dict) -> int
    Return the parsimony score for a binary rooted tree T whose leaves
    are labeled by symbols stored in Character (i.e. Character(v) is the label of leaf v).
    '''
    
    # use a dictionary to keep track of processed nodes
    # i.e., Tag(v) = 1 if the array sk(v) has been computed and Tag(v) = 0 otherwise
    Tag = {}
    
    # keep track of the parsimony score of each node
    Score = {}
    
    # for each node v in tree T
    for v in T:
        # Tag(v) ← 0
        Tag[v] = 0
        # if v is a leaf
            # Tag(v) ← 1
        if IsLeaf(T, v) == True:
            Tag[v] = 1
            if v not in Score:
                Score[v] = {}
            #  for each symbol k in the alphabet
            for k in ['A', 'C', 'G', 'T']:
                # if Character(v) = k --> sk(v) ← 0 else: sk(v) ← ∞
                if Character[v] == k:
                    Score[v][k] = 0
                else:
                    Score[v][k] = math.inf

    # while there exist ripe nodes in T
    while CountRipeNodes(T, Tag) > 0:
        # v ← a ripe node in T
        v = GetRipeNode(T, Tag)
        # Tag(v) ← 1
        Tag[v] = 1
        # identify son and duaghter of node v
        children = [i for i in T[v] if i < v]
        son, daughter = children[0], children[1]
        if v not in Score:
            Score[v] = {}
        #  for each symbol k in the alphabet
        for k in ['A', 'C', 'G', 'T']:
            # sk(v) ← minimumall symbols i {si(Daughter(v))+αi,k} + minimumall symbols j {sj(Son(v))+αj,k}
            Score[v][k] = ComputeMinimumScore(Score, son, k)[0] + ComputeMinimumScore(Score, daughter, k)[0]
    assert IsRoot(T, v) == True 
    
    # return minimum over all symbols k {sk(v)}
    # get the score and nucleotide at the root
    Minimum = min([Score[v][k] for k in Score[v]])
    nucleotide = FindCharacterMinScore(Score, v, Minimum)
    # add root nucleotide to character
    Character[v] = nucleotide
    
    return Minimum, v, Score, Character
    
    
def InferAncestralCharacter(T, Minimum, Root, Score, Character):
    '''
    (dict, int, int, dict, dict) -> dict
    Take a Tree, the minimum parsimony score of the Tree, the root label,
    a dictionary with Scores at each node and a dictionary with character for leaf nodes
    amd return a dictionary with character at all nodes minimizing the Tree score for that character

    '''

    C = {}
    C[Root] = Character[Root]
    # use breadth-first-search to find the paths between the root and each leaf
    Leafs = [i for i in T if IsLeaf(T,i) == True] 
    for i in Leafs:
        C[i] = Character[i]
        path = BFS(T, Root, i)
        for j in range(len(path)):
            # find the minimum score of node
            Minimum = min([Score[path[j]][k] for k in Score[path[j]]])
            # list all nucleotides with minimum score for node
            nucleotides = [n for n in Score[path[j]] if Score[path[j]][n] == Minimum]
            if path[j] not in C:
                # check if parent nucleotide in list
                if C[path[j-1]] in nucleotides:
                    # use parent nucleotide
                    C[path[j]] = C[path[j-1]]
                else:
                    # take first nucleotide in list
                    C[path[j]] = nucleotides[0]
    return C



def ReadTreeFromFile(PbFile):
    '''
    (file) -> dict    
    Return a dictionary representing the Tree written as adjacency list if file 
    '''
    T = {}
    i = 0
    infile = open(PbFile)
    infile.readline()
    for line in infile:
        if '-' in line:
            line = line.rstrip().split('->')
            if line[1].isalpha() == True:
                if int(line[0]) not in T:
                    T[int(line[0])] = set()
                T[int(line[0])].add(i)
                if i not in T:
                    T[i] = set()
                T[i].add(int(line[0]))
                i += 1
            else:
                if int(line[0]) not in T:
                    T[int(line[0])] = set()
                T[int(line[0])].add(int(line[1]))
                if int(line[1]) not in T:
                    T[int(line[1])] = set()
                T[int(line[1])].add(int(line[0]))
    infile.close()
    return T            
    

def GetCharacter(PbFile):
    '''
    (file) -> dict
    Return a dictionary with the leaf label: leaf sequence from the adjacency 
    in file
    '''
    T = {}
    i = 0
    infile = open(PbFile)
    infile.readline()
    for line in infile:
        if '-' in line:
            line = line.rstrip().split('->')
            if line[1].isalpha() == True:
                T[i] = line[1]
                i += 1
    infile.close()
    return T            


def SolvePb(PbFile):
    '''
    (file) -> None
    Print the minimum score of the Tree in file and an adjacency list with sequences 
    at internal and leaf nodes minimizing the core of the tree
    '''
    
    # get the leaf sequences 
    S = GetCharacter(PbFile)
    # get the tree
    T = ReadTreeFromFile(PbFile)
    # make a list of dictionaries with individual leaf nucleotides
    L = []
    for i in range(len(S[0])):
        d = {}
        for j in S:
            d[j] = S[j][i]
        L.append(d)
    # apply small parsimony to get the minimum tree score and the sequences
    # minimizing the score of the tree
    D, Seqs = [], []
    for i in range(len(L)):
        Minimum, Root, Score, L[i] = SmallParsimony(T, L[i])
        D.append(Minimum) 
        Seqs.append(InferAncestralCharacter(T, Minimum, Root, Score, L[i]))

    # built sequences by adding each character at each node to first character
    for i in range(1, len(Seqs)):
        for node in Seqs[i]:
            Seqs[0][node] = Seqs[0][node] + Seqs[i][node]
    
    print(sum(D))
    
    for i in T:
        for j in T[i]:
            print(Seqs[0][i] + '->' + Seqs[0][j] + ':' + str(HammingDistance(Seqs[0][i], Seqs[0][j])))
    
    newfile = open('SmallParsimonySolution.txt', 'w')
    newfile.write(str(sum(D)) + '\n')
    for i in T:
        for j in T[i]:
            newfile.write(Seqs[0][i] + '->' + Seqs[0][j] + ':' + str(HammingDistance(Seqs[0][i], Seqs[0][j])) + '\n')
    newfile.close()
    
    
