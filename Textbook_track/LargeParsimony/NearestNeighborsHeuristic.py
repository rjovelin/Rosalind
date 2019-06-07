# -*- coding: utf-8 -*-
"""
Created on Wed May 29 19:21:07 2019

@author: rjovelin
"""


import math
import copy



def FindTwoConnectedInternalNodes(T):
    '''
    (dict) -> tuple
    Return a tuple of connected internal nodes in tree T
    '''
    
    for i in T:
        if IsRoot(T, i) ==False and IsLeaf(T, i) == False:
            # make a list of children nodes
            L = [j for j in T[i]]
            # make a list of boolean for chidren nodes that are not root and not leafs
            M = []
            for j in L:
                if IsRoot(T, j) ==False and IsLeaf(T, j) == False:
                    M.append(True)
                else:
                    M.append(False)
            # check if at least 1 child node is not root and not leaf
            if True in M:
                k = M.index(True)
                return (i, L[k])


def ListAllInternalEdges(T):
    '''
    (dict) -> list
    Return a list with all internal nodes in tree T    
    '''
    I = []
    
    # make a list of internal nodes
    L = [i for i in T if IsRoot(T, i) ==False and IsLeaf(T, i) == False]
    for i in range(len(L)-1):
        for j in range(i+1, len(L)):
            path = BFS(T, L[i], L[j])
            if len(path) == 2:
                I.append(tuple(path))
            else:
                for k in range(len(path)-1):
                    I.append((path[k], path[k+1]))
    I = list(set(I))
    return I
    

def RootTree(T):
    '''
    (dict) -> dict
    Return a modified dictionary T with a root between 2 connected internal node
    '''
    
    # find 2 connected internal nodes 
    node1, node2 = FindTwoConnectedInternalNodes(T)
    # make a list with all nodes in T
    nodes = list(T.keys())
    # label root node
    Root = max(nodes) + 1
    # add root to tree between node 1 and 2
    T[Root] = {node1, node2}
    # add root ro node1 and node 2
    T[node1].add(Root)
    T[node2].add(Root)
    # remove link between node1 and 2
    T[node1].remove(node2)
    T[node2].remove(node1)
    
    return T

def RemoveRoot(T):
    '''
    (dict) -> dict
    Return an unrooted tree T by removing the root from T
    '''
    Root = ''
    # find the root of the tree
    for i in T:
        if IsRoot(T, i) == True:
            Root = i
            break
    assert Root != ''
    # make a list of children nodes of root
    L = list(T[Root])
    # add link bewteen children nodes    
    T[L[0]].add(L[1])
    T[L[1]].add(L[0])
    # remove link between chidren and root
    for i in L:
        T[i].remove(Root)
    # remove root from tree
    del T[Root]
    return T

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
    # identify the root of the Tree
    Root = ''
    for i in T:
        if IsRoot(T, i) == True:
            Root = i
            break
    assert Root != ''
    for i in T:
        # check if internal node
        if IsLeaf(T,i) == False and Tag[i] == 0:
            if i == Root:
                children = [j for j in T[i]]
            else:
                # get the path from the root to node i
                path = BFS(T, Root, i)
                parent = path[-2]
                # identify son and duaghter of node v
                children = [j for j in T[i] if j != parent]
            assert len(children) == 2
            childrenTags = [Tag[j] for j in children]
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
    
    Root = ''
    # identify the root of the Tree
    for i in T:
        if IsRoot(T, i) == True:
            Root = i
            break
    assert Root != ''
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
    while 0 in Tag.values():
        # v ← a ripe node in T
        v = GetRipeNode(T, Tag)
        # Tag(v) ← 1
        Tag[v] = 1
        if v == Root:
            children = [i for i in T[v]]
        else:
            # get the path from the rot to node v
            path = BFS(T, Root, v)
            parent = path[-2]
            # identify son and duaghter of node v
            children = [i for i in T[v] if i != parent]
        
        assert len(children) == 2
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
    
    Leafs = {}
    
    T = {}
    i = 0
    infile = open(PbFile)
    infile.readline()
    for line in infile:
        if '-' in line:
            line = line.rstrip().split('->')
            if line[0].isalpha() ==True:
                if line[0] not in Leafs:
                    Leafs[line[0]] = i
                    i += 1
                if Leafs[line[0]] not in T:
                    T[Leafs[line[0]]] = set()
                T[Leafs[line[0]]].add(int(line[1]))
            else:
                if line[1].isalpha() == True:
                    if line[1] not in Leafs:
                        Leafs[line[1]] = i
                        i += 1
                    if int(line[0]) not in T:
                        T[int(line[0])] = set()
                    T[int(line[0])].add(Leafs[line[1]])
                else:
                    if int(line[0]) not in T:
                        T[int(line[0])] = set()
                    T[int(line[0])].add(int(line[1]))
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


def InferFullSequences(T, Sequences):
    '''
    (file) -> None
    Print the minimum score of the Tree in file and an adjacency list with sequences 
    at internal and leaf nodes minimizing the core of the tree
    '''
    
    # root the tree
    T = RootTree(T)
    
    # make a list of dictionaries with individual leaf nucleotides
    L = []
    for i in range(len(Sequences[0])):
        d = {}
        for j in Sequences:
            d[j] = Sequences[j][i]
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
    
    # remove root from Tree
    T = RemoveRoot(T)
    
    # return Tree score and sequences at all nodes removing the root     
    toremove = [i for i in Seqs[0] if i not in T]
    for i in toremove:
        del Seqs[0][i]
    
    return sum(D), Seqs[0]
    
    

def NearestNeighborsInterchange(T, a, b):
    '''
    (dict, int, int) -> tuple
    Return a tuple with dictionaries representing nearest neighbors
    of tree T after performing nearest neighbors interchange
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

def FormatOutput(S, T, Seqs):
    '''
    (int, dict, dict) -> None
    Print required output for the grader
    '''
    
    print(S)
    for i in T:
        for j in T[i]:
            print(Seqs[i] + '->' + Seqs[j] + ':' + str(HammingDistance(Seqs[i], Seqs[j])))
    print('\n')


def LargeParsimonyHeuristic(Tree, LeafSeqs):
    '''
    (dict, dict) -> dict
    Take a tree and a dictionary with leaf sequences and return a Tree with
    mininmum parsimony score over subtrees of T through a heuristic
    
    '''
    
    # score ← ∞
    score = math.inf
    # generate an arbitrary unrooted binary tree Tree with |Strings| leaves
    
    # label the leaves of Tree by arbitrary strings from Strings
    
    
    # solve  the  Small Parsimony in an Unrooted Tree Problem for Tree
    # label the internal nodes of Tree according to a most parsimonious labeling
    # newScore ← the parsimony score of Tree
    NewScore, Sequences = InferFullSequences(Tree, LeafSeqs)
        
    #newTree ← Tree
    NewTree = Tree
    # get the sequences at all nodes for Tree
    NewLabels = Sequences
    
    while NewScore < score:
        # score ← newScore
        score = NewScore
        # Tree ← newTree
        Tree = NewTree
        for i in Tree:
            assert IsRoot(Tree, i) == False
        # for each internal edge e in Tree
        for a, b in ListAllInternalEdges(Tree):
            # make a list of nearest-neighbors
            Neighbors = NearestNeighborsInterchange(Tree, a, b)
            # for each nearest neighbor NeighborTree of Tree with respect to the edge e
            for T in Neighbors:
                # solve the Small Parsimony in an Unrooted Tree Problem for NeighborTree
                # label the internal nodes of Tree according to a most parsimonious labeling
                # neighborScore ← the minimum parsimony score of NeighborTree
                NeighborScore, Sequences = InferFullSequences(T, LeafSeqs)
                # if neighborScore < newScore: newScore ← neighborScore, newTree ← NeighborTree
                if NeighborScore < NewScore:
                    NewLabels = Sequences
                    NewScore = NeighborScore
                    NewTree = T 
        # print output if NewScore , score
        if NewScore < score:
            FormatOutput(NewScore, NewTree, NewLabels)             
    return NewTree


def SolveProblem(PbFile):
    
    # get the leaf sequences 
    C = GetCharacter(PbFile)
    # get the tree
    T = ReadTreeFromFile(PbFile)
    # infer Tree using heuristic search
    Tree = LargeParsimonyHeuristic(T, C)
    
    return Tree
    