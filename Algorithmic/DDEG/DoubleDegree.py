# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 01:38:36 2019

@author: Richard
"""

def DoubleDegree(PbFile):
    '''
    
    
    '''
    
    # extract data from file
    infile = open(PbFile)
    NumNodes, NumEdges = list(map(lambda x: int(x), infile.readline().rstrip().split()))
    # extract interactions between nodes
    interactions = infile.read().strip().split('\n')
    infile.close()
    # create a dict to store neighbors of node i D[i] = {j, k, m}
    D = {}
    # initialize D with empty sets to account for possible nodes without connection
    for i in range(1, NumNodes + 1):
        D[i] = set()
    # record interactions
    for i in interactions:
        i = i.split()
        D[int(i[0])].add(int(i[1]))
        D[int(i[1])].add(int(i[0]))
        
    # count interactions of each neighbor
    K = {}
    for i in D:
        K[i] = sum([len(D[j]) for j in D[i]])
    # get a sorted list of nodes
    nodes = list(K.keys())
    nodes.sort()
    count = [str(K[i]) for i in nodes]
    print(' '.join(count))
        
