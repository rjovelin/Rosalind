# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 11:47:49 2015

@author: Richard
"""


def kmer_to_debruijn_graph(input_file):
    '''
    (file) -> dict
    Given an integer k and a collection of kmers in input file,
    return a DeBruijn graph in the form of an adjacent list
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    
    # get k
    k = int(infile.readline().rstrip())
    
    # make a list with the sequences
    sequences = []
    sequences = []
    for line in infile:
        line = line.rstrip()
        if line != '':
            sequences.append(line)
    # close file
    infile.close()
        
    # create a dictionnary to store the left and right k-1 mer
    LRkm1 = {}
    # loop over the k-mers
    for kmer in sequences:
        # find the left and right k-1 mers
        left = kmer[:-1]
        right = kmer[1:]
        # graph is directed from left to right k-1 mer
        if left in LRkm1:
            LRkm1[left].append(right)
        else:
            LRkm1[left] = [right]
    
    return LRkm1
    
    
def debruijn_to_eulerian_path(graph):
    '''
    (dict) -> dict
    Given a dictionnary of debruijn graph in the form of an adjacency list,
    return an Eulerian path in this graph.
    '''
    
    # count the number of edges
    edges = 0
    for node in graph:
        edges += len(graph[node])
    
       
    # make a list of nodes already explored  
    explored = []
    
    # make a list to store the eulerian cycle
    cycle = []
     
    # find starting node and end nodes
    first_nodes = [node for node in graph]
    second_nodes = []
    for node in graph:
        for item in graph[node]:
            second_nodes.append(item)
    
    for node in first_nodes:
        # starting node outdegree > indegree
        if len(graph[node]) > second_nodes.count(node):
            firstnode = node
    for node in second_nodes:
        # end node indegree < outdegree
        if node in graph:
            if len(graph[node]) < second_nodes.count(node):
                lastnode = node
        else:
            lastnode = node
         
    # create a link between the lastnode and the first node to form a 
    # graph that contains an eulerian cycle
    if lastnode in graph:
        graph[lastnode].append(firstnode)
    else:
        graph[lastnode] = [firstnode]
        
    # count the number of edges
    edges = 0
    for node in graph:
        edges += len(graph[node])
    
    # form an eulerian cycle starting at the firstnode    
    # add node to cycle
    cycle.append(firstnode)
    secondnode = graph[firstnode][0]
    cycle.append(secondnode)
    explored.append(secondnode)
    # remove secondnode from graph
    graph[firstnode].remove(secondnode)
    firstnode = secondnode
    
    # iterate until the cycle has the expected number of edges
    while len(cycle) < edges:
        
        # build the cycle until its stuck and return to the first node
        while firstnode in graph and graph[firstnode] != []:
            secondnode = graph[firstnode][0]
            graph[firstnode].remove(graph[firstnode][0])
            cycle.append(secondnode)
            explored.append(secondnode)
            firstnode = secondnode
        
        # read the cycle and find the node with nodes unexplored
        # start from that node and go through cycle, populating the new cycle
        # until that node is reached. add that node also to new cycle
        # resume looking for nodes in graph until its stuck
        
        not_found = True
        new_cycle = []
        for i in range(len(cycle)):
            if not_found:
                firstnode = cycle[i]
                for node in graph[firstnode]:
                    if not node in explored:
                        firstnode = cycle[i]
                        if cycle[0] == cycle[-1]:
                            secondnode = cycle[i:-1]
                        else:
                            secondnode = cycle[i:]
                        new_cycle.extend(secondnode)
                        new_cycle.extend(cycle[:i+1])
                        cycle = new_cycle
                        not_found = False
                        break
        
    
    # transform the eulerian cycle into an eulerian path from starting node
    # to end node
    euler_path = []
    end = cycle.index(lastnode)
    euler_path.extend(cycle[end+1:-1])
    euler_path.extend(cycle[:end+1])
    
    return euler_path
    
def string_reconstruction(input_file):
    '''
    (file) -> str
    
    Given an integer k and a collection of kmers in input file forming
    the kmer composition of a string sequence, return the string sequence
    (Return any sequence if multiple sequence exist)
    '''
    
    # construct the debruijn graph from the kmer composition
    graph = kmer_to_debruijn_graph(input_file)
    
    # find the eulerian path in the debruijn graph
    eulerian_path = debruijn_to_eulerian_path(graph)
    
    # build the DNA sequence from the kmer in the eulerian path
    DNA = ''
    DNA = eulerian_path[0]
    for i in range(1, len(eulerian_path)):
        DNA += eulerian_path[i][-1]
    
    return DNA

