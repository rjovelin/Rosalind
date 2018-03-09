# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:10:19 2015

@author: Alivia
"""

def grab_sequences(input_file):
    '''
    (file) -> list
    Return a list with DNA sequences from inputf file
    '''
    # open file
    infile = open(input_file, 'r')
    # create list
    sequences= []
    for line in infile:
        line = line.rstrip()
        if line != '':
            sequences.append(line)
    # close file
    infile.close()
    return sequences



def make_debruijn_graph(sequences):
    '''
    (list) -> dict
    Return the debruijn graph of a list of kmers
    '''
    
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
    
    
    
def eulerian_path_circular_string(graph):
    '''
    (dict) -> dict
    Given a dictionnary of debruijn graph in the form of an adjacency list,
    return an Eulerian path in this graph.
    Precondition: the DNA string is cicular and the debruijn graph consists
    of a simple cycle with no repeated nodes
    '''
    
    # make a list of nodes already explored  
    explored = []
    
    # make a list to store the eulerian cycle
    cycle = []
     
    # makde lists of nodes
    nodes = [node for node in graph]
    
    
    # form an eulerian cycle starting anywhere    
    # start at first node = nodes[0]
    firstnode = nodes[0]
    secondnode = graph[firstnode][0]
    cycle.append(firstnode)
    cycle.append(secondnode)
    explored.append(secondnode)
    graph[firstnode].remove(secondnode)
    firstnode = secondnode
    
    # build the cycle until all secondnodes have been used
    while firstnode in graph and graph[firstnode] != []:
        secondnode = graph[firstnode][0]
        graph[firstnode].remove(graph[firstnode][0])
        cycle.append(secondnode)
        explored.append(secondnode)
        firstnode = secondnode
    
    return cycle    
    
def circular_string_reconstruction(input_file):
    '''
    (file) -> str
    Given a collection of all DNA k-mers taken from the same strand of a
    circular chromosome, such that their de Bruijn graph consists of
    exactly one simple cycle, return a cyclic superstring of minimal length
    containing the reads (thus corresponding to a candidate cyclic chromosome).
    '''
    
    # get sequences
    sequences = grab_sequences(input_file)
    
    # make a debruijn graph from kmers
    graph = make_debruijn_graph(sequences)
    
    # find the eulerian path in the debruijn graph
    cycle = eulerian_path_circular_string(graph)
    
    # build the DNA sequence from the kmer in the eulerian path
    DNA = ''
    DNA = cycle[0]
    for i in range(1, len(cycle)):
        DNA += cycle[i][-1]
    
    firstkmer = cycle[0]
    duplicate = firstkmer[:-1]
    stop = DNA.index(duplicate, 1)
    
    return DNA[:stop]
            
    
    
    