# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 22:59:49 2015

@author: Richard
"""

def trie_construction(input_file):
    '''
    (file) ->
    Return a trie on a collection of string patterns. 
    If Trie(Patterns) has n nodes, first label the root 
    with 1 and then label the remaining nodes with the 
    integers 2 through n in any order you like. Each edge
    of the adjacency list of Trie(Patterns) will be encoded
    by a triple: the first two members of the triple must be
    the integers labeling the initial and terminal nodes of the
    edge, respectively; the third member of the triple must be
    the symbol labeling the edge.
    '''
    
    
    
    
    
    # create a dictionnary to hold the trie
    trie = {}
    trie[1] = {}
    
    
    for seq in patterns:
        currentNode = 1
        for i in range(len(seq)):
            if seq[1] in trie[currentNode]:
                
                
            else:
                
        
    
    
    
    ####
    TRIECONSTRUCTION(Patterns)
        Trie ← a graph consisting of a single node root
        for each string Pattern in Patterns
            currentNode ← root
            for i ← 1 to |Pattern|
                if there is an outgoing edge from currentNode with label currentSymbol
                    currentNode ← ending node of this edge
                else
                    add a new node newNode to Trie
                    add a new edge from currentNode to newNode with label currentSymbol
                    currentNode ← newNode
        return Trie