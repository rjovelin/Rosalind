# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 11:12:02 2015

@author: Richard
"""



def debruijn_graph_kmer(input_file, outputfile):
    '''
    (file) -> file
    Given an collection of DNA sequence in input file,
    return a DeBruijn graph in the form of an adjacent list
    Note: DeBruijn graph of a collection of k-mers is a directed graph
    connecting k-1 mers for which the 2 k-1 mers are respectively prefix and
    suffix of a k-mer. The edges are k-2 overlap between the k-1 mers. but are
    also the k-mers
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
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
    
    # open file for writing
    newfile = open(outputfile, 'w')
       
    # write graph to outputfile
    for node in LRkm1:
        if len(LRkm1[node]) == 1:
            newfile.write(node + ' -> ' + LRkm1[node][0] + '\n')
        else:
            newfile.write(node + ' -> ')
            for seq in LRkm1[node][:-1]:
                newfile.write(seq + ', ')
            newfile.write(LRkm1[node][-1] + '\n')
    # close file after writing
    newfile.close()
    
    
   
    
    
    