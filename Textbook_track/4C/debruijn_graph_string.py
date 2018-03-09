# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 11:12:02 2015

@author: Richard
"""

def debruijn_graph_from_string(input_file, outputfile):
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
    # grab k and sequence S
    k = int(infile.readline().rstrip())
    S = infile.readline().rstrip()
    # close file
    infile.close()
    
    # make a list with the k-mers from S
    sequences = []
    # loop over S and find the k-mers, and add to sequence list
    for i in range(len(S) -k +1):
        kmer = S[i:i+k]
        sequences.append(kmer)
    
    
    # create a dict to store graph
    LRkm1 = {}
    # loop over the sequences
    for kmer in sequences:
        # get the left and right k-1 mers
        left = kmer[:-1]
        right = kmer[1:]
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