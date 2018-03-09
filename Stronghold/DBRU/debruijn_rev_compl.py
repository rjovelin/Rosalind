# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 15:41:19 2015

@author: Richard
"""

def reverse_complement(DNA):
    '''
    (str) -> str
    Return the reverse complement of string sequence DNA
    Precondition: DNA only include valid nucleotides 'A', 'T', 'C', 'G'
    '''
    
    seq = DNA.upper()
    
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rev_comp = ''
    for i in reversed(seq):
        rev_comp += complement[i]
    return rev_comp


def debruijn_graph(sequences):
    '''
    (list) -> dict
    Take a list of kmers and return the debruijn graph of the kmers as a dict
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
    

def debruijn_kmer_rev_compl(input_file, outputfile):
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
    kmers = []
    for line in infile:
        line = line.rstrip()
        if line != '':
            kmers.append(line)
    # close file
    infile.close()
    
    
    # get the reverse complement of every kmers
    kmers_rc = []
    for seq in kmers:
        kmers_rc.append(reverse_complement(seq))
    
    # take the union of the sets of kmers and rev compl
    sequences = set(kmers).union(set(kmers_rc))
    
    # build a debruijn graph with the set of kmers
    db_graph = debruijn_graph(sequences)
    
    # open file for writing
    newfile = open(outputfile, 'w')
       
    # write graph to outputfile
    for node in db_graph:
        for item in db_graph[node]:
            newfile.write('(' + node + ', ' + item + ')' + '\n')
    # close file after writing
    newfile.close()
    