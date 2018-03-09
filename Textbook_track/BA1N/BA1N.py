# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 14:04:04 2015

@author: RJovelin
"""

import itertools


def hamming_distance(seq1, seq2):
    '''
    (str, str) -> int
    Return the number of differences between seq1 and seq2
    Precondition: seq1 and seq2 have same length
    '''
    
    # set up counter
    D = 0
    # loop over seq1
    for i in range(len(seq1)):
        if seq1[i]  != seq2[i]:
            D += 1
    return D
    


def d_neighborhood(kmer, d):
    '''
    (str, int) -> list
    Return a list of strings for which the hamming distance with kmer does
    not exceed d
    
    >>> d_neighborhood('ACG', 1)
    ['CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG', 'ACA', 'ACC', 'ACT', 'ACG']
    '''
    
    alphabet = 'ACGT'    
       
    # create a set of sequence of kmer's length
    seqs = set()
    for i in itertools.product(alphabet, repeat = len(kmer)):
        seqs.add(''.join(list(i)))
    
    # create a list of neighbors
    neighbors = [i for i in seqs if hamming_distance(i, kmer) <= d]
    return neighbors
    
    
def ba1n(filename):
    infile = open(filename, 'r')
    newfile = open('ba1n_out.txt', 'w')
    kmer = infile.readline().rstrip()
    d = int(infile.readline().rstrip())
    infile.close()
    seqs = d_neighborhood(kmer, d)
    for i in seqs:
        newfile.write(i + '\n')
    newfile.close()
        
    
        
    
    
    
    