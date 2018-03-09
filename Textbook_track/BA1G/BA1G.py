# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 01:08:08 2015

@author: Richard
"""

def hamming_distance(seq1, seq2):
    '''
    (str, str) -> int
    Return the hamming distance (number of differences) between seq1 and seq2
    Precondition: seq1 and seq2 have same length
    
    >>> hamming_distance('GGGCCGTTGGT', 'GGACCGTTGAC')
    3
    '''
    D = 0
    
    # loop over seq1
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            D += 1
    return D
    
def ba1g(filename):
    infile = open(filename, 'r')
    seq1 = infile.readline().rstrip()
    seq2 = infile.readline().rstrip()
    
    print(hamming_distance(seq1, seq2))