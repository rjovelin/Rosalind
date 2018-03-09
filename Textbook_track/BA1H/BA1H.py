# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 18:52:44 2015

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
    
    
    
def approximate_matching(pattern, genome, d):
    '''
    (str, str, int) -> list
    Return the indices where pattern matches genome with at most d differences
    
    >>>approximate_matching('ATTCTGGA', 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC', 3)
    [6, 7, 26, 27, 78]
    '''
    
    # set counter
    hits = []
        
    # loop through the text
    for i in range(0, len(genome) - len(pattern) + 1):
        # grab the sequence of length pattern in text
        seq = genome[i: i + len(pattern)].upper()
        # compute hamming distance between pattern and seq
        if hamming_distance(pattern.upper(), seq.upper()) <= d:
            # add a hit
            hits.append(i)
    return hits
    
    
def ba1h(filename):
    infile = open(filename, 'r')
    pattern = infile.readline().rstrip()
    genome = infile.readline().rstrip()
    d = int(infile.readline().rstrip())
    hits = approximate_matching(pattern, genome, d)
    for i in hits:
        print(i, end = ' ')
