# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 19:25:20 2015

@author: Richard
"""

def minimum_GC_skew(genome):
    '''
    (str) -> list
    Return a list with the positions in genome minimizing the GC skew.
    The GC skew is the difference between the total number of occurrences of 'G'
    and 'C' for in Genome. 
    
    >>> minimum_GC_skew('CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG')
    [53, 97]
    '''
    
    # create lists
    minimum , skew, = [], []
        
    # loop over genome, grab prefix sequence, compute GC skew
    for i in range(len(genome)):
        # grab prefix
        prefix = genome[:i]
        # check if empty string
        if prefix == '':
            skew.append(0)
        else:
            skew.append(prefix.count('G') - prefix.count('C'))
    
    # find the minimum skew
    lowest = min(skew)
    # loop over indices of skew
    positions = [i for i in range(len(skew)) if skew[i] == lowest]
    
    return positions
    
    
def ba1f(filename):
    infile = open(filename, 'r')
    genome = infile.readline()
    infile.close()
    skew = minimum_GC_skew(genome)
    for i in skew:
        print(i, end = ' ')