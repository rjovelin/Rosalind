# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 20:56:53 2015

@author: Richard
"""

def find_pattern(pattern, genome):
    '''
    (str, str) -> list
    Return the indices at which pattern is found in genome, including indices of
    overlapping position
    
    >>> find_pattern('ATAT', 'GATATATGCATATACTT')
    1 3 9    
    '''
    # set counter
    hits = []
        
    # loop through the text
    for i in range(0, len(genome) - len(pattern) + 1):
        # grab the sequence of length pattern in text
        seq = genome[i: i + len(pattern)].upper()
        # check if seq and pattern are the same
        if pattern.upper() == seq:
            # add a hit
            hits.append(i)
    return hits