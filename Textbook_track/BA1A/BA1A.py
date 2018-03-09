# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 20:15:21 2015

@author: Richard
"""



def count_pattern(text, pattern):
    '''
    (str, str) -> int
    Count the number of times pattern occurs in text with exact matches, 
    count also overlapping positions
    
    >>> count_pattern('GCGCG', 'GCG')
    2    
    '''
    # set counter
    hits = 0
        
    # loop through the text
    for i in range(0, len(text) - len(pattern) + 1):
        # grab the sequence of length pattern in text
        seq = text[i: i + len(pattern)].upper()
        # check if seq and pattern are the same
        if pattern.upper() == seq:
            # add a hit
            hits += 1
    return hits
    
    
