# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 12:06:45 2015

@author: Richard
"""

import itertools


def organizing_strings(input_file):
    '''
    (file) -> list
    Grab the the collection C and n from the imput file
    C is a collection of at most 10 symbols defining an ordered alphabet,
    and a positive integer n  (n≤10).
    Return all strings of length n  that can be formed from the alphabet,
    ordered lexicographically
    '''
    # open file for reading
    infile = open(input_file,'r')
    C = infile.readline().rstrip().split()
    n = int(infile.readline().rstrip())
    
    # make a sequence of symbols
    seq = ''.join(C)
    
    # find all the possible substrings of length n
    subseqs = list(itertools.product(seq, repeat = n))
    
    
    # the itertools.product conviently order the strings based
    # on the order the input sequence 
    for item in subseqs:
        # convert the tuples into strings
        subseq = ''.join(item)
        print(subseq)
    
        
    
    
    
    
    
    
    
