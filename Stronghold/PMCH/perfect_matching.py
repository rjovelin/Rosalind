# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 16:14:05 2015

@author: Richard
"""

import scipy.misc as msc


def perfect_matching(input_file):
    '''
    (file) -> int
    Grab the RNA sequence S in the input file and return the total possible
    number of perfect matching of basepair edges in the bonding graph of S
    Precondition, S has the same number of 'A' and 'U', and the same number
    of 'C' and 'G'
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # grab RNA seq
    S = ''
    for line in infile:
        if line != '\n' and not line.startswith('>'):
            line = line.rstrip()
            S += line
        
    S = S.upper()
    
    # count the number of A (nA = nU), and the number of C (nC = nG)
    A = S.count('A')
    C = S.count('C')
    
    # don't use approximation when computing factorial
    print(int(msc.factorial(A, True) * msc.factorial(C, True)))
    print(S)
    
    
    
