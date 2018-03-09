# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import math


def random_string(input_file):
    '''
    (file) -> float
    Grab the sequence s and the array A of GC content from input_file
    Return an array B having the same length as A in which B[k]
    represents the common logarithm of the probability that a random string
    constructed with the GC-content found in A[k] will match s exactly    
    '''
    # open the file for reading
    infile = open(input_file, 'r')
    # grab the sequence
    s = infile.readline().rstrip().upper()
    # make a list with the GC content
    GC = [float(i) for i in infile.readline().split()]
    
    # make a list to store the log of P
    Prob = []
    
    # compute P for each GC content
    for GC_content in GC:
        # compute the probabilities of each nucleotide given GC content
        pCG = GC_content / 2
        pAT = (1- GC_content) / 2
        P = 0
        
        # compute the probability P that the random string matches sequence s
        if s[0] == 'A' or s[0] == 'T':
            P += pAT
        else:
            P += pCG
        
        for i in range(1, len(s)):
            if s[i] == 'A' or s[i] == 'T':
                P *= pAT
            elif s[i] == 'C' or s[i] == 'G':
                P *= pCG
        
        # populate list with P
        Prob.append(round(math.log10(P), 3))
        
    # print the elements of the list to screen        
    for item in Prob:
        print(item, end = '\t')
        
        
    
    