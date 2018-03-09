# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 21:54:53 2015

@author: Richard
"""

def failure_array(input_file):
    '''
    The failure array of s is an array P of length n for which P[k] is 
    the length of the longest substring s[j:k] that is equal to some
    prefix s[1:k−j+1], where j cannot equal 1 (otherwise, P[k] would always
    equal k). By convention, P[1]=0.

Given: A DNA string s (of length at most 100 kbp) in FASTA format.

Return: The failure array of s.
    
    
    '''
    
    
    
    # open file
    infile = open(input_file, 'r')
    dna = ''    
    for line in infile:
        if not line.startswith('>') and line != '':
            dna += line.rstrip()
    infile.close()
    
    # create array and intialize with 0 values
    P = [0]*len(dna) 
    # initialize k to 0
    k = 0 
    # loop over dna, starting at 2, j cannot be == 1
    for j in range(2, len(dna)): 
        
        while k > 0 and dna[k] != dna[j-1]: 
            k = P[k-1] 
        if dna[k] == dna[j-1]: 
           k += 1 
        P[j-1] = k 
    
    return ' '.join(map(str, P))
    
        