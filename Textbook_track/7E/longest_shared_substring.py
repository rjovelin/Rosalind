# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 19:01:19 2015

@author: Richard
"""

def longest_shared_substring(input_file):
    '''
    (file) -> str
    Given 2 strings Text1 and Text2 in input_file, return the longest
    substring that occurs in both Text1 and Text2.
    '''
    # open file for reading
    infile = open(input_file, 'r')
    # create a variable to hold the Text
    dna1 = infile.readline().rstrip()
    dna2 = infile.readline().rstrip()
    
    # close file
    infile.close()
    
    # set up variable   
    longest_shared = ''
    max_length = 0
    
    # go through dna1, identify the longest substring of dna1 in dna2
    for i in range(len(dna1)):
        # reduce dna1 background
        for k in range(len(dna1), 0, -1):
            motif = dna1[i:i+k]
            if motif in dna2 and len(motif) > max_length:
                longest_shared = motif
                max_length = len(motif)
    
    # go through dna2, identify the longest substring of dna2 in dna1
    for i in range(len(dna2)):
        # reduce dna1 background
        for k in range(len(dna2), 0, -1):
            motif = dna2[i:i+k]
            if motif in dna1 and len(motif) > max_length:
                longest_shared = motif
                max_length = len(motif)
    
    print(longest_shared)