# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 18:45:48 2015

@author: Richard
"""

def longest_repeat(input_file):
    '''
    (file) -> str
    Given a string Text, return the longest substring of Text appearing more than once
    '''
    # open file for reading
    infile = open(input_file, 'r')
    # create a variable to hold the Text
    dna = ''
    for line in infile:
        line = line.rstrip()
        if line != '':
            dna += line
    # close file
    infile.close()
    
    # set up variable   
    longest = ''
    max_length = 0
    
    for i in range(len(dna)):
        for k in range(len(dna), i, -1):
            motif = dna[i:i+k]
            if dna.count(motif) > 1and len(motif) > max_length:
                longest = motif
                max_length = len(motif)
    
    print(longest)
            
        
    
    