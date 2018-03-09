# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 19:31:03 2015

@author: Richard
"""

def multiple_pattern_matching(input_file):
    '''
    (file) -> list
    Given a string Text and a collection of strings Patterns from input file,
    return all indices in Text where a string Pattern appears    
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get the dna sequence
    dna = infile.readline().rstrip()
    # make a dict of patterns
    i = 0
    patterns = {}
    for line in infile:
        line = line.rstrip()
        if line != '':
            patterns[i] = line
            i += 1
    
    # close file after reading
    infile.close()
    
    # create list to store motif indices
    positions = []
    
    # loop over dna for each seq and report indices
    for seq in patterns:
        motif = patterns[seq]
        k = len(motif)
        for i in range(0, len(dna) -k+1):
            kmer = dna[i:i+k]
            if kmer == motif:
                positions.append(i)
    
    # sort list
    positions.sort()
    
    for i in positions:
        print(i, end = ' ')
            
        
    
    
    
    
    