# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 01:55:35 2015

@author: Richard
"""

def minimum_skew(input_file):
    '''
    (file) -> None
    Take the sequence S in the inout file, compute the skew (difference
    between G and C) for each prefix i for S of length i and return
    all integer(s) i minimizing Skew(Prefixi(S)) over all values
    of i (from 0 to |Genome|)
    '''
    
    # open file
    infile = open(input_file, 'r')
    S = infile.readline().rstrip()
    # close file
    infile.close()
    
    # make a list to store the skew
    skew = []
    # count the sequence from 1 to end (0 is prefix of the 1st nucleotide)
    skew.append(0)
    
    # make a sequence that will accumulate each nucleotide
    genome = ''
    for nucleotide in S:
        genome += nucleotide
        C = genome.count('C')
        G = genome.count('G')
        skew.append(G - C)
    
    # find the position i minimizing skew
    # get the lower value of skew
    minimum = min(skew)
    # count how many times the lower value is reached
    minimum_count = skew.count(minimum)
    # make a list to store the positions i of the minimum
    minimum_i = []
    # find first position of minimum
    position_i = skew.index(minimum)
    minimum_i.append(position_i)
    # substract 1 to the counter
    minimum_count -= 1
    
    # find all occurence of minimum, starting after the first position
    while minimum_count != 0:
        position_i = skew.index(minimum, position_i + 1)
        minimum_count -= 1
        minimum_i.append(position_i)
    
    if len(minimum_i) !=0:
        for item in minimum_i:
            print(item, end = ' ')
    
        
        
        
    
    
        
    
    
       
    
    
    
    