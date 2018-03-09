# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 16:07:38 2015

@author: Richard
"""

def string_composition(input_file, outputfile):
    '''
    (file) -> file
    Get the integer k and the string s from the the input file
    Return the Composition-k-(s), where the k-mers are written in lexicographic order
    '''
    
    # open filr for reading
    infile = open(input_file, 'r')
    
    # grab k and S
    k = int(infile.readline().rstrip())
    S = infile.readline().rstrip()
    
    # close file after reading
    infile.close()
    
    # create list to store the k-mer substrings
    k_mers = []
    i = 0
    while len(S[i: i+k]) == k:
        substring = S[i: i+k]
        k_mers.append(substring)
        i += 1
    
    # sort the list of k-mers to return the substrings in the lexicographic order
    k_mers.sort()
    
    # open file for writing
    newfile = open(outputfile, 'w')
    for substring in k_mers:
        newfile.write(substring + '\n')
    # close file after writing
    newfile.close()
