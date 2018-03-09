# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

def mRNA_from_protein(input_file):
    '''
    (file) -> int
    Returns the number of possible RNA sequences modulo 1000000
    corresponding to the protein sequences found in the input_file
    '''
    # open file for reading    
    infile = open(input_file, 'r')
    
    # grab the protein sequence
    protein = infile.readline().rstrip()
    
    # close file after reading
    infile.close()
    
    # define the genetic code
    reverse_genetic_code = {'A': ['GCU', 'GCC', 'GCA', 'GCG'],
                            'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
                            'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                            'K': ['AAA', 'AAG'],
                            'N': ['AAU', 'AAC'],
                            'M': ['AUG'],
                            'D': ['GAU', 'GAC'],
                            'F': ['UUU', 'UUC'],
                            'C' : ['UGU', 'UGC'],
                            'P': ['CCU', 'CCC', 'CCA', 'CCG'],
                            'Q': ['CAA', 'CAG'],
                            'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
                            'E': ['GAA', 'GAG'],
                            'T': ['ACU', 'ACC', 'ACA', 'ACG'],
                            'G': ['GGU', 'GGC', 'GGA', 'GGG'],
                            'W' : ['UGG'],
                            'H': ['CAU', 'CAC'], 
                            'Y': ['UAU', 'UAC'],
                            'I': ['AUU', 'AUC', 'AUA'], 
                            'V': ['GUU', 'GUC', 'GUA', 'GUG'],
                            '*': ['UAA', 'UGA', 'UAG']} 
    
    # count the number of possible RNA sequences
    RNAs = 0
    
    # count the number of possible RNAs corresponding to the first AA
    RNAs = len(reverse_genetic_code[protein[0]])
    
    # loop over the protein sequence, omitting the first AA     
    for i in range(1, len(protein)):
        RNAs *= len(reverse_genetic_code[protein[i]])
    # add the number of possible sequences due to the STOP codon
    RNAs *= len(reverse_genetic_code['*'])
    
    print(RNAs % 1000000)