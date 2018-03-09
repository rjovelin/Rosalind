# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 20:47:15 2015

@author: Richard
"""

def dna_rev_compl(dna):
    '''
    (str) -> str
    Return the reverse complement of a DNA sequence
    Precondition: dna has only valid nucleotides
    
    >>> dna_rev_compl('AAAACCCGGT')
    ACCGGGTTTT
    '''
    
    # create a dictionary of complementary bases
    complement_nt = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # create a string collector
    revcomp = ''
    # loop over dna:
    for i in range(len(dna)):
        # gett he complement of nucleotide in dna and add to the begining of the 
        # new string to get the reverse sequence
        revcomp = complement_nt[dna[i].upper()] + revcomp
    
    return revcomp
    
    