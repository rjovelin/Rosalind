# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 22:02:52 2015

@author: Richard
"""

# use biopython to translate DNA sequences and manipulate sequences
from Bio.Seq import Seq


def find_longest_ORF(input_file):
    '''
    (file) -> str
    
    Given a DNA string S in input file, return the longest protein string
    that can be translated from an ORF of S.
    (return any protein if multiple proteins with longest length exist)
    Note: an ORF starts with the initiation codon, and ends either with 
    a stop codon or with the end of the DNA sequence    
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get DNA sequence and create a SeqObject
    DNA = Seq(infile.readline().rstrip())
    # close file
    infile.close()
    
    # create variable for protein length
    size = 0
    # check translations of all DNA sequences
    for i in range(len(DNA)):
        # stop translation at the first stop codon
        protein = DNA[i:].translate(to_stop = True)
        # keep if peptide length > size and peptide starts with M
        if protein.startswith('M') and len(protein) > size:
            longest = protein
            size = len(protein)
    
    # check ORF on the reverse strand
    DNArc = DNA.reverse_complement()
    for i in range(len(DNArc)):
        # stop translation at the first stop codon
        protein = DNArc[i:].translate(to_stop = True)
        # keep if peptide length > size and peptide starts with M
        if protein.startswith('M') and len(protein) > size:
            longest = protein
            size = len(protein)
    
    return longest
    
    