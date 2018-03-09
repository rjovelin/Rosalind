# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 12:20:07 2015

@author: Richard
"""

# use Biopython to count nucleotides
from Bio.Seq import Seq

def nucleotide_counts(input_file):
    '''
    (file) -> list
    Use Biopython to count the number of nucleotides 'A', 'C', 'G' and 'T' 
    in string S in input file
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    S = infile.readline().rstrip()
    
    # create alphabet to loop over
    bases = 'ACGT'
    
    # create list to store the counts
    counts = []
    
    # use Biopython Seq object to count nucleotides
    # create a sequence object
    DNA = Seq(S)
    
    for nucleotide in bases:
        counts.append(DNA.count(nucleotide))
    
    # print counts to screen
    for i in counts:
        print(i, end = ' ')
    


    
    