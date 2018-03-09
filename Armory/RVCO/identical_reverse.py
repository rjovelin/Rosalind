# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 18:43:05 2015

@author: Richard
"""




# import SeqIO from biopython to create a SeqRecord object
from Bio import SeqIO

# use Biopython to compare DA sequences with their reverse complement
def identical_reverse(input_file):
    '''
    (file) -> int
    Given a collection of n DNA strings from input file in fasta format,
    return the number of given strings that match their reverse complements
    '''
    
    # parse the fasta file into a list of SeqRecords
    sequences = list(SeqIO.parse(input_file, 'fasta'))
    
    # set up counter
    total = 0
    
    # loop over SeqRecords, compare sequence and reverse complement
    for DNA in sequences:
        if DNA.seq == DNA.seq.reverse_complement():
            total += 1
    
    return total
    
    
    
    