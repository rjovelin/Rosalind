# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:06:40 2015

@author: Richard
"""

# import module Entrez to access NCBI
from Bio import Entrez
# import SeqIO to create a SeqRecord
from Bio import SeqIO


def fetch_NCBI_sequences(input_file):
    '''
    (file) -> list
    Given 2 GenBank IDs in input file, return a list of Biopython
    SeqRecord objects corresponding to the GenBank IDs
    Print the IDs and their sequence
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get the GenBank IDs 
    IDs = infile.readline().rstrip().split()
    
    # create a handle to get information on thte GenBank sequences
    handle = Entrez.efetch(db = 'nucleotide', id= IDs, rettype = 'fasta')
    # parse handle to get SeqRecords
    sequences = list(SeqIO.parse(handle, 'fasta'))
    
    for seq in sequences:
        print(seq.id, seq.seq, end = '\n')
   
# Note: for solving Rosalind problem need
# Go to Needle online alignment tool
# http://www.ebi.ac.uk/Tools/psa/emboss_needle/nucleotide.html
# set options in needle with:
# gap opening penalty = 10
# gap extension penalty  =1
# output format = pair
# matrix = DNAfull
# end gap penalty = True
# end gap open = 10
# end gap extend = 1

# return the alignment score    
    
