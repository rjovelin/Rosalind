# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:02:34 2015

@author: Richard
"""

# import ExPASy to access the Uniprot database
from Bio import ExPASy
# import Swissprot to create a Sequence object
from Bio import SwissProt

# use Biopython to get the protein sequences
def fetch_uniprot_sequences(input_file):
    '''
    (file) -> SwissProt.record
    Given 2 Uniprot IDs in input file, return a list of SwissProt.record objects
    Print each protein ID and their protein sequences
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get the protein IDs
    prot_IDs = infile.readline().rstrip().split()
    # close file
    infile.close()
    
    # create list to store the records
    prot_sequences = []
    
    # for each id in prot_ID: create  Swissprot.record object, append to a list
    for protein_name in prot_IDs:
        # create a handle to get protein information
        handle = ExPASy.get_sprot_raw(protein_name)
        # read handle with Swissprot to create a record
        protein = SwissProt.read(handle)
        # store record in list
        prot_sequences.append(protein)
        
    for protein in prot_sequences:
        print(protein.sequence)
    
# Note: for solving Rosalind problem need:
# Go to Water online alignment tool
# http://www.ebi.ac.uk/Tools/psa/emboss_water/
# set options in Water with:
# gap opening penalty = 10
# gap extension penalty = 1
# output format = pair
# matrix = BLOSUM62
# return the alignment score    