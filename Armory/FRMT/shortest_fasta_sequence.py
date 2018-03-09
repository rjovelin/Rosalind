# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 15:31:02 2015

@author: Richard
"""


# import module Entrez from Biopython to access genbank
from Bio import Entrez
# import module SeqIO to create a seq object from a sequence format
from Bio import SeqIO

def shortest_fasta_sequence(input_file):
    '''
    (file) -> 
    Given a collection of n GenBank entry IDs in input file, return the
    shortest of the strings associated with the IDs in FASTA format.
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get the list of GenBank IDs
    IDs = infile.readline().rstrip().split()
    # close file
    infile.close()    
    
    # access genbank IDs using efetch
    # specificy the database to search with db
    # specify the sequence fasta to return with rettype
    # use id to pass the ID of the sequence (or a list of IDs)    
    
    # create a handle
    handle = Entrez.efetch(db = 'nucleotide', id= IDs, rettype = 'fasta')
    
    # handle.read() returns a string in fasta format, not a SeqRecord
    # need to read the handle before parsing it
        
    # use SeqIO to create a SeqRecord, return a generator: use list to access
    # content
    sequences = list(SeqIO.parse(handle, 'fasta'))
    
    
    # find the sequence with the shortest sequence
    # assign shortest seq to first sequence
    shortest = len(sequences[0].seq)
    shortest_seq = sequences[0]
    # compare the other sequences to shortest, updating shortest sequence
    for i in range(1, len(sequences)):
        if len(sequences[i].seq) < shortest:
            shortest_seq = sequences[i]
            shortest = len(sequences[i].seq)
    
    # get the string in fasta format of the shortest sequence
    for name in IDs:
        if name in shortest_seq.id:
            prot_ID = name
    new_handle = Entrez.efetch(db = 'nucleotide', id= prot_ID, rettype = 'fasta')
    # use read method to get the string in fasta
    sequence = new_handle.read()
    
    return sequence
    
    
    
    
    
