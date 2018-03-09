# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 19:00:35 2015

@author: Richard
"""

from Bio import SeqIO


# use biopython to convert fastq file into a fasta file
def fastq_to_fasta(input_file, outputfile):
    '''
    (file) -> file
    Convert a FASTQ file into a FASTA file
    '''
    
    # convert a fastq file into a fasta file with SeqIO
    SeqIO.convert(input_file, 'fastq', outputfile, 'fasta')
    

    
    
    
    