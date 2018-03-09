# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 23:00:41 2015

@author: Richard
"""


# import SeqIO to create a SeqRecord
from Bio import SeqIO

def quality_scores(input_file):
    '''
    (file) -> int
    Given a quality threshold, along with FASTQ entries for multiple reads
    in input file, return the number of reads whose average quality is below
    the threshold.
    '''
    
    # open file for reading, create handle
    infile = open(input_file, 'r')
    # get the threshold score
    threshold = int(infile.readline().rstrip())
    # keep file open as handle for parsing
    
    # parse the reads into a Seqrecords into a list of SeqRecords
    reads = list(SeqIO.parse(infile, 'fastq'))
    
    # close file
    infile.close()
    
    # create variable counter for the low quality reads
    low_quality = 0
    
    # get the phred scores for each nucleotide, compute the mean and
    # compare to threshold
    for read in reads:
        # phred scores for each nucleotide can be obtained with
        # record.letter_annotations: dict with phred_quality: list of phred
        phred = read.letter_annotations['phred_quality']
        # compute mean by passing the list to numpy.mean
        mean = numpy.mean(phred)
        # compare to threshold, update counter
        if mean < threshold:
            low_quality += 1
    
    return low_quality
        
    
    
    
    
    
    
    
    