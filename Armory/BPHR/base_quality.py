# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:10:24 2015

@author: Richard
"""
# import SeqIO to create a SeqRecord
from Bio import SeqIO
# import numpy
import numpy as np

# use biopython to analyze quality scores
def base_quality(input_file):
    '''
    (file) -> int
    Given a quality threshold q and a collection of reads in FASTQ format,
    return the number of positions where mean base quality falls below
    given threshold q
    Precondition: all reads have the same length
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get the quality score
    q = int(infile.readline().rstrip())
    # keep the file open as a handler for SeqIO
    reads = list(SeqIO.parse(infile, 'fastq'))
    
    # close file after parsing
    infile.close()
    
    # create a counter variable
    low_quality = 0    
    
    # loop over the base positions in all reads,
    # compute the mean quality score and update the counter
    for i in range(len(reads[0].letter_annotations['phred_quality'])):
        # create list to store the phred score at each position in all reads
        scores = []
        for j in range(len(reads)):
            scores.append(reads[j].letter_annotations['phred_quality'][i])
        # compute mean score for position i
        mean_score = np.mean(scores)
        # if mean_score < threshold: add position to counter
        if mean_score < q:
            low_quality += 1
            
    return low_quality
            
        
    
    
    
    
    
    
    