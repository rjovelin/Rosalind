# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 23:57:07 2015

@author: Richard
"""

# import SeqIO to create a SeqRecord
from Bio import SeqIO

def quality_filtering(input_file):
    '''
    (file) -> int
    Given a quality threshold value q, percentage of bases p, and a set
    of FASTQ entries, return the number of reads in filtered FASTQ entries
    Filter reads for which p% of bases have a phred score below threshold
    '''
    
    # open file for reading, create handle
    infile = open(input_file, 'r')
    # get the threshold and cutoff values
    q, p = infile.readline().rstrip().split()
    q = int(q)
    p = int(p)
        
    # keep file open as handle for parsing
    # parse the reads into a Seqrecords into a list of SeqRecords
    reads = list(SeqIO.parse(infile, 'fastq'))
    
    # close file
    infile.close()
    
    # create list of reads to filter out
    filtered = []
        
    # compute the F: % of bases with score < q for each read
    for read in reads:
        # phred scores for each nucleotide can be obtained with
        # record.letter_annotations: dict with phred_quality: list of phred
        phred = read.letter_annotations['phred_quality']
        
        # compute % of bases with phred > q
        highs = 0
        for score in phred:
            if score >= q:
                highs += 1
        
        # filter reads if fraction of bases with phred > q is < p
        if (highs/len(phred)) * 100 < p:
            filtered.append(read)
                
    # return the number of reads that pass filtering
    return len(reads) - len(filtered)
    
          
