# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 18:57:09 2015

@author: Richard
"""

# import SeqIo to create SeqRecords
from Bio import SeqIO

def convert_fastq(input_file):
    '''
    (file) -> dict
    Convert a file with reads in fastq format in a dictionnary
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # skip threshold
    infile.readline()
    # create a dictionnary to store reads ID: list of features
    # fastq[read_ID] = [sequence, '+', score]
    fastq = {}
    # read file, add i to ID, use ID as key, other line in list
    for line in infile:
        line = line.rstrip()
        if line != '':
            if line.startswith('@Rosalind'):
                name = line
                fastq[name] = []
            else:
                fastq[name].append(line)
    # close file
    infile.close()
    
    # create a lit of tuple (read, features) sorted by the read index number
    reads = [(read, features) for read, features in fastq.items()]
    reads.sort()
    
    return reads


def trimming_by_quality(input_file, outputfile):
    '''
    (file) -> file
    Given a quality cut-off value q, and a collections of reads in
    fastq with Phred33 quality scores, the fastq file trimmed from
    the both ends when base quality scores are lower than q
    '''
    
    # convert fastq into a list sorted by read index
    fastq = convert_fastq(input_file)
      
    # open file for reading
    infile = open(input_file, 'r')
    # get quality threshold q
    q = int(infile.readline().rstrip())
    # keep file open as handle
    # create a SeqRecords
    sequences = list(SeqIO.parse(infile, 'fastq'))
    
    # close file
    infile.close()
    
    # add a number and @ to sequence id so that each read can be identified
#    i = 1
    for read in sequences:
        read.id = '@' + read.id
              
    # create a dict to store read_ID: [left, right trimming]
    trimming = {}
    
    # check base quality scores
    for read in sequences:
        # create variable for the phred scores
        phred = read.letter_annotations['phred_quality']
        # create variables to store number of bases to trim
        left_trim = 0
        right_trim = 0
        # count the number of bases to left trim and right trim
        j = 0
        while phred[j] < q:
            left_trim += 1
            j += 1
        for score in reversed(phred):
            if score < q:
                right_trim += 1
            else:
                break
        # store left_trim and right_trim in dict
        trimming[read.id] = [left_trim, right_trim]
    
    # trim sequences in the fastq list
    for i in range(len(fastq)):
        read = fastq[i][0]
        left_trim, right_trim = trimming[read]
        # keep original seq if left and right = 0 
        # mote: trimming with left and right = 0 is [0:0] = []
        if not (left_trim == 0 and right_trim == 0):
            # trim sequence
            fastq[i][1][0] = fastq[i][1][0][left_trim: -right_trim]
            # trim scoring format
            fastq[i][1][-1] = fastq[i][1][-1][left_trim: -right_trim]
    
    
    # open file for writing
    newfile = open(outputfile, 'w')
    for pair in fastq:
        newfile.write(pair[0] + '\n')
        for item in pair[1]:
            newfile.write(item + '\n')
    # close file
    newfile.close()
    

    
        
    
    
        
        
    
    
    
    
    
    
    