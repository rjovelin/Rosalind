# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 20:29:35 2015

@author: Richard
"""


def get_seq_k(file):
    '''
    (file) -> str, int
    Extract a sequence and integet k from file
    '''
    # extract seq and k from file
    infile = open(file, 'r')
    seq = infile.readline()
    k = int(infile.readline())
    infile.close()
    return seq, k
    


   
    



def find_most_frequent_kmer(seq, k):
    '''
    (str, int) -> list
    Print the most frequent kmers of length k in sequence
    
    >>> find_most_frequent_kmer('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4)
    CATG GCAT
    '''
    
    # create a dict to store the k-mer : count pairs
    kmers = {}
    
    # loop through the text, record all kmer of length k and their counts
    for i in range(len(seq) - k + 1):
        # grab kmer sequence
        subseq = seq[i:i+k]
        # check if kmer is in key in dict
        if subseq not in kmers:
            # add key : 1
            kmers[subseq] = 1
        else:
            # add 1
            kmers[subseq] += 1
            
    # reverse dict
    counts = {}
    for i in kmers:
        if kmers[i] in counts:
            counts[kmers[i]].append(i)
        else:
            counts[kmers[i]] = [i]
            
    # make a list of counts
    freq = [i for i in counts]
    # get maximum count
    highest = max(freq)
    for subseq in counts[highest]:
        print(subseq, end = ' ')
    