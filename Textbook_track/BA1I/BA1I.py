# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 19:57:03 2015

@author: Richard
"""

import itertools

def hamming_distance(seq1, seq2):
    '''
    (str, str) -> int
    Return the hamming distance (number of differences) between seq1 and seq2
    Precondition: seq1 and seq2 have same length
    
    >>> hamming_distance('GGGCCGTTGGT', 'GGACCGTTGAC')
    3
    '''
    D = 0
    
    # loop over seq1
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            D += 1
    return D


def most_frequent_mismatched_kmer(genome, k, d):
    '''
    (str, int, int) -> list
    Return a list with the most frequent kmers in genome with at most d mismatches
    
    >>> most_frequent_mismatched_kmer('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1)    
    ['GATG', 'ATGC', 'ATGT']
    '''
    
    # create to store kmers : count
    kmers = set()
    alphabet = 'ACTG'
    for i in itertools.combinations_with_replacement(alphabet, k):
        for j in itertools.permutations(''.join(list(i))):
            kmers.add(''.join(list(j)))
    
    # create dict to store kmer : counts
    counts = {}
    
    # loop over kmers
    for seq in kmers:
        # loop over genome
        for i in range(len(genome) - k + 1):
            # grab kmer sequence
            subseq = genome[i:i+k].upper()
            # compute distance:
            if hamming_distance(seq, subseq) <= d:
                # check if kmer in dict
                if seq in counts:
                    counts[seq] += 1
                else:
                    counts[seq] = 1
                    
    # make a reverse dict : count : list of kmers
    counts_to_kmers = {}
    for seq in counts:
        if counts[seq] in counts_to_kmers:
            counts_to_kmers[counts[seq]].append(seq)
        else:
            counts_to_kmers[counts[seq]] = [seq]
            
    # make a list of counts
    L = [i for i in counts_to_kmers]
    highest = max(L)
    return counts_to_kmers[highest]
    
    
def ba1i(filename):
    infile = open(filename, 'r')
    genome = infile.readline().rstrip()
    k , d = infile.readline().rstrip().split()
    k = int(k)
    d = int(d)
    infile.close()
    most_freq = most_frequent_mismatched_kmer(genome, k, d)
    for i in most_freq:
        print(i, end = ' ')
    