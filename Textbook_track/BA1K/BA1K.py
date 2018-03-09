# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 01:19:14 2015

@author: Richard
"""

import itertools

def frequency_array(T, k):
    '''
    (T, k) -> list
    Given an integer k, and a string T, return a list in which the i-th element
    of the list holds the number of times that the i-th k-mer (in lexicographic order)
    appears in T.
    
    >>> frequency_array('ACGCGGCTCTGAAA', 2)
    2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0
    '''
    
    alphabet = 'ACTG'
    # create all possible kmers from alphabet
    seqs = set()
    for i in itertools.combinations_with_replacement(alphabet, k):
        for j in (list(itertools.permutations(i))):
            seqs.add(''.join(list(j)))
    # make a list of kmers
    seqs = list(seqs)
    # sort list in lexicographical order
    seqs.sort()
    
    # create a dict with kmer : counts from T
    kmers = {}
    for i in range(len(T) - k + 1):
        # get kmer seq
        seq = T[i:i+k]
        # populate dict
        if seq not in kmers:
            kmers[seq] = 1
        else:
            kmers[seq] += 1
    
    for j in kmers:
        if j not in seqs:
            print(j)
    
    # creat a dict of seq: counts
    counts = {}
    for i in seqs:
        if i in kmers:
            counts[i] = kmers[i]
        else:
            counts[i] = 0
    
           
    # create a list of (kmer, count)
    freq = [(key, val) for key, val in counts.items()]
    # sort list lexicographically
    freq.sort()
    # create a list of freq
    array = [i[1] for i in freq] 
    return array
    
    
def ba1k(filename):
    infile = open(filename, 'r')
    T = infile.readline().rstrip()
    k = int(infile.readline().rstrip())
    infile.close()
        
    array = frequency_array(T, k)
    for i in array:
        print(i, end = ' ')
        
        
 
        