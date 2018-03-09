# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 08:08:15 2015

@author: Richard
"""

def kmer_composition(input_file):
    '''
    Given a DNA string s in FASTA format in input file,
    Compute the 4-mer composition of s and return an array with the 
    counts of each 4-mer in lexicographic order in s    
    '''
        
    # open file for reading
    infile = open(input_file, 'r')
    dna = ''
    for line in infile:
        line = line.rstrip()
        if line != '' and not line.startswith('>'):
            dna += line
    infile.close()
    
    # make a list to store the 4-mers
    kmers = []
    for i in range(0, len(dna) -4+1):
        kmer = dna[i:i+4]
        kmers.append(kmer)
    # sort the kmers lexicographically (alphabet = 'ATCG')
    kmers.sort()
    
    # count each kmer in dna and add counts to list
    # cannot use the count builtin because it won't count overlapping kmers
    counts = []
    already_added = []
    
    for kmer in kmers:
        total = 0
        if kmer not in already_added:
            if dna.count(kmer) !=0:
                for i in range(len(dna) -4+1):
                    motif = dna[i:i+4]
                    if kmer == motif:
                        total += 1
            
            already_added.append(kmer)
            counts.append(total)
          
    for i in counts:
        print(i, end = ' ')
    

import itertools


def kmer_composition_2(input_file):
    '''
    Given a DNA string s in FASTA format in input file,
    Compute the 4-mer composition of s and return an array with the 
    counts of each 4-mer in lexicographic order in s    
    '''
        
    # open file for reading
    infile = open(input_file, 'r')
    dna = ''
    for line in infile:
        line = line.rstrip()
        if line != '' and not line.startswith('>'):
            dna += line
    infile.close()
    
    
    # make all possible 4-mers from alphabet {'A', 'T', 'C', 'G'}
    subseq = itertools.product('ATGC', repeat = 4)
    
    # make a list to store the 4-mers
    kmers = []
    
    # add all 4-mers to the list of kmers
    for item in subseq:
        kmers.append(''.join(item))
    
    # sort the kmers lexicographically (alphabet = 'ATCG')
    kmers.sort()
    
    # count each kmer in dna and add counts to list
    # cannot use the count builtin because it won't count overlapping kmers
    counts = []
        
    for kmer in kmers:
        total = 0
        if dna.count(kmer) !=0:
            for i in range(len(dna) -4+1):
                motif = dna[i:i+4]
                if kmer == motif:
                    total += 1
            counts.append(total)
        else:
            counts.append(total)
          
    for i in counts:
        print(i, end = ' ')     