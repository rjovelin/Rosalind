# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 16:27:54 2015
@author: Richard
"""

import itertools

def find_median_string(input_file):
    '''
    (file) -> str
    Given an integer k and a collection of strings Dna in input file,
    return a k-mer Pattern that minimizes the total distance d(Pattern, Dna)
    over all k-mers Pattern. Total distance noted d(pattern, DNA) is the
    sum of the minimum Hamming distances for each DNA sequence i
    in the collection. Return only 1 k-mer if multiple k-mer exists    
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get k
    k = int(infile.readline().rstrip())
    # make a list of DNA sequences
    sequences = []
    for line in infile:
        line = line.rstrip()
        if line != '':
            sequences.append(line)
    # close file
    infile.close()
    
    # create a dictionnary to hold the pattern: total_distance pairs
    pattern_distance = {}
      
    # go through each sequence
    for i in range(len(sequences)):
        # loop over the sequence and grab a motif of length k
        for j in range(len(sequences[i])):
            patterns = itertools.permutations(sequences[i][j:j+k], k)
            # loop over the collection of pattern
            for seq in patterns:
                # join the nucleotides in seq as a string sequence
                seq = ''.join(seq)
                # initialize the dict with list that will hold all the lowest d
                # for each seq, then sum all the d
                if seq not in pattern_distance:
                    # initialize total distance D
                    pattern_distance[seq] = 0
                    # scan each seq to find its minimum hamming distance
                    for m in range(len(sequences)):
                        # create a list to hold distances, initialize for each sequence
                        distances = []
                        # go through the sequence and grab motif of length k
                        for r in range(len(sequences[m]) -k+1):
                            motif = sequences[m][r:r+k]
                            # set up distance
                            diff = 0
                            # compute d for each motif
                            for q in range(len(seq)):
                                if seq[q] != motif[q]:
                                    diff +=1
                            # store distance in list
                            distances.append(diff)
                        # find the lowest hamming distance and add to the pattern D
                        hamming_d = min(distances)
                        pattern_distance[seq] += hamming_d
                    
    # make a list of (D, pattern) sublists
    D_pattern = [[d, subseq] for subseq, d in pattern_distance.items()]
    # sort list
    D_pattern.sort()
    
    print(D_pattern[0][1])