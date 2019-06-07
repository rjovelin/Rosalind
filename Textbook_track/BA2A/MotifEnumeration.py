# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 01:53:00 2015
@author: Richard
"""

import itertools

def kd_motif(input_file):
    '''
    (file) -> list
    Get integers k and d and a collection of DNA sequences from the input file.
    Return all (k, d)-motifs in Dna.
    (k, d)-motifs are motifs x such that lenght(x) = k and every input DNA seq
    has at least 1 variant of x at a Hamming distance of at most d
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # grab k and d
    nums = infile.readline().rstrip().split()
    k = int(nums[0])
    d = int(nums[1])
    # make a list with all sequences 
    sequences = []
    for line in infile:
        line = line.rstrip()
        if line != '':
            sequences.append(line)
    # close file
    infile.close()
    
    
    # create a set to hold the (k, d)-motifs
    kmer = set()
  
    # go through each seq
    for i in range(len(sequences)):
        print(i, len(sequences))
        # go through the sequences, grab a subseq of length k
        # and find all possible motifs corresponding to that subseq
        for r in range(len(sequences[i]) -k+1):
            subsequences = sequences[i][r:r+k]
            patterns  = itertools.permutations(subsequences, k)
            # ask is each subsequence has a variant in all sequences
            for subseq in patterns:
                # no need to check if subseq alread recorded
                if ''.join(subseq) not in kmer:
                    # set up boolean
                    found_in_all = True
                    # initialise the counter
                    j = 0
                    # keep searching for a variant
                    # stop when no variant is found or when all sequences are examined
                    while j != len(sequences) and found_in_all == True:
                        # set up boolean
                        found_it = False
                        # scan sequence to see if a variant exists
                        for m in range(len(sequences[j]) -k+1):
                            # grab sequence motif and initialise diff counter
                            motif = sequences[j][m:m+k]
                            diff = 0
                            # count number of diffs between subseq and motif
                            for n in range(len(subseq)):
                                if subseq[n] != motif[n]:
                                    diff +=1
                            # check if the motif is a variant
                            if diff <= d:
                                found_it = True # motif is a variant
                                break # stop scanning the sequence
                        # after checking sequence, exit while loop if no variant                
                        if found_it == False:
                            found_in_all = False # exit while loop and grab next subseq
                        j += 1 # accumulate counter
                    if found_in_all == True and j == len(sequences):
                        kmer.add(''.join(subseq))
                        
    if len(kmer) != 0:
        for subseq in kmer:
            print(subseq, end= ' ')
    else:
        print('there are no (k, d) motif')
        print(k, d)
        print(sequences) 