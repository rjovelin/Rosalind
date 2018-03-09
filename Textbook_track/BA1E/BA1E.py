# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 23:03:13 2015

@author: Richard
"""


def is_clumped(hits, L, t):
    '''
    (list, int, int) -> bool
    Take a list of indices of matching patterns, the length L of an interval,
    and return True if at least t indices are found within that interval
    and False otherwise
    '''

    # sort the list of indices of pattern matching
    hits.sort()
    
    # no clump if number of matches < t
    if len(hits) < t:
        return False
    else:
        # check if the patterns form clumps
       # loop starting at the last match position
       # evaluate the distance between the last matcg and the first match
       # then move on to second to last position, until a clump is found or
       # until all possible intervals have been evaluated
       for i in range(len(hits)-1, -1, -1):
           for j in range(len(hits)):
               interval = hits[i] - hits[j]
               N = i - j + 1
               if N >= t and interval <= L:
                   return True
                   
    # return False if no clump are found
    return False


def find_clump(genome, k, L, t):
    '''
    (str, int, int, int) -> list
    Return a list of kmers forming a clump of t occurences within an interval L
    in genome. 
    
    >>> find_clump('CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC', 5, 75, 4)
    [CGACA, GAAGA, AATGT]    
    '''
    
    # make a list of kmers forming clump
    clump = []
    
    # make a dicts of kmers of length k in genome
    kmers = {}
    
    # loop over genome and record kmers
    for i in range(0, len(genome) - k + 1):
        # extract kmer sequence from genome
        seq = genome[i:i+k]
        # check if kmer in dict
        if seq in kmers:
            # record match
            kmers[seq].append(i)
        else:
            kmers[seq] = [i]
            
    # loop over kmers
    for seq in kmers:
        # ask if kmer form clump
        if is_clumped(kmers[seq], L, t):
            # add kmer to list
            clump.append(seq)
    
    return clump
    
def ba1e(filename):
    infile = open(filename, 'r')
    # get genome
    genome = infile.readline().rstrip()
    nums = infile.readline().rstrip().split()
    k = int(nums[0])
    L = int(nums[1])
    t = int(nums[2])
    infile.close()

    # get the list of kmers forming clumps
    clump = find_clump(genome, k, L, t)
    
    for i in clump:
        print(i, end = ' ')
        
    