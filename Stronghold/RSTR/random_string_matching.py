# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 18:09:17 2015

@author: Richard
"""

from scipy import misc

def proba_random_string(s, GC):
    '''
    (str, float) -> float
    Given a string s and a GC content GC, return the probability
    that a random string with GC content equal to GC is identical to s
    '''
    
    # compute the probabilities of each nucleotide given GC content
    pCG = GC / 2
    pAT = (1- GC) / 2
    P = 0

    # compute the probability P that the random string matches sequence s
    if s[0] == 'A' or s[0] == 'T':
        P += pAT
    else:
        P += pCG

    for i in range(1, len(s)):
        if s[i] == 'A' or s[i] == 'T':
            P *= pAT
        elif s[i] == 'C' or s[i] == 'G':
            P *= pCG

    return P


def random_string_matching(input_file):
    '''
    (file) -> float
    Given a positive integer N ≤ 100000, a number x between 0 and 1,
    and a DNA string s of length at most 10 bp. Return the probability that
    if N random DNA strings having the same length as s are constructed
    with GC-content x then at least one of the strings equals s.
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get N, s, GC
    nums = infile.readline().rstrip().split()
    N = int(nums[0])
    GC = float(nums[1])
    s = infile.readline().rstrip()
    # close file
    infile.close()
    
    # probability that a random string with GC content GC equals s
    p = proba_random_string(s, GC)
    
    # probability that no random string equals given by binomial distribution
    P0 = misc.comb(N, 0, exact = True) * p**0 * (1-p)**(N-0)
    
    # probability that at least 1 string equals s
    Q = 1 - P0
    
    return round(Q, 3)    
    
    
    
    
    
    
    
    
    