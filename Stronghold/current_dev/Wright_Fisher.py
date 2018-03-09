# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 23:38:06 2015

@author: Richard
"""

from scipy.misc import comb


def binomial(N, k, p):
    '''
    (int, int, float) -> float
    Return the Binomial probability of observing k successes in N trials
    given that the probability of a 1 sucess is p
    '''
    # the probability of k successes in n trials is given by the binomial
    # distribution P(k) = C(n,k) * p^k * q^(n-k)
    
    P = comb(N, k) * p**k * (1-p)**(N-k)
    
    return P


def genetic_drift(input_file):
    '''
    (file) -> float
    Given 4 positive integers N, m, g and k, return the probability that
    in a population of N diploid individuals initially possessing m
    copies of a dominant allele, we will observe after g generations at
    least k copies of a recessive allele. Assume the Wright-Fisher model.
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get N, m, g and k
    nums = infile.readline().rstrip().split()
    N = int(nums[0])
    m = int(nums[1])
    g = int(nums[2])
    k = int(nums[3])
    
    # compute the change in frequencies over g generations
    
#    for i in range(g):
#        P = binomial(2*N, )
    
    
    pass
    
    