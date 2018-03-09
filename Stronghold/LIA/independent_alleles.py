# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 17:46:10 2015

@author: Richard
"""
from math import factorial

def independent_alleles(input_file):
    '''
    (file) -> float
    Given 2 positive integers k and N from input file.
    In this problem, we begin with Tom, who in the 0th generation has genotype
    Aa Bb. Tom has two children in the 1st generation, each of whom has two
    children, and so on. Each organism always mates with an organism having
    genotype Aa Bb. Return the probability that at least N AaBb organisms will
    belong to the k-th generation of Tom's family tree 
    (don't count the AaBb mates at each level).    
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get k and N
    nums = infile.readline().rstrip().split()
    k = int(nums[0])
    N = int(nums[1])
    # close file
    infile.close()
    
    # the probability of having N AaBb offspring is given by the binomial
    # distribution P(k) = C(n,k) * p^k * q^(n-k)
    
    # p = probability that 1 offspring from AaBb x AaBb cross
    p = 0.25
    
    # finding n, the number of experiments in the binomial trial
    # n is the total number of offspring after k generations
    n = 2**k
    
    # finding k, the the total number of successes in the binomial trial
    # k is the number of offspring with genotype AaBb, k = N
    
    # because we want to know the probability of at least N offspring, 
    # we need to compute the probabilities of i to N-1 offpsring = P
    # proba(N) = 1-P
    
    P = 0    
    for i in range(N):
        P += (factorial(n) / (factorial(n - i) * factorial(i))) * p**i * (1-p)**(n-i)
    
    return round((1 - P), 3)
        
    
    
    
    
    
    
    
    
