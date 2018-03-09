# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 00:31:56 2015

@author: Richard
"""

import itertools



def fact(n):
    '''
    (int) -> int
    Return the factorial of n
    '''
    if n == 0 or n == 1:
        return 1
    else:
        return n * fact(n - 1)


def partial_permutations(input_file):
    '''
    (file) -> int
    Grab the positive integers n (0 < n <=100) and k (0 < k <= 10)
    from the imput file and return the total number of partial
    permutations P(n,k) modulo 1,000,000
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    line = infile.readline().rstrip().split()
    n = int(line[0])
    k = int(line[1])
    
    # close file
    infile.close()
    
#    # slow way to do it:    
#    # make a generator of with the number of permutations of size k from n
#    perm = itertools.permutations(range(1, n+1), k)
#    # count the number of items in the genrator
#    # stop counting when sentinel is reached
#    i = 0
#    total = 0
#    while i != -1:
#        total += 1
#        i = next(perm, -1)
#    return total % 1000000
    
    
    # fast way is tp use factorials and the definiton of permutations
    return int((fact(n) / fact(n - k)) % 1000000)
    
    
