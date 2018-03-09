# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 02:22:04 2015

@author: Richard
"""

# note: cannot use factorial to compute permutations and combinations for large numbers
# instead use the module scipy.misc with functions to do just that


def alternative_splicing(input_file):
    '''
    (file) -> int
    Grab the positive integers n and m from the input file
    (0 <= m <= n <= 2000) and return the sum of combinations C(n,k)
    for all k satisfying m <= k <= n modulo 1,000,000.
    In shorthand, ∑nk=m(nk).
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # grab n and m
    line = infile.readline().rstrip().split()
    n = int(line[0])
    m = int(line[1])
    # close file
    infile.close()
    
    # compute C(n, k) for k in {m, ..., n} and return sum % 1000000
    total = 0
    for k in range(m, n+1):
        total += msc.comb(n, k, exact = True)
        
    return int(total % 1000000)

if __name__ == '__main__':
    import scipy.misc as msc