# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 23:03:41 2015

@author: Richard
"""

import scipy.misc as msc

def counting_subsets(input_file):
    '''
    (file) -> int
    Grab the positive integer N (N <= 1000) from the imput_file
    Print the total number of subsets of {1,2,…,N} modulo 1,000,000
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # grab N
    N = int(infile.readline().rstrip())
    # close file
    infile.close()
    
    # set up a counter fot the total number of subsets
    all_subsets = 0
    
    # count the number of possible combinations N choose i
    for i in range(1, N + 1):
        all_subsets += msc.comb(N, i, True)
    
    # add the empty set as part of the subset of the set of N
    all_subsets += 1
    
    # get the modulo
    all_subsets = all_subsets % 1000000
    
    print(all_subsets)
    
    
    