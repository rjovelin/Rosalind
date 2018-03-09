# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 20:12:21 2015

@author: Richard
"""

def odd_sum(input_file):
    '''
    (file) -> int
    Given 2 positive integers a and b (a<b<10000) in input file,
    return the sum of all odd integers from a through b, inclusively
    '''
    

    # open file
    infile = open(input_file, 'r')
    line = infile.readline().rstrip().split()
    a = int(line[0])
    b = int(line[1])
    infile.close()
    
    total = 0
    for i in range(a, b+1):
        if i % 2 != 0:
            total += i
    return total
    
    
    
    