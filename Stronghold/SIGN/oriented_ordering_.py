# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 18:44:20 2015

@author: Richard
"""

import itertools

def oriented_ordering(input_file, outputfile):
    '''
    (file) -> file
    Grab the postive integer n (n<=6) from the input file and save
    to outputfile the total number of signed permutations of length n
    followed by a list of all such permutations (in any order)    
    '''
    
    
    # open file for reading
    infile = open(input_file, 'r')
    # grab n
    n = int(infile.readline().rstrip())
    # close file
    infile.close()
    
    # make a list of signed elements for {1, 2..n}
    elements = []
    for i in range(1, n+1):
        elements.append(i)
        elements.append(-i)
    
    # make a list of all permutations of elements of length n
    all_perm = list(itertools.permutations(elements, n))
    
    # exclude tuples with positive and negative identical int
    signed_perm = []
    check_set = set(i for i in range(1, n+1))
    for item in all_perm:
        to_check = {abs(i) for i in item}
        if to_check == check_set:
            signed_perm.append(item)
    
    # convert each tuple into a list
    for i in range(len(signed_perm)):
        signed_perm[i] = list(signed_perm[i])
    
    
    # convert elements of each tuple to string
    for item in signed_perm:
        for i in range(len(item)):
            item[i] = str(item[i])
    
    # save permutations to file
    newfile = open(outputfile, 'w')
    
    newfile.write(str(len(signed_perm)) + '\n')    
    
    # make a string with the inner list of signed_perm
    for inner_list in signed_perm:
        perm = ' '.join(inner_list)
        newfile.write(perm + '\n')
    
    # close file
    newfile.close()
    
    
    