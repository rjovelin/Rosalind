# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 09:07:44 2015

@author: Richard
"""

def approximate_pattern(input_file):
    '''
    (file) - > list
    Grab pattern P,  sequence S and integer d from the input file and 
    return all starting positions of P in S whith at most d mismatches
    '''
    
    # open file
    infile = open(input_file, 'r')
    line1 = infile.readline().rstrip()
    line2 = infile.readline().rstrip()
    line3 = infile.readline().rstrip()
    # close file
    infile.close()
    # get the values of each variable
    d = int(line3)
    if len(line1) < len(line2):
        P, S = line1, line2
    else:
        P, S = line2, line1
    
    # create a list to record the matching positions
    positions = []    
    
    
    # go through the sequence nucleotide by nucleotide,
    # scan each motif of length P and record position if d <= 3
    for i in range(len(S) - len(P) +1):
        motif = S[i: i + len(P)]
        if motif == P:
            positions.append(i)
        else:
            diff = 0
            for j in range(len(P)):
                if P[j] != motif[j]:
                    diff += 1
            if diff <= d:
                if i not in positions:
                    positions.append(i)
    
    # print the positions to screen
    if len(positions) != 0:
        for item in positions:
            print(item, end = ' ')
    else:
        print('S has no matches for P')
    
    
        
        
    
    
    