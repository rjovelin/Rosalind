# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 09:45:45 2015

@author: Richard
"""

from scipy import misc

def expected_occurences(input_file, outputfile):
    '''
    (str, float) -> float
    Given a string s and a GC content GC, return the probability
    that a random string with GC content equal to GC is identical to s
    '''
    
    
    # open file for reading
    infile = open(input_file, 'r')
    n = int(infile.readline().rstrip())
    s = infile.readline().rstrip().upper()
    A = infile.readline().rstrip().split()
    for i in range(len(A)):
        A[i] = float(A[i])
    # close file
    infile.close()
    
    # create list to store the expected number of occurences
    occurences = []
    
    # determine the probabilities of nucleotide given GC content
    for GC in A:
        # compute the probabilities of each nucleotide given GC content
        pCG = GC / 2
        pAT = (1- GC) / 2
        nucleotides = {'A':pAT, 'T': pAT, 'C':pCG, 'G':pCG}
        # compute the probability of the string s
        Ps = 1
        for i in range(len(s)):
            Ps *= nucleotides[s[i]]
        
        # compute the expected number of occurences
        # maximum occurences = n-len(s)+1
        Exp = round(Ps * (n - len(s)+1), 3)
        occurences.append(Exp)
    
    # open file for writing
    newfile = open(outputfile, 'w')
    for item in occurences[:-1]:
        newfile.write(str(item) + ' ')
    newfile.write(str(occurences[-1]) + '\n')
    # close file
    newfile.close()
    
        
            
