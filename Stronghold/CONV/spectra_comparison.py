# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 14:01:41 2015

@author: Richard
"""

def spectral_convolution(input_file):
    '''
    (file) -> tuple
    Given 2 lists of positive real numbers S1 and S2, representing 2 simplified
    spectra (length <= 200), return the largest multiplicity of the spectral convolution 
    between S1 and S2 as well as the absolute value of the number x having the
    largest multiplicity (return any value if multiple solutions exist)
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get S1 and S2
    S1 = infile.readline().rstrip().split()
    S2 = infile.readline().rstrip().split()
    
    # close file after reading
    infile.close()
    
    # convert str into float
    for i in range(len(S1)):
        S1[i] = float(S1[i])
    for i in range(len(S2)):
        S2[i] = float(S2[i])
    
    # compute the Minkowski dofference between S1 and S2
    S3 = []
    for item in S1:
        for element in S2:
            # round each difference to 5 points after decimal
            S3.append(round((item - element), 5))
    
    # make a dictionnary with the mass : count pairs
    counter = {}
    for element in S3:
        if element in counter:
            counter[element] += 1
        else:
            counter[element] = 1
    
    # make a list to store tuples (count, element)
    convolution = [(count, element) for element, count in counter.items()]
    
    # sort the convolution according to counts
    convolution.sort()
    
    # get the tuple (count, element) with largest count
    multiplicity, element = convolution[-1]
    
    print(multiplicity, element, sep ='\n')
        
    
    
    