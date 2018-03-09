# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 11:57:39 2015

@author: Richard
"""

# align sequences with ClustalX from http://www.ebi.ac.uk/Tools/msa/clustalw2/
# save the identity matrix

import numpy as np

def most_divergent_sequence(input_file):
    '''
    (file) -> str
    Return the name of the most divergent sequence using the identity matrix
    output from the ClustalW alignment tol of EBI
    '''
    
    # make a dictionnary to store the seq_name: list of identity scores
    identity = {}
    
    # open file for reading
    infile = open(input_file, 'r')
    for line in infile:
        line = line.strip()
        if line != '':
            if not line.startswith('#'):
                line = line.split()
                identity[line[1]] = line[2:]
                
    # convert str to floats
    for name in identity:
        for i in range(len(identity[name])):
            identity[name][i] = float(identity[name][i])
            
    # change identity dict to seq_name: mean identity score pairs
    for name in identity:
        identity[name] = np.mean(identity[name])
        
    # find the sequence with the lowest score
    # initiate score for comparison
    names = [name for name in identity]
    lowest_score = identity[names[0]]
    divergent = names[0]
    
    for name in identity:
        if identity[name] < lowest_score:
            lowest_score = identity[name]
            divergent = name
    
    return divergent
    
    
    
                
                
    
        
    
    
    
    
    
    
    
    
    
