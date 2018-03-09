# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:30:38 2015

@author: Richard
"""

import math

def disease_carriers(input_file):
    '''
    (file) -> list
    Given an array A in input_file for which A[k] represents the
    proportion of homozygous recessive individuals for a given gene in a 
    in a diploid population at Hardy-Weinberg equilibrium, return an array B
    having the same length as A in which B[k] represents the probability
    that a randomly chosen individual carries at least one copy of
    the recessive allele   
    '''
        
    # Hardy-Weinberg equilibrium: p2 + 2pq + q2 = 1
    
    # probability that a individual carries at least 1 recessive allele for the
    # gene of interest if (1 - probability individual is homozygous dominant)
    # Prob = 1 - (1 - q)2
    
    # open file for reading
    infile = open(input_file, 'r')
    # create list of recessive homozygous frequencies q2
    q_deux = infile.readline().rstrip().split()
    # convert str to floats
    for i in range(len(q_deux)):
        q_deux[i] = float(q_deux[i])
    # close file 
    infile.close()
    
    recessif = []
    for q2 in q_deux:
        recessif.append(1 - (1 - math.sqrt(q2))**2)
    
    # print Prob to specified format
    for P in recessif:
        print(round(P, 3), end = ' ')
    
    
    
    