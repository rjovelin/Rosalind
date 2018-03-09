# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 16:19:39 2015

@author: Richard
"""

def sex_linked_genes(input_file):
    '''
    (file) -> list
    Given an array A of length n for which A[k] represents the proportion
    of males in a population exhibiting the recessive phenotype for the 
    the k-th of n total recessive X-linked genes. 
    Return an array B of length n  in which B[k] equals the probability that a
    randomly selected female will be a carrier for the k-th gene.
    Precondition: the population is in genetic equilibrium for all n genes
    (ie. Assume Hardy-Weinberg equilibrium)
    '''
    # open file for reading
    infile = open(input_file, 'r')
    
    # create list of males proportions for the X-linked genes
    males = infile.readline().rstrip().split()
    # convert str to float
    for i in range(len(males)):
        males[i] = float(males[i])
    # close file
    infile.close()
    
    # proportions of makes with the phenotype of a recessive X-linked gene
    # gives q the frequency of the recessive allele
    
    # NB: a carrier has the recessive alelles but does not show phenotype
    # female_carrriers are heterozygous
    
    # assuming HW equilibrium, we now p and q for the entire population
    # and the frequency of females is simply 2pq
    
    # create list to store the frequency of fmeale carriers
    females = []
    for i in range(len(males)):
        females.append(round(2 * males[i] * (1 - males[i]), 3))
        
    # print to expected format
    for carrier_freq in females:
        print(carrier_freq, end = ' ')
        
        
        