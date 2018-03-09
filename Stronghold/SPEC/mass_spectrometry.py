# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 11:22:24 2015

@author: Richard
"""

def mass_spectrometry(input_file):
    '''
    (file) -> str
    Given a list L of n (n ≤ 100) of floats representing the prefix spectrum
    of a protein sequence, return a protein string of length n−1 whose
    prefix spectrum is equal to L (return any sequence if multiple exist)
    '''
       
    # mass table of each amino acid
    mass_table = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694,
                  'E': 129.04259, 'F': 147.06841, 'G': 57.02146,
                  'H': 137.05891, 'I': 113.08406, 'K': 128.09496,
                  'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
                  'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,
                  'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 
                  'W': 186.07931, 'Y': 163.06333}
                  
    # open file for reading
    infile = open(input_file, 'r')
    # create a list to hold the prefix spectrum
    weights = []
    for line in infile:
        line = line.rstrip()
        if line != '':
            weights.append(float(line))
    # close file
    infile.close()
    
    # The prefix spectrum of a protein is the collection of all
    # its prefix weights so the order of AA in the protein sequence can be
    # deduced from the weights of the prefices
    
    # check if the weights in prefix spectrum > weights of AA
    aa_weights = [weight for weight in mass_table.values()]
    heaviest_aa = max(aa_weights)
                
    # the prefix spectrum may contains weight that are very high, indicating
    # that the protein sequence is attached to something
    if heaviest_aa < weights[0]:
        additional_weight = True
    else:
        additional_weight = False
        
    # find the weights of the prefix
    peptide_weights = []
    if additional_weight:
        for i in range(1, len(weights)):
            peptide_weights.append(round(weights[i] - weights[i-1], 4))
    
    # weights of prefix reveal the order of the AA in the protein sequence:
    # weights[0] = weight_1st_aa + extra_weight
    # weights[1] = weight_1st_aa + weight_2nd_aa + extra_weight
    # weights[1] - weight[0] = weight_2nd_aa
    
    
    # initialize protein sequence
    protein  = ''
    # initialize position of the aa in protein sequence
    i = 0
    # potential issue: identical weights fr some amino acid, take any AA
    # make a reverse dictionnary mass : aa
    mass_AA = {}
    for aa in mass_table:
        if round(mass_table[aa], 4) in mass_AA:
            mass_AA[round(mass_table[aa], 4)].append(aa)
        else:
            mass_AA[round(mass_table[aa], 4)] = [aa]
    
    for weight in peptide_weights:
        for mass in mass_AA:
            if weight == mass:
                protein += mass_AA[mass][0]
    return protein
    
    
    
    
           
     
                  
                  
    