# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 22:55:41 2015

@author: Richard
"""

def peptides_with_given_mass_counter(input_file):
    '''
    (file) -> int
    Given an integer m in input file representing the integer mass of an
    unknown peptide, return the number of linear peptides with mass m
    '''
    
    # create a dict of mass for each AA
    mass_table = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
                  'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
                  'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156,
                  'Y': 163, 'W': 186}
    
    # open file for reading
    infile = open(input_file, 'r')
    # grab peptide mass
    m = int(infile.readline().rstrip())
    # close file
    infile.close()
    
    # use dynamic programming to count the number of possible peptides
    # with given mass
   
    # create a dictionnary to store the mass : N_peptide pairs
    peptides_mass = {}
    
    # initiate dictionnary
    peptides_mass[0] = 0
    
    # iterate from 1 to m inclusive
    for i in range(1, m+1):
        # create a variable to count the number of possible AAs
        # that can be added to 
    
    
    
    # get the minimum and maximum length of the peptide
    minimum_length = int(m // 186) # minimum number of the heaviest aa
    maximum_length = math.ceil(m / 57) # maximum number of the lightest
    
    
    # create string with all amino acids
    amino_acids = ''.join(mass_table.keys())
    
    # generate all possible peptides, evaluate their mass and count if m
    for i in range(minimum_length, maximum_length + 1):
        print(i)
        # mass doesn't dpend on the order of the aa
        peptides = itertools.combinations_with_replacement(amino_acids, i)
        # check all peptides
        for seq in peptides:
            # initialise the peptide's mass
            mass = 0
            # compute mass by looking up in the mass table
            for aa in seq:
                mass += mass_table[aa]
            # accumulate the peptide counter if mass = m
            if mass == m:
                peptide_m += 1
    
    return peptide_m
                  
                  
                  
                  
    
if __name__ == '__main__':
    import itertools
    import math