# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 23:50:37 2015

@author: Richard
"""

# mass table of each amino acid as global variable
mass_table = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694,
              'E': 129.04259, 'F': 147.06841, 'G': 57.02146,
              'H': 137.05891, 'I': 113.08406, 'K': 128.09496,
              'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
              'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,
              'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 
              'W': 186.07931, 'Y': 163.06333}


def reverse_mass_table(mass_table):
    '''
    (dict) -> dict
    Return the reverse dictionnary of mass_table with masses as keys and 
    amino acid as values   
    '''
    # create a reverse dictionnary to loop up masses
    mass_AA = {}
    for aa in mass_table:
        if round(mass_table[aa], 2) in mass_AA:
            mass_AA[round(mass_table[aa], 2)].append(aa)
        else:
            mass_AA[round(mass_table[aa], 2)] = [aa]
    return mass_AA


def find_lowest_spectral_value(spectrum):
    '''
    (dict) -> num
    Search the mass differences as string in the dictionnary spectrum keys
    and return the lowest mass value   
    '''
    
    # make a list to store all masses in spectrum keys
    masses = []
        
    # find lowest value
    for mass_differences in spectrum:
        # get the masses connected by a a single AA
        pair = mass_differences.split(':')
        # convert the str masses to floats
        for i in range(len(pair)):
            pair[i] = float(pair[i])
        for mass in pair:
            masses.append(mass)
    # return the lowest mass value
    return min(masses)


def real_spectra(input_file):
    '''
    (file) -> str
    Given a list L of n (n ≤ 100) of floats representing the prefix spectrum
    of a protein sequence, return a protein string of length n−1 whose
    prefix spectrum is equal to L (return any sequence if multiple exist)
    '''
       
    # open file for reading
    infile = open(input_file, 'r')
    # create a list to hold the prefix spectrum
    weights = []
    for line in infile:
        line = line.rstrip()
        if line != '':
            weights.append(round(float(line), 5))
    # close file
    infile.close()
    
    # convert str into floats
    for i in range(len(weights)):
        weights[i] = float(weights[i])
     
    # make a reverse dict to look up the AA for each mass    
    mass_AA = reverse_mass_table(mass_table)
    
    # create a link between node2 and node1 in weights if node2 > node1 
    # and node2 - node1 in mass_AA, label the edge with the AA
    
    spectrum = {}
    for i in range(len(weights)):
        for j in range(len(weights)):
            if weights[j] > weights[i]:
#                print(round(weights[j] - weights[i], 2))
                if round(weights[j] - weights[i], 2) in mass_AA:
#                    spectrum[str(weights[j]) + ':' + str(weights[i])] = round(weights[j] - weights[i], 2)
                    spectrum[str(weights[j]) + ':' + str(weights[i])] = mass_AA[round(weights[j] - weights[i], 2)][0]
    
    # create a list of peptides
    peptides = []
    
    # loop until spectrum is empty, store all potential peptides in list
    while len(spectrum) != 0:
        # find lowest value
        lowest_val = find_lowest_spectral_value(spectrum)
                
        # initiate a protein string
        protein = ''
        
        # create a list of lower masses in mass differences keys
        lower_masses = []
        for mass_differences in spectrum:
            pair = mass_differences.split(':')
            for i in range(len(pair)):
                pair[i] = float(pair[i])
            lower_masses.append(pair[1])
        
        # loop until lowest val is in the list of lower masses
        while lowest_val in lower_masses:
            # find the lowest value in mass differences, add the AA to protein
            # re-assign lowest value to connected mass
            for mass_differences in spectrum:
                # get the masses connected by a a single AA
                pair = mass_differences.split(':')
                # convert the str masses to floats
                for i in range(len(pair)):
                    pair[i] = float(pair[i])
                if lowest_val == pair[1]:
                    # update protein
                    protein += spectrum[mass_differences]
                    # re-assign lowest_val
                    lowest_val = pair[0]
                    # assign key to be removed when loop exists
                    to_remove = mass_differences
                  
        # remove last link if the end of peptode reached
        # (ie when lowest val not in lowest masses and while loop has exited)
        del spectrum[to_remove]
        
        # add protein to peptides list
        peptides.append(protein)
    
    # find longest protein
    protein = peptides[0]
    longest = len(protein)
    for peptide in peptides:
        if len(peptide) > longest:
            protein = peptide
            longest = len(peptide)
    
    return protein
                    
            

