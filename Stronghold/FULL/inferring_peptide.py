# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:34:11 2015

@author: Richard
"""


# THERE IS A BUG IN THE SCRIPT
# to pass the test: take the second peptide minus the last amino acid
# need to make sure a single peptide is returned or the reverse peptide,
# but without extra AA






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
        if round(mass_table[aa], 5) in mass_AA:
            mass_AA[round(mass_table[aa], 5)].append(aa)
        else:
            mass_AA[round(mass_table[aa], 5)] = [aa]
    return mass_AA


def find_highest_value(adjacent):
    '''
    (dict) -> num
    Return the highest value in the dictionnary of mass differences
    '''
    # find highest item value
    highest_val = 0
    for mass_differences in adjacent.keys():
        pair = mass_differences.split(':')
        pair[0] = float(pair[0])
        pair[1] = float(pair[1])
        if pair[0] > highest_val:
            highest_val = pair[0]
    
    return highest_val

def search_mass_differences(first_val, adjacent, highest_val):
    '''
    (list, dict, num) -> list, dict
    Return a list of peptide masses that are separated by the weight of a 
    single amino acid by searching the keys in the dictionnary adjacent
    starting with the highest value and until the highest value is found in the 
    list first_val, and return also the updated dictionnary by removing masses
    already used
    '''
    
    # create a list to store the mass peptides
    by_ions = []
    
    # create a list to store pairs of masses already used
    already_added = []
    
    # masses separated by the weight of a single AA are from the same ion type
    while highest_val in first_val:
        for mass_differences in adjacent.keys():
            pair = mass_differences.split(':')
            pair[0] = float(pair[0])
            pair[1] = float(pair[1])
            if highest_val == pair[0]:
                by_ions.append(pair[0])
                already_added.append(mass_differences)
                highest_val = pair[1]
    

    # add the smalles value
    for mass_differences in adjacent.keys():
        pair = mass_differences.split(':')
        pair[0] = float(pair[0])
        pair[1] = float(pair[1])
        if by_ions[-1] == pair[0]:
            by_ions.append(pair[1])
        
           
    # remove pairs of masses from mass difference dict if used
    for mass_differences in already_added:
        del adjacent[mass_differences]
        
    # return list of peptide masses and updated dict
    return by_ions, adjacent
    


def get_first_val(adjacent):
    '''
    (dict) -> list
    Return a list with the first masses of in the keys of the dictionnary of 
    mass differences adjacent
    '''
    
    # make a list of the first elements in adjacent.keys()
    first_val = []
    for mass_differences in adjacent.keys():
        pair = mass_differences.split(':')
        first_val.append(float(pair[0]))

    return first_val


    
def make_by_ion_list(highest_val, first_val, adjacent):
    '''
    (num, list, dict) -> list, list
    Return 2 lists of peptide masses corresponding to b-ions or y-ions
    '''
           
    # make list of b_ions (b_ions and y_ions lists representing the spectrum
    # of prefic and suffix, but without knowing which is the suffix or prefix)
    b_ions, adjacent = search_mass_differences(first_val, adjacent, highest_val)
    
    # find new highest value
    highest_val = find_highest_value(adjacent)
    
    # make list of y_ions (b_ions and y_ions lists representing the spectrum
    # of prefic and suffix, but without knowing which is the suffix or prefix)
    y_ions, adjacent = search_mass_differences(first_val, adjacent, highest_val)
    
    # sort lists to obtain the prefix and suffix spectrum
    b_ions.sort()
    y_ions.sort()
    
    return b_ions, y_ions
    

def inferring_peptide(input_file):
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
    
    # parent mass is first mass in weights
    parent = weights[0]
    
    
    
    # length of internal peptide is n 
    n = int((len(weights) / 2) -3)
    
    # impossible to dtermine the b-ion from the y-ion, only pairs on ions
    # create a dictionnary to store the ion pairs
    spectrum = {}
    for i in range(len(weights)):
        print(round(parent - weights[i], 5))
        if round(parent - weights[i], 5) in weights:
            spectrum[weights[i]] = round(parent - weights[i], 5)
    
    # create a list to store the ion masses
    ions = [mass for mass in spectrum]
    
    # create a reverse dict to loop up for aa corresponding to given masses
    mass_AA = reverse_mass_table(mass_table)    
    
    # create a dictionnary to store adjacent masses
    adjacent = {}
    
    for item in ions:
        for element in ions:
            # round each difference to 5 points after decimal
            if round((item - element), 5) in mass_AA:
                adjacent[str(item) + ':' + str(element)] = round((item - element), 5)
    
    # make a list of the first elements in adjacent.keys()
    first_val = get_first_val(adjacent)
    
    # find highest value
    highest_val = find_highest_value(adjacent)
    
    
    # find the prefix and suffix spectrum
    b_ions, y_ions = make_by_ion_list(highest_val, first_val, adjacent) 
    
    # get the amino masses by substracting masses in prefic and suffix spectrum
    b_weights = []
    for i in range(1, len(b_ions)):
        b_weights.append(round(b_ions[i] - b_ions[i-1], 5))
    y_weights = []
    for i in range(1, len(y_ions)):
        y_weights.append(round(y_ions[i] - y_ions[i-1], 5))
    
    # build peptides using the amino acid masses
    b_protein = ''
    for weight in b_weights:
        for mass in mass_AA:
            if weight == mass:
                b_protein += mass_AA[mass][0]
                
    y_protein = ''
    for weight in y_weights:
        for mass in mass_AA:
            if weight == mass:
                y_protein += mass_AA[mass][0]
                
    return b_protein, y_protein
        
    
    
    
    
    
      
    
 
    
    # weights of prefix reveal the order of the AA in the protein sequence:
    # weights[0] = weight_1st_aa + extra_weight
    # weights[1] = weight_1st_aa + weight_2nd_aa + extra_weight
    # weights[1] - weight[0] = weight_2nd_aa
    
    
 
