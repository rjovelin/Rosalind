# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 11:07:27 2015

@author: Richard
"""


def matching_spectrum(input_file):
    '''
    (file) -> dict
    Given a positive integer n followed by a collection of n protein strings
    s1, s2 ,..., sn and a multiset R of positive numbers (corresponding to
    the complete spectrum of some unknown protein string), return the maximum
    multiplicity of R⊖S[sk] (spectral convolution between R and each string)
    taken over all strings sk, followed by the string sk for which this
    maximum multiplicity occurs
    (return any value, string if multiple solutions exist) 
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
    n = int(infile.readline().rstrip())
    # make a list of string peptides
    S1 = []
    while n!= 0:
        S1.append(infile.readline().rstrip())
        n -= 1
    # make a list to store the complete spectrum of the unknown protein
    S2 = []
    for line in infile:
        line = line.rstrip()
        if line != '':
            S2.append(line)
            
    # close file after reading
    infile.close()
    
    # combert str to floats
    for i in range(len(S2)):
        S2[i] = float(S2[i])
        
    # createa dictionnary to store the spectrum of each peptide
    spectrum = {}
    for peptide in S1:
        spectrum[peptide] = []
        mass = 0
        for aa  in peptide:
            mass += mass_table[aa]
            spectrum[peptide].append(mass)

    # compute the Minkowski dofference between S1 and each peptide spectrum
    convolution = {}
    for peptide in spectrum:
        convolution[peptide] = []
        for item in S2:
            for element in spectrum[peptide]:
                convolution[peptide].append(round(item - element, 5))
    
    # make a dictionnary of dictionnaries with the count: mass pairs
    counter = {}
    for peptide in convolution:
        counter[peptide] = {}
        for element in convolution[peptide]:
            if element in counter[peptide]:
                counter[peptide][element] += 1
            else:
                counter[peptide][element] = 1
    
    # sort lists in the dictionnary
    spectral_convolution = {}
    for peptide in counter:
        spectral_convolution[peptide] = [(count, element) for element, count in counter[peptide].items()]
    
    for peptide in spectral_convolution:
        spectral_convolution[peptide].sort()
    
    # find the maximum multiplicity
    multiplicity = 0
    protein = ''    
    
    for peptide in spectral_convolution:
        # check the last element of the list of each peptide
        # it contains the mass difference with highest frequency (multiplicity)
        if spectral_convolution[peptide][-1][0] > multiplicity:
            multiplicity = spectral_convolution[peptide][-1][0]
            protein = peptide
    
    print(multiplicity, protein, sep = '\n')

    
