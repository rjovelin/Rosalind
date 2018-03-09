# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 17:24:45 2015

@author: Richard
"""

def theoritical_spectrum(input_file):
    '''
    (file) -> lst
    Get the peptide sequence P from the input file and return the theoritical
    spectrum of cyclic P.
    Note1: theorical spectrum of a cyclic peptide (Cyclospectrum(Peptide)),
    is the collection of all of the masses of its subpeptides, in addition
    to the mass 0 and the mass of the entire peptide.
    Note2: assume that duplicate values are tolerated
    '''
    # create a dict of mass for each AA
    mass_table = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
                  'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
                  'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156,
                  'Y': 163, 'W': 186}
    
    # open file
    infile = open(input_file, 'r')
    P = infile.readline().rstrip()
    infile.close()
    
    # make a list to store the mass of each subpeptide
    spectrum = []
    # add mass 0
    spectrum.append(0)
    # create a peptide used to loop over as it were cyclic
    PP = P + P
    # create a list to hold all subpeptides
    subpeptides = []
    # generate subpeptides up to the length of P, and cycling through PP
    for i in range(len(P)):
        for j in range(1, len(P)):
            subpeptides.append(PP[i:i+j])
            
    # add the entire peptide
    subpeptides.append(P)
    
    # compute mass for each subpeptide and add to the spectrum list
    for peptide in subpeptides:
        # initialise counter
        mass = 0
        # look up mass of each AA and accumulate peptide mass, store in list
        for AA in peptide:
            mass += mass_table[AA]
        spectrum.append(mass)
    
    # sort list    
    spectrum.sort()    
    
    # print spectrum
    for mass in spectrum:
        print(mass, end = ' ')