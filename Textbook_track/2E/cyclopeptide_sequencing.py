# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 00:49:08 2015

@author: Richard
"""


# implementing a branch-and-bound for finding the peptide of a theoritical spectrum
def cyclopeptide_sequencing(input_file, outputfile):
    '''
    (file) -> file
    Given a collection of (possibly repeated) integer masses corresponding
    to an ideal experimental spectrum in the input file, return an amino acid
    string Peptide such that Cyclospectrum(Peptide) = Spectrum.
    Return the peptide strings by the mass of the AA seperated by a dash    
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # make a list with the theoritcal spectrum
    spectrum = infile.readline().rstrip().split()
    # convert str to int
    for i in range(len(spectrum)):
        spectrum[i] = int(spectrum[i])
    # close file
    infile.close()
        
    # create a dict of mass for each AA
    mass_table = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
                  'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
                  'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156,
                  'Y': 163, 'W': 186}
    
    # make a list of AA
    # some AA have same mass, so collapse the A with same mass into a single AA
    AAs = []
    already_added = []
    for aa in mass_table.keys():
        if mass_table[aa] not in already_added:
            AAs.append(aa)
            already_added.append(mass_table[aa])
        
    # create a list to store the starting amino acids
    potential_peptides = []
    
    # find the spectrum masses found in the mass table and add to list
    # look in the restricted sets of AAs
    for aa in AAs:
        if mass_table[aa] in spectrum:
            potential_peptides.append(aa)
        
    # create a list to store the correct peptides
    sequenced_peptides = []
    
    
    # enter the branch and bound search, stop when the list is emoty
    while len(potential_peptides) != 0:
        print('potential', len(potential_peptides))
        print('sequenced', len(sequenced_peptides))
        
        # create a list to store the branching amino acid
        branch_bound = []        
        
        # expand the AA by adding a single AA for all AA in AAs
        for aa in potential_peptides:
            # loop over the AAs and add a single aa
            for amino in AAs:
                # inititialize the new aa
                branching_aa = ''
                # create the branching aa
                branching_aa = aa + amino
                # add the branching aa to branching_round
                branch_bound.append(branching_aa)
        
        # remove all source peptides for the branching round from list
        potential_peptides = []        
        
        # merge the branching peptides with the potential peptides
        potential_peptides.extend(branch_bound)
        print('potential', len(potential_peptides))
        
        # reset the branching list to store peptides after bound
        branch_bound = []
        
        # check that the AA masses and the peptide mass are in spectrum
        for peptide in potential_peptides:
            # set up boolean
            correct_peptide_mass = True
            
            # the mass of each subpeptide must be in spectrum
            for i in range(len(peptide)):
                for j in range(1, len(peptide) -i+1):
                    # grab subpeptide
                    subpeptide = peptide[i:i+j]
                    # initialise mass_peptide
                    mass_peptide = 0
                    for aa in subpeptide:
                        mass_peptide += mass_table[aa]
                    # check that mass subpeptide in spectrum
                    if mass_peptide not in spectrum:
                        # change boolean value
                        correct_peptide_mass = False
                        # exit loop
                        break
            
            # set boolean
            correct_aa_mass = True
            # evaluate the mass of each single aa in the peptide
            for aa in peptide:
                if mass_table[aa] not in spectrum:
                    # change boolean
                    correct_aa_mass = False
                    # exit loop                        
                    break
                                    
            # evaluate that aa masses and peptide mass are in spectrum
            if correct_aa_mass and correct_peptide_mass:
                # reinitialise peptide_mass
                peptide_mass = 0
                for i in range(len(peptide)):
                    peptide_mass += mass_table[peptide[i]]
                # add to sequenced list if peptide_mass is entire peptide mass
                if peptide_mass == spectrum[-1]:
                    sequenced_peptides.append(peptide)
                else:
                    # if not, do nothing, repeat branch and bound
                    branch_bound.append(peptide)
              
        # remove all peptides
        potential_peptides = []
        # add potential peptides, if any, after the bound step
        potential_peptides.extend(branch_bound)
        print('potential', len(potential_peptides))
        
 
    # open file for writing
    newfile = open(outputfile, 'w')
    # save each peptide as a string of AA mass-dash separated
    if len(sequenced_peptides) != 0:
        for peptide in sequenced_peptides:
            for aa in peptide[:-1]:
                newfile.write(str(mass_table[aa]) + '-')
            newfile.write(str(mass_table[peptide[-1]]) + ' ')
    # close file
    newfile.close()
    
    
    
    