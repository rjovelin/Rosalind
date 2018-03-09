# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 20:54:41 2015

@author: Richard
"""

# create a dict of mass for each AA
mass_table = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
             'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
             'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156,
             'Y': 163, 'W': 186}
   
   
def cyclopeptide_theoritical_spectrum(P):
    '''
    (str) -> lst
    Return the theoritical spectrum of cyclic peptide.
    '''
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

   # return spectrum
    return spectrum   
   
 
def peptide_score(peptide, experimental_spectrum):
    '''
    (str, list) -> int
    Return a score being the number of subpeptide interger masses in common
    between the experimental spectrum and theoritical spectrum of peptide
    >>> peptide_score('NQEL', [0, 99, 113, 114, 128, 227, 257, 299,
    355, 356, 370, 371, 484])
    11
    '''
    # compute the theotical spectrum
    th_spectrum = cyclopeptide_theoritical_spectrum(peptide)
    # cget the elements in common
    score = len(set(th_spectrum).intersection(experimental_spectrum))
    
    return score


def peptide_mass(peptide):
    '''
    (str) -> int
    Return the integer mass of peptide
    '''
    mass = 0
    for aa in peptide:
        mass += mass_table[aa]
    return mass
    
    
    
    
    
    
    
   
   
#def theoritical_spectrum(peptide):
#    '''
#    (str) -> list
#    Return the theoritical spectrum of a linear peptide: a list with integer
#    masses of all subpeptides in peptide
#    >>> theoritical_spectrum('PEAL')  
#    '''
#    
#    spectrum = [] 
#    
#    for i in range(len(peptide)):
#        for j in range(1, len(peptide) -i+1):
#            # grab subpeptide
#            subpeptide = peptide[i:i+j]
#            # initialise mass_peptide
#            mass_peptide = 0
#            for aa in subpeptide:
#                mass_peptide += mass_table[aa]
#            spectrum.append(mass_peptide)
#            subs.append(subpeptide)
#    
#    # add 0 to the spectrum
#    spectrum.append(0)
#        
#    # sort spectrum
#    spectrum.sort()
#    # return spectrum
#    return spectrum, subs
    


def leaderboard_cyclopeptide_sequencing(N, exp_spectrum):
    '''
    (file) -> file
    Given an integer N and a collection of mass integers from an experimental
    Spectrum, return the sequence of the leader peptide in the form of amino
    acid mass-dased separated. 
    '''
        
    # make a list of AA
    # some AA have same mass, so collapse the A with same mass into a single AA
    AAs = []
    already_added = []
    for aa in mass_table.keys():
        if mass_table[aa] not in already_added:
            AAs.append(aa)
            already_added.append(mass_table[aa])
    print(len(AAs))
    
    # make a list of with (score, peptide) sublist
    leaderboard = []
       
    # initiate the leaderboard with empty peptides
    for i in range(len(AAs)):
        leaderboard.append([0, ''])
    
    # initalize the leaderpeptide
    leaderpeptide = ''
    leaderscore = 0
    
     # enter the branch and bound search, stop when the list is emoty
    while len(leaderboard) != 0:
        print('leaderboard', len(leaderboard))
        print('leaderpeptide', len(leaderpeptide))
        
        # create a list to store the branching amino acid
        branch_bound = []        
        
        # expand the AA by adding a single AA for all AA in AAs
        for i in range(len(leaderboard)):
            # loop over the AAs and add a single aa
            for amino in AAs:
                # inititialize the new aa
                branching_aa = ''               # create the branching aa
                branching_aa = leaderboard[i][1]  + amino
                # add the branching (score, aa) to branching_round
                branch_bound.append([peptide_score(branching_aa, exp_spectrum), branching_aa])
                
        # sort the peptides by decreasing score
        branch_bound.sort()
        branch_bound.reverse()
        print('branch_bound', len(branch_bound))
        
        # cut the peptides from the branching step by score rank
        # keep the N peptides, including ties
                        
        # check if there are ties
        
        # keep all if less values than stop otherwise cut any values with score below score at stop
        if len(branch_bound) > N:
            # find the minimum score to keep
            threshold = branch_bound[N][0]
        # find if there are ties
        for i in range(N, len(branch_bound)):
            if branch_bound[i] == threshold:
                # keep looking
                continue
            else:
                # reassign stop
                N = i
                # exit loop
                break
                
        # empty leaderboard
        leaderboard = []
        # cut the list of branching peptides and add to leaderboard
        if len(branch_bound) >= N:
            leaderboard.extend(branch_bound[:N])
        else:
            leaderboard.extend(branch_bound)
        
        # create a list to remove peptides with high mass
        to_remove = []
        
        # check the score and mass of each peptide
        for i in range(len(leaderboard)):
            # check the mass and score
            if leaderboard[i][0] > leaderscore and peptide_mass(leaderboard[i][1]) == exp_spectrum[-1]:
                # reassign leaderpeptide and leaderscore
                leaderpeptide = leaderboard[i][1]
                leaderscore = leaderboard[i][0]
            else:
                # remove peptides with mass > mass peptide
                if peptide_mass(leaderboard[i][1]) > exp_spectrum[-1]:
                    to_remove.append(leaderboard[i])
        # remove any (score, peptide) pair
        for pair in to_remove:
            leaderboard.remove(pair)
        
        print('leaderboard', len(leaderboard))
        leaderboard.reverse()
    
    # set up the final peptide to be returned
    sequenced_peptide = ''    
    
    for aa in leaderpeptide[:-1]:
        sequenced_peptide += str(mass_table[aa]) + '-'
    sequenced_peptide += str(mass_table[leaderpeptide[-1]])
    
    return sequenced_peptide
    
   

#if __name__ == '__main__':
#    input_file = input('enter name of input file: ')
#    # open file
#    infile = open(input_file, 'r')
#    N = int(infile.readline().rstrip())
#    spectrum = infile.readline().rstrip().split()
#    print(N)
#    print(spectrum)
#    # close file
#    infile.close()
#    # convert str to int
#    for i in range(len(spectrum)):
#        spectrum[i] = int(spectrum[i])
#    sequenced_peptide = leaderboard_cyclopeptide_sequencing(N, spectrum)