# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 08:39:06 2015

@author: Richard
"""

def spectral_convolution(spectrum):
    '''
    (list) -> list
    Compute the spectral convolution of spectrum and return a list of elements
    in the convolution
    >>> spectral_convolution([0, 137, 186, 323])
    [137, 186, 49, 323, 186, 137]
    '''
    
    # create a dictionnary to hold the elements on spectum convolution
    convolution = {}
    # initialize i
    i = 1
    # compute the elements of spectrum convolution as absolute value
    # of each 2 elements that are different: ie: compute a matrix of differences
    while i != len(spectrum):
        # initialize value of i
        convolution[i] = []
        for k in range(i):
            convolution[i].append(abs(spectrum[i] - spectrum[k]))
        i += 1
    # create  a list to store the elements of the spectral convolution
    spectral = []
    # add all elements to list
    for i in convolution:
        for item in convolution[i]:
            # add only item > 0
            if item > 0:
                spectral.append(item)
    
    return spectral


def multiplicity_spectral_convolution(spectrum):
    '''
    (list) -> list
    Return a list of elements from the convolution of spectrum ordered in
    decreasing order of their multiplicity. Return each element the number
    of times it appears in the convolution
    >>> multiplicity_spectral_convolution([0, 137, 186, 323])
    [137 137 186 186 323 49]
    '''
    
    # compute the convolution of spectrum
    convolution = spectral_convolution(spectrum)
    
    # make a dictionnary with element : count pairs
    multiplicity = {}
    for item in convolution:
        multiplicity[item] = convolution.count(item)
    
    # create a list of tuples (count, item) for item in multiplicity
    counts = [(count, item) for item, count in multiplicity.items()]
    # reverse sort the list
    counts.sort()
    counts.reverse()
    
    # create a list to store all elements by decreasing order
    ordered = []
    for pair in counts:
        j = pair[0]
        while j != 0:
            ordered.append(pair[1])
            j -= 1
    
    return ordered
        
    
def M_most_frequent(spectrum, M):
    '''
    (list, int) -> list
    Return the M most frequent elements between 57 and 200 from the
    spectral convolution of spectrum   
    '''
    
    # get the spectral convolution, with elements ordered by their multiplicity
    convolution = multiplicity_spectral_convolution(spectrum)
    
    # remove elements that are <= 57 or >= 200
    to_remove = []
    for element in convolution:
        if element < 57:
            to_remove.append(element)
        elif element > 200:
            to_remove.append(element)
            
    for element in to_remove:
        convolution.remove(element)
    
    # make a set of elements
    element_set = set(convolution)
    
    # create a list of most frequent element
    most_frequent = []
    
    # if N element < M: take all
    if len(element_set) < M:
        for element in element_set:
            most_frequent.append(element)
    else:
        # take the M most frequent element
        # make a dict to store element : count
        element_count = {}
        for element in convolution:
            element_count[element] = convolution.count(element)
        # make list of (count,element)
        counts_elements = [(count, element) for element, count in element_count.items()]
        # sort dict
        counts_elements.sort()
        # reverse sort to get the most frequent at the begining of list
        counts_elements.reverse()
        # take the M most frequent pairs
        
        # check if there are ties
        # find the minimum count to keep
        threshold = counts_elements[M-1][0]
        # find if there are ties
        for i in range(M-1, len(counts_elements)):
            if counts_elements[i][0] == threshold:
                # keep looking
                continue
            else:
                # reassign stop
                M = i
                # exit loop
                break
        
        M_counts = counts_elements[0:M]
        for pair in M_counts:
            most_frequent.append(pair[1])
    
    return most_frequent


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
    
    # create a list from P, removing the '-' character
    P = P.split('-')     
    
    # create a peptide used to loop over as it were cyclic
    PP = P + P
    # create a list to hold all subpeptides
    subpeptides = []
    # generate subpeptides up to the length of P, and cycling through PP
    for i in range(len(P)):
        for j in range(1, len(P)):
            subpeptides.append('-'.join(PP[i:i+j]))
            
    # add the entire peptide
    subpeptides.append('-'.join(P))
    
    # compute mass for each subpeptide and add to the spectrum list
    for peptide in subpeptides:
        # initialise counter
        mass = 0
        # look up mass of each AA and accumulate peptide mass, store in list
        mass_peptide = peptide.split('-')
        for AA in mass_peptide:
            mass += int(AA)
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
    >>> peptide_score('114-128-129-113', [0, 99, 113, 114, 128, 227, 257, 299,
    355, 356, 370, 371, 484])
    11
    '''
    # compute the theotical spectrum
    th_spectrum = cyclopeptide_theoritical_spectrum(peptide)
    # cget the elements in common
    score = len(set(th_spectrum).intersection(set(experimental_spectrum)))
#    score = 0
#    common = set(th_spectrum).intersection(set(experimental_spectrum))
#    for item in common:
#        score += experimental_spectrum.count(item)
    
    return score


def peptide_mass(peptide):
    '''
    (str) -> int
    Return the integer mass of peptide
    '''
    
    # make a list of mass from peptide
    mass_peptide = peptide.split('-')
    # convert str to int
    for i in range(len(mass_peptide)):
        mass_peptide[i] = int(mass_peptide[i])
    
    # add the masses in peptide
    mass = 0
    for aa in mass_peptide:
        mass +=aa
    return mass
    
    


def convolution_leaderboard_cyclopeptide_sequencing(input_file):
    '''
    (file) -> file
    Given an integer N, an experimental spectrum of mass integers from the
    input file,  return the sequence of the leader peptide in the form of
    amino acid mass-dased separated. 
    '''
        
    # open file
    infile = open(input_file, 'r')
    M = int(infile.readline().rstrip())
    N = int(infile.readline().rstrip())
    exp_spectrum = infile.readline().rstrip().split()
    print('M', M)
    print('N', N)
    print(exp_spectrum)
    # close file
    infile.close()
    # convert str to tint
    for i in range(len(exp_spectrum)):
        exp_spectrum[i] = int(exp_spectrum[i])
     
    # add 0 if 0 not in spectrum
    if 0 not in exp_spectrum:
        exp_spectrum.append(0)
    # sort spectrum before taking the convolution
    exp_spectrum.sort()
    
    # get the most frequemtn elements from the experimental spectrum    
    frequent = M_most_frequent(exp_spectrum, M)
    
    spectral = spectral_convolution(exp_spectrum)
    
    # define an alphabet from the most frequent elements
    # directly use the masses to build the peptides, not the AAs since some
    #  AAs are not in the mass_table
    AAs = []
    for mass in frequent:
        AAs.append(str(mass))
    print(len(AAs))
    
    # sort alphabet
    AAs.sort()
    
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
                branching_aa = leaderboard[i][1] + '-' + amino
                # remove '-' in 1st position
                if branching_aa[0] == '-':
                    branching_aa = branching_aa[1:]
                
                # add the branching (score, aa) to branching_round
                branch_bound.append([peptide_score(branching_aa, exp_spectrum), branching_aa])
                
                
        # sort the peptides by decreasing score
        branch_bound.sort()
        branch_bound.reverse()
        print('ldb before cutting', len(leaderboard))
        
        # cut the peptides from the branching step by score rank
        # keep the N peptides, including ties
                        
        # check if there are ties
        
        # keep all if less values than stop otherwise cut any values with score below score at stop
        if len(branch_bound) >= N:
            # find the minimum score to keep
            threshold = branch_bound[N-1][0]
        # find if there are ties
        for i in range(N-1, len(branch_bound)):
            if branch_bound[i][0] == threshold:
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
        print('ldb before removal', len(leaderboard))
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
        
        print('ldb after removal', len(leaderboard))
        leaderboard.reverse()
    
    return leaderpeptide
    
   

