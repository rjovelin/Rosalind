# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 22:28:04 2015

@author: Richard
"""

def profile_most_probable_kmer(input_file):
    '''
    (file) -> str
    Given a string DNA sequence, an integer k, and a 4 × k matrix Profile,
    return the Profile-most probable k-mer in DNA (return any one if multiple
    solutions exists)
    '''    
    # open file
    infile = open(input_file)
    DNA = infile.readline().rstrip()
    k = int(infile.readline().rstrip())
    
    # create a dict to hold profile
    profile = {}
    # populate profile with list of position proba   
    profile['A'] = infile.readline().rstrip().split()
    profile['C'] = infile.readline().rstrip().split()
    profile['G'] = infile.readline().rstrip().split()
    profile['T'] = infile.readline().rstrip().split()
    # close file
    infile.close()
    
    # convert list elements in profile into float
    for nucleotide in profile:
        for i in range(len(profile[nucleotide])):
            profile[nucleotide][i] = float(profile[nucleotide][i])
    
    # make a dict to store each kmer : probability pairs
    kmers = {}
    # go through DNA, grab k-mer and compute probability based on profile
    for i in range(len(DNA) -k+1):
        # grab motif
        motif = DNA[i:i+k]
        # compute probability, initialize proba
        proba = 1
        for j in range(len(motif)):
            proba *= profile[motif[j]][j]
        # store kmer : proba in dict
            kmers[motif] = proba
            
    # create a list to store the (proba, kmer) 
    kmer_probs = [(prob, seq) for seq, prob in kmers.items()]
    # sort the list
    kmer_probs.sort()
    # check the highest probability and update most_probable and highest prob
    highest_proba = kmer_probs[-1][0]
    most_probable = kmer_probs[-1][1]
        
    return most_probable

        

    
    