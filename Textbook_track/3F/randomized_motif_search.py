# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 11:23:03 2015

@author: Richard
"""

import random

def motif_score(sequences):
    '''
    (dict) -> int
    Return the score of a collection of sequences called Motif
    The score is defined as the sum of the number of unpopular nucleotides
    at each position of the sequences
    Precondition: all the sequences in the collection have same length
    '''
    # initialize score
    score = 0
    
    # make a list of keys in dict
    seqnames = [i for i in sequences.keys()]
       
    # go over each nucleotide of each sequence, count the number of unpopular nucleotides and sum to get the score
    for i in range(len(sequences[seqnames[0]])):
        # find the most popular nucleotide
        # make a dict to store the nucleotides and count
        nucleotides = {}
        # populate the dict with the nucleotides in sequences at position i
        for seq in sequences:
            if sequences[seq][i] in nucleotides:
                nucleotides[sequences[seq][i]] += 1
            else:
                nucleotides[sequences[seq][i]] = 1
        # find the most popular nucleotide
        popularity = [(count, base) for base, count in nucleotides.items()]
        popularity.sort() # sort nucleotides by count
        popular = popularity[-1][1]
        # get the position score
        position_score = 0
        for base in popularity[:-1]:
            position_score += base[0]
        # add position score to the sequences score
        score += position_score
    
    return score
        
    
def motif_matrix_pseudocounts(sequences):
    '''
    (dict) -> dict
    Return the profile of sequences as a dictionnary of list with probabilities
    of nucleotides key at each position in the sequences collection.
    Probabilities are adjusted using pseudocounts to avoid zeros 
    '''
    # probabilities with pseudocounts: (s+1) / 2n
    # with s # success and n # sample size
      
    # make a list of keys in dict
    seqnames = [i for i in sequences.keys()]
    
    # create a dict to store the frequencies
    profile = {'A': [], 'T': [], 'C': [], 'G': []}
       
    # go over each nucleotide of each sequence and store frequency
    for i in range(len(sequences[seqnames[0]])):
        # get the nucleotide count and divide by the number of sequences
        # add 1 to all values to generate pseudocounts
        A, C, G, T = 1, 1, 1, 1
        # populate the list of frequencies
        for seq in sequences:
            if sequences[seq][i]  == 'A':
                A += 1
            elif sequences[seq][i] == 'C':
                C += 1
            elif sequences[seq][i] == 'G':
                G += 1
            elif sequences[seq][i] == 'T':
                T += 1
        # divide nucleotide counts by the (number of sequence +2)
#        profile['A'].append(round(A/(2*len(sequences)), 4))
#        profile['T'].append(round(T/(2*len(sequences)), 4))
#        profile['C'].append(round(C/(2*len(sequences)), 4))    
#        profile['G'].append(round(G/(2*len(sequences)), 4))    
        
        profile['A'].append(round(A/(len(sequences)+2), 4))
        profile['T'].append(round(T/(len(sequences)+2), 4))
        profile['C'].append(round(C/(len(sequences)+2), 4))    
        profile['G'].append(round(G/(len(sequences)+2), 4))    
            
    return profile



def most_probable_kmer(sequence, profile, k):
    '''
    (str, dict, int) -> str
    Return the most probable k-mer in a DNA string sequence given the profile
    '''    
    # make a dict to store each kmer : probability pairs
    kmers = {}
    # go through DNA, grab k-mer and compute probability based on profile
    for i in range(len(sequence) -k+1):
        # grab motif
        motif = sequence[i:i+k]
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
    
    # check if there are multiple most probable kmer
    total = 0
    for motif in kmers:
        if kmers[motif] == highest_proba:
            total += 1
    
    # if more than 1, find the first most_probable motif in sequence
    if total > 1:
        # go through DNA, grab k-mer and compute probability based on profile
        for i in range(len(sequence) -k+1):
            # grab motif
            motif = sequence[i:i+k]
            # compute probability, initialize proba
            proba = 1
            for j in range(len(motif)):
                proba *= profile[motif[j]][j]
            if proba == highest_proba:
                return motif
    else:
        return most_probable
        

def randomized_motif_search(input_file):
    '''
    (file) -> dict
    Given 2 positive integers k and t and a collection of string DNA sequences,
    return the motif (ie a collection of k-mer of length k and size t) that
    are profile most probable k-mer in each DNA sequence
    Uses pseudocounts to generate the profile matrix
    '''
    # open file for reading
    infile = open(input_file, 'r')
    # get k and t
    kt = infile.readline().rstrip().split()
    k = int(kt[0])
    t = int(kt[0])
    
    # make a dictionnary to hold the sequences
    DNA = {}
    j = 0
    for line in infile:
        line = line.rstrip()
        if line != '':
            DNA[j] = line
            j += 1
    
    # close after reading
    infile.close()
    
    # Form an initial BestMotifs by randomly selecting 1 kmer in each sequence
    bestmotif = {}
    for seq in DNA:
        i = random.randint(0, len(DNA[seq]) -k)
        kmer = DNA[seq][i:i+k]
        bestmotif[seq] = kmer
        
    # Compute the score of Bestmotif --> Best score
    bestscore = motif_score(bestmotif)
    
    # need to repeat the procedures multipe times because the randmized
    # algorithm can get stuck on a local optima if the search is done once
    # repeat 1000 times
    j = 1000
    while j != 0:
        print(j)
        
        # form a motif by randomly choosing kmers in each sequence
        motif = {}
        for seq in DNA:
            i = random.randint(0, len(DNA[seq]) -k)
            kmer = DNA[seq][i:i+k]
            motif[seq] = kmer
        
        # calculate score of random motif
        score = motif_score(motif)
        
        # generate a profile from random motif
        profile = motif_matrix_pseudocounts(motif)
        
        # set up a boolean
        score_improvement = True
        
        # find a most probable kmer motif, iterate until score doesn't improve
        while score_improvement:
            
            # generate a profile most probable kmer motif
            new_motif = {}
            for m in range(0, len(DNA)):
                # Use profile to select the most probable k-mer from each string
                most_probable = most_probable_kmer(DNA[m], profile, k)
                # make a motif with each most probable kmer
                new_motif[m] = most_probable
            
            # calculate score of profile most probable kmer motif  
            score_new_motif = motif_score(new_motif)
            # generate a new profile, this profile will be used to find the
            # next profile most probable kmer motif 
            profile = motif_matrix_pseudocounts(new_motif)
            
            # check if score is improving
            if score_new_motif < score:
                # score is improving
                score_improvement = True
                # use score of new motif as baseline for comparison
                score = score_new_motif
            else:
                # score is not improving
                score_improvement = False
                # retrieve the new motif and its score, compare score to bestscore
                scoremotif = score_new_motif
                motif = new_motif
          
        # check if score of motif < score bestmotif
        # eventually update the bestmotif, bestscore and profile         
        scoremotif = motif_score(motif)
        if scoremotif < bestscore:
            bestscore = scoremotif
            bestmotif = motif
            
        # update j
        j -= 1

    # output BestMotif
    for i in range(len(bestmotif)):
        print(bestmotif[i])
       
    

    
    
    
    
