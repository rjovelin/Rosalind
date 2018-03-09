# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 02:00:02 2015

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
        

def gibbs_sampler(input_file):
    '''
    (file) -> dict
    Given 3 positive integers k, t and N and a collection of string DNA
    sequences, return the motif (ie a collection of k-mer of length k and
    size t) of sequences. Uses pseudocounts to generate the profile matrix
    '''
    
    # In contrast to randomized_search in which motifs can be completely
    # changed at each iteration and bestmotif is only evaluated when motif_score
    # doesn't improve, in Gibbs sampling the motif is changed 1 kmer
    # at a time and the bestmotif is updated at each iteration of motif
    # This requires iterations to be sufficiently large to explore the sample
    # space, and also requires to start many times from random points to avoid
    # local optima
    
    # open file for reading
    infile = open(input_file, 'r')
    # get k and t
    kt = infile.readline().rstrip().split()
    k = int(kt[0])
    t = int(kt[1])
    N = int(kt[2])
    
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
    j = 200
    while j != 0:
                
        # form a motif by randomly choosing kmers in each sequence
        motif = {}
        for seq in DNA:
            i = random.randint(0, len(DNA[seq]) -k)
            kmer = DNA[seq][i:i+k]
            motif[seq] = kmer
        
        # calculate score of random motif
        score = motif_score(motif)
        
        N = 200
       
        # find a most probable kmer motif, update motif
        while N != 0:
            # randomly choose a sequence
            p = random.randint(0, len(DNA) -1)
            sequence_star = DNA[p]
            
            # remove the kmer from motif selected from that sequence 
            del motif[p]
            
            # generate a profile from random motif
            profile = motif_matrix_pseudocounts(motif)
            
            # find most probable kmer in sequence_star
            new_kmer = most_probable_kmer(sequence_star, profile, k)
                       
            # add the profile most probable kmer to motif
            motif[p] = new_kmer
                       
            # calculate score of profile most probable kmer motif  
            score_motif = motif_score(motif)
            
            # compare motif score to bestscore
            # eventually update bestscore and bestmotif
            if score_motif < bestscore:
                bestscore = score_motif
                bestmotif = motif
            

            N -= 1  
             
        # update j
        j -= 1

    # output BestMotif
    for i in range(len(bestmotif)):
        print(bestmotif[i])
       
    

    
    
    
    
