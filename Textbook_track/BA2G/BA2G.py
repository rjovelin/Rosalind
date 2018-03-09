# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 10:57:44 2016

@author: Richard
"""



import math
import random


# use this function to compute the score of a motif in a collection of sequences
def ScoreMotif(motif):
    '''
    (list) -> int
    Return the score of motif. The score is the sum of the number of nucleotides
    differing from the most abundant nucleotide at each position of the kmer
    in the collection of kmer motif. A lower score indicates stronger conservation
    Precondition: no gaps in the motif sequences
    '''
    
    # set up variable
    score = 0
    # loop over first kmer in motif 
    for i in range(len(motif[0])):
        # grab all nucleotides at this position in motif
        nucleotides = [kmer[i] for kmer in motif]
        # get the most abundant nucleotides
        # get any nucleotide if more than 1 most frequent
        mostfrequent = {}
        # get nucleotide counts
        for base in nucleotides:
            mostfrequent[base] = mostfrequent.get(base, 0) + 1
        counts = [(val, key) for key, val in mostfrequent.items()]
        counts.sort()
        mostfrequent = counts[-1][-1]
        # count the number of nucleotides different from mostfrequent
        for base in nucleotides:
            if base != mostfrequent:
                score += 1
    return score
     


def PatternProbability(pattern, profile):
    '''
    (str, dict) -> float
    Return the probability of a given kmer pattern given the profile matrix
    of a motf collection of kmers
    Precondition: the kmer pattern and the profile for each nucleotide
    in the matrix have same length
    '''
    # set up probability    
    P = 1
    # loop over pattern
    for i in range(len(pattern)):
        # grab probability of nucleotide at that position
        p = profile[pattern[i]][i]
        # update probability f pattern
        P *= p
    return P




def ProfilePseudoCounts(motif):
    '''
    (list) -> dict
    Return a dictionary with nucleotides as key and list of nucleotide frequencies
    for each position of the kmer motif. Adjust frequencies to convert
    null probabilities to low probabilities
    '''
    
    # create dict
    profile = {}
    # initialize dict
    for base in 'ACGT':
        profile[base] = []
    
    # loop over first kmer in motif 
    for i in range(len(motif[0])):
        # grab all nucleotides at this position in motif
        nucleotides = [kmer[i].upper() for kmer in motif]
        # set up variable to compute frequencies
        total = len(nucleotides) + sum([nucleotides.count(base) for base in 'ACGT']) 
        # get the frequency of each nucleotide
        for base in 'ACGT':
            # add 1 to each count to adjust 0 values
            freq = (nucleotides.count(base) + 1) / total
            profile[base].append(freq)
    return profile    


  
def GibbsSampler(dna, k, t, N):
    '''
    (list, int, int, int) -> list
    
    Return the most likely kmer motif from the collection of sequences in dna
    using a Gibbs Sampler running repeat times
    '''    
    
    # generate motif with randomly selected kmers, one in each seq in dna
    BestMotif = []
    for i in range(len(dna)):
        # select nucleotide index at random
        j = random.randint(0, len(dna[i]) -k)
        # grab kmer
        kmer = dna[i][j:j+k]
        BestMotif.append(kmer)
    
    # run the procedure with 20 random starts
    starts = 20
    for i in range(starts):
        # select a random motif
        motif = []
        for j in range(len(dna)):
            # select nucleotide index at random
            m = random.randint(0, len(dna[j]) -k)
            # grab kmer
            kmer = dna[j][m:m+k]
            motif.append(kmer)
        
        
        # iterate N times
        for j in range(N):
            # select a random sequence in dna
            m = random.randint(0, len(dna) -1)
            # grab the corresponding sequence
            RemovedSeq = dna[m]
            # construct profile with all kmers in Motif but the kmer from removed seq
            motif.remove(motif[m])
            # construct the profile using pseudocounts on the motif - kmerj
            profile = ProfilePseudoCounts(motif)
            # compute the probability of each kmer in removed seq
            patterns = [RemovedSeq[o:o+k] for o in range(len(RemovedSeq) -k + 1)]
            proba = [PatternProbability(pattern, profile) for pattern in patterns]
            # adjust the probabilities
            SumProba = sum(proba)
            for o in range(len(proba)):
                proba[o] = proba[o] / SumProba
            # get the kmer with highest probability
            # insert kmer at position j
            motif.insert(m, patterns[proba.index(max(proba))])
            # compare score alternative motif and best motif
            if ScoreMotif(motif) < ScoreMotif(BestMotif):
                BestMotif = motif
    return BestMotif
    
    
    