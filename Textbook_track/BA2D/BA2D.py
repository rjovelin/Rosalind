# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 11:51:58 2016

@author: RJovelin
"""

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


def FindMostProbableKmer(text, k, profile):
    '''
    (str, int, dict) -> str
    Find the most probable kmer pattern in text given the probability matrix profile
    '''
    
    Patterns= {}
    # loop over text, grab each kmer
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        # compute the probability of kmer given profile
        proba = PatternProbability(kmer, profile)
        # update Patterns
        Patterns[kmer] = proba
    # find the most probable kmer
    probable = [(p, kmer) for kmer, p in Patterns.items()]
    probable.sort()
    MostLikely = probable[-1][1]
    highest_proba = probable[-1][0]
    
    # if there are multiple probable kmers,
    # then return the kmer appearing first in text
    total = 0
    for kmer in Patterns:
        if Patterns[kmer] == highest_proba:
            total += 1
    if total > 1:
        # loop over seq, determine proba of each kmer and return the kmer with highest proba
        for i in range(len(text) -k + 1):
            kmer = text[i:i+k]
            if PatternProbability(kmer, profile) == highest_proba:
                return kmer
    else:
        return MostLikely
        

def GreedyMotifSearch(dna, k, t):
    '''
    (str, int, int) -> list
    Return the most likely kmer motif from the collection of sequences in dna
    '''
    
    # initialize BestMotif by taking the first kmer in each sequence of dna
    BestMotif = [seq[:k] for seq in dna]
    
    # iterate over first seq in dna, grab kmer
    for i in range(len(dna[0]) -k + 1):
        # assign kmer to motif1
        Motif = [dna[0][i:i+k]]
        for j in range(1, t):
            # generate profile from Motif
            profile = Profile(Motif)
            # find the most probable kmer in the next sequence
            MostLikely = FindMostProbableKmer(dna[j], k, profile)
            # add the most probable kmer to motif
            Motif.append(MostLikely)
        assert len(Motif) == len(dna)
        # evaulate the score of the current motif
        if ScoreMotif(Motif) < ScoreMotif(BestMotif):
            # update BestMotif is score is better (lower)
            BestMotif = Motif
    return BestMotif