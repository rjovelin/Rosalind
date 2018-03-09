# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 13:52:07 2016

@author: RJovelin
"""


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



def GreedyMotifSearchPseudoCounts(dna, k, t):
    '''
    (str, int, int) -> list
    Return the most likely kmer motif from the collection of sequences in dna
    Replace 0 probabilities with a pseudocounts in profile matrices
    '''
    # initialize BestMotif by taking the first kmer in each sequence of dna
    BestMotif = [seq[:k] for seq in dna]
    
    # iterate over first seq in dna, grab kmer
    for i in range(len(dna[0]) -k + 1):
        # assign kmer to motif1
        Motif = [dna[0][i:i+k]]
        for j in range(1, t):
            # generate profile from Motif
            profile = ProfilePseudoCounts(Motif)
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