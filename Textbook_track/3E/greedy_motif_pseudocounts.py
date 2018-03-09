# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 10:41:17 2015

@author: Richard
"""

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
    # probabilities with pseudocounts: (s+1) / (n+2)
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
        profile['A'].append(round(A/len(sequences)+2, 4))
        profile['T'].append(round(T/len(sequences)+2, 4))
        profile['C'].append(round(C/len(sequences)+2, 4))    
        profile['G'].append(round(G/len(sequences)+2, 4))    
            
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
        



def greedy_motif_search_pseudocounts(input_file):
    '''
    (file) -> dict
    Given 2 positive integers k and t and a collection of string DNA sequences,
    return the motif (ie a collection of k-mer of length k and size t) that
    are profile most probable k-mer in each DNA sequence
    Uses pseudocounts to generate the profile matrix
    '''
    pass
       
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
    
    # Principle: Construct a t x k motif matrix row-by-row. t is the number of the strings.
    # the ith row is from the ith string. A newly added ith row is selected from the ith string 
    # and the selection is guided by the previous i-1 rows that form a profile to select the ith row. 
    
    # 1. Form an initial BestMotifs by using the 1st k-mer from each of the strings.
    # BestMotif has a t x k dimension.
    # (Taking the 1st k-mer is arbitrary, but the grader is expecting you to do so.)
       
    bestmotif = {}
    for seq in DNA:
        bestmotif[seq] = DNA[seq][:k]
        
    # 2. Assign t x k as the score of BestMotif. Store it in a variable called BestScore.
    # Note:  t x k is the maximum score that a motif can have.
    # e.g. if there are t = 5 strings and k =3, the maximum score is 5 x 3 = 15.    
    bestscore = motif_score(bestmotif)
    
    # repeat 3 to 5
    for i in range(len(DNA[0]) -k+1):
        # grab the 1st kmer of the 1st string sequence
        kmer = DNA[0][i:i+k]
    
        # 3. Form a Motif by using the 1st k-mer of the 1st string. This Motif has a 1 x k dimension.
        motif = {}
        motif[0] = kmer
        
        # Compute a profile by using this 1 x k dimension Motif.
        profile = motif_matrix_pseudocounts(motif)
        
        for j in range(1, len(DNA)):
            # Use this profile to select the most probable k-mer from the 2nd string, the selected k-mer will
            # be the k-mer most similar to the k-mer of the 1st string.
            most_probable = most_probable_kmer(DNA[j], profile, k)
        
            # So, we now have two similar k-mers, one from the 1st string and the other from 2nd string.
            # These two k-mers form a new Motif, which is 2 x k dimension.
            motif[j] = most_probable
                       
            # Compute a profile by using this 2 x k dimension Motif and select the most probable
            # k-mer from the 3rd string, ..., and so on until you select the most probable k-mer from the last string.
            profile = motif_matrix_pseudocounts(motif)
        
        # Now, you get a Motif that has a t x k dimension, which has the same dimension as the BestMotif.
        # 4. Compute the score of the t x k Motif you have just got, say, scoreMotif,
        # if scoreMotif  < BestScore, then BestMotif = Motif, and BestScore = scoreMotif
        
        scoremotif = motif_score(motif)
        if scoremotif < bestscore:
            bestscore = scoremotif
            bestmotif = motif
         
        # 5. Form a NEW 1 x k Motif by using the 2nd k-mer of the 1st string. Repeat step 3 and 4,
        # until you have gone through all the k-mers of the 1st string.

    # 6. Output BestMotif.
    for i in range(len(bestmotif)):
        print(bestmotif[i])
       
    
