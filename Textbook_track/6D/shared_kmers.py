# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 22:41:52 2015

@author: Richard
"""


def kmer_positions(S, k):
    '''
    (str, int) -> dict
    Return a dictionnary with the pairs of kmer: list of indices of
    kmer in S   
    '''
    
    # make a dictionnary to store indices of kmers in each sequence
    S_pos = {}
    for i in range(len(S)-k+1):
        kmer = S[i:i+k]
        if kmer in S_pos:
            S_pos[kmer].append(i)
        else:
            S_pos[kmer] = [i]
    
    return S_pos


def reverse_complement(S):
    '''
    (str) -> str
    Return the reverse complement of sequence S
    '''
    S = S.upper()
    complements = {'A': 'T', 'T': 'A', 'C':'G', 'G':'C'}
    S_rc = ''
    for i in reversed(S):
        S_rc += complements[i]
    
    return S_rc
    

# find the shared k-mers between 2 sequences
def shared_kmers(input_file):
    '''
    (file) -> list of tuples
    Given an integer k and two string sequences in input_file, 
    return all k-mers shared by these strings, in the form of ordered
    pairs (x, y). The ordered pairs is made of the indices of occurence of the
    k-mers in string 1 and string 2    
    '''
        
    # open file for reading
    infile = open(input_file, 'r')
    # get k
    k = int(infile.readline().rstrip())
    # get sequences
    S1 = infile.readline().rstrip()
    S2 = infile.readline().rstrip()
    # close file
    infile.close()
    
    # make a dictionnary to store indices of kmers in each sequence
    S1_pos = kmer_positions(S1, k)
    S2_pos = kmer_positions(S2, k)
    
    # make dictionnaries to store the indices of the kmers in the reverse compl
    S1_rc_pos = kmer_positions(reverse_complement(S1), k)
    S2_rc_pos = kmer_positions(reverse_complement(S2), k)
    
    # create a list to store the index pairs
    ordered_pairs = []
    for kmer in S1_pos:
        if kmer in S2_pos:
            for position in S1_pos[kmer]:
                for location in S2_pos[kmer]:
                    ordered_pairs.append((position, location))
        if reverse_complement(kmer) in S2_pos:
            for position in S1_pos[kmer]:
                for location in S2_pos[reverse_complement(kmer)]:
                    ordered_pairs.append((position, location))
                
    # sort ordered pair according to the position of the kmer in sequence 1
    ordered_pairs.sort()                
                           
    for pair in ordered_pairs:
        print(pair)
            
        
    
    