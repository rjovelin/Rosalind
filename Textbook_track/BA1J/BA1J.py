# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 21:25:59 2015

@author: Richard
"""


def dna_rev_compl(dna):
    '''
    (str) -> str
    Return the reverse complement of a DNA sequence
    Precondition: dna has only valid nucleotides
    
    >>> dna_rev_compl('AAAACCCGGT')
    ACCGGGTTTT
    '''
    
    # create a dictionary of complementary bases
    complement_nt = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # create a string collector
    revcomp = ''
    # loop over dna:
    for i in range(len(dna)):
        # gett he complement of nucleotide in dna and add to the begining of the 
        # new string to get the reverse sequence
        revcomp = complement_nt[dna[i].upper()] + revcomp
    
    return revcomp
    
    
import itertools

def hamming_distance(seq1, seq2):
    '''
    (str, str) -> int
    Return the hamming distance (number of differences) between seq1 and seq2
    Precondition: seq1 and seq2 have same length
    
    >>> hamming_distance('GGGCCGTTGGT', 'GGACCGTTGAC')
    3
    '''
    D = 0
    
    # loop over seq1
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            D += 1
    return D


def build_k_mer_index(text, kmer):
    '''
    (str, int) -> dict
    Take a string text, the length of a kmer and return a dictionary with 
    strings of length kmer in text as key and a list of the kmers positions in text
    '''
    
    # initiate dict
    text_index = {}
    
    # loop ober text
    for i in range(0, len(text) - kmer + 1):
        # grab kmer string
        subseq = text[i:i+ kmer]
        # check if subseq in dict
        if subseq not in text_index:
            # add key and intiate list value
            text_index[subseq] = [i]
        elif subseq in text_index:
            # add index to list
            text_index[subseq].append(i)
            
    return text_index
    

def approximate_matches(p, T_index, T, differences, kmer):
    '''
    (str, dict, str, int, int) -> list, int
    Take pattern p, a text T and the maximum number of differences between p and 
    its occurence in T and return a tuple with a list with indices where p matches T with at
    most the number of differences and the number of kmers in p mattching an kmer index in index of T
    '''
    
    # make a set of indices matching p
    # to avoid counting matches to the same position
    # (in case different kmers are perfect matches in index)
    matches = set()
               
    # build a list of query strings
    queries = [p[i:i+kmer] for i in range(0, len(p),kmer)] 
        
    # check if the queries are matches in index
    for i in range(len(queries)):
        # check if substring is kmer in index of T
        if queries[i] in T_index:
            # for all matches of query string in T:
            for pos in T_index[queries[i]]:
                # count the number of mismatches for the
                # sequence upstream and downstream of query kmer
                upstream = ''.join(queries[0:i])
                downstream = ''.join(queries[i+1:])
                # get the sequence in T aligned to upstream sequence
                T_up = T[pos - len(upstream): pos]
                T_down = T[pos + kmer:pos + kmer + len(downstream)]
                # get the start position of upstream , ie = start position of the pattern match
                start = pos - len(upstream)
                # count the number of mismatches
                if hamming_distance(T_up, upstream) + hamming_distance(T_down, downstream) <= differences:
                    # number of mismatches is tolerable
                    # add match position to list
                    matches.add(start)
                    
    return matches





def most_frequent_mismatched_kmer_revseq(genome, k, d):
    '''
    (str, int, int) -> list
    Return a list with the most frequent kmers in genome with at most d mismatches
    
    >>> most_frequent_mismatched_kmer_revseq('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1)
    ['ATGT', 'ACAT']
    '''
    
    # create to store kmers : count
    kmers = set()
    alphabet = 'ACTG'
    for i in itertools.product(alphabet, repeat = k):
        kmers.add(''.join(list(i)))
    # create dict to store kmer : counts
    counts = {}
    print(len(kmers))
    
    genome_index = build_k_mer_index(genome, k)
    
    # loop over kmers
    for seq in kmers:
        # get the positions of seq in genome
        matches = approximate_matches(seq, genome_index, genome, d, k)
        # get rev seq
        revseq = dna_rev_compl(seq)
        revmatches = approximate_matches(revseq, genome_index, genome, d, k)
        # check if kmer in dict
        if seq in counts:
            counts[seq] += (len(matches) + len(revmatches))
        else: counts[seq] = len(matches) + len(revmatches)
                   
                       
    # make a reverse dict : count : list of kmers
    counts_to_kmers = {}
    for seq in counts:
        if counts[seq] in counts_to_kmers:
            counts_to_kmers[counts[seq]].append(seq)
        else:
            counts_to_kmers[counts[seq]] = [seq]
            
    # make a list of counts
    L = [i for i in counts_to_kmers]
    highest = max(L)
    return counts_to_kmers[highest]
    
    
def ba1j(filename):
    infile = open(filename, 'r')
    genome = infile.readline().rstrip()
    k , d = infile.readline().rstrip().split()
    k = int(k)
    d = int(d)
    infile.close()
    most_freq = most_frequent_mismatched_kmer_revseq(genome, k, d)
    for i in most_freq:
        print(i, end = ' ')