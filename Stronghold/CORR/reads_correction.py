# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 21:24:17 2015

@author: Richard
"""

def rev_complement(dna):
    '''
    (str) -> str
    Return the reverse complement of string dna
    Precondition: there are only valid nucleotides in dna
    >>> rev_complement('atcg')
    'cgat'
    >>> rev_complement('TTCGAT')
    'ATCGAA'
    '''

    dna2 = dna.upper()
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    reverse_comp_L = []
    for i in reversed(dna2):
        reverse_comp_L.append(complement[i])
    reverse_comp_dna = ''.join(reverse_comp_L)

    return reverse_comp_dna


def Hamming_D(s1, s2):
    '''
    (str, str) -> int
    Return the Hamming distance (# of differences) between 2 sequences
    s1 and s2. Precondition: s1 and s2 have the same length
    >>> Hamming_D('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT')
    7
    '''
    diffs = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            diffs += 1
    return diffs

def reads_correction(input_file,output_file):
    '''
    (file) -> file
    
    Given: A collection of up to 1000 reads of equal length (at most 50 bp)
    in FASTA format. Some of these reads were generated with a
    single-nucleotide error. For each read s in the dataset,
    one of the following applies:
    - s was correctly sequenced and appears in the dataset at least twice
    (possibly as a reverse complement);
    - s is incorrect, it appears in the dataset exactly once,
    and its Hamming distance is 1 with respect to exactly one correct read
    in the dataset (or its reverse complement).
    Return: A list of all corrections in the form "[old read]->[new read]".
    (Each correction must be a single symbol substitution,
    and you may return the corrections in any order.)
    '''
    
    # open file
    infile = open(input_file, 'r')
    # make a dict with the sequences
    fasta = {}
    for line in infile:
        line = line.rstrip()
        if line != '':
            if line.startswith('>'):
                name = line[1:]
                fasta[name] = ''
            else:
                fasta[name] += line.upper()
    
    # close file
    infile.close()
    
    # make a list with all sequence reads
    reads = [fasta[gene] for gene in fasta]
    
    # make a set of incorrect read IDs
    incorrect = set()
    for seq in fasta:
        if reads.count(fasta[seq]) == 1 and\
        reads.count(rev_complement(fasta[seq])) == 0:
            incorrect.add(fasta[seq])
    
    
    # make a set of correct reads and their reverse complement
    correct = set()
    for read in reads:
        if read not in incorrect:
            correct.add(read)
            correct.add(rev_complement(read))
    
    # make a dict of incorrect read : 1-hamming-read
    corrected = {}
    for read in incorrect:
        for correct_read in correct:
            if Hamming_D(read, correct_read) == 1:
                if read not in corrected:
                    corrected[read] = [correct_read]
                else:
                    corrected[read].append(read)
                    
    
    
    # open file for writing
    newfile = open(output_file, 'w')
        
    for read in corrected:
        if len(corrected[read]) == 1:
            newfile.write(read + '->' + corrected[read][0] + '\n')
        else:
            newfile.write(read + '->')
            for item in corrected[read][:-1]:
                newfile.write(item + ', ')
            newfile.write(corrected[read][-1] + '\n')
    newfile.close()
            
    
    
    
    
    