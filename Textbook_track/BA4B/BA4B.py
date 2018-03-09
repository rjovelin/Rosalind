# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 14:59:01 2015

@author: RJovelin
"""



def rev_complement(dna):
    '''
    (str) -> str
    Return the reverse complement of dna
    Precondition: dna has only valid nucleotides
    '''

    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_dna = ''
    for i in dna:
        rev_dna = complement[i] + rev_dna
        
    return rev_dna



def translate_rna(rna):
    '''
    (str) -> str
    
    Return the amino acid translation of RNA, without the terminak stop codon
    
    >>> translate_rna('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA')
    'MAMAPRTEINSTRING'
    '''
    
    genetic_code = {'GCU': 'A' , 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	              'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L',
                    'CUA': 'L', 'CUG': 'L', 'CGU': 'R', 'CGC': 'R',
                    'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
                    'AAA': 'K', 'AAG': 'K', 'AAU': 'N', 'AAC': 'N',
                    'AUG': 'M', 'GAU': 'D', 'GAC': 'D', 'UUU': 'F',
                    'UUC': 'F', 'UGU': 'C', 'UGC': 'C', 'CCU': 'P',
                    'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAA': 'Q',
                    'CAG': 'Q', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S',
                    'UCG': 'S', 'AGU': 'S', 'AGC': 'S', 'GAA': 'E',
                    'GAG': 'E', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T',
                    'ACG': 'T', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
                    'GGG': 'G', 'UGG': 'W', 'CAU': 'H', 'CAC': 'H',
                    'UAU': 'Y', 'UAC': 'Y',  'AUU': 'I', 'AUC': 'I',
                    'AUA': 'I', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
                    'GUG': 'V', 'UAA': '*', 'UGA': '*', 'UAG': '*'}
    
    # create protein string
    protein= ''    
    
    # loop over rna, grab codons, translate and update protein
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        protein += genetic_code[codon]
    if protein[-1] == '*':
        protein = protein[:-1]
    return protein


def find_peptide_dna(dna, peptide):
    '''
    (str, str) -> list
    Return a list of substrings of dna that encode or for which the reverse
    complement encodes peptide
    
    >>> find_peptide_dna('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA', 'MA')
    ['ATGGCC', 'GGCCAT', 'ATGGCC']
    '''
    
    # create list to store substrings
    seqs = []
    # get length of substring
    L = len(peptide) * 3
    
    # loop over dna for the given frame
    for i in range(len(dna)):
        subseq = dna[i: i+L]
        if len(subseq) % 3 == 0 and len(subseq) > 0:
            # grab sequence, convert to RNA
            subseq = subseq.replace('T', 'U')
            # take the reverse complement
            rvsubseq = rev_complement(dna[i: i+L]).replace('T', 'U')
            # check if translation is same as peptide
            if translate_rna(subseq) == peptide:
                seqs.append(subseq.replace('U', 'T'))
            elif translate_rna(rvsubseq) == peptide:
                seqs.append(subseq.replace('U', 'T'))
    return seqs
            

def ba4b(filename):
    infile = open(filename, 'r')
    dna = infile.readline().rstrip()
    peptide = infile.readline().rstrip()
    seqs = find_peptide_dna(dna, peptide)
    for i in seqs:
        print(i)