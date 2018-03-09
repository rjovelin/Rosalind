# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 14:28:54 2015

@author: RJovelin
"""




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
    

def ba4a(filename):
    infile = open(filename, 'r')
    rna = infile.readline().rstrip()
    protein = translate_rna(rna)
    print(protein)