# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 17:11:16 2015

@author: Richard
"""

# import Seq from biopython to manipulate sequences
from Bio.Seq import Seq


def which_genetic_code(input_file):
    '''
    (file) -> int
    Given a DNA string s and and a protein string translated by s,
    return the index of the genetic code variant that was used for translation
    (return any index if multiple solutions exist)
    Precondition: S starts with the first codon position, and introns/stop
    codons are ignored
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # grab DNA and protein sequences
    S = infile.readline().rstrip()
    protein = infile.readline().rstrip()
    # close file after reading
    infile.close()
    
    # create a SeqObject
    DNA = Seq(S)
           
    # create a list to store the valid codon tables
    tables = []
    
    # Note: codon tables range from 1 to 15 and tables 7 and 8 have been deleted
    
    # translate DNA with alternative codon tables
    for i in range(1, 7):
        # to_stop = True, stops translation at the first stop codon encountered
        # does not include stop codon in translation
        translation = DNA.translate(table = i, to_stop = True)
        # compare translation with protein
        if translation == protein:
            # store the index of the codon table if protein matches translation
            tables.append(i)
    for i in range(9, 15 + 1):
        # to_stop = True, stops translation at the first stop codon encountered
        # does not include stop codon in translation
        translation = DNA.translate(table = i, to_stop = True)
        # compare translation with protein
        if translation == protein:
             # store the index of the codon table if protein matches translation
            tables.append(i)
    
    # check that S can be translated into protein
    if len(tables) != 0:
        # return first index of codon table recorded
        return tables[0]
    else:
        print('no codon table can translate S in protein')
    
    
    
    
