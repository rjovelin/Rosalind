# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 15:54:46 2015

@author: Richard
"""

def hydrophobicity(input_file):
    '''
    (file) -> list
    Given a protein string s of length n, an integer l and a float k,
    return all substrings of s of length l having a hydropathy index
    higher than k.
    Note: the hydropathy index of a protein is the average of the
    hydropathy indices of each amino acid
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    protein = infile.readline().rstrip()
    nums = infile.readline().rstrip().split()
    L = int(nums[0])
    k = float(nums[1])
    # close file after reading
    infile.close()
            
    Kyte_Doolittle = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
                      'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
                      'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
                      'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
    
                      
    # create list to store the kmers with hydropath > k
    kmers = []
    
    # go through protein sequence, grab the kmer of length L, compute H-index
    # evaluate H-index and eventually store kmer in list
    for i in range(len(protein) -L+1):
        # grab kmer
        motif = protein[i:i+L]
        # compute H-index
        hydro = 0
        for aa in motif:
            hydro += Kyte_Doolittle[aa]
        # compare H-index to k
        if hydro / L > k:
            kmers.append(motif)
    # print motif in kmers
    for motif in kmers:
        print(motif)