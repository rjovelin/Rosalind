# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 21:52:03 2015

@author: Richard
"""



# use this function to find the highest global alignment score or the distance matrix
def gap_alignment(seq1, seq2, dist_matrix = 'matrix'):
    '''
    (str, str, str) -> int
    Return the edit distance (Levenshtein distance) between sequences 1 and 2
    if dist_matrix is distance and return the matrix by default.
    Uses a matrix of penalty scores and dynamic programming
    '''
    
    # create a list of nucleotides
    alphabet = ['A', 'C', 'T', 'G']
    # create a matrix of penalty scores between each nucleotides in alphabet,
    # and with empty string, penalize substitutions over gaps
    penalty =  [[0, 10, 10, 10, 5],
                [10, 0, 10, 10, 5],
                [10, 10, 0, 10, 5],
                [10, 10, 10, 0, 5],
                [5, 5, 5, 5, 5]]
    # to get the penalty score between A in seq1 -> C in seq2, use the index of nucleotides in alphabet
    # i = alphabet.index('A'); j = alphabet.index('C'); S = penalties[i][j]
    
    # to get the penalty score between C in seq1 -> A in seq2, use the index of nucleotides in alphabet
    # i = alphabet.index('C'); j = alphabet.index('A'); S = penalties[i][j]
    
    # create a distance matrix
    D = []
    
    # initiate matrix with 0s, with empty string before each sequence
    for i in range(len(seq1) + 1):
        # index i of list in D corresponds to index i in seq1
        D.append([0] * (len(seq2) + 1))
        # index j in each list of D corresponds to index j in seq2
    
    # initiate first column with distance from empty string to seq1, keep distance from empty strings 0 
    for i in range(1, len(seq1) +1):
        # grab the index in alphabet of the nucleotides in seq1
        j = alphabet.index(seq1[i-1])
        # get the penalty score between empty string and nucleotide in seq1
        S = penalty[j][-1]
        # update matrix with penalty score
        # penalty score of seq[i] = penalty score of seq[i-1] + S
        D[i][0]  = D[i-1][0] + S
    
    # initiate first column with distance from empty string to seq2, keep distance from empty strings 0
    for i in range(1, len(seq2) + 1):
        # grab the index in alphabet of the nucleotides in seq2
        j = alphabet.index(seq2[i-1])
        # get the penalty score between empty string and nucleotide in seq2
        S = penalty[-1][j]
        # update matrix with penalty score
        # penalty score of seq[i] = penalty score of seq[i-1] + S
        D[0][i] = D[0][i-1] + S

    # fill the matrix, keep distance between empty strings 0
    # loop over 1st sequence
    for i in range(1, len(seq1) + 1):
        # start at index 1 to keep distance = 0 between empty strings
        # loop over second sequence
        for j in range(1, len(seq2) + 1):
            # compute distance at position x, y in matrix based on the
            # minimum distance from (x-1, y), (x-1, y-1), (x, y-1)
            distHor = D[i][j-1] + penalty[-1][alphabet.index(seq2[j-1])]
            distVert = D[i-1][j] + penalty[alphabet.index(seq1[i-1])][-1]
            distDiag = D[i-1][j-1] + penalty[alphabet.index(seq1[i-1])][alphabet.index(seq2[j-1])]        
            # update distance at coordinate i, j with minumum distance
            D[i][j] = min(distHor, distVert, distDiag)
            
    # global alignment score is the distance in the column of tghe last row
    if dist_matrix == 'distance':
        return D[-1][-1]
    elif dist_matrix == 'matrix':
        return D
    
       
def TraceBack_LCSQ(seq1, seq2):
    '''
    (str, str, list) -> str, str
    Return the aligned sequences of seq1 and seq2 
    '''
    
    # generate the matrix of distances
    D = gap_alignment(seq1, seq2, dist_matrix = 'matrix')    
        
    # get last indices of matrix
    i = len(D) - 1
    j = len(D[-1]) - 1
    
    # create aligned sequences to be updated
    seq1_ali, seq2_ali = '', ''
    
    # start at the last indices in matrix and trace back,
    # update the indices until they reach 0 each and update aligned sequences
    while i != 0 and j != 0:
        print(i, j, seq1_ali, seq2_ali, seq1[i-1], seq2[j-1])
        # get current score
        current_score = D[i][j]
        print(i, j, current_score, current_score -  D[i][j-1], current_score - D[i-1][j-1], current_score - D[i-1][j])
        # compare nucleotides at positions i in seq1 and j in seq2
        if seq1[i-1] == seq2[j-1]:
            # the previous position in the matrix is in the upper left corner
            # update sequences
            seq1_ali = seq1[i-1] + seq1_ali
            seq2_ali = seq2[j-1] + seq2_ali
            # i and j
            i -= 1
            j -= 1
        else:
            # previous position is left or above
            # unless scores above and left > score at current position: previous position is upper left
            if D[i][j-1] <= D[i][j] and D[i][j-1] <= D[i-1][j]:
                # previous position is left
                # update sequences with gap in seq1 and nucleotide in seq2
                seq1_ali = '-' + seq1_ali
                seq2_ali = seq2[j-1] + seq2_ali
                # update j, i stays the same
                j -= 1
            elif D[i][j-1] <= D[i][j] and D[i][j-1] > D[i-1][j]:
                # previous position is above
                # update sequences with nucleotide in seq1 and gap in seq2
                seq1_ali = seq1[i-1] + seq1_ali
                seq2_ali = '-' + seq2_ali
                # update i, j stays the same
                i -= 1
            elif D[i][j-1] > D[i][j]:
                # previous position is above
                # update sequences with nucleotide in seq1 and gap in seq2
                seq1_ali = seq1[i-1] + seq1_ali
                seq2_ali = '-' + seq2_ali
                # update i, j stays the same
                i -= 1
            else:
                # previous position is upper left, substitution less costly than gap
                # update sequences
                seq1_ali = seq1[i-1] + seq1_ali
                seq2_ali = seq2[j-1] + seq2_ali
                # i and j
                i -= 1
                j -= 1
            
    return seq1_ali, seq2_ali    