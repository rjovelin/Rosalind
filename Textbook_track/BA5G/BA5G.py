# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 12:53:38 2015

@author: RJovelin
"""






# use this function to compute the edit distance between 2 sequences with dynamic programming
def editDistance(seq1, seq2):
    '''
    (str, str) -> int
    Return the edit distance (Levenshtein distance) between sequences 1 and 2
    using dynamic programming
    '''

    # initiate matrix
    D = []
    # populate matrix with lists of 0, including the empty string before seq1 and seq2
    for i in range(len(seq1) + 1):
        D.append([0] * (len(seq2)+1))
    # initiate first row and first column 
    # with distance from empty string to seq1 and seq2
    for i in range(1, len(seq1) + 1):
        D[i][0] = i
    for i in range(1, len(seq2) + 1):
        D[0][i] = i
    # fill matrix with distances
    # loop over the first sequence
    for i in range(1, len(seq1) + 1):
        # start at index 1 to keep distance = 0 between empty strings
        # loop over second sequence
        for j in range(1, len(seq2) + 1):
            # compute distance at position x, y in matrix based on the 
            # minimum distance from (i-1, j), (i-1, j-1), (i, j-1)
            distHor = D[i][j-1] + 1
            distVert = D[i-1][j] + 1
            # check if nucleotides at -1 position are the same in seq1 and seq2
            if seq1[i-1] == seq2[j-1]:
                # no penalty
                distDiag = D[i-1][j-1]
            else:
                # add penalty
                distDiag = D[i-1][j-1] + 1
            # take the minimum distance to fill the matrix
            D[i][j] = min(distHor, distVert, distDiag)
    
    # EDit distance is the distance in the last row and last column of the matrix
    return D[-1][-1]



def ba5g(rosalind_file):
    infile = open(rosalind_file, 'r')
    seq1 = infile.readline().rstrip()
    seq2 = infile.readline().rstrip()
    infile.close()
    print(editDistance(seq1, seq2))