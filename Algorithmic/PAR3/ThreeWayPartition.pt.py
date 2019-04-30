# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 22:58:45 2019

@author: Richard
"""

def ThreeWayPartition(A):
    '''
    (list) -> list
    Take a list A of integers from −10**5 to 10**5 and return an array B[1..n]
    such that it is a permutation of A and there are indices 1≤q≤r≤n such that B[i]<A[1] for all 1≤i≤q−1, B[i]=A[1] for all q≤i≤r, and B[i]>A[1] for all r+1≤i≤n.
    '''
    
    B = A[::]
    
    K = B[0]
    # keep track of position of K as it moves each time elements are insterted
    # at position 0
    j = 0
    
    for i in range(1, len(B)):
        if B[i] < K:
            # take element and insert at the begining of A
            B.insert(0, B.pop(i))
            # adjust new position of K in A
            j += 1
        # insert element where K is curently located
        elif B[i] ==K:
            B.insert(j, B.pop(i))
        
    return B



def SolvePb(PbFile):
    infile = open(PbFile)
    infile.readline()
    L = list(map(lambda x: int(x), infile.readline().strip().split()))
    A = list(map(lambda x: str(x), ThreeWayPartition(L)))
    infile.close()
    newfile = open('threeway_solution.txt', 'w')
    newfile.write(' '.join(A))
    newfile.close()
    
    