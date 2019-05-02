# -*- coding: utf-8 -*-
"""
Created on Wed May  1 22:52:51 2019

@author: Richard
"""


def ThreeWayPartition(A):
    '''
    (list) -> list
    Take a list A of integers from −10**5 to 10**5, modify A such that it is
    a permutation of A and there are indices 1≤q≤r≤n such that A[i]<A[1]
    for all 1≤i≤q−1, A[i]=A[1] for all q≤i≤r, and A[i]>A[1] for all r+1≤i≤n
    and return al 3 correspondng arrays split on q and r
    '''
    
    
    # pick first value as pivot
    K = A[0]
    # keep track of position of K as it moves each time elements are insterted
    # at position 0
    j = 0
    
    for i in range(1, len(A)):
        if A[i] < K:
            # take element and insert at the begining of A
            A.insert(0, A.pop(i))
            # adjust new position of K in A
            j += 1
        # insert element where K is curently located
        elif A[i] ==K:
            A.insert(j, A.pop(i))
            
    m = j + A.count(K)
    return A[:j] + A[j:m] + A[m:], j, m
    

def QuickSort(A):
   
    # use divide and conquer to sort
    # partition A in 3, and recursively sort the 1st and 3rd sublists
    # return list if empty
    if len(A) > 1:
        A, j, m = ThreeWayPartition(A)
        return QuickSort(A[:j]) + A[j:m] + QuickSort(A[m:])
    else:
        return A
    
    
def SolvePb(PbFile):
    infile = open(PbFile)
    infile.readline()
    L = list(map(lambda x: int(x), infile.readline().strip().split()))
    infile.close()
    L = QuickSort(L)
    L = list(map(lambda x: str(x), L))
    newfile = open('quicksort_solution.txt', 'w')
    newfile.write(' '.join(L))
    newfile.close()