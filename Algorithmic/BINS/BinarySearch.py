# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 17:46:58 2019

@author: rjovelin
"""

# The problem is to find a given set of keys in a given array.


import math

# implement binary search
def BinarySearch(A, n, T):
    '''
    (list, int, str) -> int
    Take a sorted list of size n and a return the index of T in A if present or -1
    Precondition: A is sorted
    '''
    
    # set L = 0 and R = n -1
    L, R = 0, n-1
    while L <= R:
        # search T in the array delimited by L and R
        # set m (the position of the middle element) to (L + R) / 2
        m = math.floor((L + R) / 2)
        # check if target T is on the right or the left of m
        if A[m] < T:
            # T on the right of m
            # set L and R to create a new array on the right
            L = m + 1
        elif A[m] > T:
            # T on the left of m, set L and R to create a new array on the left
            R = m - 1
        elif A[m] == T:
            # target found
            return m
    # if L > R the entire array has been searched and target is not found
    return -1



def PrintArray(PbFile):
    '''
    Given: Two positive integers n≤105 and m≤105, a sorted array A[1..n]
    of integers from −10**5 to 10**5 and a list of m integers −105≤k1,k2,…,km≤10**5.
    Return: For each ki, output an index 1≤j≤n s.t. A[j]=ki or "-1" if there is no such index.    
    '''

    infile = open(PbFile) 
    n = int(infile.readline().rstrip())
    m = int(infile.readline().rstrip())
    # make sure to convert str  --> int to get the sorted list
    # assumption of binary seqarch : arrat A is sorted
    
    A = list(map(lambda x: int(x), infile.readline().rstrip().split()))
    B = list(map(lambda x: int(x), infile.readline().rstrip().split()))

    assert n == len(A) 
    assert m == len(B)
    
    L = []
    for i in range(len(B)):
        m = BinarySearch(A, n, B[i])
        if m != -1:
            # record indices 1-base
            L.append(str(m + 1))
        else:
            L.append(str(m))
    print(' '.join(L))
        
