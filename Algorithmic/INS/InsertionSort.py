# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 00:53:09 2019

@author: Richard
"""




def InsertionSort(A):
    '''
    
    
    '''
    
    N = 0
    
    for i in range(1, len(A)):
        k = i
        while k >= 1 and A[k] < A[k-1]:
            # swap A[k] and A[k-1]
            A[k-1], A[k] = A[k], A[k-1]
            k -= 1
            # record number of swap
            N +=1

    return N




