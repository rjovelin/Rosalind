# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 09:34:25 2019

@author: Richard
"""


def TwoWayPartition(A):
    
    # take first element j of A
    # create a second array in which all elements of B up to j are < j and all elements after j are > j
    # add all elements 
    K = A[0]
    
    for i in range(1, len(A)):
        if A[i] <= K:
            A.insert(0, A.pop(i))
    return A
    




def SolvePb(PbFile):
    infile = open(PbFile)
    infile.readline()
    L = list(map(lambda x: int(x), infile.readline().rstrip().split()))
    L = list(map(lambda x: str(x), TwoWayPartition(L)))
    infile.close()
    newfile = open('par_solution.txt', 'w')
    newfile.write(' '.join(L))
    newfile.close()
    
    


#TwoWayPartition([7, 2, 5, 6, 1, 3, 9, 4, 8])



#Problem
#A partition procedure is an essential part of the Quick Sort algorithm, the subject of one of the following problems. Its main goal is to put the first element of a given array to its proper place in a sorted array. It can be implemented in linear time, by a single scan of a given array. Moreover, it is not hard to come up with an in-place algorithm.
#
#
#
#Given: A positive integer n≤105 and an array A[1..n] of integers from −105 to 105.
#
#Return: A permuted array B[1..n] such that it is a permutation of A and there is an index 1≤q≤n such that B[i]≤A[1] for all 1≤i≤q−1, B[q]=A[1], and B[i]>A[1] for all q+1≤i≤n.
#
#Sample Dataset
#9
#7 2 5 6 1 3 9 4 8
#Sample Output
#5 6 3 4 1 2 7 9 8