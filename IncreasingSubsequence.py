# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 23:25:54 2019

@author: Richard
"""



def LongestIncreasingSubsequence(X):
    """
    (list) -> list
    Returns the Longest Increasing Subsequence in the Given List/Array
    """
    
    
    N = len(X)
    
    # initialize lists or size N and N+1 respectively
      
    # M[j] — stores the index k of the smallest value X[k] such that there is
    #an increasing subsequence of length j ending at X[k] on the range k ≤ i.
    # P[k] — stores the index of the predecessor of X[k] in the longest
    # increasing subsequence ending at X[k].
    P = [0] * N
    M = [0] * (N+1)
    # L is the length of the longest increasing subsequence
    L = 0
    
    for i in range(N):
       # use binary search to find the largest positive j ≤ L such that X[M[j]] <= X[i]
       lo = 1
       hi = L
       while lo <= hi:
           mid = (lo+hi)//2
           if (X[M[mid]] < X[i]):
               lo = mid+1
           else:
               hi = mid-1
       # After searching, lo is 1 greater than the
       # length of the longest prefix of X[i]
       newL = lo
       P[i] = M[newL-1]
       M[newL] = i
 
       # The predecessor of X[i] is the last index of 
       # the subsequence of length newL-1
       if (newL > L):
           L = newL
    # Reconstruct the longest increasing subsequence
    S = []
    k = M[L]
    for i in range(L-1, -1, -1):
        S.append(X[k])
        k = P[k]
    return S[::-1]