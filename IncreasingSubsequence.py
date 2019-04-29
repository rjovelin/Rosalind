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
    



























def SolvePb(PbFile):
    infile = open(PbFile)
    infile.readline()
    X = list(map(lambda x: x.strip(), infile.readline().rstrip().split()))
    print(' '.join(LongestIncreasingSubsequence(X)))
    
    infile.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#=======================================================================
# Author: Isai Damier
# Title: Longest Decreasing Subsequence
# Project: geekviewpoint
# Package: algorithms
#
# Statement:
#   Given a sequence of numbers, find a longest decreasing subsequence.
#
#
#  Time Complexity: O(n^2)
#
# Sample Input: [5,0,3,2,1,8]
# Sample Output: [5,3,2,1]
#
# DEFINITION OF SUBSEQUENCE:
#   A sequence is a particular order in which related objects follow
#   each other (e.g. DNA, Fibonacci). A sub-sequence is a sequence
#   obtained by omitting some of the elements of a larger sequence.
#
#   SEQUENCE       SUBSEQUENCE     OMISSION
#   [3,1,2,5,4]     [1,2]            3,5,4
#   [3,1,2,5,4]     [3,1,5]          2,4
#
#   SEQUENCE       NOT SUBSEQUENCE   REASON
#   [3,1,2,5,4]     [4,2,5]           4 should follow 5
#
# STRATEGY:
#   Illustrating by finding
#   a longest decreasing subsequence of [5,0,3,2,1,8]:
#
#   - Start by finding all subsequences of size 1: [5],[0],[3],[2],[1],[8];
#     each element is its own decreasing subsequence.
#
#   - Since we already have the solutions for the size 1 subsequences,
#     we can use them in solving for the size two subsequences. For
#     instance, we already know that 5 is the smallest element of a
#     decreasing subsequence of size 1, i.e. the subsequence [5].
#     Therefore, all we need to get a subsequence of size 2 is add an
#     element smaller than 5 to [5]: [5,0], [5,3], [5,2], [5,1];
#     [3,2], [3,1], [2,1].
#
#   - Now we use the size 2 solutions to get the size 3 solutions:
#     [5,3,2], [5,3,1], [3,2,1]
#
#   - Then we use the size 3 solutions to get the size 4 solutions:
#     [5,3,2,1]. Since there are no size 5 solutions, we are done.
#
# SUMMARY:
#   Instead of directly solving the big problem, we solved a smaller
#   version and then 'copied and pasted' the solution of the subproblem
#   to find the solution to the big problem. To make the 'copy and paste'
#   part easy, we use a table (i.e. list) to track the subproblems
#   and their solutions. This strategy as a whole is called Dynamic
#   Programming. The tabling part is known as memoization, which means
#   writing memo.
#
#   To recognize whether you can use dynamic programming on a problem,
#   look for the following two traits: optimal substructures and
#   overlapping subproblems.
#
#   Optimal Substructures: the ability to 'copy and paste' the solution
#     of a subproblem plus an additional trivial amount of work so to
#     solve a larger problem. For example, we were able to use [5,3]
#     itself an optimal solution to the problem [5,0,3] to get [5,3,2]
#     as an optimal solution to the problem [5,0,3,2].
#
#   Overlapping Subproblems: Okay. So in our approach the solution grew
#     from left to right: [5] to [5,3] to [5,3,2] etc. But in reality
#     we could have solved the problem using recursion trees so that
#     for example [5,3] could be reached either from [5] or from [3].
#     That wouldn't really be a problem except we would be solving for
#     [5,3] more than once. Any time a recursive solution would lead to
#     such overlaps, the bet is dynamic programming is the way to go.
#
#          [5]                 [3]
#         / | \               / | \
#        /  |  \             /  |  \
#       /   |   \           /   |   \
#   [5,0] [5,3] [5,2]   [5,3] [3,2] [3,1]
#
# NOTE:
# Dynamic Programming = Optimal Substructures + Overlapping Subproblems
# Divide and Conquer = Optimal Substructures - Overlapping Subproblems
#   see merge sort: http://www.geekviewpoint.com/python/sorting/mergesort
#
# Alternate coding: Not really much difference here, just another code
#   that some readers will find more intuitive:



def LongestDecreasingSubSequence(A):
    '''
    
    
    '''
    
    
# STRATEGY:
#   Illustrating by finding
#   a longest decreasing subsequence of [5,0,3,2,1,8]:
#
#   - Start by finding all subsequences of size 1: [5],[0],[3],[2],[1],[8];
#     each element is its own decreasing subsequence.
#
#   - Since we already have the solutions for the size 1 subsequences,
#     we can use them in solving for the size two subsequences. For
#     instance, we already know that 5 is the smallest element of a
#     decreasing subsequence of size 1, i.e. the subsequence [5].
#     Therefore, all we need to get a subsequence of size 2 is add an
#     element smaller than 5 to [5]: [5,0], [5,3], [5,2], [5,1];
#     [3,2], [3,1], [2,1].
#
#   - Now we use the size 2 solutions to get the size 3 solutions:
#     [5,3,2], [5,3,1], [3,2,1]
#
#   - Then we use the size 3 solutions to get the size 4 solutions:
#     [5,3,2,1]. Since there are no size 5 solutions, we are done.


      # use a list to store 


### NEED TO RECORD POSITION IN THE ARRAY WHERE TO START LOOPING OVER
### ALSO LOOP OVER ALL VALUES OF D[LS]
### BUT NEED TO REMEMBER THE INDEX OF THE LAST ELEMENT OF THE SUBSEQUENCE TO START SEARCHING FROM THERE      
      
      
    # start by finding all subsequences of size 1
    Seqs = {i:[A[i]] for i in range(len(A))}
    # set the length of the longest subsequence 
    LS = 0
    newLS = 1
    # create a dict to store the subsequence of given size
    D = {}
    D[LS] = Seqs
    
    while newLS > LS:
        LS = list(D.keys())[0]    
        #print(LS, newLS)
        for i in D[LS]:
            for j in range(i+1, len(A)):
                #print('i', i, 'j', j, A[j], D[LS][i][-1])
                if A[j] < D[LS][i][-1]:
                    #print(LS, j, i, 'found a smaller value', A[j], D[LS][i][-1])
                    newLS = LS + 1
                    if newLS not in D:
                        D[newLS] = {}
                    D[newLS][j] = D[LS][i] + [A[j]]
        
        #print(D)    
        if newLS > LS:
            #print('found greater')
            #print('LS', LS, 'newLS', newLS)
            del D[LS]
        #print(D)
    return D


#    # start by finding all subsequences of size 1
#    Seqs = {i:[A[i]] for i in range(len(A))}
#    # set the length of the longest subsequence 
#    LS = 0
#    newLS = 1
#    # create a dict to store the subsequence of given size
#    D = {}
#    D[LS] = Seqs
#    
#    while newLS > LS:
#        LS = list(D.keys())[0]    
#        print(LS, newLS)
#        for i in D[LS]:
#            # use binary search to find a lower value
#            L, R = i+1, n-1
#            while L <= R:
#                m = math.floor((L + R) / 2)
#                if A[m] < D[LS][i][-1]:
#                    # 
#                    
#                    
#                # check if target T is on the right or the left of m
#        if A[m] < T:
#            # T on the right of m
#            # set L and R to create a new array on the right
#            L = m + 1
#        elif A[m] > T:
#            # T on the left of m, set L and R to create a new array on the left
#            R = m - 1
#        elif A[m] == T:
#            # target found
#            return m
#    # if L > R the entire array has been searched and target is not found
#    return -1
#
#
#
#
#
#
#            
#            
#            
#            for j in range(i+1, len(A)):
#                if A[j] < D[LS][i][-1]:
#                    print(LS, j, i, 'found a smaller value', A[j], D[LS][i][-1])
#                    newLS = LS + 1
#                    if newLS not in D:
#                        D[newLS] = {}
#                    D[newLS][j] = D[LS][i] + [A[j]]
#        #print(D)    
#        if newLS > LS:
#            #print('found greater')
#            #print('LS', LS, 'newLS', newLS)
#            del D[LS]
#        #print(D)
#    return D











def SolveOb(PbFile):
    
    infile = open(PbFile)
    infile.readline()
    L = list(map(lambda x: str(x), infile.readline().rstrip().split()))
    infile.close()
    LDS = LongestDecreasingSubSequence(L)
    print(' '.join(list(list(LDS.values())[0].values())[0]))







#    m = [1] * len( A )
#
#    for x in range(len(A)):
#        for y in range(x):
#            print(x, y, m[x], m[y], A[x], A[y])
#            if m[y] >= m[x] and A[y] > A[x]:
#                m[x]+=1
#
#    max_value = max(m)
#
#    result = []
#    for i in range(len(m)-1,-1,-1):
#        if max_value == m[i]:
#            result.append(A[i])
#            max_value-=1
#
#    result.reverse()
#    return result



#=======================================================================
  
# def LDS( A ):
#  m = [0] * len( A ) # starting with m = [1] * len( A ) is not necessary
#  for x in range( len( A ) - 2, -1, -1 ):
#    for y in range( len( A ) - 1, x, -1 ):
#      if m[x] <= m[y] and A[x] > A[y]:
#        m[x] = m[y] + 1 # or use m[x]+=1
# 
#  #===================================================================
#  # Use the following snippet or the one line below to get max_value
#  # max_value=m[0]
#  # for i in range(m):
#  #  if max_value < m[i]:
#  #    max_value = m[i]
#  #===================================================================
#  max_value = max( m )
# 
#  result = []
#  for i in range( len( m ) ):
#    if max_value == m[i]:
#      result.append( A[i] )
#      max_value -= 1
# 
#  return result    
    
    
    
    
    
    
    