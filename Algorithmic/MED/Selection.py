# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 22:17:40 2019

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
    return A[:j], A[j:m], A[m:]


def Selection(S, k):
    
    '''
    (list, int) -> int
    Take a list S and an integer k and return the kth smallest element of S
    Selection([2, 36, 5, 21, 8, 13, 11, 20, 5, 4, 1], 8) = 13
    '''
    
    # pick a random number v, always 1st element. and split S into 3 lists using ThreeWayPartition(S)
    # SL = elements lower than v, Sv = elements equal to v and SR = elements higher than v
    # always pick 1 st element
    v = S[0]
    SL, Sv, SR = ThreeWayPartition(S)
    
    
    if k <= len(SL):
        return Selection(SL,k)
    elif len(SL) < k <= len(SL) + len(Sv):
        return v
    elif k > len(SL) + len(Sv):
        return Selection(SR, k - len(SL) - len(Sv))
        
    

def SolvePb(PbFile):
    infile = open(PbFile)
    infile.readline()
    L = list(map(lambda x: int(x), infile.readline().strip().split()))
    k = int(infile.readline().rstrip())
    infile.close()
    print(Selection(L, k))
    