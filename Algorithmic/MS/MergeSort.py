# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 17:29:41 2019

@author: rjovelin
"""


def MergeSorted(x, a):
    '''
    (list, list) -> list
    Take 2 sorted lists and return a sorted list containing all elements of x and a    
    '''
    C = []
    L = len(x) + len(a)
    
    # take the minimum of each array and compare them, add the minimum to C
    # stop iteration when all elements have been added to C
    # or when no elements to compare in one of the lists
    while len(C) < L and len(a) != 0 and len(x) != 0:
        # check whether mininum element is in a or x
        if x[0] < a[0]:
            MaxFound = False
            # find the value j in x greater than minimum val of a
            # add all values in x up to j
            for i in range(len(x)):
                if x[i] > a[0]:
                    MaxFound = True
                    break
            if MaxFound == True:
                C.extend(x[:i])
                # update x
                x = x[i:]
            else:
                C. extend(x)
                x = []
        elif a[0] <= x[0]:
            MaxFound = False
            for i in range(len(a)):
                if a[i] > x[0]:
                    MaxFound = True
                    break            
            if MaxFound == True:
                C.extend(a[:i])
                a = a[i:]
            else:
                C.extend(a)
                a = []
            
    if len(a) == 0:
        # add all elements in X
        C.extend(x)
    elif len(x) == 0:
        C.extend(a)
        
    return C




def MergeSort(A):
    
    # create even-length lists of 1 element
    # if len(A) is odd, keep last element for the final round of merge
    if len(A) % 2 != 0:
        B = [[i] for i in A[:-1]]
    else:
        B = [[i] for i in A]
    # make a copy of B, this list will be modified at each iteration
    b = B[::]
    
    # make a list to store the sorted list
    C = []
    while len(C) != 1:
        # re-initiate C
        C = []    
        # loop over pairs of sorted lists, merge and collect to list
        for i in range(0, len(b)-1, 2):
            C.append(MergeSorted(b[i], b[i+1]))
        # check length is odd
        if len(C) % 2 != 0 and len(C) != 1:
            # merge 1st 2 elements to get an even list
            c = MergeSorted(C[0], C[1])
            # delete 1st 2 elements of C, keeping ordering intact
            C.pop(1)
            C. pop(0)
            C.insert(0, c)
        # start this process again on the pairs in the new list
        b = C[::]
                
    # check if len(A) is odd
    if len(A) % 2 != 0:
        # need 1 final round of merge 
        D = MergeSorted(C[0], [A[-1]])
    else:
        D = C[0]
    return D


def SolvePb(PbFile):
    infile = open(PbFile)
    infile.readline()
    L = list(map(lambda x: int(x), infile.readline().rstrip().split()))
    infile.close()
    D = MergeSort(L)
    newfile = open('mergeSort_solution.txt', 'w')
    D = list(map(lambda x: str(x), D))
    newfile.write(' '.join(D))
    newfile.close()
    
