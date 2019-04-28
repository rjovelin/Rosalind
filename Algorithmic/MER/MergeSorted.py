# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 00:37:48 2019

@author: Richard
"""



def MergeSorted(x, a):
    '''
    (list, list) -> list
    Take 2 sorted lists and return a sorted list containing all elements of x and a    
    '''
    C = []
    # make copies of arrayes
    A, X = a[::], x[::]
    # take the minimum of each array and compare them, add the minimum to C
    # stop iteration when all elements have been added to C
    # or when no elements to compare in one of the lists
    while len(C) < len(x) + len(a) and len(A) != 0 and len(X) != 0:
        if min(X) < min(A):
            C.append(X.pop(X.index(min(X))))
        elif min(A) <= min(X):
            C.append(A.pop(A.index(min(A))))
    if len(A) == 0:
        # ad all elements in X
        C.extend(X)
    elif len(X) == 0:
        C.extend(A)
        
    return C



def SolvePb(PbFile):
    infile = open(PbFile)
    infile.readline()
    A = list(map(lambda x: int(x), infile.readline().rstrip().split()))
    infile.readline()
    B = list(map(lambda x: int(x), infile.readline().rstrip().split()))
    infile.close()
    
    C = MergeSorted(A, B)
    C = list(map(lambda x: str(x), C))
    newfile = open('mer_solution.txt', 'w')
    newfile.write(' '.join(C))
    newfile.close()

