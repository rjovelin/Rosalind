# -*- coding: utf-8 -*-
"""
Created on Thu May  2 17:58:15 2019

@author: rjovelin
"""



def MajorityElement(A):
    '''
    (list, list) -> list
    Take 2 sorted lists and return a sorted list containing all elements of x and a    
    '''
    
    D = {}
    Max = 0
    T = 0
    j = ''
    for i in A:
        if i in D:
            D[i] += 1
        else: D[i] = 1
        T += 1
        if D[i] >= Max:
            Max = D[i]
            j = i
    if Max > T / 2:
        return j
    else:
        return -1

    
    
def SolvePb(PbFile):
    infile = open(PbFile)
    infile.readline()
    L = []
    for line in infile:
        if line.rstrip() != '':
            L.append(list(map(lambda x: int(x), line.rstrip().split())))
    infile.close()
    S = ' '.join([str(MajorityElement(i)) for i in L])
    print(S)
        