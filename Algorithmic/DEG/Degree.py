# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 00:56:57 2019

@author: Richard
"""



def Degree(L):
    '''
    (list) -> None
    Take a list of interctions between nodes labeled 1, 2, ... and 
    print the  degree of each node in increasing node order   
    '''
    
    
    D = {}
    for pair in L:
        pair = pair.split()
        for i in pair:
            if int(i) in D:
                D[int(i)] +=1
            else:
                D[int(i)] = 1
            
    nodes = list(D.keys())
    nodes.sort()
    degree = [str(D[i]) for i in nodes]
    print(' '.join(degree))
    


def SolvePb(pbfile):
    '''
    (file) -> None
    Take input file with node 
    interactions and print the degree of each node
    '''
    
    infile = open(pbfile)
    # skip number of nodes and edges
    infile.readline()
    L = infile.read().rstrip().split('\n')
    infile.close()
    Degree(L)
    

