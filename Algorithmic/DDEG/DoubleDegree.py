# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 01:38:36 2019

@author: Richard
"""



# use a dict to store connections:
# D[i] = {j, k, m}

# iterate through sorted node lists:
for each neighbor, count edges (len(D[j])) --> sum for i



def DoubleDegree(PbFile):
    '''
    
    
    '''
    
    
    # extract data from file
    
    NumNodes, NumEdges = list(map(lambda x: int(x), infile.readline().rstrip().split()))
    
    
    # nodes may be missing from list: add 0 connections
    