# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 00:56:57 2019

@author: Richard
"""



s = '''1 2
2 3
6 3
5 6
2 5
2 4
4 1'''


L = s.plit('\n')





def Degree(L):
    D = {}
    for pair in L:
        pair = pair.split()
        for i in pair:
            if i in D:
                D[i] +=1
            else:
                D[i] = 1
            
            








1 2
2 3
6 3
5 6
2 5
2 4
4 1