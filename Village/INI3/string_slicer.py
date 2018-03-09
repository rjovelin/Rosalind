# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 19:35:46 2015

@author: Richard
"""

def string_slicer(input_file):
    '''
    (file) -> 
    Get string s of lengh <= 200 and integers a, b, c and d
    Return the slice of s from indices a through b and c through d
    (with space in between), inclusively.
    '''
    
    # open file
    infile = open(input_file, 'r')
    s = infile.readline().rstrip()
    num = infile.readline().rstrip().split()
    a = int(num[0])
    b = int(num[1])
    c = int(num[2])
    d = int(num[3])
    infile.close()
    
    word1 = s[a:b+1]
    word2 = s[c:d+1]
    return word1 + ' ' + word2
    
    
    
