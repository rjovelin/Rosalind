# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 21:12:34 2015

@author: Richard
"""

def word_counter(input_file):
    '''
    (file) -> None
    Given a string s of length at most 10000 letters in input file,
    return how many times any word occurred in string.
    Each letter case (upper or lower) in word matters.
    Lines in output can be in any order.
    '''
    
    # open file
    infile = open(input_file, 'r')
    # grab the text as string, get rif od the new line, and split the words
    S = infile.read().rstrip().split()
    # close file
    infile.close()
    
    # make a dictionnary to hold the word : count pairs
    word_counts = {}
    for word in S:
        if word not in word_counts:
            word_counts[word] = S.count(word)
    # print to screen
    for word, count in word_counts.items():
        print(word, count)
    
    