# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 18:13:46 2015

@author: Richard
"""
import itertools


def lexicographic_string_order(input_file, outputfile):
    '''
    (file) -> list
    Given a permutation of at most 12 symbols defining an ordered alphabet A
    and a positive integer n (n≤4) in input_file, return all strings of
    length at most n formed from A, ordered lexicographically.
    Alphabet order is based on the order in which the symbols are given in A    
    '''
    
    # open file
    infile = open(input_file, 'r')
    alphabet = ''.join(infile.readline().rstrip().split())
    n = int(infile.readline().rstrip())
    infile.close()
    
    # make a list of all possible words obtained from alphabet
    words = []
    
    # find all possible words of length 1 to n
    for i in range(1, n+1):
        grammar = []
        vocabulary = itertools.product(alphabet, repeat = i)
        # add all words in vocabulary to the list
        for item in vocabulary:
            words.append(''.join(item))
            
    # create a function used to sort the strings
    lexico = lambda word: [alphabet.index(letter) for letter in word]
    
    # sort the list of words using the lexico function
    words.sort(key = lexico)
    
    # open file for writing
    newfile = open(outputfile, 'w')
    for word in words:
        newfile.write(word + '\n')
    newfile.close()
    
    
    