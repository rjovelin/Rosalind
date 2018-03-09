# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:41:58 2015

@author: Richard
"""

def occurence_index(S):
    '''
    (str) -> list of tuples 
    Return a list of tuples containing the characters of string S in their
    order of appearance in S with the index of occurence of of each character
    >>> occurence_index('TCCTCTAT')
    [('T', 1), ('C', 1), ('C', 2), ('T', 2), ('C', 3), ('T', 3), ('A', 1), ('T, 4)]
    '''
    
    # create a dictionnary to store the characters: count in which count is
    # updated each time the character is found when traversing S
    counter = {}
    
    # create list of tuples
    positions = []
    
    # loop over S
    for symbol in S:
        if symbol in counter:
            positions.append((symbol, counter[symbol] + 1))
            counter[symbol] += 1
        else:
            positions.append((symbol, 1))
            counter[symbol] = 1
    return positions
    
 
def BWmatching(input_file):
    '''
    (file) -> list
    Given a string BWT(Text), followed by a collection of strings Patterns, 
    return a list of integers, where the i-th integer corresponds to the
    number of substring matches of the i-th member of Patterns in Text.
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get the BWText
    BWtext = infile.readline().rstrip()
    # get the patterns
    patterns = infile.readline().rstrip().split()
    # close file
    infile.close()
    
    # get last column of the BW matrix
    last = [i for i in BWtext]    
    
    # get first column of the BW matrix
    first = last[:]
    first.sort()
    
    # create lists with tuples (character, index of character occurence )
    # C1: first column, C2, last column
    C1 = occurence_index(''.join(first))
    C2 = occurence_index(BWtext)
    
    # create a dictionnary to store the tuples of C1: tuples of C2 pairs
    matrix = {}
    for i in range(len(C1)):
        matrix[C1[i]] = C2[i]
    
    # create a list to store the matching counts
    matching = []
    
    # loop over the patterns:
    for pattern in patterns:
        # search the occurence of the last character in first column
        # create a list of keys to probe the matrix 
        keys = []
        # set i the index of the last character in pattern
        i = -1
        # get the keys if character in first column matches last character in pattern
        for pair in matrix:
            if pair[0] == pattern[-1]:
                keys.append(pair)
        # ipdate i with the index of the next character        
        i -= 1
        
        # loop until all the characters of pattern have been examined
        while i >= -len(pattern):
            # create a list of new keys 
            new_keys = []
            # populate new_keys if the character in last column with same index
            # as key corresponds to next character in pattern
            for pair in keys:
                if matrix[pair][0] == pattern[i]:
                    new_keys.append(matrix[pair])
            # update keys as new_keys
            keys = new_keys[:]
            # update i to get the next character
            i -= 1
        
        # add the number of characters in keys to the list of pattern matching
        # len(keys) corresponds to the number of pattern matches when loop exists
        # because at this point keys store the occurence index of the first character
        # in pattern for which the entire pattern is found in text
        matching.append(len(keys))
    
    for count in matching:
        print(count, end = ' ')
        
    
    
    
    
    