# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 17:47:41 2015

@author: Richard
"""

def BW_transform(input_file):
    '''
    (file) -> string
    Given a string Text, return the Burrows-Wheeler Transform of text
    '''
    
    
    # open file for reading
    infile = open(input_file, 'r')
    text = infile.readline().rstrip()
    # close file after reading
    infile.close()
    
    #create a list to store the different strings
    BWA = []
    
    # create string by cycling through Text
    BWA.append(text)
    
    # get the number of permutations to cycle through text
    n = len(text) -1
    
    # cycle through text by adding the last letter at the beginning
    # populate list with new string
    while n != 0:
        text = text[-1] + text[:-1]
        BWA.append(text)
        n -= 1
    
    # sort the string by lexicographic order
    BWA.sort()
    
    # the BW transform is the last column of the cyclic strings
    transform = ''
    for word in BWA:
        transform += word[-1]
    
    return transform
       
    
    
    
    
    